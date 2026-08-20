// split_polygons.h —— 用折面切割多边形组,按(来源多边形, 切割片)输出边界环
//
// 设计/验证文档:docs/polygon-split-by-folded-surface.md
// 核心链路:①建面(每个原始面 = 一个多边形,相邻多边形共享顶点/边)
//           → ②triangulate_faces + visitor 逐面标记来源(禁用连通分量记来源:
//             相邻多边形有共用边,会并成一个分量,来源信息丢失)
//           → ④PMP::split + visitor 标记传播
//           → ⑤连通分组:以"两侧来源不同的边"为 barrier 求连通分量
//             → pid = (来源多边形 × 切割片) 组号,每个多边形各自分成多片
//           → ⑥组边界环提取:face(h) 属于本组且(对侧无面或对侧异组)的半边,
//             沿 next 链成环 —— 覆盖撕开边与相邻多边形共享边两种组边界
//
// 依赖:CGAL 6.2(PMP)、GMP/MPFR(EPECK 精确核)、Boost
// 注意:全程与文档一致使用精确核 EPECK;点坐标为 K::FT 精确有理数,
//       交给下游 double 管线时逐坐标 CGAL::to_double()。

#pragma once

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/clip.h>                 // PMP::split
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>   // does_self_intersect
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>     // degenerate_faces
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/face_graph_utils.h> // Default_visitor
#include <CGAL/boost/graph/helpers.h>                          // is_valid_polygon_mesh
#include <CGAL/boost/graph/iterator.h>                         // halfedges_around_face

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace poly_split {

namespace PMP = CGAL::Polygon_mesh_processing;

// ── 输出:一片 = 来源多边形编号 + 切割组编号 + 边界环点列 ────────────────
template <class K>
struct PolygonPiece {
  std::size_t source = 0;   // 来源:②标记的原始面号(= 输入多边形号)
  std::size_t piece  = 0;   // 片号:⑤ 的组号(来源多边形 × 切割片)
  std::vector<typename K::Point_3> ring;  // 该片外边界环(带洞片会另出环,见文档⑥)
};

template <class K>
using Polygon_mesh = CGAL::Surface_mesh<typename K::Point_3>;

template <class K>
struct Split_result {
  Polygon_mesh<K> mesh;                    // 撕开后的 tm(调试/校验用)
  std::size_t component_count = 0;         // ⑤ 的组数(无洞片时应 == pieces.size())
  std::size_t unmarked_faces   = 0;        // 切割后仍为 -1 来源标记的面数(应 == 0)
  std::vector<PolygonPiece<K>> pieces;
};

// ── 三角化标记 visitor(②)─────────────────────────────────────────────
// triangulate_faces 把一个面 f 拆成若干三角形:其中一个三角形复用 f 的描述符
// (来源标记天然保留),其余新面逐个触发 after_subface_created(f_new) → 继承 f 的
// 来源标记(before_subface_creations(f_old) 记录)。钩子见 triangulate_faces.h:55-57,
// quad 快速路径(triangulate_faces.h:288-321)与 hole-filling 路径(143 行起)均回调;
// quad 路径对复用原面也回调 after_subface_created(写同值,幂等无害)。
//
// 与 corefine 不同:triangulate_faces 只作用于单一网格、钩子无 mesh 参数 → 无需 watch 守卫。
// visitor 会被按值拷贝(triangulate_faces.h:430-431);Surface_mesh::Property_map 是
// 共享底层存储的轻句柄,拷贝安全、写入互通。
template <class TagMap, class TriangleMesh>
struct Tri_source_tagger : PMP::Triangulate_faces::Default_visitor<TriangleMesh>
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  TagMap fid;
  int cur = -1;

  Tri_source_tagger() = default;
  explicit Tri_source_tagger(TagMap m) : fid(std::move(m)) {}

  void before_subface_creations(face_descriptor f_old)
  { cur = get(fid, f_old); }
  void after_subface_created(face_descriptor f_new)
  { put(fid, f_new, cur); }
};

// ── 标记传播 visitor(④)────────────────────────────────────────────────
// corefine 切分一个面时:第一个子面复用原面描述符(标记天然保留,face_graph_utils.h:614-619),
// 其余每个新面触发 after_subface_created → 继承被切面(before_subface_creations 记录)的标记。
//
// 关键:corefine 会同时细化 tm 和 splitter 两个网格,同一 visitor 对两个网格都会触发;
// 钩子的 const TriangleMesh& 参数即用于区分。必须用 watch 指针守卫,只对被标记网格操作,
// 否则会把 splitter 的面索引写进 tm 的属性表(轻则写花标记,重则数组越界断言)。
//
// visitor 会被按值拷贝(clip.h:1337);Surface_mesh::Property_map 是共享底层存储的轻句柄,
// 拷贝安全、写入互通。after_face_copy 覆盖共面复制路径(保险,通常不触发)。
template <class TagMap, class TriangleMesh>
struct Tag_propagator : PMP::Corefinement::Default_visitor<TriangleMesh>
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  TagMap fid;
  const TriangleMesh* watch = nullptr;   // 只作用于该网格(= 被标记的 tm)
  int cur = -1;

  Tag_propagator() = default;
  Tag_propagator(TagMap m, const TriangleMesh& w) : fid(std::move(m)), watch(&w) {}

  void before_subface_creations(face_descriptor f_old, const TriangleMesh& m)
  { if (&m == watch) cur = get(fid, f_old); }
  void after_subface_created(face_descriptor f_new, const TriangleMesh& m)
  { if (&m == watch) put(fid, f_new, cur); }
  void after_face_copy(face_descriptor f_old, const TriangleMesh& m1,
                       face_descriptor f_new, const TriangleMesh& /*m2*/)
  { if (&m1 == watch) put(fid, f_new, get(fid, f_old)); }
};

// ── 建 mesh:每个点列环各建一个面(形态 1/2:单多边形 / 批量不相交多边形)──
// 多边形必须简单(不自交)、互不相交。注意:此入口给每个多边形建独立顶点,
// 相邻(共边)多边形不能用 —— 共边处顶点重复、边不共享,corefine 会遇共面退化;
// 相邻多边形请手工建面共享顶点后调 split_mesh(见 tests 用例 5/6)。
template <class K>
Polygon_mesh<K> build_mesh(const std::vector<std::vector<typename K::Point_3>>& polygons)
{
  Polygon_mesh<K> m;
  std::vector<typename Polygon_mesh<K>::Vertex_index> vds;
  for (const auto& poly : polygons) {
    vds.clear();
    for (const auto& p : poly) vds.push_back(m.add_vertex(p));
    m.add_face(vds);
  }
  return m;
}

// 折面建 splitter:若干平面片(点列环)→ 建面 + 三角化。
// 前提(文档③):无自交、完全横跨被切多边形(建议向四周延伸出包围盒)。
// 需要共享边的折面(如 V 形两坡共屋脊)请自行建面共享顶点,再手动 triangulate_faces。
template <class K>
Polygon_mesh<K> build_splitter(const std::vector<std::vector<typename K::Point_3>>& patches)
{
  auto s = build_mesh<K>(patches);
  PMP::triangulate_faces(s);
  return s;
}

// ── split 前置体检(诊断"看起来正常的 mesh 也挂")──────────────────────
// PMP::split 不消费 throw_on_self_intersection(§3④):前置被违反时它不抛
// 异常、直接未定义行为(常表现为崩在 corefine 深处)。此函数逐条检查文档前置
// 与静默失效模式,返回违规清单(空 = 全部通过,可安全进入 split)。
// 检查对象语义与管线一致:
//   - tm:管线内会先三角化(②),多边形面是合法输入 → 在"三角化后的副本"上
//     体检(与 split/corefine 实际看到的网格一致);三角化后仍非三角 = 环退化;
//   - splitter:调用者必须传已三角化的(build_splitter 已保证)→ 直接要求全三角。
// 各面通过全三角检查后才做退化/自交检查(两者假设三角形面:
// degenerate_faces shape_predicates.h:272、does_self_intersect
// self_intersections.h:774)。自交最贵,放最后。
// 依据:split 内部即 corefine + split_along_edges(clip.h:1345-1350),
// 改调 PMP::corefine 躲不开这些前置,且 corefine 不撕开、撕开无公开 API。
template <class K>
std::vector<std::string> preflight_check(const Polygon_mesh<K>& tm,
                                         const Polygon_mesh<K>& splitter)
{
  typedef Polygon_mesh<K> Mesh;
  std::vector<std::string> problems;

  auto all_triangles = [](const Mesh& m) {
    for (auto f : faces(m)) {
      std::size_t d = 0;
      for (auto h : halfedges_around_face(halfedge(f, m), m)) {
        if (++d > 3) break;
      }
      if (d != 3) return false;
    }
    return true;
  };
  auto geom_checks = [&](const Mesh& m, const char* name) {
    std::vector<typename Mesh::Face_index> deg;
    PMP::degenerate_faces(m, std::back_inserter(deg));
    if (!deg.empty())
      problems.push_back(std::string(name) + ": " + std::to_string(deg.size()) +
                         " 个退化面(零面积)");
    if (PMP::does_self_intersect(m))
      problems.push_back(std::string(name) + ": 自交(split 不抛异常,直接未定义行为)");
  };

  if (!CGAL::is_valid_polygon_mesh(tm, false))
    problems.push_back("tm: 网格结构非法");
  else {
    Mesh tm_copy(tm);                       // 与管线同路径:先三角化再体检
    PMP::triangulate_faces(tm_copy);
    if (!all_triangles(tm_copy))
      problems.push_back("tm: 存在无法三角化的面(顶点共线等退化环)");
    else
      geom_checks(tm_copy, "tm");
  }

  if (!CGAL::is_valid_polygon_mesh(splitter, false))
    problems.push_back("splitter: 网格结构非法");
  else if (!all_triangles(splitter))
    problems.push_back("splitter: 存在非三角面(splitter 必须已三角化,如经 build_splitter)");
  else
    geom_checks(splitter, "splitter");

  return problems;
}

// ── 核心:标记 → 切割 → 分组 → 取环 ────────────────────────────────────
// tm 按值传入(函数内修改);splitter 传常引用,内部复制(split 会细化 splitter 副本)。
// 前提:tm 各面为简单多边形环(可非凸),每个原始面 = 一个输入多边形(相邻多边形
// 必须共享顶点/边地建入同一 tm);tm 与 splitter 同 kernel、均无自交。
template <class K>
Split_result<K> split_mesh(Polygon_mesh<K> tm, const Polygon_mesh<K>& splitter)
{
  typedef Polygon_mesh<K>               Mesh;
  typedef typename Mesh::Face_index     Face_index;
  typedef typename Mesh::Edge_index     Edge_index;
  typedef typename Mesh::Halfedge_index halfedge_descriptor;

  // ② 来源标记 + 三角化:先给每个原始面写自己的面号(= 多边形号),再用带
  //    visitor 的 triangulate_faces 三角化 —— 每个新三角形继承所在原始面的来源。
  //    (不用"先三角化再连通分量记来源":相邻多边形共边会并成一个分量。)
  auto fid = tm.add_property_map<Face_index, int>("f:source", -1).first;
  {
    std::size_t poly_id = 0;
    for (auto f : faces(tm)) put(fid, f, static_cast<int>(poly_id++));
  }
  PMP::triangulate_faces(tm,
      PMP::parameters::visitor(Tri_source_tagger<decltype(fid), Mesh>(fid)));

  // ④ 切割 + 标记传播
  //    注:两网格 split() 不消费 throw_on_self_intersection(clip.h:1325-1347 未提取、
  //    未转发给 corefine),故不传;visitor 的工厂/链式写法依据 clip.h:1346。
  Tag_propagator<decltype(fid), Mesh> tagger(fid, tm);   // 守卫:只传播 tm 的标记
  Mesh splitter_copy(splitter);
  PMP::split(tm, splitter_copy,
             CGAL::parameters::vertex_point_map(get(CGAL::vertex_point, tm))
                              .visitor(tagger));

  // 兜底:统计未被标记传播覆盖的新面(应恒为 0;非 0 说明存在未挂钩子的建面路径)
  Split_result<K> result;
  for (auto f : faces(tm))
    if (get(fid, f) == -1) ++result.unmarked_faces;

  // ⑤ 分组:两侧来源不同的边(= 相邻多边形的共享边)设为约束边,connected_components
  //    不跨约束边扩散(connected_components.h:149-152)→ 撕开分量内再按多边形隔开,
  //    pid = (来源多边形 × 切割片) 组号:每个多边形各自分成多片。
  auto pid = tm.add_property_map<Face_index, std::size_t>("f:piece", 0).first;
  auto ecm = tm.add_property_map<Edge_index, bool>("e:source_border", false).first;
  for (auto e : edges(tm)) {
    auto h  = halfedge(e, tm);
    auto f1 = face(h, tm), f2 = face(opposite(h, tm), tm);
    if (f1 != Mesh::null_face() && f2 != Mesh::null_face() &&
        get(fid, f1) != get(fid, f2))
      put(ecm, e, true);
  }
  result.component_count =
      PMP::connected_components(tm, pid, PMP::parameters::edge_is_constrained_map(ecm));

  // ⑥ 组边界环:半边 h 是组 g 的边界半边 ⟺ face(h) ∈ g 且(对侧无面[撕开边/网格
  //    外边界] 或 对侧面组号 ≠ g[相邻多边形共享边])。
  //    走法(面集边界标准旋转,不能沿 next() 链 —— next 不换面,跨不过组内共享边):
  //    从边界半边 h 的终点 v 取 g=opposite(h)(自 v 出发、指向 h 起点,必在组外),
  //    绕 v 逐步旋转 g=next(opposite(g)),直到 g 自身是组边界半边(面在组内且对侧
  //    不在组内 —— 只判"面在组内"会停在组内对角线上),它就是环上紧随 h 的边界边。
  auto vpm = get(CGAL::vertex_point, tm);
  std::vector<char> visited(tm.number_of_halfedges(), 0);
  auto group_border = [&](halfedge_descriptor hh, std::size_t grp) {
    auto fo = face(opposite(hh, tm), tm);
    return fo == Mesh::null_face() || get(pid, fo) != grp;
  };
  auto next_border = [&](halfedge_descriptor hh, std::size_t grp) {
    halfedge_descriptor g = opposite(hh, tm);      // 绕 target(hh) 旋转
    for (;;) {
      auto fg  = face(g, tm);
      auto fgo = face(opposite(g, tm), tm);
      if (fg != Mesh::null_face() && get(pid, fg) == grp &&
          (fgo == Mesh::null_face() || get(pid, fgo) != grp))
        return g;
      g = next(opposite(g, tm), tm);               // 旋转一步,仍自 target(hh) 出发
    }
  };
  for (auto h : halfedges(tm)) {
    auto f = face(h, tm);
    if (f == Mesh::null_face()) continue;
    std::size_t grp = get(pid, f);
    if (!group_border(h, grp) || visited[h.idx()]) continue;

    PolygonPiece<K> pc;
    pc.source = static_cast<std::size_t>(get(fid, f));
    pc.piece  = grp;
    auto cur = h;
    do {
      visited[cur.idx()] = 1;
      pc.ring.push_back(get(vpm, target(cur, tm)));
      cur = next_border(cur, grp);
    } while (cur != h);
    result.pieces.push_back(std::move(pc));
  }

  result.mesh = std::move(tm);
  return result;
}

// 便捷入口:多边形点列组 + 折面 → 各片(来源, 边界环)。
// 仅适用于互不相交的多边形(每个多边形一个独立面);相邻多边形见 build_mesh 注释。
template <class K>
Split_result<K> split_polygons(const std::vector<std::vector<typename K::Point_3>>& polygons,
                               const Polygon_mesh<K>& splitter)
{
  return split_mesh<K>(build_mesh<K>(polygons), splitter);
}

} // namespace poly_split
