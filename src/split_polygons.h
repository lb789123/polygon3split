// split_polygons.h —— 用折面切割多边形组,按(来源多边形, 切割片)输出边界环
//
// 设计/验证文档:docs/polygon-split-by-folded-surface.md
// 核心链路:①建面(每个原始面 = 一个多边形,相邻多边形共享顶点/边)
//           → ②三角化 + 逐面标记来源:≥5 顶点面走本文件的平面化 CDT
//             (triangulate_polygon_faces,PMP 洞填充路径对不共面环有缺陷,
//             见该函数注释);quad 交 PMP::triangulate_faces 精确路径 + visitor
//             (禁用连通分量记来源:相邻多边形有共用边,会并成一个分量)
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
#include <CGAL/number_utils.h>                                // to_double(OBJ 导出)
#include <CGAL/Constrained_Delaunay_triangulation_2.h>         // ≥5 边形平面化 CDT
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <array>
#include <cstddef>
#include <cstdlib>      // getenv(调试导出开关)
#include <filesystem>   // create_directories(调试导出目录)
#include <fstream>      // 调试 OBJ 导出
#include <queue>        // CDT 外域洪泛
#include <stdexcept>    // invalid_argument
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

// ── 单环平面化 CDT:Newell 法向投影 + 精确 2D 约束 Delaunay ─────────────
// 环顶点 → n-2 个(环下标)三角形,朝向与环一致;失败返回空。
// 失败 = 环去重后 <5 点(本函数只服务 ≥5;3/4 交还 PMP 的精确路径)、
//        Newell 法向为零(整环共线/退化)、投影后自交(constraint 相交异常)、
//        CDT 退化(维度 <2)或域内三角形投影共线。
template <class K>
std::vector<std::array<std::size_t, 3>>
planar_cdt_ring(const std::vector<typename K::Point_3>& ring)
{
  typedef CGAL::Projection_traits_3<K>                                 PT;
  typedef CGAL::Triangulation_vertex_base_with_info_2<std::size_t, PT> Vb;
  typedef CGAL::Triangulation_face_base_with_info_2<bool, PT>          Fbi;
  typedef CGAL::Constrained_triangulation_face_base_2<PT, Fbi>         Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                 TDS;
  // 环是简单环(投影后不交叉)→ 约束不会相交,可用免相交检查的快速 tag
  typedef CGAL::Constrained_Delaunay_triangulation_2<PT, TDS,
          CGAL::No_constraint_intersection_tag>                        CDT;
  typedef typename K::Vector_3                                         V3;

  const std::size_t n = ring.size();
  if (n < 5) return {};

  // Newell 法向 = Σ p_i × p_{i+1}(精确;近平面简单环 = 2×矢量面积,恒非零)
  V3 N(0, 0, 0);
  for (std::size_t i = 0; i < n; ++i) {
    const auto& a = ring[i];
    const auto& b = ring[(i + 1) % n];
    N = N + V3(a.y() * b.z() - a.z() * b.y(),
               a.z() * b.x() - a.x() * b.z(),
               a.x() * b.y() - a.y() * b.x());
  }
  if (N == CGAL::NULL_VECTOR) return {};

  PT traits(N);
  CDT cdt(traits);
  std::vector<std::pair<typename K::Point_3, std::size_t>> pts;
  pts.reserve(n);
  for (std::size_t i = 0; i < n; ++i) pts.emplace_back(ring[i], i);
  cdt.insert(pts.begin(), pts.end());
  std::vector<typename CDT::Vertex_handle> vh(n);
  for (auto v : cdt.finite_vertex_handles()) vh[v->info()] = v;
  if (cdt.number_of_vertices() != n) return {};   // 有点被合并(重复点;调用方已去重)
  try {
    for (std::size_t i = 0; i < n; ++i)
      if (vh[i] != vh[(i + 1) % n])
        cdt.insert_constraint(vh[i], vh[(i + 1) % n]);
  } catch (const typename CDT::Intersection_of_constraints_exception&) {
    return {};                                     // 投影后自交
  }
  if (cdt.dimension() != 2) return {};             // 全共线等退化

  // 外域洪泛标记:自无限面起跨非约束边扩散,余下的有限面即多边形域内
  for (auto fit = cdt.all_faces_begin(); fit != cdt.all_faces_end(); ++fit)
    fit->info() = false;
  std::queue<typename CDT::Face_handle> fq;
  fq.push(cdt.infinite_vertex()->face());
  while (!fq.empty()) {
    auto fh = fq.front();
    fq.pop();
    if (fh->info()) continue;
    fh->info() = true;
    for (int i = 0; i < 3; ++i)
      if (!cdt.is_constrained(typename CDT::Edge(fh, i)))
        fq.push(fh->neighbor(i));
  }

  // 域内三角形,朝向统一为与 Newell 法向成右手系(与环旋向一致)
  std::vector<std::array<std::size_t, 3>> tris;
  for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
    if (fit->info()) continue;                     // 外域
    const std::size_t a = fit->vertex(0)->info();
    const std::size_t b = fit->vertex(1)->info();
    const std::size_t c = fit->vertex(2)->info();
    const V3 t = CGAL::cross_product(ring[b] - ring[a], ring[c] - ring[a]);
    if (t * N == 0) return {};                     // 域内投影退化
    tris.push_back(t * N > 0 ? std::array<std::size_t, 3>{a, b, c}
                             : std::array<std::size_t, 3>{a, c, b});
  }
  if (tris.size() != n - 2) return {};             // 无洞简单环必须恰为 n-2 片
  return tris;
}

// ── ≥5 顶点面的平面化三角化(② 的前置步骤)──────────────────────────────
// 为什么不把 ≥5 顶点面直接交给 PMP::triangulate_faces:它走洞填充路径
// (triangulate_faces.h:328 → triangulate_hole_polyline),其平面 CDT 子路径的
// 共面闸门以 bbox 的 **z 跨度** 为容差基准 —— triangulate_hole.h:769-771 取
// vertex(0)-vertex(5) 的平方距离,而在 CGAL 6.2 的 Iso_cuboid_3 顶点序里这两点
// 仅相差 z(实测确认)。于是近似水平的多边形(z 跨度≈0,平面数据最常见形态)
// 容差≈0:EPECK 下哪怕 1e-17 的不共面都会被拒,然后【静默】回退到 3D Delaunay
// 洞填充(use_delaunay_triangulation 默认开)—— 那是给曲面洞补片设计的最小
// 二面角搜索,在近平面点集上三角形质量失控(对角线乱穿,视觉即"乱七八糟"),
// 大环更是 O(n³)。另实测:PMP 路径遇环上连续重复顶点会静默不三角化(非三角
// 面直接进 split,未定义行为)。
//
// 本函数做真正的"3D 多边形三角化"(对不共面环 = 投影到 Newell 平面后做平面
// 三角化,割耳思路的 Delaunay 版,三角形质量更优):对每个 ≥5 顶点的面
// (先去连续重复点),在 Newell 法向投影下建精确 2D 约束 Delaunay(PMP 平面
// 子路径同款机制,但无共面闸门),三角形朝向与环一致。3/4 顶点面原样保留:
// PMP 的 quad 路径(triangulate_faces.h:240-324)纯精确算术、无平面性前提、
// 对角线选择最优,继续交给它。输入多边形允许不共面(严格或近似)。
//
// 返回:重建后网格的面按原面分组(groups[k] = 原面 k 拆出的新面,顺序与
// 原面一致);任一环退化时返回【空 vector】且 m 原样不动(全或无)。
// 重建逐顶点复制:相邻多边形的共享顶点/边照常保留,add_face 自动缝合。
template <class K>
std::vector<std::vector<typename Polygon_mesh<K>::Face_index>>
triangulate_polygon_faces(Polygon_mesh<K>& m)
{
  typedef Polygon_mesh<K>               Mesh;
  typedef typename Mesh::Vertex_index   Vertex_index;
  typedef typename Mesh::Face_index     Face_index;
  typedef typename K::Point_3           Point_3;

  // 1. 采集原面环(顶点描述符 + 点),去连续重复点(含首尾相接处)
  auto vpm = get(CGAL::vertex_point, m);
  std::vector<std::vector<Vertex_index>> rings_v;
  std::vector<std::vector<Point_3>>      rings_p;
  for (auto f : faces(m)) {
    std::vector<Vertex_index> rv;
    std::vector<Point_3>      rp;
    for (auto h : halfedges_around_face(halfedge(f, m), m)) {
      auto v = target(h, m);
      if (rp.empty() || rp.back() != get(vpm, v)) {
        rv.push_back(v);
        rp.push_back(get(vpm, v));
      }
    }
    while (rp.size() > 1 && rp.front() == rp.back()) { rv.pop_back(); rp.pop_back(); }
    rings_v.push_back(std::move(rv));
    rings_p.push_back(std::move(rp));
  }

  // 2. ≥5 顶点环:平面化 CDT(任一失败 → 全或无,不动 m)
  std::vector<std::vector<std::array<std::size_t, 3>>> tris;
  for (const auto& rp : rings_p) {
    if (rp.size() >= 5) {
      auto t = planar_cdt_ring<K>(rp);
      if (t.empty()) return {};
      tris.push_back(std::move(t));
    } else {
      if (rp.size() < 3) return {};                // 去重后退化(<3 点不成面)
      tris.emplace_back();                         // 3/4 顶点:原样保留
    }
  }

  // 3. 重建:逐顶点复制(保留共享),逐原面加面(≥5:三角形;3/4:原环面)
  Mesh out;
  std::size_t max_vi = 1;
  for (auto v : vertices(m)) {
    const std::size_t i = v.idx() + 1;             // 去重后可能有编号空洞
    if (i > max_vi) max_vi = i;
  }
  std::vector<Vertex_index> nv(max_vi);
  for (auto v : vertices(m)) nv[v.idx()] = out.add_vertex(get(vpm, v));

  std::vector<std::vector<Face_index>> groups;
  for (std::size_t k = 0; k < rings_v.size(); ++k) {
    std::vector<Face_index> g;
    if (rings_p[k].size() >= 5) {
      g.reserve(tris[k].size());
      for (const auto& t : tris[k]) {
        auto f = out.add_face(nv[rings_v[k][t[0]].idx()],
                              nv[rings_v[k][t[1]].idx()],
                              nv[rings_v[k][t[2]].idx()]);
        if (f == Mesh::null_face()) return {};     // 非流形/退化 → 全或无
        g.push_back(f);
      }
    } else {
      std::vector<Vertex_index> fv;
      for (auto v : rings_v[k]) fv.push_back(nv[v.idx()]);
      auto f = out.add_face(fv);
      if (f == Mesh::null_face()) return {};
      g.push_back(f);
    }
    groups.push_back(std::move(g));
  }
  m = std::move(out);
  return groups;
}

// 折面建 splitter:若干平面片(点列环)→ 建面 + 三角化。
// 需要共享边的折面(如 V 形两坡共屋脊)请自行建面共享顶点,再手动
// triangulate_polygon_faces + PMP::triangulate_faces(与管线②一致)。
template <class K>
Polygon_mesh<K> build_splitter(const std::vector<std::vector<typename K::Point_3>>& patches)
{
  auto s = build_mesh<K>(patches);
  const std::size_t n_faces = s.number_of_faces();
  if (triangulate_polygon_faces<K>(s).size() != n_faces)
    throw std::invalid_argument(
        "build_splitter: 折面片无法三角化(≥5 顶点且退化,或投影后自交)");
  PMP::triangulate_faces(s);   // quad 路径收尾
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
    const std::size_t n_faces = tm_copy.number_of_faces();
    if (triangulate_polygon_faces<K>(tm_copy).size() != n_faces)
      problems.push_back("tm: 存在无法平面化三角化的 ≥5 顶点面(退化或投影自交)");
    PMP::triangulate_faces(tm_copy);        // quad 路径收尾
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

// ── 调试:分阶段 OBJ 导出(可选)────────────────────────────────────────
// split_mesh 内设置环境变量 POLYGON3SPLIT_DEBUG_DIR(输出目录)后,每次调用
// 写出一组 .obj,逐阶段人工检查(Blender / MeshLab 可直接打开):
//   <call>_00_tm_input.obj          建面后的 tm(三角化前,可能含多边形面)
//   <call>_01_tm_triangulated.obj   ② 之后(全三角,每面注释来源编号)
//   <call>_02_splitter.obj          ④ 实际消费的折面副本(已三角化)
//   <call>_03_tm_after_split.obj    ④+⑤ 之后(已撕开,每面注释 来源/片号)
//   <call>_04_pieces.obj            ⑥ 各片边界环(每片一个多边形)
// 坐标统一 CGAL::to_double;OBJ 行用 CRLF。未设置环境变量则零开销不输出。
template <class Mesh, class FaceNote>
void write_mesh_obj(const Mesh& mesh,
                    const std::string& file_path,
                    const FaceNote& face_note)
{
    std::ofstream out(file_path, std::ios::binary);
    if (!out)
    {
        return;
    }
    auto emit_line = [&out](const std::string& text)
    {
        out << text << "\r\n";
    };

    auto vpm = get(CGAL::vertex_point, mesh);
    for (auto v : vertices(mesh))
    {
        const auto& p = get(vpm, v);
        emit_line("v " + std::to_string(CGAL::to_double(p.x())) + " " +
                  std::to_string(CGAL::to_double(p.y())) + " " +
                  std::to_string(CGAL::to_double(p.z())));
    }
    for (auto f : faces(mesh))
    {
        std::string line = "f";
        for (auto h : halfedges_around_face(halfedge(f, mesh), mesh))
        {
            line += " " + std::to_string(target(h, mesh).idx() + 1);
        }
        const std::string note = face_note(f);
        if (!note.empty())
        {
            line += "  # " + note;
        }
        emit_line(line);
    }
}

// 不带每面注释的导出(仅顶点 + 面)
template <class Mesh>
void write_mesh_obj(const Mesh& mesh, const std::string& file_path)
{
    write_mesh_obj(mesh, file_path,
                   [](const typename boost::graph_traits<Mesh>::face_descriptor&)
                   {
                       return std::string();
                   });
}

// 各片边界环导出:每片 = 一个 OBJ 组 + 一个多边形面,注释标明来源/片号
template <class Result>
void write_pieces_obj(const Result& result, const std::string& file_path)
{
    std::ofstream out(file_path, std::ios::binary);
    if (!out)
    {
        return;
    }
    auto emit_line = [&out](const std::string& text)
    {
        out << text << "\r\n";
    };

    std::size_t vertex_offset = 0;
    for (const auto& pc : result.pieces)
    {
        emit_line("g source_" + std::to_string(pc.source) + "_piece_" +
                  std::to_string(pc.piece));
        for (const auto& p : pc.ring)
        {
            emit_line("v " + std::to_string(CGAL::to_double(p.x())) + " " +
                      std::to_string(CGAL::to_double(p.y())) + " " +
                      std::to_string(CGAL::to_double(p.z())));
        }
        std::string line = "f";
        for (std::size_t i = 1; i <= pc.ring.size(); ++i)
        {
            line += " " + std::to_string(vertex_offset + i);
        }
        emit_line(line + "  # ring_size=" + std::to_string(pc.ring.size()));
        vertex_offset += pc.ring.size();
    }
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

  // 调试导出开关:环境变量 POLYGON3SPLIT_DEBUG_DIR 指向输出目录,非空才写 OBJ
  const char* const DEBUG_DIR_ENV_VAR = "POLYGON3SPLIT_DEBUG_DIR";
  const char* debug_dir_env = std::getenv(DEBUG_DIR_ENV_VAR);
  const bool dump_debug_obj = (debug_dir_env != nullptr && debug_dir_env[0] != '\0');
  static std::size_t debug_call_seq = 0;
  std::string debug_prefix;
  if (dump_debug_obj)
  {
    std::filesystem::create_directories(debug_dir_env);
    debug_prefix = std::string(debug_dir_env) + "/" +
                   std::to_string(debug_call_seq++) + "_";
    write_mesh_obj(tm, debug_prefix + "00_tm_input.obj");
  }

  // ② 来源标记 + 三角化:≥5 顶点面先走平面化 CDT(triangulate_polygon_faces,
  //    见其注释:PMP 洞填充路径对不共面环会静默回退 3D Delaunay、质量失控),
  //    新三角形直接按分组写来源;剩余 quad 交给 PMP::triangulate_faces 的精确
  //    专用路径 + visitor 继承来源。
  //    (不用"先三角化再连通分量记来源":相邻多边形共边会并成一个分量。)
  const std::size_t n_orig_faces = tm.number_of_faces();
  auto groups = triangulate_polygon_faces<K>(tm);
  if (groups.size() != n_orig_faces)
    throw std::invalid_argument(
        "split_mesh: 存在无法三角化的多边形面(≥5 顶点且退化,或投影后自交)");
  auto fid = tm.add_property_map<Face_index, int>("f:source", -1).first;
  for (std::size_t poly_id = 0; poly_id < groups.size(); ++poly_id)
    for (auto f : groups[poly_id])
      put(fid, f, static_cast<int>(poly_id));
  PMP::triangulate_faces(tm,
      PMP::parameters::visitor(Tri_source_tagger<decltype(fid), Mesh>(fid)));
  if (dump_debug_obj)
  {
    write_mesh_obj(tm, debug_prefix + "01_tm_triangulated.obj", [&](Face_index f)
    {
      return "source=" + std::to_string(get(fid, f));
    });
  }

  // ④ 切割 + 标记传播
  //    注:两网格 split() 不消费 throw_on_self_intersection(clip.h:1325-1347 未提取、
  //    未转发给 corefine),故不传;visitor 的工厂/链式写法依据 clip.h:1346。
  Tag_propagator<decltype(fid), Mesh> tagger(fid, tm);   // 守卫:只传播 tm 的标记
  Mesh splitter_copy(splitter);
  if (dump_debug_obj)
  {
    write_mesh_obj(splitter_copy, debug_prefix + "02_splitter.obj");
  }
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
  if (dump_debug_obj)
  {
    write_mesh_obj(tm, debug_prefix + "03_tm_after_split.obj", [&](Face_index f)
    {
      return "source=" + std::to_string(get(fid, f)) +
             " piece=" + std::to_string(get(pid, f));
    });
  }

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

  if (dump_debug_obj)
  {
    write_pieces_obj(result, debug_prefix + "04_pieces.obj");
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
