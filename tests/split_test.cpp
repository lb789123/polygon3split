// split_test.cpp —— 验证 docs/polygon-split-by-folded-surface.md 的全部 [推断]/[未验证] 论断
// 用例 1-7 期望值全部可手算(文档第 8 节);8-10 是 ≥5 顶点多边形面三角化的回归
// (非平面 / 强非凸 / 重复连续顶点)。任一 CHECK 失败即回溯对应论断。
//
// 通过标准:全部 CHECK 通过;result.mesh 通过 is_valid_polygon_mesh;
//           component_count == pieces.size();unmarked_faces == 0。

#include "split_polygons.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/boost/graph/helpers.h>   // is_valid_polygon_mesh
#include <CGAL/number_utils.h>          // to_double

#include <array>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef poly_split::Polygon_mesh<K>     Mesh;

// ── 极简测试框架 ─────────────────────────────────────────────────────
static int g_failures = 0;
static int g_checks   = 0;
#define CHECK(cond, msg)                                       \
  do {                                                         \
    ++g_checks;                                                \
    if (!(cond)) {                                             \
      ++g_failures;                                            \
      std::cerr << "  FAIL[line " << __LINE__ << "]: " << msg << "\n"; \
    }                                                          \
  } while (0)

// ── 环比较:EPECK → double,允许旋转起点 + 反向(边界半边走向与面序相反,文档⑥)──
typedef std::array<double, 3> DPoint;

// 期望环上的 2D 点(z=0 平面用例)
static DPoint P2(double x, double y) { return { x, y, 0.0 }; }

static DPoint to_dp(const Point_3& p)
{ return { CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()) }; }

static bool close3(const DPoint& a, const DPoint& b, double tol)
{ return std::abs(a[0]-b[0])<=tol && std::abs(a[1]-b[1])<=tol && std::abs(a[2]-b[2])<=tol; }

static std::string str(const DPoint& p)
{ return "(" + std::to_string(p[0]) + ", " + std::to_string(p[1]) + ", " + std::to_string(p[2]) + ")"; }

// 共线点归一化:剥掉与前后邻点共线的顶点(容差绝对值 1e-7,用例坐标 O(1))。
// 原因:triangulate_faces 的对角线可能与切线相交(如测试 1 对角线交 x=1 于 (1,1)),
// 该交点会成为网格顶点、落在片的直线边上并出现在环里 —— 但它是几何冗余点,
// 与片多边形的角点无关。剥掉后比较,测试对三角化的对角线选择鲁棒。
static std::vector<DPoint> drop_collinear(const std::vector<DPoint>& ring)
{
  const std::size_t n = ring.size();
  std::vector<DPoint> out;
  for (std::size_t i = 0; i < n; ++i) {
    const DPoint& a = ring[(i + n - 1) % n];
    const DPoint& b = ring[i];
    const DPoint& c = ring[(i + 1) % n];
    const double cr = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]); // xy 叉积
    if (std::abs(cr) > 1e-7) out.push_back(b);
  }
  return out;
}

// act 与 exp 是否为同一环(先各自共线归一化,再循环旋转 + 允许反向)
static bool ring_match(std::vector<DPoint> exp, const std::vector<Point_3>& act, double tol = 1e-12)
{
  std::vector<DPoint> a;
  for (const auto& p : act) a.push_back(to_dp(p));
  exp = drop_collinear(exp);
  a = drop_collinear(a);
  if (exp.size() != a.size()) return false;
  const std::size_t n = exp.size();
  for (std::size_t rot = 0; rot < n; ++rot) {
    bool fwd = true, rev = true;
    for (std::size_t i = 0; i < n; ++i) {
      const DPoint fwd_p = a[(i + rot) % n];       // 正向
      const DPoint rev_p = a[(rot + n - i) % n];   // 反向
      fwd = fwd && close3(exp[i], fwd_p, tol);
      rev = rev && close3(exp[i], rev_p, tol);
    }
    if (fwd || rev) return true;
  }
  return false;
}

// 打印 + 通用校验(所有测试共用)
static void report(const char* name, poly_split::Split_result<K>& r)
{
  std::cout << "== " << name << " ==\n";
  std::cout << "  vertices=" << r.mesh.number_of_vertices()
            << " faces=" << r.mesh.number_of_faces()
            << " components=" << r.component_count
            << " pieces=" << r.pieces.size()
            << " unmarked=" << r.unmarked_faces << "\n";
  for (const auto& pc : r.pieces) {
    std::cout << "  piece: source=" << pc.source << " piece_id=" << pc.piece << " ring:";
    for (const auto& p : pc.ring) std::cout << " " << str(to_dp(p));
    std::cout << "\n";
  }
  CHECK(r.unmarked_faces == 0, "标记传播有遗漏(unmarked != 0)");
  CHECK(r.component_count == r.pieces.size(), "组数 != 边界环数(论证 5-5)");
  CHECK(CGAL::is_valid_polygon_mesh(r.mesh), "撕开后 mesh 不合法");
}

// 期望环集合匹配:每个期望环恰好匹配一个实际片(顺序无关)。
// tol:噪声用例(顶点 z≈1e-9 微离面)放宽到 1e-7。
static void check_rings(const std::vector<std::vector<DPoint>>& expected,
                        const std::vector<poly_split::PolygonPiece<K>>& pieces,
                        double tol = 1e-12)
{
  CHECK(expected.size() == pieces.size(), "片数不符");
  std::vector<bool> used(pieces.size(), false);
  for (const auto& exp : expected) {
    bool found = false;
    for (std::size_t i = 0; i < pieces.size(); ++i) {
      if (!used[i] && ring_match(exp, pieces[i].ring, tol)) { used[i] = true; found = true; break; }
    }
    CHECK(found, "未找到期望环(" + str(exp.front()) + " ...)");
  }
}

// ── 几何构造小工具 ───────────────────────────────────────────────────
static Point_3 P(double x, double y, double z) { return Point_3(x, y, z); }

// V 形折面:屋脊 C(t)=(rx,1+t,t),t∈[-3,3];两坡向 ±x 下探 u∈[0,3]。
// 两坡共享屋脊顶点(否则两分量沿屋脊棱几何重合,行为未定义)。
static Mesh make_v_splitter(double rx)
{
  const Point_3 c0(rx, -2, -3);   // C(-3)
  const Point_3 c1(rx,  4,  3);   // C(3)
  const Point_3 l2 = P(rx - 3, 4, 0);    // C(3)+3*(-1,0,-1)
  const Point_3 l3 = P(rx - 3, -2, -6);  // C(-3)+3*(-1,0,-1)
  const Point_3 r2 = P(rx + 3, -2, -6);  // C(-3)+3*(1,0,-1)
  const Point_3 r3 = P(rx + 3, 4, 0);    // C(3)+3*(1,0,-1)

  Mesh s;
  auto v0 = s.add_vertex(c0), v1 = s.add_vertex(c1);
  auto v2 = s.add_vertex(l2), v3 = s.add_vertex(l3);
  auto v4 = s.add_vertex(r2), v5 = s.add_vertex(r3);
  std::vector<Mesh::Vertex_index> left { v0, v1, v2, v3 };
  std::vector<Mesh::Vertex_index> right{ v1, v0, v4, v5 }; // 屋脊方向相反 → 两坡取向一致
  s.add_face(left);
  s.add_face(right);
  CGAL::Polygon_mesh_processing::triangulate_faces(s);
  return s;
}

// 竖直平板 x=v(y,z∈[-3,3]) —— 退化折面(平面)
static Mesh make_wall_splitter(double x)
{
  return poly_split::build_splitter<K>({
    { P(x, -3, -3), P(x, 3, -3), P(x, 3, 3), P(x, -3, 3) }
  });
}

// ── 测试 1:竖直平板切正方形(基本链路)──────────────────────────────
static void test1()
{
  auto r = poly_split::split_polygons<K>(
    { { P(0,0,0), P(2,0,0), P(2,2,0), P(0,2,0) } },
    make_wall_splitter(1));
  report("test1 plane cut", r);
  CHECK(r.pieces.size() == 2, "应得 2 片");
  // 实测(test0):4 原始 + 3 交点(含对角线交点 (1,1))+ 3 复制 = 10
  CHECK(r.mesh.number_of_vertices() == 10, "顶点应为 4+3+3=10");
  check_rings({
    { P2(0,0), P2(1,0), P2(1,2), P2(0,2) },
    { P2(1,0), P2(2,0), P2(2,2), P2(1,2) },
  }, r.pieces);
}

// ── 测试 2:V 形折面,折痕穿过多边形内部(核心用例)───────────────────
// 切割折线 A(0,1.8)→c(0.8,1)→B(1.8,2),拐点 c 在多边形内部(文档§8 已手算)。
static void test2()
{
  auto r = poly_split::split_polygons<K>(
    { { P(0,0,0), P(2,0,0), P(2,2,0), P(0,2,0) } },
    make_v_splitter(0.8));
  report("test2 V-fold cut", r);
  CHECK(r.pieces.size() == 2, "应得 2 片");
  check_rings({
    { P2(0,1.8), P2(0,2), P2(1.8,2), P2(0.8,1) },                            // 含角点 (0,2) 的 4 边形
    { P2(0,1.8), P2(0,0), P2(2,0), P2(2,2), P2(1.8,2), P2(0.8,1) },          // 其余 6 边形
  }, r.pieces);
  // 实测:4 + 交点 A/c/B 3 + 路径顶点全复制 3 = 10(内部拐点 c 也两侧各一份,
  // 与文档初版"只复制原边界端点"的推断相反 —— 已按实测修正文档 §5-3)
  CHECK(r.mesh.number_of_vertices() == 10, "顶点应为 4+3+3=10");
}

// ── 测试 3:不相交(折面远离)────────────────────────────────────────
static void test3()
{
  auto r = poly_split::split_polygons<K>(
    { { P(0,0,0), P(2,0,0), P(2,2,0), P(0,2,0) } },
    make_v_splitter(10));
  report("test3 no intersection", r);
  CHECK(r.pieces.size() == 1, "应得 1 片");
  CHECK(r.mesh.number_of_vertices() == 4, "顶点数应不变");
  check_rings({ { P2(0,0), P2(2,0), P2(2,2), P2(0,2) } }, r.pieces);
}

// ── 测试 4:批量输入(形态 2,多连通分量 + 来源标记)──────────────────
static void test4()
{
  auto r = poly_split::split_polygons<K>(
    { { P(0,0,0), P(2,0,0), P(2,2,0), P(0,2,0) },      // S1,来源 0
      { P(3,0,0), P(5,0,0), P(5,2,0), P(3,2,0) } },    // S2,来源 1(切不到,原样通过)
    make_wall_splitter(1));
  report("test4 batch input", r);
  CHECK(r.pieces.size() == 3, "应得 3 片");
  // 实测:输入 8 + S1 交点 3((1,0),(1,1) 对角线交,(1,2))+ 复制 3 = 14;S2 无变化
  CHECK(r.mesh.number_of_vertices() == 14, "顶点应为 8+3+3=14");
  check_rings({
    { P2(0,0), P2(1,0), P2(1,2), P2(0,2) },
    { P2(1,0), P2(2,0), P2(2,2), P2(1,2) },
    { P2(3,0), P2(5,0), P2(5,2), P2(3,2) },
  }, r.pieces);
  // 来源断言:S1 两片来源 0,S2 一片来源 1
  std::size_t src0 = 0, src1 = 0;
  for (const auto& pc : r.pieces) { if (pc.source == 0) ++src0; else if (pc.source == 1) ++src1; }
  CHECK(src0 == 2 && src1 == 1, "来源应为 {0,0,1}");
}

// ── 测试 5:相邻多边形(共边)+ 切线穿过共享边 ─────────────────────────
// 两个四边形:下 0≤y≤1(多边形 0)、上 1≤y≤2(多边形 1)共享边 y=1;
// 切割线 x=1 穿共享边于 (1,1,0)。每个多边形各分成左右两片 → 共 4 片,
// 来源 {0,0,1,1}(左片=下+上 经未切断的共享边左半段连通,但 barrier 按多边形隔开)。
static void test5()
{
  Mesh tm;
  std::vector<Mesh::Vertex_index> v;
  for (const auto& p : { P(0,0,0), P(2,0,0), P(2,1,0), P(0,1,0), P(2,2,0), P(0,2,0) })
    v.push_back(tm.add_vertex(p));
  // 下四边形 (0,0),(2,0),(2,1),(0,1)= 多边形 0;上四边形 (0,1),(2,1),(2,2),(0,2)
  // = 多边形 1;共享 v2-v3(同顶点 → 真共享边)
  std::vector<Mesh::Vertex_index> bottom { v[0], v[1], v[2], v[3] };
  std::vector<Mesh::Vertex_index> top    { v[3], v[2], v[4], v[5] };
  tm.add_face(bottom);
  tm.add_face(top);

  auto r = poly_split::split_mesh<K>(tm, make_wall_splitter(1));
  report("test5 adjacent polygons, cut crosses shared edge", r);
  CHECK(r.pieces.size() == 4, "应得 4 片(两个多边形 × 左右各一)");
  // 实测:输入 6 + 交点 5((1,0),(1,0.5) 下四边形对角线交,(1,1) 共享边交,
  // (1,1.5) 上四边形对角线交,(1,2))+ 复制 5 = 16 —— 共享边交点同样两侧各一份
  CHECK(r.mesh.number_of_vertices() == 16, "顶点应为 6+5+5=16");
  check_rings({
    { P2(0,0), P2(1,0), P2(1,1), P2(0,1) },   // 下左:y∈[0,1] 段上边是与上左片的共享边
    { P2(1,0), P2(2,0), P2(2,1), P2(1,1) },   // 下右
    { P2(0,1), P2(1,1), P2(1,2), P2(0,2) },   // 上左
    { P2(1,1), P2(2,1), P2(2,2), P2(1,2) },   // 上右
  }, r.pieces);
  std::size_t src0 = 0, src1 = 0;
  for (const auto& pc : r.pieces) { if (pc.source == 0) ++src0; else if (pc.source == 1) ++src1; }
  CHECK(src0 == 2 && src1 == 2, "来源应为 {0,0,1,1}");
}

// ── 测试 6:相邻多边形,切割线只切共享边一侧、两片经未切共享边连通 ──────
// 两单位方块 S1=[0,1]²(多边形 0)、S2=[1,2]×[0,1](多边形 1)共享边 x=1;
// 横墙 y=0.5 同时切两个。4 片来源 {0,0,1,1};关键:S1上 与 S2上 经未切的
// 共享边段 x=1,y∈[0.5,1] 边连通(同一撕开分量),仍按多边形分作两组 ——
// barrier 分组 + 共享边作为组边界出环的核心验证。
static void test6()
{
  Mesh tm;
  std::vector<Mesh::Vertex_index> v;
  for (const auto& p : { P(0,0,0), P(1,0,0), P(2,0,0), P(0,1,0), P(1,1,0), P(2,1,0) })
    v.push_back(tm.add_vertex(p));
  std::vector<Mesh::Vertex_index> s1 { v[0], v[1], v[4], v[3] };  // (0,0),(1,0),(1,1),(0,1)
  std::vector<Mesh::Vertex_index> s2 { v[1], v[2], v[5], v[4] };  // (1,0),(2,0),(2,1),(1,1)
  tm.add_face(s1);
  tm.add_face(s2);

  Mesh wall = poly_split::build_splitter<K>({
    { P(-3,0.5,-3), P(3,0.5,-3), P(3,0.5,3), P(-3,0.5,3) } });
  auto r = poly_split::split_mesh<K>(tm, wall);
  report("test6 adjacent squares, horizontal wall cuts both", r);
  CHECK(r.pieces.size() == 4, "应得 4 片");
  // 实测:输入 6 + 交点 5((0,0.5),(0.5,0.5) S1 对角线交,(1,0.5) 共享边交,
  // (1.5,0.5) S2 对角线交,(2,0.5))+ 复制 5 = 16
  CHECK(r.mesh.number_of_vertices() == 16, "顶点应为 6+5+5=16");
  check_rings({
    { P2(0,0), P2(1,0), P2(1,0.5), P2(0,0.5) },   // S1 下:右边 x=1,y∈[0,0.5] 是与 S2 下的共享边
    { P2(0,0.5), P2(1,0.5), P2(1,1), P2(0,1) },   // S1 上:右边 x=1,y∈[0.5,1] 是与 S2 上的共享边
    { P2(1,0), P2(2,0), P2(2,0.5), P2(1,0.5) },   // S2 下
    { P2(1,0.5), P2(2,0.5), P2(2,1), P2(1,1) },   // S2 上
  }, r.pieces);
  std::size_t src0 = 0, src1 = 0;
  for (const auto& pc : r.pieces) { if (pc.source == 0) ++src0; else if (pc.source == 1) ++src1; }
  CHECK(src0 == 2 && src1 == 2, "来源应为 {0,0,1,1}");
}

// ── 测试 7:preflight 前置体检(诊断"正常 mesh 也挂")───────────────────
// 1) 正常输入 → 无违规;2) 自交(竖直三角形穿 z=0 面);3) 退化面(共线三角形);
// 4) 非三角面。split 违反前置不抛异常(§3④),体检在进入 corefine 前给出原因。
static bool has_problem(const std::vector<std::string>& ps, const char* kw)
{
  for (const auto& s : ps) if (s.find(kw) != std::string::npos) return true;
  return false;
}

static void test7()
{
  // 1) 正常输入(测试 6 的两方块 + 横墙)
  {
    Mesh tm;
    std::vector<Mesh::Vertex_index> v;
    for (const auto& p : { P(0,0,0), P(1,0,0), P(2,0,0), P(0,1,0), P(1,1,0), P(2,1,0) })
      v.push_back(tm.add_vertex(p));
    tm.add_face(std::vector<Mesh::Vertex_index>{ v[0], v[1], v[4], v[3] });
    tm.add_face(std::vector<Mesh::Vertex_index>{ v[1], v[2], v[5], v[4] });
    Mesh wall = poly_split::build_splitter<K>({
      { P(-3,0.5,-3), P(3,0.5,-3), P(3,0.5,3), P(-3,0.5,3) } });
    auto ps = poly_split::preflight_check<K>(tm, wall);
    std::cout << "== test7 preflight ==\n  normal: " << ps.size() << " problems\n";
    for (const auto& s : ps) std::cout << "    " << s << "\n";
    CHECK(ps.empty(), "正常输入不应有违规");
  }
  // 2) 自交:两个 z=0 三角形拼成的面 + 一个 x=1 竖直三角形穿过它
  {
    Mesh tm;
    std::vector<Mesh::Vertex_index> v;
    for (const auto& p : { P(0,0,0), P(2,0,0), P(2,1,0), P(0,1,0),
                           P(1,-1,-1), P(1,2,1), P(1,0.5,2) })
      v.push_back(tm.add_vertex(p));
    tm.add_face(std::vector<Mesh::Vertex_index>{ v[0], v[1], v[2] });  // z=0 面的一半
    tm.add_face(std::vector<Mesh::Vertex_index>{ v[0], v[2], v[3] });  // 另一半(共享对角边,不算自交)
    tm.add_face(std::vector<Mesh::Vertex_index>{ v[4], v[5], v[6] });  // x=1 竖直三角形,穿过 z=0 面
    Mesh dummy = poly_split::build_splitter<K>({
      { P(-3,9,-3), P(3,9,-3), P(3,9,3), P(-3,9,3) } });  // 远离,不参与
    auto ps = poly_split::preflight_check<K>(tm, dummy);
    std::cout << "  self-intersect: " << (ps.empty() ? std::string("(none!)") : ps[0]) << "\n";
    CHECK(has_problem(ps, "自交"), "自交输入应被查出");
  }
  // 3) 退化面:三点共线的三角形
  {
    Mesh tm;
    std::vector<Mesh::Vertex_index> v;
    for (const auto& p : { P(0,0,0), P(1,0,0), P(2,0,0) })
      v.push_back(tm.add_vertex(p));
    tm.add_face(v);
    Mesh dummy = poly_split::build_splitter<K>({
      { P(-3,9,-3), P(3,9,-3), P(3,9,3), P(-3,9,3) } });
    auto ps = poly_split::preflight_check<K>(tm, dummy);
    std::cout << "  degenerate: " << (ps.empty() ? std::string("(none!)") : ps[0]) << "\n";
    CHECK(has_problem(ps, "退化"), "退化面应被查出");
  }
  // 4) 非三角 splitter:quad 面(未三角化)放在 splitter 位置才是违规
  //    (tm 位置的 quad 合法 —— 管线会先三角化,第 1 项已证)
  {
    Mesh tm;
    std::vector<Mesh::Vertex_index> v;
    for (const auto& p : { P(0,0,0), P(2,0,0), P(2,1,0), P(0,1,0) })
      v.push_back(tm.add_vertex(p));
    tm.add_face(std::vector<Mesh::Vertex_index>(v));
    Mesh quad_splitter;
    std::vector<Mesh::Vertex_index> w;
    for (const auto& p : { P(1,3,-3), P(1,3,3), P(1,-3,3), P(1,-3,-3) })
      w.push_back(quad_splitter.add_vertex(p));
    quad_splitter.add_face(std::vector<Mesh::Vertex_index>(w));
    auto ps = poly_split::preflight_check<K>(tm, quad_splitter);
    std::cout << "  non-triangle splitter: " << (ps.empty() ? std::string("(none!)") : ps[0]) << "\n";
    CHECK(has_problem(ps, "非三角"), "非三角 splitter 应被查出");
  }
}

// ── 测试 8:非平面 ≥5 边形(L 形六边形,凹点 z=1e-9)─────────────────────
// 回归:PMP::triangulate_faces 的 ≥5 顶点洞填充路径对此类输入会静默回退
// 3D Delaunay 洞填充(triangulate_hole.h:769-771 的共面容差 = bbox z 跨度²/16,
// 近水平多边形容差≈0)→ 三角形乱穿。管线②现改走 triangulate_polygon_faces
// 的平面化 CDT。切割几何与严格平面版完全一致:竖墙 x=1 切 L((0,0),(2,0),
// (2,1),(0.7,1),(0.7,2),(0,2)) 成右四边形 + 左六边形;噪声顶点出现在左片环里
// (z≈1e-9,环比较容差放宽到 1e-7)。
static void test8()
{
  const std::vector<std::vector<DPoint>> expect = {
    { P2(1,0), P2(2,0), P2(2,1), P2(1,1) },                                  // 右片
    { P2(1,0), P2(1,1), P2(0.7,1), P2(0.7,2), P2(0,2), P2(0,0) },            // 左片
  };
  {
    auto r = poly_split::split_polygons<K>(
      { { P(0,0,0), P(2,0,0), P(2,1,0), P(0.7,1,0), P(0.7,2,0), P(0,2,0) } },
      make_wall_splitter(1));
    report("test8 planar L-hexagon (>=5 vertices)", r);
    CHECK(r.pieces.size() == 2, "应得 2 片");
    check_rings(expect, r.pieces);
  }
  {
    auto r = poly_split::split_polygons<K>(
      { { P(0,0,0), P(2,0,0), P(2,1,0), P(0.7,1,1e-9), P(0.7,2,0), P(0,2,0) } },
      make_wall_splitter(1));
    report("test8 non-planar L-hexagon (reflex vertex z=1e-9)", r);
    CHECK(r.pieces.size() == 2, "应得 2 片(与平面版一致)");
    check_rings(expect, r.pieces, /*tol=*/1e-7);   // 噪声顶点 z≈1e-9 放宽
  }
}

// ── 测试 9:强非凸 C 形八边形(带噪声)切三片 ───────────────────────────
// 竖墙 x=1 过 C 形缺口:上下两横条各被切,右柱连通 → 左下/右侧连通/左上三片。
// 8 顶点 ≥5,凹角点 (2.5,2.5) 抬 z=1e-9 → 平面化 CDT 路径 + 非平面回归。
static void test9()
{
  auto r = poly_split::split_polygons<K>(
    { { P(0,0,0), P(3,0,0), P(3,3,0), P(0,3,0), P(0,2.5,0),
        P(2.5,2.5,1e-9), P(2.5,0.5,0), P(0,0.5,0) } },
    make_wall_splitter(1));
  report("test9 C-shape octagon with noise", r);
  CHECK(r.pieces.size() == 3, "应得 3 片(左下/右侧连通/左上)");
  check_rings({
    { P2(0,0), P2(1,0), P2(1,0.5), P2(0,0.5) },                              // 左下
    { P2(1,0), P2(3,0), P2(3,3), P2(1,3), P2(1,2.5),                         // 右侧连通
      P2(2.5,2.5), P2(2.5,0.5), P2(1,0.5) },
    { P2(0,2.5), P2(1,2.5), P2(1,3), P2(0,3) },                              // 左上
  }, r.pieces, /*tol=*/1e-7);
  std::size_t src0 = 0;
  for (const auto& pc : r.pieces) if (pc.source == 0) ++src0;
  CHECK(src0 == 3, "三片来源都应为 0");
}

// ── 测试 10:环上重复连续顶点(PMP 路径会静默不三角化)───────────────────
// 与测试 8 同一 L,但 (0.7,1) 重复两次:管线去重后与平面版完全一致。
static void test10()
{
  auto r = poly_split::split_polygons<K>(
    { { P(0,0,0), P(2,0,0), P(2,1,0), P(0.7,1,0), P(0.7,1,0), P(0.7,2,0), P(0,2,0) } },
    make_wall_splitter(1));
  report("test10 ring with duplicated consecutive vertex", r);
  CHECK(r.pieces.size() == 2, "应得 2 片(与测试 8 平面版一致)");
  check_rings({
    { P2(1,0), P2(2,0), P2(2,1), P2(1,1) },
    { P2(1,0), P2(1,1), P2(0.7,1), P2(0.7,2), P2(0,2), P2(0,0) },
  }, r.pieces);
}

// ── 测试 0(临时):裸 split,不带 visitor/属性映射 —— 二分定位崩溃来源 ──
static void test0_raw()
{
  std::cerr << "== test0 raw split ==\n";
  Mesh tm;
  std::vector<Mesh::Vertex_index> v;
  for (const auto& p : { P(0,0,0), P(2,0,0), P(2,2,0), P(0,2,0) })
    v.push_back(tm.add_vertex(p));
  tm.add_face(v);
  Mesh s = make_wall_splitter(1);
  std::cerr << "  triangulate tm...\n";
  CGAL::Polygon_mesh_processing::triangulate_faces(tm);
  std::cerr << "  raw PMP::split...\n";
  CGAL::Polygon_mesh_processing::split(tm, s);
  std::cerr << "  raw split done: V=" << tm.number_of_vertices()
            << " F=" << tm.number_of_faces() << "\n";
}

int main()
{
  std::cout << std::unitbuf;   // 崩溃定位:每次输出立即刷新
  std::cout << "CGAL " << CGAL_VERSION_STR << " | kernel: EPECK\n";
  test0_raw(); std::cout << "\n";
  test1(); std::cout << "\n";
  test2(); std::cout << "\n";
  test3(); std::cout << "\n";
  test4(); std::cout << "\n";
  test5(); std::cout << "\n";
  test6(); std::cout << "\n";
  test7(); std::cout << "\n";
  test8(); std::cout << "\n";
  test9(); std::cout << "\n";
  test10(); std::cout << "\n";

  std::cout << g_checks << " checks, " << g_failures << " failures\n";
  return g_failures == 0 ? 0 : 1;
}
