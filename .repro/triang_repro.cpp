// triang_repro.cpp —— 验证管线修复后(≥5 顶点面走 triangulate_polygon_faces)
// 的三角化质量;同时用 Hole_filling visitor 钩子确认 ≥5 面不再进 PMP 洞填充。
//
// 旧结论(修复前,裸 PMP::triangulate_faces):≥5 顶点面走 triangulate_hole_polyline
//      (use_2d_constrained_delaunay_triangulation=true):
//        ① 先试平面 CDT2(任一前置失败即静默放弃:法线/共面性/简单性/顶点数/constraint 相交)
//        ② 回退 use_delaunay_triangulation(默认 true)的 3D Delaunay 洞填充
//      共面闸门容差 = bbox vertex(0)-vertex(5) 平方距离/16 —— CGAL 6.2 顶点序里
//      这两点仅差 z ⇒ 近水平多边形容差≈0,1e-9 级不共面即被拒、静默回退 DT3。
//      修复:src/split_polygons.h 的 triangulate_polygon_faces(Newell 法向投影 +
//      精确 2D CDT,无共面闸门;连续重复点先去重),quad 仍交 PMP。
//
// 校验(精确 EPECK):
//        - 管线路径成功(triangulate_polygon_faces 返回组数 == 原面数)
//        - 网格结构合法 + 全三角 + 三角形数 == n-2
//        - 每个三角形法向与多边形矢量面积同向(朝向一致,无翻面)
//        - 三角形矢量面积之和 == 多边形矢量面积(精确相等 ⇔ 恰好铺满不多不少)
//        - 无自交
#include "split_polygons.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/number_utils.h>

#include <filesystem>
#include <fstream>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using K      = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;
using Mesh   = CGAL::Surface_mesh<Point_3>;
namespace PMP = CGAL::Polygon_mesh_processing;

// ── 探测 visitor:继承项目同款 Tri_source_tagger 的基类,附加 planar_phase 钩子 ──
struct Probe_state { int planar_ok = 0, planar_fail = 0; };

struct Probe_visitor : PMP::Triangulate_faces::Default_visitor<Mesh>
{
  using face_descriptor = Mesh::Face_index;
  Probe_state* st = nullptr;
  Probe_visitor() = default;
  explicit Probe_visitor(Probe_state* s) : st(s) {}

  void end_planar_phase(bool ok) const
  {
    if (st) { if (ok) ++st->planar_ok; else ++st->planar_fail; }
  }
  // 保留来源标记行为(与 split_polygons.h ② 相同的调用形态)
  Mesh::Property_map<Mesh::Face_index, int>* fid = nullptr;
  int cur = -1;
  void before_subface_creations(face_descriptor f_old)
  { if (fid) cur = (*fid)[f_old]; }
  void after_subface_created(face_descriptor f_new)
  { if (fid) (*fid)[f_new] = cur; }
};

static std::string vf(const K::FT& x) { return std::to_string(CGAL::to_double(x)); }

static void dump_obj(const Mesh& m, const std::string& path)
{
  std::ofstream out(path, std::ios::binary);
  if (!out) return;
  auto vpm = get(CGAL::vertex_point, m);
  auto line = [&out](const std::string& s) { out << s << "\r\n"; };
  for (auto v : vertices(m))
  {
    const auto& p = get(vpm, v);
    line("v " + vf(p.x()) + " " + vf(p.y()) + " " + vf(p.z()));
  }
  for (auto f : faces(m))
  {
    std::string s = "f";
    for (auto h : halfedges_around_face(halfedge(f, m), m))
      s += " " + std::to_string(target(h, m).idx() + 1);
    line(s);
  }
}

struct CaseResult
{
  std::string name;
  std::size_t n_in = 0, n_out = 0;      // 输入环点数 / 输出三角形数
  bool pipe_ok = false;                 // 管线路径成功(≥5 走平面化 CDT,无回退)
  bool structure_ok = false, all_tri = false;
  bool orientation_ok = false, area_exact = false, no_self_intersect = false;
  int planar_ok = 0, planar_fail = 0;   // PMP 洞填充路径计数(≥5 修复后应不再触发)
  bool expect_degenerate = false;        // 已知退化输入:只观测,不判失败
  std::string note;
};

static CaseResult run_case(const std::string& name,
                           const std::vector<Point_3>& ring,
                           const std::filesystem::path& out_dir,
                           bool expect_degenerate = false,
                           const std::string& note = std::string())
{
  CaseResult r;
  r.name = name;
  r.n_in = ring.size();
  r.expect_degenerate = expect_degenerate;
  r.note = note;

  Mesh m;
  std::vector<Mesh::Vertex_index> vs;
  for (const auto& p : ring) vs.push_back(m.add_vertex(p));
  m.add_face(vs);
  dump_obj(m, (out_dir / (name + "_in.obj")).string());

  // ── 管线路径(src/split_polygons.h ②):≥5 走平面化 CDT,quad 交 PMP ──
  const std::size_t n_faces0 = m.number_of_faces();
  auto groups = poly_split::triangulate_polygon_faces<K>(m);
  r.pipe_ok = (groups.size() == n_faces0);

  auto fid = m.add_property_map<Mesh::Face_index, int>("f:source", -1).first;
  for (std::size_t k = 0; k < groups.size(); ++k)
    for (auto f : groups[k]) put(fid, f, static_cast<int>(k));

  Probe_state st;
  Probe_visitor vis(&st);
  vis.fid = &fid;
  PMP::triangulate_faces(m, PMP::parameters::visitor(vis));   // quad 收尾(≥5 已是三角)
  r.planar_ok = st.planar_ok;
  r.planar_fail = st.planar_fail;

  dump_obj(m, (out_dir / (name + "_out.obj")).string());

  // 结构
  r.structure_ok = m.is_valid() && CGAL::is_valid_polygon_mesh(m, false);
  r.all_tri = true;
  r.n_out = 0;
  for (auto f : faces(m))
  {
    std::size_t d = 0;
    for (auto h : halfedges_around_face(halfedge(f, m), m)) { if (++d > 3) break; }
    if (d != 3) { r.all_tri = false; break; }
    ++r.n_out;
  }

  // 几何:矢量面积(×2 面积,即 Σ p_i × p_{i+1})
  Vector_3 Avec(0, 0, 0);   // 环的 2×矢量面积
  for (std::size_t i = 0; i < ring.size(); ++i)
  {
    const auto& a = ring[i];
    const auto& b = ring[(i + 1) % ring.size()];
    Avec = Avec + CGAL::cross_product(Vector_3(a.x(), a.y(), a.z()),
                                      Vector_3(b.x(), b.y(), b.z()));
  }
  Vector_3 Tsum(0, 0, 0);   // 三角形 2×矢量面积之和
  r.orientation_ok = true;
  auto vpm = get(CGAL::vertex_point, m);
  for (auto f : faces(m))
  {
    std::vector<Point_3> pts;
    for (auto h : halfedges_around_face(halfedge(f, m), m)) pts.push_back(get(vpm, target(h, m)));
    if (pts.size() != 3) { r.orientation_ok = false; continue; }
    const Vector_3 n = CGAL::cross_product(pts[1] - pts[0], pts[2] - pts[0]);
    if (n * Avec <= 0) r.orientation_ok = false;      // 翻面/反向
    Tsum = Tsum + n;
  }
  r.area_exact = (Tsum == Avec);
  r.no_self_intersect = !PMP::does_self_intersect(m);
  return r;
}

static void report(const CaseResult& r)
{
  const bool core = r.pipe_ok && r.structure_ok && r.all_tri && r.orientation_ok
                    && r.area_exact && r.no_self_intersect;
  const bool pass = core || r.expect_degenerate;
  std::cout << (pass ? "[PASS*]" : "[FAIL ]") << " " << r.name
            << "  n=" << r.n_in << " -> tris=" << r.n_out
            << "  pipe=" << r.pipe_ok
            << "  pmp_holefill=" << (r.planar_ok + r.planar_fail)
            << "  struct=" << r.structure_ok << " alltri=" << r.all_tri
            << " orient=" << r.orientation_ok << " area=" << r.area_exact
            << " no-self-int=" << r.no_self_intersect
            << (r.expect_degenerate ? "  (退化输入,仅观测)" : "")
            << (r.note.empty() ? "" : "  -- " + r.note)
            << "\n";
}

int main()
{
  const auto out_dir = std::filesystem::current_path().parent_path() / "out";
  std::filesystem::create_directories(out_dir);

  std::cout << "=== PMP::triangulate_faces 对 ≥5 顶点面的行为复现 ===\n\n";

  { // A:凸五边形,z=0
    std::vector<Point_3> ring{{0,0,0},{2,0,0},{3,1,0},{1.5,2,0},{0,1,0}};
    report(run_case("A_convex5_z0", ring, out_dir));
  }
  { // B:凸五边形,斜平面 x+y+z=1
    // 平面上取 2D 坐标 (u,v),p = (u, v, 1-u-v)
    auto onp = [](double u, double v) { return Point_3(u, v, 1 - u - v); };
    std::vector<Point_3> ring{onp(0,0), onp(2,0), onp(3,1), onp(1.5,2), onp(0,1)};
    report(run_case("B_convex5_slant", ring, out_dir));
  }
  { // C:L 形六边形(非凸),z=0
    std::vector<Point_3> ring{{0,0,0},{2,0,0},{2,1,0},{1,1,0},{1,2,0},{0,2,0}};
    report(run_case("C_L6_z0", ring, out_dir));
  }
  { // D:L 形六边形,竖直平面 x=y(p=(c,c,h))
    auto onv = [](double c, double h) { return Point_3(c, c, h); };
    std::vector<Point_3> ring{onv(0,0), onv(2,0), onv(2,1), onv(1,1), onv(1,2), onv(0,2)};
    report(run_case("D_L6_vertical", ring, out_dir));
  }
  { // E:L 形 + 共线冗余中点(10 顶点,三段共线)——真实数据常见
    std::vector<Point_3> ring{{0,0,0},{1,0,0},{2,0,0},{2,0.5,0},{2,1,0},
                              {1,1,0},{1,1.5,0},{1,2,0},{0.5,2,0},{0,2,0}};
    report(run_case("E_L10_collinear", ring, out_dir));
  }
  { // F:微非平面(一个点 z 抬 1e-9)——模拟 double 真实数据
    std::vector<Point_3> ring{{0,0,0},{2,0,0},{3,1,0},{1.5,2,1e-9},{0,1,0}};
    report(run_case("F_nonplanar_1e9", ring, out_dir));
  }
  { // F2:double 转换造成的非平面(0.1+0.2 != 0.3)
    const double za = 0.1 + 0.2, zb = 0.3;
    std::vector<Point_3> ring{{0,0,zb},{2,0,zb},{3,1,zb},{1.5,2,za},{0,1,zb}};
    report(run_case("F2_dbl_noise", ring, out_dir));
  }
  { // G:重复连续顶点(p4 == p5)
    std::vector<Point_3> ring{{0,0,0},{2,0,0},{2,1,0},{1,1,0},{1,1,0},{0,2,0}};
    report(run_case("G_dup_vertex", ring, out_dir));
  }
  { // H:C 形强非凸八边形,fan 三角化必然出错的那种
    std::vector<Point_3> ring{{0,0,0},{3,0,0},{3,3,0},{0,3,0},{0,2.5,0},
                              {2.5,2.5,0},{2.5,0.5,0},{0,0.5,0}};
    report(run_case("H_C8_z0", ring, out_dir));
  }
  // ── 不共面程度扫描:L 形六边形,顶点 4(凹角点 (1,1)) 抬高 delta ──
  for (double d : {1e-9, 1e-6, 1e-3, 1e-2, 0.1, 0.3, 1.0})
  {
    char nm[64];
    std::snprintf(nm, sizeof(nm), "N_lift_%g", d);
    std::vector<Point_3> ring{{0,0,0},{2,0,0},{2,1,0},{1,1,d},{1,2,0},{0,2,0}};
    report(run_case(nm, ring, out_dir));
  }
  // ── 整体弯折:V 形六边形(沿环中段整体弯起,不共面量大)──
  {
    std::vector<Point_3> ring{{0,0,0},{2,0,0},{2,1,0.4},{1,1,0.8},{1,2,0.4},{0,2,0}};
    report(run_case("N_fold_V", ring, out_dir));
  }
  // ── 真实数据规模:60 边形圆环,z 交替 ±1e-9(坐标噪声级不共面)──
  {
    std::vector<Point_3> ring;
    for (int i = 0; i < 60; ++i)
    {
      const double a = 2 * 3.14159265358979323846 * i / 60;
      const double z = (i % 2 == 0) ? 1e-9 : -1e-9;
      ring.push_back(Point_3(10 * std::cos(a), 10 * std::sin(a), z));
    }
    report(run_case("N_circle60_noise", ring, out_dir));
  }
  // ── 真实数据规模:60 边形圆环,轻微碗状弯折(z = 0.01*sin(2a))──
  {
    std::vector<Point_3> ring;
    for (int i = 0; i < 60; ++i)
    {
      const double a = 2 * 3.14159265358979323846 * i / 60;
      const double z = 0.01 * std::sin(2 * a);
      ring.push_back(Point_3(10 * std::cos(a), 10 * std::sin(a), z));
    }
    report(run_case("N_circle60_bowl", ring, out_dir));
  }
  std::cout << "\nOBJ 已导出: " << out_dir.string() << "\n";
  return 0;
}
