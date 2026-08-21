// gate_probe.cpp —— 逐闸门复刻 triangulate_hole_polyline_with_cdt 的前置检查,
// 找出到底是哪道门拒绝了 1e-9 不共面的五边形。
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/internal/Hole_filling/Triangulate_hole_polyline.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/bounding_box.h>

#include <iostream>
#include <vector>

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;

int main()
{
  // F 用例:五边形,第 4 点 z=1e-9
  std::vector<Point_3> ring{{0,0,0},{2,0,0},{3,1,0},{1.5,2,1e-9},{0,1,0}};

  std::vector<Point_3> P(ring);
  P.push_back(P.front());
  const std::size_t size = P.size() - 1;

  // ── 门 0:平均法向(逐字复刻)──
  Vector_3 avg_normal;
  {
    K::FT x = 0, y = 0, z = 0;
    int num_normals = 0;
    const Point_3& ref = P[0];
    for (std::size_t i = 1; i + 1 < size; ++i)
    {
      if (CGAL::collinear(ref, P[i], P[i+1])) continue;
      Vector_3 n = CGAL::normal(ref, P[i], P[i+1]);
      if (n.x() > 0 || (n.x() == 0 && n.y() > 0) ||
          (n.x() == 0 && n.y() == 0 && n.z() > 0))
        { x += n.x(); y += n.y(); z += n.z(); }
      else
        { x -= n.x(); y -= n.y(); z -= n.z(); }
      ++num_normals;
    }
    std::cout << "num_normals = " << num_normals << "\n";
    if (num_normals < 1) { std::cout << "门0失败: 无法计算法向\n"; return 1; }
    avg_normal = Vector_3(x / num_normals, y / num_normals, z / num_normals);
    std::cout << "avg_normal ≈ (" << CGAL::to_double(avg_normal.x())
              << ", " << CGAL::to_double(avg_normal.y())
              << ", " << CGAL::to_double(avg_normal.z()) << ")\n";
    if (avg_normal == CGAL::NULL_VECTOR) { std::cout << "门0失败: 零法向\n"; return 1; }
  }

  // ── 门 1:is_planar_2(逐字复刻调用参数)──
  {
    const K::Iso_cuboid_3 bbox = CGAL::bounding_box(P.begin(), P.end());
    for (int i = 0; i < 8; ++i)
    {
      const auto& p = bbox.vertex(i);
      std::cout << "  bbox.v" << i << " = (" << CGAL::to_double(p.x())
                << ", " << CGAL::to_double(p.y())
                << ", " << CGAL::to_double(p.z()) << ")\n";
    }
    for (int j = 1; j < 8; ++j)
      std::cout << "  sqdist(v0,v" << j << ") = "
                << CGAL::to_double(CGAL::squared_distance(bbox.vertex(0), bbox.vertex(j))) << "\n";
    K::FT max_sq = CGAL::abs(CGAL::squared_distance(bbox.vertex(0), bbox.vertex(5)));
    max_sq /= K::FT(16);
    const bool planar = CGAL::internal::is_planar_2(P, avg_normal, max_sq, K());
    std::cout << "门1 is_planar_2 (max_sq_dist=" << CGAL::to_double(max_sq)
              << ") -> " << (planar ? "通过" : "拒绝") << "\n";
    if (!planar) return 2;
  }

  // ── 门 2:投影下 is_simple_2 ──
  typedef CGAL::Projection_traits_3<K> P_traits;
  {
    const P_traits pt(avg_normal);
    const bool simple = CGAL::is_simple_2(P.begin(), P.end() - 1, pt);
    std::cout << "门2 is_simple_2 -> " << (simple ? "通过" : "拒绝") << "\n";
    if (!simple) return 3;
  }

  // ── 门 3/4:CDT 构建 + 顶点数/维度 ──
  {
    typedef CGAL::Triangulation_vertex_base_with_info_2<std::size_t, P_traits> Vb;
    typedef CGAL::Triangulation_face_base_with_info_2<bool, P_traits>          Fbi;
    typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fbi>         Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                       TDS;
    typedef CGAL::No_constraint_intersection_tag                               Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS, Itag>    CDT;
    P_traits ct(avg_normal);
    CDT cdt(ct);
    std::vector<std::pair<Point_3, std::size_t>> pts;
    for (std::size_t i = 0; i < size; ++i) pts.push_back({P[i], i});
    std::vector<CDT::Vertex_handle> vertices(size);
    cdt.insert(pts.begin(), pts.end());
    for (CDT::Vertex_handle v : cdt.finite_vertex_handles())
      vertices[v->info()] = v;
    std::cout << "门3 vertices.size()=" << vertices.size()
              << " cdt.nv=" << cdt.number_of_vertices() << " -> "
              << (vertices.size() == cdt.number_of_vertices() ? "通过" : "拒绝") << "\n";
    if (vertices.size() != cdt.number_of_vertices()) return 4;
    try
    {
      for (std::size_t i = 0; i < size; ++i)
      {
        std::size_t ip = (i + 1) % size;
        if (vertices[i] != vertices[ip])
          cdt.insert_constraint(vertices[i], vertices[ip]);
      }
      std::cout << "门4a constraint 插入 -> 通过\n";
    }
    catch (const CDT::Intersection_of_constraints_exception&)
    {
      std::cout << "门4a constraint 插入 -> 拒绝(相交异常)\n";
      return 5;
    }
    std::cout << "门4b dimension=" << cdt.dimension()
              << " nv=" << cdt.number_of_vertices()
              << " (期望 2 / " << size << ") -> "
              << ((cdt.dimension() == 2 && cdt.number_of_vertices() == size) ? "通过" : "拒绝") << "\n";

    // ── 三角形逐个过 is_valid(共线性)──
    int tris = 0, bad = 0;
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
    {
      // 只统计域内面照搬内部流程过于复杂,这里先看全量;真正的门 5 是逐三角形共线检查
      int i0 = (int)fit->vertex(0)->info(), i1 = (int)fit->vertex(1)->info(), i2 = (int)fit->vertex(2)->info();
      ++tris;
      if (CGAL::collinear(P[i0], P[i1], P[i2])) ++bad;
    }
    std::cout << "CDT 有限面总数=" << tris << " 其中 3D 共线=" << bad << "\n";
  }
  std::cout << "\n结论:以上全部通过的话,拒绝点只可能在内部流程未复刻的细节里\n";
  return 0;
}
