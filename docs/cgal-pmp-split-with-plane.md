# CGAL `PMP::split(pm, plane)` 平面切割网格详解

> 源码版本:CGAL 6.2(GPL-3.0-or-later),源码树 `D:\github\CGAL-6.2\CGAL-6.2`
> 涉及文件:
> - [`clip.h`](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h) —— `split()` 入口与 `split_along_edges()`
> - [`refine_with_plane.h`](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h) —— `refine_with_plane()` 细化实现
>
> 关联文档:[polygon-split-by-folded-surface.md](polygon-split-by-folded-surface.md) —— 在本文基础上把"平面"推广为"折面"
> (两网格版 split)的完整方案,已实现([src/split_polygons.h](../src/split_polygons.h))并通过五用例测试(2026-08-20)。

---

## 目录

1. [功能概述](#1-功能概述)
2. [函数签名与模板参数](#2-函数签名与模板参数)
3. [命名参数](#3-命名参数)
4. [整体流程](#4-整体流程)
5. [阶段一:refine_with_plane —— 求交与细化](#5-阶段一refine_with_plane--求交与细化)
6. [阶段二:split_along_edges —— 复制交线边、形成边界](#6-阶段二split_along_edges--复制交线边形成边界)
7. [典型用法](#7-典型用法)
8. [注意事项与前置条件](#8-注意事项与前置条件)
9. [源码位置索引](#9-源码位置索引)

---

## 1. 功能概述

`CGAL::Polygon_mesh_processing::split(pm, plane)` 用一个平面把输入网格 `pm` 切成**两部分,且两部分都保留**。

官方文档原文(摘自 [`clip.h:1353-1423`](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h#L1353-L1423)):

> The polygon mesh is refined with the intersection edges, and those edges are duplicated as to create a boundary, and thus separate connected components on either side of the plane.
>
> (网格被交线边所细化,随后这些边被复制以产生边界,从而把平面两侧的几何分离成各自的连通分量。)

这句话概括了算法的两个阶段:

| 阶段 | 操作 | 结果 |
|---|---|---|
| ① 细化(refine) | 把被平面穿过的面剖分,使交线成为网格中真实存在的边 | 交线边被打上标记 |
| ② 撕开(split) | 复制交线边、沿切割线断开共享顶点 | 两侧各产生一条边界,成为独立连通分量 |

若输入是封闭水密(closed & watertight)网格,输出即为两个各自封闭的网格,切口处是位置重合但拓扑独立的两组顶点/边。

### 与 `clip()` 的区别

同文件中的 `clip()` 与 `split()` 共用同一套底层机制(`refine_with_plane` + `split_along_edges`),区别在于:

- `clip()`:只保留平面一侧,丢弃另一侧(可选 `clip_volume` 补盖);
- `split()`:两侧都保留,不删除任何几何。

文档末尾的 `@see clip()` 即提示二者为同族函数。

---

## 2. 函数签名与模板参数

```cpp
template <class PolygonMesh,
          class NamedParameters = parameters::Default_named_parameters>
void split(PolygonMesh& pm,
#ifdef DOXYGEN_RUNNING
          const Plane_3& plane,      // 仅供文档显示
#else
          const typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Plane_3& plane,
#endif
          const NamedParameters& np = parameters::default_values())
```

([clip.h:1418-1423](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h#L1418-L1423))

- **`PolygonMesh`**:须为 `MutableFaceGraph`(可增删点/边/面)、`HalfedgeListGraph`、`FaceListGraph` 概念的模型,且自带 `CGAL::vertex_point_t` 属性映射。典型如 `CGAL::Surface_mesh<Kernel::Point_3>`、`CGAL::Polyhedron_3`。
- **`plane`**:平面类型由 `GetGeomTraits<...>::type::Plane_3` 推导,即**必须与网格顶点使用同一个 CGAL kernel**(例如同为 `Exact_predicates_inexact_constructions_kernel`)。
- **`#ifdef DOXYGEN_RUNNING`**:CGAL 文档惯用技巧。真实签名的平面类型是一长串嵌套 typename,文档中显示为简洁的 `Plane_3`,二者编译时等价。
- **`np`**:BGL 风格命名参数,全部可选,见下节。

---

## 3. 命名参数

| 参数 | 类型 | 默认值 | 说明 |
|---|---|---|---|
| `vertex_point_map` | `ReadWritePropertyMap`,key 为 `vertex_descriptor`,value 为 `Point_3` | `get(CGAL::vertex_point, pm)` | 顶点 → 三维坐标映射 |
| `visitor` | `PMPCorefinementVisitor` 模型 | `Corefinement::Default_visitor<PolygonMesh>` | 跟踪切割中新建的面/边/顶点。由于平面不是网格,visitor 回调中"平面一侧"的半边/面以 `null_halfedge()` / `null_face()` 填充,网格参数传 `pm` 本身 |
| `throw_on_self_intersection` | `bool` | `false` | 为 `true` 时检查交线附近的三角形组是否自交,发现则抛出 `Corefinement::Self_intersection_exception`。仅当 `pm` 是三角网格时生效 |
| `do_not_triangulate_faces` | `bool` | `false` | 输入为三角网格且此参数为 `false` 时,输出保持三角化;输入不是三角网格时该参数恒按 `true` 处理(不强制三角化) |
| `geom_traits` | `Kernel` 模型 | 由点类型经 `CGAL::Kernel_traits` 推导 | 几何 traits,必须与顶点类型兼容 |

---

## 4. 整体流程

入口实现位于 [clip.h:1424-1467](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h#L1424-L1467):

```cpp
void split(PolygonMesh& pm, const Plane_3& plane, const NamedParameters& np)
{
  // 1. 从命名参数取出 traits、点映射、并发标签
  GT traits    = ...;
  auto vpm     = ...;                     // vertex point map
  const bool throw_on_self_intersection = ...;
  bool triangulate = !do_not_triangulate_faces;

  // 2. 两个动态属性图(算法的"工作草稿")
  auto vos = get(dynamic_vertex_property_t<Oriented_side>(), pm); // 顶点在平面哪一侧
  auto ecm = get(dynamic_edge_property_t<bool>(), pm, false);     // 边是否为交线边

  // 3. 想保持三角化但输入不是纯三角网格 → 放弃三角化
  if (triangulate && !is_triangle_mesh(pm)) triangulate = false;

  // 4. 阶段一:求交 + 细化(填充 vos / ecm)
  refine_with_plane(pm, plane,
      parameters::vertex_oriented_side_map(vos)
                 .edge_is_marked_map(ecm)
                 .vertex_point_map(vpm)
                 .geom_traits(traits)
                 .do_not_triangulate_faces(!triangulate)
                 .throw_on_self_intersection(throw_on_self_intersection)
                 .concurrency_tag(Concurrency_tag())
                 .visitor(std::ref(visitor)));

  // 5. 阶段二:沿标记边把网格"撕开"
  internal::split_along_edges(pm, ecm, vpm, visitor);
}
```

流程示意:

```
        输入网格(封闭)                 阶段一 refine_with_plane
        ┌─────────────┐    求交、插点、剖分被穿过的面     ┌─────────────┐
        │   ♦♦♦♦♦     │  ────────────────────────────▶  │  ♦═══♦═══♦  │ ← 交线边已存在
        │  ♦♦   ♦♦    │         (交线成为真实边)          │  ║  ⋮  ║  │   且在 ecm 中标记
        │   ♦♦♦♦♦     │                                  │  ♦═══♦═══♦  │
        └─────────────┘                                  └──────┬──────┘
                                                                │ 阶段二 split_along_edges
                                                                ▼ 复制交线边、断开共享顶点
                                                    ┌─────────────┐
                                                    │ 上半部分(带切口边界)│
                                                    └─────────────┘
                                                    ┌─────────────┐
                                                    │ 下半部分(带切口边界)│
                                                    └─────────────┘
                                          两个连通分量,仍共存于同一个 mesh 中
```

---

## 5. 阶段一:`refine_with_plane` —— 求交与细化

实现位于 [refine_with_plane.h:241-621](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L241-L621)。它完成"把交线插入网格"的全部工作,分 5 步。

### 5.0 凸网格快速通道(split 不走此路径)

函数开头([refine_with_plane.h:295-303](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L295-L303)):若 `use_convex_specialization=true` 且调用方**未**传入 `edge_is_constrained_map` / `edge_is_marked_map` / `vertex_oriented_side_map`,则改走针对凸网格的专用快速实现。`split()` 因为传入了 `edge_is_marked_map`,**始终走通用路径**。

### 5.1 顶点分类([refine_with_plane.h:356-374](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L356-L374))

遍历所有顶点,用 `traits.oriented_side_3_object()(plane, p)` 判定每个顶点相对平面的位置,写入 `vertex_os` 映射:

- `ON_POSITIVE_SIDE` —— 平面正侧(法向一侧);
- `ON_NEGATIVE_SIDE` —— 平面负侧;
- `ON_ORIENTED_BOUNDARY` —— 恰好落在平面上(记入 `on_obnd` 列表)。

同时维护三个全局标志:`all_in`(全在负侧)、`all_out`(全在正侧)、`at_least_one_on`(至少一点在平面上)。**若网格整体在平面一侧(`all_in || all_out`),函数直接返回,不产生任何交线边**([refine_with_plane.h:422-426](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L422-L426))——此时 `split()` 后网格保持原样,仍是一个连通分量。

### 5.2 边分类([refine_with_plane.h:376-420](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L376-L420))

遍历所有边,分两类处理:

1. **穿越边(crossing edge)**:两端点分别在平面两侧(一正一负)→ 收入 `inters` 列表,稍后要在边上插交点;
2. **共面边(coplanar edge)**:两端点都恰好落在平面上 → 检查其两侧面是否**都**完全位于平面内("pure coplanar"):
   - 若是,该边属于恰好嵌在平面里的面片,不作为切割标记;
   - 若否(至少一个邻面伸出平面),该边直接标记为交线边(`put(edge_is_marked, e, true)`)。

### 5.3 可选:自交检查([refine_with_plane.h:430-447](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L430-L447))

若 `throw_on_self_intersection=true`(且是三角网格):收集所有穿越边的邻接面,去重后对其运行 `does_self_intersect()`。发现自交则抛出 `Corefinement::Self_intersection_exception`。

> 这就是文档中该参数描述的含义:只检查"交线附近的那组三角形",而非全网格,开销可控。

### 5.4 穿越边插交点([refine_with_plane.h:451-474](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L451-L474))

对每条穿越边 `e`:

```cpp
Point_3 ip = intersection_point(plane, p_src, p_tgt);  // 平面 ∩ 边段 的精确交点
h = CGAL::Euler::split_edge(h, pm);                    // 边中插入新顶点,一分为二
put(vpm, target(h, pm), ip);                           // 新顶点坐标 = 交点
put(vertex_os, target(h, pm), ON_ORIENTED_BOUNDARY);   // 新顶点归类为"在平面上"
```

同时把被该边穿过的(最多两个)面连同相应半边登记进 `splitted_faces` 映射,供下一步剖分。若原边在 `edge_is_constrained_map` 中被约束,分裂后的两条新边继承约束标记。

### 5.5 面剖分([refine_with_plane.h:479-619](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L479-L619))

先把"有顶点恰好落在平面上"的邻接面也登记进 `splitted_faces`(切线/tangency 情形:切割线只擦过顶点。若该点不是穿越点,会把同一半边登记两次以作标记,[refine_with_plane.h:493-494](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L493-L494))。

然后逐面处理:

- **跳过完全嵌在平面内的面**([refine_with_plane.h:504-514](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L504-L514));
- **两个半边的简单情形**(切割线正好横穿面,一次分割即可):
  `CGAL::Euler::split_face(h1, h2)` 生成切割边,并在 `edge_is_marked` 中标记为交线边([refine_with_plane.h:518-528](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L518-L528));
  若 `triangulate=true`,再判断分割出的面是否为三角形,不是则继续 `split_face` 补充对角线,保持全三角([refine_with_plane.h:530-554](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L530-L554));
- **多于两个半边的复杂情形**(切割折线在一个面内多次进出,常见于大 polygon 面):
  先按半边目标点的字典序排序并去除首尾重复的切线点([refine_with_plane.h:556-583](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L556-L583)),再两两配对依次 `Euler::split_face(h[i], h[i+1])`,完成整个面的重剖分([refine_with_plane.h:587-615](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L587-L615))。

阶段一结束后:

- 交线(平面 ∩ 网格)已成为网格中真实存在的边链,全部在 `ecm`(edge_is_marked_map)中标记为 `true`;
- 所有顶点的归属(`vos`)已确定;
- 断言 `is_valid_polygon_mesh(pm)` 通过([refine_with_plane.h:620](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L620))。

---

## 6. 阶段二:`split_along_edges` —— 复制交线边、形成边界

实现位于 [clip.h:525-662](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h#L525-L662)。目标:沿标记边把网格撕开,使两侧在切口处各自形成边界。半边数据结构中"撕开"="复制边 + 重连 next/prev + (必要时)复制顶点"。

### 6.1 收集与预处理([clip.h:537-569](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h#L537-L569))

- 收集所有 `ecm` 标记的边(`shared_edges`);
- 针对原本就带边界的网格:预先收集交线边端点周围的**既有边界半边**(`extra_border_hedges`),标记 `no_target_update`,并把相关顶点的出半边固定到边界半边上——这保证后续步骤不会把既有边界的信息冲掉。

### 6.2 复制边并重连([clip.h:572-617](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h#L572-L617))

对每条非边界的标记边 `e`(半边 `h`):

```cpp
user_visitor.before_edge_duplicated(h, tm);
new_hedge = halfedge(add_edge(tm), tm);   // 新增一条边(含一对半边 new_hedge/new_opp)
user_visitor.after_edge_duplicated(h, new_hedge, tm);

// new_hedge 顶替 h 在面环中的位置
set_next(new_hedge, next(h));  set_next(prev(h), new_hedge);
set_face(new_hedge, face(h));  set_halfedge(face(h), new_hedge);
set_target(new_hedge, vt);     set_target(new_opp, vs);

// h 与 new_opp 变成边界半edge(无面)
set_face(new_opp, null_face);  set_face(h, null_face);
```

效果:原来的一条内部交线边变成**两条几何重合的边界边**——一条属于正侧面的边界环,一条属于负侧面的边界环。端点是否需要复制视 `no_target_update` 而定,记入 `vertices_to_duplicate`。

### 6.3 缝合边界环([clip.h:620-634](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h#L620-L634))

新成为边界的半边(`hedges_to_update`)的 `next` 指针尚不正确。对每条这样的半边 `h`:

```cpp
candidate = opposite(prev(opposite(h)));      // 绕顶点扇面旋转
while (!is_border(candidate))                 // 直到找到下一条边界半边
  candidate = opposite(prev(candidate));
set_next(h, candidate);                       // 接上,闭合所在侧的边界环
```

即沿顶点周围的扇面(halfedge fan)"绕行",把同侧的边界半边按顺序首尾相接,两侧各形成一条闭合的切口边界环。

### 6.4 顶点复制([clip.h:636-659](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h#L636-L659))

对 6.2 中登记的需要复制的顶点(典型场景:切割线终止于网格**原有边界**的角点):

```cpp
vertex_descriptor nv = add_vertex(tm);
put(vpm, nv, get(vpm, p.second));       // 坐标相同
for (h : halfedges_around_target(p.first))
  set_target(h, nv);                    // 但把其中一侧的半边改指向新顶点
```

于是切口端点处的同一个几何点被拆成两个拓扑独立的顶点,每侧各执一个。最后一趟统一修正各边界顶点周围半边的 target,并以 `is_valid_polygon_mesh(tm)` 断言收尾([clip.h:661](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h#L661))。

阶段二结束后,平面两侧的几何只在切口处**位置重合**,拓扑上完全断开——它们成为同一 mesh 内的两个(或多个)连通分量。这就是文档所说 *"duplicated as to create a boundary, and thus separate connected components on either side of the plane"* 的确切含义。

---

## 7. 典型用法

```cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3                                     Point_3;
typedef K::Plane_3                                     Plane_3;
typedef CGAL::Surface_mesh<Point_3>                    Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main()
{
  Mesh mesh = /* 封闭水密网格 */;

  // 平面 z = d(kernel 必须与 mesh 顶点一致)
  Plane_3 plane(0, 0, 1, -d);

  // 切割:两侧都保留,结果仍在同一个 mesh 里
  PMP::split(mesh, plane);

  // 按连通分量拆成独立 mesh:parts[0]、parts[1] 即两半
  std::vector<Mesh> parts;
  PMP::split_connected_components(mesh, parts);

  // 可选:只想留某一半时,也可直接用 clip(mesh, plane, true);
  return 0;
}
```

常用变体:

```cpp
// 切割时对交线附近三角形做自交检查,失败抛异常
PMP::split(mesh, plane, PMP::parameters::throw_on_self_intersection(true));

// 输入是三角网格,但希望切割面保持为多边形(不额外三角化)
PMP::split(mesh, plane, PMP::parameters::do_not_triangulate_faces(true));
```

自定义 visitor(概念 `PMPCorefinementVisitor`,可跟踪新点/边/面的产生;注意回调中"平面一侧"以 `null_halfedge()`/`null_face()` 出现):

```cpp
struct MyVisitor : public CGAL::Polygon_mesh_processing::Corefinement::Default_visitor<Mesh>
{
  // 注意:Surface_mesh 的嵌套描述符名是 Vertex_index,没有 vertex_descriptor 嵌套名
  void new_vertex_added(std::size_t i_id, Mesh::Vertex_index v, const Mesh&)
  { std::cout << "切割产生新顶点 #" << i_id << "\n"; }
};

PMP::split(mesh, plane, PMP::parameters::visitor(MyVisitor()));
```

---

## 8. 注意事项与前置条件

1. **kernel 一致性**:平面必须与顶点映射中的点属于同一 kernel,否则编译期 `static_assert`(见 [refine_with_plane.h:290](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L290))或类型推导失败。
2. **输入不应自交**:同族函数(`clip` 等)文档要求 `!does_self_intersect(pm)`。对可靠性要求高的场合开启 `throw_on_self_intersection(true)`,注意它只对三角网格生效,且只检查交线附近的三角形。
3. **网格整体在平面一侧**:算法优雅短路,网格原样保留,`split` 后仍是一个连通分量,调用方需自行处理这种情况。
4. **退化情形**:顶点或边恰好落在平面上、面完全嵌入平面、切割线与既有边界相切/相交等,实现中均有专门处理(共面边判定、切线点双重登记、`vertices_to_duplicate` 等,见 5.2 / 5.5 / 6.4),但使用精确 kernel(如 `EPECK`)可获得更稳健的结果。
5. **切口处顶点成对**:切割线上的每个点在两侧各有一个副本,后续若做布尔运算或缝合(stitch)需留意;`PMP::stitch_borders()` 可把这样的重合边界重新缝起来。(此行为已在折面方案的两网格 split 测试中实测证实:切线路径上**全部**顶点——含内部拐点、内部共享边交点、三角化对角线交点——撕开后均两侧各一份,见关联文档 §5-3。)
6. **输出连通分量数量**:不一定是 2。平面可能同时穿过多个不相连的部件,或把一个部件切成多个分量;按实际连通分量数处理 `split_connected_components` 的结果。
7. **三角化**:输入为三角网格且未设 `do_not_triangulate_faces(true)` 时输出保持三角网格;输入非三角网格时,只有被切割穿过的面会被剖分,其余面保持原样。

---

## 9. 源码位置索引

| 内容 | 位置 |
|---|---|
| `split(pm, plane)` 文档注释 | [clip.h:1353-1423](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h#L1353-L1423) |
| `split(pm, plane)` 实现 | [clip.h:1424-1467](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h#L1424-L1467) |
| `split_along_edges()`(阶段二) | [clip.h:525-662](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h#L525-L662) |
| 复制边并重连 | [clip.h:572-617](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h#L572-L617) |
| 缝合边界环 | [clip.h:620-634](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h#L620-L634) |
| 顶点复制 | [clip.h:636-659](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h#L636-L659) |
| `refine_with_plane()` 文档注释 | [refine_with_plane.h:167-240](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L167-L240) |
| `refine_with_plane()` 实现(阶段一) | [refine_with_plane.h:241-621](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L241-L621) |
| 顶点分类 | [refine_with_plane.h:356-374](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L356-L374) |
| 边分类(穿越边/共面边) | [refine_with_plane.h:376-420](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L376-L420) |
| 自交检查 | [refine_with_plane.h:430-447](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L430-L447) |
| 穿越边插交点(`Euler::split_edge`) | [refine_with_plane.h:451-474](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L451-L474) |
| 面剖分(`Euler::split_face`) | [refine_with_plane.h:479-619](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/refine_with_plane.h#L479-L619) |
