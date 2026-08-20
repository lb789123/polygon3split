# polygon3split

用**折面**(分段线性曲面)切割多边形组,输出各切割片的 (来源编号, 切割组编号, 边界环)。支持**相邻(共享边)多边形**输入,每个多边形各自被切成多片。全程精确核 **EPECK**,基于 CGAL 6.2 `PMP::split`(两网格版)。

方案设计与逐条验证:**[docs/polygon-split-by-folded-surface.md](docs/polygon-split-by-folded-surface.md)** ·
split 内部原理(平面版):**[docs/cgal-pmp-split-with-plane.md](docs/cgal-pmp-split-with-plane.md)**

## 目录

```
src/split_polygons.h     核心实现(header-only,poly_split 命名空间)
tests/split_test.cpp     六用例验证测试(期望值全部手算)
CMakeLists.txt           手动 header-only 构建(依赖路径见下)
docs/                    设计/验证文档 ×2
```

## 用法

```cpp
#include "split_polygons.h"
typedef CGAL::Exact_predicates_exact_constructions_kernel K;

// 形态 A:互不相交的多边形点列组
auto result = poly_split::split_polygons<K>(polygons, splitter);
// polygons: std::vector<std::vector<K::Point_3>>,每个点列环一个多边形(互不相交、不自交)

// 形态 B:相邻(共享边)多边形 —— 手工建面,共享顶点,每个面 = 一个多边形
poly_split::Polygon_mesh<K> tm;
/* ... add_vertex / add_face(相邻多边形复用共享边的顶点) ... */
auto result = poly_split::split_mesh<K>(tm, splitter);

for (const auto& piece : result.pieces) {
  piece.source;   // 来自哪个输入多边形(②逐面标记,三角化时即写入)
  piece.piece;    // 切割组编号(⑤ 来源×切割片)
  piece.ring;     // 该片边界环点列(K::FT 精确坐标,含与相邻多边形的共享边段)
}
// 校验字段:result.component_count == pieces.size() 且 result.unmarked_faces == 0 为正常

// 崩溃诊断:split 违反前置(自交/退化面/splitter 未三角化)时不抛异常、直接未定义
// 行为 —— 先跑 preflight,拿到具体违规原因再进 split:
auto problems = poly_split::preflight_check<K>(tm, splitter);   // 空 = 可安全切割
```

注意:`build_mesh`/`split_polygons` 给每个多边形建独立顶点,**相邻多边形不能用**,须手工建面共享顶点(见 tests 用例 5/6)。

## 环境与构建

依赖(路径可在 CMake 里 -D 覆盖):

| 依赖 | 位置 |
|---|---|
| CGAL 6.2(发行版,含 gmp/mpfr) | `D:/github/CGAL-6.2/CGAL-6.2` |
| Boost 1.91(纯头文件使用) | `D:/github/boost_1_91_0` |
| 工具链 | CMake 4.3 + Visual Studio 18 2026(MSVC v145) |

采用**手动 header-only 模式**(`CGAL_HEADER_ONLY=1`,显式链 `auxiliary/gmp/lib/{gmp,mpfr}.lib`,运行时自动拷贝 `gmp-10.dll`/`mpfr-6.dll`)。
原因:CMake 4.x 已移除 FindBoost 模块,而 Boost 源码树无 `BoostConfig.cmake`,发行版 `CGALConfig.cmake` 的 Boost 查找会失败 —— 详见设计文档 §8。

```powershell
cmake -B build -G "Visual Studio 18 2026"
cmake --build build --config Release
.\build\Release\split_test.exe
```

## 测试状态

**2026-08-20:Debug + Release 双配置,59 项检查 0 失败。**

| # | 用例 | 验证点 |
|---|---|---|
| 1 | 竖直平板切正方形 | 基本链路;开放单面网格可用于两网格 split |
| 2 | V 形折面,折痕穿多边形内部 | 面内插点;交线 A→c→B 精确落位 |
| 3 | 折面远离不相交 | 无交优雅通过 |
| 4 | 批量输入(两分量,一切一不切) | 来源标记贯穿切割(0/0/1) |
| 5 | 相邻两四边形,切线穿共享边 | 相邻多边形各切成两片({0,0,1,1});逐面标记 + barrier 分组 |
| 6 | 相邻两方块,横墙同时切两者 | 同分量内按多边形分组;共享边作为组边界出环 |
| 7 | 前置体检 preflight | 自交/退化面/非三角 splitter 预先诊断 |

## 已知注意点(实测确立,详见设计文档)

- **来源标记必须逐面打**(三角化时带 visitor):相邻多边形有共用边,连通分量会把它们并成一个分量,来源信息丢失;`triangulate_faces` 有 visitor(`Triangulate_faces::Default_visitor`,单网格、无需守卫);
- corefine 会**同时细化两个网格**:其 visitor 钩子必须按 `const TriangleMesh&` 参数守卫,否则会把 splitter 的面索引写进 tm 的属性表(越界断言崩溃);
- **组边界环要用绕点旋转**(`g = next(opposite(g))`),不能沿 `next()` 链 —— next 不换面,跨不过组内共享边,环会被撕成碎片;
- 切割路径上**全部**顶点撕开后均两侧各一份(含内部拐点、三角化对角线交点);环上会出现共线冗余点,要纯角点环需共线归一化;
- 两网格 `split()` 不消费 `throw_on_self_intersection`(传了被静默忽略);
- Surface_mesh 嵌套描述符名是 `Vertex_index`/`Face_index`/`Halfedge_index`(没有 `*_descriptor`)。
