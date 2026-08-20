# 用折面切割多边形:方案设计与验证文档

> 目标读者:需要核对方案正确性的开发者
> 源码版本:CGAL 6.2(`D:\github\CGAL-6.2\CGAL-6.2`)
> 关联文档:[cgal-pmp-split-with-plane.md](cgal-pmp-split-with-plane.md)(`split` 内部原理)
> **实施状态(2026-08-20)**:实现 [src/split_polygons.h](../src/split_polygons.h),测试 [tests/split_test.cpp](../tests/split_test.cpp),
> MSVC v145(VS 18 2026)Debug + Release 双配置,**59 项检查全部通过**(七用例:六切割用例 + 前置体检)。

---

## 0. 可信度声明(先读这个)

本文档论断分档标注;**2026-08-20 第 8 节六个测试全部通过后,原 [推断]/[未验证] 项已全部落定**:

| 标注 | 含义 | 核验方式 |
|---|---|---|
| **[已验证]** | 在本地 CGAL 6.2 头文件中读到原文(附文件:行号) | 点击行号链接自行查看 |
| **[已测试]** | 由第 8 节测试实际运行证实(2026-08-20,Debug+Release) | 复跑 `split_test` |
| **[已修正]** | 初版推断被实测推翻,已按实测改写(附推翻证据) | 见对应测试 |

完整可编译代码在 `src/split_polygons.h`(已随测试通过),第 7 节只摘录关键部分。

## 1. 问题定义

- **输入 A**:多边形 —— 三种形态(见第 3 节 ①):单个简单多边形点列 / 一批互不相交的多边形 / **共享边的相邻多边形拼成的网格**(每个面 = 一个输入多边形);点列为 3D,近似平面;
- **输入 B**:折面 —— 若干平面片拼成的分段线性曲面(可以折叠),可三角化;
- **精度策略**:**全程精确核 EPECK**(`Exact_predicates_exact_constructions_kernel`),谓词与构造均精确,交点坐标是精确有理数,无构造误差;
- **输出**:每片 = (来源多边形编号, 切割组编号, 该片边界环点列);来源由**逐面标记(三角化时)+ 标记传播(切割时)+ 按多边形边界的连通分组**追踪(第 3 节 ②②′④⑤)。**[已测试]**(测试 4:来源 0/0/1;测试 5/6:相邻多边形来源各自追踪 0/0/1/1)。

**方案**:不自己写求交/切割算法,把折面建成三角网格 `splitter`,调用 CGAL 两网格版 `split`,再按 (来源 × 切割片) 分组取回各片 —— **每个多边形各自被切成多片**。

## 2. 流水线总览

```
多边形(单/批/相邻共边)──① add_face──▶ mesh(每个原始面 = 一个多边形)
                                     │
                                     ├──② 来源标记 + 三角化:
                                     │      先给每个原始面写面号(= 多边形号),
                                     │      再带 visitor 调 triangulate_faces
                                     │      —— 每个新三角形继承所在原始面的标记
                                     │      (不能用连通分量记来源:相邻多边形共边,
                                     │        会并成一个分量,来源信息丢失)
                                     │
折面数据 ──③ 建网格+三角化──▶ splitter ──④ PMP::split(tm, splitter,
                                     │                visitor=标记传播)──┘
                                     │
                                     ├──⑤ 分组:两侧来源不同的边(= 相邻多边形
                                     │      共享边)设为约束边,connected_components
                                     │      不跨约束边扩散 → pid = (来源多边形 ×
                                     │      切割片) 组号
                                     │
                                     └──⑥ 组边界环:face(h) 属于本组且(对侧无面
                                            [撕开边] 或 对侧异组[多边形共享边])
                                            的半边,绕点旋转成环 ──▶ 各片多边形点列
```

## 3. 逐步依据

### ⓪ 精度策略:EPECK 精确核(全程)+ 依赖现状

- kernel 统一用 `CGAL::Exact_predicates_exact_constructions_kernel`(EPECK),`tm` 与 `splitter` **必须同一点类型/同一 kernel**(两网格求交要求几何 traits 兼容,统一 EPECK 最省事)。
- 精确性:谓词与构造全精确——corefine 求出的交点坐标是精确有理数,无浮点构造误差;退化情形结果可复现,第 6 节风险 3/6 的脆弱区显著缩小。
- 输出转换:EPECK 点坐标是 `K::FT`,交给下游 double 管线时逐坐标 `CGAL::to_double()`。
- **依赖现状(2026-08-20 实测环境)**:
  - CGAL 6.2:`D:\github\CGAL-6.2\CGAL-6.2`(Windows 发行版,`auxiliary/gmp` 自带 gmp/mpfr 的头文件 + .lib + .dll);
  - Boost 1.91:`D:\github\boost_1_91_0`(源码树 + b2 构建产物;本方案仅用其头文件,PMP 用到的模块均为 header-only);
  - 构建:CMake 4.3 + MSVC v145。**注意**:CMake 4.x 已移除 FindBoost 模块,而 Boost 源码树无 `BoostConfig.cmake`,发行版 `CGALConfig.cmake` 的 Boost 查找会失败 → 实际采用**手动 header-only 模式**(直接给 include 路径、`CGAL_HEADER_ONLY=1`、显式链 gmp.lib/mpfr.lib、运行时拷贝 gmp-10.dll/mpfr-6.dll),见第 8 节构建脚本。

### ① 多边形 → mesh(单面 / 批量 / 相邻共边三种形态)

**形态 1:单个多边形 → 单面**。`Surface_mesh::add_face(vertex_range)` 接受一个顶点环生成单个面;非凸简单多边形可以直接建面(环序即面序)。

**形态 2:一批多边形 → 多分量 mesh**。每个多边形各建一个面,全部放进同一个 `Surface_mesh`,一次 split 批量切割。条件:各多边形**互不相交**,否则违反 tm 无自交前置。**[已测试]**(测试 4)。

**形态 3:相邻多边形拼成的网格(共享顶点/边)**。每个面 = 一个输入多边形,相邻多边形共享边处**必须共享顶点建面**(同一 `add_face` 序列里复用顶点描述符)。管线不变:corefine 照常求交插点(交点可能落在共享边上);**每个多边形各自追踪来源、各自切成多片**(②′⑤⑥)。**[已测试]**(测试 5/6)。
**注意**:`build_mesh`(点列入口)给每个多边形建独立顶点,共边处顶点重复、边不共享 —— **相邻多边形不能用它**,必须手工建面(见测试 5/6 的构造方式)。

**共同前置:多边形必须简单(不自交)**。两网格 split 的文档前置明确要求无自交(**[已验证]** [clip.h:1255-1256](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/clip.h))。

**来源追踪**:批量/相邻输入时必须知道每片来自哪个多边形 → 逐面标记(②)+ 切割时的标记传播(④)+ 分组(⑤)。**[已测试]**。

### ② 来源标记 + 三角化:`triangulate_faces(tm, visitor=Tri_source_tagger)`

- 依据:两网格版 `split` 的文档写明 `@param tm input triangulated surface mesh`、`@param splitter triangulated surface mesh used to split tm`。**[已验证]**(clip.h:1263-1264)
- **triangulate_faces 有 visitor**(初版文档称"无 visitor"是**错的**):named parameter `visitor`,概念 `PMPTriangulateFaceVisitor`,基类 **[已验证]** [Triangulate_faces::Default_visitor](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/triangulate_faces.h#L46-L63)(triangulate_faces.h:46-63),钩子与 corefine 同名:`before_subface_creations(f_old)` / `after_subface_created(f_new)` / `after_subface_creations()`;quad 快速路径(triangulate_faces.h:288-321)与 hole-filling 路径(143 行起)均回调,quad 路径对复用的原面也回调 `after_subface_created`(写同值,幂等无害)。**[已修正]**(由本条方案实施实测:标记传播全程 unmarked==0)
- **标记流程**:三角化**前**给每个原始面写自己的面号(= 多边形号,插入序);visitor 在 `before_subface_creations(f_old)` 记录、`after_subface_created(f_new)` 继承 → 每个三角形知道自己来自哪个多边形。
- **与 corefine visitor 的差别**:triangulate_faces 只作用于单一网格、钩子**无** mesh 参数 → 无需 watch 守卫。visitor 同样被按值拷贝(triangulate_faces.h:430-431),`Surface_mesh::Property_map` 共享存储,安全。**[已验证]**
- **对角线交点(实测补充)**:三角化引入的对角线若与切线相交(如正方形对角线交 x=1 于 (1,1)),交点会成为网格顶点,并作为**共线冗余点**出现在片的边界环上(测试 1/5/6 实测:(1,1)、(1,0.5)、(0.5,0.5) 等)。它不影响片的多边形语义;下游要"纯角点环"时需做共线归一化(测试已采用,剥掉与前后邻点共线的顶点)。**[已测试]**

### ②′ 为什么不能用连通分量记来源(设计约束)

`PMP::connected_components` 按**共享边**连通(connected_components.h:1092-1095)。**相邻多边形有共用边** → 三角化后同属一个边连通分量 → 分量号无法区分它们(实测:两个共边四边形切一刀,旧方案两片同源;新方案 4 片 {0,0,1,1},测试 5)。因此来源标记必须**逐面**在三角化时完成(②),连通分量只用于切割片分组(⑤,且带多边形边界 barrier)。

### ③ 折面 → `splitter` mesh

- 与 ①② 相同方式建面 + `triangulate_faces`。
- 三个硬性要求:
  1. **无自交**:`!does_self_intersect(tm)`、`!does_self_intersect(splitter)` 是两网格 split 的文档前置条件。**[已验证]**(clip.h:1255-1256);
  2. **完全横跨多边形**:切割线必须从多边形边界通到边界,否则只形成内部裂缝(slit),仍是一个连通分量。**[已测试]**(测试 2 折面正确横跨,得 2 片;slit 反例未单测);
  3. **共享边的折面必须共享顶点建面**(如 V 形两坡共屋脊):若两坡各自独立建面,屋脊棱上两组顶点坐标重合但拓扑分离,行为未定义。测试 2 的 splitter 用共享屋脊顶点构造。**[已测试]**
- 开放曲面可以当 splitter:两网格 split 的文档没有封闭性要求。**[已测试]**(测试 1-6 的 splitter 均为开放曲面)。

### ④ 切割:`PMP::split(tm, splitter, np_tm, np_s)`

- 签名与实现 **[已验证]**(clip.h:1307-1351),实现只有两步:corefine(把交线边写入 `ecm`)+ `internal::split_along_edges`(撕开)。
- **标记传播(来源追踪的关键)**:corefine 切分一个面时的调用模式 **[已验证]**(face_graph_utils.h:587-629):`before_subface_creations(f_old)` 带被切面调用一次 → **第一个子面复用 `f_old` 的面描述符**(标记天然保留)→ 其余每个新面触发 `after_subface_created(f_new)`。两钩子"记录/继承"即可让标记贯穿切割。**[已测试]**(测试 4/5/6:全部面标记非 -1)。
- **visitor 必须按网格守卫(初版文档缺失,实测抓出)**:corefine 会**同时细化 tm 和 splitter 两个网格**,同一 visitor 对两个网格的切分都会触发;钩子的 `const TriangleMesh&` 参数就是用来区分的。不守卫会把 splitter 的面索引写进 tm 的属性表——轻则写花标记,重则 `Properties.h` 数组越界断言崩溃。守卫方式:visitor 记 `watch` 指针,`&m == watch` 才操作。**[已修正]**(无守卫版本实测崩溃于 `_idx < data_.size()`,加守卫后全过)
- **`throw_on_self_intersection` 不要传**:两网格 `split()` 实现只提取 vertex_point_map / visitor / do_not_modify,不提取也不转发该参数(clip.h:1325-1347),传了被静默忽略。**[已修正]**(初版文档建议传它是错的)
- `do_not_modify`(np_s):不改动折面。实现里 splitter 按值传入 `split_polygons`,内部复制后切割,折面原件不受影响,故未用该参数。
- visitor 会按值拷贝(clip.h:1337);`Surface_mesh::Property_map` 是共享底层存储的轻句柄,拷贝安全、写入互通。**[已验证]**

### ⑤ 分组:`connected_components(tm, pid, edge_is_constrained_map=多边形边界)`

- **barrier 语义**:**[已验证]** [connected_components.h:149-152](../../github/CGAL-6.2/CGAL-6.2/include/CGAL/Polygon_mesh_processing/connected_components.h#L149-L152) —— 扩散前检查 `if(!get(ecmap, edge(h)))`,**约束边阻断扩散**。把"两侧面来源不同的边"(= 相邻多边形的共享边,撕开后仍共享)设为约束边 → 撕开分量内再按多边形隔开 → `pid` = **(来源多边形 × 切割片) 组号**:**每个多边形各自分成多片**。
- 撕开后的切割边一侧无面,天然隔开不同切割片(论证 5-3);共享边被 barrier 隔开不同多边形 —— 两类边界都不跨,组 = 切割片 ∩ 单个多边形。**[已测试]**(测试 5:共边两四边形 4 组;测试 6:经未切共享边连通的两片正确分作两组)
- **连通性按"共享边"定义,不是共享顶点**(**[已验证]** connected_components.h:1092-1095)。见第 5 节论证 3。
- 备选:若确需每片独立 mesh,分组之后再调 `PMP::split_connected_components`(**[已验证]** connected_components.h:1101-1103);其内部经 `copy_face_graph` 复制,自定义面属性是否带出仍未验证,故分组逻辑不依赖它。

### ⑥ 转回多边形:组边界半边 + 绕点旋转成环

- **组边界半边** `h`:`face(h) ∈ 组` 且(`face(opposite(h))` 无面[撕开边/网格外边界] 或 组号不同[相邻多边形共享边])。它覆盖两类组边界:切割弧与多边形共享边 —— **测试 5/6 的环都含共享边段**(如测试 6 中 S1 上片的环含 x=1, y∈[0.5,1] 段,即与 S2 上片的交界)。
- **走法必须"绕点旋转",不能沿 next() 链**(实测抓出的 bug):从边界半边 h 的终点 v 取 `g = opposite(h)`(自 v 出发、指向 h 的起点,必在组外),绕 v 逐步旋转 `g = next(opposite(g))`,直到 g **自身**是组边界半边(面在组内**且**对侧不在组内——只判"面在组内"会停在组内对角线上)。沿 `next()` 链走是错的:next 不换面,跨不过组内共享边,会把环撕成每面一段的碎片(实测:2 个三角形的正方形被拆成两个 2 点"环")。**[已修正]**
- 旋转步 `next(opposite(g))` 在 Surface_mesh 上对普通半边(对面内)与边界半边(沿 border 环)都正确:两者都得到仍自 v 出发的下一条半边。**[已验证]**(Surface_mesh 边界半边维护 next 链)
- 每组是盘状(论证 5-4)→ 恰 1 个环;环上顶点按序即多边形点列;方向与面顶点序一致。**[已测试]**
- 环中可能含共线冗余点(见 ②)与输入原有的共线顶点;要纯角点环需共线归一化。**[已测试]**
- 环数 > 1 的组(带洞/内岛):输入本身带洞或封闭折面切出环形片时出现,同 (source, piece) 会输出多个环;开放折面切无洞多边形不会出现。

## 4. API 核实清单

| API | 位置 | 状态 |
|---|---|---|
| `PMP::split(tm, splitter, np_tm, np_s)` | clip.h:1310 | 已验证(签名+实现:corefine+撕开 1345-1350)+ 已测试 |
| `PMP::corefine` | corefinement.h | 已验证(被 clip.h 包含;会细化两个网格,visitor 需守卫) |
| `internal::split_along_edges` | clip.h:526 | 已验证(逐行读过;internal 无公开替代 → 撕开只有 `split` 一条公开路径) |
| `PMP::triangulate_faces(pmesh, np)` + visitor | triangulate_faces.h:448-469, 519-524 | 已验证 + 已测试(来源标记继承,§3②) |
| `Triangulate_faces::Default_visitor` 钩子 | triangulate_faces.h:46-63 | 已验证 + 已测试(单网格,无需守卫) |
| `PMP::connected_components(pm, fcm, np)` + `edge_is_constrained_map` | connected_components.h:175-217 | 已验证(barrier 语义 149-152)+ 已测试(分组) |
| `PMP::degenerate_faces` | shape_predicates.h:272 | 已验证 + 已测试(preflight) |
| `PMP::does_self_intersect` | self_intersections.h:774 | 已验证 + 已测试(preflight;两者均要求三角网格) |
| `PMP::split_connected_components` | connected_components.h:1102 | 已验证(签名;未使用) |
| `Corefinement::Default_visitor` 钩子 | face_graph_utils.h:417-487 | 已验证 + 已测试(befor/after_subface + after_face_copy,需 watch 守卫) |
| `Surface_mesh::add_property_map` / `Property_map` | Surface_mesh/Properties.h | 已测试(共享句柄,拷贝互通;越界有断言) |
| ~~`CGAL::extract_boundary_cycles`~~ | boost/graph/border.h:247 | 已弃用(⑥ 改为组边界旋转成环后不再使用) |

### ⓢ 前置体检 `preflight_check(tm, splitter)`(诊断"正常 mesh 也挂")

split 内部即 corefine + 撕开(clip.h:1345-1350),**改调 `PMP::corefine` 躲不开任何前置**,且 corefine 不撕开、撕开(`internal::split_along_edges`)无公开 API。split 又不消费 `throw_on_self_intersection`(§3④)—— 前置被违反时**不抛异常、直接未定义行为**(常表现为崩在 corefine 深处)。`preflight_check` 在进入 corefine 前逐条给出违规原因,检查对象语义与管线一致:

- **tm**:管线会先三角化(②),多边形面合法 → 在"三角化后的副本"上体检;三角化后仍非三角 = 退化环;
- **splitter**:调用者必须传已三角化的(`build_splitter` 已保证);
- 全三角通过后才查退化面(`degenerate_faces`)与自交(`does_self_intersect`,最贵放最后)。

**[已测试]**(测试 7:正常输入 0 违规;自交/退化面/非三角 splitter 各被准确抓出)。

## 5. 核心论证链:为什么"按连通性转回"是正确的

1. **交线成为真实边**:`corefine` 保证 tm 与 splitter 的交线全部成为 tm 的边并写入 `ecm`。**[已验证]**(clip.h:1345-1346 传参即此意图)
2. **交线边全部被复制**:`split_along_edges` 复制每一条 `ecm==true` 的边,复制后两侧不再共享任何交线边。**[已验证]**(clip.h:572-617)
3. **切割片无共享边;相邻多边形共享边保留**:输入各多边形是盘状,被"边界到边界"的切割弧分成若干盘状片;片与片只能沿切割弧相邻,而切割弧的边已全部复制 → 不同切割片无共享边。相邻多边形的共享边**不在**切割弧上时不被复制,两侧仍共享 → ⑤ 用 barrier(约束边)隔开,组 = 切割片 ∩ 单个多边形。**[已测试]**(六个测试 `component_count == pieces.size()`;测试 5/6 相邻片经共享边连通仍正确分组)。
   - **顶点复制细节(实测修正)**:初版推断"只有落在原有边界上的端点被复制、内部拐点两侧共享"是**错的**。实测(测试 1/2/5/6 的顶点数与环内容):切割路径上的**全部**顶点(原边界端点、内部拐点 c、多边形共享边交点、三角化对角线交点)撕开后均两侧各一份。这不影响分组——连通性按**边**定义,共享与否都不合并片。**[已修正]**
4. **每组恰有一个边界环**:组 = 盘状多边形 ∩ 盘状切割片,仍盘状 → 边界是一个环 → ⑥ 对每组返回 1 个环,环上顶点依序即该片多边形。**[已测试]**
5. **来源标记贯穿全程**:原始面写面号(②)→ 三角化新面继承(`Triangulate_faces::Default_visitor` 钩子)→ corefine 切面时第一个子面复用原面描述符、其余新面经 `after_subface_created` 继承(**[已验证]** 钩子调用模式)→ visitor 按网格守卫(④)→ 撕开不移动面描述符 → 组内任一面标记即该片来源。**[已测试]**(测试 4:来源 0/0/1;测试 5/6:相邻多边形来源 {0,0,1,1};无 -1 残留)
6. **组边界环走法正确**:绕点旋转每步在有限个自 v 出发的半边内循环,遇到"面在组内且对侧不在组内"即停 —— 该半边就是环上下一段(组在 v 处占一个连续扇区,旋转从组外扇区出发必先碰到该扇区的结束边);沿环一周回到起点闭合。实测与手算环全部吻合(测试 1-6)。**[已测试]**

## 6. 风险与失效模式

| # | 情形 | 后果 | 对策 | 状态 |
|---|---|---|---|---|
| 1 | 折面未横跨多边形 | 内部裂缝(slit),仍 1 片 | 折面延伸出多边形包围盒 | 理论推断,未单测 |
| 2 | 折面自交 / 多边形自交 | 违反前置,结果不可预测 | 构造时保证 | 前置(已验证) |
| 3 | 切线恰过顶点、与边共线 | 退化情形,EPECK 下可复现 | 测试覆盖 | EPECK 已缓解 |
| 4 | 精度 | 全程 EPECK:交点精确 | 已执行 | **[已测试]** |
| 5 | 封闭折面切出内岛/环形片 | 该组 2+ 个边界环 | 同 (source, piece) 输出多环 | 未单测 |
| 6 | 折面与多边形共面重叠 | coplanar 情形,corefine 有专门处理 | 尽量避免 | 未单测 |
| 7 | 三角化对角线与切线相交 | 环上出现共线冗余点 | 下游共线归一化(测试已采用) | **[已测试]** |
| 8 | corefine visitor 不按网格守卫 | splitter 面索引写进 tm 属性表,越界断言崩溃 | `watch` 指针守卫(④) | **[已修正]** |
| 9 | 折面/多边形共享边不共享顶点 | 棱上坐标重合拓扑分离,行为未定义 | 共享顶点建面(③-3、① 形态 3) | **[已测试]**(正例) |
| 10 | 用连通分量记来源(相邻输入) | 共边多边形并成一个分量,来源丢失 | 逐面标记(②),分量仅用于分组且带 barrier(⑤) | **[已测试]**(测试 5/6) |
| 11 | 组边界环沿 next() 链走 | next 不换面,环被撕成每面一段的碎片(实测 2 点"环") | 绕点旋转 `g = next(opposite(g))`(⑥) | **[已修正]** |
| 12 | splitter 未三角化 / tm 环退化 / 输入自交 | split 不抛异常,崩在 corefine 深处难定位 | `preflight_check`(§4ⓢ)预诊断 | **[已测试]**(测试 7) |

## 7. 实现(已编译,已通过第 8 节全部测试)

完整代码:[src/split_polygons.h](../src/split_polygons.h)(`poly_split::split_polygons(polygons, splitter)` / `split_mesh(tm, splitter)` / `preflight_check(tm, splitter)`)。关键部分:

```cpp
namespace poly_split {
namespace PMP = CGAL::Polygon_mesh_processing;

// 三角化标记 visitor(②):单网格、钩子无 mesh 参数 → 无需 watch 守卫。
// quad 路径对复用原面也回调 after_subface_created(写同值,幂等)。
template <class TagMap, class TriangleMesh>
struct Tri_source_tagger : PMP::Triangulate_faces::Default_visitor<TriangleMesh>
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  TagMap fid;
  int cur = -1;
  void before_subface_creations(face_descriptor f_old)
  { cur = get(fid, f_old); }                    // 被拆面:记录
  void after_subface_created(face_descriptor f_new)
  { put(fid, f_new, cur); }                     // 新三角形:继承
};

// 切割传播 visitor(④):corefine 同时细化两个网格,必须 watch 守卫
template <class TagMap, class TriangleMesh>
struct Tag_propagator : PMP::Corefinement::Default_visitor<TriangleMesh>
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  TagMap fid;
  const TriangleMesh* watch = nullptr;
  int cur = -1;

  void before_subface_creations(face_descriptor f_old, const TriangleMesh& m)
  { if (&m == watch) cur = get(fid, f_old); }
  void after_subface_created(face_descriptor f_new, const TriangleMesh& m)
  { if (&m == watch) put(fid, f_new, cur); }
  void after_face_copy(face_descriptor f_old, const TriangleMesh& m1,
                       face_descriptor f_new, const TriangleMesh&)
  { if (&m1 == watch) put(fid, f_new, get(fid, f_old)); }
  // 第一个子面复用 f_old 描述符(face_graph_utils.h:614-619),标记天然保留
};

// 核心管线(节选,错误处理/类型别名见 src/)
template <class K>
Split_result<K> split_mesh(Polygon_mesh<K> tm, const Polygon_mesh<K>& splitter)
{
  // ② 来源标记 + 三角化:原始面写面号(= 多边形号),visitor 让新三角形继承
  auto fid = tm.add_property_map<Mesh::Face_index, int>("f:source", -1).first;
  { std::size_t i = 0; for (auto f : faces(tm)) put(fid, f, int(i++)); }
  PMP::triangulate_faces(tm,
      PMP::parameters::visitor(Tri_source_tagger<decltype(fid), Mesh>(fid)));

  // ④ 切割 + 标记传播(watch 守卫;不传 throw_on_self_intersection,§3④)
  Tag_propagator<decltype(fid), Mesh> tagger(fid, tm);
  Mesh splitter_copy(splitter);
  PMP::split(tm, splitter_copy,
             CGAL::parameters::vertex_point_map(get(CGAL::vertex_point, tm))
                              .visitor(tagger));

  // ⑤ 分组:两侧来源不同的边设为约束边 → connected_components 不跨它扩散
  //    (connected_components.h:149-152)→ pid = (来源多边形 × 切割片) 组号
  auto pid = tm.add_property_map<Mesh::Face_index, std::size_t>("f:piece", 0).first;
  auto ecm = tm.add_property_map<Mesh::Edge_index, bool>("e:source_border", false).first;
  for (auto e : edges(tm)) { /* 两侧 fid 不同 → put(ecm, e, true) */ }
  component_count = PMP::connected_components(tm, pid,
      PMP::parameters::edge_is_constrained_map(ecm));

  // ⑥ 组边界环:边界半边 = face(h)∈组 且(对侧无面或对侧异组);
  //    绕点旋转 g=next(opposite(g)) 找下一段(沿 next 链是 bug,§3⑥)
  //    每组一环,环点列 = 依次 target(边界半边)
}

// ⓢ 前置体检:返回违规清单(空 = 可安全进入 split),见 §4ⓢ
} // namespace poly_split
```

## 8. 测试计划与结果(2026-08-20 全部通过)

测试程序 [tests/split_test.cpp](../tests/split_test.cpp):`split_polygons`/`split_mesh` 六切割用例 + 前置体检用例 + 环比较(共线归一化 + 旋转/反向匹配)+ 通用校验(unmarked==0、component_count==pieces、is_valid)。**Debug + Release 均 59 项检查 0 失败。**

| # | 用例 | 关键断言(实测值) | 验证点 |
|---|---|---|---|
| 1 | 竖直平板 x=1 切正方形 | 2 片;环归一化后为两个 4 边形;V=10(4+3 交点+3 复制,含对角线交点 (1,1));F=6 | 基本链路;开放单面网格可用 |
| 2 | V 形折面切正方形(折痕穿内部) | 2 片;环恰为预测的 4 边形 {(0,1.8),(0,2),(1.8,2),(0.8,1)} 与 6 边形;V=10(内部拐点 c 也复制) | 折痕穿面内部;交线 A→c→B 精确落位 |
| 3 | 折面远离(x=10) | 1 片;环 = 原多边形;V=4 | 无交优雅通过 |
| 4 | 批量:S1 切、S2 不切 | 3 片;来源 **0/0/1**;unmarked=0;V=14(8+3+3) | 多分量输入;**来源标记贯穿切割** |
| 5 | **相邻**两四边形(共享边 y=1)+ 切线 x=1 穿共享边 | **4 片;来源 {0,0,1,1}**;每片环含共享边段;V=16(6+5 交点+5 复制) | 相邻多边形各切成两片;**逐面标记 + barrier 分组**;共享边交点复制 |
| 6 | **相邻**两方块(共享边 x=1)+ 横墙 y=0.5 同时切两者 | **4 片;来源 {0,0,1,1}**;S1上/S2上 经未切共享边连通仍分作两组,环含该共享边段;V=16(6+5+5) | barrier 隔开同分量内不同多边形;**共享边作为组边界出环** |
| 7 | 前置体检 preflight | 正常输入(quad tm + 三角墙)0 违规;自交 tm 被抓;共线三角形被抓;quad 放 splitter 位置被抓 | **split 前置违反的预先诊断**(§4ⓢ) |

与初版手算的差异(均已按实测修正):测试 1/4/5 初版顶点数(8/12/11)漏算了**三角化对角线与切线的交点**和**路径内部顶点的复制**;"只复制原边界端点"的推断被测试 1/2/5 推翻(§5-3);"triangulate_faces 无 visitor"与"组边界沿 next 链走"两处初版论断被本次实施推翻(§3②/§3⑥)。

### 构建与运行(实际使用,手动 header-only 模式)

```cmake
# CMakeLists.txt(节选,完整见仓库根)
set(CGAL_DIR  "D:/github/CGAL-6.2/CGAL-6.2"   CACHE PATH "CGAL 6.2 发行版根目录")
set(BOOST_DIR "D:/github/boost_1_91_0"        CACHE PATH "Boost 源码树")
set(GMP_DIR   "${CGAL_DIR}/auxiliary/gmp"     CACHE PATH "GMP/MPFR")
target_include_directories(split_test PRIVATE src "${CGAL_DIR}/include" "${BOOST_DIR}" "${GMP_DIR}/include")
target_compile_definitions(split_test PRIVATE CGAL_HEADER_ONLY=1)
target_link_libraries(split_test PRIVATE "${GMP_DIR}/lib/gmp.lib" "${GMP_DIR}/lib/mpfr.lib")
# MSVC 需 /bigobj /utf-8;运行时把 auxiliary/gmp/bin 下 gmp-10.dll、mpfr-6.dll 拷到输出目录
```

```powershell
cmake -B build -G "Visual Studio 18 2026"
cmake --build build --config Release
.\build\Release\split_test.exe    # 59 checks, 0 failures
```
