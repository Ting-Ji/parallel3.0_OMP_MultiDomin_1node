# AI 协作记录:多域预条件器改进

本文件是 **DeepSeek Harness(代码侧,可读写仓库)** 与 **ChatGPT(理论/复核侧,无仓库访问)** 之间的交流窗口与长期记录。

## 使用协议

1. **信息流向**
   - 代码侧(我)在「第 1 部分」追加条目,并把需要在「第 2 部分」回答的问题编号化。
   - 人类把本文件的相关段落(或整个文件)复制给 ChatGPT。
   - 人类把 ChatGPT 的回答粘贴回「第 2 部分」对应编号下。
   - 我读取回答,在文末「变更记录」登记处理结果,必要时追加新问题。
2. **编号约定**
   - 事实条目:`F-<序号>`(Fact,已逐行读代码确认)
   - 判断条目:`J-<序号>`(Judgement,推断/未实测)
   - 提问条目:`Q-<序号>`(Question,需要 ChatGPT 回答)
   - 回答条目:`A-<序号>`(Answer,人类粘贴)
   - 回答必须标明它对应哪个 `Q` 编号;新增问题从当前最大编号继续。
3. **证据规范**
   - 引用代码一律用 `相对路径:行号`(相对于本仓库根目录),不用行号区间以外的模糊描述。
   - 任何"我认为"必须标成 `J-`,不得混入 `F-`。
   - 未被任何一次运行记录过的数值,一律不许当作结论。
4. **当前范围约定**
   - 默认只讨论**多域路径**(`MultiDomainModel` + `GlobalDofMap` + `CCSRMat` + `BuildMappedLeafPreConditioner`)。
   - 单域 `GMRESPreConditioner` 只作为兼容入口与对照,不是改进对象。
   - `DomainCount()==1` 但走多域代码路径的算例,允许作为"隔离界面效应"的对照。
   - 详见 `README.md` 开头「当前工作约定:默认只考虑多域路径」。

## 索引

| 会话 | 日期 | 主题 | 状态 |
| --- | --- | --- | --- |
| S1 | 2026-02(本轮) | 多域叶块预条件器结构、缺块成因、创新路线 | 待 ChatGPT 复核 `Q-1`~`Q-7` |

---

# 第 1 部分:代码侧 → ChatGPT

## S1 背景与结论摘要

研究对象:三维时域边界元(TDBEM)常单元程序,多材料子域 + 匹配界面耦合,每时间步解块稀疏系统 `A_n x = b_n`(CCSR 存 3×3 块),用左预条件 GMRES。

### 一句话结论

`missingLocalBlocks` 的**主导成因不是分域或界面**,而是"波前带直径远小于叶块直径";分域只在"一个叶混进不同域的点"时产生缺块;界面不产生缺块,它只改变"哪一行配哪个列"。此外,预条件器只用 step 0 的矩阵构造一次后全程复用,而波前带随步数单调变宽,这可能比缺块本身更影响收敛。

## S1-F 已确认事实(逐行读代码)

### F-1 预条件器代数形式

- 位置:`DBEM1/multidomain.cpp:4186`(`BuildMappedLeafPreConditioner`)、`DBEM1/pre_gmres.cpp:822`(`GMRES_M_X`)。
- `M = blockdiag(A_L)`,`L` 为空间树(八叉树)叶节点;`A_L = [A(r, p(c))]_{r,c∈L}`,元素为 3×3 块。
- 构造期对每个叶块求稠密逆并存于 `m_PreM[leaf]`;应用期 `M⁻¹` 作用 = 按 `m_InputPID` gather → 乘叶块逆 → 按 `m_OutputPID` scatter。
- 归类:**带方程行—未知量列置换的非重叠块 Jacobi(局部近似逆)预条件器**。非 ILU,无跨叶消元,无粗空间。

### F-2 行块与列块是两套编号

- 行块 = 物理点号,恒等(`m_equationRows[node] = node`)。
- 列块 = 未知量块号,按 BCID/界面类型分配。规则:外边界 `BCID=123` 的行 → 该点未知面力块;`BCID=456` → 未知位移块;界面配对 A 侧行 → 共享位移块,B 侧行 → 共享面力块。
- `p(r) = GlobalDofMap::GetPreferredUnknownBlockForRow(r)`;`Validate()` 强制 `p` 为双射(`DBEM1/multidomain.cpp:1017-1025`)。
- 映射写入点:`DBEM1/multidomain.cpp:4303-4308`(`m_InputPID[i]=m_RePID[i]`,`m_OutputPID[i]=p(m_RePID[i])`)。

### F-3 索引约定

- `IID[p] = 3*p`,其中 `p` 是**叶内局部下标**(不是全局下标)。
- `m_RePID[g]` 为"重编号序 → 原物理点号";叶内点连续存放(`RenumberPointID`,`DBEM1/tree.cpp:431`)。
- 单域路径**不分配** `m_InputPID/m_OutputPID`;`GMRES_M_X` 以 `?:` 回落到 `m_RePID`(`DBEM1/pre_gmres.cpp:839`)。
- `GMRESPreconditionerOutputUnique()`(`DBEM1/pre_gmres.cpp:763`)先验证 `m_OutputPID` 在全体叶块上双射,成立才走 `#pragma omp parallel for` 的免原子 scatter;否则退串行。

### F-4 缺块的检测点

```cpp
if (FindCCSRBlockInRow(matrix, row, col) < 0) ++missingBlocks;
AssignCCSRBlockToLocalPreconditioner(matrix, row, col, pre.m_PreM[leaf], r, c);
```

- 位置:`DBEM1/multidomain.cpp:4323`、`DBEM1/multidomain.cpp:3912`。
- 找不到块时写入 3×3 **零块**(不跳过、不触发重组)。`FindCCSRBlockInRow` 为行内线性扫描(`DBEM1/multidomain.cpp:3748`)。
- 含义:预条件器只读取叶内条目,跨叶耦合**从不查询**;缺块即"预条件器把矩阵里没有的块当零"。

### F-5 一个 3×3 块进入 A 的 4 个必要条件

| # | 条件 | 位置 |
| --- | --- | --- |
| 1 | `NeedCal()` 返回非 0 | `DBEM1/dsquareelement.cpp:1157`、判据在 `:1184` |
| 2 | 上层 `Compute*Block` 的 `if (!flag) return 0` 未短路 | `DBEM1/multidomain.cpp:1883`、`:2489` |
| 3 | 未被阈值剔除:`IsZeroBlock9 \|\| IsNearZeroBlock9`(绝对阈值 1e-14) | `DBEM1/multidomain.cpp:1145` |
| 4 | 源单元与场单元属于**同一物理域** | 串行 `DBEM1/multidomain.cpp:2517-2563`(源/场循环都遍历同一 `domain.elementIds`) |

`NeedCal` 判据原文要点:

```cpp
double r = Dist(m_center, o_Point);          // 场单元中心 ↔ 源点(=源单元中心)
double Rmax = r + m_radius, Rmin = max(0.0, r - m_radius);
if (ndt >= Rmax / C2 || np1dt <= Rmin / C1) return 0;   // ndt = n*dt, np1dt = (n+1)*dt
```

即"两单元球壳与当前时滞时间窗无交集 ⇒ 整块不组装"(时间因果/波前截断)。

### F-6 界面的编码方式

- 界面两侧**共享同一个未知块**:`m_uRefs[nodeA]` 与 `m_uRefs[nodeB]` 指向同一 `globalBlock`(位移 `sign=+1`;面力 B 侧 `sign=normalSign`,通常 −1)。位置:`DBEM1/multidomain.cpp:831-848`。
- 界面单元的装配是"同一行同时尝试 `BCID=456` 与 `123` 两个版本",经列变换 `TransformColumnBlock9` 后各写一个块(`DBEM1/multidomain.cpp:2538-2555`)。
- **不存在** `(A域行, B域列)` 形式的跨域块。

### F-7 构造时机与 GMRES 参数

- 预条件器在时间循环**之前**、只从 step 0 的 `A` 构造一次,之后所有时间步复用:`DBEM1/multidomain.cpp:4654-4661`。
- 多域调用 `gmres(elements, A, rhs, pre, iterations, error, /*restart=*/1, x0, x, flag, iter)`:`DBEM1/multidomain.cpp:4757`。
- 收敛判据为真实残差 `||b-Ax||`;每步 `x0 = 0`(`:4751`)。

### F-8 树的实现细节

- 根:`DBEM_MULTIDOMAIN_TREE_ROOT=cuboid`(默认)走 `AssignPaddedBoundingBox(...,1.02)`,`L_k = 1.02·max(Δ_k, ε_L)`,`ε_L = 64·eps_mach·max(1, max_k Δ_k, max|x|)`(`DBEM1/geo.cpp:171`);`cube` 用最长跨度作三方向统一边长。
- 每层沿 x/y/z 各二分(8 子);落点按子节点编号 0..7 顺序,分割面上的点确定性进入编号最小的可含子节点;舍入导致无处归属时按"子框边长归一化的中心距离"选最近子节点并计入 `fallbackAssignments`(`DBEM1/tree.cpp:203-241`)。
- 两类不设点数上限的强制叶:坐标跨度低于数值分辨率 → `degenerateLeaf`;到达第 64 层 → `maxDepthLeaf`(`DBEM1/tree.cpp:161-181`)。
- 多域另有严格自检:根框有限为正、所有点在根内、父子指针、子框半边长/四分之一中心偏移、点唯一归属、叶点总数、重编号为置换;失败即终止,不降级(`DBEM1/multidomain.cpp:3979-4129`)。

### F-9 可观测性现状

- `missingLocalBlocks` 只 `printf` 到 stdout(`DBEM1/multidomain.cpp:4359`),**未写入任何 CSV / `validation_metrics.txt`**。全仓库检索确认:不存在任何一次运行留下的该数值。
- `PrintPreconditionerSelfCheck`(计算 `M⁻¹Ax` 与 `x` 的偏差)门禁为 `DBEMPolicy::solveDiagnostics && DomainCount()==1`(`DBEM1/multidomain.cpp:4666`),**对多域路径不生效**。
- `BuildBlockDiagonalPreConditioner`(`DBEM1/multidomain.cpp:4367`)已实现但主路径未调用,且按单域同编号假设编写。
- `DBEM_POLICY_VALIDATION_OUTPUT=1` 时输出 `tree_leaf_map.csv`(每点所属叶、叶中心、边长、层数)与 `validation_metrics.txt`。

## S1-J 判断条目(推断,均未实测)

| 编号 | 判断 | 依据 | 若被推翻的后果 |
| --- | --- | --- | --- |
| J-1 | 缺块主因是"波前带直径 ≪ 叶直径":叶内可查块对约 `m²`(`m`=每叶点数),step 0 只有落在大约 `C2·dt` 半径内的点对才有块;叶近似立方体时边长 `L ≈ h·m^{1/3}` | F-5 条件 1 + F-8 分叶 | 创新点"物理参数化分叶"的动机消失 |
| J-2 | 分域只在"空间树把一个叶切成跨域"时制造缺块(叶内同时含域 A、域 B 的点且非界面配对) | F-5 条件 4 + F-6 | 界面/分域章节需改写 |
| J-3 | "界面导致缺块"不存在独立路径;界面只是让列置换 `p(r)` 非平凡,并缩小矩阵规模 | F-6 + F-2 | 界面章节的核心论点失效 |
| J-4 | 列置换 `p()` 不改变"找得到/找不到":`p(r)` 正是该行装配时写块的那个列,故每个方程行在 `p(r)` 处必有块 | F-2 + F-5 条件 4 | 若存在反例,则现状可能已经错配 |
| J-5 | 退化叶/后备分配会放大问题:节点中心重合的两个单元必被判 `degenerateLeaf`,无论 `premaxleafpointnum` 多小都永远同叶,而它们 `Dist=0`,`n≥1` 时 `NeedCal` 返回 0,该叶内非对角块整片缺失 | F-8 + F-5 条件 1 | 需新增"按 `Dist=0` 特判"的分支 |
| J-6 | 若某叶的**对角块**本身缺失,该行 `I_L` 整行为零 ⇒ `inv_mat` 失败 ⇒ 整个求解硬失败(`DBEM1/multidomain.cpp:4334`),不是静默降级 | F-4 + F-5 条件 3 | — |
| J-7 | 更值得优先验证的是构造时机:预条件器用 step 0 的稀疏带状矩阵构造却全程复用,而波前带随步数单调变宽,预条件器所近似的算子与实际要解的算子随步数偏离 | F-7 + F-5 条件 1 | 若偏离可忽略,则该优先级下降 |
| J-8 | `restart=1` ⇒ 每外迭代只做 1 步内迭代并重算真实残差,形式上接近预条件 Richardson/最速下降,度量的是"`M⁻¹` 近似单位阵有多好",而非 GMRES 的谱分布信息 | F-7 | 若语义理解有误,所有迭代次数对比都要重做 |
| J-9 | 由 `NeedCal` 可得:同材料内 `r_ij ≥ C1(n+1)Δt ⇒ A_n(i,j) = 0`(严格零),即 `A_n` 的块图在距离度量下是半径 `B_n ≈ (n+1)β·h_min` 的带;定义 `Θ_n = B_n / L_leaf` 为无量纲控制参数 | F-5 条件 1 | 创新路线 A 的立论基础 |
| J-10 | 若 `diam(L) ≤ B_n`,则被丢弃的跨叶块在 `A_n` 中**严格为零**,`M⁻¹ = blockdiag(A_L⁻¹)` 在锥内**精确**(非近似) | J-9 + F-4 | 创新路线 A/B 的"精确性"论证 |

## S1-Q 需要 ChatGPT 复核的问题

### Q-1(J-1 的定量推理是否成立)

我的推理:叶内可查块对约 `m²`,但有块的仅落在半径 `B = C2·dt` 球内的点对,比例约 `(B/L)³` 量级,`L ≈ h·m^{1/3}`,故缺块比例很大。

请判断:
(a) 这个量级推理是否成立?
(b) 是否可能缺块其实主要由 J-2(跨域叶)或 J-5(退化叶)造成,而非波前带?
(c) 若 `β` 较大或叶较小时 `B` 可能接近甚至超过 `L`,此时 J-1 是否自动失效?

### Q-2(J-3 / J-4 是否有反例)

请重点检查:
(a) 界面行在 A 侧与 B 侧、面力符号 `sign`、列变换 `TransformColumnBlock9` 之后,是否可能出现"该行在 `p(r)` 处确实没有块"的情况(例如变换后被 1e-14 阈值剔除)?
(b) 如果 (a) 不成立,那么"界面不产生缺块、`p(r)` 处必有块"是否可以视为严格结论?

### Q-3(J-7 的处置策略)

用 step 0 的 `A` 构造预条件器并全程复用,在时域波前问题里通常应如何选择?可选:
(1) 每步重建;(2) 固定滞后重建;(3) 用末期矩阵构造;(4) 增量更新。
我倾向 (4),但需要你判断在"矩阵序列单调增长"这一结构下,(4) 的增量更新(Sherman–Morrison–Woodbury / 稀疏低秩更新 + 触发式重建)在文献里是否已有等价做法。

### Q-4(J-8 的确认)

`restart=1` 是否确实等价于每外迭代 1 步内迭代?在这种设置下,"`M⁻¹` 的作用次数"与"块 Jacobi 迭代"之间是什么关系?如果要做预条件器对比,`restart` 应设为什么值才符合惯例?

### Q-5(创新路线评估)

以下哪条最可能被同行认为有真正创新性?各需要补什么?

| 路线 | 内容 |
| --- | --- |
| A | 用因果锥参数化分叶:`Θ_n = B_n/L_leaf` 取代经验旋钮 `premaxleafpointnum`,并证明"锥内精确、锥外可控" |
| B | 锥外截断误差的可证界(远场块数值秩低,复用现有 ACA 设施验证) |
| C | 利用 `{B_n}` 单调性做增量预条件器(总成本 ∝ 新增块数 + 触发式重建) |
| D | 界面 Schur 补粗空间(界面共享未知量是唯一的非局部耦合通道) |

另外请判断:把 `B_n/L` 作为预条件器的设计变量,与已有的"用因果块稀疏加速组装/卷积"的工作(MOT 时域积分方程)之间,差异是否足够构成创新点?

### Q-6(实验设计)

我计划在叶内抽取循环处增加缺块分类计数(判定所需数据在构建期全部现成):

- `missingByWavefront`:对 step 0 重算 `NeedCal` / 用 `Dist` 与 `C2·dt` 比较,落在窗外;
- `missingCrossDomain`:`m_domainByElement[row] != m_domainByElement[col]`;
- `missingThreshold`:以上都不成立者;
- `leafSpansDomains`:该叶域归属取值数 > 1;
- `missingDiagonalBlocks`:`row == col` 却缺。

请评估:(a) 这个分类是否足以判定 J-1/J-2 的相对权重?(b) 还缺什么量?(c) 用 `tree_leaf_map.csv` + 域归属做离线 join 是否可行?

### Q-7(风险出口)

我识别出四个风险,并认为每个"失败"都仍有出版价值:

| 编号 | 风险 | 我的"失败出口" |
| --- | --- | --- |
| R1 | `B_n/L` 本就大,缺块少,改进空间有限 | 转为"当前预条件器足够好"的负结果 + 代价模型 |
| R2 | 远场块数值秩不衰减 | 转为"块 Jacobi 类预条件器在时域 BEM 中有本质天花板,需要 Dirichlet–Neumann 类界面方法" |
| R3 | 增量更新触发重建太频繁 | 退回"固定滞后重建 + 周期估计"的工程版本 |
| R4 | 界面粗空间收益微小 | 用占比数据给出设计准则 |

请判断这些"失败出口"是否真的可发表,以及有没有我没想到的第五个风险。

---

# 第 2 部分:ChatGPT → 代码侧

> 说明:把回答粘贴到对应编号下方。每条回答请标明对应的 `Q` 编号;若回答中引入新问题,请用 `Q-8`、`Q-9`… 继续编号,我会在下一轮回应。

## A-1(对应 Q-1)

<!-- 待粘贴 -->

## A-2(对应 Q-2)

<!-- 待粘贴 -->

## A-3(对应 Q-3)

<!-- 待粘贴 -->

## A-4(对应 Q-4)

<!-- 待粘贴 -->

## A-5(对应 Q-5)

<!-- 待粘贴 -->

## A-6(对应 Q-6)

<!-- 待粘贴 -->

## A-7(对应 Q-7)

<!-- 待粘贴 -->

---

# 第 3 部分:任务与决定

## 当前待办

| 编号 | 任务 | 状态 | 备注 |
| --- | --- | --- | --- |
| T-1 | 确认多域 `gmres` 的 `restart` 语义,评估改为 `restart=30` 的代价 | 待办 | 优先级高于任何算法改动;见 J-8 / Q-4 |
| T-2 | 实现缺块分类计数(`missingByWavefront` 等五项) | 待办 | 见 Q-6 |
| T-3 | 用 `tree_leaf_map.csv` + 域归属做离线 join,统计"跨域叶"比例 | 待办 | 零代码改动,可先做 |
| T-4 | 测量 `nnz(A_n)` 随 `n` 的增长,拟合 `B_n` | 待办 | 验证 J-9 |
| T-5 | 扫 `premaxleafpointnum`,画"迭代次数 vs `Θ_n`"主图 | 待办 | 决定创新路线 A 是否值得做 |
| T-6 | 测量远场块的 ACA 数值秩随距离的衰减 | 待办 | 为创新路线 B 提供实测指数 |

## 已确定的决定

| 日期 | 决定 | 影响 |
| --- | --- | --- |
| 2026-02 | 默认只考虑多域路径;单域仅作兼容入口与对照 | 见 `README.md` 开头工作约定;`DomainCount()==1` 走多域代码路径的算例保留为对照 |

## 变更记录

| 日期 | 变更 |
| --- | --- |
| 2026-02 | 创建本文件;登记 S1 的 F-1~F-9、J-1~J-10、Q-1~Q-7、T-1~T-6 |
