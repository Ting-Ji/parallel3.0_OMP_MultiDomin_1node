# 三维多域时域边界元（TDBEM）并行求解程序

本项目是在原有三维时域边界元程序上发展出的研究代码，用于求解各向同性线弹性体的瞬态动力响应。当前 `master` 分支的主要工作对象是：**常物理单元、多个材料子域、匹配界面耦合、3×3 块稀疏矩阵（CCSR/mCCSR）与 GMRES 迭代求解，以及面向正确性验证的诊断输出**。

程序保留了原有单域 Gauss、GMRES 和 ACA 等路径，但目前新增且持续验证的核心路径是“动力学 + 多域 + GMRES + CCSR”。本文档按当前代码实际行为说明任务内容、数学模型、离散方法、代码组织、输入输出、并行实现、已完成目标和已知限制。

## 1. 项目任务与当前范围

给定边界几何、材料参数、时间步数和随时间变化的位移或面力边界条件，程序完成以下计算：

1. 读取 8 节点二次四边形边界网格和边界条件；
2. 由全局参数计算统一时间步长，并为每个子域建立独立材料参数；
3. 将边界单元按连续编号划分到多个子域，建立显式或自动界面配对；
4. 对外边界未知量和界面共享未知量进行全局自由度编号；
5. 装配当前时刻矩阵、已知边界矩阵和各历史时滞矩阵；
6. 在每个时间步构造卷积历史项，通过预条件 GMRES 求解未知位移或面力；
7. 回填物理状态，检查线性残差和界面连续/平衡条件，并输出 Tecplot 与 CSV 结果。

### 1.1 已完成的主要目标

- 已接入可选多域参数段，并支持各子域独立的 `E`、`v`、`Rou`。
- 已将动态材料参数封装为 `DynaMat`，多域积分时通过线程局部材料上下文选择当前源域材料，避免线程间覆盖全局材料。
- 已建立多域全局自由度映射：外边界按 `BCID` 选择未知量，界面两侧共享位移块，并按法向符号共享面力块。
- 已实现 3×3 块级 CCSR 装配器，支持重复块累加、近零块过滤、有限值检查和确定性合并。
- 已实现多域当前项与历史项的 pthread 并行装配，并保留串行回退和并行/串行一致性检查。
- 已实现多域 GMRES 时间递推、叶块预条件、矩阵结构检查、残差诊断和界面传递审计。
- 已统一主要诊断文件的算例输出目录，并提供单域/多域回归脚本和一维杆解析验证脚本。
- 当前常单元版本已统一为“每个边界单元一个物理未知点”，即 `NodeNum = EleNum`。

### 1.2 尚未宣称完成的目标

- 多域与等价单域在所有网格、时间步和材料组合下的数值等价性；
- 网格与时间步收敛性、复杂三维算例的高精度基准验证；
- 非匹配界面、mortar/投影耦合和每域独立时间网格；
- 自适应时间步、分布式多节点并行和全求解流程的并行化；
- 多域路径中的分量混合边界条件、批量边界文件和 ACA/Gauss 全面接入。

因此，当前代码适合算法开发、回归测试和受控算例研究；用于正式工程结论前，应针对目标问题重新进行网格、时间步、材料和解析/实验基准验证。

## 2. 数学原理

### 2.1 三维线弹性动力学

忽略体力时，子域 `Ω_d` 内的动力平衡方程为

```math
\nabla \cdot \boldsymbol{\sigma}(\mathbf{x},t)
= \rho\,\ddot{\mathbf{u}}(\mathbf{x},t),
\qquad \mathbf{x}\in\Omega_d,
```

其中 `u` 为位移，`σ` 为应力，`ρ` 为密度。各向同性线弹性本构关系为

```math
\boldsymbol{\varepsilon}
=\frac{1}{2}\left(\nabla\mathbf{u}+\nabla\mathbf{u}^{\mathsf T}\right),
\qquad
\boldsymbol{\sigma}
=\lambda\,\mathrm{tr}(\boldsymbol{\varepsilon})\mathbf{I}
+2G\boldsymbol{\varepsilon},
```

```math
G=\frac{E}{2(1+\nu)},
\qquad
\lambda=\frac{2G\nu}{1-2\nu}.
```

边界面力定义为

```math
\mathbf{t}=\boldsymbol{\sigma}\mathbf{n},
```

其中 `n` 是子域外法向。代码中的纵波和横波速度分别为

```math
c_1=\sqrt{\frac{\lambda+2G}{\rho}}
=\sqrt{\frac{2(1-\nu)}{1-2\nu}\frac{G}{\rho}},
\qquad
c_2=\sqrt{\frac{G}{\rho}}.
```

`BuildDynaMat()` 预计算 `c1`、`c2`、`c1*dt`、`c2*dt` 及动态基本解积分所需的倒数和组合系数。

### 2.2 时域边界积分方程

直接时域边界元将域内偏微分方程转换为边界上的时空卷积方程。其示意形式为

```math
C_{ij}(\boldsymbol{\xi})u_j(\boldsymbol{\xi},t)
=\int_0^t\!\int_{\Gamma}
U_{ij}(\mathbf{x},t-\tau;\boldsymbol{\xi})t_j(\mathbf{x},\tau)
\,\mathrm d\Gamma\,\mathrm d\tau
-\int_0^t\!\int_{\Gamma}
T_{ij}(\mathbf{x},t-\tau;\boldsymbol{\xi})u_j(\mathbf{x},\tau)
\,\mathrm d\Gamma\,\mathrm d\tau .
```

`U_ij` 和 `T_ij` 分别是三维弹性动力学位移基本解和面力基本解，`C_ij` 是边界自由项。具体核函数、奇异/近奇异处理和数值积分位于 `dsquareelement.cpp`、`dsquareelementdyna.cpp` 及相关高斯积分模块中。

### 2.3 几何插值与常物理单元

当前分支需要区分两类节点：

- **几何节点**：每个边界单元由 8 节点二次四边形（serendipity）描述，用于曲面映射、雅可比、法向量和积分点坐标；
- **物理未知点**：每个单元只在局部中心 `(ξ,η)=(0,0)` 设置一个物理点，位移和面力在该单元上按常量处理。

因此当前自由度规模为

```math
N_{\text{physical node}}=N_{\text{element}},
\qquad
N_{\text{scalar equation}}=3N_{\text{element}}.
```

每个块行或块列对应一个单元中心的三个笛卡尔分量，矩阵基本存储单元是 `3×3` 子块。参数 `gap` 仍保留在历史单元与积分初始化接口中，但当前物理未知点固定在单元中心。

### 2.4 外边界条件与界面条件

多域 GMRES 路径目前只接受整块边界条件：

- `BCID = 123`：三个位移分量已知，三维面力为未知量；
- `BCID = 456`：三维面力已知，三个位移分量为未知量。

分量混合形式的 `BCID` 在旧代码中存在历史处理，但尚未完整接入当前多域自由度映射。

对两个相邻子域 A、B 的匹配界面，物理条件为位移连续和面力平衡：

```math
\mathbf{u}^{A}=\mathbf{u}^{B},
\qquad
\mathbf{t}^{A}+\mathbf{t}^{B}=\mathbf{0}.
```

代码使用同一个全局位移块表示两侧位移，并令

```math
\mathbf{t}^{B}=s_n\mathbf{t}^{A},
```

其中通常 `s_n = -1`。B 侧局部量在装配和回填时会转换到 A 侧参考坐标系。当前常单元界面每个配对单元只有一个局部物理点，界面映射固定为 `0 ↔ 0`。

### 2.5 时间步长与历史递推

程序先用参数卡中的全局材料和平均单元尺寸计算统一时间步长：

```math
\Delta t=\beta\frac{h_{\mathrm{avg}}}{c_1}.
```

随后所有子域用同一个 `dt` 构造各自的 `DynaMat`。各域的波速可以不同，但当前没有局部时间步或自适应时间步。

离散后，每一时刻求解可概括为

```math
\mathbf{A}\mathbf{x}_n
=\mathbf{K}\mathbf{y}_n
+\sum_{\ell=1}^{\min(n,N_h)}
\left(\mathbf{G}_{\ell}\widehat{\mathbf{t}}_{n,\ell}
-\mathbf{T}_{\ell}\widehat{\mathbf{u}}_{n,\ell}\right),
```

其中 `x_n` 是当前未知边界量，`y_n` 是当前已知边界量，`A` 和 `K` 分别对应代码中的当前未知矩阵和 `KnownM`，历史矩阵对应 `MGS[ℓ]` 和 `MTS[ℓ]`。历史卷积组合按可用状态构造为

```math
\widehat{\mathbf{z}}_{n,\ell}
=\mathbf{z}_{n-\ell+1}+2\mathbf{z}_{n-\ell}+\mathbf{z}_{n-\ell-1}.
```

代码中 `MTS = T2 + T1`，`MGS = U`。在边界类型、网格和固定 `dt` 不变时，当前矩阵只装配一次，各时间步重复使用；历史矩阵按时滞预先装配并在递推中复用。

### 2.6 块稀疏矩阵与 GMRES

`MultiDomainCCSRBuilder` 以 `(blockRow, blockCol)` 为键累加 `3×3` 子块，并在最终构建时生成 CCSR 数据：

- 同一位置的多个贡献进行代数累加；
- 绝对阈值 `1e-14` 过滤近零块，相对阈值 `1e-10` 用于比较和对称性判断；
- 非有限块会被拒绝或在结构诊断中报告；
- 历史 `MGS` 的每个 `3×3` 块满足对称检查时使用六分量存储，否则自动回退为完整九分量 CCSR；
- 基于有序映射和固定线程合并顺序，减少并行装配结果的非确定性。

线性方程采用预条件 GMRES 求解。单域兼容路径使用原有 `GMRESPreConditioner`；多域路径根据物理点空间树构造映射后的叶块预条件器，从全局 CCSR 中提取局部稠密块并求逆。该预条件器改善迭代性能，但不改变离散方程本身。

### 2.7 当前多域预处理器详解

#### 2.7.1 类型与作用方式

当前多域代码实际启用的是 `BuildMappedLeafPreConditioner()`。从数值线性代数角度，它可以归类为：

> 带“方程行—未知量列”置换的、按几何空间树划分的非重叠块 Jacobi（block Jacobi）/局部近似逆预处理器。

它不是 ILU/ILUT，不包含跨叶块消元，也不是以物理子域为单位构造的 Schur 补或粗网格校正。设当前时刻线性系统为

```math
\mathbf{A}\mathbf{x}=\mathbf{b}.
```

GMRES 中先计算真实残差 `r=b-Ax`，再由 `GMRES_M_X()` 计算近似逆作用 `M^{-1}r`；Arnoldi 过程中同样计算 `M^{-1}Av`。因此实际构造的是左预处理系统

```math
\mathbf{M}^{-1}\mathbf{A}\mathbf{x}
=\mathbf{M}^{-1}\mathbf{b}.
```

收敛判断仍使用未预处理的真实残差 `||b-Ax||`，不是只检查预处理残差。

#### 2.7.2 为什么多域需要额外的列映射

全局矩阵的块行按物理点/方程编号排列，但多域未知向量同时包含外边界未知量、界面共享位移和界面共享面力，未知块顺序不一定与方程行相同。代码为每个方程块行 `r` 定义首选未知块

```math
p(r)=\texttt{GetPreferredUnknownBlockForRow}(r).
```

映射规则为：

- 外边界 `BCID=123` 的行优先对应该点未知面力块；
- 外边界 `BCID=456` 的行优先对应该点未知位移块；
- 对界面配对，A 侧行优先对应共享位移块，B 侧行优先对应共享面力块。

`GlobalDofMap::Validate()` 要求 `p(r)` 对全部未知块形成一一映射。若没有这一置换，直接抽取 `A[r,c]` 的同编号对角叶块，会把部分界面方程与错误类型的未知量配对，所得局部块不再代表合理的局部方程系统。

#### 2.7.3 构造过程

预处理器只从当前未知矩阵 `A` 构造一次，之后由所有时间步复用。具体过程为：

1. 使用全部单元中心物理点的坐标建立空间树；
2. 通过 `premaxleafpointnum` 限制每个叶块最多包含的物理点数；
3. `RenumberPointID()` 得到按叶块连续存放的物理点编号 `m_RePID`；
4. 对树序位置 `i` 设置：
   - `m_InputPID[i] = m_RePID[i]`，表示原矩阵方程块行；
   - `m_OutputPID[i] = p(m_RePID[i])`，表示与该行配对的首选未知块列；
5. 对每个叶块 `L`，抽取局部稠密矩阵

```math
\mathbf{A}_L
=\left[\mathbf{A}_{r,p(c)}\right]_{r,c\in L},
```

   其中每个元素是 `3×3` CCSR 子块；全局 CCSR 中不存在的局部块按零块处理，跨叶块耦合被直接舍弃；
6. 对大小为 `3|L| × 3|L|` 的 `A_L` 原地求逆，将结果保存到 `m_PreM[leaf]`；
7. 应用预处理器时，把方程顺序的输入向量按 `m_InputPID` 收集，乘以各叶块逆矩阵，再按 `m_OutputPID` 散射到未知量顺序。

对任意方程空间向量 `v`，叶块内的作用可写为

```math
\left(\mathbf{M}^{-1}\mathbf{v}\right)_{p(c)}
=\sum_{r\in L}\left(\mathbf{A}_L^{-1}\right)_{c,r}\mathbf{v}_r,
\qquad c\in L.
```

各叶块互不重叠，且 `m_OutputPID` 经验证为全局唯一，因此 `GMRES_M_X()` 可以用 OpenMP 并行计算不同叶块，不需要对输出向量做原子累加。

#### 2.7.4 空间树根节点、递归和自检

该树仍是每层最多八个子节点的八叉树，但多域预处理器默认不再强制使用立方体根节点。设全部物理点在方向 `k` 上的包围盒跨度为

```math
\Delta_k=x_k^{\max}-x_k^{\min},\qquad
c_k=\frac{x_k^{\min}+x_k^{\max}}{2}.
```

默认 `cuboid` 模式使用轴对齐长方体，正常非退化方向的边长严格为 `L_k=1.02\Delta_k`。为了使平面、直线、全部重合点以及“大坐标、小跨度”点集仍有有限正尺寸，代码定义

```math
S=\max\left(1,\max_k\Delta_k,\max_{i,k}|x_{i,k}|\right),\qquad
\varepsilon_L=64\varepsilon_{\mathrm{mach}}S,
```

并最终取 `L_k=1.02 max(Delta_k, epsilon_L)`。输入空指针、非正点数、非有限坐标、非法填充因子或非有限几何量都会使预处理器构造立即失败。

每个非叶节点仍同时沿 `x/y/z` 二分。子节点编号沿用位编码：bit 0、1、2 分别控制 `x/y/z`，0 取负半区、1 取正半区。子节点各方向边长是父节点的一半，中心偏移是父边长的四分之一。点按 0 到 7 的顺序检查，因此分割面上的点确定性地进入编号最小的可包含子节点；若舍入误差使点未落入任何子框，则按三个方向以子框边长归一化后的中心距离选最近子节点，并累计 `TreeFallbackAssignmentCount`。

除叶点数阈值外，递归有两项强制保护：节点内全部坐标跨度均低于数值分辨率时生成 `degenerateLeaf`；达到第 64 层时生成 `maxDepthLeaf`。这两类叶节点可以超过 `premaxleafpointnum`，但都会记录原因并输出计数，防止退化点集无限递归。

树建成后、抽取局部矩阵前会自检根框、父子指针、半边长/四分之一中心偏移、点的唯一归属、叶点总数、普通叶点数上限及重编号排列。任何缺点、重复点、非法几何或错误父子关系都会终止构造。`DBEM_VALIDATION_OUTPUT=1` 时还会生成 `tree_leaf_map.csv`，记录每个点所属叶节点及该叶的中心、边长和层数；`validation_metrics.txt` 记录根框、长宽比、叶数、最大层数、退化叶、最大深度叶及后备分配等统计。

运行时可用下列开关做严格对照：

```text
DBEM_MULTIDOMAIN_TREE_ROOT=cuboid|cube
```

默认值为 `cuboid`；`cube` 使用修改前的最长包围盒跨度作为三个方向的统一边长。非法值直接报错。该开关只接入多域 GMRES 叶块预处理器；旧单域 GMRES 和 ACA 仍通过兼容入口建立立方体树。两种模式只允许改变点到叶块的归属，不改变全局矩阵、自由度映射或预处理代数形式。

界面划分不要求“同一界面的所有单元必须进入同一叶节点”。例如含 16 个单元的界面可以按空间局部性进入多个叶节点；强行把整个界面合成一叶会随界面尺寸增长而形成昂贵且更易病态的稠密块。当前几何树也不主动保证每一对界面单元位于同一叶，界面两侧强耦合可能被跨叶舍弃；若实验表明这是主要收敛瓶颈，后续应实现界面配对感知、重叠叶块或界面粗空间，而不是把整个界面无条件捆绑。

#### 2.7.5 参数、代价与预期效果

`premaxleafpointnum` 是当前预处理器唯一的主要强度参数。若单个叶块包含 `m` 个物理点，则局部矩阵阶数为 `3m`：

- 构造存储量约为 `O((3m)^2)`；
- 稠密求逆代价约为 `O((3m)^3)`；
- 每次应用代价约为 `O((3m)^2)`。

较小叶块构造和应用便宜，但只保留很局部的耦合，GMRES 迭代次数可能增加；较大叶块能保留更多近场耦合，但求逆成本、内存和病态风险迅速上升。最佳值依赖网格、界面位置、材料波阻抗差和边界条件，不能只凭单域经验固定。

它对当前问题有效的原因是：边界积分矩阵虽有全局耦合，但几何上相近的单元、同一界面附近的位移/面力通常具有较强耦合。叶块逆矩阵先消除这些局部强耦合，再由 GMRES 处理叶块之间被舍弃的长程耦合。

#### 2.7.6 当前实现的限制

- 空间树按几何位置划叶，不认识材料子域、界面配对和波阻抗；界面两侧点可能被分到不同叶块。
- 叶块不重叠，没有相邻叶块信息、界面粗空间或全局低频修正，材料反差增大时收敛可能明显恶化。
- 缺失的叶内 CCSR 块直接置零，只打印 `missingLocalBlocks` 数量，没有据此自适应扩大或重组叶块。
- `inv_mat()` 仍以固定绝对阈值判断奇异；当前构造会检查其返回值和逆矩阵有限性并在失败时终止，但尚无条件数估计、尺度自适应阈值、正则化或稳健回退。
- 当前没有运行时预处理器选择。代码中虽定义了 `BuildBlockDiagonalPreConditioner()`，多域求解主路径并未调用它；它只能视为尚未接入的候选基线。
- 只有单域分支调用了 `PrintPreconditionerSelfCheck()`；多域分支当前没有系统检查 `M^{-1}Ax`、叶块逆误差或预处理前后谱/残差变化。

因此，日志中的 `leaves`、`localBlocks` 和 `missingLocalBlocks` 只描述构造规模，不能单独说明预处理质量。评估预处理器至少还应比较 GMRES 迭代数、总求解时间、真实残差、叶块求逆失败数，以及随子域数和材料反差变化的鲁棒性。

## 3. 当前求解流程

程序的主执行过程如下：

1. 从当前工作目录读取 `BEM_DATACARD.DAT`；
2. 根据第一行算例名读取 `input/<case>.model` 和 `input/<case>.bd`，或读取对应二进制文件；
3. 计算平均单元尺寸、全局波速和统一 `dt`；
4. 读取可选多域段，连续划分子域并构造各域 `DynaMat`；
5. 初始化 8 节点几何单元，并在每个单元中心生成一个物理未知点；
6. 标记外边界/界面单元，解析界面配对和局部坐标变换；
7. 建立外边界、界面、已知量及历史量的全局块自由度映射；
8. 并行装配 `A`、`KnownM`、`MTS[ℓ]`、`MGS[ℓ]`，检查矩阵结构并构造预条件器；
9. 对 `n = 1...NStep` 形成当前已知项和历史卷积项，调用 GMRES；
10. 将解向量散射回位移/面力状态，执行界面审计和线性残差检查；
11. 输出状态、指标、诊断 CSV 和 Tecplot 文件。

主程序当前将问题类型固定为动力学；GMRES 分支还会将动态方法标志强制设为 `2`。动力学影响长度 `MaxLe` 当前在主程序中固定为 `1000`。

## 4. 代码结构与实现职责

| 文件/目录 | 主要职责 |
| --- | --- |
| `DBEM1/DBEM.cpp` | 主程序；参数卡、模型和边界条件读取；文本/二进制输入；时间步计算；求解器分派；公共验证状态输出 |
| `DBEM1/dsquareelement.h/.cpp` | 8 节点曲面单元、几何映射、动态材料参数、基本解积分、线程局部 `DynaMat` 入口 |
| `DBEM1/dsquareelementdyna.cpp` | 动态单元矩阵和时域核积分相关实现 |
| `DBEM1/multidomain.h/.cpp` | 多域输入、子域/界面模型、自由度映射、块稀疏 builder、pthread 装配、多域 GMRES 和诊断 |
| `DBEM1/mat_vec.h/.cpp` | 向量、稠密矩阵、CCSR/SymCCSR 数据结构及基础运算 |
| `DBEM1/pre_gmres.h/.cpp` | 原有单域 GMRES、预条件和单域 RHS/残差诊断 |
| `DBEM1/PthreadCCSR.h/.cpp` | 原有 CCSR pthread 计算辅助代码 |
| `DBEM1/tree.h/.cpp` | 空间树、叶节点划分和预条件器局部块组织 |
| `DBEM1/Gauss_solver.*`、`DBEM1/aca.*` | 保留的直接求解与 ACA 路径 |
| `DBEM1/output_path.h/.cpp` | 生成并管理 `output/<case>_<NStep>/` 算例目录 |
| `DBEM1/write_tec.*`、`DBEM1/postpro.*` | Tecplot 数据和后处理输出 |
| `convert_inputs_to_binary.py` | 将 `.model`、`.bd` 文本输入转换为版本化二进制输入 |
| `tools/generate_multidomain_rod.py` | 生成多域杆测试输入 |
| `tools/run_multidomain_validation.py` | 单域/多域回归、界面、残差和 RHS 对照测试 |
| `tools/run_rod120_100_analytic_validation.py` | 材料污染检查和一维杆解析解验证 |

## 5. 输入文件

程序通常从 `DBEM1` 目录运行，目录关系为：

```text
DBEM1/
├── BEM_DATACARD.DAT
├── gaussposition
├── gaussweight
├── input/
│   ├── <case>.model
│   └── <case>.bd
└── output/
```

### 5.1 参数卡公共部分

`BEM_DATACARD.DAT` 的公共字段按以下顺序逐项读取：

| 顺序 | 字段 | 含义 |
| ---: | --- | --- |
| 1 | `caseName` | 输入文件基名，不含 `.model`/`.bd` |
| 2 | `thread_num` | OpenMP/pthread 线程数 |
| 3 | `FlagDyna` | 动态离散方法标志；GMRES 路径当前会设为 `2` |
| 4 | `Beta` | 时间步尺度系数 `β` |
| 5 | `NStep` | 最大时间步编号，状态范围为 `0...NStep` |
| 6 | `gap` | 历史不连续单元/积分初始化参数 |
| 7 | `E` | 全局弹性模量，同时用于计算全局 `dt` |
| 8 | `v` | 全局泊松比 |
| 9 | `Rou` | 全局密度 |
| 10 | `amplitude` | Tecplot 变形显示放大系数 |
| 11 | `solver` | `1=Gauss`，`2=GMRES`，`3=ACA` |

GMRES（`solver=2`）继续读取：

```text
iterations
error
premaxleafpointnum
InputBinaryMode
BatchBoundaryMode
```

- `iterations`：GMRES 最大迭代控制参数；
- `error`：GMRES 收敛容差；
- `premaxleafpointnum`：预条件器每个空间树叶块的最大物理点数；
- `InputBinaryMode`：`0` 读取文本，`1` 强制读取二进制，`2` 存在二进制时优先读取；
- `BatchBoundaryMode`：`0` 单一边界文件，`1` 显式文件列表，`2` 按目录、前缀和编号范围生成文件名。

当 `BatchBoundaryMode=1` 时，后续格式为：

```text
BatchBoundaryCaseCount
boundaryFile1
boundaryFile2
...
```

当 `BatchBoundaryMode=2` 时，后续格式为：

```text
BatchBoundaryCaseCount
BatchBoundaryDir
BatchBoundaryPrefix
BatchBoundaryStartIndex
BatchBoundaryEndIndex
BatchBoundaryExt
```

模式 2 最终按 `dir + prefix + index + ext` 生成文件名；编号范围数量应与 `BatchBoundaryCaseCount` 一致。

多域 GMRES 当前不支持 `BatchBoundaryMode=1/2`，启用后会主动终止。Gauss 分支能够忽略一段兼容旧参数卡的 GMRES 参数，但新增多域求解器目前只在 `solver=2` 分支被主程序调用。

ACA（`solver=3`）参数依次为：

```text
iterations error maxleafelenum premaxleafpointnum
errorvalue errorstop Flag_ReadWrite
```

### 5.2 多域扩展段

不追加任何字段时，`ReadOptionalMultiDomainInput()` 返回未启用，程序沿旧单域路径运行。要启用多域，必须在求解器参数之后追加：

```text
DomainCount
domainId E v Rou
...                  # 必须恰好 DomainCount 行
[InterfaceCount]
[domainA eleA domainB eleB normalSign]
...
```

一个两域片段示例为：

```text
2
1 20000 0.25 8
2 30000 0.25 8
1
1 6 2 7 -1
```

该示例要求第 6 个单元属于域 1、第 7 个单元属于域 2，并将二者配为界面。

多域输入必须遵守以下规则：

- `DomainCount > 0`，且实际使用时不能多于单元数；
- 材料行不是可选回退项：解析器要求恰好读取 `DomainCount` 行，建议每个域 ID 恰好出现一次；
- 域 ID 和显式界面中的单元 ID 都支持 0 基或 1 基；只要相应段中出现 `0`，该段就按 0 基解释，因此不要混用编号体系；
- 单元按输入顺序连续、近似均匀地分给各域，余数优先分给编号较小的域；当前不能在参数卡中给出任意的“单元—域”列表；
- `InterfaceCount` 若存在，后面必须有相同数量的界面行；写入 `0` 表示显式声明无界面；
- `InterfaceCount` 若完全省略，程序只自动连接每对相邻域的“前一域最后一个单元”和“后一域第一个单元”。该规则仅适合极简测试，真实界面应显式逐单元列出；
- 同一个界面单元不能在多个 `InterfacePair` 中重复使用；界面单元必须确实属于所写子域；
- `normalSign=0` 会被转换为 `-1`，物理上通常应直接填写 `-1`。

### 5.3 `.model` 几何文件

文本格式为：

```text
PointNum
x y z                         # PointNum 行
EleNum
n1 n2 n3 n4 n5 n6 n7 n8      # EleNum 行，节点号为 1 基
InfEleNum
i1 i2 i3 i4                  # InfEleNum 行，4 个关联编号为 1 基
```

普通边界单元必须给出 8 个二次四边形几何节点。`InfEleNum` 可以为 `0`。程序读入文本后会把普通单元节点号和无限单元的 4 个关联编号转换为内部 0 基编号。

### 5.4 `.bd` 边界文件

文本格式为：

```text
BCID[0]
...
BCID[EleNum-1]
value_x value_y value_z       # step=0, element=0
...
value_x value_y value_z       # step=NStep, element=EleNum-1
```

即先给出 `EleNum` 个边界类型，再给出 `(NStep+1)×EleNum` 组三分量值。输入值按单元局部坐标读取，求解过程中会在局部坐标与全局笛卡尔坐标之间转换。

### 5.5 二进制输入

文本文件可转换为 `<name>.model.bin` 和 `<name>.bd.bin`：

```powershell
python convert_inputs_to_binary.py --model DBEM1\input\case.model
python convert_inputs_to_binary.py --bd DBEM1\input\case.bd --ele-num 1000 --nstep 100
```

二进制文件含 magic、版本、类型和端序标志。`InputBinaryMode=1` 时文件缺失或头信息不匹配会终止；`InputBinaryMode=2` 时仅在对应二进制文件存在时使用它。

## 6. 构建与运行

### 6.1 构建环境

工程是 Visual Studio C++ x64 项目，项目文件当前指定 `PlatformToolset=v142`，Release 配置启用 OpenMP，并依赖 pthread 兼容头/库。

`DBEM1/DBEM1.vcxproj` 仍包含原开发环境的绝对头文件和库目录。换机器构建前，应将 `AdditionalIncludeDirectories`、`AdditionalLibraryDirectories` 调整为本机或仓库中的有效路径，并确认已安装兼容的 MSVC v142 工具集和 x64 库。

Visual Studio 2022 的典型构建命令为：

```powershell
& 'C:\Program Files\Microsoft Visual Studio\2022\Community\MSBuild\Current\Bin\MSBuild.exe' `
  DBEM1\DBEM1.vcxproj /p:Configuration=Release /p:Platform=x64 /m:1 /v:minimal
```

默认可执行文件位置为：

```text
DBEM1\x64\Release\DBEM1.exe
```

代码中仍有部分历史编码、未使用变量和格式类型 warning；构建时应首先保证无编译和链接错误，再逐步清理 warning。

### 6.2 运行

从 `DBEM1` 目录启动，确保参数卡、输入目录、高斯点和高斯权重文件可访问：

```powershell
cd DBEM1
.\x64\Release\DBEM1.exe
```

该程序沿用历史返回码约定：正常完成时 `main()` 返回 `1`。自动化脚本因此将 `0` 和 `1` 都视为可接受的进程返回码，但仍需同时检查求解器 flag、日志和验证指标。

## 7. 并行实现

当前并行粒度包括：

- 使用 OpenMP 并行初始化单元几何和单元中心物理点；
- 多域当前项和历史项按子域遍历，并把源单元区间分配给 pthread worker；
- 每个 worker 使用独立 `MultiDomainCCSRBuilder`，不共享写同一个稀疏容器；
- 线程结束后按固定 worker 顺序合并子块，重复位置统一累加；
- 每次积分通过 `ScopedDynaMat` 设置当前线程的源域材料，作用域退出后恢复原材料指针。

这不是“全流程并行”：全局自由度构建、部分界面当前项处理、预条件器构造、GMRES 主迭代和部分后处理仍主要为串行或沿用旧实现。

多域装配相关环境变量如下：

| 环境变量 | 默认值 | 作用 |
| --- | --- | --- |
| `DBEM_MULTIDOMAIN_GMRES_PTHREAD` | 开启 | 设为 `0/false/off/no` 时强制使用串行装配 |
| `DBEM_MULTIDOMAIN_GMRES_PTHREAD_CHECK` | 关闭 | 设为真时同时计算串行结果并比较；超出混合容差会采用串行结果 |
| `DBEM_MULTIDOMAIN_TREE_ROOT` | `cuboid` | 多域叶块预处理器树根；可设为 `cube` 做旧立方体行为对照，其他值直接报错 |

pthread 创建或连接失败时，装配函数会回退到串行实现。并行检查主要验证装配等价性，不等价于最终物理解的完整正确性证明。

## 8. 输出与诊断

按当前目录设计，从 `DBEM1` 运行时，主要新输出位于：

```text
DBEM1/output/<caseName>_<NStep>/
```

| 输出 | 内容 |
| --- | --- |
| `validation_state.csv` | 每个时间步、单元中心的坐标、域、表面类型、`u/t` 分量 |
| `validation_metrics.txt` | `dt`、域材料、GMRES flag/迭代、矩阵块数、耗时、结构检查和界面自由度不变量 |
| `interface_transfer_audit.csv` | 界面两侧位移跳量和面力平衡残差 |
| `rhs_breakdown.csv` | 当前已知项、`G` 历史项、`T` 历史项、总 RHS、`Ax`、残差和解 |
| `residual_probe.csv` | 各步最大绝对残差、L2 残差、RHS 范数和相对残差 |
| `TecValueFile_1node*` | 单元中心结果的 Tecplot 输出 |

少量旧路径仍通过主程序早期构造的父目录 `output/` 写入，例如部分 `log.txt`、`dt.txt` 和历史性能文件；分析结果时应优先使用算例目录中的结构化验证文件，并留意旧输出位置。

诊断开关如下：

| 环境变量 | 默认值 | 作用 |
| --- | --- | --- |
| `DBEM_VALIDATION_OUTPUT` | 开启 | 设为 `0/false/off/no` 可关闭结构化验证输出 |
| `DBEM_RHS_DIAGNOSTIC` | 关闭 | 设为 `1/true` 启用或保留详细 RHS/残差诊断 |
| `DBEM_PARALLEL_KNOWN_RHS_CHECK` | 关闭 | 用于旧单域 GMRES 已知 RHS 并行路径的一致性检查 |

矩阵求解前还会检查零块行、零块列、非有限块、对角块和可逆对角块；界面自由度映射会检查位移/面力共享关系和全局未知块一一对应关系。任何结构检查通过都只说明代数结构满足预期，不能替代解析解或收敛性验证。

## 9. 验证方法与当前证据

### 9.1 多域回归测试

快速测试包含 1 域旧/新路径的短时间步对照和 2 域同材料测试：

```powershell
python tools\run_multidomain_validation.py `
  --exe DBEM1\x64\Release\DBEM1.exe
```

可选模式：

```powershell
# 增加 10 域测试
python tools\run_multidomain_validation.py --ten --exe DBEM1\x64\Release\DBEM1.exe

# 增加 10 域、100 域和分段材料测试
python tools\run_multidomain_validation.py --full --exe DBEM1\x64\Release\DBEM1.exe

# 只诊断单域门槛和两域 RHS 差异
python tools\run_multidomain_validation.py --diagnose-two-domain --exe DBEM1\x64\Release\DBEM1.exe

# 分别运行长方体与立方体树，比较 validation_metrics.txt 中的树统计和 GMRES 指标
python tools\run_multidomain_validation.py --tree-root-mode cuboid --exe DBEM1\x64\Release\DBEM1.exe
python tools\run_multidomain_validation.py --tree-root-mode cube --exe DBEM1\x64\Release\DBEM1.exe
```

脚本检查进程状态、GMRES flag、单域状态差异、界面连续/平衡、线性残差、矩阵结构和 RHS 分解。更改矩阵装配、自由度映射、材料上下文或时间递推后，应至少重新运行快速测试。

路径兼容注意：当前主程序把结构化诊断写入 `output/<caseName>_<NStep>/`，但当前版本的 `run_multidomain_validation.py` 仍从测试目录下的 `output/validation_metrics.txt`、`output/validation_state.csv` 等旧位置读取。使用当前源码重新构建可执行文件后，应先让脚本同步查找算例子目录，或将对应结果整理到脚本预期位置；否则脚本会因“缺少输出文件”失败，而不是因为数值求解本身失败。解析验证脚本已经通过 `case_output_dir()` 兼容两种目录结构。

### 9.2 一维杆解析验证

解析验证脚本会故意修改参数卡的全局材料，同时保持多域材料行不变，以检查积分是否真正使用域材料，然后将受载端轴向位移与一维杆解析解比较：

```powershell
python tools\run_rod120_100_analytic_validation.py `
  --repo . `
  --exe DBEM1\x64\Release\DBEM1.exe `
  --in-place `
  --timeout 900
```

仓库当前的 `DBEM1/analytic_validation.txt` 记录：

```text
material_metrics_ok=1
max_abs=8.413148e-02
max_rel=9.976766e-02
passed=1
```

这表明该次运行通过了脚本当前的材料生效检查和宽松解析阈值；最坏相对误差约为 `9.98%`，不能据此宣称已达到高精度。后续仍应补充网格加密、`Beta/dt` 收敛、多时间区间和不同材料比的系统结果。

### 9.3 建议的验证层次

1. **装配级**：pthread 与 serial 的块数、块值和 RHS 一致；
2. **代数级**：无零行/零列/非有限块，GMRES flag 正常，`||Ax-b||` 满足容差；
3. **接口级**：`uA-uB` 与 `tA+tB` 接近零；
4. **回归级**：单域新路径复现旧路径，多域均质模型与等价单域结果一致；
5. **物理级**：与解析解、半解析解、其他成熟程序或实验数据比较；
6. **收敛级**：验证空间网格、时间步和迭代容差收敛。

界面误差精确为零可能来自共享自由度的代数构造，只能证明约束按预期施加，不能单独证明边界积分核、时间卷积或整体物理解正确。

## 10. 当前限制与注意事项

- 当前是常物理单元：8 节点只用于几何插值，每个单元只有一个三分量物理未知点；不应按 8 个物理节点估计自由度。
- 多域主程序只分派到 GMRES/CCSR 路径；`DynaGaussSolverMultiDomainDense` 虽有实现入口，但当前 `main()` 的 Gauss 分支未调用它。
- 多域外边界只支持 `BCID=123/456`，尚不支持三个分量分别混合给定位移和面力。
- 子域只能按单元输入顺序连续均分，不能直接表示任意拓扑分区。
- 自动界面生成每对相邻域只配一个单元；多单元界面必须显式给出全部配对。
- 界面要求一一匹配的常单元对，同一单元不能重复出现在多个界面中；没有非匹配网格投影。
- 全部子域共享由全局材料和平均单元尺寸计算的固定 `dt`，没有局部/自适应步长和异步卷积。
- `MaxLe=1000` 仍为主程序硬编码值，不是根据当前几何自适应计算。
- 多域历史矩阵按所有有效时滞预装配并驻留内存；大规模、长时间步计算需要继续评估内存和装配复杂度。
- 当前多域预条件器固定为映射叶块方案，尚无运行时选择和系统的条件数/强缩放研究。
- Visual Studio 工程包含历史绝对路径，尚未形成开箱即用的跨机器构建配置。
- 正常返回码为 `1` 是历史约定，外部调度系统若只接受 `0` 需单独适配。
- 当前验证证明的是若干受控算例与代数不变量，不能外推为任意模型下的工程精度保证。

## 11. 后续开发优先级

建议按以下顺序继续推进：

1. 固化可重复的 Release 构建环境，移除工程绝对路径并清理关键编译 warning；
2. 将快速回归纳入每次装配、自由度和时间递推修改后的必跑流程；
3. 解决并量化均质多域与等价单域的全场差异，补齐网格和时间步收敛曲线；
4. 将任意“单元—域”映射和显式界面拓扑从参数卡或独立文件读入；
5. 扩展分量混合边界条件、多域批量边界输入及更多求解器分支；
6. 研究更稳健的块预条件、历史矩阵压缩和 GMRES/后处理并行；
7. 在数值基线稳定后，再发展非匹配界面、局部时间步和分布式多节点版本。

---

本文档描述当前 `master` 分支的实际实现。代码行为发生变化时，应同步更新输入格式、数学离散、验证结果和限制列表，尤其不要把“能够运行”“代数约束为零”和“物理精度已经验证”混为一谈。
