# 多域 TDBEM 的 GMRES-mCCSR 实施说明

本文档说明本工程中“多域 + GMRES + mCCSR”版本的完整目标流程，以及当前已经实现的工作范围。当前代码保留原有单域求解主线，并以独立多域入口逐步接入，避免破坏已有单域结果。

阶段性串行实现和验证结论见：[docs/多域串行实现说明.md](docs/多域串行实现说明.md)。

## 目标完整求解过程

目标版本把任意多个共形子域的 TDBEM 问题转化为统一的全局块稀疏线性系统，然后继续使用现有 GMRES 与 mCCSR/CCSR 存储求解。

完整流程如下：

1. 读取参数卡、模型文件和边界文件。
   - 旧输入文件不包含多域段时，仍按单域流程运行。
   - 多域输入通过可选 `DomainCount` 段启用。
   - 第一版多域假定长杆按单元顺序划分子域，复杂网格自动识别暂不做。

2. 构造几何单元和物理节点。
   - 当前 `DSquareElement` 实际为常值/单物理节点边界单元。
   - 当前实现中 `NodeNum = EleNum`，每个单元中心对应一个 3 自由度块。

3. 构造多域拓扑。
   - 每个单元被分配到一个 `Domain`。
   - 相邻子域之间生成或读取 `InterfacePair`。
   - 界面要求共形，一对界面单元一一对应。

4. 构造每个子域的材料参数。
   - 每个子域可以拥有独立的 `E / v / Rou`。
   - 内部统一转换成 `DynaMat`，包含 `C1`、`C2`、`C1Dt`、`C2Dt` 等动态积分所需系数。

5. 建立全局自由度映射。
   - 外边界沿用现有 `BCID`。
   - 内界面使用共享界面位移和单侧界面面力：
     - `u_A = u_I`
     - `u_B = u_I`
     - `t_A = t_I`
     - `t_B = -t_I`
   - 第一版按 3 自由度块编号，不拆分单个方向。

6. 多域边界积分方程装配。
   - 每个子域只在本域内部做源点与边界单元积分。
   - 不同子域之间不直接使用基本解互相积分。
   - 子域之间只通过界面自由度共享和符号关系耦合。

7. 构造全局 mCCSR/CCSR 矩阵。
   - 使用块稀疏装配器按 `(rowBlock, colBlock)` 聚合 3x3 子块。
   - 重复块必须累加，不能直接使用旧 `WmatrixCCSR::assign` 的追加逻辑。
   - 装配完成后压缩成现有 `CCSRMat`，供 `SparseMul` 和 `gmres` 使用。

8. GMRES 求解。
   - 复用现有 `gmres(..., CCSRMat& A, ...)`。
   - RHS 来自当前已知边界项和历史卷积项。
   - 求解后反写每个子域的位移和面力状态。

9. 时域递推。
   - 每个时间步组装 RHS。
   - 执行 GMRES。
   - 反写子域变量。
   - 更新历史卷积状态。
   - 输出界面连续性误差和性能数据。

## 当前已经实现的工作

当前实现已经推进到多域动态递推验证阶段：数据结构、材料局部化基础、全局块编号、块稀疏装配器、可选多域输入读取、`DomainCount = 1` 矩阵等价自检，以及 `NStep >= 1` 的多域 GMRES 历史卷积入口。

已实现内容包括：

- 新增 `DBEM1/multidomain.h` 和 `DBEM1/multidomain.cpp`。
- 新增核心结构：
  - `Domain`
  - `InterfacePair`
  - `MultiDomainModel`
  - `DomainMaterialContext`
  - `GlobalDofMap`
  - `MultiDomainCCSRBuilder`
  - `MultiDomainState`
- 给 `DSquareElement` 增加：
  - `DomainID`
  - `LocalEleID`
  - `SurfaceType`
- 新增 `BuildDynaMat(...)`，并让原有 `DSquareElement::StaticInit(...)` 调用它，保证单域材料初始化逻辑不重复。
- 新增可选多域参数读取：
  - 老参数卡没有 `DomainCount` 时，继续走原有单域路径。
  - 参数卡追加 `DomainCount` 后，启用多域模型构建。
- 新增长杆顺序均分子域能力：
  - `EleNum` 按 `DomainCount` 均分。
  - 余数单元分配给前面的子域。
  - 未显式给出界面时，自动生成相邻子域界面。
- 新增 `GlobalDofMap`：
  - 外边界按 `BCID` 决定未知量。
  - 内界面共享一套位移和一套面力。
  - 打印未知量分类统计。
- 新增 `MultiDomainCCSRBuilder`：
  - 使用 `std::map<pair<long,long>, array<double,9>>` 聚合块。
  - 支持重复 3x3 块累加。
  - 可输出 `CCSRMat` 和 dense 调试矩阵。
  - 可输出 `SymCCSRMat`，用于沿用旧 `MGS` 对称块乘法路径。
- 新增 `DynaGMRESSolverMultiDomainCCSR(...)`：
  - 装配多域 `Step0` 未知矩阵与已知矩阵。
  - 装配 `MTS[step] = T2 + T1` 和 `MGS[step] = U` 历史矩阵。
  - 使用 `MultiDomainState` 保存 `globalU/globalT/globalUnknown/globalKnown`。
  - 按旧单域公式构造 RHS：当前已知边界项 + 历史 `G*T` - 历史 `T*U`。
  - 使用现有 `gmres(..., CCSRMat& A, ...)` 逐步求解。
  - 求解后反写边界位移和面力，并输出界面误差。
- 更新 Visual Studio 工程：
  - `DBEM1/DBEM1.vcxproj`
  - `DBEM1/DBEM1.vcxproj.filters`
  - 同步更新了旧工程文件 `DBEM1/DBEM.vcxproj` 和 filters，便于查看源码。

## 当前限制

当前实现仍保持保守，已经接入多步动态递推，但需要继续按算例做数值回归和多域验证。

限制如下：

- 多域路径已移除 `NStep = 1` 限制，但仍应先从 `NStep = 2`、`NStep = 5` 小步验证。
- `DomainCount = 1` 时会先运行矩阵-向量等价自检：
  - 检查 `MTS0`、`MTS[0]`、`MTS[1]`、`MGS[1]`。
  - 不通过时会停止进入多域 GMRES 求解。
- 多域路径当前只支持整块边界条件：
  - `BCID = 123`：给定位移，未知面力。
  - `BCID = 456`：给定面力，未知位移。
  - `126 / 153 / 156 / 423 / 426 / 453` 这类逐方向混合边界条件暂不支持。
- 多域路径暂不支持 GMRES batch boundary mode。
- 多域路径暂不使用旧的几何树预处理器。
  - 当前先使用恒等块预处理器，保证编号语义清楚。
  - 后续再接块对角或多域专用预处理器。
- 当前界面条件依赖共形单元一一配对。
- 第一版多域按当前常值单元实现，不支持 8 节点物理自由度版本。

## 多域输入格式

旧参数卡保持不变。要启用多域，在原 GMRES 参数后追加多域段。

基本格式：

```txt
DomainCount
domainId E v Rou
...
InterfaceCount
domainA eleA domainB eleB normalSign
...
```

说明：

- `DomainCount` 不存在时，按单域运行。
- `domainId` 支持从 0 或 1 开始编号，但同一段内应保持一致。
- `eleA / eleB` 支持从 0 或 1 开始编号，但同一段内应保持一致。
- `InterfaceCount` 可省略；省略时自动按相邻子域生成界面。
- `normalSign` 建议写 `-1`，表示 B 侧面力相对 A 侧取反。

示例：2 子域同材料，界面自动生成。

```txt
2
1 20000.0 0.0 8.0
2 20000.0 0.0 8.0
```

示例：2 子域显式界面。

```txt
2
1 20000.0 0.0 8.0
2 20000.0 0.0 8.0
1
1 50 2 51 -1
```

## 推荐验证顺序

后续开发和验证建议按以下顺序推进：

1. 单域旧路径回归。
2. 多域 `DomainCount = 1`，先通过矩阵-向量等价自检。
3. 多域 `DomainCount = 1`，`NStep = 2` 和 `NStep = 5`，与旧单域路径结果对比。
4. `rod_2domain`，相同材料，`NStep = 2/5`。
5. `rod_10domain`，相同材料，检查界面编号、历史矩阵非零块和界面误差。
6. `rod_100domain`，相同材料，检查任意多子域拓扑和 GMRES 收敛。
7. 分段材料长杆，检查每域 `C1/C2/C1Dt/C2Dt` 输出和界面误差。

## 后续实施重点

下一阶段建议优先完成：

1. 准备 `rod_2domain / rod_10domain / rod_100domain` 输入算例并跑数值回归。
2. 对比 `DomainCount = 1` 新旧路径在 `NStep = 2/5` 的位移、面力、GMRES 迭代步和残差。
3. 细化界面误差输出：

```txt
max |u_left - u_right|
max |t_left + t_right|
```

4. 多域块对角预处理器。
5. 整理性能表：矩阵装配时间、GMRES 时间、非零块数、平均迭代步。
