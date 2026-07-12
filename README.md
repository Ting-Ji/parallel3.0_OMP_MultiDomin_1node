# 多域 TDBEM / GMRES-mCCSR 项目说明

本项目是在原有三维时域边界元程序基础上，接入多域材料、多域界面耦合、GMRES 块稀疏求解和验证诊断的一版研究代码。当前重点是多域 GMRES 并行求解路径；旧的单域、Gauss、ACA 路径尽量保持兼容。

阶段性串行实现说明见 [docs/多域串行实现说明.md](docs/多域串行实现说明.md)。

## 当前实现状态

已实现并进入验证范围的内容：

- 多域输入读取：参数卡可选追加 `DomainCount`、多域材料行和显式界面行。
- 多域模型构建：按单元顺序划分 `Domain`，支持显式界面，也支持相邻子域自动生成界面。
- 多域材料参数：每个域使用参数卡多域材料行的 `E / v / Rou` 构造 `DynaMat`。
- 线程局部材料上下文：`DSquareElement::ActiveDynaMat()` 在多域装配时读取当前 worker 域材料，避免并行线程互相覆盖全局静态材料。
- 多域自由度映射：外边界按 `BCID` 确定未知量，界面共享位移和反号面力。
- 多域 mCCSR/CCSR 装配：用块稀疏 builder 聚合 3x3 子块，重复块累加。
- 多域 GMRES 动态递推：装配 Step0/current 矩阵、历史 `MTS/MGS` 矩阵，并按时间步构造 RHS 求解。
- pthread 装配路径：Step0/current 和 history 矩阵支持按域 pthread 并行装配。
- pthread/serial check：可用环境变量开启并行装配和串行 fallback 的一致性检查。
- 混合容差策略：近零块、矩阵差异和 3x3 对称判断使用绝对下限 + 相对尺度，降低材料尺度变化带来的误判。
- 验证输出：支持 `validation_state.csv`、`validation_metrics.txt`、界面传递审计和 RHS 分解诊断。
- `rod120_100` 解析验证脚本：污染全局静态材料行，保持多域材料行不变，用一维杆解析解检查多域材料是否真正生效。

尚未实现的内容：

- 真正的自适应 `dt`。
- 非匹配时间步的多域接口 mortar/投影耦合。
- 每域独立时间网格下的历史卷积重构。
- 混合方向边界条件在多域 GMRES 路径中的完整支持。

当前程序仍使用全局 `Beta` 先计算一个固定 `dt`，再用这个 `dt` 为各域材料构造 `DynaMat`。多域材料已经局部化，但时间步长控制仍不是自适应的。

## 主要源码

- [DBEM1/DBEM.cpp](DBEM1/DBEM.cpp)：主程序、参数卡读取、模型/边界读取、求解器分派、验证输出。
- [DBEM1/dsquareelement.h](DBEM1/dsquareelement.h)、[DBEM1/dsquareelement.cpp](DBEM1/dsquareelement.cpp)：边界单元、动态基本解、`DynaMat` 构造、线程局部材料入口。
- [DBEM1/multidomain.h](DBEM1/multidomain.h)、[DBEM1/multidomain.cpp](DBEM1/multidomain.cpp)：多域输入、拓扑、自由度映射、块稀疏装配、多域 GMRES 求解和诊断。
- [DBEM1/output_path.h](DBEM1/output_path.h)、[DBEM1/output_path.cpp](DBEM1/output_path.cpp)：按算例组织输出目录。
- [tools/run_rod120_100_analytic_validation.py](tools/run_rod120_100_analytic_validation.py)：`rod120_100` 多域材料污染测试和解析解验证。

## 构建

推荐使用 Visual Studio 2022 的 MSBuild：

```powershell
& 'C:\Program Files\Microsoft Visual Studio\2022\Community\MSBuild\Current\Bin\MSBuild.exe' DBEM1\DBEM1.vcxproj /p:Configuration=Release /p:Platform=x64 /m:1 /v:minimal
```

输出程序位于：

```text
DBEM1\x64\Release\DBEM1.exe
```

项目中存在若干历史 warning，例如编码、未使用变量、`printf` 类型不匹配等；当前构建验证以无编译错误为准。

## 运行

程序运行目录通常为 `DBEM1`，需要当前目录下存在 `BEM_DATACARD.DAT`，并能从 `input/` 找到参数卡中指定的 `.model` 和 `.bd` 文件。

```powershell
cd DBEM1
.\x64\Release\DBEM1.exe
```

常用验证运行方式：

```powershell
python tools\run_rod120_100_analytic_validation.py --repo . --exe DBEM1\x64\Release\DBEM1.exe --in-place --timeout 900
```

该脚本会临时写入 `DBEM1/BEM_DATACARD.DAT`，运行结束后恢复原文件。

## 多域参数卡格式

旧参数卡不追加多域段时，程序按单域路径运行。启用多域时，在原 GMRES 参数后追加：

```text
DomainCount
domainId E v Rou
...
InterfaceCount
domainA eleA domainB eleB normalSign
...
```

说明：

- `DomainCount` 为子域数量。
- 每个子域建议显式给出材料行，材料参数为 `E / v / Rou`。
- 当前代码兼容缺失材料行并回退全局静态材料；多域验证和正式算例不建议依赖该回退。
- `domainId` 支持 0 基或 1 基编号，但同一段内应保持一致。
- `InterfaceCount` 可省略；省略时按相邻子域自动生成界面。
- 显式界面行中 `normalSign=-1` 表示 B 侧面力相对 A 侧取反。

`rod120_100` 当前示例为 5 个均质子域：

```text
5
1 20 0 20
2 20 0 20
3 20 0 20
4 20 0 20
5 20 0 20
16
1 5 2 25 -1
...
```

## 多域求解技术要点

### 材料局部化

`BuildDynaMat(v, G, Rou, dt)` 统一构造动态材料参数，包括：

- `C1 / C2`
- `Dt / Dt2 / Dt_2 / Dt_3`
- `C1Dt / C2Dt`
- 动态基本解需要的倒数和组合系数

多域路径中，每个 `Domain` 持有自己的 `DynaMat`。动态积分、波前判定、Step0/current 装配、history 装配通过 `DSquareElement::ActiveDynaMat()` 获取当前域材料。

并行装配中使用：

```cpp
thread_local const DynaMat* DSquareElement::m_ThreadDMat;
```

worker 进入域装配时设置线程局部材料，退出时恢复，避免写全局 `DSquareElement::m_DMat` 造成线程间污染。

### 自由度映射

当前多域 GMRES 路径按常值单元处理，每个边界单元中心对应一个 3 自由度块。

外边界：

- `BCID = 123`：给定位移，未知面力。
- `BCID = 456`：给定面力，未知位移。

界面：

- A、B 两侧共享界面位移块。
- A、B 两侧共享一套面力块，B 侧按 `normalSign` 取符号。
- 当前只支持共形界面的一一单元配对。

### 矩阵装配

`MultiDomainCCSRBuilder` 负责块级聚合：

- key 为 `(rowBlock, colBlock)`。
- value 为 3x3 子块。
- 重复块累加。
- 可输出 `CCSRMat`。
- 对 `MGS` 可尝试输出 `SymCCSRMat`；若 3x3 块不满足相对容差下的对称性，则回退为 full CCSR。

近零块和对称判断使用混合容差：

```text
threshold = absTol + relTol * scale
```

当前默认：

```text
absTol = 1.0e-14
relTol = 1.0e-10
```

### 多域 GMRES 时间递推

多域 GMRES 求解器入口：

```cpp
DynaGMRESSolverMultiDomainCCSR(...)
```

主要流程：

1. 构造 `GlobalDofMap`。
2. 装配 Step0/current 未知矩阵和已知矩阵。
3. 装配历史矩阵 `MTS[step] = T2 + T1` 和 `MGS[step] = U`。
4. 对每个时间步构造 RHS：
   - 当前已知边界项。
   - 历史 `G*T`。
   - 历史 `T*U`。
5. 调用现有 `gmres(..., CCSRMat& A, ...)`。
6. 反写 `BoundaryValue`，输出界面误差和诊断数据。

### pthread 并行装配

Step0/current 和 history 装配支持按域切分 pthread worker。常用环境变量：

```text
DBEM_MULTIDOMAIN_GMRES_PTHREAD_CHECK=1
```

开启后会比较 pthread 装配和串行 fallback 的差异；若差异超过混合容差，会回退串行矩阵。

关闭 pthread 装配可设置：

```text
DBEM_MULTIDOMAIN_GMRES_PTHREAD=0
```

## 验证与诊断输出

开启验证输出：

```text
DBEM_VALIDATION_OUTPUT=1
```

主要输出文件位于算例输出目录下：

- `validation_state.csv`：每步、每单元代表点的位移/面力状态。
- `validation_metrics.txt`：求解器、时间步、材料、界面误差、多域矩阵和性能指标。
- `interface_transfer_audit.csv`：界面两侧位移连续与面力平衡审计。
- `rhs_breakdown_*.csv`：RHS 当前项、历史项、总 RHS 和残差诊断。

`validation_metrics.txt` 中关键字段：

- `DomainMaterial[*]`：各域 `C1/C2/C1Dt/C2Dt`。
- `MultiDomainFinalFlag`：GMRES 最终标志。
- `MultiDomainSolvedSteps`：已求解步数。
- `InterfaceOverallMaxAbsU`：界面位移最大不连续。
- `InterfaceOverallMaxAbsTBalance`：界面面力平衡最大误差。
- `MultiDomainMatrixStructureValid`：矩阵结构诊断。

## 解析验证

`tools/run_rod120_100_analytic_validation.py` 用于验证多域材料是否真正生效：

1. 使用 `rod120_100` 作为算例。
2. 临时污染全局静态材料行，例如把全局 `E` 或 `Rou` 改成错误值。
3. 保持参数卡多域材料行不变。
4. 运行多域 GMRES。
5. 读取 `validation_state.csv`。
6. 用多域材料行计算一维杆解析解。
7. 输出最大绝对误差、最大相对误差和材料日志检查结果。

示例：

```powershell
python tools\run_rod120_100_analytic_validation.py --repo . --exe DBEM1\x64\Release\DBEM1.exe --in-place --timeout 900
```

通过条件：

- 数值解与解析解误差在脚本阈值内。
- `validation_metrics.txt` 中 `DomainMaterial[*]` 与多域材料行一致。
- 运行日志中的 `MultiDomain GMRES material domain=*` 与多域材料行一致。

## 当前限制

- 自适应 `dt` 尚未实现。
- 当前 `dt` 仍由全局静态材料行和 `Beta` 先计算得到；多域材料用于各域 `DynaMat` 和积分/装配，但时间步长控制尚未局部化。
- 多域路径只覆盖 GMRES 并行求解链路，不重构旧 Gauss、ACA、单域 GMRES。
- 当前多域 GMRES 只支持整块边界条件 `123` 和 `456`。
- 多域路径暂不支持 batch boundary mode。
- 当前界面要求共形、一一配对。
- 当前常值单元路径只使用每单元一个代表物理节点，不是完整 8 节点自由度版本。
- 非匹配时间步、多速率自适应时间网格、时间 mortar 接口耦合仍停留在方案讨论阶段。

## 后续方向

下一阶段如果继续做时间步控制，建议不要采用“所有域统一子步数”的方案，因为它等价于全局缩小 `dt`，无法体现多域自适应收益。

真正的目标应是：

- 每域按材料波速和局部误差选择自己的 `dt_d`。
- 历史卷积按真实物理时间区间重构。
- 接口用隐式时间 mortar 或等价投影方法耦合非匹配时间点。
- 输出宏步状态时，对每域局部时间网格做插值。

这部分尚未实现，不能作为当前数值结果的通过标准。
