# RHS 构造 CPU 并行实现与验收

实现分支：`codex/rhs-parallel`。新路径由 `DBEM_RHS_CPU_MODE=parallel` 选择，
未设置或设为 `legacy` 时保留原 RHS 路径。GMRES 的算法、接口和默认后端没有改变。
本轮物理验收显式设置 GMRES 为 `omp`，全部使用 40 线程。

## 实现与所有权

`multidomain.cpp` 的求解会话显式持有 `RhsCpuWorkspace`。初始化时按原遍历顺序
建立已知边界目标块到源节点/符号的映射，重复目标保留最后一次赋值，未映射目标置零。
KnownM 和各历史矩阵按累计非零块数加行数划分连续区间，区间在整个时间步循环中复用。

扩展原有常驻 pthread 历史乘法线程池，每个时间步提交一次完整 RHS 任务：
已知边界收集 → KnownM 乘法 → 按滞后项顺序构造 U/T 卷积向量、G 加、T 减 →
有限值检查及 x0/x 清零。阶段间使用可复用条件变量屏障，每个输出位置只有一个写入者。
RHS 返回后所有工作线程休眠，随后 GMRES 使用它已有的 OpenMP 线程组。

卷积每个元素从零开始，仍按原来的 1、2、1 顺序累加有效历史步。
G/T 每行先按原顺序算完矩阵乘局部结果，再对 RHS 加减；不跨滞后项归约。
CCSR/SymCCSR 的块内和行内运算表达式保持原样。工作线程不分配内存、不写文件、
不从屏障间提前退出；非有限值在全部线程完成任务后由主线程报告。
线程初始化失败沿用现有显式串行回退策略，默认停止；部分创建成功的线程会先被回收。
工作区在求解结束释放，不使用全局地址缓存。

## 构建与算子测试

从 WSL 在项目根目录调用现有 Windows Python/MSVC：

```sh
WIN_PY=/mnt/c/Users/33435/AppData/Local/Programs/Python/Python312/python.exe
"$WIN_PY" tools/build_rhs_parallel.py
(cd tmp/rhs_parallel/build && ./rhs_parallel_test.exe && ./gmres_openmp_test.exe)
"$WIN_PY" tools/run_rhs_parallel_acceptance.py --audit-only
"$WIN_PY" tools/run_rhs_parallel_acceptance.py --detach
```

独立目录 `tmp/rhs_parallel/build` 使用与旧基线相同的 MSVC 编译参数，包括
`/O2 /openmp /fp:precise`、独立真实残差复核和关闭 Tecplot 输出。完整场量 CSV 仍输出。
`build_manifest.json` 保存参数、源码及可执行文件哈希。

算子测试调用实际工作区和线程池，以原串行卷积、原 SparseMul 和原向量累加为独立参照，
通过 `memcmp` 检查逐值/位模式一致。覆盖符号、重复/未映射目标、索引越界、
不均匀稀疏行和空行、行数小于线程数的任务分区、CCSR/SymCCSR、初始及截断历史步、
G/T 顺序、诊断向量、反复提交、清零、NaN 拒绝及拒绝后再次提交、部分线程创建失败、
重新初始化和重复关闭，以及模式选择。GMRES 原算子/求解器测试同时保留。
测试输出中的模拟线程创建失败是预期的故障注入，最终 PASS 才表示测试通过。

当前构建 SHA256：`cd72bccf91ddbe5ad43f1c8b4773afac36327b58237039c3c9e1d8f673e3843b`。
RHS 与 GMRES 两套测试均通过，日志在构建目录 `rhs_test.log`、`gmres_test.log`。

## 计时口径

`multidomain.rhs` 保留原有外层口径；一次性 `rhs.cpu.workspace_setup` 同时包含在
`multidomain.state_init`、`run.solve` 和 `run.total` 中，不能重复相加。

`DBEM_RHS_CPU_PROFILE=1` 启用细粒度计时：

| 记录 | 含义 |
|---|---|
| rhs.cpu.known_gather | 已知边界收集 |
| rhs.cpu.known_spmv | KnownM 乘法 |
| rhs.cpu.convolution | U/T 卷积向量构造 |
| rhs.cpu.g_spmv_accumulate | G 乘法及融合累加 |
| rhs.cpu.t_spmv_accumulate | T 乘法及融合累加 |
| rhs.cpu.check_zero | 有限值检查、初值清零 |
| rhs.cpu.worker_wait | 工作线程阶段屏障等待 |
| rhs.cpu.task_wall | 主线程提交至完成的墙钟时间，逐步记录 |
| rhs.cpu.dispatch_completion | task_wall 减去该步最长工作线程任务时长的非负余量 |

每个工作线程在内存中累计，求解结束由主线程写各阶段的最大线程累计时间。
这些最大值和等待时间存在重叠，不可相加解释为完整阶段墙钟；dispatch_completion 也是
调度/完成开销的估计，不能精确分离操作系统调度。诊断向量额外工作仅在原诊断开启时存在。
正式验收关闭 RHS 和 GMRES 细粒度计时；不会额外重跑诊断物理算例。

## 已有基线复用与数值门槛

旧证据只读：`tmp/gmres_openmp/acceptance_v2` 中各模型、步数的 `omp_*_0`。
旧二进制 SHA256：`41d545684ffe41d326848ffcdfec559d4467426bb07248acf7e784cdcf6d9611`。
它的 GMRES 已并行，RHS 使用本轮改造前路径，作为正确的增量性能参照。
不使用旧 `legacy` GMRES 或额外 repeat 结果，不重新启动旧基线求解。

队列启动前审计全部八项：编译参数、未修改的非 RHS 源码哈希、输入清单和实际输入文件、
40 线程、旧 OMP/repeat=0、完整场量文件存在、逐步残差及计时证据。
`baseline_match.json` 保存匹配结果和原始输入 SHA256。二进制改变是本轮预期行为。
如证据缺失或不匹配则停止并报告，不补跑基线。

门控顺序：杆件单/双域各 24 步，TSV01/04/09 各 1000 步；
1000/2000 域立方体各 100 步；1000 域立方体 1000 步。
TSV 隔离卡全局使用铜 E=110、nu=.34、rho=8.96，域材料不变；立方体保留原材料。
每个新版本模型/步数组合只运行一次。队列使用排他锁，不自动重跑未完成或失败目录。

独立原始 SparseMul 逐步核对真实残差，必须满足原参数卡容差。
完整位移及牵引分别要求全时程相对 L2 <=1e-6，最大绝对差异 <=基线全局峰值的1e-5，
零场下限1e-30。保留全部场量，逐行身份、有限性及迭代数均核查。
旧结果无逐步 RHS 指纹，因此不宣称物理时程 RHS 本身逐值一致。
任何门槛失败停止后续算例并保留证据，不放宽标准或自动重算。

RHS、GMRES、计算阶段及全程均按单次时间比报告，另报告工作区准备、峰值采样内存和
迭代数变化。跨轮单次对照不能排除机器状态差异；不取平均或中位数。
两个大规模 100 步 RHS 加速比均达到1.5为工程目标，不达标仍如实保留结果及显式选项。
不自动修改默认后端。历史物理精度告警不属于本轮已经解决的问题。

## 实时证据

完整物理验收已启动，结果随队列更新，尚未宣称全部通过。

- `tmp/rhs_parallel/acceptance/acceptance_report.md`：实时数值/性能结果。
- 同目录 `acceptance_queue.json`、`small.json`、`large-short.json`、`large-full.json`。
- 各运行目录 `result.json`、`input_manifest.json`、`run.log` 和完整输出。

上述路径均在新目录内；旧验收结果和输入没有被修改。
