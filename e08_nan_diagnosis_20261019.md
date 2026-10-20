# E08 非有限值排查（2026-10-19）

本次“单域路径”明确指新多域框架的DomainCount=1。此前把旧入口的失败归于新框架单域，范围不准确，现据代码及输出更正。

|项目|用户手动运行|失败重跑|
|---|---|---|
|输入卡|有域数1、材料行、0个接口|只有原16项，没有域配置段|
|MultiDomain.IsActive()|true|false|
|入口|DynaGMRESSolverMultiDomainCCSR|DynaGMRESSolverNewCCSRNewPth|
|残差标记|multidomain|single|
|奇异初始化|串行|OpenMP并行|

## 新框架域数1的结果

已扫描用户手动结果的1000步残差及2,016,000条逐分量诊断。解、Ax、残差、已知右端项、历史项和总右端项均有限；最大相对残差9.9928697224e-10，末步约2.32417e-10。

归档目录没有validation_state.csv，所以本次确认求解诊断有限，不据此证明状态导出及后处理全部完成。timing.txt的时间戳属于另一次尝试，不能与这份残差拼接计算完整耗时。

E08原672单元、材料、时间步长，缩短至前2步并显式启用新框架1域：修复前后各以1/4/8/64线程运行，各线程矩阵及状态逐字节一致，修复前后亦逐字节一致，最大相对残差约9.01435e-10。

## NaN来自保留旧入口

失败缓存MTS0第395号单元的自作用3x3块有9个NaN。首步RHS有限，但解和残差已为NaN，第二步历史项被污染。maxAbsResidual=0由NaN比较为假造成，不能表示收敛。

与此前本地旧入口缓存相比，非自作用块全部逐字节一致；当前矩阵238个自作用块、已知项矩阵5个自作用块存在明显差异。

旧入口并行PthPhyInit会共写四组类静态数组：m_OnEleN、m_OnEleJ、m_OnEleR、m_StaticTij。线程互相覆盖距离、法向、Jacobian及静态积分结果，导致错误有限值或NaN。新框架当前串行初始化，未触发该竞争，符合用户观察。该竞争早于普通波前积分后端改造存在。

96单元对照试验复现旧问题，最大分量差约1.486；四数组改为thread_local后，1/4/8/64线程各重复4次，均逐字节匹配串行参考。此前本地旧入口1000步“完成且有限”不能证明这些自作用系数正确，不能作为新框架单域验收结果。

## 修复和验证

四数组改为线程私有；不改奇异公式、积分参数及时间支撑。旧矩阵缓存修订号升至3。增加矩阵、GMRES范数、诊断与最终状态输出的非有限值检查，阻止NaN传入历史。

注入NaN的RHS、预条件块、矩阵、状态输出均被拒绝；正常小系统求解通过。旧入口E08前2步在本地和服务器均以1/4/8/64线程通过，矩阵和状态逐字节一致。独立稀疏LU首步残差4.19e-16，GMRES与其相对解差约7.54e-10。

源码已同步并独立编译五个文件：DBEM.cpp、dsquareelement.cpp、dsquareelement.h、pre_gmres.cpp、multidomain.cpp。安装保护发现用户默认EXE已变化，阻止了替换；默认程序仍为用户手动版本0EE2E3A9512C486ED883CBF373AB08F0A18DDB263A4FBBC3DA3B3227ACAA7072。手动卡片与原结果未改。

后续精度、性能、长期稳定性以新框架域数1为验收目标。残差有限、解有限与长期稳定需要分别检查。

证据位于validation_runs/e08_nan_diagnosis_20261019/：matrix_audit.json、manual_new_framework_audit.json、new_framework_comparison.json、independent_first_step.json及before/after控制试验。

后续收尾及停止记录见 [legacy_single_completion_20261019.md](legacy_single_completion_20261019.md)。用户要求停止后未继续数值测试；E08完整验证未完成，默认EXE未替换。本次临时输出已清理，已完成检查摘要保留。
