# 波前积分修复：远程应用记录（2026-10-19）

已通过SSH连接tdbem-server，核验服务器四个文件与修复前版本一致后，备份并应用近零积分判据及保守筛选修复。服务器为Windows，使用现有MSVC v142、Release/x64配置，独立bin/obj完整编译，验证后更新默认程序。

- 安装时间：2026-10-19T20:08:03.0438908+08:00。
- 项目：`D:\wtj\3D_MultiDomain_TDBEM`。
- 默认程序：`D:\wtj\3D_MultiDomain_TDBEM\x64\Release\DBEM1.exe`。
- 新程序SHA256：`099DCA679007633E21E309BCD2BF07F92553D5F8BFA0293CC2A353BA1FCA1778`。
- 原程序SHA256：`979701F1CBBE45E3D00847421E9A4FD8982A3D4AA7601F1CFC7A3BDA9A5D2FCA`。
- 源码备份：`D:\wtj\3D_MultiDomain_TDBEM\backups\wavefront_budget_fix_20261019_200100\source_backup`。
- 默认EXE/PDB备份：`D:\wtj\3D_MultiDomain_TDBEM\backups\wavefront_budget_fix_20261019_200100\binary_backup`。

同步文件为DBEM1/wavefront_quadrature.cpp、wavefront_quadrature.h、solver_policy.h、pre_gmres.cpp。每块200万次求值上限、相对容差系数和时间解析核保持原值；绝对项采用整元特征尺度保护，旧矩阵缓存的积分修订号升为2。

服务器验证：原失败T1块（单元22、lag65）用203次求值完成；378个波前擦边、时间步扰动及尺度缩放块通过独立参考校验，最大实际误差/允许误差约0.00264884；正式生产版本完成34单元8步短程求解，输出170条有限状态记录。

本次服务器验收范围为上述局部积分和8步求解；原E08的1000步完整回归为此前本地隔离运行结果，不能混作服务器本次重跑。历史实验结果和原图保留。

远程详细记录：`D:\wtj\3D_MultiDomain_TDBEM\backups\wavefront_budget_fix_20261019_200100`下的deployment.json与verification_records.zip，后者保存验证结果及编译日志。本地副本位于`validation_runs/wavefront_budget_20261019/remote_sync_20261019_200100`。长期稳定性不由本次安装或短程验证判定。


## 临时文件清理

2026-10-19 20:11（北京时间）已删除本次临时编译、测试和传输目录，回收约62.73 MiB。源码、旧EXE/PDB备份及34项验证记录已归档到`D:\wtj\3D_MultiDomain_TDBEM\backups\wavefront_budget_fix_20261019_200100`；归档逐项通过SHA256校验。默认EXE及四个源码文件哈希保持不变。
