# 本轮研究实现与关键证据

保存日期：2026-11-03。此目录为独立分支的精简证据副本；不改动原研究目录。

- 核块、串并行与短时程回归：kernel_acceptance.json。
- 几何294组、1764项标量积分及解析复核：geometry_*。
- 完整几何覆盖证明及TikZ图：几何操作的正确性与完备性.tex。
- 三组已完成大算例旧/新与反射解析解的物理误差：large_analytic_*、large_peak_timing.json；2000域800步按用户要求停止，不列为完成。
- 长杆候选及实际未达目标的结果：rod/。初次切向加载已被原目录标为无效轴向对照，此处保存核验正确轴向加载后的结果。

核块正确并不保证完整解更准确。1000域1000步端部L2为旧13.2274%、新15.5167%；小杆轴向L2为旧13.3525%、新11.5103%。不将本实现标为已优于原方法。

复现实际核块回归：Windows Python运行tools/test_wavefront_integration.py。复现小杆：Windows Python运行本目录rod/run_small_rod.py --axial（使用安装的MSVC；仅在独立输出目录构建）。大算例解析后处理脚本tools/compare_wavefront_large_analytic.py依赖本地保留的大CSV，不会重跑求解器。

可执行文件、完整大CSV、构建中间文件和所有原实验结果保留在本地tmp及wavefront_topic忽略目录，没有放入Git。桌面组会TeX项目也保留原位置。
