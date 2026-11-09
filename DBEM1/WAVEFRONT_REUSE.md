# 精确积分路径的几何复用

本改进建立在622e2ef的现有波前后端上。该版本已将同一角度的射线裁剪和三角函数移出径向循环；本次没有重复实现或把已有收益算作新增加速。

新增两层线程私有复用：

1. 同一不可变面元几何及源点的整片保守距离/面积/Jacobian界限，在Classify、Integrate与完整参数片评估间复用。此结果不依赖材料和时间窗；变化后的材料尺度与核仍按当前调用计算。
2. 同一源点—面元、波速和物理时间窗的平面波前分区，在U(n)、T1(n)、T2(n+1)之间复用。T2(n)属于前一窗，不可与T1(n)混用。每种核独立运行5/9规则和自适应误差控制。

每线程每层最多两项。分区仅缓存不超过4096个初始片的计划，超出时照常构造和积分、不保留计划。几何所有者以shared_ptr保留，避免释放后指针复用误命中；GeoInit创建新几何时自然失效。源点按精确浮点位模式匹配，分区还匹配C1/C2、窗端点及全部波前。没有核积分结果缓存，也没有跨核共享误差估计、叶子队列或计算预算。

`DBEM_WAVEFRONT_REUSE_PARTITIONS`默认为1，设为0可关闭上述新增复用作严格对照。Statistics新增partitionBuilds/partitionCacheHits及boundsBuilds/boundsCacheHits。修改了Result中的统计布局，必须完整重编译，不混用旧对象文件。

## 构建与验证

Windows Python（可从WSL调用）运行：

```
python tools/validate_wavefront_reuse.py --dependency-root <含include和lib的项目目录> --output <独立输出目录> --build
python tools/smoke_wavefront_reuse.py --dependency-root <项目目录> --output <同一输出目录>
```

独立构建器build_gmres_openmp.py支持`--reuse-partitions 0/1`和`--dependency-root`。完整时程驱动沿用项目已有小杆及moderate_v2模板（需要保留对应研究数据目录），交替先后顺序运行缓存开/关，固定输入、dt、线程、核规则和真实残差门槛。

核块测试：210个独立参考块，重复调用、U/T1/对齐T2、源点/波速/时间变化、同波速不同密度/刚度、GeoInit重建、预算失败、容差变化、曲面路径和双线程隔离。缓存开关输出全部数值、自适应层数与核求值次数逐字节比较。GMRES/RHS既有回归分别运行。

完整时程：82单元56步三组配对；328单元112步两组配对。位移和面力完整CSV及输入摘要必须完全相同；全部当步真实残差通过1e-10。短时程另比较单域、双域同材、双域异材的1/2线程和缓存开/关。

本次工作区为codex/wavefront-reuse-20261104，原主工作区分支保持原状。完整实测结果另见本工作区测试证据报告；一次统计中的工作量降幅不能直接当作完整求解加速比。新增复用保持原物理解，不修正已知的时空离散误差。
