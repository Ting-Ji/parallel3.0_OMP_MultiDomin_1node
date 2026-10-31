# GMRES 可选 CUDA 后端

2026-10-28：已经实现并在本机真实 GPU 上验证。CPU 仍为默认后端；原生产参数卡、输入和既有 EXE 未被替换。该目录的 DLL、可选 EXE 和试验输出独立存放于 `build/`、`validation/`，不自动提交或上传。

## 当前 CPU 并行到底在哪里

- `DBEM1/pre_gmres.cpp`：预条件器构造中的部分叶块工作、求逆及满足独立输出条件的 `GMRES_M_X` 使用 OpenMP。
- `DBEM1/multidomain.cpp`：历史项 SpMV 使用 pthread 工作线程池。这属于时间步右端项构造，不是 GMRES 的 Arnoldi 内循环。
- `DBEM1/mat_vec.cpp` 的 `SparseMul(CCSRMat&, ...)` 是串行循环；当前 CCSR GMRES 调用的就是它。
- Arnoldi 正交化、向量点积和更新等在原 GMRES 中主要仍为串行循环。因此当前实现属于部分 CPU 并行，不能称为整套 GMRES 已充分并行。

## 本次实现范围

1. 在 CCSR GMRES 入口增加可选后端。将原 3×3 块矩阵转为等价标量 CSR，保留 FP64。
2. 将现有叶块逆预条件器转为等价稀疏线性算子，保留 `IID`、输入/输出节点映射和重叠输出的累加；不是换成对角预条件器来降低计算量。
3. GPU 上执行 A*x、P*x、点积、范数、缩放、AXPY 和 Krylov 基组合。小 Hessenberg 系统、Givens 旋转及迭代调度保留在 CPU。
4. 多域求解器（包含走该框架的单域算例）用作用域对象管理 GPU 会话，A/P 上传一次，在多个时间步复用。旧 CCSR 入口也可调用 GPU，但每次调用建立临时会话，尚未优化跨时间步复用。
5. 每个候选解使用未预条件化的真实残差判停；返回后再用原 CPU `SparseMul` 复核。不因估计残差较小就接受错误解。
6. GPU 使用标准重启左预条件 GMRES：大向量驻留显存，CPU 上回代小上三角系统。原 CPU 实现用增量 W 更新及末步矩阵求逆；数学子空间与预条件方向相同，浮点运算顺序不同，不要求逐位一致。

本次没有迁移矩阵装配、历史卷积、预条件器构造或完整后处理。普通密集 `Wmatrix`、普通 `CSRMat` 的旧 GMRES 重载仍是 CPU 路径；切换范围明确限定为 CCSR。

## 构建与运行

本机工具：Windows CUDA 12.6、Visual Studio 2022 x64、Windows Python 3.12；GPU 为 GTX 1050 Ti 4 GB。无须在 WSL 安装第二套 CUDA。编译辅助程序调用 Windows 专属的 nvcc/MSVC，其余本地操作可继续使用 WSL。

在 WSL 中进入项目后运行：

```bash
cd "$(wslpath 'D:\博二下学期\parallel3.0_OMP_MultiDomin_1node')"
WIN_PY=/mnt/c/Users/33435/AppData/Local/Programs/Python/Python312/python.exe
"$WIN_PY" "$(wslpath -w "$PWD/gpu/build_windows.py")" --solver --benchmark
"$WIN_PY" "$(wslpath -w "$PWD/gpu/test_cuda.py")"
"$WIN_PY" "$(wslpath -w "$PWD/gpu/run_validation.py")"
```

`build_windows.py` 自动定位现有 Visual Studio，不安装 SDK、不修改系统环境。默认 `--arch sm_61`；其他机器可显式指定目标。基准和完整求解器编译使用原项目源文件，输出为 `gpu/build/DBEM1_cuda_optional.exe`、`benchmark_project.exe` 和 `dbem_cuda_gmres.dll`。

选好自己的参数卡后，可手动运行：

```bash
"$WIN_PY" "$(wslpath -w "$PWD/gpu/run_solver.py")" --backend cuda
# 同一可选 EXE 的 CPU 路径
"$WIN_PY" "$(wslpath -w "$PWD/gpu/run_solver.py")" --backend cpu
```

默认工作目录为项目 `DBEM1`，也可通过 `--workdir` 指定独立算例目录。启动器不改参数卡；求解器按原逻辑生成输出。原程序成功返回码为1，启动器保留该约定。

后端由 `DBEM_GMRES_BACKEND=cpu|cuda` 显式选择，未设置时为 CPU。`DBEM_CUDA_GMRES_DLL` 可指定插件路径，默认查找 EXE 同目录插件，支持中文路径。直接启动 EXE 时，CUDA `bin` 应在该 Windows 进程 PATH 中。CPU 模式不加载 CUDA 插件，验证时已用不存在的插件路径检查这一点。

## 精度、生命周期和失败行为

- 只使用双精度，不启用 FP32/TF32 或混合精度替代。
- 保留原矩阵和原预条件器；矩阵、预条件器在 GPU Scope 内必须保持不变。新算例创建新作用域，不以可能复用的裸指针永久缓存数值。
- 收敛检查为 `||b-Ax||₂ / ||b||₂ <= tol`，使用调用方容差。CPU 最终复核不通过则不能报告成功。零右端项、精确初值、奇异/停滞及非有限输入有独立检查。
- Krylov 大向量均驻留显存，每步仍有小标量/小系数数组传输与同步。因此小系统可能明显更慢。
- A 与 P 使用独立的 cuSPARSE 外部工作区，零字节工作区保留空指针。算子串联检查覆盖 `P*(A*x)`，防止只验证两个独立乘法而遗漏工作区复用问题。
- 显式矩阵/工作区分配限制在可用显存预算内；超预算、缺 DLL、无设备或 CUDA 错误会明确报错，不自动改精度或悄悄退回 CPU。
- 打印的 `device_MiB` 仅统计显式分配，不包含驱动和 NVIDIA 库上下文的全部显存。

## 验证与性能记录

详见同目录 `RESULTS.md`。核心检查包括：

- DLL 级数值测试及 A/P/PA 算子检查、不同初值、重启、多个右端项、零右端项、非有限输入和奇异系统。
- 项目原 CCSR CPU GMRES 与 GPU 的合成矩阵对照，含非恒定预条件器及非平凡输入/输出映射。
- 14面元单域和20面元两域物理长杆，各24步，4 CPU线程；两域含一个人工接口。设置 beta=0.8，触发多次迭代和跨域波传播，而非仅检查首步局部响应。
- 上述试验输入、原始输出、日志和 JSON 保存在 `gpu/validation/`，运行脚本只清理自己拥有的试验输出，不操作生产输入输出。

本次比较以当前项目 CPU 实现为基线，未与充分优化的多线程 CPU SpMV/MKL 比较；合成系统加速不代表完整边界元模型的端到端加速，更不代表更改了物理或时间离散精度。

NVIDIA 接口依据：[cuSPARSE 12.6 SpMV](https://docs.nvidia.com/cuda/archive/12.6.0/cusparse/generic-api/generic-api-functions.html)、[cuBLAS 12.6](https://docs.nvidia.com/cuda/archive/12.6.0/cublas/index.html)。

## 实际 TSV 扩展验证

2026-10-29 的 TSV 01/04/09 数值对照见 [TSV_RESULTS.md](TSV_RESULTS.md)。本地参数卡统一为 40 线程，只替换全局铜材料参数 E=110、nu=0.34、rho=8.96，各域硅/铜材料和原 1000 步边界历程保留。生产参数卡不改动；脚本将输入复制到隔离目录并记录 SHA256。

复现时通过 Windows Python 调用 `gpu/run_tsv_validation.py --cases 4 1 9 --threads 40 --output <新的隔离目录>`。每个后端顺序运行，包含完整场量核对；已有结果只在源码/二进制身份匹配时复用。求解器子进程无窗口启动并关闭标准输入，避免启动时依赖控制台状态。

本轮存在其他 GPU 负载，计时只保留为原始记录。用户要求 TSV 09 结束后停止，TSV 16 和空闲性能复测未启动。线性收敛、CPU/GPU 全时程差异和物理/长期稳定性分别评价，不把 flag=0 当作全部验证通过。

随后用户确认背景 GPU 负载已消除，重新完成 TSV 04 的 40 线程、1000 步 CPU/GPU 对照，见 [TSV_IDLE_RESULTS.md](TSV_IDLE_RESULTS.md)。GMRES 为 CPU 7.713 秒、GPU 7.793 秒，未加速；GPU 侧整段利用率采样中位数为 0%，大量时间仍在 CPU 和输出。逐后端的完整状态 CSV 与上轮哈希完全一致，后期差异仍可复现。按最新指令，在本组比较结束后停止，本轮没有启动 TSV 09/16。
