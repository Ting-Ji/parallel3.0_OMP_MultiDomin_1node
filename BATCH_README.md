# 连续计算

在项目根目录编辑 `cases.txt`，每行写一个模型名（不带扩展名）。同名模型位于不同子目录时，写相对于 `DBEM1/input` 的路径，例如 `cube1000_million/cube_md1000_face169`。模型必须有同名的 `.model`、`.bd` 和 `.DATACARD`。空行和以 `#` 开头的行会被跳过。

先检查输入：

```text
python run_batch.py cases.txt --dry-run
```

在 Windows Python 中开始连续计算（计算程序是 Windows EXE）：

```text
python run_batch.py cases.txt
```

默认 EXE 为 `DBEM1/x64/Release/DBEM1.exe`；需要更换时使用 `--exe 路径`。清单和 EXE 的相对路径均以项目根目录为基准。脚本逐项替换 `DBEM1/BEM_DATACARD.DAT` 并等待本次计算结束；完成后恢复原参数卡。原有求解结果仍由 DBEM 自行保存。每批的控制台日志、公共 `log.txt` 副本、参数卡备份和汇总 CSV 位于项目根目录 `output/batch_时间/`。

若强制关闭脚本导致参数卡未恢复，先确认 DBEM 进程已停止，再按 `output/.dbem_batch.lock` 中的 `backup` 路径恢复 `DBEM1/BEM_DATACARD.DAT`，最后删除该锁文件。同名、同时间步数的算例重复运行时，DBEM 原有结果可能被覆盖。

## C 版本

`run_batch.c` 与 `cases.txt` 位于项目根目录。在 **Visual Studio x64 开发者命令提示符**中，进入项目根目录并编译：

```text
cl /nologo /W4 /O2 /Fe:x64\Release\run_batch.exe run_batch.c
```

然后运行（可从其他工作目录启动；程序会从 `run_batch.exe` 位置向上找到项目根目录）：

```text
x64\Release\run_batch.exe cases.txt --dry-run
x64\Release\run_batch.exe cases.txt
```

C 版默认调用项目根目录的 `x64\Release\DBEM1.exe`，也支持 `--exe 路径`。DBEM 的计算输出会实时显示在批处理程序的控制台中；当前算例结束后才会开始下一个。失败继续，最后恢复原参数卡。每个算例的 DBEM 日志副本、备份和 `summary.csv` 保存在项目根目录的 `output/batch_时间/`。C 版使用 Windows API，因此需在 Windows 上编译和运行。

复制到另一台电脑时，保留项目根目录中的 `cases.txt`、`DBEM1` 目录以及 `x64/Release` 目录的相对位置。至少需要 `x64/Release/run_batch.exe`、`x64/Release/DBEM1.exe`、`DBEM1/input` 中列出的三件套、`DBEM1/BEM_DATACARD.DAT`、`DBEM1/gaussposition` 和 `DBEM1/gaussweight`。双击新版 `run_batch.exe` 时窗口会停留，显示缺少的具体文件。也可以在命令提示符中先运行 `x64\Release\run_batch.exe cases.txt --dry-run`。

批处理程序本身只依赖 Windows 系统库；被调用的 `DBEM1.exe` 还需要 x64 的 `pthreadVC2.dll` 以及 Visual C++ 运行库（包括 OpenMP 的 `VCOMP140.DLL`）。新版程序会在计算前检查这些依赖并报告缺项。`pthreadVC2.dll` 可放在 `DBEM1.exe` 同目录；Visual C++ 运行库应在目标电脑安装对应的 x64 版本。
