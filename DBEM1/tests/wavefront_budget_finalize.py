"""Finalize the E08 evidence once the unchanged full solve has finished."""
from pathlib import Path
import json,re,shutil,zipfile,os,hashlib
ROOT=Path(__file__).resolve().parents[2]
work=ROOT/'validation_runs/wavefront_budget_20261019/e08_full_20261019_173023'
r=json.loads((work/'summary.json').read_text())
case=Path(r["case"])
for ext,expected in r["input_sha256"].items():
 path=case/"BEM_DATACARD.DAT" if ext=="DATACARD" else next(case.joinpath("input").glob("*."+ext))
 assert hashlib.sha256(path.read_bytes()).hexdigest()==expected,(ext,path)
assert hashlib.sha256(case.joinpath("DBEM1.exe").read_bytes()).hexdigest()==r["executable_sha256"]
assert r.get('passed') and r['last_state_step']==1000,r
message=("<!-- E08_RESULT_BEGIN -->\n"
 "本地原 E08（672单元、64线程）已完成全部1000阶历史矩阵装配及1000步时程计算。"
 "返回码1，最后输出步1000，共%d条状态记录，全部位移/面力为有限数值；没有再次出现核函数求值上限错误。"
 "本机隔离运行耗时约%.1f秒，期间有其他回归占用CPU，不作为服务器性能基准。"
 "model、bd、DATACARD哈希与原输入一致。"
 "完整记录见 validation_runs/wavefront_budget_20261019/e08_full_20261019_173023/summary.json。"
 "这一结果证明计算完成，不代表长期稳定性已通过。\n<!-- E08_RESULT_END -->")%(r['state_rows'],r['seconds'])
report=ROOT/'wavefront_budget_fix_20261019.md'
s=report.read_text(encoding='utf-8')
pattern=r'<!-- E08_RESULT_BEGIN -->.*?<!-- E08_RESULT_END -->'
assert len(re.findall(pattern,s,re.S))==1
report.write_text(re.sub(pattern,lambda m:message,s,flags=re.S),encoding='utf-8')
bundle=ROOT/'validation_runs/wavefront_budget_20261019/repair_bundle'
shutil.copy2(report,bundle/'README.md')
archive=ROOT/'wavefront_budget_fix_20261019.zip';temporary=archive.with_suffix('.zip.tmp')
with zipfile.ZipFile(temporary,'w',zipfile.ZIP_DEFLATED) as z:
 for f in bundle.rglob('*'):
  if f.is_file():z.write(f,f.relative_to(bundle))
os.replace(temporary,archive)
print('PASS finalized E08 report and source bundle')
