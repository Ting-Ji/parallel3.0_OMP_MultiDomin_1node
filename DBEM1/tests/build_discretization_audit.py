"""Build isolated source snapshots for discretization/integration audit (Windows Python)."""
from pathlib import Path
import argparse, json, re, shutil, subprocess, xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[2]
WORK=ROOT/'validation_runs/discretization_audit_20261019'
EXTRA_DEFINES=[]
NS={'ms':'http://schemas.microsoft.com/developer/msbuild/2003'}
KERNEL=r"""
static int RunDiscretizationKernelProbe() {
    FILE* f=0; fopen_s(&f,"kernel_probe.csv","w");
    if(!f) return 2;
    fprintf(f,"case,nu,dt,rx,ry,rz,n,k,U,T1,T2,I,J\n");
    DSquareElement element;
    for(int q=0;q<12;++q) {
        double nu=(q%3==0 ? 0.0 : q%3==1 ? 0.25 : 0.45);
        double dt=(q%2 ? 0.0005 : 0.002);
        DSquareElement::m_DMat=BuildDynaMat(nu,20000.0/(2*(1+nu)),8.0,dt);
        double v[3]={.043+.007*q,.061-.002*q,.089+.003*q};
        double r=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        double R[6]={0,r,r*r,r*r*r,r*r*r*r,r*r*r*r*r};
        double RI[4]={0,v[0],v[1],v[2]};
        double N[3]={.3,.4,sqrt(.75)};
        int maxlag=(int)ceil(r/(DSquareElement::m_DMat.C2*dt))+2;
        for(int n=0;n<=maxlag;++n) {
            double U[9],T1[9],T2[9],is,js;
            element.GetDynaUij(U,n,dt,R,RI);
            element.GetDynaT1ij(T1,n,dt,R,RI,N);
            element.GetDynaT2ij(T2,n,dt,R,RI,N);
            element.ISeries(n,dt,R,is); element.JSeries(n,dt,R,js);
            for(int k=0;k<9;++k)
                fprintf(f,"%d,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    q,nu,dt,v[0],v[1],v[2],n,k,U[k],T1[k],T2[k],is,js);
        }
    }
    fclose(f); return 0;
}
"""
DUMP=r"""
    if (getenv("DBEM_AUDIT_DUMP")) {
        auto dumpMatrix = [](const char* name, const auto& matrix) {
            FILE* f=0; const std::string path=DBEMOutputPath(name);
            if(fopen_s(&f,path.c_str(),"wb")!=0 || !f)
                throw std::runtime_error("audit matrix open failed");
            const bool ok=matrix.WriteBinary(f); fclose(f);
            if(!ok) throw std::runtime_error("audit matrix write failed");
        };
        dumpMatrix("audit_A.bin",A); dumpMatrix("audit_K.bin",KnownM);
        FILE* f=0; fopen_s(&f,DBEMOutputPath("audit_history.csv").c_str(),"w");
        fprintf(f,"lag,symmetric\n");
        for(long lag=1;lag<=MaxN;++lag) {
            char name[64]; sprintf_s(name,"audit_T_%04ld.bin",lag);
            dumpMatrix(name,MTS[lag]);
            sprintf_s(name,"audit_G_%04ld.bin",lag);
            if(MGSUseSym[lag]) dumpMatrix(name,MGS[lag]);
            else dumpMatrix(name,MGSFull[lag]);
            fprintf(f,"%ld,%d\n",lag,(int)MGSUseSym[lag]);
        }
        fclose(f);
        fopen_s(&f,DBEMOutputPath("audit_geometry.csv").c_str(),"w");
        fprintf(f,"element,node,domain,bcid");
        for(int k=0;k<9;++k) fprintf(f,",t%d",k);
        fprintf(f,"\n");
        for(long e=0;e<model.elementCount;++e) {
            const long node=elements[e].m_nodeID[0];
            fprintf(f,"%ld,%ld,%d,%d",e,node,elements[e].DomainID,elements[e].BCID);
            for(int k=0;k<9;++k) fprintf(f,",%.17g",DSquareElement::m_transmat[node][k]);
            fprintf(f,"\n");
        }
        fclose(f);
        if(getenv("DBEM_AUDIT_DUMP_ONLY")) return 1;
    }
"""
def build(name):
    changes={'base':{},'pw12':{'NUMPW':12},'pw24':{'NUMPW':24},
             'sin40':{'SINNUMPW':40},'accurate':{'NUMPW':12,'SINNUMPW':40,'GAUSSPOINT':8,'SINGAUSSPOINT':8}}.get(name, {'NUMPW':int(name[2:])} if re.fullmatch(r'pw[1-9][0-9]*',name) else None)
    if re.fullmatch(r'g[1-9][0-9]*',name):
        changes={'NUMPW':6,'GAUSSPOINTPW':int(name[1:])}
        assert 1<=changes['GAUSSPOINTPW']<=40,'Stored Gauss rules support at most 40 points per axis'
    assert changes is not None,name
    folder=WORK/name;src=folder/'src';src.mkdir(parents=True,exist_ok=True)
    for p in (ROOT/'DBEM1').iterdir():
        if p.is_file() and p.suffix.lower() in ('.cpp','.h','.hpp','.inl'):
            shutil.copy2(p,src/p.name)
    p=src/'precompiler.h';data=p.read_bytes()
    for key,value in changes.items():
        data,n=re.subn(rb'(const int '+key.encode()+rb'\s*=\s*)\d+',lambda m:m[1]+str(value).encode(),data)
        assert n==1,(key,n)
    p.write_bytes(data)
    p=src/'DBEM.cpp';data=p.read_bytes()
    needle=b'int main()\r\n{'
    assert data.count(needle)==1
    data=data.replace(needle,KERNEL.encode()+b'\r\n'+needle+b'\r\n if(getenv("DBEM_AUDIT_KERNEL_ONLY")) return RunDiscretizationKernelProbe();')
    p.write_bytes(data)
    p=src/'multidomain.cpp';data=p.read_bytes()
    needle=b'\tassemblyPool.reset();'
    assert data.count(needle)==1
    p.write_bytes(data.replace(needle,needle+DUMP.encode()))
    ET.register_namespace('',NS['ms'])
    tree=ET.parse(ROOT/'DBEM1/DBEM1.vcxproj')
    project=tree.getroot()
    for node in project.findall('.//ms:ItemGroup/ms:ClCompile',NS):
        node.set('Include',str(src/node.attrib['Include']))
    for group in project.findall('.//ms:ItemDefinitionGroup',NS):
        for node in group.findall('.//ms:AdditionalIncludeDirectories',NS):
            node.text=str(src)+';'+str(ROOT/'include')+';%(AdditionalIncludeDirectories)'
        for node in group.findall('.//ms:AdditionalLibraryDirectories',NS):
            node.text=str(ROOT/'lib/x64')+';%(AdditionalLibraryDirectories)'
        for node in group.findall('ms:ClCompile',NS):
            defs=node.find('ms:PreprocessorDefinitions',NS)
            if defs is None:defs=ET.SubElement(node,'{'+NS['ms']+'}PreprocessorDefinitions')
            defs.text=';'.join(EXTRA_DEFINES)+';_CRT_SECURE_NO_WARNINGS;DBEM_POLICY_VALIDATION_OUTPUT=1;DBEM_POLICY_ASSEMBLY_CHECK=0;DBEM_POLICY_SOLVE_DIAGNOSTICS=1;DBEM_POLICY_TECPLOT_OUTPUT=0;'+(defs.text or '')
    project_file=folder/(name+'.vcxproj')
    tree.write(project_file,encoding='utf-8',xml_declaration=True)
    for helper in ('gaussposition','gaussweight'):shutil.copy2(ROOT/'DBEM1'/helper,folder/helper)
    ms='C:/Program Files/Microsoft Visual Studio/2022/Community/MSBuild/Current/Bin/MSBuild.exe'
    done=subprocess.run([ms,str(project_file),'/p:Configuration=Release','/p:Platform=x64','/m:1','/v:quiet','/clp:ErrorsOnly'],cwd=folder,capture_output=True)
    (folder/'build.log').write_bytes(done.stdout+done.stderr)
    if done.returncode:raise RuntimeError((done.stdout+done.stderr).decode(errors='replace'))
    (folder/'changes.json').write_text(json.dumps(changes,indent=2))
    print(folder/'x64/Release'/(name+'.exe'),flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('variants',nargs='*',default=['base','pw12','sin40'])
    for name in parser.parse_args().variants:build(name)
