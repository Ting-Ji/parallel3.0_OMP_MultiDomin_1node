#include "../DBEM1/pre_gmres.h"
#include "../DBEM1/gpu_gmres_bridge.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>

static double seconds() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
static void run(int blocks,bool dense) {
    const int n=3*blocks;
    CCSRMat a; a.initial(blocks,dense?blocks*blocks:3*blocks-2);
    int at=0; a.row_ptr[0]=0;
    for(int i=0;i<blocks;++i) {
        int begin=dense?0:(std::max)(0,i-1),end=dense?blocks:(std::min)(blocks,i+2);
        for(int j=begin;j<end;++j) {
            a.col_index[at]=j;
            for(int r=0;r<3;++r) for(int c=0;c<3;++c) {
                double value=i==j?(r==c?4.+r*.1:.05):
                    (dense?.5/std::sqrt(double(n)):-.8)*std::sin((3*i+r+1)*12.9898+(3*j+c+1)*78.233);
                a.value[size_t(at)*9+r*3+c]=value;
            }
            ++at;
        }
        a.row_ptr[i+1]=at;
    }
    PreConditioner p;p.m_LeafCount=blocks;
    p.m_LeafPointCount=new long[blocks];p.m_LeafBeginID=new long[blocks];p.m_RePID=new long[blocks];p.IID=new long[blocks];
    p.m_InputPID=new long[blocks];p.m_OutputPID=new long[blocks];p.m_PreM=new Wmatrix[blocks];
    for(int i=0;i<blocks;++i) {
        p.m_LeafPointCount[i]=1;p.m_LeafBeginID[i]=i;p.m_RePID[i]=i;p.IID[i]=3*i;
        // Deliberately non-identity maps test project gather/scatter handling.
        p.m_InputPID[i]=blocks-1-i;p.m_OutputPID[i]=blocks-1-i;
        p.m_PreM[i].initial(3,3);p.m_PreM[i]=0.;
        for(int c=1;c<=3;++c) p.m_PreM[i](c,c)=1/(4.+(c-1)*.1);
    }
    auto pa=DbemGpu::PackMatrix(a),pp=DbemGpu::PackPreconditioner(p,n);
    Wvector probe(n),ap(n),mp(n);
    for(int i=1;i<=n;++i)probe[i]=std::sin(i*.013);
    SparseMul(a,probe,ap);GMRES_M_X(p,probe,mp);
    for(int r=0;r<n;++r) {
        double av=0,pv=0;
        for(int k=pa.row[r];k<pa.row[r+1];++k)av+=pa.val[k]*probe.b[pa.col[k]];
        for(int k=pp.row[r];k<pp.row[r+1];++k)pv+=pp.val[k]*probe.b[pp.col[k]];
        if(std::abs(av-ap.b[r])>1e-11 || std::abs(pv-mp.b[r])>1e-12)throw std::runtime_error("CSR conversion changed A/P");
    }
    Wvector b(n),x0(n,0),xc(n),xg(n),known(n),check(n);
    double setup=seconds();_putenv_s("DBEM_GMRES_BACKEND","cuda");
    DbemGpu::Scope scope(a,p);setup=seconds()-setup;
    for(int rep=0;rep<4;++rep) {
        for(int i=1;i<=n;++i)known[i]=std::sin(i*.013+rep*.1);
        SparseMul(a,known,b);
        int cf=0,gf=0,ci[2]={},gi[2]={};
        _putenv_s("DBEM_GMRES_BACKEND","cpu");double cpu=seconds();
        gmres(nullptr,a,b,p,60,1e-10,4,x0,xc,cf,ci);cpu=seconds()-cpu;
        _putenv_s("DBEM_GMRES_BACKEND","cuda");double gpu=seconds();
        gmres(nullptr,a,b,p,60,1e-10,4,x0,xg,gf,gi);gpu=seconds()-gpu;
        SparseMul(a,xg,check);double e=0,bn=0,xd=0,xn=0;
        for(int i=1;i<=n;++i) {e+=(check[i]-b[i])*(check[i]-b[i]);bn+=b[i]*b[i];xd+=(xg[i]-xc[i])*(xg[i]-xc[i]);xn+=xc[i]*xc[i];}
        double residual=std::sqrt(e/bn),difference=std::sqrt(xd/xn);
        if(cf || gf || residual>1e-10 || difference>1e-8)throw std::runtime_error("CPU/GPU GMRES regression failed");
        printf("BENCH,{\"n\":%d,\"dense\":%s,\"repeat\":%d,\"cpu\":%.9g,\"gpu\":%.9g,\"setup\":%.9g,\"residual\":%.9g,\"solution_difference\":%.9g,\"cpu_inner\":%d,\"gpu_inner\":%d}\n",n,dense?"true":"false",rep,cpu,gpu,setup,residual,difference,ci[1],gi[1]);
    }
}
int main() {
    try {setvbuf(stdout,nullptr,_IONBF,0);thread_num=4;run(64,false);run(256,true);run(1024,true);return 0;}
    catch(const std::exception& e) {fprintf(stderr,"GPU benchmark failed: %s\n",e.what());return 2;}
}
