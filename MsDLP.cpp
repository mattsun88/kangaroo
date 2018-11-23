#include <math.h>
#include <NTL/ZZ.h>
#include <time.h>
#include <vector>
#include <chrono>
#include <fstream>

#include "randwalk.h"

#define P_BITLEN 160
#define L_xNUM 10//how many x 
#define LOOP 100

using namespace NTL;
using namespace std;

void set_initstate(ZZ& P,ZZ& p,ZZ& g,long p_bitlen);
void MYDLP(ZZ P,ZZ p,ZZ g,ZZ *x,ZZ *y,ZZ w,int &tmp);
Randwalk MakeTame(ZZ P,ZZ p,ZZ g,ZZ w,ZZ a[],ZZ e[],ZZ &Rete);

int main(){
  ZZ P,p,g,m_x[L_xNUM+1],y[L_xNUM+1];//y=g^x mod P, g^p=1 mod P
  ZZ x[L_xNUM+1];

  std::chrono::system_clock::time_point  start, end;
  ofstream ofs;

  set_initstate(P,p,g,P_BITLEN);

  ZZ seed;
  time_t t = time(NULL);
  seed = (long)t;
  SetSeed(seed);

  ZZ w;//x range(0<=x<w)
  w = to_ZZ(pow(2,20));//2^20  
  start = std::chrono::system_clock::now();
  int fail=0;
  int count=0;
  int sum=0;
  int tmp=0;
  int failnum=0;
  for(;count<LOOP;count++){
    puts("---------------------------------------------------");
    for(int i=1;i<=L_xNUM;i++){
      x[i] = RandomBnd(w-1)+1;
      y[i] = PowerMod(g,x[i],P);
      m_x[i] = to_ZZ(0);
    }

    MYDLP(P,p,g,m_x,y,w,tmp);
    sum +=tmp;

    for(int i=1;i<=L_xNUM;i++){
      if(x[i]!=m_x[i]){
	cout<<"fail"<<endl;
	failnum+=1;
      }
    }
  }
  end = std::chrono::system_clock::now();
  cout<<"avg count="<<(double)sum/(double)(LOOP * L_xNUM)<<endl;
  cout<<"fail="<<failnum<<endl;
  //cout<<"P="<<P<<endl;
  //cout<<"p="<<p<<endl;
  //cout<<"g="<<g<<endl;
  //cout<<"w="<<w<<endl;
  
  double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
  cout<<"avg time="<<elapsed/LOOP<<"milli sec"<<endl;
  
  ofs.open( "10.txt" , ios::out | ios::app);
  ofs << elapsed/LOOP<< endl;
  ofs.close();
  

}

void set_initstate(ZZ& P,ZZ& p,ZZ& g,long p_bitlen){
  ZZ s;
  s = to_ZZ((unsigned int)time(NULL)) ;
  SetSeed(s);

  p = GenGermainPrime_ZZ(p_bitlen,80);
  P = 2*p+1;
  do{g = RandomBnd(P);}
  while(PowerMod(g,2,P)==1 || PowerMod(g,p,P)==1);
  g = PowerMod(g,2,P);
}

void MYDLP(ZZ P,ZZ p,ZZ g,ZZ *x,ZZ *y,ZZ w,int &tmp){

  ZZ m,N;
  m = sqrt(L_xNUM*to_long(w));
  N = m/(2*L_xNUM);

  ZZ a[L_xNUM+1];
  ZZ b[L_xNUM+1],b0[L_xNUM+1];
  ZZ d[L_xNUM+1],e[L_xNUM+1],en,u[L_xNUM+1];
  ZZ max_e;
  Randwalk walk(g,m,P);

  //tame系列を作成しランダムウォーク関数を受け取る
  walk = MakeTame(P,p,g,w,a,e,max_e);

  cout<<"wild"<<endl;
  //wild
  int count=0;
  int loop=0;
  for(int k=1;k<=L_xNUM;k++){
    //(W1)
    b[k]=y[k];
    b0[k]=b[k];
    u[k]=0;
    d[k]=u[k];
    count=0;
  rset2: 
    //(W2)
    if(count>5&&loop==0){ //5回失敗が続いたらtame系列を作り直す
      walk=MakeTame(P,p,g,w,a,e,max_e);
      loop++;
      puts("OneMore!");
    }
    for(int i=1;;i++){
      //(a)
      d[k]=AddMod(d[k],walk.get_r(b[k]),p);
      b[k]=walk.get_R(b[k]);
      
      //(b)
      for(int h=1;h<=L_xNUM;h++){
	if(b[k]==a[h]){
	  x[k]= SubMod(e[h],d[k],p);
	  /*
	    cout<<"a = "<<a[h][to_int(N)]<<endl;
	    cout<<"b = "<<b[k][i]<<endl;
	    cout<<"e = "<<e[h]<<endl;
	    cout<<"d = "<<d[k]<<endl;
	    cout<<"h = "<<h<<endl;
	  */
	  tmp+=count;
	}
      }
      if(x[k]!=0)
	break;
      //(c)
      if(d[k]>max_e){ //Tame系列の値を通り越したら失敗。値を1増やしてやり直し
	count++;
	b[k] = MulMod(b0[k],g,P);
	b0[k] = b[k];
	u[k]=u[k]+1;
	d[k]=u[k];
	goto rset2;
      }
      if(count>10){ //10回失敗したら探索打ち切りで失敗を出力
	cout<<"count>10..."<<endl;
	x[k]=to_ZZ(0);
	tmp+=10;
	break;
      }
    }
  }
}

Randwalk MakeTame(ZZ P,ZZ p,ZZ g,ZZ w,ZZ a[],ZZ e[],ZZ &Rete){

  ZZ m,N;
  m = sqrt(L_xNUM*to_long(w));
  N = m/(2*L_xNUM);

  ZZ an ,a0[L_xNUM+1];
  ZZ en,f[L_xNUM+1];
  ZZ max_e;
  Randwalk walk(g,m,P);
  cout<<"Make Tames..."<<endl;
  //walk.R_print();
  a[1]=PowerMod(g,w,P);
  a0[1]=a[1];
  e[1]=w;
  f[1]=w; 
  //(0<=i<N-1)
  for(int i=0;i<N;i++){
    e[1] = AddMod(e[1],walk.get_r(a[1]),p);
    a[1]=walk.get_R(a[1]);
  }
  max_e = e[1];
  //(2<=k<=L)
  for(int k=2;k<=L_xNUM;k++){
    //T1
    a[k]=MulMod(a0[k-1],g,P);
    a0[k]=a[k];
    f[k]=f[k-1]+1;
    e[k]=f[k];
  rset1: 
    //T2
    for(int i=0;i<N;i++){
      //(a)
      e[k]=AddMod(e[k],walk.get_r(a[k]),p);
      a[k]=walk.get_R(a[k]);
      //(b)
      for(int j=1;j<k;j++){
	if(a[k]==a[j]){
	  a[k] = MulMod(a0[k],g,P);
	  a0[k]= a[k];
	  f[k] = f[k]+1;
	  e[k] = f[k];
	  goto rset1;
	}
      }
      
    }
    //T3
    an = a[k];
    en = e[k];
    for(int i=0;i<N;i++){
      if(en >max_e ){break;}
      //(a+)
      en=AddMod(en,walk.get_r(an),p);
      an=walk.get_R(an);
      //(b+)
      for(int j=1;j<k;j++){
	if(an==a[j]){
	  a[k] = MulMod(a0[k],g,P);
	  a0[k]= a[k];
	  f[k] = f[k]+1;
	  e[k] = f[k];
	  goto rset1;
	}
      }

    }
    if(max_e<e[k]){max_e=e[k];}
  }
  Rete = max_e;
  return walk;
}
