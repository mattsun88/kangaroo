#include <math.h>
#include <NTL/ZZ.h>
#include <time.h>
#include <vector>
#include <chrono>
#include <fstream>
#include "randwalk.h"
#include "HashTable.h"

#define P_BITLEN 160
//#define L_xNUM 10//how many x 
#define LOOP 100
#define LIMIT 20

int L_xNUM =  1;
using namespace NTL;
using namespace std;

void set_initstate(ZZ& P,ZZ& p,ZZ& g,long p_bitlen);
void MYDLP(ZZ P,ZZ p,ZZ g,ZZ *x,ZZ *y,ZZ w,int &tmp);
Randwalk MakeTame(ZZ P,ZZ p,ZZ g,ZZ w,HashTable &tame,ZZ &Rete);
Randwalk MakeTame2(ZZ P,ZZ p,ZZ g,ZZ w,ZZ *a,ZZ *e,ZZ &Rete);
void TameCheck(ZZ P,ZZ p,ZZ g,ZZ w,ZZ *a,ZZ *e);
int main(){

  while(L_xNUM<=20){
    cout<<"L="<<L_xNUM<<endl;
    ZZ P,p,g,m_x[L_xNUM+1],y[L_xNUM+1];//y=g^x mod P, g^p=1 mod P
    ZZ x[L_xNUM+1];
    //HashTable htable;
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
      //puts("---------------------------------------------------");

      //HashTable htable;
      for(int i=1;i<=L_xNUM;i++){
	x[i] = RandomBnd(w-1)+1;//xを格納
	y[i] = PowerMod(g,x[i],P);
	m_x[i] = to_ZZ(0);//アルゴリズムの回答が返ってくる配列
      }

      MYDLP(P,p,g,m_x,y,w,tmp);
      sum +=tmp;

      for(int i=1;i<=L_xNUM;i++){//正誤チェック
	if(x[i]!=m_x[i]){
	  if(m_x[i] != to_ZZ(0)){//アルゴリズムの回答が間違っている場合は確認のため出力
	  
	    cout<<x[i]<<endl;
	    cout<<m_x[i]<<endl;
	    puts("");
	  
	  }
	  //cout<<"fail"<<endl;
	  failnum+=1;
	}
      }
    }
    end = std::chrono::system_clock::now();
    //cout<<"avg count="<<(double)sum/(double)(LOOP * L_xNUM)<<endl;
    cout<<"fail="<<failnum<<endl;
    //cout<<"P="<<P<<endl;
    //cout<<"p="<<p<<endl;
    //cout<<"g="<<g<<endl;
    //cout<<"w="<<w<<endl;
  
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    cout<<"avg time="<<elapsed/LOOP<<"milli sec"<<endl;
    ofs.open( "2_19_160_1001_f.txt" , ios::out | ios::app);
    ofs << failnum<< endl;
    ofs.close();  
    ofs.open( "2_19_160_1001.txt" , ios::out | ios::app);
    ofs << elapsed/LOOP<< endl;
    ofs.close();
    L_xNUM++;
  }
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
  N = m/2*(L_xNUM);

  ZZ a[L_xNUM+1];
  ZZ b[L_xNUM+1],b0[L_xNUM+1];
  ZZ d[L_xNUM+1],e[L_xNUM+1],en,u[L_xNUM+1];
  ZZ max_e;
  Randwalk walk(g,m,P);
  HashTable tame;
  vector<vector<ZZ>> TmpTames;
  //tame系列を作成しランダムウォーク関数を受け取る
  walk = MakeTame(P,p,g,w,tame,max_e);
  /*
  if(L_xNUM > 5){//Lが多い時はハッシュテーブルに格納
    walk = MakeTame(P,p,g,w,tame,max_e);
  }
  else{
    walk = MakeTame2(P,p,g,w,a,e,max_e);
  }
  */
  //walk = MakeTame2(P,p,g,w,a,e,max_e);
  //TameCheck(P,p,g,w,a,e);
  //puts("");
  //cout<<"wild"<<endl;
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
    //(W2
    for(int i=1;;i++){
      //(a)
      d[k]=AddMod(d[k],walk.get_r(b[k]),p);
      b[k]=walk.get_R(b[k]);
      
      //(b)
      //if(L_xNUM > 0 ){
	ZZ tmpnum;
	tame.Search(b[k],tmpnum);
	if(tmpnum!=0){
	  x[k]= SubMod(tmpnum,d[k],p);
	  tmp+=count;	  
	}
      //}
      /*
      else{
	for(int h=1;h<=L_xNUM;h++){
	  if(b[k]==a[h]){
	    x[k]= SubMod(e[h],d[k],p);
	    tmp+=count;
	  }
	}
      }
      */
      if(x[k]!=0){
	for(int i = 0; i < count;i++){
	  //xがわかればそれ以前の失敗した系列の指数部がわかるのでTameカンガルーとして追加
	  /*
	    cout<<PowerMod(g,TmpTames[i][1] +  x[k],P)<<endl;
	    cout<<TmpTames[i][0]<<endl;
	  */
	  tame.Add( TmpTames[i][0] ,AddMod(TmpTames[i][1],x[k],p) );
	}
	break;
      }
      //(c)
      if(d[k]>max_e){ //Tame系列の値を通り越したら失敗。値を1増やしてやり直し
	//cout<<"here"<<endl;
	count++;

	//そのときの値と指数部を保存しておく
	vector<ZZ> tmp(2); 
	tmp[0]=b[k];tmp[1]=d[k];
	TmpTames.push_back(tmp);

	b[k] = MulMod(b0[k],g,P);
	b0[k] = b[k];
	u[k]=u[k]+1;
	d[k]=u[k];
	goto rset2;
      }
      if(count>LIMIT ){ //10回失敗したら探索打ち切りで失敗を出力
	//cout<<"count>20..."<<endl;
	x[k]=to_ZZ(0);
	//tmp+=10;
	break;
      }
    }
  }
}

Randwalk MakeTame(ZZ P,ZZ p,ZZ g,ZZ w,HashTable &tame,ZZ &Rete){

  ZZ m,N;
  m = sqrt(L_xNUM*to_long(w));
  N = m/(2*L_xNUM);
  
  ZZ a[L_xNUM+1],e[L_xNUM+1];
  ZZ an ,a0[L_xNUM+1];
  ZZ en,f[L_xNUM+1];
  ZZ max_e;
  Randwalk walk(g,m,P);
  // cout<<"Make Tames..."<<endl;
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
  tame.Add(a[1],e[1]);
  max_e = e[1];
  //(2<=k<=L)
  for(int k=2;k<=L_xNUM;k++){
    //T1
    a[k]=PowerMod(g,w+((k-1)*20),P);
    a0[k]=a[k];
    f[k]=w+((k-1)*20);
    e[k]=f[k];
  rset1: 
    //T2
    for(int i=0;i<N;i++){
      //(a)
      e[k]=AddMod(e[k],walk.get_r(a[k]),p);
      a[k]=walk.get_R(a[k]);
      //(b)
      ZZ tmp;
      tame.Search(a[k],tmp);
      if(tmp!=0){
	  a[k] = MulMod(a0[k],g,P);
	  a0[k]= a[k];
	  f[k] = f[k]+1;
	  e[k] = f[k];
	  //cout<<"here"<<endl;
	  goto rset1;
      }
      
    }

    
    //T3
    an = a[k];
    en = e[k];
    for(int i=0;i<N;i++){
      if(en >max_e ){break;}
      //cout<<"--------------------------------------------------------------------------------------------"<<endl;
      //(a+)
      en=AddMod(en,walk.get_r(an),p);
      an=walk.get_R(an);
      //(b+)
      ZZ tmp;
      tame.Search(an,tmp);
      if(tmp!=0){
	  a[k] = MulMod(a0[k],g,P);
	  a0[k]= a[k];
	  f[k] = f[k]+1;
	  e[k] = f[k];
	  //cout<<"here?"<<endl;
	  goto rset1;
      }

    }
    //cout<<a[k]<<endl;
    //cout<<e[k]<<endl;
    tame.Add(a[k],e[k]);
    if(max_e<e[k]){max_e=e[k];}
  }
  Rete = max_e;
  return walk;
}

Randwalk MakeTame2(ZZ P,ZZ p,ZZ g,ZZ w,ZZ a[],ZZ e[],ZZ &Rete){

  ZZ m,N;
  m = sqrt(L_xNUM*to_long(w));
  N = m/(2*L_xNUM);

  ZZ an ,a0[L_xNUM+1];
  ZZ en,f[L_xNUM+1];
  ZZ max_e;
  Randwalk walk(g,m,P);
  //cout<<"Make Tames..."<<endl;
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
    a[k]=PowerMod(g,w+((k-1)*20),P);
    a0[k]=a[k];
    f[k]=w+((k-1)*20);
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
	//cout<<"-----------------------------------------------------------"<<endl;
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

void TameCheck(ZZ P,ZZ p,ZZ g,ZZ w,ZZ a[],ZZ e[]){
  //ZZ sorted_a[L_xNUM+1],sorted_e[L_xNUM+1];
  ZZ tmpmin;
  ZZ tmpe,tmpa;
  int minnum;
  for(int i=1;i<=L_xNUM;i++){
    tmpmin=e[i];
    minnum=i;
    for(int k=i+1;k<=L_xNUM;k++){
      if(e[k]<tmpmin){
	minnum=k;
      }
    }
    tmpa=a[i];tmpe=e[i];
    a[i]=a[minnum];e[i]=e[minnum];
    a[minnum]=tmpa;e[minnum]=tmpe;
  }
  /*
  for(int i=1;i<=L_xNUM;i++){
    cout<<e[i]<<endl;
  }
  */
}
