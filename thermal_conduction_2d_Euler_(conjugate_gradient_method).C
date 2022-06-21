#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include "nrutil.h"
#include <stdio.h>
FILE  *  initgr(int i,int j,int k,int l,char * s);
void closegr(FILE  *g1);
void linbcg(unsigned long n, double b[], double x[], int itol, double tol,
	    int itmax, int *iter, double *err);

using namespace std;

const double L=1.,dt=.005,tmax=.1;//dt=.00005
const int N=10;
const double gy3=sqrt(3.),t1=0,t2=0,t3=100,T0=200,T2=100,T3=30,Q=10,D=10,A=10;
int NN;

double sa[4*N*N];
unsigned long ija[4*N*N];

double xx(int i,int j){
    double x;
    x=L/N/2.*(i+j-1.);
    return x;
}

double yy(int i,int j){
    double y;
    y=L/N*gy3/2.*(j-1.+1./3.*(1.+((i+1)%2)));
    return y;
}

int nn(int i,int j){
    int n;
    n=i+(2*N+1-j)*(j-1);
    return n;
}


void ij(int n, int *i, int *j){
  double x;
  x=(double)(N+1);
  x=x-sqrt(x*x-(double) n+2.-2.*x);
  *j=(int)x;
  *i=n-(2*N+1-*j)*(*j-1);
}

void szomszedok(int n, int *nsz){

    for(int i=1;i<=3;i++) nsz[i]=0;    
    int i,j;
    ij(n,&i,&j);


    if(i%2==1){

	if(j>1) nsz[1]=nn(i+1,j-1);
	if(i>1) nsz[2]=n-1;
	if(i<2*N+1-2*j) nsz[3]=n+1;
    }else{
	nsz[3]=nn(i-1,j+1);
	nsz[1]=n-1;
	nsz[2]=n+1;
    } 
}

void ritka_matrix_indexelese(){
    int nsz[4];
    unsigned long ii=NN+1;
    for(int n=1;n<=NN;n++){
	szomszedok(n,nsz);
	ija[n]=0;
	for(int j=1;j<=3;j++){
	    if(nsz[j]>0){
		ii++;
		if(ija[n]==0) ija[n]=ii;
		ija[ii]=nsz[j];
	    }
	}
    }
    ija[NN+1]=ii+1;
}



void showdomain(){
    FILE *g2; 
    char cim[40];
    sprintf(cim,"The domain");
    g2=initgr(600,600,700,0,cim);
    fprintf(g2,"set xrange[0.:%f]\n",L);
    fprintf(g2,"set yrange[0.:%f]\n",L);
    fprintf(g2,"set multiplot\n",L);
    for(int j=1;j<=N;j++){
	for(int i=1;i<=2*N+1-2*j;i++){
	    double x=xx(i,j);
	    double y=yy(i,j);
	    int n=nn(i,j);
	    fprintf(g2,"set label \"%d\" at %21.6f, %21.6f center tc lt 1\n",n,x,y);
	}
    }
    fprintf(g2,"plot '-' w l lc 3\n");
    for(int j=1;j<N;j++){
	fprintf(g2,"%f21.6 %f21.6\n",L/(double)N*.5*(double)j,(double)j*L/(double)N*gy3/2.);
	fprintf(g2,"%f21.6 %f21.6\n",L-L/(double)N*.5*(double)j,(double)j*L/(double)N*gy3/2.);
	fprintf(g2,"\n");
    }     
   for(int j=0;j<N;j++){
       fprintf(g2,"%f21.6 %f21.6\n",L/(double)N*(double)j,0.);
       fprintf(g2,"%f21.6 %f21.6\n",.5*L+.5*L/(double)N*(double)j,gy3*.5*L-gy3*.5*L/(double)N*(double)j);
       fprintf(g2,"\n");   
       fprintf(g2,"%f21.6 %f21.6\n",L-L/(double)N*(double)j,0.);
       fprintf(g2,"%f21.6 %f21.6\n",.5*L-.5*L/(double)N*(double)j,gy3*.5*L-gy3*.5*L/(double)N*(double)j);
       fprintf(g2,"\n"); 
   }
    fprintf(g2,"e \n");
    closegr(g2);

}

void kezdeti_feltetelek(double * T){
    double x1=0,y1=0,x2=L,y2=0,x3=.5*L,y3=.5*gy3*L,x,y,a,b,u;
    a=(y2-y1)*(t3-t1)-(t2-t1)*(y3-y1);
    b=(t2-t1)*(x3-x1)-(x2-x1)*(t3-t1);
    u=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);
    for(int j=1;j<=N;j++){
	for(int i=1;i<=2*N+1-2*j;i++){
	    x=xx(i,j);
	    y=yy(i,j);
	    int n=nn(i,j);
	    T[n]=t1-a*(x-x1)/u-b*(y-y1)/u;
	}
    }
}



void ritka_matrix_elemei(){

    int n;
    double h=L/N;
    double d=4*D/h/h;


    for(int j=1;j<=N;j++){
	for(int i=1;i<=2*N+1-2*j;i++){
	    n=nn(i,j); 
    if(i>1&& i < 2*N+1-2*j&&(j>1||i%2==0)){
	    sa[n]=1+3*d*dt;
	    for(int k=0;k<3;k++) sa[ija[n]+k]=-d*dt; 
	    /*if(i%2==0){
		n1=nn(i-1,j+1);
	    }else{
		n1=nn(i+1,j-1);
	    }
	    dT[n]=d*(T[n+1]+T[n-1]+T[n1]-3*T[n]);*/
	}
    if(i==1&&j>1&&j<N){
	sa[n]=1+dt*d*(2.+4./(1+3*gy3*D/A/h));
	sa[ija[n]]=-dt*d*(1.+.5/(1+3*gy3*D/A/h));
	sa[ija[n]+1]=-dt*d*(1.+.5/(1+3*gy3*D/A/h));
/*	n1=nn(i+1,j-1);
	dT[n]=d*(T[n+1]+T[n1]-2*T[n])+2*d/(1+3*gy3*D/A/h)*(3./2.*T0-2.*T[n]+1./4.*(T[n+1]+T[n1]));*/
    }

    if(i==2*N+1-2*j&&j>1&&j<N){
	sa[n]=1+dt*d*6.;
	sa[ija[n]]=-dt*d*1.5;
	sa[ija[n]+1]=-dt*d*1.5;	
/*	n1=nn(i+1,j-1);
	dT[n]=d*(T[n-1]+T[n1]-2*T[n])+d*2*(1./4.*(T[n-1]+T[n1])-2*T[n]+3./2.*(T2*(N-j+.5)/double(N)+T3*(j-.5)/double(N)));*/
    }  

    if(i>1&&i<2*N+1-2*j&&j==1&&i%2==1){
	sa[n]=1+dt*d*2.;
	sa[ija[n]]=-dt*d;
	sa[ija[n]+1]=-dt*d;	
	//dT[n]=d*(T[n+1]+T[n-1]-2*T[n])+Q/h*4./gy3;
    }
    if(i==1&&j==1){
	sa[n]=1+dt*d*(1.+2./(1+2*gy3*D/A/h));
	sa[ija[n]]=-dt*d;
	//dT[n]=d*(T[n+1]-T[n])+Q/h*4./gy3+2*d/(1+2*gy3*D/A/h)*(T0-T[n]);
    }

    if(i==2*N+1-2*j&&j==1){
	sa[n]=1+dt*d*3.;
	sa[ija[n]]=-dt*d;
	//dT[n]=d*(T[n-1]-T[n])+Q/h*4./gy3+2*d*(T2*(N-j+.5)/double(N)+T3*(j-.5)/double(N)-T[n]);
    }  

  if(i==1&&j==N){
      sa[n]=1+dt*d*(3.+2./(1+2*gy3*D/A/h));
      sa[ija[n]]=-dt*d;
      /*n1=nn(i+1,j-1);
      dT[n]=d*(T[n1]-T[n])+2*d/(1+2*gy3*D/A/h)*(T0-T[n]);
      dT[n]+=2*d*(T2*(N-j+.5)/double(N)+T3*(j-.5)/double(N)-T[n]);*/
  }

	}
    }

}

void implicit_Euler(double * T){

  double *b,*x;
  b=dvector(1,NN);
  x=dvector(1,NN);

  int n;
  double h=L/N;
  double d=4*D/h/h;
    for(int j=1;j<=N;j++){
	for(int i=1;i<=2*N+1-2*j;i++){
	    n=nn(i,j); 
	if(i>1&& i < 2*N+1-2*j&&(j>1||i%2==0)){
	    b[n]=T[n];
	    /*if(i%2==0){
		n1=nn(i-1,j+1);
	    }else{
		n1=nn(i+1,j-1);
	    }
	    dT[n]=d*(T[n+1]+T[n-1]+T[n1]-3*T[n]);*/
	}
    if(i==1&&j>1&&j<N){
	b[n]=T[n]+dt*2*d/(1+3*gy3*D/A/h)*(3./2.*T0);
/*	n1=nn(i+1,j-1);
	dT[n]=d*(T[n+1]+T[n1]-2*T[n])+2*d/(1+3*gy3*D/A/h)*(3./2.*T0-2.*T[n]+1./4.*(T[n+1]+T[n1]));*/
    }

    if(i==2*N+1-2*j&&j>1&&j<N){
	b[n]=T[n]+dt*d*2*(3./2.*(T2*(N-j+.5)/double(N)+T3*(j-.5)/double(N)));
/*	n1=nn(i+1,j-1);
	dT[n]=d*(T[n-1]+T[n1]-2*T[n])+d*2*(1./4.*(T[n-1]+T[n1])-2*T[n]+3./2.*(T2*(N-j+.5)/double(N)+T3*(j-.5)/double(N)));*/
    }  

    if(i>1&&i<2*N+1-2*j&&j==1&&i%2==1){
	b[n]=T[n]+dt*Q/h*4./gy3;
		
	//dT[n]=d*(T[n+1]+T[n-1]-2*T[n])+Q/h*4./gy3;
    }
    if(i==1&&j==1){
	b[n]=T[n]+dt*(Q/h*4./gy3+2*d/(1+2*gy3*D/A/h)*T0);
	//dT[n]=d*(T[n+1]-T[n])+Q/h*4./gy3+2*d/(1+2*gy3*D/A/h)*(T0-T[n]);
    }

    if(i==2*N+1-2*j&&j==1){
	b[n]=T[n]+dt*(Q/h*4./gy3+2*d*(T2*(N-j+.5)/double(N)+T3*(j-.5)/double(N)));
	//dT[n]=d*(T[n-1]-T[n])+Q/h*4./gy3+2*d*(T2*(N-j+.5)/double(N)+T3*(j-.5)/double(N)-T[n]);
    }  

  if(i==1&&j==N){
      b[n]=T[n]+dt*(2*d/(1+2*gy3*D/A/h)*T0+2*d*(T2*(N-j+.5)/double(N)+T3*(j-.5)/double(N)));
      /*n1=nn(i+1,j-1);
      dT[n]=d*(T[n1]-T[n])+2*d/(1+2*gy3*D/A/h)*(T0-T[n]);
      dT[n]+=2*d*(T2*(N-j+.5)/double(N)+T3*(j-.5)/double(N)-T[n]);*/
  }

	}
    }


  int iter;
  double err;
  linbcg(NN, b, x, 1, 1.e-6, 10000, &iter, &err);

  for(int i=1;i<=NN;i++) T[i]=x[i];

  free_dvector(b,1,NN);
  free_dvector(x,1,NN);

}


int main(){
    NN=nn(1,N); 
    double T[NN+1];

    showdomain();

/*    int ix,iy;

    int nsz[4];
    for(int i=1;i<=NN;i++){
	ij(i,&ix,&iy);
	cout<<setfill(' ')<<setw(5)<<i;
	cout<<setfill(' ')<<setw(5)<<ix;
	cout<<setfill(' ')<<setw(5)<<iy<<"|";
	szomszedok(i,nsz);
	for(int j=1;j<=3;j++) cout<<setfill(' ')<<setw(5)<<nsz[j];
	cout<<endl;
    }
    exit(0);*/

    FILE * graphics,*g2; 
    char cim[40];
    sprintf(cim,"Heat diffusion on a triangular domain");
    graphics=initgr(600,600,0,0,cim);
    fprintf(graphics,"set xrange[0.:%f]\n",L);
    fprintf(graphics,"set yrange[0.:%f]\n",L);
   
    
    kezdeti_feltetelek(T);

    ritka_matrix_indexelese();
    ritka_matrix_elemei();

    for(double t=0;t<tmax;t+=dt){
	implicit_Euler(T);
	fprintf(graphics,"splot '-' u 1:2:3 w l lc 2\n");
	double x,y;
	for(int j=1;j<=N;j++){
	    for(int i=1;i<=2*N+1-2*j;i++){
		x=xx(i,j);
		y=yy(i,j);
		int n=nn(i,j);
		fprintf(graphics,"%f21.6 %f21.6 %f21.6\n",x,y,T[n]);
	    }
	    for(int i=2*N+2-2*j;i<=2*N;i++){ 
		x+=L/N/100.;
		fprintf(graphics,"%f21.6 %f21.6 %f21.6\n",x,y,(T2*(N-j+.5)/double(N)+T3*(j-.5)/double(N)));
	    }
	    fprintf(graphics," \n");
	}
	fprintf(graphics,"e \n");
    }
    cout<<"The end.\n"<<endl;
    int k;
    cin>>k;
    closegr(graphics);

}
