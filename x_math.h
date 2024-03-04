#ifndef X_MATH_H
#define X_MATH_H

#define NAMEDIF abcdefghijklmnopqrstuvwxyz012345678901234567890

const double MY_PI = 4 * atan(1.0);
/*Enumeration of filter passband types*/
enum band_type{LOWPASS=1, HIGHPASS, BANDPASS, BANDSTOP};
/*Enumeration of window types to which FIR filters are added*/
enum win_type{RECTANG=1,BARTLETT,HANNING,HAMMING,BLACKMAN};
const double INFS=0.0;
const double INF=1.0/INFS;
const int MAXIT = 100;
const double EPSS = 6.0e-8;


class COMPLEX{
public:
	double r; //real part
	double i; //imaginary part

	COMPLEX(double r0=0,double i0=0);

	COMPLEX operator+(const COMPLEX& a) const;
	COMPLEX operator+(const double& a) const;
	void operator+=(const COMPLEX& a);
	void operator+=(const double& a);

	COMPLEX operator-(const COMPLEX& a) const;
	COMPLEX operator-(const double& a) const;
	void operator-=(const COMPLEX& a);
	void operator-=(const double& a);
	COMPLEX operator-();

	COMPLEX operator*(const COMPLEX& a) const;
	COMPLEX operator*(const double& a) const;
	void operator*=(const COMPLEX& a);
	void operator*=(const double& a);

	COMPLEX operator/(const COMPLEX& a) const;
	COMPLEX operator/(const double& a) const;	
	void operator/=(const COMPLEX& a);
	void operator/=(const double& a);
};

COMPLEX operator+(double a, const COMPLEX& b);
COMPLEX operator-(double a, const COMPLEX& b);
COMPLEX operator*(double a, const COMPLEX& b);
COMPLEX operator/(double a, const COMPLEX& b);

double abs(const COMPLEX& a);

/*The square root of (-INF) is inconsistent with matlab, it should be INFi, but the result here is INF+INFi*/
COMPLEX sqrt(const COMPLEX z);

/*Find the complex polynomial a[0]+a[1]s+a[2]s^2+... +a[n-1]s^(n-1) complex number. Similar to matlab polyval*/
COMPLEX polyval(COMPLEX a[],int n,COMPLEX s);
COMPLEX polyval(double a[],int n,COMPLEX s);
double polyval(double a[],int n,double s);

/*Polynomial multiplication is convolution, calculate y[L]=X[N]@h[M],L=M+ n-1,y[] without clearing 0, matlab*/
void conv(double y[],double x[],int N,double h[],int M);
void conv(COMPLEX y[],COMPLEX x[],int N,COMPLEX h[],int M);

/*The coefficients of the corresponding polynomial are obtained from a set of roots,
 where a[] is the output coefficient and b[] is the root of the input.
 a[0]+a[1]s+...+a[N]s^N=(s-b[0])(s-b[1])...(s-b[N-1])*/
void polycoef(COMPLEX a[],COMPLEX b[],int N);
void polycoef(double a[],double b[],int N);

/*For the n-1 complex polynomial a[0]+a[1]S+... Plus a[n-1]s to the n-1, take the m derivative*/
void polyder(COMPLEX a[],int n,int m);
/*For n-1 degree real polynomial a[0]+a[1]x+... Plus a[n-1]x to the n-1, take the m derivative*/
void polyder(double a[],int n,int m);

/*Discrete Fourier transform.
idata[] stores the result of the transformation, 
and data[] is the array to be transformed. Both arrays have length n.
isign=1 for the forward transformation and isign=-1 for the inverse transformation.*/
void dft(COMPLEX idata[],COMPLEX data[],int n,int isign);

/*Fast Fourier transform.
data[] is both an input array and an output array, and the array length n must be a power of two.
isign=1 for the forward transformation and isign=-1 for the inverse transformation.*/
void fft(COMPLEX data[],int n,int isign);

/*x=x@h.*/
double* cv(double *x,int N,double *h,int M);

/*here :a[]=a[]+b[]*/
double* pAdd(double* a,double* b,int M);

/*y[]=x0A<L-1>+x1A<L-2>B<1>+...+x[L-2]A<1>B<L-2>+x[L-1]B<L>
here:  k=(M>N)?M:N;
	   K=(K-1)*(L-1)+1;*/
double* pReplace(double *y,int K,
	double* x,int L,double* a,int M,double* b,int N);
int freqImage(double *e,int K,double *f,int L,
		 double *a,int M,double *b,int N,
		 double *c,int P,double *d,int Q);
void cpoly_coef(COMPLEX b[],int n,COMPLEX a[]);
void btwNormal(double d[],int N);
/*Calculate N of IIR filter*/
int  btwOrder(int bandType, double db1, double db2, double fs, double f1, double f2, ...);
int btwAf(double b[][2],
				  int order);
double btwC23(double c[2][3],
			   int bandType,int N,double db1,double fs,double f1,double f2,...);
void btwAf2Df(double H[][2][5],int L,
			  double b[][2],double c[2][3]);
double btw20lgHs(double f,double *p,int *N,double *fc);
double btw20lgHs(double f,double fc,double p[],int N);
double btw20lgHz(double f,double *fs,double H[][2][5],int *L);
double btw20lgHz(double f,double fs,double H[][2][5],int L);
void firWin(int winType,double w[],int N);
void firDesgin(double b[],int L,int bandType,int winType,double fs,double fL,...);
double firAbs(double f,double fs,
			   double b[],int L);
double firPhase(double f,double fs,
			   double b[],int L);
void polyLag(COMPLEX a[],int m,COMPLEX *x,double r_err,int polish);


/////////////���ֺ���/////////////////////////////////////

//func��[a,b]���仮�ֳ�2^N������������,����ֵ���Ϊy��������Ϊdouble;
//�β�x������������ֻ���������Ǹ��������Ա�����
//�÷��磺double y;_Trapz(y,t,atan2(t,7.12),0,4,11);
#define _iTrapz(y,a,b,N,func,x) \
{double x,sum##NAMEDEF,del##NAMEDEF,it##NAMEDEF=1.0;int i##NAMEDEF,j##NAMEDEF;\
x=a,y=func,x=b,y=y+func,y=0.5*(b-a)*y;\
for(i##NAMEDEF=2;i##NAMEDEF<=N;i##NAMEDEF++)\
{del##NAMEDEF=(b-a)/it##NAMEDEF,x=a+0.5*del##NAMEDEF,sum##NAMEDEF=0.0,j##NAMEDEF=1;\
for(;j##NAMEDEF<=it##NAMEDEF;j##NAMEDEF++,x+=del##NAMEDEF)sum##NAMEDEF+=func;\
y=0.5*(y+sum##NAMEDEF*del##NAMEDEF);it##NAMEDEF=it##NAMEDEF+it##NAMEDEF;}\
}


//double y�������
//����ֵ��ȷ��6λ��Ч���֡�
#define _iTrap(y,a,b,func,x) \
{double oldy##NAMEDEF;int kk##NAMEDEF;\
_iTrapz(oldy##NAMEDEF,a,b,1,func,x);\
_iTrapz(y,a,b,2,func,x);\
for(kk##NAMEDEF=3;(fabs(y-oldy##NAMEDEF)>=fabs((1.0e-6)*y))&&kk##NAMEDEF<20;kk##NAMEDEF++)\
{oldy##NAMEDEF=y;_iTrapz(y,a,b,kk##NAMEDEF,func,x);}}

double iTrapz(double a,double b,int N,double (*func)(double));
double iTrap(double a,double b,double (*func)(double));
double wal(double t,int n);
double Wal(double t,int n[1]);
double* WAL(int N);

//22:double coef�������
#define _WalCoef(coef,k,T,fx,t) \
{if(T<=0.0||k<0){printf("_WalCoef()��T����k�����������");_getch();exit(0);}\
else _iTrap(coef,0.0,T,fx*wal(t/T,k),t);\
coef=coef/T;}

double WalSer(double t,double *T,double a[],int n[1]);
double WalSer(double t,double T,double a[],int n);

//double coef�������
#define _fsinCoef(coef,k,T,fx,x) \
{double w##NAMEDEF=MY_PI*2.0/T;\
if(T<=0.0||k<0){printf("fsinCoef()��T����k�����������");_getch();exit(0);}\
else {_iTrap(coef,0.0,T,fx*sin(w##NAMEDEF*x*k),x);\
coef=coef*2/T;}}


#define _fcosCoef(coef,k,T,fx,x) \
{double w##NAMEDEF=MY_PI*2.0/T;\
if(T<=0.0||k<0){printf("fsinCoef()��T����k�����������");_getch();exit(0);}\
else {_iTrap(coef,0.0,T,fx*cos(w##NAMEDEF*x*k),x);\
coef=coef*2/T;}\
if(k==0)coef=coef/2;}

double FouSer(double t,double *T,
					 double a[],double b[],int *N);
double FouSer(double t,double T,
					 double a[],double b[],int N);

//��������
#define _iTConv(s,ti,a,N,func1,func2,x) \
{double x,s1,s2=0,del,it=1.0;int i,j;\
x=a,s1=func1,x=ti-a,s1=s1*func2,x=ti-a,s2=func1,x=a,s2=s2*func2,s=s1+s2,s=0.5*(ti-2*a)*s;\
for(i=2;i<=N;i++)\
{del=(ti-2*a)/it,x=a+0.5*del,s2=0.0,j=1;\
for(;j<=it;j++,x+=del){s1=func1;x=ti-x;s2+=func2*s1;x=ti-x;}\
s=0.5*(s+s2*del);it=it+it;}\
}

//��֪�����������ķ�����㣺t>=a��f1(x)=0,f2(x)=0ʱ������������Ϊ[a,t-a]
#define _iConv(y,ti,a,fx1,fx2,x) \
{double oldy;int kk;\
_iTConv(oldy,ti,a,1,fx1,fx2,x);\
_iTConv(y,ti,a,2,fx1,fx2,x);\
for(kk=3;(fabs(y-oldy)>=fabs((1.0e-6)*y))&&kk<20;kk++)\
{oldy=y;_iTConv(y,ti,a,kk,fx1,fx2,x);}}


#define _fTConv(y,t,a,b,N,fx1,fx2,x) \
{double x,s1,s2=0,del,it=1.0;int i,j;\
x=a,s1=fx1,s2=fx2,x=b,y=s1*fx2+s2*fx1,y=0.5*(b-a)*y;\
for(i=2;i<=N;i++)\
{del=(b-a)/it,x=a+0.5*del,s2=0.0,j=1;\
for(;j<=it;j++,x+=del){s1=fx1;x=t-x;s2+=fx2*s1;x=t-x;}\
y=0.5*(y+s2*del);it=it+it;}\
}

//��֪��������������һ�����ķ������䣺ֻ��[a��b]��f1(x)!=0,����f1(x)=0��
#define _fConv(y,t,a,b,fx1,fx2,x) \
{double oldy;int kk;\
_fTConv(oldy,t,a,b,1,fx1,fx2,x);\
_fTConv(y,t,a,b,2,fx1,fx2,x);\
for(kk=3;(fabs(y-oldy)>=fabs((1.0e-6)*y))&&kk<20;kk++)\
{oldy=y;_fTConv(y,t,a,b,kk,fx1,fx2,x);}}

#endif
