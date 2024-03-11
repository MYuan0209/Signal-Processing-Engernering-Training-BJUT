#include<bits/stdc++.h>
#include<windows.h>
#include<conio.h>
#include"../x_math.h"
#include"../x_graph.h"
using namespace std;

double _DMFRV_=0;//DOUBLE MACRO FUNCTION RETURN Value,全局变量,用于放置宏函数的返回值。


double func1(double t) {
	return t >= 0 ? sin(t) : 0;
}


double func2(double t) {
	return (t >= 0) ? exp(-t / 2) : 0;
}


double func3(double t) {
	if (t > -1 && t < 1) {
		return 1;
	}
	else if (t > 4 && t < 5) {
		return 3;
	}
	else {
		return 0;
	}
}


double func4(double t) {
	return (t > -5 && t < 5) ? t / 3 : 0;
}


//方法1：用最原始的矩形积分法，将做卷积的两函数在函数f1CONVf2()实现。
//观察该法的精度和速度，分析优缺点。
double f1CONVf2(double t) {
	double sum = 0, a = 0, dt = 0.01;
	int N;
	t = t - a * 2;
	N = (int)(t / dt);
	for(int k = 0;k<N;k++) {
		sum += func1(a + k * dt) * func2(a + (N - k) * dt);
	}
	return sum*dt;
}


//方法2：用积分宏函数_iTrapz(),改进方法1,注意宏参数ehtao的调用。
//将函数func2修改为：double func2(double t){return (-1<=t&&t<=3)?1:0;}比较修改前后的差异，找出原因。
//观察积分速度变化，总结该法的优缺点。
#define _plotConvf1(t,a,ehtao,tao)\
{double t;\
t=win2.x1;\
_iTrapz(_DMFRV_,(a),(t-(a)),10,ehtao,tao);\
moveto2(t,_DMFRV_);\
for(; t<win2.x2; t=t+win2.dxPerPixel){_iTrapz(_DMFRV_,(a),(t-(a)),10,ehtao,tao); lineto2(t,_DMFRV_);}\
t=win2.x2; _iTrapz(_DMFRV_,(a),(t-(a)),10,ehtao,tao); lineto2(t,_DMFRV_);}


//方法3：改进方法1，用宏来实现任意两个带参数函数的卷积，总结优缺点。
#define convf(t,a,dt,f1,f2,thita)\
{double sum=0,f1t,f2t;int N,k=0;\
if(dt<=0.0){printf("convf()函数的dt参数设置错误!\n"),exit(1);};\
double thita=(t)-2*(a);N=(int)((thita)/(dt));\
for(;k<N;k++){thita=a+k*(dt),f1t=f1,thita=a+(N-k)*(dt),f2t=f2,sum+=f1t*f2t;}\
_DMFRV_=sum*(dt);}

//编写专门针对convf()的绘图函数。
#define _plotConvf2(a,dt,f1,f2,thita)\
{double t;\
t=win2.x1;\
convf(t,a,dt,f1,f2,thita);\
moveto2(t,_DMFRV_);\
	for(;t<win2.x2;t=t+win2.dxPerPixel){convf(t,a,dt,f1,f2,thita);lineto2(t,_DMFRV_);}\
t=win2.x2;convf(t,a,dt,f1,f2,thita);lineto2(t,_DMFRV_);}

//方法4：把两函数的非0值区域采样后变成两数组，放入一个ehg结构中。观察运算速度，总结优缺点。
//注意取值区间[a,b]的设置，注意宏的输入参数名字不要与内部变量名i,N重名。 

struct ehg{int N;double *e,*h,a,b,dt;};

//在外面做struct ehg y;
#define ehArray(y,t1,t2,delt,f1t,f2t,t)\
{int i,N;\
double t;\
if((dt)<=0.0){printf("ehArray()的宏参数dt输入错误\n");exit(1);}\
N=(int)((t2-(t1))/(delt));\
y.e=new double[N];\
y.h=new double[N];\
y.N=N;y.a=t1;y.b=t2;y.dt=delt;\
for(i=0;i<N;i++)\
   {t=t1+i*(dt);y.e[i]=f1t;y.h[i]=f2t;}\
}

double convf2(double t,struct ehg g) {
	double sum=0;int N,k=0;
	t=t-2*g.a;
	N=(int)(t/g.dt);
	if(N>g.N){return 0;}
	for(k=0;k<N;k++)
		sum+=g.e[k]*g.h[N-k];
	return sum*g.dt;
}

int main() {
	extern WIN2 win2;
	double a = 0, dt=0.01;

	window2((WCHAR*)"工程训练5,卷积运算", -10,5,10,-5, (char*)"t", (char*)"f");

	//绘制两卷积的函数：
	HPEN pen1=CreatePen(PS_SOLID,3,RED),oldpen=(HPEN)SelectObject(win2.hdc,pen1);
	_plotxy2(t,func1(t));_getch();
	_plotxy2(t,func2(t));_getch();

	//方法1：
	HPEN pen2=CreatePen(PS_SOLID,1,BLUE);SelectObject(win2.hdc,pen2);
	_plotxy2(t,f1CONVf2(t));_getch();

	//方法2：
	HPEN pen3=CreatePen(PS_SOLID,2,GREEN);SelectObject(win2.hdc,pen3);
	_plotConvf1(t,-4,func1(tao)*func2((t-(tao))),tao);_getch();

	//方法3：
	HPEN pen4=CreatePen(PS_SOLID,2,YELLOW);SelectObject(win2.hdc,pen4);
	_plotConvf2(a,dt,func1(thita),func2(thita),thita);_getch();

	//方法4：
	HPEN pen5=CreatePen(PS_SOLID,3,BLACK);SelectObject(win2.hdc,pen5);
	struct ehg g;
	ehArray(g, 0, 10,dt,func1(t),func2(t),t);
	_plotxy2(t,convf2(t,g));
	_getch();

	delete[] g.e;
	delete[] g.h;

	a = -5;
	frame2(win2.xstr, win2.ystr);
	SelectObject(win2.hdc, pen1);
	_plotxy2(t, func3(t));
	_getch();
	_plotxy2(t, func4(t));
	_getch();
	SelectObject(win2.hdc, pen2);
	_plotConvf2(a,dt,func3(thita),func4(thita),thita);
	_getch();
	SelectObject(win2.hdc, pen3);
	_plotConvf1(t, a,func3(tao)*func4((t-(tao))),tao);
	_getch();
	SelectObject(win2.hdc, pen4);
	ehArray(g, a, 20, dt,func3(t),func4(t),t);
	_plotxy2(t,convf2(t,g));
	delete[] g.e;
	delete[] g.h;

	SelectObject(win2.hdc,oldpen);//恢复原来画笔.
	DeleteObject(pen1);//删除新建画笔.
	DeleteObject(pen2);//删除新建画笔.
	DeleteObject(pen3);//删除新建画笔.
	DeleteObject(pen4);//删除新建画笔.
	DeleteObject(pen5);//删除新建画笔.
	return 0;
}