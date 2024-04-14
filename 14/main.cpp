#include<bits/stdc++.h>
#include<windows.h>
#include<conio.h>
#include"../x_math.h"
#include"../x_graph.h"
using namespace std;


extern WIN2 win2;
extern WIN3 win3;


void IIR(int bandType, double db1, double db2, double fs, double f1, double f2, double f3, double f4) {
	int order, L;
	double b[10][2], c[2][3], H[10][2][5];

	order = btwOrder(bandType, db1, db2, fs, f1, f2, f3, f4);

	L = btwAf(b, order);

	btwC23(c, bandType, order, db1, fs, f1, f2, f3, f4);

	btwAf2Df(H, L, b, c);

	//任务1：用printf()格式化把系统函数的系数H数组显示出来，以后理解数组的设计含义。
	//任务2：结合训练13中的数字滤波器制作，用设计的滤波器，对一信号（可以是方波、三角波等）进行抽样滤波。

	window2((WCHAR*)"lowpass",-1.,5.,900.,-60., (char*)"hz", (char*)"db");
	xy2(BLUE);
	plotxy2(RED,2,f,btw20lgHz(f,fs,H,L));

	//绘制通阻带参考线：
	if(bandType==LOWPASS||bandType==HIGHPASS) {
		line2(0,-db1,win2.x2,-db1);
		line2(0,-db2,win2.x2,-db2);
		line2(f1,0,f1,win2.y1);
		line2(f2,0,f2,win2.y1);
	}
	else{
		line2(0,-db1,win2.x2,-db1);
		line2(0,-db2,win2.x2,-db2);
		line2(f1,0,f1,win2.y1);
		line2(f2,0,f2,win2.y1);
		line2(f3,0,f3,win2.y1);
		line2(f4,0,f4,win2.y1);
	}

	_getch();
}


double func1(int i,double X[]){
   if(i>=0)return X[i];
   else return 0;}


double filter(int j,double X[] ,double Y[]  ,double H[][5]){
    int i;
	int result=0;
	for(i=0;i<5;i++){ result+=H[0][i]*func1(j-i,X);
	}
	for(i=1;i<5;i++){ result+=H[1][i]*func1(j-i,Y);}
	return result;
}

int main()
{
	IIR(1,3.0,18.0,2000.,200.,300.,400.,500.);
	IIR(2,3.0,18.0,2000.,200.,300.,400.,500.);
	IIR(3,3.0,18.0,2000.,200.,300.,400.,500.);
	IIR(4,3.0,18.0,2000.,200.,300.,400.,500.);
	_getch();

	
	return 0;
}
