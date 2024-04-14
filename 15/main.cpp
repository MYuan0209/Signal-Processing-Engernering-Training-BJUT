#include<bits/stdc++.h>
#include<windows.h>
#include<conio.h>
#include"../x_math.h"
#include"../x_graph.h"
using namespace std;


extern WIN2 win2;
extern WIN3 win3;


//根据卷积公式编写卷积函数
//x[N]卷h[M],结果长度:L=M+N-1;
//由于:0<=k<=M+N-2,0<=n-k<=M-1,0<=k<=N-1
//所以:K>=0,K>=n-M+1;
//	 k<=N-1,K<=n; 
void convol(double x[], int N, double h[], int M, double y[], int L) {
	int n, k, kk;
	for(n = 0; n < L; n++){
		for(k=(n>=M)?n-M+1:0,kk=(N>n)?n:N-1,y[n]=0;k<=kk;k++) {
			y[n]+=x[k]*h[n-k];
		}
	}
}


int main()
{
	int N = 5, M = 3, L = 10, i;
	double x[8] = {1, 2, 3, 4, 5}, h[8] = {3, 2, 1}, y[20];
	COMPLEX f[8], r[8];
	//复数数组f[8]={1,0,2,0,3,0,4,0,5,0};
	//复数数组r[8]={3,0,2,0,1,0,0,0};
	for(i = 0; i < 5; i++) {
		f[i] = COMPLEX(i + 1, 0);
	}
	for(i = 5; i < 8; i++) {
		f[i] = COMPLEX(0, 0);
	}
	for(i = 0; i < 3; i++) {
		r[i] = COMPLEX(3 - i, 0);
	}
	for(i = 5; i < 8; i++) {
		r[i] = COMPLEX(0, 0);
	}

	/*第1步：用卷积函数直接做卷积*/
	printf("(1)直接卷积 \n");
	convol(h, M, x, N, y, L);
	for(i = 0; i < L; i++) {
		printf("y[%d]=%f\n", i, y[i]);
	}
	_getch();

	/*第2步：用fft做卷积,与直接法的结果比较验证fft卷积结果是否正确性。*/
	printf("(2) 用fft做卷积\n");
	for(i = M + N - 1, L = 1; i != 0; ) {
		i = i >> 1;
		L = L << 1;
	}
	fft(f, L, 1);
	fft(r, L, 1);
	for(i = 0; i < L; i++) {
		f[i] = f[i] * r[i];
	}
	fft(f, L, -1);
	for(i = 0; i < L; i++) {
		printf("f[%d]=%f+j%f,\n", i, f[i].r, f[i].i);
	}
	_getch();

	/*第3步  用FFT实现一个长的数据序列的卷积。*/
	printf("(3) 15个信号数据代表长序列，用一次FFT实现数据序列的卷积\n");
	//信号复数组 x[32]={1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9,0,10,0,11,0,12,0,13,0,14,0,15,0};
	//系统复数组 h[32]={3,0,2,0,1,0,};
	//信号数据实际长度为15，系统数据实际长度为3，
	//15+3-1=17，所以要做32的FFT，故两数组长度都声明为32。
	COMPLEX x2[32], h2[32];
	for(i = 0; i < 15; i++) {
		x2[i] = COMPLEX(i+1, 0);
	}
	for(i=15;i<32;i++)x2[i]=COMPLEX(0,0);
	for(i=0;i<3;i++)h2[i]=COMPLEX(3-i,0);
	for(i=3;i<32;i++)h2[i]=COMPLEX(0,0);

	fft(x2,32,1);
	fft(h2,32,1);
	for(i=0;i<32;i++)x2[i]=x2[i]*h2[i];
	fft(x2,32,-1);
	for(i=0;i<20;i++)printf("x2[%d]=%f+j%f,\n",i,x2[i].r,x2[i].i);
	_getch();


	/*第4步：叠接舍去法实现长序列的数据卷积。*/
	//以第3步的卷积结果为正确结果，检验下面叠接舍去法做出的结果。

	printf("(4) 15个数据用叠接舍去法做出的结果，检查与第3步结果是否一致\n");
	//长序列信号：x2[15+3+(8-(3+15)%5)]={[0,0,0,0,0,0,(1,0,2,0],3,0,<4,0,5,0),6,0,7,0>,
	//                                 8,0,9,0,10,0,11,0,12,0,13,0,14,0,15,0};
	//代表可做8个FFT的系统h2[8]={3,0,2,0,1,0}
	//每次分组所取信号数据暂存在y[8];
	COMPLEX x3[15 + 3 + (8 - (3 + 15) % 5)], h3[8], g[8];
	//x3[]的长度为15+3+(8-(3+15)%5)=23，
	//x3分组取数如下：
	//第1次取：0,0,0,0,0,0,1,0,2,0,|3,0,4,0, 5,0,这8个数据
	//第2次取：3,0,4,0,5,0,6,0,7,0,|8,0,9,0,10,0,这8个数据
	//第3次取：8,0,9,0,10,0,11,0,12,0,|13,0,14,0,15,0,这8个数据
	//第4次取：13,0,14,0,15,0,0,0,0,0,|0,0,0,0,0,0,这8个数据
	for(i=0;i<23;i++)x3[i]=COMPLEX(0,0);
	for(i=3;i<18;i++)x3[i]=COMPLEX(i-2,0);
	for(i=0;i<3;i++)h3[i]=COMPLEX(3-i,0);
	for(i=3;i<8;i++)h3[i]=COMPLEX(0,0);

	int j;
	fft(h3,8,1);
	for(j=0;j<=23-8;j=j+5){
		for(i=0;i<8;i++)g[i]=x3[i+j];
		fft(g,8,1);
		for(i=0;i<8;i++)g[i]=g[i]*h3[i];
		fft(g,8,-1);
		//思考：为什么从i=3开始？
		for(i=3;i<8;i++)printf("g[%d]=%f+%fj\n",i,g[i].r,g[i].i);
	}
	_getch();
	/*
	COMPLEX x4[15 + 3 + (8 - (3 + 15) % 6)], h4[8], y4[8];
	for (i = 0; i < 24; ++i) {
		x4[i] = COMPLEX(0, 0);
	}
	for (i = 2; i < 17; ++i) {
		x4[i] = COMPLEX(i - 1, 0);
	}
	for (i = 0; i < 3; ++i) {
		h4[i] = COMPLEX(3 - i, 0);
	}
	for (i = 3; i < 8; ++i) {
		h4[i] = COMPLEX(0, 0);
	}

	fft(h4, 8, 1);
	for (int j = 0; j < 24 - 8; j = j + 6) {
		for (i = 0; i < 8; ++i) {
			y4[i] = x4[i + j];
		}
		fft(y4, 8, 1);
		for (i = 0; i < 8; ++i) {
			y4[i] = y4[i] * h4[i];
		}
		fft(y4, 8, -1);
		for (i = 2; i < 8; ++i) {
			printf("g[%d]=%f+%fj\n", i, y4[i].r, y4[i].i);
		}
	}*/

	/*第5步：自己完成！实现叠接相加法，结果与上面正确卷积结果互相验证*/
	//问题：分组数多少？前后叠接项是多少？一次FFT得到的非0数据长度是多少？
	printf("(5) 实现叠接相加法，结果与上面正确卷积结果互相验证\n");
	COMPLEX x4[16], h4[8], g4[8], y4[32];
	for (i = 0; i < 15; ++i) {
		x4[i] = COMPLEX(i + 1, 0);
	}
	for (; i < 16; ++i) {
		x4[i] = COMPLEX(0, 0);
	}
	for (i = 0; i < 3; ++i) {
		h4[i] = COMPLEX(3 - i, 0);
	}
	for (; i < 8; ++i) {
		h4[i] = COMPLEX(0, 0);
	}
	for (i = 0; i < 32; ++i) {
		y4[i] = COMPLEX(0, 0);
	}

	fft(h4, 8, 1);
	for (int j = 0; j < 16 / 4; ++j) {
		for (i = 0; i < 4; ++i) {
			g4[i] = COMPLEX(0, 0);
		}
		for (i = 4; i < 8; ++i) {
			g4[i] = x4[j * 8 + i - 4];
		}
		fft(g4, 8, 1);
		for (i = 0; i < 8; ++i) {
			g4[i] = g4[i] * h4[i];
		}
		fft(g4, 8, -1);
		for (i = 0; i < 8; ++i) {
			printf("g[%d]=%f+%fj\n", i, g4[i].r, g4[i].i);
		}
	}
	return 0;
}
