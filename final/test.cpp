#include<bits/stdc++.h>
#include<windows.h>
#include<conio.h>
#include"../x_math.h"
#include"../x_graph.h"
#include"../func.h"
using namespace std;


extern WIN2 win2;
extern WIN3 win3;


void sampling (double *x, double (*f)(double, double), int N) {
	int T = 1;
	for (int i = 0; i < N; ++i) {
		x[i] = f(1.0 * i / N, T) + 0.05 * sin(2 * MY_PI * 450 * i / N);
	}
}


int main(int argc, char* argv[])
{
	int N = 100;
    double temp_max;
    double *rectang = new double[N];
    for (int i = 0; i < 3 * N / 4; ++i) {
        rectang[i] = 1;
    }
    for (int i = 3 * N / 4; i < N; ++i) {
        rectang[i] = 0;
    }
    double *hanning = new double[N];
    for (int i = 0; i < 3 * N / 4; ++i) {
        hanning[i] = 0.5 * (1.0 - cos(2 * i * MY_PI /(3 * N / 4 - 1)));
    }
    for (int i = 3 * N / 4; i < N; ++i) {
        hanning[i] = 0;
    }
    double *hamming = new double[N];
    for (int i = 0; i < 3 * N / 4; ++i) {
        hamming[i] = 0.54 - 0.46 * cos(2 * i*  MY_PI / (3 * N / 4 - 1));
    }
    for (int i = 3 * N / 4; i < N; ++i) {
        hamming[i] = 0;
    }

    window2((WCHAR*)"IIR Digital Fitter", 0, 0, N, 1.2, (char*)"n", (char*)"y[n]");
	xy2(BLACK);
	plotgri2(WHITE, GREEN, rectang, N);
	_getch();
    plotgri2(WHITE, RED, hanning, N);
	_getch();
    plotgri2(WHITE, BLUE, hamming, N);
	_getch();

    COMPLEX *Rectang = new COMPLEX[N];
    COMPLEX *result_rec = new COMPLEX[N];
    for (int i = 0; i < N; ++i) {
        Rectang[i].r = rectang[i];
        Rectang[i].i = 0;
    }
    dft(result_rec, Rectang, N, 1);
    for (int i = 0; i < N; ++i) {
        rectang[i] = abs(result_rec[i]);
    }
    temp_max = *max_element(rectang, rectang + N);
    for (int i = 0; i < N; ++i) {
        rectang[i] = 20 * log10(rectang[i] / temp_max);
    }
    COMPLEX *Hanning = new COMPLEX[N];
    COMPLEX *result_han = new COMPLEX[N];
    for (int i = 0; i < N; ++i) {
        Hanning[i].r = hanning[i];
        Hanning[i].i = 0;
    }
    dft(result_han, Hanning, N, 1);
    for (int i = 0; i < N; ++i) {
        hanning[i] = abs(result_han[i]);
    }
    temp_max = *max_element(hanning, hanning + N);
    for (int i = 0; i < N; ++i) {
        hanning[i] = 20 * log10(hanning[i] / temp_max);
    }
    COMPLEX *Hamming = new COMPLEX[N];
    COMPLEX *result_ham = new COMPLEX[N];
    for (int i = 0; i < N; ++i) {
        Hamming[i].r = hamming[i];
        Hamming[i].i = 0;
    }
    dft(result_ham, Hamming, N, 1);
    for (int i = 0; i < N; ++i) {
        hamming[i] = abs(result_ham[i]);
    }
    temp_max = *max_element(hamming, hamming + N);
    for (int i = 0; i < N; ++i) {
        hamming[i] = 20 * log10(hamming[i] / temp_max);
    }

    window2((WCHAR*)"frequency spectrum", 0, -200, N / 2, 0, (char*)"k", (char*)"W[k]");
	xy2(BLACK);
	plot_frequency_spectrum(BLUE, rectang, N / 2);
	_getch();
    frame2((char*)"k", (char*)"W[k]");
    xy2(BLACK);
    plot_frequency_spectrum(BLUE, hanning, N / 2);
	_getch();
     frame2((char*)"k", (char*)"W[k]");
    xy2(BLACK);
    plot_frequency_spectrum(BLUE, hamming, N / 2);
	_getch();

    int M = 50+1, fs = 2000, fc = fs / 4;
	double *b_rec = new double[M];
	firDesgin(b_rec, M, 1, RECTANG, 1.0 * fs, 1.0 * fc);
    double *b_han = new double[M];
	firDesgin(b_han, M, 1, HANNING, 1.0 * fs, 1.0 * fc);
    double *b_ham = new double[M];
	firDesgin(b_ham, M, 1, HAMMING, 1.0 * fs, 1.0 * fc);

    window2((WCHAR*)"LOWPASS", 0, -150, fs / 2, 5, (char*)"hz", (char*)"db");
	xy2(BLACK);
    line2(fc, win2.y2, fc, win2.y1);
	plotxy2(GREEN, 2, f, 20 * log10(firAbs(f, fs, b_rec, M)));
    _getch();
    plotxy2(RED, 2, f, 20 * log10(firAbs(f, fs, b_han, M)));
    _getch();
    plotxy2(BLUE, 2, f, 20 * log10(firAbs(f, fs, b_ham, M)));
	_getch();
	return 0;
}