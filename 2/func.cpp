#include<bits/stdc++.h>
#include"func.h"
using namespace std;


double ramp(double t) {
    return t > 0 ? t : 0;
}


double sawtoothwave(double t, double T) {
    if (T <= 0) {
        cout << "sawtoothwave() input error" << endl;
        exit(1);
    }
    double result = fmod(t, T);
    if (t < 0 && result != 0) {
        result += T;
    }
    return result;
}


double sinc(double t, double w, double W) {
    if (w == 0.0 || t == 0.0) {
        return cos(W * t);
    }
    else {
        return sin(w * t) * cos(W * t) / w * t;
    }
}
