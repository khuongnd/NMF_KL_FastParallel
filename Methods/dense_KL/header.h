#ifndef H_HEADER_H
#define H_HEADER_H
#include <cmath>
#include <cstring>
#include <algorithm>
#include "pthread.h"
#include "math.h"
#include "mex.h"
#include <time.h>
#include <cmath>
#include <stdlib.h>
#include <ctime>
#define MAX_SUPPORTED_THREADS 100


#define For(i,n) for(int i=0; i<n; i++)
#define ForI(i, a, b) for (int i=a; i<b; i++)
#define ForC(it, v) for(map<int, double>::iterator it=v.begin(); it!=v.end(); it++)

//#undefing max
//#define max(a, b) (a > b ?  a : b)
#define abs(x) (x >= 0 ? x : -x)
#define eps 1e-9

using namespace std;

#define printfFnc(...) { mexEvalString("disp('')"); mexPrintf(__VA_ARGS__); mexEvalString("disp('')");}

double getTime(){
    return clock();
}
double getTime(double begin){
    return (clock() - begin)/CLOCKS_PER_SEC;
}

void printArray(double** a, int m, int n){
    For(i, m){
        For(j, n) mexPrintf("%.4f ", a[i][j]);
        mexPrintf("\n");
    }
}
#endif