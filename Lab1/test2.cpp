#include <iostream>
#include <windows.h>
using namespace std ;

inline void normal(int *a,int n){
    int sum = 0 ;
    for (int i = 0 ; i < n;i++){
        sum += a[i] ;
    }
    return ;
}

inline void chain(int *a,int n){
    int sum1 = 0 ,sum2 = 0 ;
    for (int i = 0 ; i < n ; i+=2){
        sum1 += a[i] ;
        sum2 += a[i+1] ;
    }
    int sum = sum1 + sum2 ;
    return ;
}

inline void recursion(int *a ,int n){
    if (n == 1) return ;
    else{
        for (int i = 0 ; i < n/2 ; i++)
            a[i] += a[n - i - 1] ;
        n = n/2 ;
        recursion(a,n) ;
    }
    return ;
}
inline void twoloop(int *a ,int n){
    for (int m = n ; m > 1; m /= 2){
        for (int i = 0 ; i < m/2 ; i++){
            a[i] = a[2*i] + a[2*i+1] ;
        }
    }
    return ;
}

int main(){
    int n ;
    n = 5000 ;
    int a[n] ;
    for (int i = 0; i < n ; i++){
        a[i] = 2*i ;
    }
    long long head, tail, freq;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    normal(a,n) ;
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    cout<<"time: "<< (tail - head) * 1000.0 / freq << "ms" << endl ;

    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    chain(a,n) ;
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    cout<<"time: "<< (tail - head) * 1000.0 / freq << "ms" << endl ;

    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    recursion(a,n) ;
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    cout<<"time: "<< (tail - head) * 1000.0 / freq << "ms" << endl ;

    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    twoloop(a,n) ;
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    cout<<"time: "<< (tail - head) * 1000.0 / freq << "ms" << endl ;
    return 0 ;
}