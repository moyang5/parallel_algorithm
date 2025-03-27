#include <iostream>
#include <windows.h>
using namespace std;
// inline void normal(int **b,int *a,int n){
//     int sum[n] ;
//     for (int i = 0 ; i < n ;i++){
//         sum[i] = 0 ;
//         for (int j = 0 ; j < n ; j++){
//             sum[i] += b[j][i] * a[j] ;        
//         }
        
//     }
//     return ;

// }

inline void cache(int **b ,int *a, int n){
    int sum[n] ;
    for (int i = 0 ;i < n ;i++){
        sum[i] = 0 ;
    }
    for (int j = 0 ; j < n ;j++){
        for (int i = 0 ; i < n ; i++){
            sum[i] += b[j][i] * a[j] ;        
        }         
    }
}



int main(){
    int n = 5000 ;
    int **b = new int*[n];
    for (int i = 0; i < n ;i++){
        b[i] = new int[n] ;
        for (int j = 0 ; j < n ; j++){
            b[i][j] = i + j;
        }
    }
    int a[n];
    for(int i = 0 ; i< n ;i++){
        a[i] = i ;
    }
    long long head, tail, freq;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
    // QueryPerformanceCounter((LARGE_INTEGER *)&head);
    // normal( b, a, n) ;
    // QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    // cout<<"time: "<< (tail - head) * 1000.0 / freq << "ms" << endl ;
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    cache(b,a,n) ;
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    cout<<"time: "<< (tail - head) * 1000.0 / freq << "ms" << endl ;
    return 0 ;
}