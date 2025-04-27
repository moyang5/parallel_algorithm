#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <algorithm>
// 可以自行添加需要的头文件

const int g = 3;//原根


void fRead(int *a, int *b, int *n, int *p, int input_id){
    // 数据输入函数
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strin = str1 + str2 + ".in";
    char data_path[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), data_path);
    data_path[strin.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    fin>>*n>>*p;
    for (int i = 0; i < *n; i++){
        fin>>a[i];
    }
    for (int i = 0; i < *n; i++){   
        fin>>b[i];
    }
}

void fCheck(int *ab, int n, int input_id){
    // 判断多项式乘法结果是否正确
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    for (int i = 0; i < n * 2 - 1; i++){
        int x;
        fin>>x;
        if(x != ab[i]){
            std::cout<<"多项式乘法结果错误"<<std::endl;
            return;
        }
    }
    std::cout<<"多项式乘法结果正确"<<std::endl;
    return;
}

void fWrite(int *ab, int n, int input_id){
    // 数据输出函数, 可以用来输出最终结果, 也可用于调试时输出中间数组
    std::string str1 = "files/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++){
        fout<<ab[i]<<'\n';
    }
}


void poly_multiply(int *a, int *b, int *ab, int n, int p){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            ab[i+j]=(1ll * a[i] * b[j] % p + ab[i+j]) % p;
        }
    }
}

//快速幂算法
long long quick_pow(long long a, long long b, long long p) {
    long long res = 1;
    a %= p;
    while(b) {
        if(b & 1) res = res * a % p;
        a = a * a % p;
        b >>= 1;
    }
    return res;
}

// 计算模逆元
int inv(int x ,int p) {
    return quick_pow(x, p - 2 ,p);
}

//蝴蝶变换
void bit_reverse(int *a, int n) {
    for (int i = 0, j = 0; i < n; i++) {
        if(i < j) std::swap(a[i], a[j]);
        for(int k = n >> 1; (j ^= k) < k; k >>= 1);
    }
}

    void NTT(int *a, int len, int p, int opt = 1) { 
        bit_reverse(a, len);
        
        for(int mid = 1; mid < len; mid <<= 1) {
            long long Wn = quick_pow(g, (p - 1) / (mid << 1), p);
            if(opt == -1) Wn = inv(Wn, p);
            
            for(int j = 0; j < len; j += (mid << 1)) {
                long long w = 1;
                for(int k = 0; k < mid; k++, w = w * Wn % p) {
                    long long x = a[j + k], y = w * a[j + k + mid] % p;
                    a[j + k] = (x + y) % p;
                    a[j + k + mid] = (x - y + p) % p;
                }
            }
        }
        
        if(opt == -1) {
            long long inv_len = inv(len, p);
            for(int i = 0; i < len; i++) {
                a[i] = a[i] * inv_len % p;
            }
        }
    }

void poly_multiply_NTT(int *a, int *b, int *ab, int n, int p) {
    int len = 1;
    while (len < 2 * n - 1) len <<= 1;
    
    int *fa = new int[len]();
    int *fb = new int[len]();
    std::copy(a, a + n, fa);
    std::copy(b, b + n, fb);
    
    NTT(fa, len, p); 
    NTT(fb, len, p);
    
    for (int i = 0; i < len; i++) {
        fa[i] = 1LL * fa[i] * fb[i] % p;
    }
    
    NTT(fa, len, p, -1); 
    
    std::copy(fa, fa + 2*n-1, ab);
    
    delete[] fa;
    delete[] fb;
}


int a[300000], b[300000], ab[300000];
int main(int argc, char *argv[])
{
    
    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 263882790666241
    // 第四个模数超过了整型表示范围, 如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求, 如果要自行探索大模数 NTT, 请在完成前三个模数的基础代码及优化后实现大模数 NTT
    // 输入文件共五个, 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久, 推荐调试正确性时只使用输入文件 1
    int test_begin = 0;
    int test_end = 1;
    for(int i = test_begin; i <= test_end; ++i){
        long double ans = 0;
        int n_, p_;
        fRead(a, b, &n_, &p_, i);
        memset(ab,0,sizeof(ab));
        auto Start = std::chrono::high_resolution_clock::now();
        // TODO : 将 poly_multiply 函数替换成你写的 ntt
        poly_multiply_NTT(a, b, ab, n_, p_);
        auto End = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double,std::ratio<1,1000>>elapsed = End - Start;
        ans += elapsed.count();
        fCheck(ab, n_, i);
        std::cout<<"average latency for n = "<<n_<<" p = "<<p_<<" : "<<ans<<" (us) "<<std::endl;
        // 可以使用 fWrite 函数将 ab 的输出结果打印到 files 文件夹下
        // 禁止使用 cout 一次性输出大量文件内容
        fWrite(ab, n_, i);
    }
    return 0;
}
