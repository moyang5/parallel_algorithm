#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <algorithm>
#include <arm_neon.h>
#include <cstdint>
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
uint64_t quick_pow(uint64_t a, uint64_t b, uint64_t p) {
    uint64_t res = 1;
    a %= p;
    while(b) {
        if(b & 1) res = (res * a) % p;
        a = (a * a) % p;
        b >>= 1;
    }
    return res;
}

// 计算模逆元
int inv(int x ,int p) {
    return quick_pow(x, p - 2 ,p);
}

//蝴蝶变换
void bit_reverse(uint64_t *a, int n) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j >= bit; bit >>= 1) j -= bit;
        j += bit;
        if (i < j) std::swap(a[i], a[j]);
    }
}

// 蒙哥马利参数结构体

struct Montgomery {
    uint64_t N;
    uint64_t R;         // R=2^64 mod N
    uint64_t N_inv;     // N^{-1} mod R
    uint64_t R2modN;    // R² mod N

    Montgomery(uint64_t mod) : N(mod) {
        if (mod % 2 == 0) {
            std::cerr << "Modulus must be odd" << std::endl;
            exit(1);
        }

        // 正确计算R = 2^64 mod N
        R = 1;
        for (int i = 0; i < 64; ++i) {
            R = (R << 1) % N;
        }

        // 计算R² mod N
        R2modN = (__uint128_t)R * R % N;

        // 牛顿迭代法计算N的模逆元
        N_inv = 1;
        uint64_t t = 1ULL << 63;
        for (int i = 0; i < 5; ++i) {
            N_inv *= 2 - N * N_inv;
            t >>= 1;
        }
        N_inv = -N_inv;
    }

    uint64_t REDC(__uint128_t T) const {
        uint64_t m = ((uint64_t)T * N_inv);
        uint64_t t = (T + (__uint128_t)m * N) >> 64;
        return t >= N ? t - N : t;
    }

    uint64_t to_mont(uint64_t x) const {
        return REDC((__uint128_t)x * R2modN);
    }

    uint64_t from_mont(uint64_t x) const {
        return REDC(x);
    }

    uint64_t mul(uint64_t a, uint64_t b) const {
        return REDC((__uint128_t)a * b);
    }
};

void NTT(uint64_t* a, int len, const Montgomery& mont, int opt = 1) {
    bit_reverse(a, len);

    // 转换为蒙哥马利形式
    for (int i = 0; i < len; ++i) {
        a[i] = mont.to_mont(a[i] % mont.N);
    }

    for (int mid = 1; mid < len; mid <<= 1) {
        uint64_t exp = (mont.N - 1) / (mid << 1);
        if (opt == -1) exp = mont.N - 1 - exp;
        uint64_t Wn = quick_pow(g, exp, mont.N);
        Wn = mont.to_mont(Wn);

        for (int j = 0; j < len; j += (mid << 1)) {
            uint64_t w = mont.to_mont(1);

            for (int k = 0; k < mid; ++k) {
                uint64_t x = a[j + k];
                uint64_t y = mont.mul(a[j + k + mid], w);

                a[j + k] = (x + y) % mont.N;
                a[j + k + mid] = (x - y + mont.N) % mont.N;

                w = mont.mul(w, Wn);
            }
        }
    }

    if (opt == -1) {
        uint64_t inv_len = quick_pow(len, mont.N-2, mont.N); // 费马小定理
        for (int i = 0; i < len; ++i) {
            a[i] = mont.mul(a[i], inv_len);
            a[i] = mont.from_mont(a[i]);
            a[i] %= mont.N;
        }
    }
}

void poly_multiply_NTT(int* a, int* b, int* ab, int n, int p) {
    int len = 1;
    while (len < 2*n) len <<= 1;

    uint64_t* fa = new uint64_t[len]();
    uint64_t* fb = new uint64_t[len]();

    // 数据预处理
    for (int i = 0; i < n; ++i) {
        fa[i] = (uint64_t)((a[i] % p + p) % p);
        fb[i] = (uint64_t)((b[i] % p + p) % p);
    }

    Montgomery mont(p);

    NTT(fa, len, mont, 1);
    NTT(fb, len, mont, 1);

    for (int i = 0; i < len; ++i) {
        fa[i] = mont.mul(fa[i], fb[i]);
    }

    NTT(fa, len, mont, -1);

    // 结果处理
    for (int i = 0; i < 2*n-1; ++i) {
        ab[i] = (int)(fa[i] % p);
    }

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
    int test_end = 3;
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
