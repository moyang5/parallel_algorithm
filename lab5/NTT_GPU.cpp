#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include <chrono>
#include <random>

#define CUDA_CHECK_ERROR(call) \
    do { \
        cudaError_t err = (call); \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA Error at " << __FILE__ << ":" << __LINE__ \
                      << " - " << cudaGetErrorString(err) << std::endl; \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

const int g = 3; // 原根

// 快速幂函数
__host__ __device__ long long quick_pow(long long a, long long b, long long p) {
    long long res = 1;
    a %= p;
    while (b) {
        if (b & 1) res = res * a % p;
        a = a * a % p;
        b >>= 1;
    }
    return res;
}

// 模逆元计算
__host__ __device__ long long mod_inv(long long x, long long p) {
    return quick_pow(x, p - 2, p);
}

// 位反转置换内核
__global__ void bit_reverse_kernel(int* a, const int* rev, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n && idx < rev[idx]) {
        int temp = a[idx];
        a[idx] = a[rev[idx]];
        a[rev[idx]] = temp;
    }
}

// NTT蝴蝶操作内核 - 修改为两层循环划分
__global__ void ntt_butterfly_kernel(int* a, int len, int p, int mid, const int* w, 
                                    int start_group, int end_group) {
    int group = blockIdx.x * blockDim.x + threadIdx.x + start_group;
    if (group >= end_group) return;
    
    for (int k = 0; k < mid; k++) {
        int j = group * 2 * mid;
        if (j + k < len && j + k + mid < len) {
            long long x = a[j + k];
            long long y = (long long)a[j + k + mid] * w[k] % p;
            a[j + k] = (x + y) % p;
            a[j + k + mid] = (x - y + p) % p;
        }
    }
}

// 缩放内核
__global__ void scale_kernel(int* a, int len, long long inv_len, int p) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < len) {
        a[idx] = (long long)a[idx] * inv_len % p;
    }
}

// GPU加速的NTT
void NTT_GPU(std::vector<int>& a, int p, int opt = 1) {
    int len = a.size();
    int* d_a = nullptr;
    CUDA_CHECK_ERROR(cudaMalloc(&d_a, len * sizeof(int)));
    CUDA_CHECK_ERROR(cudaMemcpy(d_a, a.data(), len * sizeof(int), cudaMemcpyHostToDevice));

    // 预计算位反转数组
    std::vector<int> h_rev(len);
    int* d_rev = nullptr;
    h_rev[0] = 0;
    for (int i = 1; i < len; i++) {
        h_rev[i] = (h_rev[i >> 1] >> 1) | ((i & 1) ? (len >> 1) : 0);
    }
    CUDA_CHECK_ERROR(cudaMalloc(&d_rev, len * sizeof(int)));
    CUDA_CHECK_ERROR(cudaMemcpy(d_rev, h_rev.data(), len * sizeof(int), cudaMemcpyHostToDevice));

    // 执行位反转置换
    int threads = 256;
    int blocks = (len + threads - 1) / threads;
    bit_reverse_kernel<<<blocks, threads>>>(d_a, d_rev, len);
    CUDA_CHECK_ERROR(cudaGetLastError());
    CUDA_CHECK_ERROR(cudaDeviceSynchronize());
    CUDA_CHECK_ERROR(cudaFree(d_rev));

    // 预分配旋转因子设备内存
    int* d_w = nullptr;
    CUDA_CHECK_ERROR(cudaMalloc(&d_w, (len/2) * sizeof(int)));
    
    // 迭代计算NTT
    for (int mid = 1; mid < len; mid <<= 1) {
        int step = mid << 1;
        // 计算主单位根
        long long w_n = quick_pow(g, (p-1)/step, p);
        if (opt == -1) w_n = mod_inv(w_n, p);
        
        // 计算当前层的旋转因子数组
        std::vector<int> h_w(mid);
        for (int k = 0; k < mid; k++) {
            h_w[k] = quick_pow(w_n, k, p);
        }
        
        // 复制旋转因子到设备
        CUDA_CHECK_ERROR(cudaMemcpy(d_w, h_w.data(), mid * sizeof(int), cudaMemcpyHostToDevice));

        // 计算组数和划分策略
        int num_groups = len / step;
        int groups_per_block = 16; // 每组块处理16个组
        int num_blocks = (num_groups + groups_per_block - 1) / groups_per_block;
        threads = 256;
        
        // 执行蝴蝶操作 - 使用两层循环划分
        ntt_butterfly_kernel<<<num_blocks, threads>>>(d_a, len, p, mid, d_w, 
                                                   0, num_groups);
        CUDA_CHECK_ERROR(cudaGetLastError());
        CUDA_CHECK_ERROR(cudaDeviceSynchronize());
    }

    // 逆变换：乘以长度逆元
    if (opt == -1) {
        long long inv_len = mod_inv(len, p);
        scale_kernel<<<(len + 255)/256, 256>>>(d_a, len, inv_len, p);
        CUDA_CHECK_ERROR(cudaGetLastError());
        CUDA_CHECK_ERROR(cudaDeviceSynchronize());
    }

    // 将结果拷贝回主机
    CUDA_CHECK_ERROR(cudaMemcpy(a.data(), d_a, len * sizeof(int), cudaMemcpyDeviceToHost));

    // 释放设备内存
    CUDA_CHECK_ERROR(cudaFree(d_w));
    CUDA_CHECK_ERROR(cudaFree(d_a));
}

// 多项式乘法
void poly_multiply_NTT(const std::vector<int>& a, const std::vector<int>& b, std::vector<int>& ab, int p) {
    int n = a.size();
    int len = 1;
    while (len < 2 * n - 1) len <<= 1;

    std::vector<int> fa(len, 0);
    std::vector<int> fb(len, 0);
    std::copy(a.begin(), a.end(), fa.begin());
    std::copy(b.begin(), b.end(), fb.begin());

    NTT_GPU(fa, p, 1);
    NTT_GPU(fb, p, 1);

    for (int i = 0; i < len; i++) {
        fa[i] = 1LL * fa[i] * fb[i] % p;
    }

    NTT_GPU(fa, p, -1);
    ab.assign(fa.begin(), fa.begin() + 2 * n - 1);
}

// 生成随机多项式
std::vector<int> generate_random_poly(int n, int p) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, p - 1);

    std::vector<int> poly(n);
    for (int i = 0; i < n; ++i) {
        poly[i] = dis(gen);
    }
    return poly;
}

int main() {
    const int n = 28734;  
    const int p = 998244353;

    // 生成随机多项式
    std::cout << "生成随机多项式..." << std::endl;
    std::vector<int> a = generate_random_poly(n, p);
    std::vector<int> b = generate_random_poly(n, p);
    std::vector<int> ab_ntt;

    // 测试NTT算法
    std::cout << "运行NTT算法..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    poly_multiply_NTT(a, b, ab_ntt, p);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed_ntt = end - start;
    
    std::cout << "NTT算法耗时: " << elapsed_ntt.count() << " ms" << std::endl;
    std::cout << "多项式长度: " << n << std::endl;
    std::cout << "结果长度: " << ab_ntt.size() << std::endl;

    // 输出前10个结果作为示例
    std::cout << "前10个结果: ";
    for (int i = 0; i < 10 && i < ab_ntt.size(); ++i) {
        std::cout << ab_ntt[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}