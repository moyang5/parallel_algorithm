#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sys/time.h>
#include <algorithm>
#include <pthread.h>
#include <vector>
#include <mpi.h>

const int g = 3;

void fRead(int *a, int *b, int *n, int *p, int input_id) {
    std::string strin = "/nttdata/" + std::to_string(input_id) + ".in";
    std::ifstream fin(strin);
    fin >> *n >> *p;
    for (int i = 0; i < *n; i++) fin >> a[i];
    for (int i = 0; i < *n; i++) fin >> b[i];
}

void fCheck(int *ab, int n, int input_id) {
    std::string strout = "/nttdata/" + std::to_string(input_id) + ".out";
    std::ifstream fin(strout);
    for (int i = 0; i < n * 2 - 1; i++) {
        int x;
        fin >> x;
        if (x != ab[i]) {
            std::cout << "多项式乘法结果错误" << std::endl;
            return;
        }
    }
    std::cout << "多项式乘法结果正确" << std::endl;
}

void fWrite(int *ab, int n, int input_id) {
    std::string strout = "files/" + std::to_string(input_id) + ".out";
    std::ofstream fout(strout);
    for (int i = 0; i < n * 2 - 1; i++) fout << ab[i] << '\n';
}

long long quick_pow(long long a, long long b, long long p) {
    long long res = 1;
    a %= p;
    while (b) {
        if (b & 1) res = res * a % p;
        a = a * a % p;
        b >>= 1;
    }
    return res;
}

long long inv(long long x, long long p) {
    return quick_pow(x, p - 2, p);
}

void bit_reverse(int *a, int n) {
    for (int i = 0, j = 0; i < n; i++) {
        if (i < j) std::swap(a[i], a[j]);
        for (int k = n >> 1; (j ^= k) < k; k >>= 1);
    }
}

void NTT(int* array, int length, int modulus, int opt = 1) {
    bit_reverse(array, length);
    for (int current_size = 1; current_size < length; current_size <<= 1) {
        long long root = quick_pow(g, (modulus - 1) / (current_size << 1), modulus);
        if (opt == -1) root = inv(root, modulus);

        for (int i = 0; i < length; i += 2 * current_size) {
            long long w = 1;
            for (int j = 0; j < current_size; ++j) {
                int idx_even = i + j;
                int idx_odd = i + j + current_size;
                long long even = array[idx_even];
                long long odd = array[idx_odd] * w % modulus;

                array[idx_even] = (even + odd) % modulus;
                array[idx_odd] = (even - odd + modulus) % modulus;
                w = w * root % modulus;
            }
        }
    }

    if (opt == -1) {
        long long inv_length = inv(length, modulus);
        for (int i = 0; i < length; ++i)
            array[i] = static_cast<long long>(array[i]) * inv_length % modulus;
    }
}

void crt_merge(const std::vector<std::vector<int>>& ab_mod, int* ab, int n, int p_target) {
    const int mods[] = {998244353, 1004535809, 469762049};
    const __int128 m1 = mods[0], m2 = mods[1], m3 = mods[2];
    const __int128 m12 = m1 * m2;

    const __int128 inv_m1 = inv(m1 % m2, m2);
    const __int128 m12_mod_m3 = (m1 % m3) * (m2 % m3) % m3;
    const __int128 inv_m12 = inv(m12_mod_m3, m3);

    for (int i = 0; i < 2 * n - 1; ++i) {
        int a1 = ab_mod[0][i], a2 = ab_mod[1][i], a3 = ab_mod[2][i];
        __int128 diff = (a2 - a1) % m2;
        if (diff < 0) diff += m2;
        __int128 k = (diff * inv_m1) % m2;
        __int128 x12 = a1 + k * m1;
        x12 = (x12 % m12 + m12) % m12;

        __int128 x12_mod_m3 = x12 % m3;
        if (x12_mod_m3 < 0) x12_mod_m3 += m3;
        __int128 diff3 = (a3 - x12_mod_m3) % m3;
        if (diff3 < 0) diff3 += m3;
        k = (diff3 * inv_m12) % m3;
        __int128 x = x12 + k * m12;

        x = (x % (m12 * m3) + m12 * m3) % (m12 * m3);
        ab[i] = (x % p_target + p_target) % p_target;
    }
}

void poly_multiply_CRT(int *a, int *b, int *ab, int n, int p_target) {
    const int mods[] = {998244353, 1004535809, 469762049};
    const int num_mods = 3;
    int length = 1;
    while (length < 2 * n - 1) length <<= 1;

    std::vector<std::vector<int>> ab_mod(num_mods, std::vector<int>(length, 0));
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size < 4) {
        if (rank == 0) {
            std::cerr << "需要至少4个进程（1个主进程+3个工作进程）" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (rank == 0) {
        // 主进程：分发任务
        for (int m = 1; m <= 3; m++) {
            MPI_Send(&n, 1, MPI_INT, m, 0, MPI_COMM_WORLD);
            MPI_Send(&length, 1, MPI_INT, m, 1, MPI_COMM_WORLD);
            MPI_Send(&mods[m-1], 1, MPI_INT, m, 2, MPI_COMM_WORLD);
            MPI_Send(a, n, MPI_INT, m, 3, MPI_COMM_WORLD);
            MPI_Send(b, n, MPI_INT, m, 4, MPI_COMM_WORLD);
        }

        // 接收结果
        for (int m = 1; m <= 3; m++) {
            MPI_Recv(ab_mod[m-1].data(), 2*n-1, MPI_INT, m, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // 合并结果
        crt_merge(ab_mod, ab, n, p_target);
    } 
    else if (rank >= 1 && rank <= 3) {
        // 工作进程：接收数据并计算
        int local_n, local_length, modulus;
        MPI_Recv(&local_n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&local_length, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&modulus, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        std::vector<int> local_a(local_length, 0);
        std::vector<int> local_b(local_length, 0);
        std::vector<int> local_ab(local_length, 0);

        MPI_Recv(local_a.data(), local_n, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(local_b.data(), local_n, MPI_INT, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // 计算NTT
        NTT(local_a.data(), local_length, modulus, 1);
        NTT(local_b.data(), local_length, modulus, 1);
        
        for (int i = 0; i < local_length; ++i) {
            local_ab[i] = static_cast<long long>(local_a[i]) * local_b[i] % modulus;
        }
        
        NTT(local_ab.data(), local_length, modulus, -1);

        // 发送结果回主进程
        MPI_Send(local_ab.data(), 2*local_n-1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    // 确保所有进程同步
    MPI_Barrier(MPI_COMM_WORLD);
}

int a[300000], b[300000], ab[300000];
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int test_begin = 0, test_end = 4;
    for (int i = test_begin; i <= test_end; ++i) {
        long double ans = 0;
        int n_ = 0, p_ = 0;
        
        // 只有主进程读取输入
        if (rank == 0) {
            fRead(a, b, &n_, &p_, i);
            memset(ab, 0, sizeof(ab));
        }

        // 广播必要参数（所有进程都必须参与）
        MPI_Bcast(&n_, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&p_, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(a, n_, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(b, n_, MPI_INT, 0, MPI_COMM_WORLD);

        auto Start = std::chrono::high_resolution_clock::now();
        
        poly_multiply_CRT(a, b, ab, n_, p_);
        
        if (rank == 0) {
            auto End = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::ratio<1, 1000>> elapsed = End - Start;
            ans += elapsed.count();
            fCheck(ab, n_, i);
            std::cout << "average latency for n = " << n_ << " p = " << p_ << " : " << ans << " (us) " << std::endl;
            fWrite(ab, n_, i);
        }
    }
    
    MPI_Finalize();
    return 0;
}
