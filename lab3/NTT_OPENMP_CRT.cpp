#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <omp.h>
#include <vector>

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

void parallel_NTT(int* array, int length, int modulus, int opt = 1) {
    bit_reverse(array, length);

    for (int current_size = 1; current_size < length; current_size <<= 1) {
        int num_blocks = length / (current_size << 1);
        if (num_blocks == 0) break;

        long long root = quick_pow(g, (modulus - 1) / (current_size << 1), modulus);
        if (opt == -1) root = inv(root, modulus);

        #pragma omp parallel for schedule(static)
        for (int block = 0; block < num_blocks; ++block) {
            int i = block * (current_size << 1);
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
        #pragma omp parallel for
        for (int i = 0; i < length; ++i)
            array[i] = static_cast<long long>(array[i]) * inv_length % modulus;
    }
}

void crt_worker(int m, int *a, int *b, int n, int length, const int* mods, std::vector<std::vector<int>> &ab_mod) {
    int p = mods[m];
    std::vector<int> fa(length, 0), fb(length, 0);
    std::copy(a, a + n, fa.begin());
    std::copy(b, b + n, fb.begin());

    parallel_NTT(fa.data(), length, p, 1);
    parallel_NTT(fb.data(), length, p, 1);

    #pragma omp parallel for
    for (int i = 0; i < length; ++i)
        fa[i] = static_cast<long long>(fa[i]) * fb[i] % p;

    parallel_NTT(fa.data(), length, p, -1);

    for (int i = 0; i < 2 * n - 1; ++i)
        ab_mod[m][i] = fa[i];
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

    #pragma omp parallel for
    for (int m = 0; m < num_mods; ++m)
        crt_worker(m, a, b, n, length, mods, ab_mod);

    crt_merge(ab_mod, ab, n, p_target);
}

int a[300000], b[300000], ab[300000];

int main(int argc, char *argv[]) {
    int test_begin = 0, test_end = 4;

    for (int i = test_begin; i <= test_end; ++i) {
        long double ans = 0;
        int n_, p_;
        fRead(a, b, &n_, &p_, i);
        memset(ab, 0, sizeof(ab));

        auto Start = std::chrono::high_resolution_clock::now();
        poly_multiply_CRT(a, b, ab, n_, p_);
        auto End = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::ratio<1, 1000>> elapsed = End - Start;
        ans += elapsed.count();
        fCheck(ab, n_, i);
        std::cout << "average latency for n = " << n_ << " p = " << p_ << " : " << ans << " (us) " << std::endl;
        fWrite(ab, n_, i);
    }
    return 0;
}
