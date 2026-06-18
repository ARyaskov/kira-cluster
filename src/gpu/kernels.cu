#include <cuda_runtime.h>
#include <stdint.h>
#include <vector>

extern "C" {

struct CudaAlignmentResult {
    int32_t score;
    uint32_t aligned_len;
    uint32_t matches;
};

__global__ void kernel_hamming_filter(
    const uint8_t* seq_a,
    const uint32_t* a_off,
    const uint32_t* a_len,
    const uint8_t* seq_b,
    const uint32_t* b_off,
    const uint32_t* b_len,
    uint32_t n_pairs,
    uint32_t* out_mismatches
) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_pairs) return;

    uint32_t la = a_len[i];
    uint32_t lb = b_len[i];
    uint32_t n = la < lb ? la : lb;
    const uint8_t* a = seq_a + a_off[i];
    const uint8_t* b = seq_b + b_off[i];

    uint32_t mm = 0;
    for (uint32_t k = 0; k < n; ++k) {
        mm += (a[k] != b[k]);
    }
    out_mismatches[i] = mm;
}

__global__ void kernel_ungapped_filter(
    const uint8_t* seq_a,
    const uint32_t* a_off,
    const uint32_t* a_len,
    const uint8_t* seq_b,
    const uint32_t* b_off,
    const uint32_t* b_len,
    uint32_t n_pairs,
    CudaAlignmentResult* out_results
) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_pairs) return;

    uint32_t la = a_len[i];
    uint32_t lb = b_len[i];
    uint32_t n = la < lb ? la : lb;
    const uint8_t* a = seq_a + a_off[i];
    const uint8_t* b = seq_b + b_off[i];

    uint32_t matches = 0;
    for (uint32_t k = 0; k < n; ++k) {
        matches += (a[k] == b[k]);
    }

    out_results[i].score = (int32_t)matches;
    out_results[i].aligned_len = n;
    out_results[i].matches = matches;
}

__global__ void kernel_sw_gotoh(
    const uint8_t* seq_a,
    const uint32_t* a_off,
    const uint32_t* a_len,
    const uint8_t* seq_b,
    const uint32_t* b_off,
    const uint32_t* b_len,
    uint32_t n_pairs,
    int16_t match_score,
    int16_t mismatch_score,
    int16_t gap_open,
    int16_t gap_extend,
    CudaAlignmentResult* out_results
) {
    const int NEG_INF = -1000000000;
    const int MAX_COLS = 1024;

    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_pairs) return;

    uint32_t la = a_len[i];
    uint32_t lb = b_len[i];
    if (lb + 1 > (uint32_t)MAX_COLS) lb = MAX_COLS - 1;
    if (la > 2048) la = 2048;

    const uint8_t* a = seq_a + a_off[i];
    const uint8_t* b = seq_b + b_off[i];

    int prev_h[MAX_COLS];
    int curr_h[MAX_COLS];
    int prev_f[MAX_COLS];
    int curr_f[MAX_COLS];
    uint32_t prev_l[MAX_COLS], curr_l[MAX_COLS], prev_m[MAX_COLS], curr_m[MAX_COLS];
    uint32_t prev_fl[MAX_COLS], curr_fl[MAX_COLS], prev_fm[MAX_COLS], curr_fm[MAX_COLS];

    for (uint32_t j = 0; j <= lb; ++j) {
        prev_h[j] = 0;
        prev_l[j] = 0;
        prev_m[j] = 0;
        prev_f[j] = NEG_INF;
        prev_fl[j] = 0;
        prev_fm[j] = 0;
    }

    int best_s = 0;
    uint32_t best_l = 0;
    uint32_t best_m = 0;

    for (uint32_t r = 1; r <= la; ++r) {
        curr_h[0] = 0;
        curr_l[0] = 0;
        curr_m[0] = 0;
        curr_f[0] = NEG_INF;
        curr_fl[0] = 0;
        curr_fm[0] = 0;

        int e_s = NEG_INF;
        uint32_t e_l = 0, e_m = 0;

        for (uint32_t c = 1; c <= lb; ++c) {
            int sub = (a[r - 1] == b[c - 1]) ? match_score : mismatch_score;

            int diag_s = prev_h[c - 1] + sub;
            uint32_t diag_l = prev_l[c - 1] + 1;
            uint32_t diag_m = prev_m[c - 1] + ((a[r - 1] == b[c - 1]) ? 1 : 0);

            int e_open_s = curr_h[c - 1] + gap_open + gap_extend;
            uint32_t e_open_l = curr_l[c - 1] + 1;
            uint32_t e_open_m = curr_m[c - 1];

            int e_ext_s = e_s + gap_extend;
            uint32_t e_ext_l = e_l + 1;
            uint32_t e_ext_m = e_m;

            if (e_open_s > e_ext_s || (e_open_s == e_ext_s && e_open_m >= e_ext_m)) {
                e_s = e_open_s;
                e_l = e_open_l;
                e_m = e_open_m;
            } else {
                e_s = e_ext_s;
                e_l = e_ext_l;
                e_m = e_ext_m;
            }

            int f_open_s = prev_h[c] + gap_open + gap_extend;
            uint32_t f_open_l = prev_l[c] + 1;
            uint32_t f_open_m = prev_m[c];

            int f_ext_s = prev_f[c] + gap_extend;
            uint32_t f_ext_l = prev_fl[c] + 1;
            uint32_t f_ext_m = prev_fm[c];

            if (f_open_s > f_ext_s || (f_open_s == f_ext_s && f_open_m >= f_ext_m)) {
                curr_f[c] = f_open_s;
                curr_fl[c] = f_open_l;
                curr_fm[c] = f_open_m;
            } else {
                curr_f[c] = f_ext_s;
                curr_fl[c] = f_ext_l;
                curr_fm[c] = f_ext_m;
            }

            int best_cell_s = diag_s;
            uint32_t best_cell_l = diag_l;
            uint32_t best_cell_m = diag_m;

            if (e_s > best_cell_s || (e_s == best_cell_s && e_m > best_cell_m)) {
                best_cell_s = e_s;
                best_cell_l = e_l;
                best_cell_m = e_m;
            }
            if (curr_f[c] > best_cell_s || (curr_f[c] == best_cell_s && curr_fm[c] > best_cell_m)) {
                best_cell_s = curr_f[c];
                best_cell_l = curr_fl[c];
                best_cell_m = curr_fm[c];
            }

            if (best_cell_s <= 0) {
                curr_h[c] = 0;
                curr_l[c] = 0;
                curr_m[c] = 0;
            } else {
                curr_h[c] = best_cell_s;
                curr_l[c] = best_cell_l;
                curr_m[c] = best_cell_m;
            }

            if (curr_h[c] > best_s || (curr_h[c] == best_s && curr_m[c] > best_m)) {
                best_s = curr_h[c];
                best_l = curr_l[c];
                best_m = curr_m[c];
            }
        }

        for (uint32_t c = 0; c <= lb; ++c) {
            prev_h[c] = curr_h[c];
            prev_l[c] = curr_l[c];
            prev_m[c] = curr_m[c];
            prev_f[c] = curr_f[c];
            prev_fl[c] = curr_fl[c];
            prev_fm[c] = curr_fm[c];
        }
    }

    out_results[i].score = best_s;
    out_results[i].aligned_len = best_l;
    out_results[i].matches = best_m;
}

int kira_cuda_available() {
    int count = 0;
    cudaError_t st = cudaGetDeviceCount(&count);
    if (st != cudaSuccess || count <= 0) {
        return 0;
    }
    return 1;
}

int kira_cuda_sw_batch(
    const uint8_t* const* a_ptrs,
    const uint32_t* a_lens,
    const uint8_t* const* b_ptrs,
    const uint32_t* b_lens,
    uint32_t n_pairs,
    int16_t match_score,
    int16_t mismatch_score,
    int16_t gap_open,
    int16_t gap_extend,
    CudaAlignmentResult* out_results
) {
    if (n_pairs == 0) return 0;

    std::vector<uint32_t> h_a_off(n_pairs), h_b_off(n_pairs);
    uint32_t a_total = 0, b_total = 0;
    for (uint32_t i = 0; i < n_pairs; ++i) {
        h_a_off[i] = a_total;
        h_b_off[i] = b_total;
        a_total += a_lens[i];
        b_total += b_lens[i];
    }

    std::vector<uint8_t> h_a(a_total), h_b(b_total);
    for (uint32_t i = 0; i < n_pairs; ++i) {
        memcpy(h_a.data() + h_a_off[i], a_ptrs[i], a_lens[i]);
        memcpy(h_b.data() + h_b_off[i], b_ptrs[i], b_lens[i]);
    }

    uint8_t *d_a = nullptr, *d_b = nullptr;
    uint32_t *d_a_off = nullptr, *d_b_off = nullptr, *d_a_len = nullptr, *d_b_len = nullptr;
    CudaAlignmentResult* d_out = nullptr;

    cudaError_t st = cudaSuccess;
    st = cudaMalloc(&d_a, a_total);
    if (st != cudaSuccess) return (int)st;
    st = cudaMalloc(&d_b, b_total);
    if (st != cudaSuccess) return (int)st;
    st = cudaMalloc(&d_a_off, sizeof(uint32_t) * n_pairs);
    if (st != cudaSuccess) return (int)st;
    st = cudaMalloc(&d_b_off, sizeof(uint32_t) * n_pairs);
    if (st != cudaSuccess) return (int)st;
    st = cudaMalloc(&d_a_len, sizeof(uint32_t) * n_pairs);
    if (st != cudaSuccess) return (int)st;
    st = cudaMalloc(&d_b_len, sizeof(uint32_t) * n_pairs);
    if (st != cudaSuccess) return (int)st;
    st = cudaMalloc(&d_out, sizeof(CudaAlignmentResult) * n_pairs);
    if (st != cudaSuccess) return (int)st;

    cudaMemcpy(d_a, h_a.data(), a_total, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, h_b.data(), b_total, cudaMemcpyHostToDevice);
    cudaMemcpy(d_a_off, h_a_off.data(), sizeof(uint32_t) * n_pairs, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b_off, h_b_off.data(), sizeof(uint32_t) * n_pairs, cudaMemcpyHostToDevice);
    cudaMemcpy(d_a_len, a_lens, sizeof(uint32_t) * n_pairs, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b_len, b_lens, sizeof(uint32_t) * n_pairs, cudaMemcpyHostToDevice);

    const int block = 64;
    const int grid = (n_pairs + block - 1) / block;
    kernel_sw_gotoh<<<grid, block>>>(
        d_a, d_a_off, d_a_len, d_b, d_b_off, d_b_len, n_pairs,
        match_score, mismatch_score, gap_open, gap_extend, d_out
    );

    st = cudaDeviceSynchronize();
    if (st == cudaSuccess) {
        st = cudaMemcpy(out_results, d_out, sizeof(CudaAlignmentResult) * n_pairs, cudaMemcpyDeviceToHost);
    }

    cudaFree(d_out);
    cudaFree(d_b_len);
    cudaFree(d_a_len);
    cudaFree(d_b_off);
    cudaFree(d_a_off);
    cudaFree(d_b);
    cudaFree(d_a);

    return (int)st;
}

}  // extern "C"
