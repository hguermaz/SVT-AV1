#include <immintrin.h>
#include "EbHighbdIntraPrediction_SSE2.h"
#include "EbDefinitions.h"
#include "aom_dsp_rtcd.h"

// =============================================================================

// DC RELATED PRED

// Handle number of elements: up to 64.
static INLINE __m128i dc_sum_large(const __m256i src) {
    const __m128i s_lo = _mm256_extracti128_si256(src, 0);
    const __m128i s_hi = _mm256_extracti128_si256(src, 1);
    __m128i sum, sum_hi;
    sum = _mm_add_epi16(s_lo, s_hi);
    sum_hi = _mm_srli_si128(sum, 8);
    sum = _mm_add_epi16(sum, sum_hi);
    // Unpack to avoid 12-bit overflow.
    sum = _mm_unpacklo_epi16(sum, _mm_setzero_si128());

    return dc_sum_4x32bit(sum);
}

static INLINE int dc_common_predictor_32xh_kernel_avx512(uint16_t *dst,
    const ptrdiff_t stride, const int32_t h, const __m512i dc) {
    for (int32_t i = 0; i < h; i++) {
    _mm512_store_si512((__m512i *)dst, dc);
    dst += stride;
    }
}

static INLINE void dc_common_predictor_32xh(uint16_t *const dst,
    const ptrdiff_t stride, const int32_t h, const __m128i dc) {
    const __m512i expected_dc = _mm512_broadcastw_epi16(dc);
    dc_common_predictor_32xh_kernel_avx512(dst, stride, h, expected_dc);
}

static INLINE void dc_common_predictor_64xh_kernel_avx512(uint16_t *dst,
    const ptrdiff_t stride, const int32_t h, const __m512i dc) {
    for (int32_t i = 0; i < h; i++) {
        _mm512_store_si512((__m512i *)(dst + 0x00), dc);
        _mm512_store_si512((__m512i *)(dst + 0x20), dc);
        dst += stride;
    }
}

static INLINE void dc_common_predictor_64xh(uint16_t *const dst,
    const ptrdiff_t stride, const int32_t h, const __m128i dc) {
    const __m512i expected_dc = _mm512_broadcastw_epi16(dc);
    dc_common_predictor_64xh_kernel_avx512(dst, stride, h, expected_dc);
}



static INLINE __m128i dc_sum_16(const uint16_t *const src) {
    const __m256i s = _mm256_loadu_si256((const __m256i *) src);
    const __m128i s_lo = _mm256_extracti128_si256(s, 0);
    const __m128i s_hi = _mm256_extracti128_si256(s, 1);
    const __m128i sum = _mm_add_epi16(s_lo, s_hi);
    return dc_sum_8x16bit(sum);
}

static INLINE __m128i dc_sum_32(const uint16_t *const src) {
    /*const __m256i s0 = _mm256_loadu_si256((const __m256i *)(src + 0x00));
    const __m256i s1 = _mm256_loadu_si256((const __m256i *)(src + 0x10));
    const __m256i sum = _mm256_add_epi16(s0, s1); uncomment only if faster than avx512 code*/

    const __m512i s32 = _mm512_loadu_si512((const __m512i *) src);
    const __m256i s0 = _mm512_extracti64x4_epi64(s32, 0);
    const __m256i s1 = _mm512_extracti64x4_epi64(s32, 1);
    const __m256i sum = _mm256_add_epi16(s0, s1);
    return dc_sum_large(sum);
}

static INLINE __m128i dc_sum_64(const uint16_t *const src) {
    const __m512i s0 = _mm512_loadu_si512((const __m512i *)(src + 0x00));
    const __m512i s1 = _mm512_loadu_si512((const __m512i *)(src + 0x20));
    const __m512i s01 = _mm512_add_epi16(s0, s1);

    const __m256i s2 = _mm512_extracti64x4_epi64(s01, 0);
    const __m256i s3 = _mm512_extracti64x4_epi64(s01, 1);

    const __m256i sum = _mm256_add_epi16(s2, s3);
    return dc_sum_large(sum);
}

// 32xN

void aom_highbd_dc_left_predictor_32x8_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    const __m128i round = _mm_cvtsi32_si128(4);
    __m128i sum;
    (void)above;

    sum = dc_sum_8(left);
    sum = _mm_add_epi16(sum, round);
    sum = _mm_srli_epi16(sum, 3);
    dc_common_predictor_32xh(dst, stride, 8, sum);
}

void aom_highbd_dc_left_predictor_32x16_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    const __m128i round = _mm_cvtsi32_si128(8);
    __m128i sum;
    (void)above;

    sum = dc_sum_16(left);
    sum = _mm_add_epi16(sum, round);
    sum = _mm_srli_epi16(sum, 4);
    dc_common_predictor_32xh(dst, stride, 16, sum);
}

void aom_highbd_dc_left_predictor_32x32_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    const __m128i round = _mm_cvtsi32_si128(16);
    __m128i sum;
    (void)above;

    sum = dc_sum_32(left);
    sum = _mm_add_epi32(sum, round);
    sum = _mm_srli_epi32(sum, 5);
    dc_common_predictor_32xh(dst, stride, 32, sum);
}

void aom_highbd_dc_left_predictor_32x64_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    const __m128i round = _mm_cvtsi32_si128(32);
    __m128i sum;
    (void)above;

    sum = dc_sum_64(left);
    sum = _mm_add_epi32(sum, round);
    sum = _mm_srli_epi32(sum, 6);
    dc_common_predictor_32xh(dst, stride, 64, sum);
}

// 64xN

void aom_highbd_dc_left_predictor_64x16_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    const __m128i round = _mm_cvtsi32_si128(8);
    __m128i sum;
    (void)above;

    sum = dc_sum_16(left);
    sum = _mm_add_epi16(sum, round);
    sum = _mm_srli_epi16(sum, 4);
    dc_common_predictor_64xh(dst, stride, 16, sum);
}

void aom_highbd_dc_left_predictor_64x32_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    const __m128i round = _mm_cvtsi32_si128(16);
    __m128i sum;
    (void)above;

    sum = dc_sum_32(left);
    sum = _mm_add_epi32(sum, round);
    sum = _mm_srli_epi32(sum, 5);
    dc_common_predictor_64xh(dst, stride, 32, sum);
}

void aom_highbd_dc_left_predictor_64x64_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    const __m128i round = _mm_cvtsi32_si128(32);
    __m128i sum;
    (void)above;

    sum = dc_sum_64(left);
    sum = _mm_add_epi32(sum, round);
    sum = _mm_srli_epi32(sum, 6);
    dc_common_predictor_64xh(dst, stride, 64, sum);
}