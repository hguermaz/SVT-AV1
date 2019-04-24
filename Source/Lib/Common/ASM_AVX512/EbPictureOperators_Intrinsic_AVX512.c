/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"
#include <immintrin.h>
#include "EbPictureOperators_AVX2.h"
#include "EbPictureOperators_SSE2.h"
#include "EbMemory_AVX2.h"

static INLINE int32_t Hadd32_AVX512_INTRIN(const __m512i src) {
    const __m256i src_L = _mm512_extracti64x4_epi64(src, 0);
    const __m256i src_H = _mm512_extracti64x4_epi64(src, 1);
    const __m256i sum = _mm256_add_epi32(src_L, src_H);

    return Hadd32_AVX2_INTRIN(sum);
}

static INLINE void Distortion_AVX512_INTRIN(const __m256i input,
    const __m256i recon, __m512i *const sum) {
    const __m512i in = _mm512_cvtepu8_epi16(input);
    const __m512i re = _mm512_cvtepu8_epi16(recon);
    const __m512i diff = _mm512_sub_epi16(in, re);
    const __m512i dist = _mm512_madd_epi16(diff, diff);
    *sum = _mm512_add_epi32(*sum, dist);
}

#if 0
// Slightly slower than AVX2 version for small area_height. Disabled.
uint64_t SpatialFullDistortionKernel16xN_AVX512_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *recon,
    uint32_t   recon_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    int32_t row_count = area_height;
    __m512i sum = _mm512_setzero_si512();

    (void)area_width;

    do
    {
        const __m128i in0 = _mm_loadu_si128((__m128i *)(input + 0 * input_stride));
        const __m128i in1 = _mm_loadu_si128((__m128i *)(input + 1 * input_stride));
        const __m128i re0 = _mm_loadu_si128((__m128i *)(recon + 0 * recon_stride));
        const __m128i re1 = _mm_loadu_si128((__m128i *)(recon + 1 * recon_stride));
        const __m256i in = _mm256_setr_m128i(in0, in1);
        const __m256i re = _mm256_setr_m128i(re0, re1);
        Distortion_AVX512_INTRIN(in, re, &sum);
        input += 2 * input_stride;
        recon += 2 * recon_stride;
        row_count -= 2;
    } while (row_count);

    return Hadd32_AVX512_INTRIN(sum);
}
#endif

static INLINE void SpatialFullDistortionKernel64_AVX512_INTRIN(
    const uint8_t *const input, const uint8_t *const recon, __m512i *const sum)
{
    const __m512i in = _mm512_loadu_si512((__m512i *)input);
    const __m512i re = _mm512_loadu_si512((__m512i *)recon);
    const __m512i max = _mm512_max_epu8(in, re);
    const __m512i min = _mm512_min_epu8(in, re);
    const __m512i diff = _mm512_sub_epi8(max, min);
    const __m512i diff_L = _mm512_unpacklo_epi8(diff, _mm512_setzero_si512());
    const __m512i diff_H = _mm512_unpackhi_epi8(diff, _mm512_setzero_si512());
    const __m512i dist_L = _mm512_madd_epi16(diff_L, diff_L);
    const __m512i dist_H = _mm512_madd_epi16(diff_H, diff_H);
    const __m512i dist = _mm512_add_epi32(dist_L, dist_H);
    *sum = _mm512_add_epi32(*sum, dist);
}

uint64_t spatial_full_distortion_kernel32x_n_avx512_intrin(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *recon,
    uint32_t   recon_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    int32_t row_count = area_height;
    __m512i sum = _mm512_setzero_si512();

    (void)area_width;

    do
    {
        const __m256i in = _mm256_loadu_si256((__m256i *)input);
        const __m256i re = _mm256_loadu_si256((__m256i *)recon);
        Distortion_AVX512_INTRIN(in, re, &sum);
        input += input_stride;
        recon += recon_stride;
    } while (--row_count);

    return Hadd32_AVX512_INTRIN(sum);
}

uint64_t spatial_full_distortion_kernel64x_n_avx512_intrin(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *recon,
    uint32_t   recon_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    int32_t row_count = area_height;
    __m512i sum = _mm512_setzero_si512();

    (void)area_width;

    do
    {
        SpatialFullDistortionKernel64_AVX512_INTRIN(input, recon, &sum);
        input += input_stride;
        recon += recon_stride;
    } while (--row_count);

    return Hadd32_AVX512_INTRIN(sum);
}

uint64_t spatial_full_distortion_kernel128x_n_avx512_intrin(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *recon,
    uint32_t   recon_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    int32_t row_count = area_height;
    __m512i sum = _mm512_setzero_si512();

    (void)area_width;

    do
    {
        SpatialFullDistortionKernel64_AVX512_INTRIN(input + 0 * 64, recon + 0 * 64, &sum);
        SpatialFullDistortionKernel64_AVX512_INTRIN(input + 1 * 64, recon + 1 * 64, &sum);
        input += input_stride;
        recon += recon_stride;
    } while (--row_count);

    return Hadd32_AVX512_INTRIN(sum);
}
