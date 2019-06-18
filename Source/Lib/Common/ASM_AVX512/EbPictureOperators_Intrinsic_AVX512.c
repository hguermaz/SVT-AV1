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

void full_distortion_kernel_cbf_zero32_bits_avx512(
    int32_t  *coeff,
    uint32_t   coeff_stride,
    int32_t  *recon_coeff,
    uint32_t   recon_coeff_stride,
    uint64_t   distortion_result[DIST_CALC_TOTAL],
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t rowCount, col_count;
    __m512i sum = _mm512_setzero_si512();
    __m256i s1, s2;
    __m128i s3, s4;

    if (area_width == 4) { //when row size is less than 8 pixels
        rowCount = area_height / 2;
        do {
            int32_t *coeffTemp = coeff;
             __m256i x0;
             __m512i y0, z0;

             __m128i in0 = _mm_loadu_si128((__m128i *)(coeffTemp));
             __m128i in1 = _mm_loadu_si128((__m128i *)(coeffTemp + coeff_stride));
             x0 = _mm256_setr_m128i(in0, in1);

             y0 = _mm512_cvtepi32_epi64(x0);
             z0 = _mm512_mul_epi32(y0, y0);
             sum = _mm512_add_epi64(sum, z0);
             coeff += (2 * coeff_stride);
             recon_coeff += (2 * coeff_stride);
             rowCount -= 1;
        } while (rowCount > 0);
    }
    else {
        rowCount = area_height;
        do {
            int32_t *coeffTemp = coeff;

            col_count = area_width / 8;
            do {
                __m256i x0;
                __m512i y0, z0;
                x0 = _mm256_loadu_si256((__m256i *)(coeffTemp));
                coeffTemp += 8;
                y0 = _mm512_cvtepi32_epi64(x0);
                z0 = _mm512_mul_epi32(y0, y0);
                sum = _mm512_add_epi64(sum, z0);
            } while (--col_count);

            coeff += coeff_stride;
            recon_coeff += coeff_stride;
            rowCount -= 1;
        } while (rowCount > 0);
    }

    s1 = _mm512_extracti64x4_epi64(sum, 0);
    s2 = _mm512_extracti64x4_epi64(sum, 1);
    s2 = _mm256_add_epi64(s1, s2);

    s3 = _mm256_extracti128_si256(s2, 0);
    s4 = _mm256_extracti128_si256(s2, 1);

    s3 = _mm_add_epi64(s3, s4);
    s4 = _mm_shuffle_epi32(s3, 0x4e);
    s4 = _mm_add_epi64(s3, s4);

    _mm_storeu_si128((__m128i *)distortion_result, s4);

    (void)recon_coeff_stride;
}

void full_distortion_kernel32_bits_avx512(
    int32_t  *coeff,
    uint32_t   coeff_stride,
    int32_t  *recon_coeff,
    uint32_t   recon_coeff_stride,
    uint64_t   distortion_result[DIST_CALC_TOTAL],
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t rowCount, col_count;
    __m512i sum1 = _mm512_setzero_si512();
    __m512i sum2 = _mm512_setzero_si512();
    __m128i temp1, temp2, temp3;
    __m256i s0, s1;

    if (area_width == 4) { //when number of pixels in row is less than 8
        rowCount = area_height / 2;
        do {
            int32_t *coeffTemp = coeff;
            int32_t *reconCoeffTemp = recon_coeff;

            __m256i x0, y0;
            __m512i x, y, z;

            __m128i in0 = _mm_loadu_si128((__m128i *)(coeffTemp));
            __m128i in1 = _mm_loadu_si128((__m128i *)(coeffTemp + coeff_stride));
            x0 = _mm256_setr_m128i(in0, in1);
            __m128i re0 = _mm_loadu_si128((__m128i *)(reconCoeffTemp));
            __m128i re1 = _mm_loadu_si128((__m128i *)(reconCoeffTemp + recon_coeff_stride));
            y0 = _mm256_setr_m128i(re0, re1);

            x = _mm512_cvtepi32_epi64(x0);
            y = _mm512_cvtepi32_epi64(y0);
            z = _mm512_mul_epi32(x, x);
            sum2 = _mm512_add_epi64(sum2, z);
            x = _mm512_sub_epi64(x, y);
            x = _mm512_mul_epi32(x, x);
            sum1 = _mm512_add_epi64(sum1, x);

            coeff += (2 * coeff_stride);
            recon_coeff += (2 * recon_coeff_stride);
            rowCount -= 1;
        } while (rowCount > 0);
    }
    else {
        rowCount = area_height;
        do {
            int32_t *coeffTemp = coeff;
            int32_t *reconCoeffTemp = recon_coeff;

            col_count = area_width / 8;
            do {
                __m256i x0, y0;
                __m512i x, y, z;
                x0 = _mm256_loadu_si256((__m256i *)(coeffTemp));
                y0 = _mm256_loadu_si256((__m256i *)(reconCoeffTemp));
                x = _mm512_cvtepi32_epi64(x0);
                y = _mm512_cvtepi32_epi64(y0);
                z = _mm512_mul_epi32(x, x);
                sum2 = _mm512_add_epi64(sum2, z);
                x = _mm512_sub_epi64(x, y);
                x = _mm512_mul_epi32(x, x);
                sum1 = _mm512_add_epi64(sum1, x);
                coeffTemp += 8;
                reconCoeffTemp += 8;
            } while (--col_count);

            coeff += coeff_stride;
            recon_coeff += recon_coeff_stride;
            rowCount -= 1;
        } while (rowCount > 0);
    }

    s0 = _mm512_extracti64x4_epi64(sum1, 0);
    s1 = _mm512_extracti64x4_epi64(sum1, 1);
    s0 = _mm256_add_epi64(s0, s1);
    
    temp1 = _mm256_extracti128_si256(s0, 0);
    temp2 = _mm256_extracti128_si256(s0, 1);
    temp1 = _mm_add_epi64(temp1, temp2);
    temp2 = _mm_shuffle_epi32(temp1, 0x4e);
    temp3 = _mm_add_epi64(temp1, temp2);

    s0 = _mm512_extracti64x4_epi64(sum2, 0);
    s1 = _mm512_extracti64x4_epi64(sum2, 1);
    s0 = _mm256_add_epi64(s0, s1);

    temp1 = _mm256_extracti128_si256(s0, 0);
    temp2 = _mm256_extracti128_si256(s0, 1);
    temp1 = _mm_add_epi64(temp1, temp2);
    temp2 = _mm_shuffle_epi32(temp1, 0x4e);
    temp1 = _mm_add_epi64(temp1, temp2);

    temp1 = _mm_unpacklo_epi64(temp3, temp1);
    _mm_storeu_si128((__m128i *)distortion_result, temp1);
}