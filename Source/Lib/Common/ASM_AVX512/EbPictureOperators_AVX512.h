/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureOperators_AVX512
#define EbPictureOperators_AVX512

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif
    void full_distortion_kernel_cbf_zero32_bits_avx512(
        int32_t  *coeff,
        uint32_t  coeff_stride,
        int32_t  *recon_coeff,
        uint32_t  recon_coeff_stride,
        uint64_t  distortion_result[DIST_CALC_TOTAL],
        uint32_t  area_width,
        uint32_t  area_height);

    void full_distortion_kernel32_bits_avx512(
        int32_t  *coeff,
        uint32_t  coeff_stride,
        int32_t  *recon_coeff,
        uint32_t  recon_coeff_stride,
        uint64_t  distortion_result[DIST_CALC_TOTAL],
        uint32_t  area_width,
        uint32_t  area_height);

    uint64_t spatial_full_distortion_kernel32x_n_avx512_intrin(
        uint8_t   *input,
        uint32_t   input_stride,
        uint8_t   *recon,
        uint32_t   recon_stride,
        uint32_t   area_width,
        uint32_t   area_height);

    uint64_t spatial_full_distortion_kernel64x_n_avx512_intrin(
        uint8_t   *input,
        uint32_t   input_stride,
        uint8_t   *recon,
        uint32_t   recon_stride,
        uint32_t   area_width,
        uint32_t   area_height);

    uint64_t spatial_full_distortion_kernel128x_n_avx512_intrin(
        uint8_t   *input,
        uint32_t   input_stride,
        uint8_t   *recon,
        uint32_t   recon_stride,
        uint32_t   area_width,
        uint32_t   area_height);

#ifdef __cplusplus
}
#endif
#endif // EbPictureOperators_AVX512
