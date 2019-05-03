#ifndef AV1_COMMON_X86_AV1_INV_TXFM_AVX512_H_
#define AV1_COMMON_X86_AV1_INV_TXFM_AVX512_H_

#include <immintrin.h>

#ifdef __cplusplus
extern "C" {
#endif
#define NewSqrt2Bits ((int32_t)12)
// 2^12 * sqrt(2)
static const int32_t NewSqrt2 = 5793;
// 2^12 / sqrt(2)
static const int32_t NewInvSqrt2 = 2896;
#ifdef __cplusplus
}
#endif

#endif  // AV1_COMMON_X86_AV1_INV_TXFM_AVX2_H_
