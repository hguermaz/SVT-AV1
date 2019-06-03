/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "gtest/gtest.h"
#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbComputeSAD_AVX2.h"
#include "EbComputeSAD_AVX512.h"
#include "EbComputeSAD_C.h"
#include "EbComputeSAD_SSE4_1.h"
#include "EbMeSadCalculation.h"
#include "EbMotionEstimation.h"
#include "EbUnitTest.h"
#include "EbUnitTestUtility.h"

static const int num_test = 28;
static const int num_sad = 6;

struct DistInfo {
    uint32_t width;
    uint32_t height;
};

const struct DistInfo opt_sad_size_info[num_sad] = {
    { 64, 16 },
    { 64, 32 },
    { 64, 64 },
    { 64, 128 },
    { 128, 64 },
    { 128, 128 }
};

typedef uint32_t(*aom_sad_fn_t)(const uint8_t *a, int a_stride,
    const uint8_t *b, int b_stride);

typedef uint32_t(*aom_sad_avg_fn_t)(const uint8_t *a, int a_stride,
    const uint8_t *b, int b_stride,
    const uint8_t *second_pred);

typedef void(*aom_sad_multi_d_fn_t)(const uint8_t *a, int a_stride,
    const uint8_t *const b_array[],
    int b_stride, uint32_t *sad_array);

aom_sad_fn_t aom_sad_func_ptr_array[num_sad][3] = {
    { aom_sad64x16_c, aom_sad64x16_avx2, aom_sad64x16_avx512 },
    { aom_sad64x32_c, aom_sad64x32_avx2, aom_sad64x32_avx512 },
    { aom_sad64x64_c, aom_sad64x64_avx2, aom_sad64x64_avx512 },
    { aom_sad64x128_c, aom_sad64x128_avx2, aom_sad64x128_avx512 },
    { aom_sad128x64_c, aom_sad128x64_avx2, aom_sad128x64_avx512 },
    { aom_sad128x128_c, aom_sad128x128_avx2, aom_sad128x128_avx512 }
};

aom_sad_multi_d_fn_t aom_sad_4d_func_ptr_array[num_sad][3] = {
    { aom_sad64x16x4d_c, aom_sad64x16x4d_avx2, aom_sad64x16x4d_avx2 },
    { aom_sad64x32x4d_c, aom_sad64x32x4d_avx2, aom_sad64x32x4d_avx2 },
    { aom_sad64x64x4d_c, aom_sad64x64x4d_avx2, aom_sad64x64x4d_avx2 },
    { aom_sad64x128x4d_c, aom_sad64x128x4d_avx2, aom_sad64x128x4d_avx2 },
    { aom_sad128x64x4d_c, aom_sad128x64x4d_avx2, aom_sad128x64x4d_avx512 },
    { aom_sad128x128x4d_c, aom_sad128x128x4d_avx2, aom_sad128x128x4d_avx512 }
};

static void init_data_sadMxN(uint8_t **src_ptr, uint32_t *src_stride, uint8_t **ref_ptr,
    uint32_t *ref_stride)
{
    *src_stride = eb_create_random_aligned_stride(MAX_SB_SIZE, 64);
    *ref_stride = eb_create_random_aligned_stride(MAX_SB_SIZE, 64);
    *src_ptr = (uint8_t*)malloc(sizeof(**src_ptr) * MAX_SB_SIZE * *src_stride);
    *ref_ptr = (uint8_t*)malloc(sizeof(**ref_ptr) * MAX_SB_SIZE * *ref_stride);
    eb_buf_random_u8(*src_ptr, MAX_SB_SIZE * *src_stride);
    eb_buf_random_u8(*ref_ptr, MAX_SB_SIZE * *ref_stride);
}

static void init_data_sadMxNx4d(uint8_t **src_ptr, uint32_t *src_stride, uint8_t *ref_ptr[4],
    uint32_t *ref_stride)
{
    *src_stride = eb_create_random_aligned_stride(MAX_SB_SIZE, 64);
    *ref_stride = eb_create_random_aligned_stride(MAX_SB_SIZE, 64);
    *src_ptr = (uint8_t*)malloc(sizeof(**src_ptr) * MAX_SB_SIZE * *src_stride);
    ref_ptr[0] = (uint8_t*)malloc(sizeof(**ref_ptr) * (MAX_SB_SIZE + 3) * *ref_stride);
    eb_buf_random_u8(*src_ptr, MAX_SB_SIZE * *src_stride);
    eb_buf_random_u8(ref_ptr[0], (MAX_SB_SIZE + 3) * *ref_stride);
    ref_ptr[1] = ref_ptr[0] + *ref_stride;
    ref_ptr[2] = ref_ptr[1] + *ref_stride;
    ref_ptr[3] = ref_ptr[2] + *ref_stride;
}

static void uninit_data(uint8_t *src_ptr, uint8_t *ref_ptr)
{
    free(src_ptr);
    free(ref_ptr);
}

TEST(MotionEstimation, sadMxN)
{
    uint8_t *src_ptr, *ref_ptr;
    uint32_t src_stride, ref_stride;

    for (int i = 0; i < 10; i++) {
        init_data_sadMxN(&src_ptr, &src_stride, &ref_ptr, &ref_stride);

        for (int j = 0; j < num_sad; j++) {
            for (int k = 1; k < 3; k++) {
                const uint32_t sad_org = aom_sad_func_ptr_array[j][0](src_ptr, src_stride, ref_ptr, ref_stride);
                const uint32_t sad_opt = aom_sad_func_ptr_array[j][k](src_ptr, src_stride, ref_ptr, ref_stride);

                EXPECT_EQ(sad_org, sad_opt);
            }
        }

        uninit_data(src_ptr, ref_ptr);
    }
}

TEST(MotionEstimation, sadMxNx4d)
{
    uint8_t *src_ptr, *ref_ptr[4];
    uint32_t src_stride, ref_stride;
    uint32_t sad_array_org[4], sad_array_opt[4];

    for (int i = 0; i < 10; i++) {
        init_data_sadMxNx4d(&src_ptr, &src_stride, ref_ptr, &ref_stride);

        for (int j = 0; j < num_sad; j++) {
            for (int k = 1; k < 3; k++) {
                eb_buf_random_u32(sad_array_opt, 4);
                aom_sad_4d_func_ptr_array[j][0](src_ptr, src_stride, ref_ptr, ref_stride, sad_array_org);
                aom_sad_4d_func_ptr_array[j][k](src_ptr, src_stride, ref_ptr, ref_stride, sad_array_opt);

                for (int l = 0; l < 4; l++) {
                    EXPECT_EQ(sad_array_org[l], sad_array_opt[l]);
                }
            }
        }

        uninit_data(src_ptr, ref_ptr[0]);
    }
}

// ===================================

static const int max_ref_stride = 512;

static const struct DistInfo sad_loop_size_info[num_test] = {
    { 4,   2 },
    { 4,   4 },
    { 4,   8 },
    { 4,  16 },
    { 8,   2 },
    { 8,   4 },
    { 8,   8 },
    { 8,  16 },
    { 8,  32 },
    { 16,  4 },
    { 16,  8 },
    { 16, 12 },
    { 16, 16 },
    { 16, 32 },
    { 16, 64 },
    { 24, 16 },
    { 24, 32 },
    { 32,  8 },
    { 32, 16 },
    { 32, 24 },
    { 32, 32 },
    { 32, 64 },
    { 48, 32 },
    { 48, 64 },
    { 64, 16 },
    { 64, 32 },
    { 64, 48 },
    { 64, 64 },
};

static void init_data_sad_loop_kernel(uint8_t *const src, const int32_t size_of_src, uint8_t *const ref,
    const int32_t size_of_ref, const int idx)
{
    if (!idx) {
        memset(src, 0, size_of_src);
        memset(ref, 0, size_of_ref);
    }
    else if (1 == idx) {
        eb_buf_random_u8_to_255(src, size_of_src);
        eb_buf_random_u8_to_255(ref, size_of_ref);
    }
    else if (2 == idx) {
        memset(src, 0, size_of_src);
        eb_buf_random_u8_to_255(ref, size_of_ref);
    }
    else if (3 == idx) {
        eb_buf_random_u8_to_255(src, size_of_src);
        memset(ref, 0, size_of_ref);
    }
    else if (4 == idx) {
        eb_buf_random_u8_to_0_or_255(src, size_of_src);
        eb_buf_random_u8_to_0_or_255(ref, size_of_ref);
    }
    else if (!(idx % 4)) {
        eb_buf_random_u8_to_small(src, size_of_src);
        eb_buf_random_u8_to_large(ref, size_of_ref);
    }
    else if (1 == (idx % 4)) {
        eb_buf_random_u8_to_small_or_large(src, size_of_src);
        eb_buf_random_u8_to_small_or_large(ref, size_of_ref);
    }
    else if (2 == (idx % 4)) {
        const uint32_t range = 16;
        uint8_t val;
        eb_buf_random_u8(&val, 1);
        eb_buf_random_u8_to_near_value(src, size_of_src, val, range);
        eb_buf_random_u8_to_near_value(ref, size_of_ref, val, range);
    }
    else {
        eb_buf_random_u8(src, size_of_src);
        eb_buf_random_u8(ref, size_of_ref);
    }
}

TEST(MotionEstimation, sad_loop_kernel)
{
    const uint32_t src_stride = 256;
    const uint32_t ref_stride = max_ref_stride;
    uint8_t src[MAX_SB_SIZE * src_stride], ref[200 * max_ref_stride];
    uint32_t ref_stride_raw = ref_stride;
    uint64_t best_sad_org, best_sad_opt;
    int16_t x_search_center_org = 0, x_search_center_opt = 0;
    int16_t y_search_center_org = 0, y_search_center_opt = 0;
    const int16_t search_area_height = 64;

    for (int i = 0; i < 10; i++) {
        init_data_sad_loop_kernel(src, sizeof(src), ref, sizeof(ref), i);

        for (int j = 0; j < num_test; j++) {
            const uint32_t width = sad_loop_size_info[j].width;
            const uint32_t height = sad_loop_size_info[j].height;
            for (int16_t search_area_width = 1; search_area_width <= 32; search_area_width++) {
                sad_loop_kernel(src, src_stride, ref, ref_stride, height, width, &best_sad_org, &x_search_center_org, &y_search_center_org, ref_stride_raw, search_area_width, search_area_height);
                sad_loop_kernel_avx512_intrin(src, src_stride, ref, ref_stride, height, width, &best_sad_opt, &x_search_center_opt, &y_search_center_opt, ref_stride_raw, search_area_width, search_area_height);

                if ((best_sad_org != best_sad_opt) || (x_search_center_org != x_search_center_opt) || (y_search_center_org != y_search_center_opt)) {
                    printf("[%d, %d], search_area_width = %d\n", width, height, search_area_width);
                }

                EXPECT_EQ(best_sad_org, best_sad_opt);
                EXPECT_EQ(x_search_center_org, x_search_center_opt);
                EXPECT_EQ(y_search_center_org, y_search_center_opt);
            }
        }
    }
}

#if 0
TEST(MotionEstimation, sad_loop_kernel_speed)
{
    const uint32_t src_stride = 256;
    const uint32_t ref_stride = max_ref_stride;
    uint8_t src[MAX_SB_SIZE * src_stride], ref[200 * max_ref_stride];
    uint32_t ref_stride_raw = ref_stride;
    uint64_t best_sad_org, best_sad_opt;
    int16_t x_search_center_org = 0, x_search_center_opt = 0;
    int16_t y_search_center_org = 0, y_search_center_opt = 0;
    int16_t search_area_width = 64, search_area_height = 64;
    double time_c, time_o;
    uint64_t start_time_seconds, start_time_useconds;
    uint64_t middle_time_seconds, middle_time_useconds;
    uint64_t finish_time_seconds, finish_time_useconds;

    printf("%40s", "sad_loop_kernel\n");

    eb_buf_random_u8(src, sizeof(src));
    eb_buf_random_u8(ref, sizeof(ref));

    for (int j = 0; j < num_test; j++) {
        const uint32_t width = sad_loop_size_info[j].width;
        const uint32_t height = sad_loop_size_info[j].height;
        const uint64_t num_loop = 10000000 / (width + height);

        EbStartTime(&start_time_seconds, &start_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            sad_loop_kernel_avx2_intrin(src, src_stride, ref, ref_stride, height, width, &best_sad_org, &x_search_center_org, &y_search_center_org, ref_stride_raw, search_area_width, search_area_height);
        }

        EbStartTime(&middle_time_seconds, &middle_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            sad_loop_kernel_avx512_intrin(src, src_stride, ref, ref_stride, height, width, &best_sad_opt, &x_search_center_opt, &y_search_center_opt, ref_stride_raw, search_area_width, search_area_height);
        }

        EbStartTime(&finish_time_seconds, &finish_time_useconds);
        EbComputeOverallElapsedTimeMs(start_time_seconds, start_time_useconds,
            middle_time_seconds, middle_time_useconds, &time_c);
        EbComputeOverallElapsedTimeMs(middle_time_seconds, middle_time_useconds,
            finish_time_seconds, finish_time_useconds, &time_o);

        EXPECT_EQ(best_sad_org, best_sad_opt);
        EXPECT_EQ(x_search_center_org, x_search_center_opt);
        EXPECT_EQ(y_search_center_org, y_search_center_opt);

        printf("Average Nanoseconds per Function Call\n");
        printf("    sad_loop_kernel_%2dx%2d_AVX2()   : %6.2f\n",
            width, height, 1000000 * time_c / num_loop);
        printf("    sad_loop_kernel_%2dx%2d_AVX512() : %6.2f   (Comparison: "
            "%5.2fx)\n", width, height, 1000000 * time_o / num_loop, time_c / time_o);
    }
}
#endif

// ===================================

static void init_context(MeContext **const context_ptr, uint8_t *const src, uint8_t *const ref, const uint32_t  src_stride,
const uint32_t  ref_stride, const uint32_t listIndex, const uint32_t refPicIndex, const uint32_t asm_type)
{
    const uint32_t number_of_sb_quad = 1;
    EbErrorType return_error = me_context_ctor(context_ptr);
    (void)return_error;

    (*context_ptr)->sb_src_ptr = src;
    (*context_ptr)->integer_buffer_ptr[listIndex][0] = ref + 200 * ref_stride + 199;
    (*context_ptr)->sb_src_stride = src_stride;
    (*context_ptr)->interpolated_full_stride[listIndex][0] = ref_stride;

    initialize_buffer32bits_func_ptr_array[asm_type]((*context_ptr)->p_sb_best_sad[listIndex][refPicIndex], MAX_ME_PU_COUNT, 1, MAX_SAD_VALUE);
    initialize_buffer32bits_func_ptr_array[asm_type]((*context_ptr)->p_sb_best_mv[listIndex][refPicIndex], MAX_ME_PU_COUNT, 1, MAX_SAD_VALUE);

    (*context_ptr)->p_best_sad64x64 = &((*context_ptr)->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x64]);
    (*context_ptr)->p_best_sad32x32 = &((*context_ptr)->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_0]);
    (*context_ptr)->p_best_sad16x16 = &((*context_ptr)->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_0]);
    (*context_ptr)->p_best_sad8x8 = &((*context_ptr)->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_0]);
    (*context_ptr)->p_best_sad64x32 = &((*context_ptr)->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x32_0]);
    (*context_ptr)->p_best_sad32x16 = &((*context_ptr)->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x16_0]);
    (*context_ptr)->p_best_sad16x8 = &((*context_ptr)->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x8_0]);
    (*context_ptr)->p_best_sad32x64 = &((*context_ptr)->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x64_0]);
    (*context_ptr)->p_best_sad16x32 = &((*context_ptr)->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x32_0]);
    (*context_ptr)->p_best_sad8x16 = &((*context_ptr)->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x16_0]);
    (*context_ptr)->p_best_sad32x8 = &((*context_ptr)->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x8_0]);
    (*context_ptr)->p_best_sad8x32 = &((*context_ptr)->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x32_0]);
    (*context_ptr)->p_best_sad64x16 = &((*context_ptr)->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x16_0]);
    (*context_ptr)->p_best_sad16x64 = &((*context_ptr)->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x64_0]);

    (*context_ptr)->p_best_mv64x64 = &((*context_ptr)->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x64]);
    (*context_ptr)->p_best_mv32x32 = &((*context_ptr)->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_0]);
    (*context_ptr)->p_best_mv16x16 = &((*context_ptr)->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_0]);
    (*context_ptr)->p_best_mv8x8 = &((*context_ptr)->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_0]);
    (*context_ptr)->p_best_mv64x32 = &((*context_ptr)->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x32_0]);
    (*context_ptr)->p_best_mv32x16 = &((*context_ptr)->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x16_0]);
    (*context_ptr)->p_best_mv16x8 = &((*context_ptr)->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x8_0]);
    (*context_ptr)->p_best_mv32x64 = &((*context_ptr)->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x64_0]);
    (*context_ptr)->p_best_mv16x32 = &((*context_ptr)->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x32_0]);
    (*context_ptr)->p_best_mv8x16 = &((*context_ptr)->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x16_0]);
    (*context_ptr)->p_best_mv32x8 = &((*context_ptr)->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x8_0]);
    (*context_ptr)->p_best_mv8x32 = &((*context_ptr)->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x32_0]);
    (*context_ptr)->p_best_mv64x16 = &((*context_ptr)->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x16_0]);
    (*context_ptr)->p_best_mv16x64 = &((*context_ptr)->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x64_0]);

#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH

    (*context_ptr)->p_best_ssd64x64 = &((*context_ptr)->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x64]);
    (*context_ptr)->p_best_ssd32x32 = &((*context_ptr)->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_0]);
    (*context_ptr)->p_best_ssd16x16 = &((*context_ptr)->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_0]);
    (*context_ptr)->p_best_ssd8x8 = &((*context_ptr)->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_0]);
    (*context_ptr)->p_best_ssd64x32 = &((*context_ptr)->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x32_0]);
    (*context_ptr)->p_best_ssd32x16 = &((*context_ptr)->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x16_0]);
    (*context_ptr)->p_best_ssd16x8 = &((*context_ptr)->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x8_0]);
    (*context_ptr)->p_best_ssd32x64 = &((*context_ptr)->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x64_0]);
    (*context_ptr)->p_best_ssd16x32 = &((*context_ptr)->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x32_0]);
    (*context_ptr)->p_best_ssd8x16 = &((*context_ptr)->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x16_0]);
    (*context_ptr)->p_best_ssd32x8 = &((*context_ptr)->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x8_0]);
    (*context_ptr)->p_best_ssd8x32 = &((*context_ptr)->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x32_0]);
    (*context_ptr)->p_best_ssd64x16 = &((*context_ptr)->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x16_0]);
    (*context_ptr)->p_best_ssd16x64 = &((*context_ptr)->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x64_0]);
#endif
}

TEST(MotionEstimation, open_loop_me)
{
    MeContext *context_ptr_org;         // input parameter, ME context Ptr, used to get SB Ptr
    MeContext *context_ptr_opt;         // input parameter, ME context Ptr, used to get SB Ptr
    const uint32_t listIndex = 0;         // input parameter, reference list index
    const uint32_t ref_pic_index = 0;
    const int16_t x_search_area_origin = 0;
    const int16_t y_search_area_origin = 0;
    const uint32_t search_area_width = 23;
    const uint32_t search_area_height = 23;
    const uint32_t refPicIndex = 0;
    uint32_t local_memory_map_index = 0;
    uint64_t local_total_lib_memory;
    uint8_t src[200 * 200];
    uint8_t ref[800 * max_ref_stride];
    const uint32_t src_stride = 200;
    const uint32_t ref_stride = max_ref_stride;

    // Allocate Memory
    local_total_lib_memory = sizeof(EbMemoryMapEntry) * MAX_NUM_PTR;
    memory_map = (EbMemoryMapEntry*)malloc(sizeof(EbMemoryMapEntry) * MAX_NUM_PTR);
    memory_map_index = &local_memory_map_index;
    total_lib_memory = &local_total_lib_memory;

    for (int skip = 0; skip <= 1; skip++) {
        for (int i = 0; i < 10; i++) {
            init_context(&context_ptr_org, src, ref, src_stride, ref_stride, listIndex, refPicIndex, ASM_NON_AVX2);
            init_context(&context_ptr_opt, src, ref, src_stride, ref_stride, listIndex, refPicIndex, ASM_AVX2);
            eb_buf_random_u8(src, sizeof(src));
            eb_buf_random_u8(ref, sizeof(ref));
            context_ptr_org->me_search_method = context_ptr_opt->me_search_method = skip;

            open_loop_me_fullpel_search_sblock(context_ptr_org, listIndex,
#if MRP_ME
                ref_pic_index,
#endif
                x_search_area_origin, y_search_area_origin, search_area_width, search_area_height, ASM_NON_AVX2);
            open_loop_me_fullpel_search_sblock(context_ptr_opt, listIndex,
#if MRP_ME
                ref_pic_index,
#endif
                x_search_area_origin, y_search_area_origin, search_area_width, search_area_height, ASM_AVX2);

            bool result_best_sad = eb_buf_compare_u32(context_ptr_org->p_sb_best_sad[listIndex][refPicIndex], context_ptr_opt->p_sb_best_sad[listIndex][refPicIndex], MAX_ME_PU_COUNT);
            bool result_best_mv = eb_buf_compare_u32(context_ptr_org->p_sb_best_mv[listIndex][refPicIndex], context_ptr_opt->p_sb_best_mv[listIndex][refPicIndex], MAX_ME_PU_COUNT);
            EXPECT_EQ(result_best_sad, true);
            EXPECT_EQ(result_best_mv, true);
        }

        for (int i = 0; i < 10; i++) {
            init_context(&context_ptr_org, src, ref, src_stride, ref_stride, listIndex, refPicIndex, ASM_NON_AVX2);
            init_context(&context_ptr_opt, src, ref, src_stride, ref_stride, listIndex, refPicIndex, ASM_AVX2);
            eb_buf_random_u8(src, sizeof(src));
            eb_buf_random_u8(ref, sizeof(ref));
            context_ptr_org->me_search_method = context_ptr_opt->me_search_method = skip;

            FullPelSearch_LCU(context_ptr_org, listIndex,
#if MRP_ME
                ref_pic_index,
#endif
                x_search_area_origin, y_search_area_origin, search_area_width, search_area_height, ASM_NON_AVX2);
            FullPelSearch_LCU(context_ptr_opt, listIndex,
#if MRP_ME
                ref_pic_index,
#endif
                x_search_area_origin, y_search_area_origin, search_area_width, search_area_height, ASM_AVX2);

            bool result_best_sad = eb_buf_compare_u32(context_ptr_org->p_sb_best_sad[listIndex][refPicIndex], context_ptr_opt->p_sb_best_sad[listIndex][refPicIndex], MAX_ME_PU_COUNT);
            bool result_best_mv = eb_buf_compare_u32(context_ptr_org->p_sb_best_mv[listIndex][refPicIndex], context_ptr_opt->p_sb_best_mv[listIndex][refPicIndex], MAX_ME_PU_COUNT);
            EXPECT_EQ(result_best_sad, true);
            EXPECT_EQ(result_best_mv, true);
        }
    }

    free(memory_map);
}
