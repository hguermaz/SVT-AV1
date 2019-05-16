/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "gtest/gtest.h"
#include "EbDefinitions.h"
#include "EbComputeSAD_AVX2.h"
#include "EbComputeSAD_SSE4_1.h"
#include "EbMeSadCalculation.h"
#include "EbMotionEstimation.h"
#include "EbUnitTest.h"
#include "EbUnitTestUtility.h"

static const int max_ref_stride = 800;

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
