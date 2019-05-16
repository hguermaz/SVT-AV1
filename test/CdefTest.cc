/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "gtest/gtest.h"
#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbCdef.h"
#include "EbUnitTestUtility.h"

static const int num_size = 4;
static const int max_cdef_count = MI_SIZE_128X128 * MI_SIZE_128X128;

const BlockSize size_info[num_size] = {
    BLOCK_4X4,
    BLOCK_4X8,
    BLOCK_8X4,
    BLOCK_8X8
};

const uint32_t max_cdef_counts[num_size] = {
    MI_SIZE_128X128 * MI_SIZE_128X128,
    MI_SIZE_128X128 * MI_SIZE_128X128 / 2,
    MI_SIZE_128X128 * MI_SIZE_128X128 / 2,
    MI_SIZE_128X128 * MI_SIZE_128X128 / 4
};

static void init_data(uint16_t **src, uint16_t **dst, int32_t *const dstride,
    cdef_list *const dlist, const int idx, const uint32_t bd)
{
    *dstride = eb_create_random_aligned_stride(2 * MAX_SB_SIZE, 64);
    *src = (uint16_t*)malloc(sizeof(**src) * MAX_SB_SIZE * MAX_SB_SIZE);
    *dst = (uint16_t*)malloc(sizeof(**dst) * 2 * MAX_SB_SIZE * *dstride);
    if (!idx) {
        eb_buf_random_u16_to_bd(*src, MAX_SB_SIZE * MAX_SB_SIZE, bd);
        eb_buf_random_u16_to_bd(*dst, 2 * MAX_SB_SIZE * *dstride, bd);
    }
    else if (1 == idx) {
        eb_buf_random_u16_to_bd(*src, MAX_SB_SIZE * MAX_SB_SIZE, bd);
        memset(*dst, 0, sizeof(**dst) * 2 * MAX_SB_SIZE * *dstride);
    }
    else if (2 == idx) {
        memset(*src, 0, sizeof(**src) * MAX_SB_SIZE * MAX_SB_SIZE);
        eb_buf_random_u16_to_bd(*dst, 2 * MAX_SB_SIZE * *dstride, bd);
    }
    else if (3 == idx) {
        eb_buf_random_u16_to_0_or_bd(*src, MAX_SB_SIZE * MAX_SB_SIZE, bd);
        eb_buf_random_u16_to_0_or_bd(*dst, 2 * MAX_SB_SIZE * *dstride, bd);
    }
    else {
        eb_buf_random_u16_with_bd(*src, MAX_SB_SIZE * MAX_SB_SIZE, bd);
        eb_buf_random_u16_with_bd(*dst, 2 * MAX_SB_SIZE * *dstride, bd);
    }

    for (int i = 0; i < max_cdef_count; i++) {
        eb_buf_random_u8_with_max(&dlist[i].bx, 1, MAX_SB_SIZE >> 3);
        eb_buf_random_u8_with_max(&dlist[i].by, 1, MAX_SB_SIZE >> 3);
    }
}

static void uninit_data(uint16_t *src, uint16_t *dst)
{
    free(src);
    free(dst);
}

TEST(Cdef, compute_cdef_dist_test)
{
    uint16_t *src, *dst;
    int32_t dstride;
    cdef_list dlist[max_cdef_count];

    for (int i = 0; i < 10; i++) {
        for (uint32_t bd = 8; bd <= 12; bd++) {
            const int32_t coeff_shift = bd - 8;
            init_data(&src, &dst, &dstride, dlist, i, bd);

            for (int32_t pli = 0; pli < 2; pli++) {
                for (int j = 0; j < num_size; j++) {
                    const BlockSize bsize = size_info[j];

                    for (int cdef_count = 0; cdef_count <= max_cdef_counts[j]; cdef_count++) {
                        const uint64_t dist_org = compute_cdef_dist_c(dst, dstride, src, dlist, cdef_count, bsize, coeff_shift, pli);
                        const uint64_t dist_opt = compute_cdef_dist_avx2(dst, dstride, src, dlist, cdef_count, bsize, coeff_shift, pli);

                        EXPECT_EQ(dist_org, dist_opt);
                    }
                }
            }

            uninit_data(src, dst);
        }
    }
}
