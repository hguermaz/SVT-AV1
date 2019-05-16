/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "gtest/gtest.h"
#include "EbDefinitions.h"
#include "EbEncHandle.h"
#include "EbPictureOperators.h"
#include "EbUnitTestUtility.h"

static const int num_test = 27;

struct DistInfo {
    uint32_t width;
    uint32_t height;
};

const struct DistInfo size_info[num_test] = {
    { 4,  4},
    { 4,  8},
    { 4, 16},
    { 8,  4},
    { 8,  8},
    { 8, 16},
    { 8, 32},
    {16,  4},
    {16,  8},
    {16, 12},
    {16, 16},
    {16, 32},
    {16, 64},
    {32,  8},
    {32, 16},
    {32, 24},
    {32, 32},
    {32, 64},
    {32, 128},
    {64, 16},
    {64, 32},
    {64, 48},
    {64, 64},
    {64, 128},
    {128, 32},
    {128, 64},
    {128, 128}
};

static void init_data(uint8_t **input, uint32_t *input_stride, uint8_t **recon,
    uint32_t *recon_stride)
{
    *input_stride = eb_create_random_aligned_stride(MAX_SB_SIZE, 64);
    *recon_stride = eb_create_random_aligned_stride(MAX_SB_SIZE, 64);
    *input = (uint8_t*)malloc(sizeof(**input) * MAX_SB_SIZE * *input_stride);
    *recon = (uint8_t*)malloc(sizeof(**recon) * MAX_SB_SIZE * *recon_stride);
    eb_buf_random_u8(*input, MAX_SB_SIZE * *input_stride);
    eb_buf_random_u8(*recon, MAX_SB_SIZE * *recon_stride);
}

static void uninit_data(uint8_t *input, uint8_t *recon)
{
    free(input);
    free(recon);
}

TEST(PictureOperators, spatial_full_distortion)
{
    uint8_t *input, *recon;
    uint32_t input_stride, recon_stride;

    for (int i = 0; i < 10; i++) {
        init_data(&input, &input_stride, &recon, &recon_stride);

        for (int a_type = 0; a_type < ASM_TYPE_TOTAL; a_type++) {
            const EbAsm asm_type = (EbAsm)a_type;

            for (int j = 0; j < num_test; j++) {
                const uint32_t area_width  = size_info[j].width;
                const uint32_t area_height = size_info[j].height;
                const uint64_t dist_org = spatial_full_distortion_kernel(input, input_stride, recon, recon_stride, area_width, area_height);
                const uint64_t dist_opt = spatial_full_distortion_kernel_func_ptr_array[asm_type][Log2f(area_width) - 2](input, input_stride, recon, recon_stride, area_width, area_height);

                EXPECT_EQ(dist_org, dist_opt);
            }
        }

        uninit_data(input, recon);
    }
}

#if 0
TEST(PictureOperators, spatial_full_distortion_speed)
{
    const int bd = 12;
    uint8_t *input, *recon;
    uint32_t input_stride, recon_stride;
    uint32_t area_width, area_height;
    bool result = EB_TRUE;
    double time_c, time_o;
    uint64_t start_time_seconds, start_time_useconds;
    uint64_t middle_time_seconds, middle_time_useconds;
    uint64_t finish_time_seconds, finish_time_useconds;

    printf("%40s", "spatial_full_distortion\n");

    init_data(&input, &input_stride, &recon, &recon_stride);

    for (int j = 0; j < num_test; j++) {
        const uint32_t area_width = size_info[j].width;
        const uint32_t area_height = size_info[j].height;
        const uint64_t num_loop = 1000000000 / (area_width + area_height);
        uint64_t dist_org, dist_opt;

        EbStartTime(&start_time_seconds, &start_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            dist_org = spatial_full_distortion_kernel_func_ptr_array[ASM_AVX2][Log2f(area_width) - 2](input, input_stride, recon, recon_stride, area_width, area_height);
        }

        EbStartTime(&middle_time_seconds, &middle_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            dist_opt = spatial_full_distortion_kernel_func_ptr_array[ASM_AVX512][Log2f(area_width) - 2](input, input_stride, recon, recon_stride, area_width, area_height);
        }

        EbStartTime(&finish_time_seconds, &finish_time_useconds);
        EbComputeOverallElapsedTimeMs(start_time_seconds, start_time_useconds,
            middle_time_seconds, middle_time_useconds, &time_c);
        EbComputeOverallElapsedTimeMs(middle_time_seconds, middle_time_useconds,
            finish_time_seconds, finish_time_useconds, &time_o);

        const bool result = (dist_org == dist_opt);

        EXPECT_EQ(result, true);

        printf("Average Nanoseconds per Function Call\n");
        printf("    SpatialFullDistortionKernel_%2dx%2d_SSE2()   : %6.2f\n",
            area_width, area_height, 1000000 * time_c / num_loop);
        printf("    SpatialFullDistortionKernel_%2dx%2d_AVX2() : %6.2f   (Comparison: "
            "%5.2fx)\n", area_width, area_height, 1000000 * time_o / num_loop, time_c / time_o);
    }
}
#endif
