/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

/***************************************
* Includes
***************************************/
#include "EbRateDistortionCost.h"
#include "aom_dsp_rtcd.h"

#include <assert.h>

#define AV1_COST_PRECISION          0
#define MV_COST_WEIGHT              108

block_size GetBlockSize(uint8_t cu_size) {
    return (cu_size == 64 ? BLOCK_64X64 : cu_size == 32 ? BLOCK_32X32 : cu_size == 16 ? BLOCK_16X16 : cu_size == 8 ? BLOCK_8X8 : BLOCK_4X4);
}

#if ICOPY
int av1_allow_intrabc(const Av1Common *const cm);
int32_t is_chroma_reference(int32_t mi_row, int32_t mi_col, block_size bsize,
#else
static INLINE int32_t is_chroma_reference(int32_t mi_row, int32_t mi_col, block_size bsize,
#endif
    int32_t subsampling_x, int32_t subsampling_y) {
    const int32_t bw = mi_size_wide[bsize];
    const int32_t bh = mi_size_high[bsize];
    int32_t ref_pos = ((mi_row & 0x01) || !(bh & 0x01) || !subsampling_y) &&
        ((mi_col & 0x01) || !(bw & 0x01) || !subsampling_x);
    return ref_pos;
}


uint8_t av1_drl_ctx(const CandidateMv *ref_mv_stack,
    int32_t ref_idx) {

    if (ref_mv_stack[ref_idx].weight >= REF_CAT_LEVEL &&
        ref_mv_stack[ref_idx + 1].weight >= REF_CAT_LEVEL)
        return 0;

    if (ref_mv_stack[ref_idx].weight >= REF_CAT_LEVEL &&
        ref_mv_stack[ref_idx + 1].weight < REF_CAT_LEVEL)
        return 1;

    if (ref_mv_stack[ref_idx].weight < REF_CAT_LEVEL &&
        ref_mv_stack[ref_idx + 1].weight < REF_CAT_LEVEL)
        return 2;

    return 0;
}

/* Symbols for coding which components are zero jointly */
//#define MV_JOINTS 4
//typedef enum {
//    MV_JOINT_ZERO = 0,   /* Zero vector */
//    MV_JOINT_HNZVZ = 1,  /* Vert zero, hor nonzero */
//    MV_JOINT_HZVNZ = 2,  /* Hor zero, vert nonzero */
//    MV_JOINT_HNZVNZ = 3, /* Both components nonzero */
//} MV_JOINT_TYPE;



MV_JOINT_TYPE av1_get_mv_joint(const MV *mv) {
    if (mv->row == 0) {
        return mv->col == 0 ? MV_JOINT_ZERO : MV_JOINT_HNZVZ;
    }
    else {
        return mv->col == 0 ? MV_JOINT_HZVNZ : MV_JOINT_HNZVNZ;
    }
}
int32_t mv_cost(const MV *mv, const int32_t *joint_cost,
    int32_t *const comp_cost[2]) {

    int32_t jnC = av1_get_mv_joint(mv);
    int32_t res =
        joint_cost[jnC] + comp_cost[0][mv->row] +
        comp_cost[1][mv->col];

    return res;
}

int32_t av1_mv_bit_cost(const MV *mv, const MV *ref, const int32_t *mvjcost,
    int32_t *mvcost[2], int32_t weight) {
    const MV diff = { mv->row - ref->row, mv->col - ref->col };
    return ROUND_POWER_OF_TWO(mv_cost(&diff, mvjcost, mvcost) * weight, 7);
}

/////////////////////////////COEFFICIENT CALCULATION //////////////////////////////////////////////
static INLINE int32_t get_golomb_cost(int32_t abs_qc) {
    if (abs_qc >= 1 + NUM_BASE_LEVELS + COEFF_BASE_RANGE) {
        const int32_t r = abs_qc - COEFF_BASE_RANGE - NUM_BASE_LEVELS;
        const int32_t length = get_msb(r) + 1;
        return av1_cost_literal(2 * length - 1);
    }
    return 0;
}

static INLINE int32_t get_txb_bwl(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_wide_log2[tx_size];
}

static INLINE int32_t get_txb_wide(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_wide[tx_size];
}

static INLINE int32_t get_txb_high(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_high[tx_size];
}
static INLINE uint8_t *set_levels(uint8_t *const levels_buf, const int32_t width) {
    return levels_buf + TX_PAD_TOP * (width + TX_PAD_HOR);
}

void av1_txb_init_levels_c(
    const tran_low_t *const coeff,
    const int32_t width,
    const int32_t height,
    uint8_t *const levels) {
    const int32_t stride = width + TX_PAD_HOR;
    uint8_t *ls = levels;

    memset(levels - TX_PAD_TOP * stride, 0,
        sizeof(*levels) * TX_PAD_TOP * stride);
    memset(levels + stride * height, 0,
        sizeof(*levels) * (TX_PAD_BOTTOM * stride + TX_PAD_END));

    for (int32_t i = 0; i < height; i++) {
        for (int32_t j = 0; j < width; j++) {

            *ls++ = (uint8_t)clamp(abs(coeff[i * width + j]), 0, INT8_MAX);

        }
        for (int32_t j = 0; j < TX_PAD_HOR; j++) {
            *ls++ = 0;
        }
    }
}


static const PredictionMode fimode_to_intradir[FILTER_INTRA_MODES] = {
    DC_PRED, V_PRED, H_PRED, D157_PRED, DC_PRED
};
// TODO(angiebird): use this function whenever it's possible
int32_t Av1TransformTypeRateEstimation(
#if CABAC_UP
    uint8_t        allow_update_cdf,
    FRAME_CONTEXT *fc,
#endif
    struct ModeDecisionCandidateBuffer_s    *candidate_buffer_ptr,
    EbBool                                  is_inter,
    EbBool                                  useFilterIntraFlag,
    TxSize                                  transform_size,
    TxType                                  transform_type,
    EbBool                                  reduced_tx_set_used)
{


    uint8_t filterIntraMode = 0; // NM- hardcoded to zero for the moment until we support different intra filtering modes.
    const TxSize square_tx_size = txsize_sqr_map[transform_size];

    //const MbModeInfo *mbmi = &xd->mi[0]->mbmi;
    //const int32_t is_inter = is_inter_block(mbmi);

    if (get_ext_tx_types(transform_size, is_inter, reduced_tx_set_used) > 1  /*&&    !xd->lossless[xd->mi[0]->mbmi.segment_id]  WE ARE NOT LOSSLESS*/) {

        const int32_t ext_tx_set = get_ext_tx_set(transform_size, is_inter, reduced_tx_set_used);
        if (is_inter) {
            if (ext_tx_set > 0)
#if CABAC_UP
            {
                if (allow_update_cdf) {

                    const TxSetType tx_set_type =
                        get_ext_tx_set_type(transform_size, is_inter, reduced_tx_set_used);

                    update_cdf(fc->inter_ext_tx_cdf[ext_tx_set][square_tx_size],
                        av1_ext_tx_ind[tx_set_type][transform_type],
                        av1_num_ext_tx_set[tx_set_type]);
                }
                return candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->interTxTypeFacBits[ext_tx_set][square_tx_size][transform_type];
            }
#else
                return candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->interTxTypeFacBits[ext_tx_set][square_tx_size][transform_type];
#endif
        }
        else {
            if (ext_tx_set > 0) {
                PredictionMode intra_dir;
                if (useFilterIntraFlag)
                    intra_dir = fimode_to_intradir[filterIntraMode];
                else
                    intra_dir = candidate_buffer_ptr->candidate_ptr->pred_mode;
                ASSERT(intra_dir < INTRA_MODES);
#if CABAC_UP
                const TxSetType tx_set_type =
                    get_ext_tx_set_type(transform_size, is_inter, reduced_tx_set_used);

                if (allow_update_cdf) {
                    update_cdf(
                        fc->intra_ext_tx_cdf[ext_tx_set][square_tx_size][intra_dir],
                        av1_ext_tx_ind[tx_set_type][transform_type],
                        av1_num_ext_tx_set[tx_set_type]);
                }
#endif
                return candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->intraTxTypeFacBits[ext_tx_set][square_tx_size][intra_dir][transform_type];
            }
        }
    }
    return 0;
}

const int16_t k_eob_group_start[12] = { 0, 1, 2, 3, 5, 9, 17, 33, 65, 129, 257, 513 };
const int16_t k_eob_offset_bits[12] = { 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

static const int8_t eob_to_pos_small[33] = {
    0, 1, 2,                                        // 0-2
    3, 3,                                           // 3-4
    4, 4, 4, 4,                                     // 5-8
    5, 5, 5, 5, 5, 5, 5, 5,                         // 9-16
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6  // 17-32
};

static const int8_t eob_to_pos_large[17] = {
    6,                               // place holder
    7,                               // 33-64
    8, 8,                           // 65-128
    9, 9, 9, 9,                   // 129-256
    10, 10, 10, 10, 10, 10, 10, 10,  // 257-512
    11                               // 513-
};

static INLINE int32_t get_eob_pos_token(const int32_t eob, int32_t *const extra) {
    int32_t t;

    if (eob < 33) {
        t = eob_to_pos_small[eob];
    }
    else {
        const int32_t e = AOMMIN((eob - 1) >> 5, 16);
        t = eob_to_pos_large[e];
    }

    *extra = eob - k_eob_group_start[t];

    return t;
}
#if CABAC_UP
#define TX_SIZE TxSize
static INLINE TX_SIZE get_txsize_entropy_ctx(TX_SIZE txsize) {
    return (TX_SIZE)((txsize_sqr_map[txsize] + txsize_sqr_up_map[txsize] + 1) >>
        1);
}
void av1_update_eob_context(int eob, TX_SIZE tx_size, TX_CLASS tx_class,
    PLANE_TYPE plane, FRAME_CONTEXT *ec_ctx,
    uint8_t allow_update_cdf) {
    int eob_extra;
    const int eob_pt = get_eob_pos_token(eob, &eob_extra);
    TX_SIZE txs_ctx = get_txsize_entropy_ctx(tx_size);

    const int eob_multi_size = txsize_log2_minus4[tx_size];
    const int eob_multi_ctx = (tx_class == TX_CLASS_2D) ? 0 : 1;


    switch (eob_multi_size) {
    case 0:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi16[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf)
            update_cdf(ec_ctx->eob_flag_cdf16[plane][eob_multi_ctx], eob_pt - 1, 5);
        break;
    case 1:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi32[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf)
            update_cdf(ec_ctx->eob_flag_cdf32[plane][eob_multi_ctx], eob_pt - 1, 6);
        break;
    case 2:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi64[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf)
            update_cdf(ec_ctx->eob_flag_cdf64[plane][eob_multi_ctx], eob_pt - 1, 7);
        break;
    case 3:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi128[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf) {
            update_cdf(ec_ctx->eob_flag_cdf128[plane][eob_multi_ctx], eob_pt - 1,
                8);
        }
        break;
    case 4:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi256[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf) {
            update_cdf(ec_ctx->eob_flag_cdf256[plane][eob_multi_ctx], eob_pt - 1,
                9);
        }
        break;
    case 5:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi512[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf) {
            update_cdf(ec_ctx->eob_flag_cdf512[plane][eob_multi_ctx], eob_pt - 1,
                10);
        }
        break;
    case 6:
    default:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi1024[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf) {
            update_cdf(ec_ctx->eob_flag_cdf1024[plane][eob_multi_ctx], eob_pt - 1,
                11);
        }
        break;
    }

    if (k_eob_offset_bits[eob_pt] > 0) {
        int eob_ctx = eob_pt - 3;
        int eob_shift = k_eob_offset_bits[eob_pt] - 1;
        int bit = (eob_extra & (1 << eob_shift)) ? 1 : 0;
#if CONFIG_ENTROPY_STATS
        counts->eob_extra[cdf_idx][txs_ctx][plane][eob_pt][bit]++;
#endif  // CONFIG_ENTROPY_STATS
        if (allow_update_cdf)
            update_cdf(ec_ctx->eob_extra_cdf[txs_ctx][plane][eob_ctx], bit, 2);
    }
}
#endif
static int32_t get_eob_cost(int32_t eob, const LV_MAP_EOB_COST *txb_eob_costs,
    const LV_MAP_COEFF_COST *txb_costs, TxType tx_type) {
    int32_t eob_extra;
    const int32_t eob_pt = get_eob_pos_token(eob, &eob_extra);
    int32_t eob_cost = 0;
    const int32_t eob_multi_ctx = (tx_type_to_class[tx_type] == TX_CLASS_2D) ? 0 : 1;
    eob_cost = txb_eob_costs->eob_cost[eob_multi_ctx][eob_pt - 1];

    if (k_eob_offset_bits[eob_pt] > 0) {
        const int32_t eob_shift = k_eob_offset_bits[eob_pt] - 1;
        const int32_t bit = (eob_extra & (1 << eob_shift)) ? 1 : 0;
        eob_cost += txb_costs->eob_extra_cost[eob_pt][bit];
        const int32_t offset_bits = k_eob_offset_bits[eob_pt];
        if (offset_bits > 1) eob_cost += av1_cost_literal(offset_bits - 1);
    }
    return eob_cost;
}

// The ctx offset table when TX is TX_CLASS_2D.
// TX col and row indices are clamped to 4.
const int8_t av1_nz_map_ctx_offset[TX_SIZES_ALL][5][5] = {
    // TX_4X4
    { { 0, 1, 6, 6, 0 },
    { 1, 6, 6, 21, 0 },
    { 6, 6, 21, 21, 0 },
    { 6, 21, 21, 21, 0 },
    { 0, 0, 0, 0, 0 } },
    // TX_8X8
    { { 0, 1, 6, 6, 21 },
    { 1, 6, 6, 21, 21 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_16X16
    { { 0, 1, 6, 6, 21 },
    { 1, 6, 6, 21, 21 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_32X32
    { { 0, 1, 6, 6, 21 },
    { 1, 6, 6, 21, 21 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_64X64
    { { 0, 1, 6, 6, 21 },
    { 1, 6, 6, 21, 21 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_4X8
    { { 0, 11, 11, 11, 0 },
    { 11, 11, 11, 11, 0 },
    { 6, 6, 21, 21, 0 },
    { 6, 21, 21, 21, 0 },
    { 21, 21, 21, 21, 0 } },
    // TX_8X4
    { { 0, 16, 6, 6, 21 },
    { 16, 16, 6, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 0, 0, 0, 0, 0 } },
    // TX_8X16
    { { 0, 11, 11, 11, 11 },
    { 11, 11, 11, 11, 11 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_16X8
    { { 0, 16, 6, 6, 21 },
    { 16, 16, 6, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 } },
    // TX_16X32
    { { 0, 11, 11, 11, 11 },
    { 11, 11, 11, 11, 11 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_32X16
    { { 0, 16, 6, 6, 21 },
    { 16, 16, 6, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 } },
    // TX_32X64
    { { 0, 11, 11, 11, 11 },
    { 11, 11, 11, 11, 11 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_64X32
    { { 0, 16, 6, 6, 21 },
    { 16, 16, 6, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 } },
    // TX_4X16
    { { 0, 11, 11, 11, 0 },
    { 11, 11, 11, 11, 0 },
    { 6, 6, 21, 21, 0 },
    { 6, 21, 21, 21, 0 },
    { 21, 21, 21, 21, 0 } },
    // TX_16X4
    { { 0, 16, 6, 6, 21 },
    { 16, 16, 6, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 0, 0, 0, 0, 0 } },
    // TX_8X32
    { { 0, 11, 11, 11, 11 },
    { 11, 11, 11, 11, 11 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_32X8
    { { 0, 16, 6, 6, 21 },
    { 16, 16, 6, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 } },
    // TX_16X64
    { { 0, 11, 11, 11, 11 },
    { 11, 11, 11, 11, 11 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_64X16
    { { 0, 16, 6, 6, 21 },
    { 16, 16, 6, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 } }
};

static INLINE int32_t get_br_ctx(const uint8_t *const levels,
    const int32_t c,  // raster order
    const int32_t bwl, const TxType tx_type) {
    const int32_t row = c >> bwl;
    const int32_t col = c - (row << bwl);
    const int32_t stride = (1 << bwl) + TX_PAD_HOR;
    const TX_CLASS tx_class = tx_type_to_class[tx_type];
    const int32_t pos = row * stride + col;
    int32_t mag = levels[pos + 1];
    mag += levels[pos + stride];
    switch (tx_class) {
    case TX_CLASS_2D:
        mag += levels[pos + stride + 1];
        mag = AOMMIN((mag + 1) >> 1, 6);
        if (c == 0) return mag;
        if ((row < 2) && (col < 2)) return mag + 7;
        break;
    case TX_CLASS_HORIZ:
        mag += levels[pos + 2];
        mag = AOMMIN((mag + 1) >> 1, 6);
        if (c == 0) return mag;
        if (col == 0) return mag + 7;
        break;
    case TX_CLASS_VERT:
        mag += levels[pos + (stride << 1)];
        mag = AOMMIN((mag + 1) >> 1, 6);
        if (c == 0) return mag;
        if (row == 0) return mag + 7;
        break;
    default: break;
    }

    return mag + 14;
}

static INLINE int32_t av1_cost_skip_txb(
#if CABAC_UP
    uint8_t        allow_update_cdf,
    FRAME_CONTEXT *ec_ctx,
#endif
    struct ModeDecisionCandidateBuffer_s    *candidate_buffer_ptr,
    TxSize                                  transform_size,
    PLANE_TYPE                               plane_type,
    int16_t                                   txb_skip_ctx)
{
    const TxSize txs_ctx = (TxSize)((txsize_sqr_map[transform_size] + txsize_sqr_up_map[transform_size] + 1) >> 1);
    const LV_MAP_COEFF_COST *const coeff_costs = &candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->coeffFacBits[txs_ctx][plane_type];
   
#if CABAC_UP   
    if (allow_update_cdf) {
        update_cdf(ec_ctx->txb_skip_cdf[txs_ctx][txb_skip_ctx], 1, 2);
    }
#endif 
    return coeff_costs->txb_skip_cost[txb_skip_ctx][1];
}
// Note: don't call this function when eob is 0.
uint64_t av1_cost_coeffs_txb(
#if CABAC_UP
    uint8_t        allow_update_cdf,
    FRAME_CONTEXT *ec_ctx,
#endif
    struct ModeDecisionCandidateBuffer_s    *candidate_buffer_ptr,
    const tran_low_t                        *const qcoeff,
    uint16_t                                   eob,
    PLANE_TYPE                               plane_type,
    TxSize                                  transform_size,
    /*const uint32_t                             area_size,
    const uint32_t                             stride,*/
    int16_t                                   txb_skip_ctx,
    int16_t                                   dc_sign_ctx,
    EbBool                                  reducedTransformSetFlag)

{

    //Note: there is a different version of this function in AOM that seems to be efficient as its name is:
    //warehouse_efficients_txb

    const TxSize txs_ctx = (TxSize)((txsize_sqr_map[transform_size] + txsize_sqr_up_map[transform_size] + 1) >> 1);
    const TxType transform_type = candidate_buffer_ptr->candidate_ptr->transform_type[plane_type];
    const TX_CLASS tx_class = tx_type_to_class[transform_type];
    int32_t c, cost;
    const int32_t bwl = get_txb_bwl(transform_size);
    const int32_t width = get_txb_wide(transform_size);
    const int32_t height = get_txb_high(transform_size);
    const SCAN_ORDER *const scan_order = &av1_scan_orders[transform_size][transform_type]; // get_scan(tx_size, tx_type);
    const int16_t *const scan = scan_order->scan;
    uint8_t levels_buf[TX_PAD_2D];
    uint8_t *const levels = set_levels(levels_buf, width);
    DECLARE_ALIGNED(16, int8_t, coeff_contexts[MAX_TX_SQUARE]);
    const LV_MAP_COEFF_COST *const coeff_costs = &candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->coeffFacBits[txs_ctx][plane_type];

    const int32_t eob_multi_size = txsize_log2_minus4[transform_size];
    const LV_MAP_EOB_COST *const eobBits = &candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->eobFracBits[eob_multi_size][plane_type];
    // eob must be greater than 0 here.
    assert(eob > 0);
    cost = coeff_costs->txb_skip_cost[txb_skip_ctx][0];

#if CABAC_UP   
    if (allow_update_cdf) {
        update_cdf(ec_ctx->txb_skip_cdf[txs_ctx][txb_skip_ctx], eob == 0, 2);
    }
#endif 
    av1_txb_init_levels(qcoeff, width, height, levels); // NM - Needs to be optimized - to be combined with the quantisation.


    // Transform type bit estimation
    cost += plane_type > PLANE_TYPE_Y ? 0 :
        Av1TransformTypeRateEstimation(
#if CABAC_UP
            allow_update_cdf,
            ec_ctx,
#endif
            candidate_buffer_ptr,
            candidate_buffer_ptr->candidate_ptr->type == INTER_MODE ? EB_TRUE : EB_FALSE,
            EB_FALSE, // NM - Hardcoded to false for the moment until we support the intra filtering
            transform_size,
            transform_type,
            reducedTransformSetFlag);

    // Transform ebo bit estimation
    int32_t eob_cost = get_eob_cost(eob, eobBits, coeff_costs, transform_type);
    cost += eob_cost;
#if CABAC_UP
    if (allow_update_cdf)
        av1_update_eob_context(eob, transform_size, tx_class,
            plane_type, ec_ctx, allow_update_cdf);
#endif
    // Transform non-zero coeff bit estimation
    av1_get_nz_map_contexts(
        levels,
        scan,
        eob,
        transform_size,
        tx_class,
        coeff_contexts); // NM - Assembly version is available in AOM

#if CABAC_UP
    if (allow_update_cdf)
    {
        for (int c = eob - 1; c >= 0; --c) {
            const int pos = scan[c];
            const int coeff_ctx = coeff_contexts[pos];
            const tran_low_t v = qcoeff[pos];
            const tran_low_t level = abs(v);

            if (allow_update_cdf) {
                if (c == eob - 1) {
                    assert(coeff_ctx < 4);
                    update_cdf(
                        ec_ctx->coeff_base_eob_cdf[txs_ctx][plane_type][coeff_ctx],
                        AOMMIN(level, 3) - 1, 3);
                }
                else {
                    update_cdf(ec_ctx->coeff_base_cdf[txs_ctx][plane_type][coeff_ctx],
                        AOMMIN(level, 3), 4);
                }
            }

            {
                if (c == eob - 1) {
                    assert(coeff_ctx < 4);
#if CONFIG_ENTROPY_STATS
                    ++td->counts->coeff_base_eob_multi[cdf_idx][txsize_ctx][plane_type]
                        [coeff_ctx][AOMMIN(level, 3) - 1];
                }
                else {
                    ++td->counts->coeff_base_multi[cdf_idx][txsize_ctx][plane_type]
                        [coeff_ctx][AOMMIN(level, 3)];
#endif
                }
            }

            if (level > NUM_BASE_LEVELS) {
                const int base_range = level - 1 - NUM_BASE_LEVELS;
                const int br_ctx = get_br_ctx(levels, pos, bwl, tx_class);

                for (int idx = 0; idx < COEFF_BASE_RANGE; idx += BR_CDF_SIZE - 1) {

                    const int k = AOMMIN(base_range - idx, BR_CDF_SIZE - 1);
                    if (allow_update_cdf) {
                        update_cdf(ec_ctx->coeff_br_cdf[AOMMIN(txs_ctx, TX_32X32)]
                            [plane_type][br_ctx],
                            k, BR_CDF_SIZE);
                    }
                    for (int lps = 0; lps < BR_CDF_SIZE - 1; lps++) {
#if CONFIG_ENTROPY_STATS
                        ++td->counts->coeff_lps[AOMMIN(txsize_ctx, TX_32X32)][plane_type][lps]
                            [br_ctx][lps == k];
#endif  // CONFIG_ENTROPY_STATS
                        if (lps == k) break;
                    }
#if CONFIG_ENTROPY_STATS
                    ++td->counts->coeff_lps_multi[cdf_idx][AOMMIN(txsize_ctx, TX_32X32)]
                        [plane_type][br_ctx][k];
#endif
                    if (k < BR_CDF_SIZE - 1) break;
                }
            }
        }

        if (qcoeff[0] != 0) {
            const int dc_sign = (qcoeff[0] < 0) ? 1 : 0;
            if (allow_update_cdf)
                update_cdf(ec_ctx->dc_sign_cdf[plane_type][dc_sign_ctx], dc_sign, 2);
        }


        //TODO: CHKN  for 128x128 where we need more than one TXb, we need to update the txb_context(dc_sign+skip_ctx) in a Txb basis.

        return 0;
    }

#endif

    for (c = eob - 1; c >= 0; --c) {

        const int32_t pos = scan[c];
        const tran_low_t v = qcoeff[pos];
        const int32_t is_nz = (v != 0);
        const int32_t level = abs(v);
        const int32_t coeff_ctx = coeff_contexts[pos];

        if (c == eob - 1) {
            ASSERT((AOMMIN(level, 3) - 1) >= 0);
            cost += coeff_costs->base_eob_cost[coeff_ctx][AOMMIN(level, 3) - 1];
        }
        else {
            cost += coeff_costs->base_cost[coeff_ctx][AOMMIN(level, 3)];
        }
        if (is_nz) {
            const int32_t sign = (v < 0) ? 1 : 0;
            // sign bit cost
            if (c == 0) {
                cost += coeff_costs->dc_sign_cost[dc_sign_ctx][sign];
            }
            else {
                cost += av1_cost_literal(1);
            }
            if (level > NUM_BASE_LEVELS) {
                int32_t ctx;
                ctx = get_br_ctx(levels, pos, bwl, transform_type);

                const int32_t base_range = level - 1 - NUM_BASE_LEVELS;
                if (base_range < COEFF_BASE_RANGE) {
                    cost += coeff_costs->lps_cost[ctx][base_range];
                }
                else {
                    cost += coeff_costs->lps_cost[ctx][COEFF_BASE_RANGE];
                }

                if (level >= 1 + NUM_BASE_LEVELS + COEFF_BASE_RANGE) {
                    cost += get_golomb_cost(level);
                }
            }
        }
    }
    return cost;
}
/*static*/ void model_rd_from_sse(
    block_size bsize,
    int16_t quantizer,
    //const AV1_COMP *const cpi,
    //const MacroBlockD *const xd,
    //block_size bsize,
    //int32_t plane,
    uint64_t sse,
    uint32_t *rate,
    uint64_t *dist);

#if REST_FAST_RATE_EST
uint64_t av1_intra_fast_cost(
    CodingUnit_t            *cu_ptr,
    ModeDecisionCandidate_t *candidate_ptr,
    uint32_t                 qp,
    uint64_t                 luma_distortion,
    uint64_t                 chroma_distortion,
    uint64_t                 lambda,
#if USE_SSE_FL
    EbBool                   use_ssd,
#endif
    PictureControlSet_t     *picture_control_set_ptr,
    CandidateMv             *ref_mv_stack,
    const BlockGeom         *blk_geom,
    uint32_t                 miRow,
    uint32_t                 miCol,
    uint32_t                 left_neighbor_mode,
    uint32_t                 top_neighbor_mode)
#else
//static INLINE int32_t av1_get_skip_mode_context(const MacroBlockD *xd) {
//    const MbModeInfo *const above_mi = xd->above_mbmi;
//    const MbModeInfo *const left_mi = xd->left_mbmi;
//    const int32_t above_skip_mode = above_mi ? above_mi->skip_mode : 0;
//    const int32_t left_skip_mode = left_mi ? left_mi->skip_mode : 0;
//    return above_skip_mode + left_skip_mode;
//}
//
//static INLINE int32_t av1_get_skip_context(const MacroBlockD *xd) {
//    const MbModeInfo *const above_mi = xd->above_mbmi;
//    const MbModeInfo *const left_mi = xd->left_mbmi;
//    const int32_t above_skip = above_mi ? above_mi->skip : 0;
//    const int32_t left_skip = left_mi ? left_mi->skip : 0;
//    return above_skip + left_skip;
//}
/*********************************************************************************
* av1_intra_fast_cost function is used to estimate the cost of an intra candidate mode
* for fast mode decisoion module in Intra or inter frame.
* Chroma cost is excluded from fast cost functions. Only the fast_chroma_rate is stored
* for future use in full loop
*
*   @param *cu_ptr(input)
*       cu_ptr is the pointer of the target CU.
*   @param *candidate_buffer_ptr(input)
*       chromaBufferPtr is the buffer pointer of the candidate luma mode.
*   @param qp(input)
*       qp is the quantizer parameter.
*   @param luma_distortion (input)
*       luma_distortion is the intra condidate luma distortion.
*   @param lambda(input)
*       lambda is the Lagrange multiplier
**********************************************************************************/
EbErrorType av1_intra_fast_cost(
    struct ModeDecisionContext_s            *context_ptr,
    CodingUnit_t                            *cu_ptr,
#if REST_FAST_RATE_EST
    ModeDecisionCandidate_t                 *candidate_ptr,
#else
    struct ModeDecisionCandidateBuffer_s    *candidate_buffer_ptr,
#endif
    uint32_t                                  qp,
    uint64_t                                  luma_distortion,
    uint64_t                                  chroma_distortion,
    uint64_t                                  lambda,
    PictureControlSet_t                     *picture_control_set_ptr)
#endif
{
#if !REST_FAST_RATE_EST
    EbErrorType return_error = EB_ErrorNone;

    (void)qp;
    (void)picture_control_set_ptr;
#endif
    UNUSED(qp);
    UNUSED(ref_mv_stack);
    UNUSED(miRow);
    UNUSED(miCol);
    UNUSED(left_neighbor_mode);
    UNUSED(top_neighbor_mode);

#if ICOPY 
   
    if (av1_allow_intrabc(picture_control_set_ptr->parent_pcs_ptr->av1_cm) && candidate_ptr->use_intrabc) {

        uint64_t lumaSad = (LUMA_WEIGHT * luma_distortion) << AV1_COST_PRECISION;
        uint64_t chromaSad = chroma_distortion << AV1_COST_PRECISION;
        uint64_t totalDistortion = lumaSad + chromaSad;

        uint64_t rate = 0;

        EbReflist refListIdx = 0;
        int16_t predRefX = candidate_ptr->motion_vector_pred_x[refListIdx];
        int16_t predRefY = candidate_ptr->motion_vector_pred_y[refListIdx];
        int16_t mvRefX = candidate_ptr->motionVector_x_L0;
        int16_t mvRefY = candidate_ptr->motionVector_y_L0;
        MV mv;
        mv.row = mvRefY;
        mv.col = mvRefX;
        MV ref_mv;
        ref_mv.row = predRefY;
        ref_mv.col = predRefX;

        int *dvcost[2] = { (int *)&candidate_ptr->md_rate_estimation_ptr->dv_cost[0][MV_MAX],
                           (int *)&candidate_ptr->md_rate_estimation_ptr->dv_cost[1][MV_MAX] };

        int32_t mvRate = av1_mv_bit_cost(
            &mv,
            &ref_mv,
            candidate_ptr->md_rate_estimation_ptr->dv_joint_cost,
            dvcost, MV_COST_WEIGHT_SUB);

        rate = mvRate + candidate_ptr->md_rate_estimation_ptr->intrabcFacBits[candidate_ptr->use_intrabc];

        candidate_ptr->fast_luma_rate = rate;
        candidate_ptr->fast_chroma_rate = 0;

        lumaSad = (LUMA_WEIGHT * luma_distortion) << AV1_COST_PRECISION;
        chromaSad = chroma_distortion << AV1_COST_PRECISION;
        totalDistortion = lumaSad + chromaSad;

       
        return(RDCOST(lambda, rate, totalDistortion));

    }
    else {
#endif

    EbBool isMonochromeFlag = EB_FALSE; // NM - isMonochromeFlag is harcoded to false.
#if REST_FAST_RATE_EST
    EbBool isCflAllowed = (blk_geom->bwidth <= 32 && blk_geom->bheight <= 32) ? 1 : 0;
#else
    EbBool isCflAllowed = (context_ptr->blk_geom->bwidth <= 32 &&
        context_ptr->blk_geom->bheight <= 32) ? 1 : 0;
#endif

    uint8_t   subSamplingX = 1; // NM - subsampling_x is harcoded to 1 for 420 chroma sampling.
    uint8_t   subSamplingY = 1; // NM - subsampling_y is harcoded to 1 for 420 chroma sampling.
#if !REST_FAST_RATE_EST
    block_size cuSizeIndex = context_ptr->blk_geom->bsize;
    uint32_t miRow = context_ptr->cu_origin_y >> MI_SIZE_LOG2;
    uint32_t miCol = context_ptr->cu_origin_x >> MI_SIZE_LOG2;
#endif
    // In fast loop CFL alphas are not know yet. The chroma mode bits are calculated based on DC Mode, and if CFL is the winner compared to CFL, ChromaBits are updated
#if REST_FAST_RATE_EST
    uint32_t chroma_mode = candidate_ptr->intra_chroma_mode == UV_CFL_PRED ? UV_DC_PRED : candidate_ptr->intra_chroma_mode;
#else
    uint32_t chroma_mode = candidate_buffer_ptr->candidate_ptr->intra_chroma_mode == UV_CFL_PRED ? UV_DC_PRED : candidate_buffer_ptr->candidate_ptr->intra_chroma_mode;
#endif

    // Number of bits for each synatax element
    uint64_t intraModeBitsNum = 0;
    uint64_t intraLumaModeBitsNum = 0;
    uint64_t intraLumaAngModeBitsNum = 0;
    uint64_t intraChromaModeBitsNum = 0;
    uint64_t intraChromaAngModeBitsNum = 0;
    uint64_t skipModeRate = 0;
    uint8_t  skipModeCtx = cu_ptr->skip_flag_context; // NM - Harcoded to 1 until the skip_mode context is added.
#if REST_FAST_RATE_EST
    PredictionMode intra_mode = (PredictionMode)candidate_ptr->pred_mode;
#else
    PredictionMode intra_mode = (PredictionMode)candidate_buffer_ptr->candidate_ptr->pred_mode;
#endif
    // Luma and chroma rate
    uint32_t rate;
    uint32_t lumaRate = 0;
    uint32_t chromaRate = 0;
    uint64_t lumaSad, chromaSad;

    // Luma and chroma distortion
    uint64_t totalDistortion;
#if !REST_FAST_RATE_EST
    uint32_t left_neighbor_mode = context_ptr->intra_luma_left_mode;
    uint32_t top_neighbor_mode = context_ptr->intra_luma_top_mode;
#endif
    const int32_t AboveCtx = intra_mode_context[top_neighbor_mode];
    const int32_t LeftCtx = intra_mode_context[left_neighbor_mode];
#if REST_FAST_RATE_EST
    intraModeBitsNum = picture_control_set_ptr->slice_type != I_SLICE ? (uint64_t)candidate_ptr->md_rate_estimation_ptr->mbModeFacBits[size_group_lookup[blk_geom->bsize]][intra_mode] : ZERO_COST;
    skipModeRate = picture_control_set_ptr->slice_type != I_SLICE ? (uint64_t)candidate_ptr->md_rate_estimation_ptr->skipModeFacBits[skipModeCtx][0] : ZERO_COST;
#else
    intraModeBitsNum = picture_control_set_ptr->slice_type != I_SLICE ? (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->mbModeFacBits[size_group_lookup[cuSizeIndex]][intra_mode] : ZERO_COST;
    skipModeRate = picture_control_set_ptr->slice_type != I_SLICE ? (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->skipModeFacBits[skipModeCtx][0] : ZERO_COST;
#endif

    // Estimate luma nominal intra mode bits
#if REST_FAST_RATE_EST
    intraLumaModeBitsNum = picture_control_set_ptr->slice_type == I_SLICE ? (uint64_t)candidate_ptr->md_rate_estimation_ptr->yModeFacBits[AboveCtx][LeftCtx][intra_mode] : ZERO_COST;
#else
    intraLumaModeBitsNum = picture_control_set_ptr->slice_type == I_SLICE ? (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->yModeFacBits[AboveCtx][LeftCtx][intra_mode] : ZERO_COST;
#endif
    // Estimate luma angular mode bits
#if REST_FAST_RATE_EST
    if (candidate_ptr->is_directional_mode_flag && candidate_ptr->use_angle_delta) {
#else
    if (candidate_buffer_ptr->candidate_ptr->is_directional_mode_flag && candidate_buffer_ptr->candidate_ptr->use_angle_delta) {
#endif
        ASSERT((intra_mode - V_PRED) < 8);
        ASSERT((intra_mode - V_PRED) >= 0);
#if REST_FAST_RATE_EST
        intraLumaAngModeBitsNum = candidate_ptr->md_rate_estimation_ptr->angleDeltaFacBits[intra_mode - V_PRED][MAX_ANGLE_DELTA + candidate_ptr->angle_delta[PLANE_TYPE_Y]];
#else
        intraLumaAngModeBitsNum = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->angleDeltaFacBits[intra_mode - V_PRED][MAX_ANGLE_DELTA + candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y]];
#endif
    }


    // NM- Harcoded assuming luma mode is equal to chroma mode
    //if (!cm->seq_params.monochrome &&
    //    is_chroma_reference(mi_row, mi_col, bsize, xd->plane[1].subsampling_x,
    //    xd->plane[1].subsampling_y)) {
    //    mbmi->uv_mode =
    //        read_intra_mode_uv(ec_ctx, r, is_cfl_allowed(xd), mbmi->mode);
    //    if (mbmi->uv_mode == UV_CFL_PRED) {
    //        mbmi->cfl_alpha_idx =
    //            read_cfl_alphas(xd->tile_ctx, r, &mbmi->cfl_alpha_signs);
    //        xd->cfl.store_y = 1;
    //    }
    //    else {
    //        xd->cfl.store_y = 0;
    //    }
    //    mbmi->angle_delta[PLANE_TYPE_UV] =
    //        use_angle_delta && av1_is_directional_mode(get_uv_mode(mbmi->uv_mode))
    //        ? read_angle_delta(r,
    //        ec_ctx->angle_delta_cdf[mbmi->uv_mode - V_PRED])
    //        : 0;
    //}
    //else {
    //    // Avoid decoding angle_info if there is is no chroma prediction
    //    mbmi->uv_mode = UV_DC_PRED;
    //    xd->cfl.is_chroma_reference = 0;
    //    xd->cfl.store_y = 1;
    //}
#if REST_FAST_RATE_EST
    if (blk_geom->has_uv) {
        if (!isMonochromeFlag && is_chroma_reference(miRow, miCol, blk_geom->bsize, subSamplingX, subSamplingY)) {
#else
    if (context_ptr->blk_geom->has_uv) {
        if (!isMonochromeFlag && is_chroma_reference(miRow, miCol, cuSizeIndex, subSamplingX, subSamplingY)) {
#endif
            // Estimate luma nominal intra mode bits
#if REST_FAST_RATE_EST
            intraChromaModeBitsNum = (uint64_t)candidate_ptr->md_rate_estimation_ptr->intraUVmodeFacBits[isCflAllowed][intra_mode][chroma_mode];
#else
            intraChromaModeBitsNum = (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->intraUVmodeFacBits[isCflAllowed][intra_mode][chroma_mode];
#endif
            // Estimate luma angular mode bits
#if REST_FAST_RATE_EST
            if (candidate_ptr->is_directional_chroma_mode_flag && candidate_ptr->use_angle_delta) {
                intraChromaAngModeBitsNum = candidate_ptr->md_rate_estimation_ptr->angleDeltaFacBits[chroma_mode - V_PRED][MAX_ANGLE_DELTA + candidate_ptr->angle_delta[PLANE_TYPE_UV]];
#else
            if (candidate_buffer_ptr->candidate_ptr->is_directional_chroma_mode_flag && candidate_buffer_ptr->candidate_ptr->use_angle_delta) {
                intraChromaAngModeBitsNum = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->angleDeltaFacBits[chroma_mode - V_PRED][MAX_ANGLE_DELTA + candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_UV]];
#endif
            }
        }
    }

#if REST_FAST_RATE_EST
    uint32_t isInterRate = picture_control_set_ptr->slice_type != I_SLICE ? candidate_ptr->md_rate_estimation_ptr->intraInterFacBits[cu_ptr->is_inter_ctx][0] : 0;
#else
    uint32_t isInterRate = picture_control_set_ptr->slice_type != I_SLICE ? candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->intraInterFacBits[cu_ptr->is_inter_ctx][0] : 0;
#endif
    lumaRate = intraModeBitsNum + skipModeRate + intraLumaModeBitsNum + intraLumaAngModeBitsNum + isInterRate;
#if ICOPY
    if (av1_allow_intrabc(picture_control_set_ptr->parent_pcs_ptr->av1_cm))
        lumaRate += candidate_ptr->md_rate_estimation_ptr->intrabcFacBits[candidate_ptr->use_intrabc];
#endif

    chromaRate = intraChromaModeBitsNum + intraChromaAngModeBitsNum;

    // Keep the Fast Luma and Chroma rate for future use
#if REST_FAST_RATE_EST
    candidate_ptr->fast_luma_rate = lumaRate;
    candidate_ptr->fast_chroma_rate = chromaRate;
#else
    candidate_buffer_ptr->candidate_ptr->fast_luma_rate = lumaRate;
    candidate_buffer_ptr->candidate_ptr->fast_chroma_rate = chromaRate;
#endif
#if USE_SSE_FL // cost
    if (use_ssd) {

        int32_t current_q_index = MAX(0, MIN(QINDEX_RANGE - 1, picture_control_set_ptr->parent_pcs_ptr->base_qindex));
        Dequants *const dequants = &picture_control_set_ptr->parent_pcs_ptr->deq;

        int16_t quantizer = dequants->y_dequant_Q3[current_q_index][1];
        rate = 0;
        model_rd_from_sse(
            blk_geom->bsize,
            quantizer,
            luma_distortion,
            &rate,
            &lumaSad);
        lumaRate += rate;
        totalDistortion = lumaSad;

        rate = 0;
        model_rd_from_sse(
            blk_geom->bsize_uv,
            quantizer,
            chroma_distortion,
            &chromaRate,
            &chromaSad);
        chromaRate += rate;
        totalDistortion += chromaSad;

        rate = lumaRate + chromaRate;

        return(RDCOST(lambda, rate, totalDistortion));
    }
    else {
#endif
        lumaSad = (LUMA_WEIGHT * luma_distortion) << AV1_COST_PRECISION;
        chromaSad = chroma_distortion << AV1_COST_PRECISION;
        totalDistortion = lumaSad + chromaSad;

        rate = lumaRate + chromaRate;

        // Assign fast cost
#if REST_FAST_RATE_EST
        return(RDCOST(lambda, rate, totalDistortion));
#else
        *(candidate_buffer_ptr->fast_cost_ptr) = RDCOST(lambda, rate, totalDistortion);

        return return_error;
#endif
#if USE_SSE_FL
    }
#endif
#if ICOPY
    }
#endif
}

//extern INLINE int32_t have_newmv_in_inter_mode(PredictionMode mode);
static INLINE int32_t have_newmv_in_inter_mode(PredictionMode mode) {
    return (mode == NEWMV || mode == NEW_NEWMV || mode == NEAREST_NEWMV ||
        mode == NEW_NEARESTMV || mode == NEAR_NEWMV || mode == NEW_NEARMV);
}

extern void av1_set_ref_frame(MvReferenceFrame *rf,
    int8_t ref_frame_type);

// This function encodes the reference frame
uint64_t EstimateRefFramesNumBits(
    PictureControlSet_t                    *picture_control_set_ptr,
#if REST_FAST_RATE_EST
    ModeDecisionCandidate_t                *candidate_ptr,
#else
    ModeDecisionCandidateBuffer_t          *candidate_buffer_ptr,
#endif
    CodingUnit_t                           *cu_ptr,
    uint32_t                                 bwidth,
    uint32_t                                 bheight,
    uint8_t                                  ref_frame_type,
    EbBool                                is_compound)
{

    uint64_t refRateBits = 0;
    uint64_t refRateA = 0;
    uint64_t refRateB = 0;
    uint64_t refRateC = 0;
    uint64_t refRateD = 0;
    uint64_t refRateE = 0;
    uint64_t refRateF = 0;
    uint64_t refRateG = 0;
    uint64_t refRateH = 0;
    uint64_t refRateI = 0;
    uint64_t refRateJ = 0;
    uint64_t refRateK = 0;
    uint64_t refRateL = 0;
    uint64_t refRateM = 0;

    // If segment level coding of this signal is disabled...
    // or the segment allows multiple reference frame options
    /*if (segfeature_active(&cm->seg, segment_id, SEG_LVL_REF_FRAME)) {
    assert(!is_compound);
    assert(mbmi->ref_frame[0] ==
    get_segdata(&cm->seg, segment_id, SEG_LVL_REF_FRAME));
    }
    else if (segfeature_active(&cm->seg, segment_id, SEG_LVL_SKIP) ||
    segfeature_active(&cm->seg, segment_id, SEG_LVL_GLOBALMV)) {
    assert(!is_compound);
    assert(mbmi->ref_frame[0] == LAST_FRAME);
    }
    else*/ {
    // does the feature use compound prediction or not
    // (if not specified at the frame/segment level)
        if (picture_control_set_ptr->parent_pcs_ptr->reference_mode == REFERENCE_MODE_SELECT) {
            if (MIN(bwidth, bheight) >= 8) {
                int32_t context = 0;
                context = cu_ptr->reference_mode_context;
                assert(context >= 0 && context < 5);
#if REST_FAST_RATE_EST
                refRateA = candidate_ptr->md_rate_estimation_ptr->compInterFacBits[context][is_compound];
#else
                refRateA = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->compInterFacBits[context][is_compound];
#endif

            }
        }
        else {
            assert((!is_compound) == (picture_control_set_ptr->parent_pcs_ptr->reference_mode == SINGLE_REFERENCE));
        }
        int32_t context = 0;
        if (is_compound) {
            const COMP_REFERENCE_TYPE comp_ref_type = /*has_uni_comp_refs(mbmi)
                                                      ? UNIDIR_COMP_REFERENCE
                                                      : */BIDIR_COMP_REFERENCE;
            MvReferenceFrame refType[2];
            av1_set_ref_frame(refType, ref_frame_type);


            context = cu_ptr->compoud_reference_type_context;
            assert(context >= 0 && context < 5);
#if REST_FAST_RATE_EST
            refRateB = candidate_ptr->md_rate_estimation_ptr->compRefTypeFacBits[context][comp_ref_type];
#else
            refRateB = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->compRefTypeFacBits[context][comp_ref_type];
#endif

            if (comp_ref_type == UNIDIR_COMP_REFERENCE) {
                printf("ERROR[AN]: UNIDIR_COMP_REFERENCE not supported\n");
                //const int32_t bit = mbmi->ref_frame[0] == BWDREF_FRAME;
                //WRITE_REF_BIT(bit, uni_comp_ref_p);

                //if (!bit) {
                //    assert(mbmi->ref_frame[0] == LAST_FRAME);
                //    const int32_t bit1 = mbmi->ref_frame[1] == LAST3_FRAME ||
                //        mbmi->ref_frame[1] == GOLDEN_FRAME;
                //    WRITE_REF_BIT(bit1, uni_comp_ref_p1);
                //    if (bit1) {
                //        const int32_t bit2 = mbmi->ref_frame[1] == GOLDEN_FRAME;
                //        WRITE_REF_BIT(bit2, uni_comp_ref_p2);
                //    }
                //}
                //else {
                //    assert(mbmi->ref_frame[1] == ALTREF_FRAME);
                //}

                //return;
            }

            assert(comp_ref_type == BIDIR_COMP_REFERENCE);

            const int32_t bit = (refType[0] == GOLDEN_FRAME ||
                refType[0] == LAST3_FRAME);

            context = av1_get_pred_context_comp_ref_p(cu_ptr->av1xd);
            assert(context >= 0 && context < 3);
#if REST_FAST_RATE_EST
            refRateC = candidate_ptr->md_rate_estimation_ptr->compRefFacBits[context][0][bit];
#else
            refRateC = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->compRefFacBits[context][0][bit];
#endif
            //            WRITE_REF_BIT(bit, comp_ref_p);

            if (!bit) {
                const int32_t bit1 = (refType[0] == LAST2_FRAME);
                context = av1_get_pred_context_comp_ref_p1(cu_ptr->av1xd);
                /*aom_write_symbol(ecWriter, bit1, frameContext->comp_ref_cdf[context][1],
                    2);*/
                assert(context >= 0 && context < 3);
#if REST_FAST_RATE_EST
                refRateD = candidate_ptr->md_rate_estimation_ptr->compRefFacBits[context][1][bit1];
#else
                refRateD = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->compRefFacBits[context][1][bit1];
#endif

                //WRITE_REF_BIT(bit1, comp_ref_p1);
            }
            else {
                const int32_t bit2 = (refType[0] == GOLDEN_FRAME);
                context = av1_get_pred_context_comp_ref_p2(cu_ptr->av1xd);
                /*aom_write_symbol(ecWriter, bit2, frameContext->comp_ref_cdf[context][2],
                    2);*/
                assert(context >= 0 && context < 3);
#if REST_FAST_RATE_EST
                refRateE = candidate_ptr->md_rate_estimation_ptr->compRefFacBits[context][2][bit2];
#else
                refRateE = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->compRefFacBits[context][2][bit2];
#endif

                //WRITE_REF_BIT(bit2, comp_ref_p2);
            }

            const int32_t bit_bwd = (refType[1] == ALTREF_FRAME);
            context = av1_get_pred_context_comp_bwdref_p(cu_ptr->av1xd);
            /*aom_write_symbol(ecWriter, bit_bwd, frameContext->comp_bwdref_cdf[context][0],
                2);*/
            assert(context >= 0 && context < 3);
#if REST_FAST_RATE_EST
            refRateF = candidate_ptr->md_rate_estimation_ptr->compBwdRefFacBits[context][0][bit_bwd];
#else
            refRateF = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->compBwdRefFacBits[context][0][bit_bwd];
#endif
            //WRITE_REF_BIT(bit_bwd, comp_bwdref_p);


            if (!bit_bwd) {
                context = av1_get_pred_context_comp_bwdref_p1(cu_ptr->av1xd);
                /*aom_write_symbol(ecWriter, refType[1] == ALTREF2_FRAME, frameContext->comp_bwdref_cdf[context][1],
                    2);*/
                assert(context >= 0 && context < 3);
#if REST_FAST_RATE_EST
                refRateG = candidate_ptr->md_rate_estimation_ptr->compBwdRefFacBits[context][1][refType[1] == ALTREF2_FRAME];
#else
                refRateG = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->compBwdRefFacBits[context][1][refType[1] == ALTREF2_FRAME];
#endif
                //WRITE_REF_BIT(mbmi->ref_frame[1] == ALTREF2_FRAME, comp_bwdref_p1);

            }

        }
        else {
            const int32_t bit0 = (ref_frame_type <= ALTREF_FRAME &&
                ref_frame_type >= BWDREF_FRAME);//0

            context = av1_get_pred_context_single_ref_p1(cu_ptr->av1xd);
            /*aom_write_symbol(ecWriter, bit0, frameContext->single_ref_cdf[context][0],
                2);*/
            assert(context >= 0 && context < 3);
#if REST_FAST_RATE_EST
            refRateH = candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][0][bit0];
#else
            refRateH = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][0][bit0];
#endif
            //WRITE_REF_BIT(bit0, single_ref_p1);


            if (bit0) {
                const int32_t bit1 = (ref_frame_type == ALTREF_FRAME);
                context = av1_get_pred_context_single_ref_p2(cu_ptr->av1xd);
                assert(context >= 0 && context < 3);
                /*aom_write_symbol(ecWriter, bit1, frameContext->single_ref_cdf[context][1],
                    2);*/
#if REST_FAST_RATE_EST
                refRateI = candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][1][bit1];
#else
                refRateI = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][1][bit1];
#endif
                //WRITE_REF_BIT(bit1, single_ref_p2);


                if (!bit1) {
                    context = av1_get_pred_context_single_ref_p6(cu_ptr->av1xd);
                    /*aom_write_symbol(ecWriter, cu_ptr->prediction_unit_array[0].ref_frame_type == ALTREF2_FRAME, frameContext->single_ref_cdf[context][5],
                        2);*/
                    assert(context >= 0 && context < 3);
#if REST_FAST_RATE_EST
                    refRateJ = candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][5][ref_frame_type == ALTREF2_FRAME];
#else
                    refRateJ = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][5][ref_frame_type == ALTREF2_FRAME];
#endif
                    //WRITE_REF_BIT(mbmi->ref_frame[0] == ALTREF2_FRAME, single_ref_p6);

                }
            }
            else {
                const int32_t bit2 = (ref_frame_type == LAST3_FRAME ||
                    ref_frame_type == GOLDEN_FRAME); //0
                context = av1_get_pred_context_single_ref_p3(cu_ptr->av1xd);
                /*aom_write_symbol(ecWriter, bit2, frameContext->single_ref_cdf[context][2],
                    2);*/
                assert(context >= 0 && context < 3);
#if REST_FAST_RATE_EST
                refRateK = candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][2][bit2];
#else
                refRateK = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][2][bit2];
#endif
                //WRITE_REF_BIT(bit2, single_ref_p3);

                if (!bit2) {
                    const int32_t bit3 = (ref_frame_type != LAST_FRAME); //0;
                    context = av1_get_pred_context_single_ref_p4(cu_ptr->av1xd);
                    assert(context >= 0 && context < 3);
                    /*aom_write_symbol(ecWriter, bit3, frameContext->single_ref_cdf[context][3],
                        2);*/
#if REST_FAST_RATE_EST
                    refRateL = candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][3][bit3];
#else
                    refRateL = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][3][bit3];
#endif
                    //WRITE_REF_BIT(bit3, single_ref_p4);

                }
                else {
                    const int32_t bit4 = (ref_frame_type != LAST3_FRAME);
                    context = av1_get_pred_context_single_ref_p5(cu_ptr->av1xd);
                    /*aom_write_symbol(ecWriter, bit4, frameContext->single_ref_cdf[context][4],
                        2);*/
                    assert(context >= 0 && context < 3);
#if REST_FAST_RATE_EST
                    refRateM = candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][4][bit4];
#else
                    refRateM = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][4][bit4];
#endif
                    //WRITE_REF_BIT(bit4, single_ref_p5);

                }
            }
        }
    }

    refRateBits = refRateA + refRateB + refRateC + refRateD + refRateE + refRateF + refRateG + refRateH + refRateI + refRateJ + refRateK + refRateL + refRateM;
    return refRateBits;
}
//extern INLINE int16_t Av1ModeContextAnalyzer(const int16_t *const mode_context, const MvReferenceFrame *const rf);

extern  int8_t av1_ref_frame_type(const MvReferenceFrame *const rf);
uint16_t compound_mode_ctx_map_2[3][COMP_NEWMV_CTXS] = {
   { 0, 1, 1, 1, 1 },
   { 1, 2, 3, 4, 4 },
   { 4, 4, 5, 6, 7 },
};
static INLINE int16_t Av1ModeContextAnalyzer(
    const int16_t *const mode_context, const MvReferenceFrame *const rf) {
    const int8_t ref_frame = av1_ref_frame_type(rf);

    if (rf[1] <= INTRA_FRAME) return mode_context[ref_frame];

    const int16_t newmv_ctx = mode_context[ref_frame] & NEWMV_CTX_MASK;
    const int16_t refmv_ctx =
        (mode_context[ref_frame] >> REFMV_OFFSET) & REFMV_CTX_MASK;
    ASSERT((refmv_ctx >> 1) < 3);
    const int16_t comp_ctx = compound_mode_ctx_map_2[refmv_ctx >> 1][AOMMIN(
        newmv_ctx, COMP_NEWMV_CTXS - 1)];
    return comp_ctx;
}

#if REST_FAST_RATE_EST
uint64_t av1_inter_fast_cost(  
    CodingUnit_t            *cu_ptr,
    ModeDecisionCandidate_t *candidate_ptr,
    uint32_t                 qp,
    uint64_t                 luma_distortion,
    uint64_t                 chroma_distortion,
    uint64_t                 lambda,
#if USE_SSE_FL
    EbBool                   use_ssd,
#endif
    PictureControlSet_t     *picture_control_set_ptr,
    CandidateMv             *ref_mv_stack,
    const BlockGeom         *blk_geom,
    uint32_t                 miRow,
    uint32_t                 miCol,
    uint32_t                 left_neighbor_mode,
    uint32_t                 top_neighbor_mode)
#else
/*********************************************************************************
* av1_inter_fast_cost function is used to estimate the cost of an inter candidate mode
* for fast mode decisoion module in Inter frame.
*
*   @param *cu_ptr(input)
*       cu_ptr is the pointer of the target CU.
*   @param *candidate_buffer_ptr(input)
*       chromaBufferPtr is the buffer pointer of the candidate luma mode.
*   @param qp(input)
*       qp is the quantizer parameter.
*   @param luma_distortion (input)
*       luma_distortion is the inter condidate luma distortion.
*   @param lambda(input)
*       lambda is the Lagrange multiplier
**********************************************************************************/
EbErrorType av1_inter_fast_cost(
    struct ModeDecisionContext_s           *context_ptr,
    CodingUnit_t                           *cu_ptr,
    ModeDecisionCandidateBuffer_t          *candidate_buffer_ptr,
    uint32_t                                  qp,
    uint64_t                                  luma_distortion,
    uint64_t                                  chroma_distortion,
    uint64_t                                  lambda,
    PictureControlSet_t                    *picture_control_set_ptr)
#endif
{
    UNUSED(top_neighbor_mode);
    UNUSED(left_neighbor_mode);
    UNUSED(miCol);
    UNUSED(miRow);

#if !REST_FAST_RATE_EST    
    EbErrorType  return_error = EB_ErrorNone;

    ModeDecisionCandidate_t *candidate_ptr = candidate_buffer_ptr->candidate_ptr;
#endif
    // Luma rate
    uint32_t           lumaRate = 0;
    uint32_t           chromaRate = 0;
    uint64_t           mvRate = 0;
    uint64_t           skipModeRate;
    // Luma and chroma distortion
    uint64_t           lumaSad;
    uint64_t             chromaSad;
    uint64_t           totalDistortion;

    uint32_t           rate;

    int16_t           predRefX;
    int16_t           predRefY;
    int16_t           mvRefX;
    int16_t           mvRefY;

    EbReflist       refListIdx;

    (void)qp;

    PredictionMode inter_mode = (PredictionMode)candidate_ptr->pred_mode;

    uint64_t interModeBitsNum = 0;

    uint8_t skipModeCtx = cu_ptr->skip_flag_context;
    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, candidate_ptr->ref_frame_type);
    uint32_t modeCtx = Av1ModeContextAnalyzer(cu_ptr->inter_mode_ctx, rf);
#if REST_FAST_RATE_EST
    skipModeRate = candidate_ptr->md_rate_estimation_ptr->skipModeFacBits[skipModeCtx][0];
#else
    skipModeRate = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->skipModeFacBits[skipModeCtx][0];
#endif
    uint64_t referencePictureBitsNum = 0;

    //Reference Type and Mode Bit estimation

    referencePictureBitsNum = EstimateRefFramesNumBits(
        picture_control_set_ptr,
#if REST_FAST_RATE_EST
        candidate_ptr,
#else
        candidate_buffer_ptr,
#endif
        cu_ptr,
#if REST_FAST_RATE_EST
        blk_geom->bwidth,
        blk_geom->bheight,
#else
        context_ptr->blk_geom->bwidth,
        context_ptr->blk_geom->bheight,
#endif
        candidate_ptr->ref_frame_type,
        candidate_ptr->is_compound);


    if (candidate_ptr->is_compound) {
#if REST_FAST_RATE_EST
        interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->interCompoundModeFacBits[modeCtx][INTER_COMPOUND_OFFSET(inter_mode)];
#else
        interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->interCompoundModeFacBits[modeCtx][INTER_COMPOUND_OFFSET(inter_mode)];
#endif
    }
    else {
        //uint32_t newmv_ctx = modeCtx & NEWMV_CTX_MASK;
        //interModeBitsNum = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->newMvModeFacBits[mode_ctx][0];

        int16_t newmv_ctx = modeCtx & NEWMV_CTX_MASK;
        //aom_write_symbol(ecWriter, mode != NEWMV, frameContext->newmv_cdf[newmv_ctx], 2);
#if REST_FAST_RATE_EST
        interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->newMvModeFacBits[newmv_ctx][inter_mode != NEWMV];
#else
        interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->newMvModeFacBits[newmv_ctx][inter_mode != NEWMV];
#endif
        if (inter_mode != NEWMV) {
            const int16_t zeromvCtx = (modeCtx >> GLOBALMV_OFFSET) & GLOBALMV_CTX_MASK;
            //aom_write_symbol(ecWriter, mode != GLOBALMV, frameContext->zeromv_cdf[zeromvCtx], 2);
#if REST_FAST_RATE_EST
            interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->zeroMvModeFacBits[zeromvCtx][inter_mode != GLOBALMV];
#else
            interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->zeroMvModeFacBits[zeromvCtx][inter_mode != GLOBALMV];
#endif
            if (inter_mode != GLOBALMV) {
                int16_t refmvCtx = (modeCtx >> REFMV_OFFSET) & REFMV_CTX_MASK;
                /*aom_write_symbol(ecWriter, mode != NEARESTMV, frameContext->refmv_cdf[refmv_ctx], 2);*/
#if REST_FAST_RATE_EST
                interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->refMvModeFacBits[refmvCtx][inter_mode != NEARESTMV];
#else
                interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->refMvModeFacBits[refmvCtx][inter_mode != NEARESTMV];
#endif
            }
        }
    }
    if (inter_mode == NEWMV || inter_mode == NEW_NEWMV || have_nearmv_in_inter_mode(inter_mode)) {

        //drLIdex cost estimation
        const int32_t new_mv = inter_mode == NEWMV || inter_mode == NEW_NEWMV;
        if (new_mv) {
            int32_t idx;
            for (idx = 0; idx < 2; ++idx) {
                if (cu_ptr->av1xd->ref_mv_count[candidate_ptr->ref_frame_type] > idx + 1) {
#if REST_FAST_RATE_EST                    
                    uint8_t drl1Ctx =
                        av1_drl_ctx(ref_mv_stack, idx);
                    interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->drlModeFacBits[drl1Ctx][candidate_ptr->drl_index != idx];
#else
                    uint8_t drl1Ctx =
                        av1_drl_ctx(context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[candidate_ptr->ref_frame_type], idx);
                    interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->drlModeFacBits[drl1Ctx][candidate_ptr->drl_index != idx];
#endif
                    if (candidate_ptr->drl_index == idx) break;
                }
            }
        }

        if (have_nearmv_in_inter_mode(inter_mode)) {
            int32_t idx;
            // TODO(jingning): Temporary solution to compensate the NEARESTMV offset.
            for (idx = 1; idx < 3; ++idx) {
                if (cu_ptr->av1xd->ref_mv_count[candidate_ptr->ref_frame_type] > idx + 1) {
#if REST_FAST_RATE_EST
                    uint8_t drl_ctx =
                        av1_drl_ctx(ref_mv_stack, idx);
                    interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->drlModeFacBits[drl_ctx][candidate_ptr->drl_index != (idx - 1)];
#else
                    uint8_t drl_ctx =
                        av1_drl_ctx(context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[candidate_ptr->ref_frame_type], idx);
                    interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->drlModeFacBits[drl_ctx][candidate_ptr->drl_index != (idx - 1)];
#endif

                    if (candidate_ptr->drl_index == (idx - 1)) break;
                }
            }
        }

    }

    if (have_newmv_in_inter_mode(inter_mode)) {
        if (candidate_ptr->is_compound) {

            mvRate = 0;

            if (inter_mode == NEW_NEWMV) {


                for (refListIdx = 0; refListIdx < 2; ++refListIdx) {

                    predRefX = candidate_ptr->motion_vector_pred_x[refListIdx];
                    predRefY = candidate_ptr->motion_vector_pred_y[refListIdx];
                    mvRefX = refListIdx == REF_LIST_1 ? candidate_ptr->motionVector_x_L1 : candidate_ptr->motionVector_x_L0;
                    mvRefY = refListIdx == REF_LIST_1 ? candidate_ptr->motionVector_y_L1 : candidate_ptr->motionVector_y_L0;


                    MV mv;
                    mv.row = mvRefY;
                    mv.col = mvRefX;

                    MV ref_mv;
                    ref_mv.row = predRefY;
                    ref_mv.col = predRefX;

                    mvRate += av1_mv_bit_cost(
                        &mv,
                        &ref_mv,
#if REST_FAST_RATE_EST
                        candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                        candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
#else
                        candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                        candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
#endif
                        MV_COST_WEIGHT);
                }

            }
            else if (inter_mode == NEAREST_NEWMV || inter_mode == NEAR_NEWMV) {

                predRefX = candidate_ptr->motion_vector_pred_x[REF_LIST_1];
                predRefY = candidate_ptr->motion_vector_pred_y[REF_LIST_1];
                mvRefX = candidate_ptr->motionVector_x_L1;
                mvRefY = candidate_ptr->motionVector_y_L1;


                MV mv;
                mv.row = mvRefY;
                mv.col = mvRefX;

                MV ref_mv;
                ref_mv.row = predRefY;
                ref_mv.col = predRefX;

                mvRate += av1_mv_bit_cost(
                    &mv,
                    &ref_mv,
#if REST_FAST_RATE_EST
                    candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                    candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
#else
                    candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                    candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
#endif
                    MV_COST_WEIGHT);

            }
            else {
                assert(inter_mode == NEW_NEARESTMV || inter_mode == NEW_NEARMV);

                predRefX = candidate_ptr->motion_vector_pred_x[REF_LIST_0];
                predRefY = candidate_ptr->motion_vector_pred_y[REF_LIST_0];
                mvRefX = candidate_ptr->motionVector_x_L0;
                mvRefY = candidate_ptr->motionVector_y_L0;

                MV mv;
                mv.row = mvRefY;
                mv.col = mvRefX;

                MV ref_mv;
                ref_mv.row = predRefY;
                ref_mv.col = predRefX;

                mvRate += av1_mv_bit_cost(
                    &mv,
                    &ref_mv,
#if REST_FAST_RATE_EST
                    candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                    candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
#else
                    candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                    candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
#endif
                    MV_COST_WEIGHT);


            }
        }
        else {

            refListIdx = candidate_ptr->prediction_direction[0] == 0 ? 0 : 1;

            predRefX = candidate_ptr->motion_vector_pred_x[refListIdx];
            predRefY = candidate_ptr->motion_vector_pred_y[refListIdx];

            mvRefX = refListIdx == 0 ? candidate_ptr->motionVector_x_L0 : candidate_ptr->motionVector_x_L1;
            mvRefY = refListIdx == 0 ? candidate_ptr->motionVector_y_L0 : candidate_ptr->motionVector_y_L1;

            MV mv;
            mv.row = mvRefY;
            mv.col = mvRefX;

            MV ref_mv;
            ref_mv.row = predRefY;
            ref_mv.col = predRefX;

            mvRate = av1_mv_bit_cost(
                &mv,
                &ref_mv,
#if REST_FAST_RATE_EST
                candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
#else
                candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
#endif
                MV_COST_WEIGHT);
        }
    }

    // NM - To be added when the intrainter mode is adopted
    //  read_interintra_mode(is_compound)

    EbBool is_inter = inter_mode >= SINGLE_INTER_MODE_START && inter_mode < SINGLE_INTER_MODE_END;
    if (is_inter
        && picture_control_set_ptr->parent_pcs_ptr->switchable_motion_mode
        && rf[1] != INTRA_FRAME)
    {

#if REST_FAST_RATE_EST
        MOTION_MODE motion_mode_rd = candidate_ptr->motion_mode;
        block_size bsize = blk_geom->bsize;
        cu_ptr->prediction_unit_array[0].num_proj_ref = candidate_ptr->num_proj_ref;
#else
        MOTION_MODE motion_mode_rd = candidate_buffer_ptr->candidate_ptr->motion_mode;
        block_size bsize = context_ptr->blk_geom->bsize;

        cu_ptr->prediction_unit_array[0].num_proj_ref = candidate_buffer_ptr->candidate_ptr->num_proj_ref;
#endif
        MOTION_MODE last_motion_mode_allowed = motion_mode_allowed(
            picture_control_set_ptr,
            cu_ptr,
            bsize,
            rf[0],
            rf[1],
            inter_mode);

        switch (last_motion_mode_allowed) {
        case SIMPLE_TRANSLATION: break;
        case OBMC_CAUSAL:
            assert(motion_mode_rd == SIMPLE_TRANSLATION); // TODO: remove when OBMC added
#if REST_FAST_RATE_EST
            interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->motionModeFacBits1[bsize][motion_mode_rd];
#else
            interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->motionModeFacBits1[bsize][motion_mode_rd];
#endif
            break;
        default:
#if REST_FAST_RATE_EST
            interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->motionModeFacBits[bsize][motion_mode_rd];
#else
            interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->motionModeFacBits[bsize][motion_mode_rd];
#endif
        }
    }

    // NM - To be added when the overlappable mode is adopted
    //    read_compound_type(is_compound)
    // NM - To be added when switchable filter is adopted
    //    if (interpolation_filter == SWITCHABLE) {
    //        for (dir = 0; dir < (enable_dual_filter ? 2 : 1); dir++) {
    //            if (needs_interp_filter()) {
    //                interp_filter[dir]    S()
    //            }
    //            else {
    //                interp_filter[dir] = EIGHTTAP
    //            }
    //        }
    //        if (!enable_dual_filter)
    //            interp_filter[1] = interp_filter[0]
    //    }
    //    else {
    //        for (dir = 0; dir < 2; dir++)
    //            interp_filter[dir] = interpolation_filter
    //    }
#if REST_FAST_RATE_EST
    uint32_t isInterRate = candidate_ptr->md_rate_estimation_ptr->intraInterFacBits[cu_ptr->is_inter_ctx][1];
#else
    uint32_t isInterRate = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->intraInterFacBits[cu_ptr->is_inter_ctx][1];
#endif
    lumaRate = referencePictureBitsNum + skipModeRate + interModeBitsNum + mvRate + isInterRate;


    //chromaRate = intraChromaModeBitsNum + intraChromaAngModeBitsNum;

    // Keep the Fast Luma and Chroma rate for future use
#if REST_FAST_RATE_EST
    candidate_ptr->fast_luma_rate = lumaRate;
    candidate_ptr->fast_chroma_rate = chromaRate;
#else 
    candidate_buffer_ptr->candidate_ptr->fast_luma_rate = lumaRate;
    candidate_buffer_ptr->candidate_ptr->fast_chroma_rate = chromaRate;
#endif

#if USE_SSE_FL // cost
    if (use_ssd) {

        int32_t current_q_index = MAX(0, MIN(QINDEX_RANGE - 1, picture_control_set_ptr->parent_pcs_ptr->base_qindex));
        Dequants *const dequants = &picture_control_set_ptr->parent_pcs_ptr->deq;

        int16_t quantizer = dequants->y_dequant_Q3[current_q_index][1];
        rate = 0;
        model_rd_from_sse(
            blk_geom->bsize,
            quantizer,
            luma_distortion,
            &rate,
            &lumaSad);
        lumaRate += rate;
        totalDistortion = lumaSad;

        rate = 0;
        model_rd_from_sse(
            blk_geom->bsize_uv,
            quantizer,
            chroma_distortion,
            &chromaRate,
            &chromaSad);
        chromaRate += rate;
        totalDistortion += chromaSad;

        rate = lumaRate + chromaRate;

        if (candidate_ptr->merge_flag) {
            uint64_t skipModeRate = candidate_ptr->md_rate_estimation_ptr->skipModeFacBits[skipModeCtx][1];
            if (skipModeRate < rate) {
                return(RDCOST(lambda, skipModeRate, totalDistortion));
            }
        }
        return(RDCOST(lambda, rate, totalDistortion));
    }
    else {
#endif
        lumaSad = (LUMA_WEIGHT * luma_distortion) << AV1_COST_PRECISION;
        chromaSad = chroma_distortion << AV1_COST_PRECISION;
        totalDistortion = lumaSad + chromaSad;
#if REST_FAST_RATE_EST
        if (blk_geom->has_uv == 0 && chromaSad != 0) {
#else
        if (context_ptr->blk_geom->has_uv == 0 && chromaSad != 0) {
#endif
            printf("av1_inter_fast_cost: Chroma error");
        }


        rate = lumaRate + chromaRate;


        // Assign fast cost
#if REST_FAST_RATE_EST
        if (candidate_ptr->merge_flag) {
            uint64_t skipModeRate = candidate_ptr->md_rate_estimation_ptr->skipModeFacBits[skipModeCtx][1];
            if (skipModeRate < rate) {
                return(RDCOST(lambda, skipModeRate, totalDistortion));
            }
        }
        return(RDCOST(lambda, rate, totalDistortion));

#else
        *(candidate_buffer_ptr->fast_cost_ptr) = RDCOST(lambda, rate, totalDistortion);

        if (candidate_buffer_ptr->candidate_ptr->merge_flag) {
            uint64_t skipModeRate = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->skipModeFacBits[skipModeCtx][1];
            if (skipModeRate < rate) {

                *(candidate_buffer_ptr->fast_cost_ptr) = RDCOST(lambda, skipModeRate, totalDistortion);
            }

        }


        return return_error;
#endif
#if USE_SSE_FL
    }
#endif
}


EbErrorType Av1TuEstimateCoeffBits(
#if CABAC_UP
    uint8_t        allow_update_cdf,
    FRAME_CONTEXT *ec_ctx,
#endif
    PictureControlSet_t                    *picture_control_set_ptr,
    struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
    CodingUnit_t                           *cu_ptr,
    uint32_t                                  tuOriginIndex,
    uint32_t                                  tuChromaOriginIndex,
    EntropyCoder_t                         *entropy_coder_ptr,
    EbPictureBufferDesc_t                  *coeff_buffer_sb,
    uint32_t                                 yEob,
    uint32_t                                 cbEob,
    uint32_t                                 crEob,
    uint64_t                                 *y_tu_coeff_bits,
    uint64_t                                 *cb_tu_coeff_bits,
    uint64_t                                 *cr_tu_coeff_bits,
    TxSize                                 txsize,
    TxSize                                 txsize_uv,
    COMPONENT_TYPE                          component_type,
    EbAsm                                  asm_type)
{
    (void)asm_type;
    (void)entropy_coder_ptr;
    EbErrorType return_error = EB_ErrorNone;


    int32_t *coeff_buffer;


    int16_t  luma_txb_skip_context = cu_ptr->luma_txb_skip_context;
    int16_t  luma_dc_sign_context = cu_ptr->luma_dc_sign_context;
    int16_t  cb_txb_skip_context = cu_ptr->cb_txb_skip_context;
    int16_t  cb_dc_sign_context = cu_ptr->cb_dc_sign_context;
    int16_t  cr_txb_skip_context = cu_ptr->cr_txb_skip_context;
    int16_t  cr_dc_sign_context = cu_ptr->cr_dc_sign_context;


    EbBool reducedTransformSetFlag = picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used ? EB_TRUE : EB_FALSE;

    //Estimate the rate of the transform type and coefficient for Luma

    if (component_type == COMPONENT_LUMA || component_type == COMPONENT_ALL) {
        if (yEob) {
            coeff_buffer = (int32_t*)&coeff_buffer_sb->buffer_y[tuOriginIndex * sizeof(int32_t)];

            *y_tu_coeff_bits = av1_cost_coeffs_txb(
#if CABAC_UP
                allow_update_cdf,
                ec_ctx,
#endif
                candidate_buffer_ptr,
                coeff_buffer,
                (uint16_t)yEob,
                PLANE_TYPE_Y,
                txsize,
                luma_txb_skip_context,
                luma_dc_sign_context,
                reducedTransformSetFlag);
        }
        else {
            *y_tu_coeff_bits = av1_cost_skip_txb(
#if CABAC_UP
                allow_update_cdf,
                ec_ctx,
#endif
                candidate_buffer_ptr,
                txsize,
                PLANE_TYPE_Y,
                luma_txb_skip_context);
        }
    }
    //Estimate the rate of the transform type and coefficient for chroma Cb

    if (component_type == COMPONENT_CHROMA_CB || component_type == COMPONENT_CHROMA || component_type == COMPONENT_ALL) {

        if (cbEob) {

            coeff_buffer = (int32_t*)&coeff_buffer_sb->bufferCb[tuChromaOriginIndex * sizeof(int32_t)];


            *cb_tu_coeff_bits = av1_cost_coeffs_txb(
#if CABAC_UP
                allow_update_cdf,
                ec_ctx,
#endif
                candidate_buffer_ptr,
                coeff_buffer,
                (uint16_t)cbEob,
                PLANE_TYPE_UV,
                txsize_uv,
                cb_txb_skip_context,
                cb_dc_sign_context,
                reducedTransformSetFlag);

        }
        else {
            *cb_tu_coeff_bits = av1_cost_skip_txb(
#if CABAC_UP
                allow_update_cdf,
                ec_ctx,
#endif
                candidate_buffer_ptr,
                txsize_uv,
                PLANE_TYPE_UV,
                cb_txb_skip_context);
        }
    }

    if (component_type == COMPONENT_CHROMA_CR || component_type == COMPONENT_CHROMA || component_type == COMPONENT_ALL) {

        //Estimate the rate of the transform type and coefficient for chroma Cr
        if (crEob) {

            coeff_buffer = (int32_t*)&coeff_buffer_sb->bufferCr[tuChromaOriginIndex * sizeof(int32_t)];

            *cr_tu_coeff_bits = av1_cost_coeffs_txb(
#if CABAC_UP
                allow_update_cdf,
                ec_ctx,
#endif
                candidate_buffer_ptr,
                coeff_buffer,
                (uint16_t)crEob,
                PLANE_TYPE_UV,
                txsize_uv,
                cr_txb_skip_context,
                cr_dc_sign_context,
                reducedTransformSetFlag);

        }
        else {
            *cr_tu_coeff_bits = av1_cost_skip_txb(
#if CABAC_UP
                allow_update_cdf,
                ec_ctx,
#endif
                candidate_buffer_ptr,
                txsize_uv,
                PLANE_TYPE_UV,
                cr_txb_skip_context);
        }
    }


    return return_error;
}
/*********************************************************************************
* av1_intra_full_cost function is used to estimate the cost of an intra candidate mode
* for full mode decisoion module.
*
*   @param *cu_ptr(input)
*       cu_ptr is the pointer of the target CU.
*   @param *candidate_buffer_ptr(input)
*       chromaBufferPtr is the buffer pointer of the candidate luma mode.
*   @param qp(input)
*       qp is the quantizer parameter.
*   @param luma_distortion (input)
*       luma_distortion is the intra condidate luma distortion.
*   @param lambda(input)
*       lambda is the Lagrange multiplier
**********************************************************************************/
EbErrorType Av1FullCost(
    PictureControlSet_t                    *picture_control_set_ptr,
    ModeDecisionContext_t                  *context_ptr,
    struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
    CodingUnit_t                           *cu_ptr,
    uint64_t                               *y_distortion,
    uint64_t                               *cb_distortion,
    uint64_t                               *cr_distortion,
    uint64_t                                lambda,
    uint64_t                               *y_coeff_bits,
    uint64_t                               *cb_coeff_bits,
    uint64_t                               *cr_coeff_bits,
    block_size                               bsize)
{
    UNUSED(picture_control_set_ptr);
    UNUSED(bsize);
    UNUSED(cu_ptr);
    EbErrorType return_error = EB_ErrorNone;

    // Luma and chroma rate
    uint64_t lumaRate = 0;
    uint64_t chromaRate = 0;
    uint64_t coeffRate = 0;

    // Luma and chroma SSE
    uint64_t luma_sse;
    uint64_t chromaSse;
    uint64_t totalDistortion;
    uint64_t rate;
    
    //Estimate the rate of the transform type and coefficient for Luma
    // Add fast rate to get the total rate of the subject mode
    lumaRate += candidate_buffer_ptr->candidate_ptr->fast_luma_rate;
    chromaRate += candidate_buffer_ptr->candidate_ptr->fast_chroma_rate;

    // For CFL, costs of alphas are not computed in fast loop, since they are computed in the full loop. The rate costs are added to the full loop.
    // In fast loop CFL alphas are not know yet. The chroma mode bits are calculated based on DC Mode, and if CFL is the winner compared to CFL, ChromaBits are updated in Full loop
    if (context_ptr->blk_geom->has_uv) {

        if (candidate_buffer_ptr->candidate_ptr->type == INTRA_MODE && candidate_buffer_ptr->candidate_ptr->intra_chroma_mode == UV_CFL_PRED) {

            EbBool isCflAllowed = (context_ptr->blk_geom->bwidth <= 32 &&
                context_ptr->blk_geom->bheight <= 32) ? 1 : 0;

            chromaRate += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->cflAlphaFacBits[candidate_buffer_ptr->candidate_ptr->cfl_alpha_signs][CFL_PRED_U][CFL_IDX_U(candidate_buffer_ptr->candidate_ptr->cfl_alpha_idx)] +
                candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->cflAlphaFacBits[candidate_buffer_ptr->candidate_ptr->cfl_alpha_signs][CFL_PRED_V][CFL_IDX_V(candidate_buffer_ptr->candidate_ptr->cfl_alpha_idx)];

            chromaRate += (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->intraUVmodeFacBits[isCflAllowed][candidate_buffer_ptr->candidate_ptr->intra_luma_mode][UV_CFL_PRED];
            chromaRate -= (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->intraUVmodeFacBits[isCflAllowed][candidate_buffer_ptr->candidate_ptr->intra_luma_mode][UV_DC_PRED];
        }
    }

    // Coeff rate
    coeffRate = (*y_coeff_bits + *cb_coeff_bits + *cr_coeff_bits);
    luma_sse = y_distortion[0];
    chromaSse = cb_distortion[0] + cr_distortion[0];

    // *Note - As of Oct 2011, the JCT-VC uses the PSNR forumula
    //  PSNR = (LUMA_WEIGHT * PSNRy + PSNRu + PSNRv) / (2+LUMA_WEIGHT)
    luma_sse = LUMA_WEIGHT * (luma_sse << AV1_COST_PRECISION);

    // *Note - As in JCTVC-G1102, the JCT-VC uses the Mode Decision forumula where the chromaSse has been weighted
    //  CostMode = (luma_sse + wchroma * chromaSse) + lambdaSse * rateMode
    //chromaSse = (((chromaSse * ChromaWeightFactorLd[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT); // Low delay and Random access have the same value of chroma weight

    chromaSse = (chromaSse << AV1_COST_PRECISION);

    totalDistortion = luma_sse + chromaSse;

    rate = lumaRate + chromaRate + coeffRate;
    // Assign full cost
    *(candidate_buffer_ptr->full_cost_ptr) = RDCOST(lambda, rate, totalDistortion);

    candidate_buffer_ptr->full_lambda_rate = *candidate_buffer_ptr->full_cost_ptr - totalDistortion;
    coeffRate = *y_coeff_bits;
    candidate_buffer_ptr->full_cost_luma = RDCOST(lambda, lumaRate + *y_coeff_bits, luma_sse);

    return return_error;
}

/*********************************************************************************
* merge_skip_full_cost function is used to estimate the cost of an AMVPSkip candidate
* mode for full mode decisoion module.
*
*   @param *cu_ptr(input)
*       cu_ptr is the pointer of the target CU.
*   @param *candidate_buffer_ptr(input)
*       chromaBufferPtr is the buffer pointer of the candidate luma mode.
*   @param qp(input)
*       qp is the quantizer parameter.
*   @param luma_distortion (input)
*       luma_distortion is the inter condidate luma distortion.
*   @param lambda(input)
*       lambda is the Lagrange multiplier
**********************************************************************************/
EbErrorType  Av1MergeSkipFullCost(
    PictureControlSet_t                    *picture_control_set_ptr,
    ModeDecisionContext_t                  *context_ptr,
    struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
    CodingUnit_t                           *cu_ptr,
    uint64_t                               *y_distortion,
    uint64_t                               *cb_distortion,
    uint64_t                               *cr_distortion,
    uint64_t                                lambda,
    uint64_t                               *y_coeff_bits,
    uint64_t                               *cb_coeff_bits,
    uint64_t                               *cr_coeff_bits,
    block_size                               bsize)
{
    UNUSED(bsize);
    UNUSED(context_ptr);
    UNUSED(picture_control_set_ptr);

    EbErrorType  return_error = EB_ErrorNone;
    uint64_t skipModeCtx = cu_ptr->skip_flag_context;
    uint64_t mergeRate = 0;
    uint64_t skipRate = 0;
    // Merge
    //uint64_t mergeChromaRate;
    uint64_t mergeDistortion;
    uint64_t merge_cost;
    //uint64_t mergeLumaCost;
    uint64_t mergeLumaSse;
    uint64_t mergeChromaSse;
    uint64_t coeffRate;
    //uint64_t lumaCoeffRate;

    // SKIP
    uint64_t skipDistortion;
    uint64_t skip_cost;
    //uint64_t skipLumaCost;

    // Luma and chroma transform size shift for the distortion
    uint64_t skipLumaSse;
    uint64_t skipChromaSse;

    uint64_t skipModeRate = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->skipModeFacBits[skipModeCtx][1];

    // Coeff rate
    coeffRate = (*y_coeff_bits + *cb_coeff_bits + *cr_coeff_bits);

    // Compute Merge Cost
    mergeLumaSse = y_distortion[0] << AV1_COST_PRECISION;
    mergeChromaSse = (cb_distortion[0] + cr_distortion[0]) << AV1_COST_PRECISION;

    skipLumaSse = y_distortion[1] << AV1_COST_PRECISION;
    skipChromaSse = (cb_distortion[1] + cr_distortion[1]) << AV1_COST_PRECISION;

    // *Note - As in JCTVC-G1102, the JCT-VC uses the Mode Decision forumula where the chromaSse has been weighted
    //  CostMode = (luma_sse + wchroma * chromaSse) + lambdaSse * rateMode

    //if (picture_control_set_ptr->parent_pcs_ptr->pred_structure == EB_PRED_RANDOM_ACCESS) {
    //    // Random Access
    //    if (picture_control_set_ptr->temporal_layer_index == 0) {
    //        mergeChromaSse = (((mergeChromaSse * ChromaWeightFactorRa[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else if (picture_control_set_ptr->temporal_layer_index < 3) {
    //        mergeChromaSse = (((mergeChromaSse * ChromaWeightFactorRaQpScalingL1[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else {
    //        mergeChromaSse = (((mergeChromaSse * ChromaWeightFactorRaQpScalingL3[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //}
    //else {
    //    // Low delay
    //    if (picture_control_set_ptr->temporal_layer_index == 0) {
    //        mergeChromaSse = (((mergeChromaSse * ChromaWeightFactorLd[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else {
    //        mergeChromaSse = (((mergeChromaSse * ChromaWeightFactorLdQpScaling[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //}

    // Add fast rate to get the total rate of the subject mode
    mergeRate += candidate_buffer_ptr->candidate_ptr->fast_luma_rate;
    mergeRate += candidate_buffer_ptr->candidate_ptr->fast_chroma_rate;


    mergeRate += coeffRate;

    mergeDistortion = (mergeLumaSse + mergeChromaSse);

    //merge_cost = mergeDistortion + (((lambda * coeffRate + lambda * mergeLumaRate + lambda_chroma * mergeChromaRate) + MD_OFFSET) >> MD_SHIFT);

    merge_cost = RDCOST(lambda, mergeRate, mergeDistortion);
    // mergeLumaCost = mergeLumaSse    + (((lambda * lumaCoeffRate + lambda * mergeLumaRate) + MD_OFFSET) >> MD_SHIFT);


    // *Note - As in JCTVC-G1102, the JCT-VC uses the Mode Decision forumula where the chromaSse has been weighted
    //  CostMode = (luma_sse + wchroma * chromaSse) + lambdaSse * rateMode

    //if (picture_control_set_ptr->parent_pcs_ptr->pred_structure == EB_PRED_RANDOM_ACCESS) {

    //    if (picture_control_set_ptr->temporal_layer_index == 0) {
    //        skipChromaSse = (((skipChromaSse * ChromaWeightFactorRa[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else if (picture_control_set_ptr->temporal_layer_index < 3) {
    //        skipChromaSse = (((skipChromaSse * ChromaWeightFactorRaQpScalingL1[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else {
    //        skipChromaSse = (((skipChromaSse * ChromaWeightFactorRaQpScalingL3[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //}
    //else {
    //    // Low Delay
    //    if (picture_control_set_ptr->temporal_layer_index == 0) {
    //        skipChromaSse = (((skipChromaSse * ChromaWeightFactorLd[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else {
    //        skipChromaSse = (((skipChromaSse * ChromaWeightFactorLdQpScaling[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //}

    skipDistortion = skipLumaSse + skipChromaSse;
    skipRate = skipModeRate;
    skip_cost = RDCOST(lambda, skipRate, skipDistortion);

    // Assigne full cost
    *candidate_buffer_ptr->full_cost_ptr = (skip_cost <= merge_cost) ? skip_cost : merge_cost;

    uint64_t tempDistortion;
    tempDistortion = (skip_cost <= merge_cost) ? skipDistortion : mergeDistortion;
    candidate_buffer_ptr->full_lambda_rate = *candidate_buffer_ptr->full_cost_ptr - tempDistortion;
    *candidate_buffer_ptr->full_cost_merge_ptr = merge_cost;
    *candidate_buffer_ptr->full_cost_skip_ptr = skip_cost;
    // Assigne merge flag
    candidate_buffer_ptr->candidate_ptr->merge_flag = EB_TRUE;
    // Assigne skip flag

    candidate_buffer_ptr->candidate_ptr->skip_flag = (skip_cost <= merge_cost) ? EB_TRUE : EB_FALSE;

    //CHKN:  skip_flag context is not accurate as MD does not keep skip info in sync with EncDec.



    return return_error;
}
/*********************************************************************************
* av1_intra_full_cost function is used to estimate the cost of an intra candidate mode
* for full mode decisoion module.
*
*   @param *cu_ptr(input)
*       cu_ptr is the pointer of the target CU.
*   @param *candidate_buffer_ptr(input)
*       chromaBufferPtr is the buffer pointer of the candidate luma mode.
*   @param qp(input)
*       qp is the quantizer parameter.
*   @param luma_distortion (input)
*       luma_distortion is the intra condidate luma distortion.
*   @param lambda(input)
*       lambda is the Lagrange multiplier
**********************************************************************************/
EbErrorType av1_intra_full_cost(
    PictureControlSet_t                    *picture_control_set_ptr,
    ModeDecisionContext_t                  *context_ptr,
    struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
    CodingUnit_t                           *cu_ptr,
    uint64_t                                 *y_distortion,
    uint64_t                                 *cb_distortion,
    uint64_t                                 *cr_distortion,
    uint64_t                                  lambda,
    uint64_t                                 *y_coeff_bits,
    uint64_t                                 *cb_coeff_bits,
    uint64_t                                 *cr_coeff_bits,
    block_size                              bsize)


{
    EbErrorType return_error = EB_ErrorNone;


    Av1FullCost(
        picture_control_set_ptr,
        context_ptr,
        candidate_buffer_ptr,
        cu_ptr,
        y_distortion,
        cb_distortion,
        cr_distortion,
        lambda,
        y_coeff_bits,
        cb_coeff_bits,
        cr_coeff_bits,
        bsize);



    return return_error;
}

/*********************************************************************************
* av1_inter_full_cost function is used to estimate the cost of an inter candidate mode
* for full mode decisoion module in inter frames.
*
*   @param *cu_ptr(input)
*       cu_ptr is the pointer of the target CU.
*   @param *candidate_buffer_ptr(input)
*       chromaBufferPtr is the buffer pointer of the candidate luma mode.
*   @param qp(input)
*       qp is the quantizer parameter.
*   @param luma_distortion (input)
*       luma_distortion is the inter condidate luma distortion.
*   @param lambda(input)
*       lambda is the Lagrange multiplier
**********************************************************************************/
EbErrorType av1_inter_full_cost(
    PictureControlSet_t                    *picture_control_set_ptr,
    ModeDecisionContext_t                  *context_ptr,
    struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
    CodingUnit_t                           *cu_ptr,
    uint64_t                                 *y_distortion,
    uint64_t                                 *cb_distortion,
    uint64_t                                 *cr_distortion,
    uint64_t                                  lambda,
    uint64_t                                 *y_coeff_bits,
    uint64_t                                 *cb_coeff_bits,
    uint64_t                                 *cr_coeff_bits,
    block_size                              bsize
)
{
    EbErrorType  return_error = EB_ErrorNone;

    if (candidate_buffer_ptr->candidate_ptr->merge_flag == EB_TRUE) {

        Av1MergeSkipFullCost(
            picture_control_set_ptr,
            context_ptr,
            candidate_buffer_ptr,
            cu_ptr,
            y_distortion,
            cb_distortion,
            cr_distortion,
            lambda,
            y_coeff_bits,
            cb_coeff_bits,
            cr_coeff_bits,
            bsize);
    }
    else {

        Av1FullCost(
            picture_control_set_ptr,
            context_ptr,
            candidate_buffer_ptr,
            cu_ptr,
            y_distortion,
            cb_distortion,
            cr_distortion,
            lambda,
            y_coeff_bits,
            cb_coeff_bits,
            cr_coeff_bits,
            bsize);
    }
    return return_error;
}

/************************************************************
* Coding Loop Context Generation
************************************************************/
void coding_loop_context_generation(
    ModeDecisionContext_t      *context_ptr,
    CodingUnit_t               *cu_ptr,
    uint32_t                      cu_origin_x,
    uint32_t                      cu_origin_y,
    uint32_t                      sb_sz,

    NeighborArrayUnit_t        *skip_coeff_neighbor_array,
    NeighborArrayUnit_t        *luma_dc_sign_level_coeff_neighbor_array,
    NeighborArrayUnit_t        *cb_dc_sign_level_coeff_neighbor_array,
    NeighborArrayUnit_t        *cr_dc_sign_level_coeff_neighbor_array,
    NeighborArrayUnit_t        *inter_pred_dir_neighbor_array,
    NeighborArrayUnit_t        *ref_frame_type_neighbor_array,

    NeighborArrayUnit_t        *intra_luma_mode_neighbor_array,
    NeighborArrayUnit_t        *skip_flag_neighbor_array,
    NeighborArrayUnit_t        *mode_type_neighbor_array,
    NeighborArrayUnit_t        *leaf_depth_neighbor_array,
    NeighborArrayUnit_t       *leaf_partition_neighbor_array)
{
    (void)sb_sz;
    uint32_t modeTypeLeftNeighborIndex = get_neighbor_array_unit_left_index(
        mode_type_neighbor_array,
        cu_origin_y);
    uint32_t modeTypeTopNeighborIndex = get_neighbor_array_unit_top_index(
        mode_type_neighbor_array,
        cu_origin_x);
    uint32_t leafDepthLeftNeighborIndex = get_neighbor_array_unit_left_index(
        leaf_depth_neighbor_array,
        cu_origin_y);
    uint32_t leafDepthTopNeighborIndex = get_neighbor_array_unit_top_index(
        leaf_depth_neighbor_array,
        cu_origin_x);
    uint32_t skipFlagLeftNeighborIndex = get_neighbor_array_unit_left_index(
        skip_flag_neighbor_array,
        cu_origin_y);
    uint32_t skipFlagTopNeighborIndex = get_neighbor_array_unit_top_index(
        skip_flag_neighbor_array,
        cu_origin_x);
    uint32_t intraLumaModeLeftNeighborIndex = get_neighbor_array_unit_left_index(
        intra_luma_mode_neighbor_array,
        cu_origin_y);
    uint32_t intraLumaModeTopNeighborIndex = get_neighbor_array_unit_top_index(
        intra_luma_mode_neighbor_array,
        cu_origin_x);

    uint32_t partition_left_neighbor_index = get_neighbor_array_unit_left_index(
        leaf_partition_neighbor_array,
        cu_origin_y);
    uint32_t partition_above_neighbor_index = get_neighbor_array_unit_top_index(
        leaf_partition_neighbor_array,
        cu_origin_x);

    // Intra Luma Neighbor Modes

    cu_ptr->prediction_unit_array->intra_luma_left_mode = (uint32_t)(
        (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != INTRA_MODE) ? (uint32_t)DC_PRED :
        intra_luma_mode_neighbor_array->leftArray[intraLumaModeLeftNeighborIndex]);

    cu_ptr->prediction_unit_array->intra_luma_top_mode = (uint32_t)(
        (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != INTRA_MODE) ? (uint32_t)DC_PRED :
        intra_luma_mode_neighbor_array->topArray[intraLumaModeTopNeighborIndex]);

    int32_t contextIndex;
    if (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != (uint8_t)INVALID_MODE && mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != (uint8_t)INVALID_MODE) {
        contextIndex = (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == (uint8_t)INTRA_MODE && mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] == (uint8_t)INTRA_MODE) ? 3 :
            (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == (uint8_t)INTRA_MODE || mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] == (uint8_t)INTRA_MODE) ? 1 : 0;

    }
    else  if (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != (uint8_t)INVALID_MODE) {
        contextIndex = (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == (uint8_t)INTRA_MODE) ? 2 : 0;
    }
    else if (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != (uint8_t)INVALID_MODE) {
        contextIndex = (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] == (uint8_t)INTRA_MODE) ? 2 : 0;
    }
    else {
        contextIndex = 0;
    }

    cu_ptr->is_inter_ctx = contextIndex;
    //  if(cu_ptr->is_inter_ctx!=0) //
    //      printf("ctx:%i \n",cu_ptr->is_inter_ctx);

      //   Top Intra Mode Neighbor Array instead of a Full
      // Skip Flag Context
    cu_ptr->skip_flag_context =
        (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == (uint8_t)INVALID_MODE) ? 0 :
        (skip_flag_neighbor_array->leftArray[skipFlagLeftNeighborIndex] == EB_TRUE) ? 1 : 0;
    cu_ptr->skip_flag_context +=
        (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] == (uint8_t)INVALID_MODE) ? 0 :
        (skip_flag_neighbor_array->topArray[skipFlagTopNeighborIndex] == EB_TRUE) ? 1 : 0;

    // Split Flag Context (neighbor info)
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].left_neighbor_mode = (uint32_t)(
        (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != INTRA_MODE) ? (uint32_t)DC_PRED :
        intra_luma_mode_neighbor_array->leftArray[intraLumaModeLeftNeighborIndex]);
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].left_neighbor_depth = leaf_depth_neighbor_array->leftArray[leafDepthLeftNeighborIndex];
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].top_neighbor_mode = (uint32_t)(
        (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != INTRA_MODE) ? (uint32_t)DC_PRED :
        intra_luma_mode_neighbor_array->topArray[intraLumaModeTopNeighborIndex]);
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].top_neighbor_depth = leaf_depth_neighbor_array->topArray[leafDepthTopNeighborIndex];


    // Generate Partition context
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].above_neighbor_partition = (((PartitionContext*)leaf_partition_neighbor_array->topArray)[partition_above_neighbor_index].above == (int8_t)INVALID_NEIGHBOR_DATA) ?
        0 : ((PartitionContext*)leaf_partition_neighbor_array->topArray)[partition_above_neighbor_index].above;

    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].left_neighbor_partition = (((PartitionContext*)leaf_partition_neighbor_array->leftArray)[partition_left_neighbor_index].left == (int8_t)INVALID_NEIGHBOR_DATA) ?
        0 : ((PartitionContext*)leaf_partition_neighbor_array->leftArray)[partition_left_neighbor_index].left;

    // Skip Coeff AV1 Context
    uint32_t skipCoeffLeftNeighborIndex = get_neighbor_array_unit_left_index(
        skip_coeff_neighbor_array,
        cu_origin_y);
    uint32_t skipCoeffTopNeighborIndex = get_neighbor_array_unit_top_index(
        skip_coeff_neighbor_array,
        cu_origin_x);


    cu_ptr->skip_coeff_context =
        (skip_coeff_neighbor_array->leftArray[skipCoeffLeftNeighborIndex] == (uint8_t)INVALID_NEIGHBOR_DATA) ? 0 :
        (skip_coeff_neighbor_array->leftArray[skipCoeffLeftNeighborIndex]) ? 1 : 0;


    cu_ptr->skip_coeff_context +=
        (skip_coeff_neighbor_array->topArray[skipCoeffTopNeighborIndex] == (uint8_t)INVALID_NEIGHBOR_DATA) ? 0 :
        (skip_coeff_neighbor_array->topArray[skipCoeffTopNeighborIndex]) ? 1 : 0;

    // Skip and Dc sign context generation

    block_size plane_bsize = context_ptr->blk_geom->bsize;

    cu_ptr->luma_txb_skip_context = 0;
    cu_ptr->luma_dc_sign_context = 0;
    cu_ptr->cb_txb_skip_context = 0;
    cu_ptr->cb_dc_sign_context = 0;
    cu_ptr->cr_txb_skip_context = 0;
    cu_ptr->cr_dc_sign_context = 0;

    int32_t txb_count = context_ptr->blk_geom->txb_count;
    int32_t txb_itr = 0;
    for (txb_itr = 0; txb_itr < txb_count; txb_itr++) {


        GetTxbCtx(                  //SB128_TODO move inside Full loop
            COMPONENT_LUMA,
            luma_dc_sign_level_coeff_neighbor_array,
            cu_origin_x,
            cu_origin_y,
            plane_bsize,
            context_ptr->blk_geom->txsize[txb_itr],
            &cu_ptr->luma_txb_skip_context,
            &cu_ptr->luma_dc_sign_context);


#if CHROMA_BLIND
        if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level == CHROMA_MODE_0) {
#else
        if (context_ptr->blk_geom->has_uv) {
#endif
            GetTxbCtx(
                COMPONENT_CHROMA,
                cb_dc_sign_level_coeff_neighbor_array,
                context_ptr->round_origin_x >> 1,
                context_ptr->round_origin_y >> 1,
                context_ptr->blk_geom->bsize_uv,
                context_ptr->blk_geom->txsize_uv[txb_itr],
                &cu_ptr->cb_txb_skip_context,
                &cu_ptr->cb_dc_sign_context);
            GetTxbCtx(
                COMPONENT_CHROMA,
                cr_dc_sign_level_coeff_neighbor_array,
                context_ptr->round_origin_x >> 1,
                context_ptr->round_origin_y >> 1,
                context_ptr->blk_geom->bsize_uv,
                context_ptr->blk_geom->txsize_uv[txb_itr],
                &cu_ptr->cr_txb_skip_context,
                &cu_ptr->cr_dc_sign_context);
        }
    }

    // Generate reference mode context

    cu_ptr->reference_mode_context = (uint8_t)Av1GetReferenceModeContext(
        cu_origin_x,
        cu_origin_y,
        mode_type_neighbor_array,
        inter_pred_dir_neighbor_array);

    cu_ptr->compoud_reference_type_context = (uint8_t)Av1GetCompReferenceTypeContext(
        cu_origin_x,
        cu_origin_y,
        mode_type_neighbor_array,
        inter_pred_dir_neighbor_array);

    //Collect Neighbor ref cout
    Av1CollectNeighborsRefCounts(
        cu_ptr,
        cu_origin_x,
        cu_origin_y,
        mode_type_neighbor_array,
        inter_pred_dir_neighbor_array,
        ref_frame_type_neighbor_array);

    return;
}

/********************************************
* tu_calc_cost
*   computes TU Cost and generetes TU Cbf
********************************************/
EbErrorType av1_tu_calc_cost(
    ModeDecisionCandidate_t *candidate_ptr,                        // input parameter, prediction result Ptr
    int16_t                   txb_skip_ctx,
    uint32_t                   tu_index,                             // input parameter, TU index inside the CU
    uint32_t                   y_count_non_zero_coeffs,                 // input parameter, number of non zero Y quantized coefficients
    uint32_t                   cb_count_non_zero_coeffs,                // input parameter, number of non zero cb quantized coefficients
    uint32_t                   cr_count_non_zero_coeffs,                // input parameter, number of non zero cr quantized coefficients
    uint64_t                   y_tu_distortion[DIST_CALC_TOTAL],      // input parameter, Y distortion for both Normal and Cbf zero modes
    uint64_t                   cb_tu_distortion[DIST_CALC_TOTAL],     // input parameter, Cb distortion for both Normal and Cbf zero modes
    uint64_t                   cr_tu_distortion[DIST_CALC_TOTAL],     // input parameter, Cr distortion for both Normal and Cbf zero modes
    COMPONENT_TYPE           component_type,
    uint64_t                  *y_tu_coeff_bits,                        // input parameter, Y quantized coefficients rate
    uint64_t                  *cb_tu_coeff_bits,                       // input parameter, Cb quantized coefficients rate
    uint64_t                  *cr_tu_coeff_bits,                       // input parameter, Cr quantized coefficients rate
    TxSize                  txsize,
    uint64_t                   lambda)                              // input parameter, lambda for Luma

{
    (void)cr_tu_coeff_bits;
    (void)cb_tu_coeff_bits;
    (void)cr_tu_distortion;
    (void)cb_tu_distortion;
    EbErrorType return_error = EB_ErrorNone;
    // Non Zero coeff mode variables
    uint64_t y_nonzero_coeff_distortion = y_tu_distortion[DIST_CALC_RESIDUAL];
    uint64_t y_nonzero_coeff_rate;

    uint64_t y_nonzero_coeff_cost = 0;

    // Zero Cbf mode variables
    uint64_t y_zero_coeff_distortion = y_tu_distortion[DIST_CALC_PREDICTION];

    uint64_t y_zero_coeff_luma_flag_bits_num = 0;

    uint64_t y_zero_coeff_rate;

    uint64_t y_zero_coeff_cost = 0;
    if (component_type == COMPONENT_LUMA || component_type == COMPONENT_ALL) {

        // Non Zero Distortion
        // *Note - As of Oct 2011, the JCT-VC uses the PSNR forumula
        //  PSNR = (LUMA_WEIGHT * PSNRy + PSNRu + PSNRv) / (2+LUMA_WEIGHT)
        y_nonzero_coeff_distortion = LUMA_WEIGHT * (y_nonzero_coeff_distortion << AV1_COST_PRECISION);

        // Zero distortion
        // *Note - As of Oct 2011, the JCT-VC uses the PSNR forumula
        //  PSNR = (LUMA_WEIGHT * PSNRy + PSNRu + PSNRv) / (2+LUMA_WEIGHT)
        y_zero_coeff_distortion = LUMA_WEIGHT * (y_zero_coeff_distortion << AV1_COST_PRECISION);



        // **Compute Rate

        // Esimate Cbf's Bits

        const TxSize txs_ctx = (TxSize)((txsize_sqr_map[txsize] + txsize_sqr_up_map[txsize] + 1) >> 1);
        const LV_MAP_COEFF_COST *const coeff_costs = &candidate_ptr->md_rate_estimation_ptr->coeffFacBits[txs_ctx][0];

        y_zero_coeff_luma_flag_bits_num = coeff_costs->txb_skip_cost[txb_skip_ctx][1];

        y_nonzero_coeff_rate = *y_tu_coeff_bits; // yNonZeroCbfLumaFlagBitsNum is already calculated inside y_tu_coeff_bits

        y_zero_coeff_rate = y_zero_coeff_luma_flag_bits_num;

#if CBF_ZERO_OFF
        if (1) {
#else
        if (candidate_ptr->type == INTRA_MODE) {
#endif
            y_zero_coeff_cost = 0xFFFFFFFFFFFFFFFFull;

        }
        else {

            y_zero_coeff_cost = RDCOST(lambda, y_zero_coeff_rate, y_zero_coeff_distortion);
        }

        // **Compute Cost
        y_nonzero_coeff_cost = RDCOST(lambda, y_nonzero_coeff_rate, y_nonzero_coeff_distortion);

        candidate_ptr->y_has_coeff |= (((y_count_non_zero_coeffs != 0) && (y_nonzero_coeff_cost < y_zero_coeff_cost)) << tu_index);
        *y_tu_coeff_bits = (y_nonzero_coeff_cost < y_zero_coeff_cost) ? *y_tu_coeff_bits : 0;
        y_tu_distortion[DIST_CALC_RESIDUAL] = (y_nonzero_coeff_cost < y_zero_coeff_cost) ? y_tu_distortion[DIST_CALC_RESIDUAL] : y_tu_distortion[DIST_CALC_PREDICTION];
        }
    if (component_type == COMPONENT_CHROMA_CB || component_type == COMPONENT_CHROMA || component_type == COMPONENT_ALL) {

        candidate_ptr->u_has_coeff |= ((cb_count_non_zero_coeffs != 0) << tu_index);
    }
    if (component_type == COMPONENT_CHROMA_CR || component_type == COMPONENT_CHROMA || component_type == COMPONENT_ALL) {

        candidate_ptr->v_has_coeff |= ((cr_count_non_zero_coeffs != 0) << tu_index);
    }

    return return_error;
    }

/********************************************
* tu_calc_cost
*   computes TU Cost and generetes TU Cbf
********************************************/

EbErrorType av1_tu_calc_cost_luma(

    int16_t                   txb_skip_ctx,
    ModeDecisionCandidate_t *candidate_ptr,                        // input parameter, prediction result Ptr
    uint32_t                   tu_index,                             // input parameter, TU index inside the CU
    TxSize                  tx_size,
    uint32_t                   y_count_non_zero_coeffs,                 // input parameter, number of non zero Y quantized coefficients
    uint64_t                   y_tu_distortion[DIST_CALC_TOTAL],      // input parameter, Y distortion for both Normal and Cbf zero modes
    uint64_t                  *y_tu_coeff_bits,                        // input parameter, Y quantized coefficients rate
    uint64_t                  *y_full_cost,
    uint64_t                   lambda)                              // input parameter, lambda for Luma

{

    EbErrorType return_error = EB_ErrorNone;

    // Non Zero Cbf mode variables
    uint64_t yNonZeroCbfDistortion = y_tu_distortion[DIST_CALC_RESIDUAL];

    uint64_t yNonZeroCbfRate;

    uint64_t yNonZeroCbfCost = 0;

    // Zero Cbf mode variables
    uint64_t yZeroCbfDistortion = y_tu_distortion[DIST_CALC_PREDICTION];

    uint64_t yZeroCbfLumaFlagBitsNum = 0;

    uint64_t yZeroCbfRate;

    uint64_t yZeroCbfCost = 0;

    // **Compute distortion
    // Non Zero Distortion
    // *Note - As of Oct 2011, the JCT-VC uses the PSNR forumula
    //  PSNR = (LUMA_WEIGHT * PSNRy + PSNRu + PSNRv) / (2+LUMA_WEIGHT)
    yNonZeroCbfDistortion = LUMA_WEIGHT * (yNonZeroCbfDistortion << AV1_COST_PRECISION);

    // Zero distortion
    // *Note - As of Oct 2011, the JCT-VC uses the PSNR forumula
    //  PSNR = (LUMA_WEIGHT * PSNRy + PSNRu + PSNRv) / (2+LUMA_WEIGHT)
    yZeroCbfDistortion = LUMA_WEIGHT * (yZeroCbfDistortion << AV1_COST_PRECISION);

    // **Compute Rate

    // Esimate Cbf's Bits

    const TxSize txs_ctx = (TxSize)((txsize_sqr_map[tx_size] + txsize_sqr_up_map[tx_size] + 1) >> 1);
    const LV_MAP_COEFF_COST *const coeff_costs = &candidate_ptr->md_rate_estimation_ptr->coeffFacBits[txs_ctx][0];

    yZeroCbfLumaFlagBitsNum = coeff_costs->txb_skip_cost[txb_skip_ctx][1];

    yNonZeroCbfRate = *y_tu_coeff_bits; // yNonZeroCbfLumaFlagBitsNum is already calculated inside y_tu_coeff_bits

    yZeroCbfRate = yZeroCbfLumaFlagBitsNum;

#if CBF_ZERO_OFF
    if (1) {
#else
    if (candidate_ptr->type == INTRA_MODE) {
#endif
        yZeroCbfCost = 0xFFFFFFFFFFFFFFFFull;

    }
    else {

        yZeroCbfCost = RDCOST(lambda, yZeroCbfRate, yZeroCbfDistortion);
    }

    // **Compute Cost
    yNonZeroCbfCost = RDCOST(lambda, yNonZeroCbfRate, yNonZeroCbfDistortion);
    candidate_ptr->y_has_coeff |= ((y_count_non_zero_coeffs != 0) << tu_index);
    *y_tu_coeff_bits = (yNonZeroCbfCost < yZeroCbfCost) ? *y_tu_coeff_bits : 0;
    y_tu_distortion[DIST_CALC_RESIDUAL] = (yNonZeroCbfCost < yZeroCbfCost) ? y_tu_distortion[DIST_CALC_RESIDUAL] : y_tu_distortion[DIST_CALC_PREDICTION];

    *y_full_cost = MIN(yNonZeroCbfCost, yZeroCbfCost);

    return return_error;
    }

static INLINE int32_t partition_cdf_length(block_size bsize) {
    if (bsize <= BLOCK_8X8)
        return PARTITION_TYPES;
    else if (bsize == BLOCK_128X128)
        return EXT_PARTITION_TYPES - 2;
    else
        return EXT_PARTITION_TYPES;
}

static int32_t cdf_element_prob(const int32_t *cdf,
    size_t element) {
    assert(cdf != NULL);
    return (element > 0 ? cdf[element - 1] : CDF_PROB_TOP) - cdf[element];
}
static void partition_gather_horz_alike(int32_t *out,
    block_size bsize,
    const int32_t *const in) {
    out[0] = CDF_PROB_TOP;
    out[0] -= cdf_element_prob(in, PARTITION_HORZ);
    out[0] -= cdf_element_prob(in, PARTITION_SPLIT);
    out[0] -= cdf_element_prob(in, PARTITION_HORZ_A);
    out[0] -= cdf_element_prob(in, PARTITION_HORZ_B);
    out[0] -= cdf_element_prob(in, PARTITION_VERT_A);
    if (bsize != BLOCK_128X128) out[0] -= cdf_element_prob(in, PARTITION_HORZ_4);
    out[0] = AOM_ICDF(out[0]);
    out[1] = AOM_ICDF(CDF_PROB_TOP);
}

static void partition_gather_vert_alike(int32_t *out,
    block_size bsize,
    const int32_t *const in) {
    out[0] = CDF_PROB_TOP;
    out[0] -= cdf_element_prob(in, PARTITION_VERT);
    out[0] -= cdf_element_prob(in, PARTITION_SPLIT);
    out[0] -= cdf_element_prob(in, PARTITION_HORZ_A);
    out[0] -= cdf_element_prob(in, PARTITION_VERT_A);
    out[0] -= cdf_element_prob(in, PARTITION_VERT_B);
    if (bsize != BLOCK_128X128) out[0] -= cdf_element_prob(in, PARTITION_VERT_4);
    out[0] = AOM_ICDF(out[0]);
    out[1] = AOM_ICDF(CDF_PROB_TOP);
}

//static INLINE int32_t partition_plane_context(const MacroBlockD *xd, int32_t mi_row,
//    int32_t mi_col, block_size bsize) {
//    const PARTITION_CONTEXT *above_ctx = xd->above_seg_context + mi_col;
//    const PARTITION_CONTEXT *left_ctx =
//        xd->left_seg_context + (mi_row & MAX_MIB_MASK);
//    // Minimum partition point is 8x8. Offset the bsl accordingly.
//    const int32_t bsl = mi_size_wide_log2[bsize] - mi_size_wide_log2[BLOCK_8X8];
//    int32_t above = (*above_ctx >> bsl) & 1, left = (*left_ctx >> bsl) & 1;
//
//    assert(mi_size_wide_log2[bsize] == mi_size_high_log2[bsize]);
//    assert(bsl >= 0);
//
//    return (left * 2 + above) + bsl * PARTITION_PLOFFSET;
//}


/*********************************************************************************
* split_flag_rate function is used to generate the Split rate
*
*   @param *cu_ptr(input)
*       cu_ptr is the pointer of the target CU.
*   @param split_flag(input)
*       split_flag is the split flag value.
*   @param split_rate(output)
*       split_rate contains rate.
*   @param lambda(input)
*       lambda is the Lagrange multiplier
*   @param md_rate_estimation_ptr(input)
*       md_rate_estimation_ptr is pointer to MD rate Estimation Tables
**********************************************************************************/
EbErrorType av1_split_flag_rate(
    SequenceControlSet_t                  *sequence_control_set_ptr,
    ModeDecisionContext_t                  *context_ptr,
    CodingUnit_t                           *cu_ptr,
    uint32_t                                  leaf_index,
    PartitionType                          partitionType,
    uint64_t                                 *split_rate,
    uint64_t                                  lambda,
    MdRateEstimationContext_t              *md_rate_estimation_ptr,
    uint32_t                                  tb_max_depth)
{
    (void)tb_max_depth;
    (void)leaf_index;

    const BlockGeom          *blk_geom = get_blk_geom_mds(cu_ptr->mds_idx);
    EbErrorType return_error = EB_ErrorNone;

    uint32_t cu_origin_x = context_ptr->sb_origin_x + blk_geom->origin_x;
    uint32_t cu_origin_y = context_ptr->sb_origin_y + blk_geom->origin_y;

    PartitionType p = partitionType;

    uint32_t cu_depth = blk_geom->depth;
    UNUSED(cu_depth);
    block_size bsize = blk_geom->bsize;
    ASSERT(bsize<BlockSizeS_ALL);
    const int32_t is_partition_point = blk_geom->bsize >= BLOCK_8X8;


    if (is_partition_point) {

        const int32_t hbs = (mi_size_wide[bsize] << 2) >> 1;
        const int32_t hasRows = (cu_origin_y + hbs) < sequence_control_set_ptr->luma_height;
        const int32_t hasCols = (cu_origin_x + hbs) < sequence_control_set_ptr->luma_width;

        uint32_t contextIndex = 0;

        const PARTITION_CONTEXT left_ctx = context_ptr->md_local_cu_unit[cu_ptr->mds_idx].left_neighbor_partition == (int8_t)(INVALID_NEIGHBOR_DATA) ? 0 : context_ptr->md_local_cu_unit[cu_ptr->mds_idx].left_neighbor_partition;
        const PARTITION_CONTEXT above_ctx = context_ptr->md_local_cu_unit[cu_ptr->mds_idx].above_neighbor_partition == (int8_t)(INVALID_NEIGHBOR_DATA) ? 0 : context_ptr->md_local_cu_unit[cu_ptr->mds_idx].above_neighbor_partition;

        const int32_t bsl = mi_size_wide_log2[bsize] - mi_size_wide_log2[BLOCK_8X8];

        int32_t above = (above_ctx >> bsl) & 1, left = (left_ctx >> bsl) & 1;

        assert(mi_size_wide_log2[bsize] == mi_size_high_log2[bsize]);
        assert(bsl >= 0);

        contextIndex = (left * 2 + above) + bsl * PARTITION_PLOFFSET;

        if (hasRows && hasCols) {

            *split_rate = (uint64_t)md_rate_estimation_ptr->partitionFacBits[partition_cdf_length(bsize)][partitionType];

        }
        else if (!hasRows && hasCols) {
            int32_t cdf[2];
            partition_gather_vert_alike(cdf, bsize, md_rate_estimation_ptr->partitionFacBits[contextIndex]);
            *split_rate = (uint64_t)md_rate_estimation_ptr->partitionFacBits[partition_cdf_length(bsize)][partitionType];

            *split_rate = (uint64_t)cdf[p == PARTITION_SPLIT];
        }
        else {
            int32_t cdf[2];
            partition_gather_horz_alike(cdf, bsize, md_rate_estimation_ptr->partitionFacBits[contextIndex]);
            *split_rate = (uint64_t)cdf[p == PARTITION_SPLIT];
        }
    }
    else {
        *split_rate = (uint64_t)md_rate_estimation_ptr->partitionFacBits[0][partitionType];
    }

    *split_rate = RDCOST(lambda, *split_rate, 0);

    return return_error;
}

/********************************************
* tu_calc_cost
*   Computes TU Cost and generetes TU Cbf
*   at the level of the encode pass
********************************************/
EbErrorType av1_encode_tu_calc_cost(
    EncDecContext_t          *context_ptr,
    uint32_t                   *count_non_zero_coeffs,
    uint64_t                    y_tu_distortion[DIST_CALC_TOTAL],
    uint64_t                   *y_tu_coeff_bits,
    uint32_t                    component_mask
)
{
    CodingUnit_t              *cu_ptr = context_ptr->cu_ptr;
    uint32_t                     tu_index = context_ptr->txb_itr;
    MdRateEstimationContext_t *md_rate_estimation_ptr = context_ptr->md_rate_estimation_ptr;
    uint64_t                     lambda = context_ptr->full_lambda;
    uint32_t                     y_count_non_zero_coeffs = count_non_zero_coeffs[0];
    uint32_t                     cb_count_non_zero_coeffs = count_non_zero_coeffs[1];
    uint32_t                     cr_count_non_zero_coeffs = count_non_zero_coeffs[2];

    EbErrorType return_error = EB_ErrorNone;

    // Non Zero Cbf mode variables
    uint64_t yNonZeroCbfDistortion = y_tu_distortion[DIST_CALC_RESIDUAL];

    uint64_t yNonZeroCbfRate;

    uint64_t yNonZeroCbfCost = 0;

    // Zero Cbf mode variables
    uint64_t yZeroCbfDistortion = y_tu_distortion[DIST_CALC_PREDICTION];

    uint64_t yZeroCbfLumaFlagBitsNum = 0;

    uint64_t yZeroCbfRate;

    uint64_t yZeroCbfCost = 0;
    int16_t  txb_skip_ctx = cu_ptr->luma_txb_skip_context;

    // **Compute distortion
    if (component_mask == PICTURE_BUFFER_DESC_LUMA_MASK || component_mask == PICTURE_BUFFER_DESC_FULL_MASK) {
        // Non Zero Distortion
        // *Note - As of Oct 2011, the JCT-VC uses the PSNR forumula
        //  PSNR = (LUMA_WEIGHT * PSNRy + PSNRu + PSNRv) / (2+LUMA_WEIGHT)
        yNonZeroCbfDistortion = LUMA_WEIGHT * (yNonZeroCbfDistortion << AV1_COST_PRECISION);


        // Zero distortion
        // *Note - As of Oct 2011, the JCT-VC uses the PSNR forumula
        //  PSNR = (LUMA_WEIGHT * PSNRy + PSNRu + PSNRv) / (2+LUMA_WEIGHT)
        yZeroCbfDistortion = LUMA_WEIGHT * (yZeroCbfDistortion << AV1_COST_PRECISION);
        TxSize    txSize = context_ptr->blk_geom->txsize[context_ptr->txb_itr];

        const TxSize txs_ctx = (TxSize)((txsize_sqr_map[txSize] + txsize_sqr_up_map[txSize] + 1) >> 1);
        const LV_MAP_COEFF_COST *const coeff_costs = &md_rate_estimation_ptr->coeffFacBits[txs_ctx][0];

        yZeroCbfLumaFlagBitsNum = coeff_costs->txb_skip_cost[txb_skip_ctx][1];

        yNonZeroCbfRate = *y_tu_coeff_bits; // yNonZeroCbfLumaFlagBitsNum is already calculated inside y_tu_coeff_bits

        yZeroCbfRate = yZeroCbfLumaFlagBitsNum;
#if ENABLE_EOB_ZERO_CHECK
        TransformUnit_t       *txb_ptr = &cu_ptr->transform_unit_array[context_ptr->txb_itr];
        if (txb_ptr->transform_type[PLANE_TYPE_Y] != DCT_DCT) {
#else
#if CBF_ZERO_OFF || TX_TYPE_FIX
        if (1) {
#else
        if (cu_ptr->prediction_mode_flag == INTRA_MODE) {
#endif
#endif
            yZeroCbfCost = 0xFFFFFFFFFFFFFFFFull;

        }
        else {

            yZeroCbfCost = RDCOST(lambda, yZeroCbfRate, yZeroCbfDistortion);
        }

        // **Compute Cost
        yNonZeroCbfCost = RDCOST(lambda, yNonZeroCbfRate, yNonZeroCbfDistortion);
        cu_ptr->transform_unit_array[tu_index].y_has_coeff = ((y_count_non_zero_coeffs != 0) && (yNonZeroCbfCost < yZeroCbfCost)) ? EB_TRUE : EB_FALSE;
        *y_tu_coeff_bits = (yNonZeroCbfCost < yZeroCbfCost) ? *y_tu_coeff_bits : 0;
        y_tu_distortion[DIST_CALC_RESIDUAL] = (yNonZeroCbfCost < yZeroCbfCost) ? y_tu_distortion[DIST_CALC_RESIDUAL] : y_tu_distortion[DIST_CALC_PREDICTION];

        }
    else {
        cu_ptr->transform_unit_array[tu_index].y_has_coeff = EB_FALSE;
    }
    cu_ptr->transform_unit_array[tu_index].u_has_coeff = cb_count_non_zero_coeffs != 0 ? EB_TRUE : EB_FALSE;
    cu_ptr->transform_unit_array[tu_index].v_has_coeff = cr_count_non_zero_coeffs != 0 ? EB_TRUE : EB_FALSE;

    return return_error;
    }


uint64_t GetPMCost(
    uint64_t                   lambda,
    uint64_t                   tuDistortion,
    uint64_t                   y_tu_coeff_bits
)
{

    uint64_t yNonZeroCbfDistortion = LUMA_WEIGHT * (tuDistortion << COST_PRECISION);
    uint64_t yNonZeroCbfRate = (y_tu_coeff_bits);
    uint64_t yNonZeroCbfCost = yNonZeroCbfDistortion + (((lambda       * yNonZeroCbfRate) + MD_OFFSET) >> MD_SHIFT);

    return yNonZeroCbfCost;
}
