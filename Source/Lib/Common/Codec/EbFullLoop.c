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

#include "EbDefinitions.h"
#include "EbModeDecisionProcess.h"
#include "EbTransforms.h"
#include "EbFullLoop.h"
#include "EbRateDistortionCost.h"
#include "aom_dsp_rtcd.h"
#ifdef __GNUC__
#define LIKELY(v) __builtin_expect(v, 1)
#define UNLIKELY(v) __builtin_expect(v, 0)
#else
#define LIKELY(v) (v)
#define UNLIKELY(v) (v)
#endif
static PartitionType from_shape_to_part[] =
{
PARTITION_NONE,
PARTITION_HORZ,
PARTITION_VERT,
PARTITION_HORZ_A,
PARTITION_HORZ_B,
PARTITION_VERT_A,
PARTITION_VERT_B,
PARTITION_HORZ_4,
PARTITION_VERT_4,
PARTITION_SPLIT

};
void quantize_b_helper_c_II(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
    int32_t skip_block, const int16_t *zbin_ptr,
    const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr,
    uint16_t *eob_ptr, const int16_t *scan,
    const int16_t *iscan, const qm_val_t *qm_ptr,
    const qm_val_t *iqm_ptr, const int32_t log_scale) {
    const int32_t zbins[2] = { ROUND_POWER_OF_TWO(zbin_ptr[0], log_scale),
                           ROUND_POWER_OF_TWO(zbin_ptr[1], log_scale) };
    const int32_t nzbins[2] = { zbins[0] * -1, zbins[1] * -1 };
    int32_t i, non_zero_count = (int32_t)n_coeffs, eob = -1;
    (void)iscan;

    memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
    memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

    if (!skip_block) {
        // Pre-scan pass
        for (i = (int32_t)n_coeffs - 1; i >= 0; i--) {
            const int32_t rc = scan[i];
            const qm_val_t wt = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
            const int32_t coeff = coeff_ptr[rc] * wt;

            if (coeff < (zbins[rc != 0] * (1 << AOM_QM_BITS)) &&
                coeff >(nzbins[rc != 0] * (1 << AOM_QM_BITS)))
                non_zero_count--;
            else
                break;
        }

        // Quantization pass: All coefficients with index >= zero_flag are
        // skippable. Note: zero_flag can be zero.
        for (i = 0; i < non_zero_count; i++) {
            const int32_t rc = scan[i];
            const int32_t coeff = coeff_ptr[rc];
            const int32_t coeff_sign = (coeff >> 31);
            const int32_t abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
            int32_t tmp32;

            const qm_val_t wt = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
            if (abs_coeff * wt >= (zbins[rc != 0] << AOM_QM_BITS)) {
                int64_t tmp =
                    clamp(abs_coeff + ROUND_POWER_OF_TWO(round_ptr[rc != 0], log_scale),
                        INT16_MIN, INT16_MAX);
                tmp *= wt;
                tmp32 = (int32_t)(((((tmp * quant_ptr[rc != 0]) >> 16) + tmp) *
                    quant_shift_ptr[rc != 0]) >>
                    (16 - log_scale + AOM_QM_BITS));  // quantization
                qcoeff_ptr[rc] = (tmp32 ^ coeff_sign) - coeff_sign;
                const int32_t iwt = iqm_ptr != NULL ? iqm_ptr[rc] : (1 << AOM_QM_BITS);
                const int32_t dequant =
                    (dequant_ptr[rc != 0] * iwt + (1 << (AOM_QM_BITS - 1))) >>
                    AOM_QM_BITS;
                const tran_low_t abs_dqcoeff = (tmp32 * dequant) >> log_scale;
                dqcoeff_ptr[rc] = (tran_low_t)((abs_dqcoeff ^ coeff_sign) - coeff_sign);

                if (tmp32) eob = i;
            }
        }
    }
    *eob_ptr = (uint16_t)(eob + 1);
}
void aom_quantize_b_c_II(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
    int32_t skip_block, const int16_t *zbin_ptr,
    const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr,
    uint16_t *eob_ptr, const int16_t *scan,
    const int16_t *iscan) {
    quantize_b_helper_c_II(coeff_ptr, n_coeffs, skip_block, zbin_ptr, round_ptr,
        quant_ptr, quant_shift_ptr, qcoeff_ptr, dqcoeff_ptr,
        dequant_ptr, eob_ptr, scan, iscan, NULL, NULL, 0);
}

void aom_quantize_b_32x32_c_II(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
    int32_t skip_block, const int16_t *zbin_ptr,
    const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr,
    tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
    const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan) {
    quantize_b_helper_c_II(coeff_ptr, n_coeffs, skip_block, zbin_ptr, round_ptr,
        quant_ptr, quant_shift_ptr, qcoeff_ptr, dqcoeff_ptr,
        dequant_ptr, eob_ptr, scan, iscan, NULL, NULL, 1);
}

void aom_quantize_b_64x64_c_II(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
    int32_t skip_block, const int16_t *zbin_ptr,
    const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr,
    tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
    const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan) {
    quantize_b_helper_c_II(coeff_ptr, n_coeffs, skip_block, zbin_ptr, round_ptr,
        quant_ptr, quant_shift_ptr, qcoeff_ptr, dqcoeff_ptr,
        dequant_ptr, eob_ptr, scan, iscan, NULL, NULL, 2);
}

void quantize_b_helper_c(
    const tran_low_t *coeff_ptr,
    int32_t stride,
#
    int32_t width,
    int32_t height,

    intptr_t n_coeffs,
    int32_t skip_block,
    const int16_t *zbin_ptr,
    const int16_t *round_ptr,
    const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr,
    tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr,
    const int16_t *dequant_ptr,
    uint16_t *eob_ptr,
    const int16_t *scan,
    const int16_t *iscan,
    const qm_val_t *qm_ptr,
    const qm_val_t *iqm_ptr,
    const int32_t log_scale)
{

    const int32_t zbins[2] = {
        ROUND_POWER_OF_TWO(zbin_ptr[0], log_scale),
        ROUND_POWER_OF_TWO(zbin_ptr[1], log_scale)
    };
    const int32_t nzbins[2] = { zbins[0] * -1, zbins[1] * -1 };
    int32_t i, non_zero_count = (int32_t)n_coeffs, eob = -1;
    (void)iscan;

    // Nader quantisation
    for (int32_t x = 0; x < height; x++) {
        memset(qcoeff_ptr + (x * stride), 0, width /*n_coeffs*/ * sizeof(*qcoeff_ptr));
        memset(dqcoeff_ptr + (x * stride), 0, width /*n_coeffs*/ * sizeof(*dqcoeff_ptr));
    }


    if (!skip_block) {
        // Pre-scan pass
        for (i = (int32_t)n_coeffs - 1; i >= 0; i--) {
            const int32_t mapRc = scan[i];

            const int32_t rc = ((mapRc / MIN(32, height))  * stride) + (mapRc % MIN(32, width));

            const qm_val_t wt = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
            const int32_t coeff = coeff_ptr[rc] * wt;


            ////if (mapRc != NewTab[rc])
            //printf("%d\n", coeff);

            if (coeff < (zbins[rc != 0] * (1 << AOM_QM_BITS)) &&
                coeff >(nzbins[rc != 0] * (1 << AOM_QM_BITS)))
                non_zero_count--;
            else
                break;
        }
        // Quantization pass: All coefficients with index >= zero_flag are
        // skippable. Note: zero_flag can be zero.
        for (i = 0; i < non_zero_count; i++) {
            const int32_t mapRc = scan[i];

            const int32_t rc = ((mapRc / MIN(32, height))  * stride) + (mapRc % MIN(32, width));
            const int32_t coeff = coeff_ptr[rc];
            const int32_t coeff_sign = (coeff >> 31);
            const int32_t abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
            int32_t tmp32;

            const qm_val_t wt = qm_ptr != NULL ? qm_ptr[mapRc] : (1 << AOM_QM_BITS);

            if (abs_coeff * wt >= (zbins[rc != 0] << AOM_QM_BITS)) {

                int64_t tmp = clamp(abs_coeff + ROUND_POWER_OF_TWO(round_ptr[rc != 0], log_scale), INT16_MIN, INT16_MAX);

                tmp *= wt;

                tmp32 = (int32_t)(((((tmp * quant_ptr[rc != 0]) >> 16) + tmp) *    quant_shift_ptr[rc != 0]) >> (16 - log_scale + AOM_QM_BITS));  // quantization

                qcoeff_ptr[rc] = (tmp32 ^ coeff_sign) - coeff_sign;

                const int32_t iwt = iqm_ptr != NULL ? iqm_ptr[mapRc] : (1 << AOM_QM_BITS);

                const int32_t dequant = (dequant_ptr[rc != 0] * iwt + (1 << (AOM_QM_BITS - 1))) >> AOM_QM_BITS;

                dqcoeff_ptr[rc] = qcoeff_ptr[rc] * dequant / (1 << log_scale);

                if (tmp32) eob = i;
            }
        }
    }


    *eob_ptr = (uint16_t)(eob + 1);
}
void highbd_quantize_b_helper_c(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block,
    const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan, const qm_val_t *qm_ptr,
    const qm_val_t *iqm_ptr, const int32_t log_scale) {
    int32_t i, eob = -1;
    const int32_t zbins[2] = { ROUND_POWER_OF_TWO(zbin_ptr[0], log_scale),
        ROUND_POWER_OF_TWO(zbin_ptr[1], log_scale) };
    const int32_t nzbins[2] = { zbins[0] * -1, zbins[1] * -1 };
    int32_t dequant;
    int32_t idx_arr[4096];
    (void)iscan;
    int32_t idx = 0;

    memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
    memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

    if (!skip_block) {
        // Pre-scan pass
        for (i = 0; i < n_coeffs; i++) {
            const int32_t rc = scan[i];
            const qm_val_t wt = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
            const int32_t coeff = coeff_ptr[rc] * wt;

            // If the coefficient is out of the base ZBIN range, keep it for
            // quantization.
            if (coeff >= (zbins[rc != 0] * (1 << AOM_QM_BITS)) ||
                coeff <= (nzbins[rc != 0] * (1 << AOM_QM_BITS)))
                idx_arr[idx++] = i;
        }

        // Quantization pass: only process the coefficients selected in
        // pre-scan pass. Note: idx can be zero.
        for (i = 0; i < idx; i++) {
            const int32_t rc = scan[idx_arr[i]];
            const int32_t coeff = coeff_ptr[rc];
            const int32_t coeff_sign = (coeff >> 31);
            const qm_val_t wt = qm_ptr != NULL ? qm_ptr[rc] : (1 << AOM_QM_BITS);
            const qm_val_t iwt = iqm_ptr != NULL ? iqm_ptr[rc] : (1 << AOM_QM_BITS);
            const int32_t abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
            const int64_t tmp1 =
                abs_coeff + ROUND_POWER_OF_TWO(round_ptr[rc != 0], log_scale);
            const int64_t tmpw = tmp1 * wt;
            const int64_t tmp2 = ((tmpw * quant_ptr[rc != 0]) >> 16) + tmpw;
            const int32_t abs_qcoeff = (int32_t)((tmp2 * quant_shift_ptr[rc != 0]) >>
                (16 - log_scale + AOM_QM_BITS));
            qcoeff_ptr[rc] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
            dequant = (dequant_ptr[rc != 0] * iwt + (1 << (AOM_QM_BITS - 1))) >>
                AOM_QM_BITS;
            const tran_low_t abs_dqcoeff = (abs_qcoeff * dequant) >> log_scale;
            dqcoeff_ptr[rc] = (tran_low_t)((abs_dqcoeff ^ coeff_sign) - coeff_sign);
            if (abs_qcoeff) eob = idx_arr[i];
        }
    }
    *eob_ptr = (uint16_t)(eob + 1);
}

void aom_highbd_quantize_b_c(const tran_low_t *coeff_ptr, intptr_t n_coeffs,
    int32_t skip_block, const int16_t *zbin_ptr,
    const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr,
    tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr,
    const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan) {
    highbd_quantize_b_helper_c(coeff_ptr, n_coeffs, skip_block, zbin_ptr,
        round_ptr, quant_ptr, quant_shift_ptr, qcoeff_ptr,
        dqcoeff_ptr, dequant_ptr, eob_ptr, scan, iscan,
        NULL, NULL, 0);
}

void aom_highbd_quantize_b_32x32_c(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block,
    const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan) {
    highbd_quantize_b_helper_c(coeff_ptr, n_coeffs, skip_block, zbin_ptr,
        round_ptr, quant_ptr, quant_shift_ptr, qcoeff_ptr,
        dqcoeff_ptr, dequant_ptr, eob_ptr, scan, iscan,
        NULL, NULL, 1);
}

void aom_highbd_quantize_b_64x64_c(
    const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block,
    const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr,
    const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr,
    const int16_t *scan, const int16_t *iscan) {
    highbd_quantize_b_helper_c(coeff_ptr, n_coeffs, skip_block, zbin_ptr,
        round_ptr, quant_ptr, quant_shift_ptr, qcoeff_ptr,
        dqcoeff_ptr, dequant_ptr, eob_ptr, scan, iscan,
        NULL, NULL, 2);
}

void av1_highbd_quantize_b_facade(const tran_low_t *coeff_ptr,
    intptr_t n_coeffs, const MacroblockPlane *p,
    tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, uint16_t *eob_ptr,
    const SCAN_ORDER *sc,
    const QUANT_PARAM *qparam) {
    // obsolete skip_block
    const int32_t skip_block = 0;
    const qm_val_t *qm_ptr = qparam->qmatrix;
    const qm_val_t *iqm_ptr = qparam->iqmatrix;
    if (qm_ptr != NULL && iqm_ptr != NULL) {
        highbd_quantize_b_helper_c(coeff_ptr, n_coeffs, skip_block, p->zbin_QTX,
            p->round_QTX, p->quant_QTX, p->quant_shift_QTX,
            qcoeff_ptr, dqcoeff_ptr, p->dequant_QTX, eob_ptr,
            sc->scan, sc->iscan, qm_ptr, iqm_ptr,
            qparam->log_scale);
    }
    else {
        switch (qparam->log_scale) {
        case 0:
            if (LIKELY(n_coeffs >= 8)) {

                aom_highbd_quantize_b(coeff_ptr, n_coeffs, skip_block, p->zbin_QTX,
                    p->round_QTX, p->quant_QTX, p->quant_shift_QTX,
                    qcoeff_ptr, dqcoeff_ptr, p->dequant_QTX,
                    eob_ptr, sc->scan, sc->iscan);

            }
            else {

                aom_highbd_quantize_b_c(coeff_ptr, n_coeffs, skip_block, p->zbin_QTX,
                    p->round_QTX, p->quant_QTX,
                    p->quant_shift_QTX, qcoeff_ptr, dqcoeff_ptr,
                    p->dequant_QTX, eob_ptr, sc->scan, sc->iscan);
            }
            break;
        case 1:
            aom_highbd_quantize_b_32x32(
                coeff_ptr, n_coeffs, skip_block, p->zbin_QTX, p->round_QTX,
                p->quant_QTX, p->quant_shift_QTX, qcoeff_ptr, dqcoeff_ptr,
                p->dequant_QTX, eob_ptr, sc->scan, sc->iscan);
            break;
        case 2:
            aom_highbd_quantize_b_64x64(
                coeff_ptr, n_coeffs, skip_block, p->zbin_QTX, p->round_QTX,
                p->quant_QTX, p->quant_shift_QTX, qcoeff_ptr, dqcoeff_ptr,
                p->dequant_QTX, eob_ptr, sc->scan, sc->iscan);
            break;
        default: assert(0);
        }
    }
}

/*
static INLINE void highbd_quantize_dc(
    const tran_low_t *coeff_ptr, int32_t n_coeffs, int32_t skip_block,
    const int16_t *round_ptr, const int16_t quant, tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr, const int16_t dequant_ptr, uint16_t *eob_ptr,
    const qm_val_t *qm_ptr, const qm_val_t *iqm_ptr, const int32_t log_scale) {
    int32_t eob = -1;

    memset(qcoeff_ptr, 0, n_coeffs * sizeof(*qcoeff_ptr));
    memset(dqcoeff_ptr, 0, n_coeffs * sizeof(*dqcoeff_ptr));

    if (!skip_block) {
        const qm_val_t wt = qm_ptr != NULL ? qm_ptr[0] : (1 << AOM_QM_BITS);
        const qm_val_t iwt = iqm_ptr != NULL ? iqm_ptr[0] : (1 << AOM_QM_BITS);
        const int32_t coeff = coeff_ptr[0];
        const int32_t coeff_sign = (coeff >> 31);
        const int32_t abs_coeff = (coeff ^ coeff_sign) - coeff_sign;
        const int64_t tmp = abs_coeff + ROUND_POWER_OF_TWO(round_ptr[0], log_scale);
        const int64_t tmpw = tmp * wt;
        const int32_t abs_qcoeff =
            (int32_t)((tmpw * quant) >> (16 - log_scale + AOM_QM_BITS));
        qcoeff_ptr[0] = (tran_low_t)((abs_qcoeff ^ coeff_sign) - coeff_sign);
        const int32_t dequant =
            (dequant_ptr * iwt + (1 << (AOM_QM_BITS - 1))) >> AOM_QM_BITS;

        const tran_low_t abs_dqcoeff = (abs_qcoeff * dequant) >> log_scale;
        dqcoeff_ptr[0] = (tran_low_t)((abs_dqcoeff ^ coeff_sign) - coeff_sign);
        if (abs_qcoeff) eob = 0;
    }
    *eob_ptr = (uint16_t)(eob + 1);
}
*/

void av1_quantize_b_facade_II(
    const tran_low_t *coeff_ptr,
    int32_t stride,
    int32_t                width,
    int32_t                height,
    intptr_t n_coeffs,

    const MacroblockPlane *p,
    tran_low_t *qcoeff_ptr,
    tran_low_t *dqcoeff_ptr,
    uint16_t *eob_ptr,
    const SCAN_ORDER *sc,
    const QUANT_PARAM *qparam)
{
    // obsolete skip_block
    const int32_t skip_block = 0;
    const qm_val_t *qm_ptr = qparam->qmatrix;
    const qm_val_t *iqm_ptr = qparam->iqmatrix;
    if (qm_ptr != NULL && iqm_ptr != NULL) {
        quantize_b_helper_c(
            coeff_ptr,
            stride,
            width,
            height,
            n_coeffs,
            skip_block,
            p->zbin_QTX,
            p->round_QTX,
            p->quant_QTX,
            p->quant_shift_QTX,
            qcoeff_ptr,
            dqcoeff_ptr,
            p->dequant_QTX,
            eob_ptr,
            sc->scan,
            sc->iscan,
            qm_ptr,
            iqm_ptr,
            qparam->log_scale);
    }
    else {
        switch (qparam->log_scale) {
        case 0:
            aom_quantize_b(coeff_ptr, n_coeffs, skip_block, p->zbin_QTX,
                p->round_QTX, p->quant_QTX, p->quant_shift_QTX,
                qcoeff_ptr, dqcoeff_ptr, p->dequant_QTX, eob_ptr,
                sc->scan, sc->iscan);

            break;
        case 1:

            aom_quantize_b_32x32(coeff_ptr, n_coeffs, skip_block, p->zbin_QTX,
                p->round_QTX, p->quant_QTX, p->quant_shift_QTX,
                qcoeff_ptr, dqcoeff_ptr, p->dequant_QTX, eob_ptr,
                sc->scan, sc->iscan);

            break;
        case 2:

            aom_quantize_b_64x64(coeff_ptr, n_coeffs, skip_block, p->zbin_QTX,
                p->round_QTX, p->quant_QTX, p->quant_shift_QTX,
                qcoeff_ptr, dqcoeff_ptr, p->dequant_QTX, eob_ptr,
                sc->scan, sc->iscan);

            break;
        default: assert(0);
        }
    }
}

/*********************************************************************
* UnifiedQuantizeInvQuantize
*
*  Unified Quant +iQuant
*********************************************************************/
void av1_quantize_inv_quantize_ii(
    PictureControlSet_t  *picture_control_set_ptr,
    int32_t               *coeff,
    const uint32_t          coeff_stride,
    int32_t               *quant_coeff,
    int32_t               *recon_coeff,
    uint32_t                qp,
    uint32_t                width,
    uint32_t                height,
    TxSize              transform_size,
    uint16_t                *eob,
    //  MacroblockPlane      candidate_plane,
    EbAsm                asm_type,
    uint32_t                *y_count_non_zero_coeffs,
#if !PF_N2_SUPPORT
    EbPfMode              pf_mode,
#endif
    uint8_t                 enable_contouring_qc_update_flag,
    uint32_t                component_type,
    uint32_t                bit_increment,

    TxType               tx_type,
    EbBool               clean_sparse_coeff_flag)
{
    (void)clean_sparse_coeff_flag;
    (void)enable_contouring_qc_update_flag;
#if !PF_N2_SUPPORT
    (void)pf_mode;
#endif
    (void)asm_type;
#if !ADD_DELTA_QP_SUPPORT
    (void) qp;
#endif
    MacroblockPlane      candidate_plane ;

    //    EB_SLICE          slice_type = picture_control_set_ptr->slice_type;
    //    uint32_t            temporal_layer_index = picture_control_set_ptr->temporal_layer_index;

    // Use a flat matrix (i.e. no weighting) for 1D and Identity transforms
    //const qm_val_t *qmatrix =
    //    IS_2D_TRANSFORM(tx_type) ? pd->seg_qmatrix[seg_id][qm_tx_size]
    //    : cm->gqmatrix[NUM_QM_LEVELS - 1][0][qm_tx_size];
    //const qm_val_t *iqmatrix =
    //    IS_2D_TRANSFORM(tx_type)
    //    ? pd->seg_iqmatrix[seg_id][qm_tx_size]
    //    : cm->giqmatrix[NUM_QM_LEVELS - 1][0][qm_tx_size];

    const qm_val_t *qMatrix = picture_control_set_ptr->parent_pcs_ptr->gqmatrix[NUM_QM_LEVELS - 1][0][transform_size];
    const qm_val_t *iqMatrix = picture_control_set_ptr->parent_pcs_ptr->giqmatrix[NUM_QM_LEVELS - 1][0][transform_size];
#if ADD_DELTA_QP_SUPPORT
    uint32_t qIndex = qp;
#else
    uint32_t qIndex = picture_control_set_ptr->parent_pcs_ptr->base_qindex;
#endif
    if (bit_increment == 0) {
        if (component_type == COMPONENT_LUMA) {
            candidate_plane.quant_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.y_quant[qIndex];
            candidate_plane.quant_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.y_quant_fp[qIndex];
            candidate_plane.round_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.y_round_fp[qIndex];
            candidate_plane.quant_shift_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.y_quant_shift[qIndex];
            candidate_plane.zbin_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.y_zbin[qIndex];
            candidate_plane.round_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.y_round[qIndex];
            candidate_plane.dequant_QTX = picture_control_set_ptr->parent_pcs_ptr->deqMd.y_dequant_QTX[qIndex];
        }

        if (component_type == COMPONENT_CHROMA_CB) {
            candidate_plane.quant_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.u_quant[qIndex];
            candidate_plane.quant_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.u_quant_fp[qIndex];
            candidate_plane.round_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.u_round_fp[qIndex];
            candidate_plane.quant_shift_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.u_quant_shift[qIndex];
            candidate_plane.zbin_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.u_zbin[qIndex];
            candidate_plane.round_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.u_round[qIndex];
            candidate_plane.dequant_QTX = picture_control_set_ptr->parent_pcs_ptr->deqMd.u_dequant_QTX[qIndex];

        }

        if (component_type == COMPONENT_CHROMA_CR) {
            candidate_plane.quant_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.v_quant[qIndex];
            candidate_plane.quant_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.v_quant_fp[qIndex];
            candidate_plane.round_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.v_round_fp[qIndex];
            candidate_plane.quant_shift_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.v_quant_shift[qIndex];
            candidate_plane.zbin_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.v_zbin[qIndex];
            candidate_plane.round_QTX = picture_control_set_ptr->parent_pcs_ptr->quantsMd.v_round[qIndex];
            candidate_plane.dequant_QTX = picture_control_set_ptr->parent_pcs_ptr->deqMd.v_dequant_QTX[qIndex];

        }

    }
    else {
        if (component_type == COMPONENT_LUMA) {
            candidate_plane.quant_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_quant[qIndex];
            candidate_plane.quant_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_quant_fp[qIndex];
            candidate_plane.round_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_round_fp[qIndex];
            candidate_plane.quant_shift_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_quant_shift[qIndex];
            candidate_plane.zbin_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_zbin[qIndex];
            candidate_plane.round_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.y_round[qIndex];
            candidate_plane.dequant_QTX = picture_control_set_ptr->parent_pcs_ptr->deq.y_dequant_QTX[qIndex];
        }

        if (component_type == COMPONENT_CHROMA_CB) {
            candidate_plane.quant_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_quant[qIndex];
            candidate_plane.quant_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_quant_fp[qIndex];
            candidate_plane.round_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_round_fp[qIndex];
            candidate_plane.quant_shift_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_quant_shift[qIndex];
            candidate_plane.zbin_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_zbin[qIndex];
            candidate_plane.round_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.u_round[qIndex];
            candidate_plane.dequant_QTX = picture_control_set_ptr->parent_pcs_ptr->deq.u_dequant_QTX[qIndex];

        }

        if (component_type == COMPONENT_CHROMA_CR) {
            candidate_plane.quant_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_quant[qIndex];
            candidate_plane.quant_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_quant_fp[qIndex];
            candidate_plane.round_fp_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_round_fp[qIndex];
            candidate_plane.quant_shift_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_quant_shift[qIndex];
            candidate_plane.zbin_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_zbin[qIndex];
            candidate_plane.round_QTX = picture_control_set_ptr->parent_pcs_ptr->quants.v_round[qIndex];
            candidate_plane.dequant_QTX = picture_control_set_ptr->parent_pcs_ptr->deq.v_dequant_QTX[qIndex];
        }
    }

    const SCAN_ORDER *const scan_order = &av1_scan_orders[transform_size][tx_type];  //get_scan(tx_size, tx_type);

    const int32_t n_coeffs = av1_get_max_eob(transform_size);

    QUANT_PARAM qparam;

    qparam.log_scale = av1_get_tx_scale(transform_size);
    qparam.tx_size = transform_size;
    qparam.qmatrix = qMatrix;
    qparam.iqmatrix = iqMatrix;

    if (bit_increment)
        av1_highbd_quantize_b_facade(
        (tran_low_t*)coeff,
            n_coeffs,
            &candidate_plane,
            quant_coeff,
            (tran_low_t*)recon_coeff,
            eob,
            scan_order,
            &qparam);
    else
        av1_quantize_b_facade_II(
        (tran_low_t*)coeff,
            coeff_stride,
            width,
            height,
            n_coeffs,
            &candidate_plane,
            quant_coeff,
            (tran_low_t*)recon_coeff,
            eob,
            scan_order,
            &qparam);

    *y_count_non_zero_coeffs = *eob;

    }
void av1_quantize_inv_quantize(
    PictureControlSet_t  *picture_control_set_ptr,
    int32_t               *coeff,
    const uint32_t          coeff_stride,
    int32_t               *quant_coeff,
    int32_t               *recon_coeff,
    uint32_t                qp,
    uint32_t                width,
    uint32_t                height,
    TxSize               txsize,
    uint16_t                *eob,
    MacroblockPlane      candidate_plane,
    EbAsm                asm_type,
    uint32_t                *y_count_non_zero_coeffs,
#if !PF_N2_SUPPORT
    EbPfMode              pf_mode,
#endif
    uint8_t                 enable_contouring_qc_update_flag,
    uint32_t                component_type,
    uint32_t                bit_increment,
    TxType               tx_type,
    EbBool               clean_sparse_coeff_flag)
{
    (void)coeff_stride;
    (void)candidate_plane;
    (void)enable_contouring_qc_update_flag;
#if !PF_N2_SUPPORT
    (void)pf_mode;
#endif
    //Note: Transformed, Quantized, iQuantized coeff are stored in 1D fashion. 64x64 is hence in the first 32x32 corner.

    uint32_t i;

    for (i = 0; i < height; i++)
    {
        memset(quant_coeff + i * width, 0, width * sizeof(int32_t));
        memset(recon_coeff + i * width, 0, width * sizeof(int32_t));
    }




    av1_quantize_inv_quantize_ii(
        picture_control_set_ptr,
        coeff,
        0,
        quant_coeff,
        recon_coeff,
        qp,
        width,
        height,
        txsize,
        &eob[0],
        asm_type,
        y_count_non_zero_coeffs,
#if !PF_N2_SUPPORT
        0,
#endif
        0,
        component_type,
        bit_increment,
        tx_type,
        clean_sparse_coeff_flag);



}

/****************************************
 ************  Full loop ****************
****************************************/
void ProductFullLoop(
    ModeDecisionCandidateBuffer_t  *candidateBuffer,
    ModeDecisionContext_t          *context_ptr,
    PictureControlSet_t            *picture_control_set_ptr,
    uint32_t                          qp,
    uint32_t                          *y_count_non_zero_coeffs,
    uint64_t                         *y_coeff_bits,
    uint64_t                         *y_full_distortion)
{
    uint32_t                       tuOriginIndex;
    uint64_t                      y_full_cost;
    SequenceControlSet_t        *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    EbAsm                         asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;

    EbBool                      clean_sparse_coeff_flag = EB_FALSE;

    //    uint32_t   currentTuIndex,tuIt;
    uint64_t   y_tu_coeff_bits;
    uint64_t   tuFullDistortion[3][DIST_CALC_TOTAL];
    context_ptr->three_quad_energy = 0;
    uint32_t  txb_1d_offset = 0;
    uint32_t txb_itr = 0;
    for (txb_itr = 0; txb_itr < context_ptr->blk_geom->txb_count; txb_itr++)
    {
        uint16_t tx_org_x = context_ptr->blk_geom->tx_org_x[txb_itr];
        uint16_t tx_org_y = context_ptr->blk_geom->tx_org_y[txb_itr];
        tuOriginIndex = tx_org_x + (tx_org_y * candidateBuffer->residual_ptr->stride_y);
        y_tu_coeff_bits = 0;

        // Y: T Q iQ
        av1_estimate_transform(
            &(((int16_t*)candidateBuffer->residual_ptr->buffer_y)[tuOriginIndex]),
            candidateBuffer->residual_ptr->stride_y,
            &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->buffer_y)[txb_1d_offset]),
            NOT_USED_VALUE,
            context_ptr->blk_geom->txsize[txb_itr],
            &context_ptr->three_quad_energy,
            context_ptr->transform_inner_array_ptr,
            0,
            candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y],
            asm_type,
            PLANE_TYPE_Y,
#if PF_N2_SUPPORT
            DEFAULT_SHAPE);
#else
            context_ptr->pf_md_mode);
#endif
        av1_quantize_inv_quantize(
            picture_control_set_ptr,
            &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->buffer_y)[txb_1d_offset]),
            NOT_USED_VALUE,
            &(((int32_t*)candidateBuffer->residualQuantCoeffPtr->buffer_y)[txb_1d_offset]),
            &(((int32_t*)candidateBuffer->reconCoeffPtr->buffer_y)[txb_1d_offset]),
            qp,
            context_ptr->blk_geom->tx_width[txb_itr],
            context_ptr->blk_geom->tx_height[txb_itr],
            context_ptr->blk_geom->txsize[txb_itr],
            &candidateBuffer->candidate_ptr->eob[0][txb_itr],
            candidateBuffer->candidate_ptr->candidate_plane[0],
            asm_type,
            &(y_count_non_zero_coeffs[txb_itr]),
#if !PF_N2_SUPPORT
            context_ptr->pf_md_mode,
#endif
            0,
            COMPONENT_LUMA,
            BIT_INCREMENT_8BIT,
            candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y],
            clean_sparse_coeff_flag);

        candidateBuffer->candidate_ptr->quantized_dc[0] = (((int32_t*)candidateBuffer->residualQuantCoeffPtr->buffer_y)[txb_1d_offset]);

#if SPATIAL_SSE
        if (context_ptr->spatial_sse_full_loop) {
            EbPictureBufferDesc_t          *input_picture_ptr = picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr;
            uint32_t input_tu_origin_index = (context_ptr->sb_origin_x + tx_org_x + input_picture_ptr->origin_x) + ((context_ptr->sb_origin_y + tx_org_y + input_picture_ptr->origin_y) * input_picture_ptr->stride_y);
            uint32_t y_has_coeff           = y_count_non_zero_coeffs[txb_itr] > 0;

            if (y_has_coeff) {
                (void)context_ptr;
                uint8_t     *pred_buffer = &(candidateBuffer->prediction_ptr->buffer_y[tuOriginIndex]);
                uint8_t     *rec_buffer = &(candidateBuffer->recon_ptr->buffer_y[tuOriginIndex]);

                uint32_t j;

                for (j = 0; j < context_ptr->blk_geom->tx_height[txb_itr]; j++)
                    memcpy(rec_buffer + j * candidateBuffer->recon_ptr->stride_y, pred_buffer + j * candidateBuffer->prediction_ptr->stride_y, context_ptr->blk_geom->tx_width[txb_itr]);

                av1_inv_transform_recon8bit(

                    &(((int32_t*)candidateBuffer->reconCoeffPtr->buffer_y)[txb_1d_offset]),
                    rec_buffer,
                    candidateBuffer->recon_ptr->stride_y,
                    context_ptr->blk_geom->txsize[txb_itr],
                    candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y],
                    PLANE_TYPE_Y,
                    (uint16_t)candidateBuffer->candidate_ptr->eob[0][txb_itr]);

            }
            else {

                picture_copy8_bit(
                    candidateBuffer->prediction_ptr,
                    tuOriginIndex,
                    0,
                    candidateBuffer->recon_ptr,
                    tuOriginIndex,
                    0,
                    context_ptr->blk_geom->tx_width[txb_itr],
                    context_ptr->blk_geom->tx_height[txb_itr],
                    0,
                    0,
                    PICTURE_BUFFER_DESC_Y_FLAG,
                    asm_type);
            }

            tuFullDistortion[0][DIST_CALC_PREDICTION] = spatial_full_distortion_kernel_func_ptr_array[asm_type][context_ptr->blk_geom->bwidth == context_ptr->blk_geom->bheight ? Log2f(context_ptr->blk_geom->bwidth) - 2 : 6](
                input_picture_ptr->buffer_y + input_tu_origin_index,
                input_picture_ptr->stride_y,
                candidateBuffer->prediction_ptr->buffer_y + tuOriginIndex,
                candidateBuffer->prediction_ptr->stride_y,
                context_ptr->blk_geom->tx_width[txb_itr],
                context_ptr->blk_geom->tx_height[txb_itr]);

            tuFullDistortion[0][DIST_CALC_RESIDUAL] = spatial_full_distortion_kernel_func_ptr_array[asm_type][context_ptr->blk_geom->bwidth == context_ptr->blk_geom->bheight ? Log2f(context_ptr->blk_geom->bwidth) - 2 : 6](
                input_picture_ptr->buffer_y + input_tu_origin_index,
                input_picture_ptr->stride_y,
                &(((int8_t*)candidateBuffer->recon_ptr->buffer_y)[tuOriginIndex]),
                candidateBuffer->recon_ptr->stride_y,
                context_ptr->blk_geom->tx_width[txb_itr],
                context_ptr->blk_geom->tx_height[txb_itr]);

            tuFullDistortion[0][DIST_CALC_PREDICTION] <<= 4;
            tuFullDistortion[0][DIST_CALC_RESIDUAL] <<= 4;
        }
        else {
            // LUMA DISTORTION
            picture_full_distortion32_bits(
                context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr,
                txb_1d_offset,
                0,
                candidateBuffer->reconCoeffPtr,
                txb_1d_offset,
                0,
                context_ptr->blk_geom->tx_width[txb_itr],
                context_ptr->blk_geom->tx_height[txb_itr],
                NOT_USED_VALUE,
                NOT_USED_VALUE,
                tuFullDistortion[0],
                NOT_USED_VALUE,
                NOT_USED_VALUE,
                y_count_non_zero_coeffs[txb_itr],
                0,
                0,
                COMPONENT_LUMA,
                asm_type);


            tuFullDistortion[0][DIST_CALC_RESIDUAL] += context_ptr->three_quad_energy;
            tuFullDistortion[0][DIST_CALC_PREDICTION] += context_ptr->three_quad_energy;
            //assert(context_ptr->three_quad_energy == 0 && context_ptr->cu_stats->size < 64);
            TxSize    tx_size = context_ptr->blk_geom->txsize[0];
            int32_t shift = (MAX_TX_SCALE - av1_get_tx_scale(tx_size)) * 2;
            tuFullDistortion[0][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_RESIDUAL], shift);
            tuFullDistortion[0][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_PREDICTION], shift);
        }
#else

        // LUMA DISTORTION
        picture_full_distortion32_bits(
            context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr,
            txb_1d_offset,
            0,
            candidateBuffer->reconCoeffPtr,
            txb_1d_offset,
            0,
            context_ptr->blk_geom->tx_width[txb_itr],
            context_ptr->blk_geom->tx_height[txb_itr],
            NOT_USED_VALUE,
            NOT_USED_VALUE,
            tuFullDistortion[0],
            NOT_USED_VALUE,
            NOT_USED_VALUE,
            y_count_non_zero_coeffs[txb_itr],
            0,
            0,
            COMPONENT_LUMA,
            asm_type);


        tuFullDistortion[0][DIST_CALC_RESIDUAL] += context_ptr->three_quad_energy;
        tuFullDistortion[0][DIST_CALC_PREDICTION] += context_ptr->three_quad_energy;
        //assert(context_ptr->three_quad_energy == 0 && context_ptr->cu_stats->size < 64);
        TxSize    txSize = context_ptr->blk_geom->txsize[0];
        int32_t shift = (MAX_TX_SCALE - av1_get_tx_scale(txSize)) * 2;
        tuFullDistortion[0][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_RESIDUAL], shift);
        tuFullDistortion[0][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_PREDICTION], shift);
#endif

        //LUMA-ONLY
        Av1TuEstimateCoeffBits(
#if CABAC_UP
            0,//allow_update_cdf,
            NULL,//FRAME_CONTEXT *ec_ctx,
#endif
            picture_control_set_ptr,
            candidateBuffer,
            context_ptr->cu_ptr,
            txb_1d_offset,
            0,
            context_ptr->coeff_est_entropy_coder_ptr,
            candidateBuffer->residualQuantCoeffPtr,
            y_count_non_zero_coeffs[txb_itr],
            0,
            0,
            &y_tu_coeff_bits,
            &y_tu_coeff_bits,
            &y_tu_coeff_bits,
            context_ptr->blk_geom->txsize[0],
            context_ptr->blk_geom->txsize_uv[0],
            COMPONENT_LUMA,
            asm_type);



        //TODO: fix cbf decision
        av1_tu_calc_cost_luma(
            context_ptr->cu_ptr->luma_txb_skip_context,//this should be updated here.
            candidateBuffer->candidate_ptr,
            txb_itr,
            context_ptr->blk_geom->txsize[0],
            y_count_non_zero_coeffs[txb_itr],
            tuFullDistortion[0],      //gets updated inside based on cbf decision
            &y_tu_coeff_bits,            //gets updated inside based on cbf decision
            &y_full_cost,
            context_ptr->full_lambda);


        (*y_coeff_bits) += y_tu_coeff_bits;

        y_full_distortion[DIST_CALC_RESIDUAL] += tuFullDistortion[0][DIST_CALC_RESIDUAL];
        y_full_distortion[DIST_CALC_PREDICTION] += tuFullDistortion[0][DIST_CALC_PREDICTION];
        txb_1d_offset += context_ptr->blk_geom->tx_width[txb_itr] * context_ptr->blk_geom->tx_height[txb_itr];



    }
}
// T1
uint8_t allowed_tx_set_a[TX_SIZES_ALL][TX_TYPES] = {
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    1,    0,    1},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    1,    0,    1},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    1,    0,    1,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0} };

uint8_t allowed_tx_set_b[TX_SIZES_ALL][TX_TYPES] = {
{1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    1,    1,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0},
{1,    1,    1,    1,    0,    0,    0,    0,    0,    0,    1,    1,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    1,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{0,    0,    1,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    1,    0,    1,    0,    0,    0,    0,    0,    1,    1,    1,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
{1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0}
};

void ProductFullLoopTxSearch(
    ModeDecisionCandidateBuffer_t  *candidateBuffer,
    ModeDecisionContext_t          *context_ptr,
    PictureControlSet_t            *picture_control_set_ptr)
{
    uint32_t                       tuOriginIndex;
    SequenceControlSet_t          *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    EbAsm                          asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;
    EbBool                         clean_sparse_coeff_flag = EB_FALSE;
    uint64_t                       y_tu_coeff_bits;
    uint64_t                       tuFullDistortion[3][DIST_CALC_TOTAL];
    int32_t                        plane = 0;
    const int32_t                  is_inter = (candidateBuffer->candidate_ptr->type == INTER_MODE || candidateBuffer->candidate_ptr->use_intrabc) ? EB_TRUE : EB_FALSE;
    uint64_t                       bestFullCost = UINT64_MAX;
    uint64_t                       y_full_cost = MAX_CU_COST;
    uint32_t                       yCountNonZeroCoeffsTemp;
    TxType                         txk_start = DCT_DCT;
    TxType                         txk_end = TX_TYPES;
    TxType                         tx_type;
    int32_t                        txb_itr = 0;
    TxSize                         txSize = context_ptr->blk_geom->txsize[txb_itr];
    const TxSetType                tx_set_type =
        get_ext_tx_set_type(txSize, is_inter, picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);


    int32_t allowed_tx_mask[TX_TYPES] = { 0 };  // 1: allow; 0: skip.
    int32_t allowed_tx_num = 0;
    TxType uv_tx_type = DCT_DCT;

    for (int32_t tx_type_index = txk_start; tx_type_index < txk_end; ++tx_type_index) {
        tx_type = (TxType)tx_type_index;
        allowed_tx_mask[tx_type] = 1;
        if (plane == 0) {
            if (allowed_tx_mask[tx_type]) {
                const TxType ref_tx_type = ((!av1_ext_tx_used[tx_set_type][tx_type]) || txsize_sqr_up_map[txSize] > TX_32X32) ? DCT_DCT : tx_type;
                if (tx_type != ref_tx_type) {

                    allowed_tx_mask[tx_type] = 0;
                }
            }
        }

        allowed_tx_num += allowed_tx_mask[tx_type];
    }
    // Need to have at least one transform type allowed.
    if (allowed_tx_num == 0) {
        allowed_tx_mask[plane ? uv_tx_type : DCT_DCT] = 1;
    }
    TxType best_tx_type = DCT_DCT;
    for (int32_t tx_type_index = txk_start; tx_type_index < txk_end; ++tx_type_index) {
        tx_type = (TxType)tx_type_index;
        if (!allowed_tx_mask[tx_type]) continue;
        if (picture_control_set_ptr->parent_pcs_ptr->tx_search_reduced_set)
            if (!allowed_tx_set_a[txSize][tx_type]) continue;

        context_ptr->three_quad_energy = 0;
        uint32_t txb_itr = 0;
        for (txb_itr = 0; txb_itr < context_ptr->blk_geom->txb_count; txb_itr++)


        {
            tuOriginIndex = context_ptr->blk_geom->origin_x + (context_ptr->blk_geom->origin_y * candidateBuffer->residual_ptr->stride_y);
            y_tu_coeff_bits = 0;

            // Y: T Q iQ
            av1_estimate_transform(
                &(((int16_t*)candidateBuffer->residual_ptr->buffer_y)[tuOriginIndex]),
                candidateBuffer->residual_ptr->stride_y,

                &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->buffer_y)[tuOriginIndex]),
                NOT_USED_VALUE,
                context_ptr->blk_geom->txsize[txb_itr],
                &context_ptr->three_quad_energy,
                context_ptr->transform_inner_array_ptr,
                0,
                tx_type,
                asm_type,
                PLANE_TYPE_Y,
                context_ptr->pf_md_mode);

            av1_quantize_inv_quantize(
                picture_control_set_ptr,
                &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->buffer_y)[tuOriginIndex]),
                NOT_USED_VALUE,
                &(((int32_t*)candidateBuffer->residualQuantCoeffPtr->buffer_y)[tuOriginIndex]),
                &(((int32_t*)candidateBuffer->reconCoeffPtr->buffer_y)[tuOriginIndex]),
                context_ptr->cu_ptr->qp,
                context_ptr->blk_geom->bwidth,
                context_ptr->blk_geom->bheight,
                context_ptr->blk_geom->txsize[txb_itr],
                &candidateBuffer->candidate_ptr->eob[0][txb_itr],
                candidateBuffer->candidate_ptr->candidate_plane[0],
                asm_type,
                &yCountNonZeroCoeffsTemp,
#if !PF_N2_SUPPORT
                context_ptr->pf_md_mode,
#endif
                0,
                COMPONENT_LUMA,
                BIT_INCREMENT_8BIT,
                tx_type,
                clean_sparse_coeff_flag);

            candidateBuffer->candidate_ptr->quantized_dc[0] = (((int32_t*)candidateBuffer->residualQuantCoeffPtr->buffer_y)[tuOriginIndex]);


            //tx_type not equal to DCT_DCT and no coeff is not an acceptable option in AV1.
            if (yCountNonZeroCoeffsTemp == 0 && tx_type != DCT_DCT) {

                continue;
            }



            // LUMA DISTORTION
            picture_full_distortion32_bits(
                context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr,
                tuOriginIndex,
                0,
                candidateBuffer->reconCoeffPtr,
                tuOriginIndex,
                0,
                context_ptr->blk_geom->bwidth,
                context_ptr->blk_geom->bheight,
                context_ptr->blk_geom->bwidth_uv,
                context_ptr->blk_geom->bheight_uv,
                tuFullDistortion[0],
                tuFullDistortion[0],
                tuFullDistortion[0],
                yCountNonZeroCoeffsTemp,
                0,
                0,
                COMPONENT_LUMA,
                asm_type);

            tuFullDistortion[0][DIST_CALC_RESIDUAL] += context_ptr->three_quad_energy;
            tuFullDistortion[0][DIST_CALC_PREDICTION] += context_ptr->three_quad_energy;

            int32_t shift = (MAX_TX_SCALE - av1_get_tx_scale(txSize)) * 2;
            tuFullDistortion[0][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_RESIDUAL], shift);
            tuFullDistortion[0][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_PREDICTION], shift);
            candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = tx_type;
            //LUMA-ONLY
            Av1TuEstimateCoeffBits(
#if CABAC_UP
                0,//allow_update_cdf,
                NULL,//FRAME_CONTEXT *ec_ctx,
#endif
                picture_control_set_ptr,
                candidateBuffer,
                context_ptr->cu_ptr,
                tuOriginIndex,
                0,
                context_ptr->coeff_est_entropy_coder_ptr,
                candidateBuffer->residualQuantCoeffPtr,
                yCountNonZeroCoeffsTemp,
                0,
                0,
                &y_tu_coeff_bits,
                &y_tu_coeff_bits,
                &y_tu_coeff_bits,
                context_ptr->blk_geom->txsize[txb_itr],
                context_ptr->blk_geom->txsize_uv[txb_itr],
                COMPONENT_LUMA,
                asm_type);

            av1_tu_calc_cost_luma(
                context_ptr->cu_ptr->luma_txb_skip_context,
                candidateBuffer->candidate_ptr,
                txb_itr,
                context_ptr->blk_geom->txsize[txb_itr],
                yCountNonZeroCoeffsTemp,
                tuFullDistortion[0],
                &y_tu_coeff_bits,
                &y_full_cost,
                context_ptr->full_lambda);

        }

        if (y_full_cost < bestFullCost) {
            bestFullCost = y_full_cost;
            best_tx_type = tx_type;

        }

        //if (cpi->sf.adaptive_txb_search_level) {
        //    if ((best_rd - (best_rd >> cpi->sf.adaptive_txb_search_level)) >
        //        ref_best_rd) {
        //        break;
        //    }
        //}
        //// Skip transform type search when we found the block has been quantized to
        //// all zero and at the same time, it has better rdcost than doing transform.
        //if (cpi->sf.tx_type_search.skip_tx_search && !best_eob) break;


    }
    candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = best_tx_type;
    // For Inter blocks, transform type of chroma follows luma transfrom type
    if (is_inter)
        candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV] = candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y];
}

void encode_pass_tx_search(
    PictureControlSet_t            *picture_control_set_ptr,
    EncDecContext_t                *context_ptr,
    LargestCodingUnit_t            *sb_ptr,
    uint32_t                       cbQp,
    EbPictureBufferDesc_t          *coeffSamplesTB,          
    EbPictureBufferDesc_t          *residual16bit,           
    EbPictureBufferDesc_t          *transform16bit,          
    EbPictureBufferDesc_t          *inverse_quant_buffer,
    int16_t                        *transformScratchBuffer,
    EbAsm                          asm_type,
    uint32_t                       *count_non_zero_coeffs,
    uint32_t                       component_mask,
    uint32_t                       use_delta_qp,
    uint32_t                       dZoffset,
    uint16_t                       *eob,
    MacroblockPlane                *candidate_plane){

    (void)dZoffset;
    (void)use_delta_qp;
    (void)cbQp;
    UNUSED(count_non_zero_coeffs);
    UNUSED(component_mask);

    CodingUnit_t          *cu_ptr = context_ptr->cu_ptr;
    TransformUnit_t       *txb_ptr = &cu_ptr->transform_unit_array[context_ptr->txb_itr];
    uint32_t               qp = cu_ptr->qp;
    const uint32_t         scratchLumaOffset = context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr] + context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr] * SB_STRIDE_Y;
    const uint32_t         coeff1dOffset = context_ptr->coded_area_sb;

    EbBool                 clean_sparse_coeff_flag = EB_FALSE;
    uint64_t               y_tu_coeff_bits;
    uint64_t               tuFullDistortion[3][DIST_CALC_TOTAL];
    const int32_t          is_inter = context_ptr->is_inter;
    uint64_t               bestFullCost = UINT64_MAX;
    uint64_t               y_full_cost;
    uint32_t               yCountNonZeroCoeffsTemp;
    TxType                 txk_start = DCT_DCT;
    TxType                 txk_end = TX_TYPES;
    TxType                 tx_type;
    TxSize                 txSize = context_ptr->blk_geom->txsize[context_ptr->txb_itr];
    const TxSetType        tx_set_type =
        get_ext_tx_set_type(txSize, is_inter, picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);

    TxType best_tx_type = DCT_DCT;

    for (int32_t tx_type_index = txk_start; tx_type_index < txk_end; ++tx_type_index) {
        tx_type = (TxType)tx_type_index;

        if(picture_control_set_ptr->parent_pcs_ptr->tx_search_reduced_set)
            if (!allowed_tx_set_a[txSize][tx_type]) continue;

        const int32_t eset = get_ext_tx_set(txSize, is_inter, picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);
        // eset == 0 should correspond to a set with only DCT_DCT and there
        // is no need to send the tx_type
        if (eset <= 0) continue;
        if (av1_ext_tx_used[tx_set_type][tx_type] == 0) continue;

        context_ptr->three_quad_energy = 0;

        y_tu_coeff_bits = 0;


        av1_estimate_transform(
            ((int16_t*)residual16bit->buffer_y) + scratchLumaOffset,
            residual16bit->stride_y,
            ((tran_low_t*)transform16bit->buffer_y) + coeff1dOffset,
            NOT_USED_VALUE,
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            &context_ptr->three_quad_energy,
            transformScratchBuffer,
            BIT_INCREMENT_8BIT,
            tx_type,
            asm_type,
            PLANE_TYPE_Y,
#if PF_N2_SUPPORT
            DEFAULT_SHAPE);
#else
            context_ptr->trans_coeff_shape_luma);
#endif

        av1_quantize_inv_quantize(
            sb_ptr->picture_control_set_ptr,
            ((tran_low_t*)transform16bit->buffer_y) + coeff1dOffset,
            NOT_USED_VALUE,
            ((int32_t*)coeffSamplesTB->buffer_y) + coeff1dOffset,
            ((int32_t*)inverse_quant_buffer->buffer_y) + coeff1dOffset,
            qp,
            context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height[context_ptr->txb_itr],
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            &eob[0],
            candidate_plane[0],
            asm_type,
            &yCountNonZeroCoeffsTemp,
#if !PF_N2_SUPPORT
            0,
#endif
            0,
            COMPONENT_LUMA,
            BIT_INCREMENT_8BIT,
            tx_type,
            clean_sparse_coeff_flag);

        //tx_type not equal to DCT_DCT and no coeff is not an acceptable option in AV1.
        if (yCountNonZeroCoeffsTemp == 0 && tx_type != DCT_DCT) {
            continue;
        }


        // LUMA DISTORTION
        picture_full_distortion32_bits(
            transform16bit,
            coeff1dOffset,
            0,
            inverse_quant_buffer,
            coeff1dOffset,
            0,
            context_ptr->blk_geom->bwidth,
            context_ptr->blk_geom->bheight,
            context_ptr->blk_geom->bwidth_uv,
            context_ptr->blk_geom->bheight_uv,
            tuFullDistortion[0],
            tuFullDistortion[0],
            tuFullDistortion[0],
            yCountNonZeroCoeffsTemp,
            0,
            0,
            COMPONENT_LUMA,
            asm_type);

        tuFullDistortion[0][DIST_CALC_RESIDUAL] += context_ptr->three_quad_energy;
        tuFullDistortion[0][DIST_CALC_PREDICTION] += context_ptr->three_quad_energy;

        int32_t shift = (MAX_TX_SCALE - av1_get_tx_scale(txSize)) * 2;
        tuFullDistortion[0][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_RESIDUAL], shift);
        tuFullDistortion[0][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_PREDICTION], shift);
        txb_ptr->transform_type[PLANE_TYPE_Y] = tx_type;

        //LUMA-ONLY

        ModeDecisionCandidateBuffer_t         **candidateBufferPtrArrayBase = context_ptr->md_context->candidate_buffer_ptr_array;
        ModeDecisionCandidateBuffer_t         **candidate_buffer_ptr_array = &(candidateBufferPtrArrayBase[0]);
        ModeDecisionCandidateBuffer_t          *candidateBuffer;

        // Set the Candidate Buffer
        candidateBuffer = candidate_buffer_ptr_array[0];
        // Rate estimation function uses the values from CandidatePtr. The right values are copied from cu_ptr to CandidatePtr
        candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_Y];
        candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_UV];
        EntropyCoder_t  *coeff_est_entropy_coder_ptr = picture_control_set_ptr->coeff_est_entropy_coder_ptr;
        candidateBuffer->candidate_ptr->type = cu_ptr->prediction_mode_flag;
        candidateBuffer->candidate_ptr->pred_mode = cu_ptr->pred_mode;

        const uint32_t coeff1dOffset = context_ptr->coded_area_sb;

        Av1TuEstimateCoeffBits(
#if CABAC_UP
            0,//allow_update_cdf,
            NULL,//FRAME_CONTEXT *ec_ctx,
#endif
            picture_control_set_ptr,
            candidateBuffer,
            context_ptr->cu_ptr,
            coeff1dOffset,
            0,
            coeff_est_entropy_coder_ptr,
            coeffSamplesTB,
            yCountNonZeroCoeffsTemp,
            0,
            0,
            &y_tu_coeff_bits,
            &y_tu_coeff_bits,
            &y_tu_coeff_bits,
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
            COMPONENT_LUMA,
            asm_type);

        av1_tu_calc_cost_luma(
            context_ptr->cu_ptr->luma_txb_skip_context,
            candidateBuffer->candidate_ptr,
            context_ptr->txb_itr,
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            yCountNonZeroCoeffsTemp,
            tuFullDistortion[0],
            &y_tu_coeff_bits,
            &y_full_cost,
            context_ptr->full_lambda);

        if (y_full_cost < bestFullCost) {
            bestFullCost = y_full_cost;
            best_tx_type = tx_type;
        }


    }

    txb_ptr->transform_type[PLANE_TYPE_Y] = best_tx_type;

    // For Inter blocks, transform type of chroma follows luma transfrom type
    if (is_inter)
        txb_ptr->transform_type[PLANE_TYPE_UV] = txb_ptr->transform_type[PLANE_TYPE_Y];

}

void encode_pass_tx_search_hbd(
    PictureControlSet_t            *picture_control_set_ptr,
    EncDecContext_t                *context_ptr,
    LargestCodingUnit_t            *sb_ptr,
    uint32_t                       cbQp,
    EbPictureBufferDesc_t          *coeffSamplesTB,
    EbPictureBufferDesc_t          *residual16bit,
    EbPictureBufferDesc_t          *transform16bit,
    EbPictureBufferDesc_t          *inverse_quant_buffer,
    int16_t                        *transformScratchBuffer,
    EbAsm                          asm_type,
    uint32_t                       *count_non_zero_coeffs,
    uint32_t                       component_mask,
    uint32_t                       use_delta_qp,
    uint32_t                       dZoffset,
    uint16_t                       *eob,
    MacroblockPlane                *candidate_plane){

    (void)dZoffset;
    (void)use_delta_qp;
    (void)cbQp;
    UNUSED(component_mask);
    UNUSED(count_non_zero_coeffs);

    CodingUnit_t    *cu_ptr               = context_ptr->cu_ptr;
    TransformUnit_t *txb_ptr              = &cu_ptr->transform_unit_array[context_ptr->txb_itr];
    uint32_t         qp                   = cu_ptr->qp;
    const uint32_t   scratchLumaOffset    = context_ptr->blk_geom->origin_x + context_ptr->blk_geom->origin_y * SB_STRIDE_Y;
    const uint32_t   coeff1dOffset        = context_ptr->coded_area_sb;
    EbBool           clean_sparse_coeff_flag = EB_FALSE;

    //Update QP for Quant
    qp += QP_BD_OFFSET;
    uint64_t                    y_tu_coeff_bits;
    uint64_t                    tuFullDistortion[3][DIST_CALC_TOTAL];
    const int32_t               is_inter = context_ptr->is_inter;
    uint64_t                    bestFullCost = UINT64_MAX;
    uint64_t                    y_full_cost;
    uint32_t                    yCountNonZeroCoeffsTemp;
    TxType                      txk_start = DCT_DCT;
    TxType                      txk_end = TX_TYPES;
    TxType                      tx_type;
    TxSize                      txSize = context_ptr->blk_geom->txsize[context_ptr->txb_itr];
    const TxSetType             tx_set_type =
        get_ext_tx_set_type(txSize, is_inter, picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);

    TxType best_tx_type = DCT_DCT;

    for (int32_t tx_type_index = txk_start; tx_type_index < txk_end; ++tx_type_index) {
        tx_type = (TxType)tx_type_index;
        ////if (!allowed_tx_mask[tx_type]) continue;
        if (picture_control_set_ptr->parent_pcs_ptr->tx_search_reduced_set)
            if (!allowed_tx_set_a[txSize][tx_type]) continue;

        const int32_t eset = get_ext_tx_set(txSize, is_inter, picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);
        // eset == 0 should correspond to a set with only DCT_DCT and there
        // is no need to send the tx_type
        if (eset <= 0) continue;
        if (av1_ext_tx_used[tx_set_type][tx_type] == 0) continue;

        context_ptr->three_quad_energy = 0;


        y_tu_coeff_bits = 0;


        av1_estimate_transform(
            ((int16_t*)residual16bit->buffer_y) + scratchLumaOffset,
            residual16bit->stride_y,
            ((tran_low_t*)transform16bit->buffer_y) + coeff1dOffset,
            NOT_USED_VALUE,
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            &context_ptr->three_quad_energy,
            transformScratchBuffer,
            BIT_INCREMENT_10BIT,
            tx_type,
            asm_type,
            PLANE_TYPE_Y,
#if PF_N2_SUPPORT
            DEFAULT_SHAPE);
#else
            context_ptr->trans_coeff_shape_luma);
#endif
        av1_quantize_inv_quantize(
            sb_ptr->picture_control_set_ptr,
            ((int32_t*)transform16bit->buffer_y) + coeff1dOffset,
            NOT_USED_VALUE,
            ((int32_t*)coeffSamplesTB->buffer_y) + coeff1dOffset,
            ((int32_t*)inverse_quant_buffer->buffer_y) + coeff1dOffset,
            qp,
            context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height[context_ptr->txb_itr],
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            &eob[0],
            candidate_plane[0],
            asm_type,
            &yCountNonZeroCoeffsTemp,
#if !PF_N2_SUPPORT
            0,
#endif
            0,
            COMPONENT_LUMA,
            BIT_INCREMENT_10BIT,
            tx_type,
            clean_sparse_coeff_flag);

        //tx_type not equal to DCT_DCT and no coeff is not an acceptable option in AV1.
        if (yCountNonZeroCoeffsTemp == 0 && tx_type != DCT_DCT) {
            continue;
        }


        // LUMA DISTORTION
        picture_full_distortion32_bits(
            transform16bit,
            coeff1dOffset,
            0,
            inverse_quant_buffer,
            coeff1dOffset,
            0,
            context_ptr->blk_geom->bwidth,
            context_ptr->blk_geom->bheight,
            context_ptr->blk_geom->bwidth_uv,
            context_ptr->blk_geom->bheight_uv,
            tuFullDistortion[0],
            tuFullDistortion[0],
            tuFullDistortion[0],
            yCountNonZeroCoeffsTemp,
            0,
            0,
            COMPONENT_LUMA,
            asm_type);

        tuFullDistortion[0][DIST_CALC_RESIDUAL] += context_ptr->three_quad_energy;
        tuFullDistortion[0][DIST_CALC_PREDICTION] += context_ptr->three_quad_energy;

        int32_t shift = (MAX_TX_SCALE - av1_get_tx_scale(txSize)) * 2;
        tuFullDistortion[0][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_RESIDUAL], shift);
        tuFullDistortion[0][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(tuFullDistortion[0][DIST_CALC_PREDICTION], shift);
        txb_ptr->transform_type[PLANE_TYPE_Y] = tx_type;

        //LUMA-ONLY

        ModeDecisionCandidateBuffer_t         **candidateBufferPtrArrayBase = context_ptr->md_context->candidate_buffer_ptr_array;
        ModeDecisionCandidateBuffer_t         **candidate_buffer_ptr_array = &(candidateBufferPtrArrayBase[0]);
        ModeDecisionCandidateBuffer_t          *candidateBuffer;

        // Set the Candidate Buffer
        candidateBuffer = candidate_buffer_ptr_array[0];
        // Rate estimation function uses the values from CandidatePtr. The right values are copied from cu_ptr to CandidatePtr
        candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_Y];
        candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_UV];
        EntropyCoder_t  *coeff_est_entropy_coder_ptr = picture_control_set_ptr->coeff_est_entropy_coder_ptr;
        candidateBuffer->candidate_ptr->type = cu_ptr->prediction_mode_flag;
        candidateBuffer->candidate_ptr->pred_mode = cu_ptr->pred_mode;

        const uint32_t coeff1dOffset = context_ptr->coded_area_sb;

        Av1TuEstimateCoeffBits(
#if CABAC_UP
            0,//allow_update_cdf,
            NULL,//FRAME_CONTEXT *ec_ctx,
#endif
            picture_control_set_ptr,
            candidateBuffer,
            context_ptr->cu_ptr,
            coeff1dOffset,
            0,
            coeff_est_entropy_coder_ptr,
            coeffSamplesTB,
            yCountNonZeroCoeffsTemp,
            0,
            0,
            &y_tu_coeff_bits,
            &y_tu_coeff_bits,
            &y_tu_coeff_bits,
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
            COMPONENT_LUMA,
            asm_type);

        av1_tu_calc_cost_luma(
            context_ptr->cu_ptr->luma_txb_skip_context,
            candidateBuffer->candidate_ptr,
            context_ptr->txb_itr,
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            yCountNonZeroCoeffsTemp,
            tuFullDistortion[0],
            &y_tu_coeff_bits,
            &y_full_cost,
            context_ptr->full_lambda);

        if (y_full_cost < bestFullCost) {
            bestFullCost = y_full_cost;
            best_tx_type = tx_type;
        }


    }


    txb_ptr->transform_type[PLANE_TYPE_Y] = best_tx_type;

    // For Inter blocks, transform type of chroma follows luma transfrom type
    if (is_inter)
        txb_ptr->transform_type[PLANE_TYPE_UV] = txb_ptr->transform_type[PLANE_TYPE_Y];

}
/****************************************
 ************  Full loop ****************
****************************************/
void FullLoop_R(
    LargestCodingUnit_t            *sb_ptr,
    ModeDecisionCandidateBuffer_t  *candidateBuffer,
    ModeDecisionContext_t          *context_ptr,
    EbPictureBufferDesc_t          *input_picture_ptr,
    PictureControlSet_t            *picture_control_set_ptr,
    uint32_t                          component_mask,
    uint32_t                          cbQp,
    uint32_t                          crQp,
    uint32_t                          *cb_count_non_zero_coeffs,
    uint32_t                          *cr_count_non_zero_coeffs)
{
    (void)sb_ptr;
    (void)crQp;
    (void)input_picture_ptr;
    int16_t                *chromaResidualPtr;
    uint32_t                 tuOriginIndex;
    UNUSED(tuOriginIndex);
    uint32_t                 tuCbOriginIndex;
    uint32_t                 tuCrOriginIndex;
    uint32_t                 tuCount;
    uint32_t                 txb_itr;
    uint32_t                 txb_origin_x;
    uint32_t                 txb_origin_y;

    // EbPictureBufferDesc_t         * tuTransCoeffTmpPtr;
     //EbPictureBufferDesc_t         * tuQuantCoeffTmpPtr;

    SequenceControlSet_t    *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    EbAsm     asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;

    EbBool clean_sparse_coeff_flag = EB_FALSE;
    context_ptr->three_quad_energy = 0;

    tuCount = context_ptr->blk_geom->txb_count;
    uint32_t  txb_1d_offset = 0;

    txb_itr = 0;
    do {

        txb_origin_x = context_ptr->blk_geom->tx_org_x[txb_itr];
        txb_origin_y = context_ptr->blk_geom->tx_org_y[txb_itr];


        // NADER - TU
        tuOriginIndex = txb_origin_x + txb_origin_y * candidateBuffer->residualQuantCoeffPtr->stride_y;
        tuCbOriginIndex = (((txb_origin_x >> 3) << 3) + (((txb_origin_y >> 3) << 3) * candidateBuffer->residualQuantCoeffPtr->strideCb)) >> 1;
        tuCrOriginIndex = (((txb_origin_x >> 3) << 3) + (((txb_origin_y >> 3) << 3) * candidateBuffer->residualQuantCoeffPtr->strideCr)) >> 1;

        //    This function replaces the previous Intra Chroma mode if the LM fast
            //    cost is better.
            //    *Note - this might require that we have inv transform in the loop
#if !PF_N2_SUPPORT
        EbPfMode    correctedPFMode = PF_OFF;
#endif
        if (component_mask & PICTURE_BUFFER_DESC_Cb_FLAG) {
            // Configure the Chroma Residual Ptr

            chromaResidualPtr = //(candidateBuffer->candidate_ptr->type  == INTRA_MODE )?
                  //&(((int16_t*) candidateBuffer->intraChromaResidualPtr->bufferCb)[tuChromaOriginIndex]):
                &(((int16_t*)candidateBuffer->residual_ptr->bufferCb)[tuCbOriginIndex]);


            // Cb Transform
            av1_estimate_transform(
                chromaResidualPtr,
                candidateBuffer->residual_ptr->strideCb,
                &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->bufferCb)[txb_1d_offset]),
                NOT_USED_VALUE,
                context_ptr->blk_geom->txsize_uv[txb_itr],
                &context_ptr->three_quad_energy,
                context_ptr->transform_inner_array_ptr,
                0,
                candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV],
                asm_type,
                PLANE_TYPE_UV,
#if PF_N2_SUPPORT
                DEFAULT_SHAPE);
#else
                correctedPFMode);
#endif
            av1_quantize_inv_quantize(
                picture_control_set_ptr,
                &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->bufferCb)[txb_1d_offset]),
                NOT_USED_VALUE,
                &(((int32_t*)candidateBuffer->residualQuantCoeffPtr->bufferCb)[txb_1d_offset]),
                &(((int32_t*)candidateBuffer->reconCoeffPtr->bufferCb)[txb_1d_offset]),
                cbQp,
                context_ptr->blk_geom->tx_width_uv[txb_itr],
                context_ptr->blk_geom->tx_height_uv[txb_itr],
                context_ptr->blk_geom->txsize_uv[txb_itr],
                &candidateBuffer->candidate_ptr->eob[1][txb_itr],
                candidateBuffer->candidate_ptr->candidate_plane[1],
                asm_type,
                &(cb_count_non_zero_coeffs[txb_itr]),
#if !PF_N2_SUPPORT
                context_ptr->pf_md_mode,
#endif
                0,
                COMPONENT_CHROMA_CB,
                BIT_INCREMENT_8BIT,
                candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV],
                clean_sparse_coeff_flag);
            candidateBuffer->candidate_ptr->quantized_dc[1] = (((int32_t*)candidateBuffer->residualQuantCoeffPtr->bufferCb)[txb_1d_offset]);
#if SPATIAL_SSE
            if (context_ptr->spatial_sse_full_loop) {
                uint32_t cb_has_coeff = cb_count_non_zero_coeffs[txb_itr] > 0;

                if (cb_has_coeff) {
                    uint8_t     *pred_buffer = &(candidateBuffer->prediction_ptr->bufferCb[tuCbOriginIndex]);
                    uint8_t     *rec_buffer = &(candidateBuffer->recon_ptr->bufferCb[tuCbOriginIndex]);

                    uint32_t j;

                    for (j = 0; j < context_ptr->blk_geom->tx_height_uv[txb_itr]; j++)
                        memcpy(rec_buffer + j * candidateBuffer->recon_ptr->strideCb, pred_buffer + j * candidateBuffer->prediction_ptr->strideCb, context_ptr->blk_geom->tx_width_uv[txb_itr]);

                    av1_inv_transform_recon8bit(
                        &(((int32_t*)candidateBuffer->reconCoeffPtr->bufferCb)[txb_1d_offset]),
                        rec_buffer,
                        candidateBuffer->recon_ptr->strideCb,
                        context_ptr->blk_geom->txsize_uv[txb_itr],
                        candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV],
                        PLANE_TYPE_UV,
                        (uint16_t)candidateBuffer->candidate_ptr->eob[1][txb_itr]);

                }
                else {

                    picture_copy8_bit(
                        candidateBuffer->prediction_ptr,
                        0,
                        tuCbOriginIndex,
                        candidateBuffer->recon_ptr,
                        0,
                        tuCbOriginIndex,
                        0,
                        0,
                        context_ptr->blk_geom->tx_width_uv[txb_itr],
                        context_ptr->blk_geom->tx_height_uv[txb_itr],
                        PICTURE_BUFFER_DESC_Cb_FLAG,
                        asm_type);
                }
            }
#endif        
        }


        if (component_mask & PICTURE_BUFFER_DESC_Cr_FLAG) {
            // Configure the Chroma Residual Ptr

            chromaResidualPtr = //(candidateBuffer->candidate_ptr->type  == INTRA_MODE )?
                //&(((int16_t*) candidateBuffer->intraChromaResidualPtr->bufferCr)[tuChromaOriginIndex]):
                &(((int16_t*)candidateBuffer->residual_ptr->bufferCr)[tuCrOriginIndex]);

            // Cr Transform
            av1_estimate_transform(
                chromaResidualPtr,
                candidateBuffer->residual_ptr->strideCr,
                &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->bufferCr)[txb_1d_offset]),
                NOT_USED_VALUE,
                context_ptr->blk_geom->txsize_uv[txb_itr],
                &context_ptr->three_quad_energy,
                context_ptr->transform_inner_array_ptr,
                0,
                candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV],
                asm_type,
                PLANE_TYPE_UV,
#if PF_N2_SUPPORT
                DEFAULT_SHAPE);
#else
                correctedPFMode);
#endif

            av1_quantize_inv_quantize(
                picture_control_set_ptr,
                &(((int32_t*)context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr->bufferCr)[txb_1d_offset]),
                NOT_USED_VALUE,
                &(((int32_t*)candidateBuffer->residualQuantCoeffPtr->bufferCr)[txb_1d_offset]),
                &(((int32_t*)candidateBuffer->reconCoeffPtr->bufferCr)[txb_1d_offset]),
                cbQp,
                context_ptr->blk_geom->tx_width_uv[txb_itr],
                context_ptr->blk_geom->tx_height_uv[txb_itr],
                context_ptr->blk_geom->txsize_uv[txb_itr],
                &candidateBuffer->candidate_ptr->eob[2][txb_itr],
                candidateBuffer->candidate_ptr->candidate_plane[2],
                asm_type,
                &(cr_count_non_zero_coeffs[txb_itr]),
#if !PF_N2_SUPPORT
                context_ptr->pf_md_mode,
#endif
                0,
                COMPONENT_CHROMA_CR,
                BIT_INCREMENT_8BIT,
                candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV],
                clean_sparse_coeff_flag);
            candidateBuffer->candidate_ptr->quantized_dc[2] = (((int32_t*)candidateBuffer->residualQuantCoeffPtr->bufferCr)[txb_1d_offset]);
#if SPATIAL_SSE
            if (context_ptr->spatial_sse_full_loop) {
                uint32_t cr_has_coeff = cr_count_non_zero_coeffs[txb_itr] > 0;

                if (cr_has_coeff) {
                    uint8_t     *pred_buffer = &(candidateBuffer->prediction_ptr->bufferCr[tuCbOriginIndex]);
                    uint8_t     *rec_buffer = &(candidateBuffer->recon_ptr->bufferCr[tuCbOriginIndex]);

                    uint32_t j;

                    for (j = 0; j < context_ptr->blk_geom->tx_height_uv[txb_itr]; j++)
                        memcpy(rec_buffer + j * candidateBuffer->recon_ptr->strideCr, pred_buffer + j * candidateBuffer->prediction_ptr->strideCr, context_ptr->blk_geom->tx_width_uv[txb_itr]);

                    av1_inv_transform_recon8bit(
                        &(((int32_t*)candidateBuffer->reconCoeffPtr->bufferCr)[txb_1d_offset]),
                        rec_buffer,
                        candidateBuffer->recon_ptr->strideCr,
                        context_ptr->blk_geom->txsize_uv[txb_itr],
                        candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV],
                        PLANE_TYPE_UV,
                        (uint16_t)candidateBuffer->candidate_ptr->eob[2][txb_itr]);

                }
                else {

                    picture_copy8_bit(
                        candidateBuffer->prediction_ptr,
                        0,
                        tuCbOriginIndex,
                        candidateBuffer->recon_ptr,
                        0,
                        tuCbOriginIndex,
                        0,
                        0,
                        context_ptr->blk_geom->tx_width_uv[txb_itr],
                        context_ptr->blk_geom->tx_height_uv[txb_itr],
                        PICTURE_BUFFER_DESC_Cr_FLAG,
                        asm_type);
                }
            }
#endif        
        
        }

        txb_1d_offset += context_ptr->blk_geom->tx_width_uv[txb_itr] * context_ptr->blk_geom->tx_height_uv[txb_itr];


        ++txb_itr;

    } while (txb_itr < tuCount);

}

//****************************************
// ************ CuFullDistortionFastTuMode ****************
//****************************************/
void CuFullDistortionFastTuMode_R(
    LargestCodingUnit_t            *sb_ptr,
    ModeDecisionCandidateBuffer_t  *candidateBuffer,
    ModeDecisionContext_t            *context_ptr,
    ModeDecisionCandidate_t           *candidate_ptr,
    PictureControlSet_t            *picture_control_set_ptr,
    uint64_t                          cbFullDistortion[DIST_CALC_TOTAL],
    uint64_t                          crFullDistortion[DIST_CALC_TOTAL],
    uint32_t                          count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU],
    COMPONENT_TYPE                  component_type,
    uint64_t                         *cb_coeff_bits,
    uint64_t                         *cr_coeff_bits,
#if SPATIAL_SSE
    EbBool                           is_full_loop,
#endif
    EbAsm                            asm_type)
{
    (void)sb_ptr;

    uint64_t                          y_tu_coeff_bits;
    uint64_t                          cb_tu_coeff_bits;
    uint64_t                          cr_tu_coeff_bits;
    uint32_t                          tuOriginIndex;
    uint32_t                          txb_origin_x;
    uint32_t                          txb_origin_y;
    uint32_t                          currentTuIndex;
    int32_t                          chromaShift;
    uint32_t                          tuChromaOriginIndex;
    uint64_t                          tuFullDistortion[3][DIST_CALC_TOTAL];
    EbPictureBufferDesc_t          *transform_buffer;
    uint32_t                          tuTotalCount;
    uint32_t                          txb_itr = 0;
    //    SequenceControlSet_t           *sequence_control_set_ptr=((SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr);

    tuTotalCount = context_ptr->blk_geom->txb_count;
    currentTuIndex = 0;
    transform_buffer = context_ptr->trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr;


    uint32_t  txb_1d_offset = 0;
    candidate_ptr->u_has_coeff = 0;
    candidate_ptr->v_has_coeff = 0;

    do {

        txb_origin_x = context_ptr->blk_geom->tx_org_x[txb_itr];
        txb_origin_y = context_ptr->blk_geom->tx_org_y[txb_itr];



        tuOriginIndex = txb_origin_x + txb_origin_y * candidateBuffer->residualQuantCoeffPtr->stride_y;
        tuChromaOriginIndex = txb_1d_offset;
        // Reset the Bit Costs
        y_tu_coeff_bits = 0;
        cb_tu_coeff_bits = 0;
        cr_tu_coeff_bits = 0;

        if (component_type == COMPONENT_CHROMA_CB || component_type == COMPONENT_CHROMA_CR || component_type == COMPONENT_CHROMA || component_type == COMPONENT_ALL) {
            uint32_t countNonZeroCoeffsAll[3];
            countNonZeroCoeffsAll[0] = count_non_zero_coeffs[0][currentTuIndex];
            countNonZeroCoeffsAll[1] = count_non_zero_coeffs[1][currentTuIndex];
            countNonZeroCoeffsAll[2] = count_non_zero_coeffs[2][currentTuIndex];

#if SPATIAL_SSE
            if (is_full_loop && context_ptr->spatial_sse_full_loop) {
                EbPictureBufferDesc_t          *input_picture_ptr = picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr;
                uint32_t input_chroma_tu_origin_index = (((context_ptr->sb_origin_y + ((txb_origin_y >> 3) << 3)) >> 1) + (input_picture_ptr->origin_y >> 1)) * input_picture_ptr->strideCb + (((context_ptr->sb_origin_x + ((txb_origin_x >> 3) << 3)) >> 1) + (input_picture_ptr->origin_x >> 1));
                uint32_t tu_uv_origin_index = (((txb_origin_x >> 3) << 3) + (((txb_origin_y >> 3) << 3) * candidateBuffer->residualQuantCoeffPtr->strideCb)) >> 1;


                tuFullDistortion[1][DIST_CALC_PREDICTION] = spatial_full_distortion_kernel_func_ptr_array[asm_type][context_ptr->blk_geom->bwidth_uv == context_ptr->blk_geom->bheight_uv ? Log2f(context_ptr->blk_geom->bwidth_uv) - 2 : 6](
                    input_picture_ptr->bufferCb + input_chroma_tu_origin_index,
                    input_picture_ptr->strideCb,
                    candidateBuffer->prediction_ptr->bufferCb + tu_uv_origin_index,
                    candidateBuffer->prediction_ptr->strideCb,
                    context_ptr->blk_geom->tx_width_uv[txb_itr],
                    context_ptr->blk_geom->tx_height_uv[txb_itr]);

                tuFullDistortion[1][DIST_CALC_RESIDUAL] = spatial_full_distortion_kernel_func_ptr_array[asm_type][context_ptr->blk_geom->bwidth_uv == context_ptr->blk_geom->bheight_uv ? Log2f(context_ptr->blk_geom->bwidth_uv) - 2 : 6](
                    input_picture_ptr->bufferCb + input_chroma_tu_origin_index,
                    input_picture_ptr->strideCb,
                    &(((int8_t*)candidateBuffer->recon_ptr->bufferCb)[tu_uv_origin_index]),
                    candidateBuffer->recon_ptr->strideCb,
                    context_ptr->blk_geom->tx_width_uv[txb_itr],
                    context_ptr->blk_geom->tx_height_uv[txb_itr]);

                tuFullDistortion[2][DIST_CALC_PREDICTION] = spatial_full_distortion_kernel_func_ptr_array[asm_type][context_ptr->blk_geom->bwidth_uv == context_ptr->blk_geom->bheight_uv ? Log2f(context_ptr->blk_geom->bwidth_uv) - 2 : 6](
                    input_picture_ptr->bufferCr + input_chroma_tu_origin_index,
                    input_picture_ptr->strideCr,
                    candidateBuffer->prediction_ptr->bufferCr + tu_uv_origin_index,
                    candidateBuffer->prediction_ptr->strideCr,
                    context_ptr->blk_geom->tx_width_uv[txb_itr],
                    context_ptr->blk_geom->tx_height_uv[txb_itr]);

                tuFullDistortion[2][DIST_CALC_RESIDUAL] = spatial_full_distortion_kernel_func_ptr_array[asm_type][context_ptr->blk_geom->bwidth_uv == context_ptr->blk_geom->bheight_uv ? Log2f(context_ptr->blk_geom->bwidth_uv) - 2 : 6](
                    input_picture_ptr->bufferCr + input_chroma_tu_origin_index,
                    input_picture_ptr->strideCr,
                    &(((int8_t*)candidateBuffer->recon_ptr->bufferCr)[tu_uv_origin_index]),
                    candidateBuffer->recon_ptr->strideCr,
                    context_ptr->blk_geom->tx_width_uv[txb_itr],
                    context_ptr->blk_geom->tx_height_uv[txb_itr]);

                tuFullDistortion[1][DIST_CALC_PREDICTION]   <<= 4;
                tuFullDistortion[1][DIST_CALC_RESIDUAL]     <<= 4;
                tuFullDistortion[2][DIST_CALC_PREDICTION]   <<= 4;
                tuFullDistortion[2][DIST_CALC_RESIDUAL]     <<= 4;
            }
            else {
                // *Full Distortion (SSE)
                // *Note - there are known issues with how this distortion metric is currently
                //    calculated.  The amount of scaling between the two arrays is not
                //    equivalent.


                picture_full_distortion32_bits(
                    transform_buffer,
                    NOT_USED_VALUE,
                    tuChromaOriginIndex,
                    candidateBuffer->reconCoeffPtr,
                    NOT_USED_VALUE,
                    tuChromaOriginIndex,
                    NOT_USED_VALUE,
                    NOT_USED_VALUE,
                    context_ptr->blk_geom->tx_width_uv[txb_itr],
                    context_ptr->blk_geom->tx_height_uv[txb_itr],
                    tuFullDistortion[0],
                    tuFullDistortion[1],
                    tuFullDistortion[2],
                    countNonZeroCoeffsAll[0],
                    countNonZeroCoeffsAll[1],
                    countNonZeroCoeffsAll[2],
                    component_type,
                    asm_type);

                TxSize    txSize = context_ptr->blk_geom->txsize_uv[txb_itr];
                chromaShift = (MAX_TX_SCALE - av1_get_tx_scale(txSize)) * 2;
                tuFullDistortion[1][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(tuFullDistortion[1][DIST_CALC_RESIDUAL], chromaShift);
                tuFullDistortion[1][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(tuFullDistortion[1][DIST_CALC_PREDICTION], chromaShift);
                tuFullDistortion[2][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(tuFullDistortion[2][DIST_CALC_RESIDUAL], chromaShift);
                tuFullDistortion[2][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(tuFullDistortion[2][DIST_CALC_PREDICTION], chromaShift);


            }
#else
            // *Full Distortion (SSE)
            // *Note - there are known issues with how this distortion metric is currently
            //    calculated.  The amount of scaling between the two arrays is not
            //    equivalent.


            picture_full_distortion32_bits(
                transform_buffer,
                NOT_USED_VALUE,
                tuChromaOriginIndex,
                candidateBuffer->reconCoeffPtr,
                NOT_USED_VALUE,
                tuChromaOriginIndex,
                NOT_USED_VALUE,
                NOT_USED_VALUE,
                context_ptr->blk_geom->tx_width_uv[txb_itr],
                context_ptr->blk_geom->tx_height_uv[txb_itr],
                tuFullDistortion[0],
                tuFullDistortion[1],
                tuFullDistortion[2],
                countNonZeroCoeffsAll[0],
                countNonZeroCoeffsAll[1],
                countNonZeroCoeffsAll[2],
                component_type,
                asm_type);

            TxSize    txSize = context_ptr->blk_geom->txsize_uv[txb_itr];
            chromaShift = (MAX_TX_SCALE - av1_get_tx_scale(txSize)) * 2;
            tuFullDistortion[1][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(tuFullDistortion[1][DIST_CALC_RESIDUAL], chromaShift);
            tuFullDistortion[1][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(tuFullDistortion[1][DIST_CALC_PREDICTION], chromaShift);
            tuFullDistortion[2][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(tuFullDistortion[2][DIST_CALC_RESIDUAL], chromaShift);
            tuFullDistortion[2][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(tuFullDistortion[2][DIST_CALC_PREDICTION], chromaShift);


#endif


            //CHROMA-ONLY
            Av1TuEstimateCoeffBits(
#if CABAC_UP
                0,//allow_update_cdf,
                NULL,//FRAME_CONTEXT *ec_ctx,
#endif
                picture_control_set_ptr,
                candidateBuffer,
                context_ptr->cu_ptr,
                tuOriginIndex,
                tuChromaOriginIndex,
                context_ptr->coeff_est_entropy_coder_ptr,
                candidateBuffer->residualQuantCoeffPtr,
                count_non_zero_coeffs[0][currentTuIndex],
                count_non_zero_coeffs[1][currentTuIndex],
                count_non_zero_coeffs[2][currentTuIndex],
                &y_tu_coeff_bits,
                &cb_tu_coeff_bits,
                &cr_tu_coeff_bits,
                context_ptr->blk_geom->txsize[txb_itr],
                context_ptr->blk_geom->txsize_uv[txb_itr],

                component_type,
                asm_type);


            // OMK Useless ? We don't calculate Chroma CBF here
            av1_tu_calc_cost(
                candidate_ptr,
                context_ptr->cu_ptr->luma_txb_skip_context,
                currentTuIndex,
                count_non_zero_coeffs[0][currentTuIndex],
                count_non_zero_coeffs[1][currentTuIndex],
                count_non_zero_coeffs[2][currentTuIndex],
                tuFullDistortion[0],
                tuFullDistortion[1],
                tuFullDistortion[2],
                component_type,
                &y_tu_coeff_bits,
                &cb_tu_coeff_bits,
                &cr_tu_coeff_bits,
                context_ptr->blk_geom->txsize[txb_itr],
                context_ptr->full_lambda);



            *cb_coeff_bits += cb_tu_coeff_bits;
            *cr_coeff_bits += cr_tu_coeff_bits;
            cbFullDistortion[DIST_CALC_RESIDUAL] += tuFullDistortion[1][DIST_CALC_RESIDUAL];
            crFullDistortion[DIST_CALC_RESIDUAL] += tuFullDistortion[2][DIST_CALC_RESIDUAL];
            cbFullDistortion[DIST_CALC_PREDICTION] += tuFullDistortion[1][DIST_CALC_PREDICTION];
            crFullDistortion[DIST_CALC_PREDICTION] += tuFullDistortion[2][DIST_CALC_PREDICTION];

        }


        txb_1d_offset += context_ptr->blk_geom->tx_width_uv[txb_itr] * context_ptr->blk_geom->tx_height_uv[txb_itr];
        currentTuIndex++;


        ++txb_itr;


    } while (txb_itr < tuTotalCount);
}

#if IMPROVE_1D_INTER_DEPTH_DECISION
/***************************************
 * Check merge_block algorithm
 ***************************************/

EbBool merge_1D_inter_block(
    ModeDecisionContext_t    *context_ptr,
    uint32_t                  sq_idx,
    uint32_t                  nsq_idx) {
    EbBool merge_blocks = EB_FALSE;
    CodingUnit_t  *parent_cu_ptr = &context_ptr->md_cu_arr_nsq[sq_idx];
    CodingUnit_t  *child_cu_ptr = &context_ptr->md_cu_arr_nsq[nsq_idx];
    int parent_diriction = parent_cu_ptr->prediction_unit_array[0].inter_pred_direction_index;
    int parent_mv_l0 = parent_cu_ptr->prediction_unit_array[0].mv[REF_LIST_0].mvUnion;
    int parent_mv_l1 = parent_cu_ptr->prediction_unit_array[0].mv[REF_LIST_1].mvUnion;
    int parent_eob = parent_cu_ptr->block_has_coeff;
    int child_0_diriction = child_cu_ptr->prediction_unit_array[0].inter_pred_direction_index;
    int child_0_mv_l0 = child_cu_ptr->prediction_unit_array[0].mv[REF_LIST_0].mvUnion;
    int child_0_mv_l1 = child_cu_ptr->prediction_unit_array[0].mv[REF_LIST_1].mvUnion;
    int child_eob = child_cu_ptr->block_has_coeff;
    if (parent_diriction == child_0_diriction && child_eob == 0) {
        switch (parent_diriction) {
        case UNI_PRED_LIST_0:
            if (parent_mv_l0 == child_0_mv_l0) {
                merge_blocks = EB_TRUE;
            }
            break;
        case UNI_PRED_LIST_1:
            if (parent_mv_l1 == child_0_mv_l1) {
                merge_blocks = EB_TRUE;
            }
            break;
        case BI_PRED:
            if (parent_mv_l0 == child_0_mv_l0 &&
                parent_mv_l1 == child_0_mv_l1) {
                merge_blocks = EB_TRUE;
            }
            break;
        default:
            merge_blocks = EB_FALSE;
            break;
        }
    }
    return merge_blocks;
}
#endif
void  d1_non_square_block_decision(
    ModeDecisionContext_t               *context_ptr
)
{
    //compute total cost for the whole block partition
    uint64_t tot_cost = 0;
    uint32_t first_blk_idx = context_ptr->cu_ptr->mds_idx - (context_ptr->blk_geom->totns - 1);//index of first block in this partition
    uint32_t blk_it;
#if IMPROVE_1D_INTER_DEPTH_DECISION
    uint32_t merge_block_cnt = 0;
    EbBool merge_block_flag = EB_FALSE;
#endif
    for (blk_it = 0; blk_it < context_ptr->blk_geom->totns; blk_it++)
    {
        tot_cost += context_ptr->md_local_cu_unit[first_blk_idx + blk_it].cost;
#if IMPROVE_1D_INTER_DEPTH_DECISION
        merge_block_cnt += merge_1D_inter_block(context_ptr, context_ptr->blk_geom->sqi_mds, first_blk_idx + blk_it);
#endif
    }
#if SPLIT_RATE_FIX
    if (context_ptr->blk_geom->bsize > BLOCK_4X4) {
        uint64_t split_cost = 0;
        uint32_t parent_depth_idx_mds = context_ptr->blk_geom->sqi_mds;
        av1_split_flag_rate(
            context_ptr->sb_ptr->picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr,
            context_ptr,
            &context_ptr->md_cu_arr_nsq[parent_depth_idx_mds],
            0,
            from_shape_to_part[context_ptr->blk_geom->shape],
            &split_cost,
            context_ptr->full_lambda,
            context_ptr->md_rate_estimation_ptr,
            context_ptr->sb_ptr->picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr->max_sb_depth);

        tot_cost += split_cost;
    }
#endif

#if IMPROVE_1D_INTER_DEPTH_DECISION
    if (merge_block_cnt == context_ptr->blk_geom->totns) merge_block_flag = EB_TRUE;

    if (context_ptr->blk_geom->shape == PART_N || (tot_cost < context_ptr->md_local_cu_unit[context_ptr->blk_geom->sqi_mds].cost && merge_block_flag == EB_FALSE))
#else
    if (context_ptr->blk_geom->shape == PART_N || tot_cost < context_ptr->md_local_cu_unit[context_ptr->blk_geom->sqi_mds].cost)
#endif
    {
        //store best partition cost in parent square
        context_ptr->md_local_cu_unit[context_ptr->blk_geom->sqi_mds].cost = tot_cost;
        context_ptr->md_cu_arr_nsq[context_ptr->blk_geom->sqi_mds].part = from_shape_to_part[context_ptr->blk_geom->shape];
        context_ptr->md_cu_arr_nsq[context_ptr->blk_geom->sqi_mds].best_d1_blk = first_blk_idx;
    }


}


/// compute the cost of curr depth, and the depth above
void   compute_depth_costs(
    ModeDecisionContext_t    *context_ptr,
    SequenceControlSet_t     *sequence_control_set_ptr,
    uint32_t                  curr_depth_mds,
    uint32_t                  above_depth_mds,
    uint32_t                  step,
    uint64_t                 *above_depth_cost,
    uint64_t                 *curr_depth_cost)
{
    uint64_t       above_non_split_rate = 0;
    uint64_t       above_split_rate = 0;


    /*
    ___________
    |     |     |
    |blk0 |blk1 |
    |-----|-----|
    |blk2 |blk3 |
    |_____|_____|
    */
    // current depth blocks
    uint32_t       curr_depth_blk0_mds = curr_depth_mds - 3 * step;
    uint32_t       curr_depth_blk1_mds = curr_depth_mds - 2 * step;
    uint32_t       curr_depth_blk2_mds = curr_depth_mds - 1 * step;
    uint32_t       curr_depth_blk3_mds = curr_depth_mds;

    // Rate of not spliting the current depth (Depth != 4) in case the children were omitted by MDC
    uint64_t       curr_non_split_rate_blk0 = 0;
    uint64_t       curr_non_split_rate_blk1 = 0;
    uint64_t       curr_non_split_rate_blk2 = 0;
    uint64_t       curr_non_split_rate_blk3 = 0;


    context_ptr->md_local_cu_unit[above_depth_mds].left_neighbor_mode = context_ptr->md_local_cu_unit[curr_depth_blk0_mds].left_neighbor_mode;
    context_ptr->md_local_cu_unit[above_depth_mds].left_neighbor_depth = context_ptr->md_local_cu_unit[curr_depth_blk0_mds].left_neighbor_depth;
    context_ptr->md_local_cu_unit[above_depth_mds].top_neighbor_mode = context_ptr->md_local_cu_unit[curr_depth_blk0_mds].top_neighbor_mode;
    context_ptr->md_local_cu_unit[above_depth_mds].top_neighbor_depth = context_ptr->md_local_cu_unit[curr_depth_blk0_mds].top_neighbor_depth;
    context_ptr->md_local_cu_unit[above_depth_mds].left_neighbor_partition = context_ptr->md_local_cu_unit[curr_depth_blk0_mds].left_neighbor_partition;
    context_ptr->md_local_cu_unit[above_depth_mds].above_neighbor_partition = context_ptr->md_local_cu_unit[curr_depth_blk0_mds].above_neighbor_partition;

    // Compute above depth  cost
    if (context_ptr->md_local_cu_unit[above_depth_mds].tested_cu_flag == EB_TRUE)
    {
#if !SPLIT_RATE_FIX
        av1_split_flag_rate(
            sequence_control_set_ptr,
            context_ptr,
            &context_ptr->md_cu_arr_nsq[above_depth_mds],
            0,
            PARTITION_NONE,//shouldn't this be final partition for above depth?
            &above_non_split_rate,
            context_ptr->full_lambda,
            context_ptr->md_rate_estimation_ptr,
            sequence_control_set_ptr->max_sb_depth);
#endif
        *above_depth_cost = context_ptr->md_local_cu_unit[above_depth_mds].cost + above_non_split_rate;
    }
    else {
        *above_depth_cost = MAX_MODE_COST;
    }

    // Compute curr depth  cost
    av1_split_flag_rate(
        sequence_control_set_ptr,
        context_ptr,
        &context_ptr->md_cu_arr_nsq[above_depth_mds],
        0,
        PARTITION_SPLIT,
        &above_split_rate,
        context_ptr->full_lambda,
        context_ptr->md_rate_estimation_ptr,
        sequence_control_set_ptr->max_sb_depth);

    if (context_ptr->blk_geom->bsize > BLOCK_4X4) {

        if (context_ptr->md_cu_arr_nsq[curr_depth_blk0_mds].mdc_split_flag == 0)
            av1_split_flag_rate(
                sequence_control_set_ptr,
                context_ptr,
                &context_ptr->md_cu_arr_nsq[curr_depth_blk0_mds],
                0,
                PARTITION_NONE,
                &curr_non_split_rate_blk0,
                context_ptr->full_lambda,
                context_ptr->md_rate_estimation_ptr,
                sequence_control_set_ptr->max_sb_depth);

        if (context_ptr->md_cu_arr_nsq[curr_depth_blk1_mds].mdc_split_flag == 0)
            av1_split_flag_rate(
                sequence_control_set_ptr,
                context_ptr,
                &context_ptr->md_cu_arr_nsq[curr_depth_blk1_mds],
                0,
                PARTITION_NONE,
                &curr_non_split_rate_blk1,
                context_ptr->full_lambda,
                context_ptr->md_rate_estimation_ptr,
                sequence_control_set_ptr->max_sb_depth);

        if (context_ptr->md_cu_arr_nsq[curr_depth_blk2_mds].mdc_split_flag == 0)
            av1_split_flag_rate(
                sequence_control_set_ptr,
                context_ptr,
                &context_ptr->md_cu_arr_nsq[curr_depth_blk2_mds],
                0,
                PARTITION_NONE,
                &curr_non_split_rate_blk2,
                context_ptr->full_lambda,
                context_ptr->md_rate_estimation_ptr,
                sequence_control_set_ptr->max_sb_depth);

        if (context_ptr->md_cu_arr_nsq[curr_depth_blk3_mds].mdc_split_flag == 0)
            av1_split_flag_rate(
                sequence_control_set_ptr,
                context_ptr,
                &context_ptr->md_cu_arr_nsq[curr_depth_blk3_mds],
                0,
                PARTITION_NONE,
                &curr_non_split_rate_blk3,
                context_ptr->full_lambda,
                context_ptr->md_rate_estimation_ptr,
                sequence_control_set_ptr->max_sb_depth);
    }
    //curr_non_split_rate_344 = splitflag_mdc_344 || 4x4 ? 0 : compute; 


    *curr_depth_cost =
        context_ptr->md_local_cu_unit[curr_depth_mds].cost + curr_non_split_rate_blk3 +
        context_ptr->md_local_cu_unit[curr_depth_mds - 1 * step].cost + curr_non_split_rate_blk2 +
        context_ptr->md_local_cu_unit[curr_depth_mds - 2 * step].cost + curr_non_split_rate_blk1 +
        context_ptr->md_local_cu_unit[curr_depth_mds - 3 * step].cost + curr_non_split_rate_blk0 +
        above_split_rate;

}

uint32_t d2_inter_depth_block_decision(
    ModeDecisionContext_t          *context_ptr,
    uint32_t                        blk_mds,
    LargestCodingUnit_t            *tbPtr,
    uint32_t                          lcuAddr,
    uint32_t                          tbOriginX,
    uint32_t                          tbOriginY,
    uint64_t                          full_lambda,
    MdRateEstimationContext_t      *md_rate_estimation_ptr,
    PictureControlSet_t            *picture_control_set_ptr)
{
    UNUSED(tbPtr);
    UNUSED(lcuAddr);
    UNUSED(tbOriginX);
    UNUSED(tbOriginY);
    UNUSED(full_lambda);
    UNUSED(md_rate_estimation_ptr);

    uint32_t                  lastCuIndex, d0_idx_mds, d1_idx_mds, d2_idx_mds, top_left_idx_mds;
    UNUSED(top_left_idx_mds);
    UNUSED(d2_idx_mds);
    UNUSED(d1_idx_mds);
    UNUSED(d0_idx_mds);
    uint64_t                    parent_depth_cost = 0, current_depth_cost = 0;
    SequenceControlSet_t     *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    EbBool                    lastDepthFlag;
    const BlockGeom          * blk_geom;

    lastDepthFlag = context_ptr->md_cu_arr_nsq[blk_mds].split_flag == EB_FALSE ? EB_TRUE : EB_FALSE;
    d1_idx_mds = blk_mds;
    d2_idx_mds = blk_mds;
    lastCuIndex = blk_mds;
    blk_geom = get_blk_geom_mds(blk_mds);
    uint32_t    parent_depth_idx_mds = blk_mds;
    uint32_t    current_depth_idx_mds = blk_mds;

    if (lastDepthFlag) {
        while (blk_geom->is_last_quadrant) {

            //get parent idx
            parent_depth_idx_mds = current_depth_idx_mds - parent_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][blk_geom->depth];
            if (picture_control_set_ptr->slice_type == I_SLICE && parent_depth_idx_mds == 0) {
                parent_depth_cost = MAX_MODE_COST;
            }
            else {

                compute_depth_costs(context_ptr, sequence_control_set_ptr, current_depth_idx_mds, parent_depth_idx_mds, ns_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][blk_geom->depth], &parent_depth_cost, &current_depth_cost);
            }
            if (parent_depth_cost <= current_depth_cost) {
                context_ptr->md_cu_arr_nsq[parent_depth_idx_mds].split_flag = EB_FALSE;
                context_ptr->md_local_cu_unit[parent_depth_idx_mds].cost = parent_depth_cost;
                lastCuIndex = parent_depth_idx_mds;
            }
            else {
                context_ptr->md_local_cu_unit[parent_depth_idx_mds].cost = current_depth_cost;
                context_ptr->md_cu_arr_nsq[parent_depth_idx_mds].part = PARTITION_SPLIT;
            }

            //setup next parent inter depth
            blk_geom = get_blk_geom_mds(parent_depth_idx_mds);
            current_depth_idx_mds = parent_depth_idx_mds;
        }
    }


    return lastCuIndex;
}



