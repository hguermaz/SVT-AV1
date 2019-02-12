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
#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbRateControlProcess.h"
#include "EbSystemResourceManager.h"
#include "EbSequenceControlSet.h"
#include "EbPictureControlSet.h"
#include "EbUtility.h"
#include "EbErrorCodes.h"

#include "EbRateControlResults.h"
#include "EbRateControlTasks.h"
#include "RateControlModel.h"


static uint8_t QP_OFFSET_LAYER_ARRAY[MAX_TEMPORAL_LAYERS] =
{
    1, 2, 4, 5, 6, 7
};

/*****************************
* Internal Typedefs
*****************************/
void RateControlLayerReset(
    RateControlLayerContext_t   *rateControlLayerPtr,
    PictureControlSet_t         *picture_control_set_ptr,
    RateControlContext_t        *rateControlContextPtr,
    uint32_t                       pictureAreaInPixel,
    EbBool                      wasUsed)
{

    SequenceControlSet_t *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
    uint32_t                sliceNum;
    uint32_t                temporal_layer_index;
    uint64_t                totalFrameInInterval;
    uint64_t                sumBitsPerSw = 0;

    rateControlLayerPtr->target_bit_rate = picture_control_set_ptr->parent_pcs_ptr->target_bit_rate*RATE_PERCENTAGE_LAYER_ARRAY[sequence_control_set_ptr->static_config.hierarchical_levels][rateControlLayerPtr->temporalIndex] / 100;
    // update this based on temporal layers
    rateControlLayerPtr->frame_rate = sequence_control_set_ptr->frame_rate;

    totalFrameInInterval = sequence_control_set_ptr->static_config.intra_period_length + 1;

    if (sequence_control_set_ptr->static_config.look_ahead_distance != 0 && sequence_control_set_ptr->intra_period_length != -1) {
        if (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0) {
            totalFrameInInterval = 0;
            for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++) {
                rateControlContextPtr->frames_in_interval[temporal_layer_index] = picture_control_set_ptr->parent_pcs_ptr->frames_in_interval[temporal_layer_index];
                totalFrameInInterval += picture_control_set_ptr->parent_pcs_ptr->frames_in_interval[temporal_layer_index];
                sumBitsPerSw += picture_control_set_ptr->parent_pcs_ptr->bits_per_sw_per_layer[temporal_layer_index];
            }
#if ADAPTIVE_PERCENTAGE
            rateControlLayerPtr->target_bit_rate = picture_control_set_ptr->parent_pcs_ptr->target_bit_rate* picture_control_set_ptr->parent_pcs_ptr->bits_per_sw_per_layer[rateControlLayerPtr->temporalIndex] / sumBitsPerSw;
#endif
        }
    }


    if (rateControlLayerPtr->temporalIndex == 0) {
        rateControlLayerPtr->coeffAveragingWeight1 = 5;
        rateControlLayerPtr->frame_rate = rateControlLayerPtr->frame_rate >> 5;
    }
    else if (rateControlLayerPtr->temporalIndex == 1) {
        rateControlLayerPtr->coeffAveragingWeight1 = 5;
        rateControlLayerPtr->frame_rate = rateControlLayerPtr->frame_rate >> 5;
    }
    else if (rateControlLayerPtr->temporalIndex == 2) {
        rateControlLayerPtr->coeffAveragingWeight1 = 5;
        rateControlLayerPtr->frame_rate = rateControlLayerPtr->frame_rate >> 4;
    }
    else if (rateControlLayerPtr->temporalIndex == 3) {
        rateControlLayerPtr->coeffAveragingWeight1 = 5;
        rateControlLayerPtr->frame_rate = rateControlLayerPtr->frame_rate >> 3;
    }
    else if (rateControlLayerPtr->temporalIndex == 4) {
        rateControlLayerPtr->coeffAveragingWeight1 = 5;
        rateControlLayerPtr->frame_rate = rateControlLayerPtr->frame_rate >> 2;
    }
    else if (rateControlLayerPtr->temporalIndex == 5) {
        rateControlLayerPtr->coeffAveragingWeight1 = 3;
        rateControlLayerPtr->frame_rate = rateControlLayerPtr->frame_rate >> 1;
    }
    if (sequence_control_set_ptr->static_config.intra_period_length != -1) {
        rateControlLayerPtr->frame_rate = sequence_control_set_ptr->frame_rate * rateControlContextPtr->frames_in_interval[rateControlLayerPtr->temporalIndex] / totalFrameInInterval;
    }

    rateControlLayerPtr->coeffAveragingWeight2 = 16 - rateControlLayerPtr->coeffAveragingWeight1;
    if (rateControlLayerPtr->frame_rate == 0) { // no frame in that layer
        rateControlLayerPtr->frame_rate = 1 << RC_PRECISION;
    }
    rateControlLayerPtr->channelBitRate = (((rateControlLayerPtr->target_bit_rate << (2 * RC_PRECISION)) / rateControlLayerPtr->frame_rate) + RC_PRECISION_OFFSET) >> RC_PRECISION;
    rateControlLayerPtr->channelBitRate = (uint64_t)MAX((int64_t)1, (int64_t)rateControlLayerPtr->channelBitRate);
    rateControlLayerPtr->ecBitConstraint = rateControlLayerPtr->channelBitRate;


    // This is only for the initial frame, because the feedback is from packetization now and all of these are considered
    // considering the bits for slice header
    // *Note - only one-slice-per picture is supported for UHD
    sliceNum = 1;

    rateControlLayerPtr->ecBitConstraint -= SLICE_HEADER_BITS_NUM * sliceNum;

    rateControlLayerPtr->ecBitConstraint = MAX(1, rateControlLayerPtr->ecBitConstraint);

    rateControlLayerPtr->previousBitConstraint = rateControlLayerPtr->channelBitRate;
    rateControlLayerPtr->bitConstraint = rateControlLayerPtr->channelBitRate;
    rateControlLayerPtr->difTotalAndEcBits = 0;

    rateControlLayerPtr->frameSameSADMinQpCount = 0;
    rateControlLayerPtr->maxQp = picture_control_set_ptr->picture_qp;

    rateControlLayerPtr->alpha = 1 << (RC_PRECISION - 1);
    {
        if (!wasUsed) {


            rateControlLayerPtr->sameSADCount = 0;

            rateControlLayerPtr->kCoeff = 3 << RC_PRECISION;
            rateControlLayerPtr->previousKCoeff = 3 << RC_PRECISION;

            rateControlLayerPtr->cCoeff = (rateControlLayerPtr->channelBitRate << (2 * RC_PRECISION)) / pictureAreaInPixel / CCOEFF_INIT_FACT;
            rateControlLayerPtr->previousCCoeff = (rateControlLayerPtr->channelBitRate << (2 * RC_PRECISION)) / pictureAreaInPixel / CCOEFF_INIT_FACT;

            // These are for handling Pred structure 2, when for higher temporal layer, frames can arrive in different orders
            // They should be modifed in a way that gets these from previous layers
            rateControlLayerPtr->previousFrameQp = 32;
            rateControlLayerPtr->previousFrameBitActual = 1200;
            rateControlLayerPtr->previousFrameQuantizedCoeffBitActual = 1000;
            rateControlLayerPtr->previousFrameSadMe = 10000000;
            rateControlLayerPtr->previousFrameQp = picture_control_set_ptr->picture_qp;
            rateControlLayerPtr->deltaQpFraction = 0;
            rateControlLayerPtr->previousFrameAverageQp = picture_control_set_ptr->picture_qp;
            rateControlLayerPtr->previousCalculatedFrameQp = picture_control_set_ptr->picture_qp;
            rateControlLayerPtr->calculatedFrameQp = picture_control_set_ptr->picture_qp;
            rateControlLayerPtr->criticalStates = 0;
        }
        else {
            rateControlLayerPtr->sameSADCount = 0;
            rateControlLayerPtr->criticalStates = 0;
        }
    }
}


void RateControlLayerResetPart2(
    RateControlLayerContext_t   *rateControlLayerPtr,
    PictureControlSet_t         *picture_control_set_ptr)
{

    // update this based on temporal layers

    rateControlLayerPtr->maxQp = (uint32_t)CLIP3(0, 63, (int32_t)(picture_control_set_ptr->picture_qp + QP_OFFSET_LAYER_ARRAY[rateControlLayerPtr->temporalIndex]));;
    {

        // These are for handling Pred structure 2, when for higher temporal layer, frames can arrive in different orders
        // They should be modifed in a way that gets these from previous layers
        rateControlLayerPtr->previousFrameQp = rateControlLayerPtr->maxQp;
        rateControlLayerPtr->previousFrameAverageQp = rateControlLayerPtr->maxQp;
        rateControlLayerPtr->previousCalculatedFrameQp = rateControlLayerPtr->maxQp;
        rateControlLayerPtr->calculatedFrameQp = rateControlLayerPtr->maxQp;
    }
}

EbErrorType HighLevelRateControlContextCtor(
    HighLevelRateControlContext_t   **entryDblPtr) {

    HighLevelRateControlContext_t *entryPtr;
    EB_MALLOC(HighLevelRateControlContext_t*, entryPtr, sizeof(HighLevelRateControlContext_t), EB_N_PTR);
    *entryDblPtr = entryPtr;

    return EB_ErrorNone;
}


EbErrorType RateControlLayerContextCtor(
    RateControlLayerContext_t   **entryDblPtr) {

    RateControlLayerContext_t *entryPtr;
    EB_MALLOC(RateControlLayerContext_t*, entryPtr, sizeof(RateControlLayerContext_t), EB_N_PTR);

    *entryDblPtr = entryPtr;

    entryPtr->firstFrame = 1;
    entryPtr->firstNonIntraFrame = 1;
    entryPtr->feedbackArrived = EB_FALSE;

    return EB_ErrorNone;
}



EbErrorType RateControlIntervalParamContextCtor(
    RateControlIntervalParamContext_t   **entryDblPtr) {

    uint32_t temporalIndex;
    EbErrorType return_error = EB_ErrorNone;
    RateControlIntervalParamContext_t *entryPtr;
    EB_MALLOC(RateControlIntervalParamContext_t*, entryPtr, sizeof(RateControlIntervalParamContext_t), EB_N_PTR);

    *entryDblPtr = entryPtr;

    entryPtr->inUse = EB_FALSE;
    entryPtr->wasUsed = EB_FALSE;
    entryPtr->lastGop = EB_FALSE;
    entryPtr->processedFramesNumber = 0;
    EB_MALLOC(RateControlLayerContext_t**, entryPtr->rateControlLayerArray, sizeof(RateControlLayerContext_t*)*EB_MAX_TEMPORAL_LAYERS, EB_N_PTR);

    for (temporalIndex = 0; temporalIndex < EB_MAX_TEMPORAL_LAYERS; temporalIndex++) {
        return_error = RateControlLayerContextCtor(&entryPtr->rateControlLayerArray[temporalIndex]);
        entryPtr->rateControlLayerArray[temporalIndex]->temporalIndex = temporalIndex;
        entryPtr->rateControlLayerArray[temporalIndex]->frame_rate = 1 << RC_PRECISION;
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    entryPtr->min_target_rate_assigned = EB_FALSE;

    entryPtr->intraFramesQp = 0;
    entryPtr->nextGopIntraFrameQp = 0;
    entryPtr->firstPicPredBits = 0;
    entryPtr->firstPicActualBits = 0;
    entryPtr->firstPicPredQp = 0;
    entryPtr->firstPicActualQp = 0;
    entryPtr->firstPicActualQpAssigned = EB_FALSE;
    entryPtr->scene_change_in_gop = EB_FALSE;
    entryPtr->extraApBitRatioI = 0;

    return EB_ErrorNone;
}

EbErrorType RateControlCodedFramesStatsContextCtor(
    CodedFramesStatsEntry_t   **entryDblPtr,
    uint64_t                      picture_number) {

    CodedFramesStatsEntry_t *entryPtr;
    EB_MALLOC(CodedFramesStatsEntry_t*, entryPtr, sizeof(CodedFramesStatsEntry_t), EB_N_PTR);

    *entryDblPtr = entryPtr;

    entryPtr->picture_number = picture_number;
    entryPtr->frameTotalBitActual = -1;

    return EB_ErrorNone;
}


EbErrorType RateControlContextCtor(
    RateControlContext_t   **context_dbl_ptr,
    EbFifo_t                *rateControlInputTasksFifoPtr,
    EbFifo_t                *rateControlOutputResultsFifoPtr,
    int32_t                   intra_period_length)
{
    uint32_t temporalIndex;
    uint32_t intervalIndex;

#if OVERSHOOT_STAT_PRINT
    uint32_t pictureIndex;
#endif

    EbErrorType return_error = EB_ErrorNone;
    RateControlContext_t *context_ptr;
    EB_MALLOC(RateControlContext_t*, context_ptr, sizeof(RateControlContext_t), EB_N_PTR);

    *context_dbl_ptr = context_ptr;

    context_ptr->rateControlInputTasksFifoPtr = rateControlInputTasksFifoPtr;
    context_ptr->rateControlOutputResultsFifoPtr = rateControlOutputResultsFifoPtr;

    // High level RC
    return_error = HighLevelRateControlContextCtor(
        &context_ptr->highLevelRateControlPtr);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    for (temporalIndex = 0; temporalIndex < EB_MAX_TEMPORAL_LAYERS; temporalIndex++) {
        context_ptr->frames_in_interval[temporalIndex] = 0;
    }

    EB_MALLOC(RateControlIntervalParamContext_t**, context_ptr->rateControlParamQueue, sizeof(RateControlIntervalParamContext_t*)*PARALLEL_GOP_MAX_NUMBER, EB_N_PTR);

    context_ptr->rateControlParamQueueHeadIndex = 0;
    for (intervalIndex = 0; intervalIndex < PARALLEL_GOP_MAX_NUMBER; intervalIndex++) {
        return_error = RateControlIntervalParamContextCtor(
            &context_ptr->rateControlParamQueue[intervalIndex]);
        context_ptr->rateControlParamQueue[intervalIndex]->firstPoc = (intervalIndex*(uint32_t)(intra_period_length + 1));
        context_ptr->rateControlParamQueue[intervalIndex]->lastPoc = ((intervalIndex + 1)*(uint32_t)(intra_period_length + 1)) - 1;
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

#if OVERSHOOT_STAT_PRINT
    context_ptr->codedFramesStatQueueHeadIndex = 0;
    context_ptr->codedFramesStatQueueTailIndex = 0;
    EB_MALLOC(CodedFramesStatsEntry_t**, context_ptr->codedFramesStatQueue, sizeof(CodedFramesStatsEntry_t*)*CODED_FRAMES_STAT_QUEUE_MAX_DEPTH, EB_N_PTR);

    for (pictureIndex = 0; pictureIndex < CODED_FRAMES_STAT_QUEUE_MAX_DEPTH; ++pictureIndex) {
        return_error = RateControlCodedFramesStatsContextCtor(
            &context_ptr->codedFramesStatQueue[pictureIndex],
            pictureIndex);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }
    context_ptr->maxBitActualPerSw = 0;
    context_ptr->maxBitActualPerGop = 0;
#endif

    context_ptr->baseLayerFramesAvgQp = 0;
    context_ptr->baseLayerIntraFramesAvgQp = 0;


    context_ptr->intraCoefRate = 4;
    context_ptr->extraBits = 0;
    context_ptr->extraBitsGen = 0;
    context_ptr->maxRateAdjustDeltaQP = 0;


    return EB_ErrorNone;
}
void HighLevelRcInputPictureMode2(
    PictureParentControlSet_t         *picture_control_set_ptr,
    SequenceControlSet_t              *sequence_control_set_ptr,
    EncodeContext_t                   *encode_context_ptr,
    RateControlContext_t              *context_ptr,
    HighLevelRateControlContext_t     *highLevelRateControlPtr)
{

    EbBool                             end_of_sequence_flag = EB_TRUE;

    HlRateControlHistogramEntry_t      *hlRateControlHistogramPtrTemp;
    // Queue variables
    uint32_t                             queueEntryIndexTemp;
    uint32_t                             queueEntryIndexTemp2;
    uint32_t                             queueEntryIndexHeadTemp;


    uint64_t                              minLaBitDistance;
    uint32_t                              selectedRefQpTableIndex;
    uint32_t                              selectedRefQp;
#if RC_UPDATE_TARGET_RATE
    uint32_t                              selectedOrgRefQp;
#endif
    uint32_t                                previous_selected_ref_qp = encode_context_ptr->previous_selected_ref_qp;
    uint64_t                                max_coded_poc = encode_context_ptr->max_coded_poc;
    uint32_t                                max_coded_poc_selected_ref_qp = encode_context_ptr->max_coded_poc_selected_ref_qp;


    uint32_t                              refQpIndex;
    uint32_t                              refQpIndexTemp;
    uint32_t                              refQpTableIndex;

    uint32_t                              areaInPixel;
    uint32_t                              numOfFullLcus;
    uint32_t                              qpSearchMin;
    uint32_t                              qpSearchMax;
    int32_t                              qpStep = 1;
    EbBool                             bestQpFound;
    uint32_t                              temporal_layer_index;
    EbBool                             tables_updated;

    uint64_t                              bitConstraintPerSw = 0;

    RateControlTables_t                    *rateControlTablesPtr;
    EB_Bit_Number                        *sadBitsArrayPtr;
    EB_Bit_Number                        *intraSadBitsArrayPtr;
    uint32_t                               pred_bits_ref_qp;

    for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++) {
        picture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_index] = 0;
    }
    picture_control_set_ptr->total_bits_per_gop = 0;

    areaInPixel = sequence_control_set_ptr->luma_width * sequence_control_set_ptr->luma_height;;

    EbBlockOnMutex(sequence_control_set_ptr->encode_context_ptr->rate_table_update_mutex);

    tables_updated = sequence_control_set_ptr->encode_context_ptr->rate_control_tables_array_updated;
    picture_control_set_ptr->percentage_updated = EB_FALSE;

    if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {

        // Increamenting the head of the hl_rate_control_historgram_queue and clean up the entores
        hlRateControlHistogramPtrTemp = (encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]);
        while ((hlRateControlHistogramPtrTemp->lifeCount == 0) && hlRateControlHistogramPtrTemp->passedToHlrc) {

            EbBlockOnMutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
            // Reset the Reorder Queue Entry
            hlRateControlHistogramPtrTemp->picture_number += INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH;
            hlRateControlHistogramPtrTemp->lifeCount = -1;
            hlRateControlHistogramPtrTemp->passedToHlrc = EB_FALSE;
            hlRateControlHistogramPtrTemp->isCoded = EB_FALSE;
            hlRateControlHistogramPtrTemp->totalNumBitsCoded = 0;

            // Increment the Reorder Queue head Ptr
            encode_context_ptr->hl_rate_control_historgram_queue_head_index =
                (encode_context_ptr->hl_rate_control_historgram_queue_head_index == HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->hl_rate_control_historgram_queue_head_index + 1;
            EbReleaseMutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
            hlRateControlHistogramPtrTemp = encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index];

        }
        // For the case that number of frames in the sliding window is less than size of the look ahead or intra Refresh. i.e. end of sequence
        if ((picture_control_set_ptr->frames_in_sw < MIN(sequence_control_set_ptr->static_config.look_ahead_distance + 1, (uint32_t)sequence_control_set_ptr->intra_period_length + 1))) {

            selectedRefQp = max_coded_poc_selected_ref_qp;

            // Update the QP for the sliding window based on the status of RC
            if ((context_ptr->extraBitsGen > (int64_t)(context_ptr->virtualBufferSize << 3))) {
                selectedRefQp = (uint32_t)MAX((int32_t)selectedRefQp - 2, 0);
            }
            else if ((context_ptr->extraBitsGen > (int64_t)(context_ptr->virtualBufferSize << 2))) {
                selectedRefQp = (uint32_t)MAX((int32_t)selectedRefQp - 1, 0);
            }
            if ((context_ptr->extraBitsGen < -(int64_t)(context_ptr->virtualBufferSize << 2))) {
                selectedRefQp += 2;
            }
            else if ((context_ptr->extraBitsGen < -(int64_t)(context_ptr->virtualBufferSize << 1))) {
                selectedRefQp += 1;
            }

            if ((picture_control_set_ptr->frames_in_sw < (uint32_t)(sequence_control_set_ptr->intra_period_length + 1)) &&
                (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0)) {
                selectedRefQp = (uint32_t)CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    selectedRefQp + 1);
            }

            queueEntryIndexHeadTemp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
            queueEntryIndexHeadTemp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
            queueEntryIndexHeadTemp = (queueEntryIndexHeadTemp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                queueEntryIndexHeadTemp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                queueEntryIndexHeadTemp;

            queueEntryIndexTemp = queueEntryIndexHeadTemp;
            {

                hlRateControlHistogramPtrTemp = (encode_context_ptr->hl_rate_control_historgram_queue[queueEntryIndexTemp]);
                refQpIndexTemp = selectedRefQp + QP_OFFSET_LAYER_ARRAY[hlRateControlHistogramPtrTemp->temporal_layer_index];
                refQpIndexTemp = (uint32_t)CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    refQpIndexTemp);

                if (hlRateControlHistogramPtrTemp->slice_type == I_SLICE) {
                    refQpIndexTemp = (uint32_t)MAX((int32_t)refQpIndexTemp + RC_INTRA_QP_OFFSET, 0);
                }

                hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = 0;
                rateControlTablesPtr = &encode_context_ptr->rate_control_tables_array[refQpIndexTemp];
                sadBitsArrayPtr = rateControlTablesPtr->sadBitsArray[hlRateControlHistogramPtrTemp->temporal_layer_index];
                intraSadBitsArrayPtr = rateControlTablesPtr->intraSadBitsArray[hlRateControlHistogramPtrTemp->temporal_layer_index];
                pred_bits_ref_qp = 0;
                numOfFullLcus = 0;

                if (hlRateControlHistogramPtrTemp->slice_type == I_SLICE) {
                    // Loop over block in the frame and calculated the predicted bits at reg QP
                    {
                        unsigned i;
                        uint32_t accum = 0;
                        for (i = 0; i < NUMBER_OF_INTRA_SAD_INTERVALS; ++i)
                        {
                            accum += (uint32_t)(hlRateControlHistogramPtrTemp->ois_distortion_histogram[i] * intraSadBitsArrayPtr[i]);
                        }

                        pred_bits_ref_qp = accum;
                        numOfFullLcus = hlRateControlHistogramPtrTemp->full_sb_count;
                    }
                    hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] += pred_bits_ref_qp;
                }

                else {
                    {
                        unsigned i;
                        uint32_t accum = 0;
                        for (i = 0; i < NUMBER_OF_SAD_INTERVALS; ++i)
                        {
                            accum += (uint32_t)(hlRateControlHistogramPtrTemp->me_distortion_histogram[i] * sadBitsArrayPtr[i]);
                        }

                        pred_bits_ref_qp = accum;
                        numOfFullLcus = hlRateControlHistogramPtrTemp->full_sb_count;

                    }
                    hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] += pred_bits_ref_qp;
                }

                // Scale for in complete
                //  pred_bits_ref_qp is normalized based on the area because of the LCUs at the picture boundries
                hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] * (uint64_t)areaInPixel / (numOfFullLcus << 12);

                // Store the pred_bits_ref_qp for the first frame in the window to PCS
                picture_control_set_ptr->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];

            }
        }
        else {
            // Loop over the QPs and find the best QP
            minLaBitDistance = MAX_UNSIGNED_VALUE;
            qpSearchMin = (uint8_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                sequence_control_set_ptr->static_config.max_qp_allowed,
                (uint32_t)MAX((int32_t)sequence_control_set_ptr->qp - 20, 0));

            qpSearchMax = (uint8_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                sequence_control_set_ptr->static_config.max_qp_allowed,
                sequence_control_set_ptr->qp + 20);

            for (refQpTableIndex = qpSearchMin; refQpTableIndex < qpSearchMax; refQpTableIndex++) {
                highLevelRateControlPtr->predBitsRefQpPerSw[refQpTableIndex] = 0;
            }

            bitConstraintPerSw = highLevelRateControlPtr->bitConstraintPerSw * picture_control_set_ptr->frames_in_sw / (sequence_control_set_ptr->static_config.look_ahead_distance + 1);

            // Update the target rate for the sliding window based on the status of RC
            if ((context_ptr->extraBitsGen > (int64_t)(context_ptr->virtualBufferSize * 10))) {
                bitConstraintPerSw = bitConstraintPerSw * 130 / 100;
            }
            else if ((context_ptr->extraBitsGen > (int64_t)(context_ptr->virtualBufferSize << 3))) {
                bitConstraintPerSw = bitConstraintPerSw * 120 / 100;
            }
            else if ((context_ptr->extraBitsGen > (int64_t)(context_ptr->virtualBufferSize << 2))) {
                bitConstraintPerSw = bitConstraintPerSw * 110 / 100;
            }
            if ((context_ptr->extraBitsGen < -(int64_t)(context_ptr->virtualBufferSize << 3))) {
                bitConstraintPerSw = bitConstraintPerSw * 80 / 100;
            }
            else if ((context_ptr->extraBitsGen < -(int64_t)(context_ptr->virtualBufferSize << 2))) {
                bitConstraintPerSw = bitConstraintPerSw * 90 / 100;
            }

            // Loop over proper QPs and find the Predicted bits for that QP. Find the QP with the closest total predicted rate to target bits for the sliding window.
            previous_selected_ref_qp = CLIP3(
                qpSearchMin,
                qpSearchMax,
                previous_selected_ref_qp);
            refQpTableIndex = previous_selected_ref_qp;
            selectedRefQpTableIndex = refQpTableIndex;
            selectedRefQp = refQpListTable[selectedRefQpTableIndex];
            bestQpFound = EB_FALSE;
            while (refQpTableIndex >= qpSearchMin && refQpTableIndex <= qpSearchMax && !bestQpFound) {

                refQpIndex = CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    refQpListTable[refQpTableIndex]);
                highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] = 0;

                // Finding the predicted bits for each frame in the sliding window at the reference Qp(s)
                queueEntryIndexHeadTemp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
                queueEntryIndexHeadTemp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                queueEntryIndexHeadTemp = (queueEntryIndexHeadTemp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                    queueEntryIndexHeadTemp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                    queueEntryIndexHeadTemp;

                queueEntryIndexTemp = queueEntryIndexHeadTemp;
                // This is set to false, so the last frame would go inside the loop
                end_of_sequence_flag = EB_FALSE;

                while (!end_of_sequence_flag &&
                    queueEntryIndexTemp <= queueEntryIndexHeadTemp + sequence_control_set_ptr->static_config.look_ahead_distance) {

                    queueEntryIndexTemp2 = (queueEntryIndexTemp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? queueEntryIndexTemp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH : queueEntryIndexTemp;
                    hlRateControlHistogramPtrTemp = (encode_context_ptr->hl_rate_control_historgram_queue[queueEntryIndexTemp2]);

                    refQpIndexTemp = refQpIndex + QP_OFFSET_LAYER_ARRAY[hlRateControlHistogramPtrTemp->temporal_layer_index];
                    refQpIndexTemp = (uint32_t)CLIP3(
                        sequence_control_set_ptr->static_config.min_qp_allowed,
                        sequence_control_set_ptr->static_config.max_qp_allowed,
                        refQpIndexTemp);

                    if (hlRateControlHistogramPtrTemp->slice_type == I_SLICE) {
                        refQpIndexTemp = (uint32_t)MAX((int32_t)refQpIndexTemp + RC_INTRA_QP_OFFSET, 0);
                    }

                    hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = 0;

                    if (refQpTableIndex == previous_selected_ref_qp) {
                        hlRateControlHistogramPtrTemp->lifeCount--;
                    }
                    if (hlRateControlHistogramPtrTemp->isCoded) {
                        // If the frame is already coded, use the actual number of bits
                        hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->totalNumBitsCoded;
                    }
                    else {
                        rateControlTablesPtr = &encode_context_ptr->rate_control_tables_array[refQpIndexTemp];
                        sadBitsArrayPtr = rateControlTablesPtr->sadBitsArray[hlRateControlHistogramPtrTemp->temporal_layer_index];
                        intraSadBitsArrayPtr = rateControlTablesPtr->intraSadBitsArray[0];
                        pred_bits_ref_qp = 0;
                        numOfFullLcus = 0;

                        if (hlRateControlHistogramPtrTemp->slice_type == I_SLICE) {
                            // Loop over block in the frame and calculated the predicted bits at reg QP
                            unsigned i;
                            uint32_t accum = 0;
                            for (i = 0; i < NUMBER_OF_INTRA_SAD_INTERVALS; ++i)
                            {
                                accum += (uint32_t)(hlRateControlHistogramPtrTemp->ois_distortion_histogram[i] * intraSadBitsArrayPtr[i]);
                            }

                            pred_bits_ref_qp = accum;
                            numOfFullLcus = hlRateControlHistogramPtrTemp->full_sb_count;
                            hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] += pred_bits_ref_qp;
                        }
                        else {
                            unsigned i;
                            uint32_t accum = 0;
                            uint32_t accumIntra = 0;
                            for (i = 0; i < NUMBER_OF_SAD_INTERVALS; ++i)
                            {
                                accum += (uint32_t)(hlRateControlHistogramPtrTemp->me_distortion_histogram[i] * sadBitsArrayPtr[i]);
                                accumIntra += (uint32_t)(hlRateControlHistogramPtrTemp->ois_distortion_histogram[i] * intraSadBitsArrayPtr[i]);

                            }
                            if (accum > accumIntra * 3)
                                pred_bits_ref_qp = accumIntra;
                            else
                                pred_bits_ref_qp = accum;
                            numOfFullLcus = hlRateControlHistogramPtrTemp->full_sb_count;
                            hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] += pred_bits_ref_qp;
                        }

                        // Scale for in complete LCSs
                        //  pred_bits_ref_qp is normalized based on the area because of the LCUs at the picture boundries
                        hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] * (uint64_t)areaInPixel / (numOfFullLcus << 12);

                    }
                    highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] += hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];
                    // Store the pred_bits_ref_qp for the first frame in the window to PCS
                    if (queueEntryIndexHeadTemp == queueEntryIndexTemp2)
                        picture_control_set_ptr->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];

                    end_of_sequence_flag = hlRateControlHistogramPtrTemp->end_of_sequence_flag;
                    queueEntryIndexTemp++;
                }

                if (minLaBitDistance >= (uint64_t)ABS((int64_t)highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] - (int64_t)bitConstraintPerSw)) {
                    minLaBitDistance = (uint64_t)ABS((int64_t)highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] - (int64_t)bitConstraintPerSw);
                    selectedRefQpTableIndex = refQpTableIndex;
                    selectedRefQp = refQpIndex;
                }
                else {
                    bestQpFound = EB_TRUE;
                }

                if (refQpTableIndex == previous_selected_ref_qp) {
                    if (highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] > bitConstraintPerSw) {
                        qpStep = +1;
                    }
                    else {
                        qpStep = -1;
                    }
                }
                refQpTableIndex = (uint32_t)(refQpTableIndex + qpStep);

            }
        }

#if RC_UPDATE_TARGET_RATE
        selectedOrgRefQp = selectedRefQp;
        if (sequence_control_set_ptr->intra_period_length != -1 && picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0 &&
            (int32_t)picture_control_set_ptr->frames_in_sw > sequence_control_set_ptr->intra_period_length) {
            if (picture_control_set_ptr->picture_number > 0) {
                picture_control_set_ptr->intra_selected_org_qp = (uint8_t)selectedRefQp;
            }
            else {
                selectedOrgRefQp = selectedRefQp + 1;
                selectedRefQp = selectedRefQp + 1;
            }
            refQpIndex = selectedRefQp;
            highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] = 0;

            if (highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] == 0) {

                // Finding the predicted bits for each frame in the sliding window at the reference Qp(s)
                //queueEntryIndexTemp = encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                queueEntryIndexHeadTemp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
                queueEntryIndexHeadTemp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                queueEntryIndexHeadTemp = (queueEntryIndexHeadTemp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                    queueEntryIndexHeadTemp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                    queueEntryIndexHeadTemp;

                queueEntryIndexTemp = queueEntryIndexHeadTemp;

                // This is set to false, so the last frame would go inside the loop
                end_of_sequence_flag = EB_FALSE;

                while (!end_of_sequence_flag &&
                    //queueEntryIndexTemp <= encode_context_ptr->hl_rate_control_historgram_queue_head_index+sequence_control_set_ptr->static_config.look_ahead_distance){
                    queueEntryIndexTemp <= queueEntryIndexHeadTemp + sequence_control_set_ptr->static_config.look_ahead_distance) {

                    queueEntryIndexTemp2 = (queueEntryIndexTemp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? queueEntryIndexTemp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH : queueEntryIndexTemp;
                    hlRateControlHistogramPtrTemp = (encode_context_ptr->hl_rate_control_historgram_queue[queueEntryIndexTemp2]);


                    refQpIndexTemp = refQpIndex + QP_OFFSET_LAYER_ARRAY[hlRateControlHistogramPtrTemp->temporal_layer_index];
                    refQpIndexTemp = (uint32_t)CLIP3(
                        sequence_control_set_ptr->static_config.min_qp_allowed,
                        sequence_control_set_ptr->static_config.max_qp_allowed,
                        refQpIndexTemp);

                    if (hlRateControlHistogramPtrTemp->slice_type == I_SLICE) {
                        refQpIndexTemp = (uint32_t)MAX((int32_t)refQpIndexTemp + RC_INTRA_QP_OFFSET, 0);
                    }

                    hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = 0;

                    if (hlRateControlHistogramPtrTemp->isCoded) {
                        hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->totalNumBitsCoded;
                    }
                    else {
                        rateControlTablesPtr = &encode_context_ptr->rate_control_tables_array[refQpIndexTemp];
                        sadBitsArrayPtr = rateControlTablesPtr->sadBitsArray[hlRateControlHistogramPtrTemp->temporal_layer_index];
                        intraSadBitsArrayPtr = rateControlTablesPtr->intraSadBitsArray[hlRateControlHistogramPtrTemp->temporal_layer_index];
                        pred_bits_ref_qp = 0;

                        numOfFullLcus = 0;

                        if (hlRateControlHistogramPtrTemp->slice_type == I_SLICE) {
                            // Loop over block in the frame and calculated the predicted bits at reg QP

                            {
                                unsigned i;
                                uint32_t accum = 0;
                                for (i = 0; i < NUMBER_OF_INTRA_SAD_INTERVALS; ++i)
                                {
                                    accum += (uint32_t)(hlRateControlHistogramPtrTemp->ois_distortion_histogram[i] * intraSadBitsArrayPtr[i]);
                                }

                                pred_bits_ref_qp = accum;
                                numOfFullLcus = hlRateControlHistogramPtrTemp->full_sb_count;
                            }
                            hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] += pred_bits_ref_qp;
                        }

                        else {
                            unsigned i;
                            uint32_t accum = 0;
                            uint32_t accumIntra = 0;
                            for (i = 0; i < NUMBER_OF_SAD_INTERVALS; ++i)
                            {
                                accum += (uint32_t)(hlRateControlHistogramPtrTemp->me_distortion_histogram[i] * sadBitsArrayPtr[i]);
                                accumIntra += (uint32_t)(hlRateControlHistogramPtrTemp->ois_distortion_histogram[i] * intraSadBitsArrayPtr[i]);

                            }
                            if (accum > accumIntra * 3)
                                pred_bits_ref_qp = accumIntra;
                            else
                                pred_bits_ref_qp = accum;
                            numOfFullLcus = hlRateControlHistogramPtrTemp->full_sb_count;
                            hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] += pred_bits_ref_qp;
                        }

                        // Scale for in complete
                        //  pred_bits_ref_qp is normalized based on the area because of the LCUs at the picture boundries
                        hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] * (uint64_t)areaInPixel / (numOfFullLcus << 12);

                    }
                    highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] += hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];
                    // Store the pred_bits_ref_qp for the first frame in the window to PCS
                    //  if(encode_context_ptr->hl_rate_control_historgram_queue_head_index == queueEntryIndexTemp2)
                    if (queueEntryIndexHeadTemp == queueEntryIndexTemp2)
                        picture_control_set_ptr->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];

                    end_of_sequence_flag = hlRateControlHistogramPtrTemp->end_of_sequence_flag;
                    queueEntryIndexTemp++;
                }
            }
        }
#endif
        picture_control_set_ptr->tables_updated = tables_updated;
        EbBool expensiveISlice = EB_FALSE;
        // Looping over the window to find the percentage of bit allocation in each layer
        if ((sequence_control_set_ptr->intra_period_length != -1) &&
            ((int32_t)picture_control_set_ptr->frames_in_sw > sequence_control_set_ptr->intra_period_length) &&
            ((int32_t)picture_control_set_ptr->frames_in_sw > sequence_control_set_ptr->intra_period_length)) {
            uint64_t iSliceBits = 0;

            if (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0) {

                queueEntryIndexHeadTemp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
                queueEntryIndexHeadTemp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                queueEntryIndexHeadTemp = (queueEntryIndexHeadTemp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                    queueEntryIndexHeadTemp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                    queueEntryIndexHeadTemp;

                queueEntryIndexTemp = queueEntryIndexHeadTemp;

                // This is set to false, so the last frame would go inside the loop
                end_of_sequence_flag = EB_FALSE;

                while (!end_of_sequence_flag &&
                    queueEntryIndexTemp <= queueEntryIndexHeadTemp + sequence_control_set_ptr->intra_period_length) {

                    queueEntryIndexTemp2 = (queueEntryIndexTemp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? queueEntryIndexTemp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH : queueEntryIndexTemp;
                    hlRateControlHistogramPtrTemp = (encode_context_ptr->hl_rate_control_historgram_queue[queueEntryIndexTemp2]);

                    refQpIndexTemp = selectedRefQp + QP_OFFSET_LAYER_ARRAY[hlRateControlHistogramPtrTemp->temporal_layer_index];
                    refQpIndexTemp = (uint32_t)CLIP3(
                        sequence_control_set_ptr->static_config.min_qp_allowed,
                        sequence_control_set_ptr->static_config.max_qp_allowed,
                        refQpIndexTemp);

                    if (hlRateControlHistogramPtrTemp->slice_type == I_SLICE) {
                        refQpIndexTemp = (uint32_t)MAX((int32_t)refQpIndexTemp + RC_INTRA_QP_OFFSET, 0);
                    }
                    if (queueEntryIndexTemp == queueEntryIndexHeadTemp) {
                        iSliceBits = hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];
                    }
                    picture_control_set_ptr->total_bits_per_gop += hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];
                    picture_control_set_ptr->bits_per_sw_per_layer[hlRateControlHistogramPtrTemp->temporal_layer_index] += hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];
                    picture_control_set_ptr->percentage_updated = EB_TRUE;

                    end_of_sequence_flag = hlRateControlHistogramPtrTemp->end_of_sequence_flag;
                    queueEntryIndexTemp++;
                }
                if (iSliceBits * 100 > 85 * picture_control_set_ptr->total_bits_per_gop) {
                    expensiveISlice = EB_TRUE;
                }
                if (picture_control_set_ptr->total_bits_per_gop == 0) {
                    for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++) {
                        picture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_index] = RATE_PERCENTAGE_LAYER_ARRAY[sequence_control_set_ptr->static_config.hierarchical_levels][temporal_layer_index];
                    }
                }
            }
        }
        else {
            for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++) {
                picture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_index] = RATE_PERCENTAGE_LAYER_ARRAY[sequence_control_set_ptr->static_config.hierarchical_levels][temporal_layer_index];
            }
        }
        if (expensiveISlice) {
            if (tables_updated) {
                selectedRefQp = (uint32_t)MAX((int32_t)selectedRefQp - 1, 0);
            }
            else {
                selectedRefQp = (uint32_t)MAX((int32_t)selectedRefQp - 3, 0);
            }
            selectedRefQp = (uint32_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                sequence_control_set_ptr->static_config.max_qp_allowed,
                selectedRefQp);
        }
        // Set the QP
        previous_selected_ref_qp = selectedRefQp;
        if (picture_control_set_ptr->picture_number > max_coded_poc && picture_control_set_ptr->temporal_layer_index < 2 && !picture_control_set_ptr->end_of_sequence_region) {

            max_coded_poc = picture_control_set_ptr->picture_number;
            max_coded_poc_selected_ref_qp = previous_selected_ref_qp;
            encode_context_ptr->previous_selected_ref_qp = previous_selected_ref_qp;
            encode_context_ptr->max_coded_poc = max_coded_poc;
            encode_context_ptr->max_coded_poc_selected_ref_qp = max_coded_poc_selected_ref_qp;

        }
        picture_control_set_ptr->best_pred_qp = (uint8_t)CLIP3(
            sequence_control_set_ptr->static_config.min_qp_allowed,
            sequence_control_set_ptr->static_config.max_qp_allowed,
            selectedRefQp + QP_OFFSET_LAYER_ARRAY[picture_control_set_ptr->temporal_layer_index]);
        if (picture_control_set_ptr->slice_type == I_SLICE) {
            picture_control_set_ptr->best_pred_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->best_pred_qp + RC_INTRA_QP_OFFSET, 0);
        }
#if RC_UPDATE_TARGET_RATE
        if (picture_control_set_ptr->picture_number == 0) {
            highLevelRateControlPtr->prevIntraSelectedRefQp = selectedRefQp;
            highLevelRateControlPtr->prevIntraOrgSelectedRefQp = selectedRefQp;
        }
        if (sequence_control_set_ptr->intra_period_length != -1) {
            if (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0) {
                highLevelRateControlPtr->prevIntraSelectedRefQp = selectedRefQp;
                highLevelRateControlPtr->prevIntraOrgSelectedRefQp = selectedOrgRefQp;
            }
        }
#endif
        picture_control_set_ptr->target_bits_best_pred_qp = picture_control_set_ptr->pred_bits_ref_qp[picture_control_set_ptr->best_pred_qp];
        //if (picture_control_set_ptr->slice_type == 2)
        // {
        //printf("\nTID: %d\t", picture_control_set_ptr->temporal_layer_index);
        //printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
        //    picture_control_set_ptr->picture_number,
        //    picture_control_set_ptr->best_pred_qp,
        //    (int32_t)picture_control_set_ptr->target_bits_best_pred_qp,
        //    (int32_t)highLevelRateControlPtr->predBitsRefQpPerSw[selectedRefQp - 1],
        //    (int32_t)highLevelRateControlPtr->predBitsRefQpPerSw[selectedRefQp],
        //    (int32_t)highLevelRateControlPtr->predBitsRefQpPerSw[selectedRefQp + 1],
        //    (int32_t)highLevelRateControlPtr->bitConstraintPerSw,
        //    (int32_t)bitConstraintPerSw,
        //    (int32_t)highLevelRateControlPtr->virtualBufferLevel);
        //}
    }
    EbReleaseMutex(sequence_control_set_ptr->encode_context_ptr->rate_table_update_mutex);
}

#if ADD_DELTA_QP_SUPPORT ||  NEW_QPS

static const uint8_t quantizer_to_qindex[] = {
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48,
    52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100,
    104, 108, 112, 116, 120, 124, 128, 132, 136, 140, 144, 148, 152,
    156, 160, 164, 168, 172, 176, 180, 184, 188, 192, 196, 200, 204,
    208, 212, 216, 220, 224, 228, 232, 236, 240, 244, 249, 255,
};
#endif

#if NEW_QPS
#define MAX_Q_INDEX 255
#define MIN_Q_INDEX 0

extern int16_t av1_ac_quant_Q3(int32_t qindex, int32_t delta, aom_bit_depth_t bit_depth);
// These functions use formulaic calculations to make playing with the
// quantizer tables easier. If necessary they can be replaced by lookup
// tables if and when things settle down in the experimental bitstream

double av1_convert_qindex_to_q(int32_t qindex, aom_bit_depth_t bit_depth) {
    // Convert the index to a real Q value (scaled down to match old Q values)
    switch (bit_depth) {
    case AOM_BITS_8: return av1_ac_quant_Q3(qindex, 0, bit_depth) / 4.0;
    case AOM_BITS_10: return av1_ac_quant_Q3(qindex, 0, bit_depth) / 16.0;
    case AOM_BITS_12: return av1_ac_quant_Q3(qindex, 0, bit_depth) / 64.0;
    default:
        assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12");
        return -1.0;
    }
}
int32_t av1_compute_qdelta(double qstart, double qtarget,
    aom_bit_depth_t bit_depth) {
    int32_t start_index = MAX_Q_INDEX;
    int32_t target_index = MAX_Q_INDEX;
    int32_t i;

    // Convert the average q value to an index.
    for (i = MIN_Q_INDEX; i < MAX_Q_INDEX; ++i) {
        start_index = i;
        if (av1_convert_qindex_to_q(i, bit_depth) >= qstart) break;
    }

    // Convert the q target to an index
    for (i = MIN_Q_INDEX; i < MAX_Q_INDEX; ++i) {
        target_index = i;
        if (av1_convert_qindex_to_q(i, bit_depth) >= qtarget) break;
    }

    return target_index - start_index;
}
#endif
#if 0 //NEW_PRED
#define STATIC_MOTION_THRESH 95

enum {
    INTER_NORMAL = 0,
    INTER_LOW = 1,
    INTER_HIGH = 2,
    GF_ARF_LOW = 3,
    GF_ARF_STD = 4,
    KF_STD = 5,
    RATE_FACTOR_LEVELS = 6
} UENUM1BYTE(RATE_FACTOR_LEVEL);

enum {
    KF_UPDATE = 0,
    LF_UPDATE = 1,
    GF_UPDATE = 2,
    ARF_UPDATE = 3,
    OVERLAY_UPDATE = 4,
    BRF_UPDATE = 5,            // Backward Reference Frame
    LAST_BIPRED_UPDATE = 6,    // Last Bi-predictive Frame
    BIPRED_UPDATE = 7,         // Bi-predictive Frame, but not the last one
    INTNL_OVERLAY_UPDATE = 8,  // Internal Overlay Frame
    INTNL_ARF_UPDATE = 9,      // Internal Altref Frame (candidate for ALTREF2)
    FRAME_UPDATE_TYPES = 10
} UENUM1BYTE(FRAME_UPDATE_TYPE);

// that are not marked as coded with 0,0 motion in the first pass.
#define STATIC_KF_GROUP_THRESH 99

static int rc_pick_q_and_bounds_two_pass(
    PictureControlSet_t         *picture_control_set_ptr,
    //const AV1_COMP *cpi, 
    //int width,
    //int height, 
    int *bottom_index,
    int *top_index, 
    int *arf_q) {

    // Config
    SequenceControlSet_t        *sequence_control_set_ptr = picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr;
    //picture_control_set_ptr->parent_pcs_ptr->av1_cm
    const Av1Common  *const cm = &picture_control_set_ptr->parent_pcs_ptr->av1_cm;
    //const RATE_CONTROL *const rc = &cpi->rc;
    //const AV1EncoderConfig *const oxcf = &cpi->oxcf;
    //const GF_GROUP *gf_group = &cpi->twopass.gf_group;
 //   const int cq_level = get_active_cq_level(rc, oxcf, frame_is_intra_only(cm),
 //       cm->superres_scale_denominator);
    int active_best_quality;
    int active_worst_quality = cpi->twopass.active_worst_quality;
    int q;
    int *inter_minq;
    int is_src_frame_alt_ref, refresh_golden_frame, refresh_alt_ref_frame, new_bwdref_update_rule;


    const int bit_depth = sequence_control_set_ptr->static_config.encoder_bit_depth;//cm->seq_params.bit_depth;
    ASSIGN_MINQ_TABLE(bit_depth, inter_minq);

    const int is_intrl_arf_boost = 0;
        //gf_group->update_type[gf_group->index] == INTNL_ARF_UPDATE;

    if (frame_is_intra_only(cm)) {
#if 0
        if (rc->frames_to_key == 1 && oxcf->rc_mode == AOM_Q) {
            // If the next frame is also a key frame or the current frame is the
            // only frame in the sequence in AOM_Q mode, just use the cq_level
            // as q.
            active_best_quality = cq_level;
            active_worst_quality = cq_level;
        }
        else if (rc->this_key_frame_forced) {
            // Handle the special case for key frames forced when we have reached
            // the maximum key frame interval. Here force the Q to a range
            // based on the ambient Q to reduce the risk of popping.
            double last_boosted_q;
            int delta_qindex;
            int qindex;

            if (cpi->twopass.last_kfgroup_zeromotion_pct >= STATIC_MOTION_THRESH) {
                qindex = AOMMIN(rc->last_kf_qindex, rc->last_boosted_qindex);
                active_best_quality = qindex;
                last_boosted_q = av1_convert_qindex_to_q(qindex, bit_depth);
                delta_qindex = av1_compute_qdelta(rc, last_boosted_q,
                    last_boosted_q * 1.25, bit_depth);
                active_worst_quality =
                    AOMMIN(qindex + delta_qindex, active_worst_quality);
            }
            else {
                // Increase the boost if the forced keyframe is a forward reference.
                // These numbers were derived empirically.
                const double boost_factor = cpi->oxcf.fwd_kf_enabled ? 0.25 : 0.50;
                qindex = rc->last_boosted_qindex;
                last_boosted_q = av1_convert_qindex_to_q(qindex, bit_depth);
                delta_qindex = av1_compute_qdelta(
                    rc, last_boosted_q, last_boosted_q * boost_factor, bit_depth);
                active_best_quality = AOMMAX(qindex + delta_qindex, rc->best_quality);
            }
        }
        else
#endif
        {
            // Not forced keyframe.
            double q_adj_factor = 1.0;
            double q_val;

            cpi->rc.worst_quality = MAXQ;
            cpi->rc.best_quality = MINQ;

            // Hsan: cross multiplication to derive kf_boost from non_moving_average_score; kf_boost range is [kf_low,kf_high], and non_moving_average_score range [NON_MOVING_SCORE_0,NON_MOVING_SCORE_3]
            cpi->rc.kf_boost = (((NON_MOVING_SCORE_3 - picture_control_set_ptr->parent_pcs_ptr->non_moving_index_average)  * (kf_high - kf_low)) / NON_MOVING_SCORE_3) + kf_low;


            // Baseline value derived from cpi->active_worst_quality and kf boost.
            active_best_quality =
                get_kf_active_quality(rc, active_worst_quality, bit_depth);
            if (picture_control_set_ptr->parent_pcs_ptr->kf_zeromotion_pct >= STATIC_KF_GROUP_THRESH) {
            //if (cpi->twopass.kf_zeromotion_pct >= STATIC_KF_GROUP_THRESH) {
                active_best_quality /= 3;
            }

            // Allow somewhat lower kf minq with small image formats.
            if ((cm->width * cm->height) <= (352 * 288)) {
                q_adj_factor -= 0.25;
            }

            // Make a further adjustment based on the kf zero motion measure.
            q_adj_factor += 0.05 - (0.001 * (double)picture_control_set_ptr->parent_pcs_ptr->kf_zeromotion_pct/*(double)cpi->twopass.kf_zeromotion_pct*/);

            // Convert the adjustment factor to a qindex delta
            // on active_best_quality.
            q_val = av1_convert_qindex_to_q(active_best_quality, bit_depth);
            active_best_quality +=
                av1_compute_qdelta(q_val, q_val * q_adj_factor, bit_depth);
        }
    }
#if 0 // disable for now
    else if (!is_src_frame_alt_ref &&
        (refresh_golden_frame || is_intrl_arf_boost ||
            refresh_alt_ref_frame)) {
        // Use the lower of active_worst_quality and recent
        // average Q as basis for GF/ARF best Q limit unless last frame was
        // a key frame.
        if (rc->frames_since_key > 1 &&
            rc->avg_frame_qindex[INTER_FRAME] < active_worst_quality) {
            q = rc->avg_frame_qindex[INTER_FRAME];
        }
        else {
            q = active_worst_quality;
        }
        // For constrained quality dont allow Q less than the cq level
#if 0
        if (oxcf->rc_mode == AOM_CQ) {
            if (q < cq_level) q = cq_level;

            active_best_quality = get_gf_active_quality(rc, q, bit_depth);

            // Constrained quality use slightly lower active best.
            active_best_quality = active_best_quality * 15 / 16;

            if (gf_group->update_type[gf_group->index] == ARF_UPDATE ||
                (is_intrl_arf_boost && !cpi->new_bwdref_update_rule)) {
                if (gf_group->update_type[gf_group->index] == ARF_UPDATE) {
                    const int min_boost = get_gf_high_motion_quality(q, bit_depth);
                    const int boost = min_boost - active_best_quality;

                    active_best_quality = min_boost - (int)(boost * rc->arf_boost_factor);
                }
                *arf_q = active_best_quality;
            }
            else if (cpi->new_bwdref_update_rule && is_intrl_arf_boost) {
                assert(rc->arf_q >= 0);  // Ensure it is set to a valid value.
                active_best_quality = rc->arf_q;
                int this_height = gf_group_pyramid_level(cpi);
                while (this_height < gf_group->pyramid_height) {
                    active_best_quality = (active_best_quality + cq_level + 1) / 2;
                    ++this_height;
                }
            }
        }
        else if (oxcf->rc_mode == AOM_Q)
#endif
        {
            if (!refresh_alt_ref_frame && !is_intrl_arf_boost) {
                active_best_quality = cq_level;
            }
            else {
                if (gf_group->update_type[gf_group->index] == ARF_UPDATE) {
                    active_best_quality = get_gf_active_quality(rc, q, bit_depth);
                    *arf_q = active_best_quality;
                    const int min_boost = get_gf_high_motion_quality(q, bit_depth);
                    const int boost = min_boost - active_best_quality;

                    active_best_quality = min_boost - (int)(boost * rc->arf_boost_factor);
                }
                else {
                    assert(rc->arf_q >= 0);  // Ensure it is set to a valid value.
                    active_best_quality = rc->arf_q;
                }
                if (new_bwdref_update_rule && is_intrl_arf_boost) {
                    int this_height = gf_group_pyramid_level(cpi);
                    while (this_height < gf_group->pyramid_height) {
                        active_best_quality = (active_best_quality + cq_level + 1) / 2;
                        ++this_height;
                    }
                }
                else {
                    // Modify best quality for second level arfs. For mode AOM_Q this
                    // becomes the baseline frame q.
                    if (gf_group->rf_level[gf_group->index] == GF_ARF_LOW)
                        active_best_quality = (active_best_quality + cq_level + 1) / 2;
                }
            }
        }
#if 0
        else {
            active_best_quality = get_gf_active_quality(rc, q, bit_depth);
            const int min_boost = get_gf_high_motion_quality(q, bit_depth);
            const int boost = min_boost - active_best_quality;

            active_best_quality = min_boost - (int)(boost * rc->arf_boost_factor);
            if (cpi->new_bwdref_update_rule && is_intrl_arf_boost) {
                int this_height = gf_group_pyramid_level(cpi);
                while (this_height < gf_group->pyramid_height) {
                    active_best_quality =
                        (active_best_quality + active_worst_quality + 1) / 2;
                    ++this_height;
                }
            }
        }
#endif
    }
    else {
#if 0
        if (oxcf->rc_mode == AOM_Q)
#endif
        {
            active_best_quality = cq_level;
        }
#if 0
        else {
            active_best_quality = inter_minq[active_worst_quality];

            // For the constrained quality mode we don't want
            // q to fall below the cq level.
            if ((oxcf->rc_mode == AOM_CQ) && (active_best_quality < cq_level)) {
                active_best_quality = cq_level;
            }
        }
#endif
    }
#endif // disable for now
#if 1 // AMIR_PRINT
    printf("\t%d\t%d\t%d\t%d\t%d\t",
        is_src_frame_alt_ref,
        refresh_golden_frame,
        is_intrl_arf_boost,
        refresh_alt_ref_frame,
        new_bwdref_update_rule);
    fflush(stdout);
#endif
#if 0
    // Extension to max or min Q if undershoot or overshoot is outside
    // the permitted range.
    if ((cpi->oxcf.rc_mode != AOM_Q) &&
        (cpi->twopass.gf_zeromotion_pct < VLOW_MOTION_THRESHOLD)) {
        if (frame_is_intra_only(cm) ||
            (!rc->is_src_frame_alt_ref &&
            (cpi->refresh_golden_frame || is_intrl_arf_boost ||
                cpi->refresh_alt_ref_frame))) {
            active_best_quality -=
                (cpi->twopass.extend_minq + cpi->twopass.extend_minq_fast);
            active_worst_quality += (cpi->twopass.extend_maxq / 2);
        }
        else {
            active_best_quality -=
                (cpi->twopass.extend_minq + cpi->twopass.extend_minq_fast) / 2;
            active_worst_quality += cpi->twopass.extend_maxq;
        }
    }

    aom_clear_system_state();
#endif
    // Static forced key frames Q restrictions dealt with elsewhere.
    if (!(frame_is_intra_only(cm)) || !rc->this_key_frame_forced ||
        (cpi->twopass.last_kfgroup_zeromotion_pct < STATIC_MOTION_THRESH)) {
        int qdelta = av1_frame_type_qdelta(cpi, gf_group->rf_level[gf_group->index],
            active_worst_quality);
        active_worst_quality =
            AOMMAX(active_worst_quality + qdelta, active_best_quality);
    }
#if 0
    // Modify active_best_quality for downscaled normal frames.
    if (av1_frame_scaled(cm) && !frame_is_kf_gf_arf(cpi)) {
        int qdelta = av1_compute_qdelta_by_rate(
            rc, cm->current_frame.frame_type, active_best_quality, 2.0, bit_depth);
        active_best_quality =
            AOMMAX(active_best_quality + qdelta, rc->best_quality);
    }
#endif
    active_best_quality =
        clamp(active_best_quality, rc->best_quality, rc->worst_quality);
    active_worst_quality =
        clamp(active_worst_quality, active_best_quality, rc->worst_quality);
#if 0
    if (oxcf->rc_mode == AOM_Q ||
        (frame_is_intra_only(cm) && !rc->this_key_frame_forced &&
            cpi->twopass.kf_zeromotion_pct >= STATIC_KF_GROUP_THRESH))
#endif
    {
        q = active_best_quality;
        // Special case code to try and match quality with forced key frames.
    }
#if 0
    else if (frame_is_intra_only(cm) && rc->this_key_frame_forced) {
        // If static since last kf use better of last boosted and last kf q.
        if (cpi->twopass.last_kfgroup_zeromotion_pct >= STATIC_MOTION_THRESH) {
            q = AOMMIN(rc->last_kf_qindex, rc->last_boosted_qindex);
        }
        else {
            q = AOMMIN(rc->last_boosted_qindex,
                (active_best_quality + active_worst_quality) / 2);
        }
    }
    else {
        q = av1_rc_regulate_q(cpi, rc->this_frame_target, active_best_quality,
            active_worst_quality, width, height);
        if (q > active_worst_quality) {
            // Special case when we are targeting the max allowed rate.
            if (rc->this_frame_target >= rc->max_frame_bandwidth)
                active_worst_quality = q;
            else
                q = active_worst_quality;
        }
    }
#endif
    clamp(q, active_best_quality, active_worst_quality);

    *top_index = active_worst_quality;
    *bottom_index = active_best_quality;

    assert(*top_index <= rc->worst_quality && *top_index >= rc->best_quality);
    assert(*bottom_index <= rc->worst_quality &&
        *bottom_index >= rc->best_quality);
    assert(q <= rc->worst_quality && q >= rc->best_quality);
    return q;
}
#if 0
int av1_rc_pick_q_and_bounds(AV1_COMP *cpi, int width, int height,
    int *bottom_index, int *top_index) {
    int q;
    if (cpi->oxcf.pass == 0) {
        if (cpi->oxcf.rc_mode == AOM_CBR)
            q = rc_pick_q_and_bounds_one_pass_cbr(cpi, width, height, bottom_index,
                top_index);
        else
            q = rc_pick_q_and_bounds_one_pass_vbr(cpi, width, height, bottom_index,
                top_index);
    }
    else {
        assert(cpi->oxcf.pass == 2 && "invalid encode pass");

        GF_GROUP *gf_group = &cpi->twopass.gf_group;
        int arf_q = -1;  // Initialize to invalid value, for sanity check later.

        q = rc_pick_q_and_bounds_two_pass(cpi, width, height, bottom_index,
            top_index, &arf_q);

        if (gf_group->update_type[gf_group->index] == ARF_UPDATE) {
            cpi->rc.arf_q = arf_q;
        }
    }

    return q;
}
#endif


#endif
void* RateControlKernel(void *input_ptr)
{
    // Context
    RateControlContext_t        *context_ptr = (RateControlContext_t*)input_ptr;
    // EncodeContext_t             *encode_context_ptr;

    RateControlIntervalParamContext_t *rateControlParamPtr;

    RateControlIntervalParamContext_t *prevGopRateControlParamPtr;
    RateControlIntervalParamContext_t *nextGopRateControlParamPtr;

    PictureControlSet_t         *picture_control_set_ptr;
    PictureParentControlSet_t   *parentPictureControlSetPtr;

    // Config
    SequenceControlSet_t        *sequence_control_set_ptr;

    // Input
    EbObjectWrapper_t           *rateControlTasksWrapperPtr;
    RateControlTasks_t          *rateControlTasksPtr;

    // Output
    EbObjectWrapper_t           *rateControlResultsWrapperPtr;
    RateControlResults_t        *rateControlResultsPtr;

    // SB Loop variables
    LargestCodingUnit_t         *sb_ptr;
    uint32_t                       lcuCodingOrder;
    uint64_t                       totalNumberOfFbFrames = 0;

    RATE_CONTROL_TASKTYPES       taskType;
    EbRateControlModel          *rc_model_ptr;

    rate_control_model_ctor(&rc_model_ptr);

    for (;;) {

        // Get RateControl Task
        EbGetFullObject(
            context_ptr->rateControlInputTasksFifoPtr,
            &rateControlTasksWrapperPtr);

        rateControlTasksPtr = (RateControlTasks_t*)rateControlTasksWrapperPtr->objectPtr;
        taskType = rateControlTasksPtr->taskType;

        // Modify these for different temporal layers later
        switch (taskType) {

        case RC_PICTURE_MANAGER_RESULT:

            picture_control_set_ptr = (PictureControlSet_t*)rateControlTasksPtr->pictureControlSetWrapperPtr->objectPtr;
            sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;

            // High level RC
            if (picture_control_set_ptr->picture_number == 0) {

                rate_control_model_init(rc_model_ptr, sequence_control_set_ptr);
                context_ptr->highLevelRateControlPtr->target_bit_rate = sequence_control_set_ptr->static_config.target_bit_rate;
                context_ptr->highLevelRateControlPtr->frame_rate = sequence_control_set_ptr->frame_rate;
                context_ptr->highLevelRateControlPtr->channelBitRatePerFrame = (uint64_t)MAX((int64_t)1, (int64_t)((context_ptr->highLevelRateControlPtr->target_bit_rate << RC_PRECISION) / context_ptr->highLevelRateControlPtr->frame_rate));

                context_ptr->highLevelRateControlPtr->channelBitRatePerSw = context_ptr->highLevelRateControlPtr->channelBitRatePerFrame * (sequence_control_set_ptr->static_config.look_ahead_distance + 1);
                context_ptr->highLevelRateControlPtr->bitConstraintPerSw = context_ptr->highLevelRateControlPtr->channelBitRatePerSw;

#if RC_UPDATE_TARGET_RATE
                context_ptr->highLevelRateControlPtr->previousUpdatedBitConstraintPerSw = context_ptr->highLevelRateControlPtr->channelBitRatePerSw;
#endif

                int32_t totalFrameInInterval = sequence_control_set_ptr->intra_period_length;
                uint32_t gopPeriod = (1 << picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels);
                context_ptr->frame_rate = sequence_control_set_ptr->frame_rate;
                while (totalFrameInInterval >= 0) {
                    if (totalFrameInInterval % (gopPeriod) == 0)
                        context_ptr->frames_in_interval[0] ++;
                    else if (totalFrameInInterval % (gopPeriod >> 1) == 0)
                        context_ptr->frames_in_interval[1] ++;
                    else if (totalFrameInInterval % (gopPeriod >> 2) == 0)
                        context_ptr->frames_in_interval[2] ++;
                    else if (totalFrameInInterval % (gopPeriod >> 3) == 0)
                        context_ptr->frames_in_interval[3] ++;
                    else if (totalFrameInInterval % (gopPeriod >> 4) == 0)
                        context_ptr->frames_in_interval[4] ++;
                    else if (totalFrameInInterval % (gopPeriod >> 5) == 0)
                        context_ptr->frames_in_interval[5] ++;
                    totalFrameInInterval--;
                }
                context_ptr->virtualBufferSize = (((uint64_t)sequence_control_set_ptr->static_config.target_bit_rate * 3) << RC_PRECISION) / (context_ptr->frame_rate);
                context_ptr->rateAveragePeriodinFrames = (uint64_t)sequence_control_set_ptr->static_config.intra_period_length + 1;
                context_ptr->virtualBufferLevelInitialValue = context_ptr->virtualBufferSize >> 1;
                context_ptr->virtualBufferLevel = context_ptr->virtualBufferSize >> 1;
                context_ptr->previousVirtualBufferLevel = context_ptr->virtualBufferSize >> 1;
                context_ptr->vbFillThreshold1 = (context_ptr->virtualBufferSize * 6) >> 3;
                context_ptr->vbFillThreshold2 = (context_ptr->virtualBufferSize << 3) >> 3;
                context_ptr->baseLayerFramesAvgQp = sequence_control_set_ptr->qp;
                context_ptr->baseLayerIntraFramesAvgQp = sequence_control_set_ptr->qp;
            }
            if (sequence_control_set_ptr->static_config.rate_control_mode)
            {
                picture_control_set_ptr->parent_pcs_ptr->intra_selected_org_qp = 0;
                HighLevelRcInputPictureMode2(
                    picture_control_set_ptr->parent_pcs_ptr,
                    sequence_control_set_ptr,
                    sequence_control_set_ptr->encode_context_ptr,
                    context_ptr,
                    context_ptr->highLevelRateControlPtr);


            }

            // Frame level RC
            if (sequence_control_set_ptr->intra_period_length == -1 || sequence_control_set_ptr->static_config.rate_control_mode == 0) {
                rateControlParamPtr = context_ptr->rateControlParamQueue[0];
                prevGopRateControlParamPtr = context_ptr->rateControlParamQueue[0];
                nextGopRateControlParamPtr = context_ptr->rateControlParamQueue[0];
            }
            else {
                uint32_t intervalIndexTemp = 0;
                EbBool intervalFound = EB_FALSE;
                while ((intervalIndexTemp < PARALLEL_GOP_MAX_NUMBER) && !intervalFound) {

                    if (picture_control_set_ptr->picture_number >= context_ptr->rateControlParamQueue[intervalIndexTemp]->firstPoc &&
                        picture_control_set_ptr->picture_number <= context_ptr->rateControlParamQueue[intervalIndexTemp]->lastPoc) {
                        intervalFound = EB_TRUE;
                    }
                    else {
                        intervalIndexTemp++;
                    }
                }
                CHECK_REPORT_ERROR(
                    intervalIndexTemp != PARALLEL_GOP_MAX_NUMBER,
                    sequence_control_set_ptr->encode_context_ptr->app_callback_ptr,
                    EB_ENC_RC_ERROR2);

                rateControlParamPtr = context_ptr->rateControlParamQueue[intervalIndexTemp];

                prevGopRateControlParamPtr = (intervalIndexTemp == 0) ?
                    context_ptr->rateControlParamQueue[PARALLEL_GOP_MAX_NUMBER - 1] :
                    context_ptr->rateControlParamQueue[intervalIndexTemp - 1];
                nextGopRateControlParamPtr = (intervalIndexTemp == PARALLEL_GOP_MAX_NUMBER - 1) ?
                    context_ptr->rateControlParamQueue[0] :
                    context_ptr->rateControlParamQueue[intervalIndexTemp + 1];
            }

            if (sequence_control_set_ptr->static_config.rate_control_mode == 0) {
                // if RC mode is 0,  fixed QP is used
                // QP scaling based on POC number for Flat IPPP structure
#if NEW_QPS
                picture_control_set_ptr->parent_pcs_ptr->base_qindex = quantizer_to_qindex[picture_control_set_ptr->picture_qp];
#endif
                if (sequence_control_set_ptr->static_config.enable_qp_scaling_flag && picture_control_set_ptr->parent_pcs_ptr->qp_on_the_fly == EB_FALSE) {
#if NEW_QPS
                    const int32_t qindex = quantizer_to_qindex[(uint8_t)sequence_control_set_ptr->qp];
                    const double q_val = av1_convert_qindex_to_q(qindex, (aom_bit_depth_t)sequence_control_set_ptr->static_config.encoder_bit_depth);
                    if (picture_control_set_ptr->slice_type == I_SLICE) {

                        const int32_t delta_qindex = av1_compute_qdelta(
                            q_val,
                            q_val * 0.25,
                            (aom_bit_depth_t)sequence_control_set_ptr->static_config.encoder_bit_depth);
                        picture_control_set_ptr->parent_pcs_ptr->base_qindex =
                            (uint8_t)CLIP3(
                            (int32_t)quantizer_to_qindex[sequence_control_set_ptr->static_config.min_qp_allowed],
                                (int32_t)quantizer_to_qindex[sequence_control_set_ptr->static_config.max_qp_allowed],
                                (int32_t)(qindex + delta_qindex));

                    }
                    else {
                        // const double delta_rate[FIXED_GF_INTERVAL] = { 0.50, 1.0, 0.85, 1.0,
                        //     0.70, 1.0, 0.85, 1.0 };
#if NEW_PRED                    
                        double delta_rate_new[2][6] = 
                                { { 0.40, 0.7, 0.85, 1.0, 1.0, 1.0 },
                                { 0.35, 0.6, 0.8,  0.9, 1.0, 1.0 } };

#else
                       const double delta_rate_new[6] = { 0.40, 0.7, 0.85, 1.0, 1.0, 1.0 };

#endif
                        const int32_t delta_qindex = av1_compute_qdelta(
                            q_val,
#if NEW_PRED
                            q_val * delta_rate_new[picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels == 4][picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index],

#else
                            q_val * delta_rate_new[picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index],
#endif     
                            (aom_bit_depth_t)sequence_control_set_ptr->static_config.encoder_bit_depth);

                        picture_control_set_ptr->parent_pcs_ptr->base_qindex =
                            (uint8_t)CLIP3(
                            (int32_t)quantizer_to_qindex[sequence_control_set_ptr->static_config.min_qp_allowed],
                                (int32_t)quantizer_to_qindex[sequence_control_set_ptr->static_config.max_qp_allowed],
                                (int32_t)(qindex + delta_qindex));

                    }
#endif

                }

                else if (picture_control_set_ptr->parent_pcs_ptr->qp_on_the_fly == EB_TRUE) {

                    picture_control_set_ptr->picture_qp = (uint8_t)CLIP3((int32_t)sequence_control_set_ptr->static_config.min_qp_allowed, (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed, picture_control_set_ptr->parent_pcs_ptr->picture_qp);
#if NEW_QPS
                    picture_control_set_ptr->parent_pcs_ptr->base_qindex = quantizer_to_qindex[picture_control_set_ptr->picture_qp];
#endif
                }

                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp;

            }
            else {
                picture_control_set_ptr->picture_qp = rate_control_get_quantizer(rc_model_ptr, picture_control_set_ptr->parent_pcs_ptr);

                if (picture_control_set_ptr->picture_number == rateControlParamPtr->firstPoc && picture_control_set_ptr->picture_number != 0 && !prevGopRateControlParamPtr->scene_change_in_gop) {
                    int16_t deltaApQp = (int16_t)prevGopRateControlParamPtr->firstPicActualQp - (int16_t)prevGopRateControlParamPtr->firstPicPredQp;
                    int64_t extraApBitRatio = (prevGopRateControlParamPtr->firstPicPredBits != 0) ?
                        (((int64_t)prevGopRateControlParamPtr->firstPicActualBits - (int64_t)prevGopRateControlParamPtr->firstPicPredBits) * 100) / ((int64_t)prevGopRateControlParamPtr->firstPicPredBits) :
                        0;
                    extraApBitRatio += (int64_t)deltaApQp * 15;
                    if (extraApBitRatio > 200) {
                        picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + 3;
                    }
                    else if (extraApBitRatio > 100) {
                        picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + 2;
                    }
                    else if (extraApBitRatio > 50) {
                        picture_control_set_ptr->picture_qp++;
                    }
                }

                if (picture_control_set_ptr->picture_number == rateControlParamPtr->firstPoc && picture_control_set_ptr->picture_number != 0) {
                    uint8_t qpIncAllowed = 3;
                    uint8_t qpDecAllowed = 4;
                    if (picture_control_set_ptr->parent_pcs_ptr->intra_selected_org_qp + 10 <= prevGopRateControlParamPtr->firstPicActualQp)
                    {
                        qpDecAllowed = (uint8_t)(prevGopRateControlParamPtr->firstPicActualQp - picture_control_set_ptr->parent_pcs_ptr->intra_selected_org_qp) >> 1;
                    }

                    if (picture_control_set_ptr->parent_pcs_ptr->intra_selected_org_qp >= prevGopRateControlParamPtr->firstPicActualQp + 10)
                    {
                        qpIncAllowed = (uint8_t)(picture_control_set_ptr->parent_pcs_ptr->intra_selected_org_qp - prevGopRateControlParamPtr->firstPicActualQp) * 2 / 3;
                        if (prevGopRateControlParamPtr->firstPicActualQp <= 15)
                            qpIncAllowed += 5;
                        else if (prevGopRateControlParamPtr->firstPicActualQp <= 20)
                            qpIncAllowed += 4;
                        else if (prevGopRateControlParamPtr->firstPicActualQp <= 25)
                            qpIncAllowed += 3;
                    }
                    else if (prevGopRateControlParamPtr->scene_change_in_gop) {
                        qpIncAllowed = 5;
                    }
                    if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
                        qpIncAllowed += 2;
                        qpDecAllowed += 4;
                    }
                    picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                        (uint32_t)MAX((int32_t)prevGopRateControlParamPtr->firstPicActualQp - (int32_t)qpDecAllowed, 0),
                        (uint32_t)prevGopRateControlParamPtr->firstPicActualQp + qpIncAllowed,
                        picture_control_set_ptr->picture_qp);
                }

                // Scene change
                if (picture_control_set_ptr->slice_type == I_SLICE && picture_control_set_ptr->picture_number != rateControlParamPtr->firstPoc) {
                    if (nextGopRateControlParamPtr->firstPicActualQpAssigned) {

                        picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                            (uint32_t)MAX((int32_t)nextGopRateControlParamPtr->firstPicActualQp - (int32_t)1, 0),
                            (uint32_t)nextGopRateControlParamPtr->firstPicActualQp + 8,
                            picture_control_set_ptr->picture_qp);
                    }
                    else {
                        if (rateControlParamPtr->firstPicActualQp < 20) {
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                                (uint32_t)MAX((int32_t)rateControlParamPtr->firstPicActualQp - (int32_t)4, 0),
                                (uint32_t)rateControlParamPtr->firstPicActualQp + 10,
                                picture_control_set_ptr->picture_qp);
                        }
                        else {
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                                (uint32_t)MAX((int32_t)rateControlParamPtr->firstPicActualQp - (int32_t)4, 0),
                                (uint32_t)rateControlParamPtr->firstPicActualQp + 8,
                                picture_control_set_ptr->picture_qp);

                        }

                    }
                }
                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    picture_control_set_ptr->picture_qp);
#if NEW_QPS
                picture_control_set_ptr->parent_pcs_ptr->base_qindex = quantizer_to_qindex[picture_control_set_ptr->picture_qp];
#endif
            }

            picture_control_set_ptr->parent_pcs_ptr->picture_qp = picture_control_set_ptr->picture_qp;
            if (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0 && sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
                context_ptr->baseLayerFramesAvgQp = (3 * context_ptr->baseLayerFramesAvgQp + picture_control_set_ptr->picture_qp + 2) >> 2;
            }
            if (picture_control_set_ptr->slice_type == I_SLICE) {
                if (picture_control_set_ptr->picture_number == rateControlParamPtr->firstPoc) {
                    rateControlParamPtr->firstPicPredQp = (uint16_t)picture_control_set_ptr->parent_pcs_ptr->best_pred_qp;
                    rateControlParamPtr->firstPicActualQp = (uint16_t)picture_control_set_ptr->picture_qp;
                    rateControlParamPtr->scene_change_in_gop = picture_control_set_ptr->parent_pcs_ptr->scene_change_in_gop;
                    rateControlParamPtr->firstPicActualQpAssigned = EB_TRUE;
                }
                {
                    if (picture_control_set_ptr->picture_number == rateControlParamPtr->firstPoc) {
                        if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
                            context_ptr->baseLayerIntraFramesAvgQp = (3 * context_ptr->baseLayerIntraFramesAvgQp + picture_control_set_ptr->picture_qp + 2) >> 2;
                        }
                    }

                    if (picture_control_set_ptr->picture_number == rateControlParamPtr->firstPoc) {
                        rateControlParamPtr->intraFramesQp = picture_control_set_ptr->picture_qp;
                        rateControlParamPtr->nextGopIntraFrameQp = picture_control_set_ptr->picture_qp;

                    }
                }
            }
            picture_control_set_ptr->parent_pcs_ptr->average_qp = 0;
            for (lcuCodingOrder = 0; lcuCodingOrder < sequence_control_set_ptr->sb_tot_cnt; ++lcuCodingOrder) {

                sb_ptr = picture_control_set_ptr->sb_ptr_array[lcuCodingOrder];
#if ADD_DELTA_QP_SUPPORT

                sb_ptr->qp = quantizer_to_qindex[picture_control_set_ptr->picture_qp];
#else
                sb_ptr->qp = (uint8_t)picture_control_set_ptr->picture_qp;
#endif
                picture_control_set_ptr->parent_pcs_ptr->average_qp += sb_ptr->qp;
            }

            // Get Empty Rate Control Results Buffer
            EbGetEmptyObject(
                context_ptr->rateControlOutputResultsFifoPtr,
                &rateControlResultsWrapperPtr);
            rateControlResultsPtr = (RateControlResults_t*)rateControlResultsWrapperPtr->objectPtr;
            rateControlResultsPtr->pictureControlSetWrapperPtr = rateControlTasksPtr->pictureControlSetWrapperPtr;

            // Post Full Rate Control Results
            EbPostFullObject(rateControlResultsWrapperPtr);

            // Release Rate Control Tasks
            EbReleaseObject(rateControlTasksWrapperPtr);

            break;

        case RC_PACKETIZATION_FEEDBACK_RESULT:

            parentPictureControlSetPtr = (PictureParentControlSet_t*)rateControlTasksPtr->pictureControlSetWrapperPtr->objectPtr;
            sequence_control_set_ptr = (SequenceControlSet_t*)parentPictureControlSetPtr->sequence_control_set_wrapper_ptr->objectPtr;

            if (sequence_control_set_ptr->static_config.rate_control_mode) {
                rate_control_update_model(rc_model_ptr, parentPictureControlSetPtr);
            }

            // Frame level RC
            if (sequence_control_set_ptr->intra_period_length == -1 || sequence_control_set_ptr->static_config.rate_control_mode == 0) {
                rateControlParamPtr = context_ptr->rateControlParamQueue[0];
                prevGopRateControlParamPtr = context_ptr->rateControlParamQueue[0];
                if (parentPictureControlSetPtr->slice_type == I_SLICE) {

                    if (parentPictureControlSetPtr->total_num_bits > MAX_BITS_PER_FRAME) {
                        context_ptr->maxRateAdjustDeltaQP++;
                    }
                    else if (context_ptr->maxRateAdjustDeltaQP > 0 && parentPictureControlSetPtr->total_num_bits < MAX_BITS_PER_FRAME * 85 / 100) {
                        context_ptr->maxRateAdjustDeltaQP--;
                    }
                    context_ptr->maxRateAdjustDeltaQP = CLIP3(0, 63, context_ptr->maxRateAdjustDeltaQP);
                    context_ptr->maxRateAdjustDeltaQP = 0;
                }
            }
            else {
                uint32_t intervalIndexTemp = 0;
                EbBool intervalFound = EB_FALSE;
                while ((intervalIndexTemp < PARALLEL_GOP_MAX_NUMBER) && !intervalFound) {

                    if (parentPictureControlSetPtr->picture_number >= context_ptr->rateControlParamQueue[intervalIndexTemp]->firstPoc &&
                        parentPictureControlSetPtr->picture_number <= context_ptr->rateControlParamQueue[intervalIndexTemp]->lastPoc) {
                        intervalFound = EB_TRUE;
                    }
                    else {
                        intervalIndexTemp++;
                    }
                }
                CHECK_REPORT_ERROR(
                    intervalIndexTemp != PARALLEL_GOP_MAX_NUMBER,
                    sequence_control_set_ptr->encode_context_ptr->app_callback_ptr,
                    EB_ENC_RC_ERROR2);

                rateControlParamPtr = context_ptr->rateControlParamQueue[intervalIndexTemp];

                prevGopRateControlParamPtr = (intervalIndexTemp == 0) ?
                    context_ptr->rateControlParamQueue[PARALLEL_GOP_MAX_NUMBER - 1] :
                    context_ptr->rateControlParamQueue[intervalIndexTemp - 1];

            }
            if (sequence_control_set_ptr->static_config.rate_control_mode != 0) {

                context_ptr->previousVirtualBufferLevel = context_ptr->virtualBufferLevel;

                context_ptr->virtualBufferLevel =
                    (int64_t)context_ptr->previousVirtualBufferLevel +
                    (int64_t)parentPictureControlSetPtr->total_num_bits - (int64_t)context_ptr->highLevelRateControlPtr->channelBitRatePerFrame;

            }

            // Queue variables
#if OVERSHOOT_STAT_PRINT
            if (sequence_control_set_ptr->intra_period_length != -1) {

                int32_t                       queueEntryIndex;
                uint32_t                       queueEntryIndexTemp;
                uint32_t                       queueEntryIndexTemp2;
                CodedFramesStatsEntry_t     *queueEntryPtr;
                EbBool                      moveSlideWondowFlag = EB_TRUE;
                EbBool                      end_of_sequence_flag = EB_TRUE;
                uint32_t                       frames_in_sw;

                // Determine offset from the Head Ptr
                queueEntryIndex = (int32_t)(parentPictureControlSetPtr->picture_number - context_ptr->codedFramesStatQueue[context_ptr->codedFramesStatQueueHeadIndex]->picture_number);
                queueEntryIndex += context_ptr->codedFramesStatQueueHeadIndex;
                queueEntryIndex = (queueEntryIndex > CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? queueEntryIndex - CODED_FRAMES_STAT_QUEUE_MAX_DEPTH : queueEntryIndex;
                queueEntryPtr = context_ptr->codedFramesStatQueue[queueEntryIndex];

                queueEntryPtr->frameTotalBitActual = (uint64_t)parentPictureControlSetPtr->total_num_bits;
                queueEntryPtr->picture_number = parentPictureControlSetPtr->picture_number;
                queueEntryPtr->end_of_sequence_flag = parentPictureControlSetPtr->end_of_sequence_flag;
                context_ptr->rateAveragePeriodinFrames = (uint64_t)sequence_control_set_ptr->static_config.intra_period_length + 1;

                //printf("\n0_POC: %d\n",
                //    queueEntryPtr->picture_number);
                moveSlideWondowFlag = EB_TRUE;
                while (moveSlideWondowFlag) {
                    //  printf("\n1_POC: %d\n",
                    //      queueEntryPtr->picture_number);
                      // Check if the sliding window condition is valid
                    queueEntryIndexTemp = context_ptr->codedFramesStatQueueHeadIndex;
                    if (context_ptr->codedFramesStatQueue[queueEntryIndexTemp]->frameTotalBitActual != -1) {
                        end_of_sequence_flag = context_ptr->codedFramesStatQueue[queueEntryIndexTemp]->end_of_sequence_flag;
                    }
                    else {
                        end_of_sequence_flag = EB_FALSE;
                    }
                    while (moveSlideWondowFlag && !end_of_sequence_flag &&
                        queueEntryIndexTemp < context_ptr->codedFramesStatQueueHeadIndex + context_ptr->rateAveragePeriodinFrames) {
                        // printf("\n2_POC: %d\n",
                        //     queueEntryPtr->picture_number);

                        queueEntryIndexTemp2 = (queueEntryIndexTemp > CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? queueEntryIndexTemp - CODED_FRAMES_STAT_QUEUE_MAX_DEPTH : queueEntryIndexTemp;

                        moveSlideWondowFlag = (EbBool)(moveSlideWondowFlag && (context_ptr->codedFramesStatQueue[queueEntryIndexTemp2]->frameTotalBitActual != -1));

                        if (context_ptr->codedFramesStatQueue[queueEntryIndexTemp2]->frameTotalBitActual != -1) {
                            // check if it is the last frame. If we have reached the last frame, we would output the buffered frames in the Queue.
                            end_of_sequence_flag = context_ptr->codedFramesStatQueue[queueEntryIndexTemp]->end_of_sequence_flag;
                        }
                        else {
                            end_of_sequence_flag = EB_FALSE;
                        }
                        queueEntryIndexTemp =
                            (queueEntryIndexTemp == CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? 0 : queueEntryIndexTemp + 1;

                    }

                    if (moveSlideWondowFlag) {
                        //get a new entry spot
                        queueEntryPtr = (context_ptr->codedFramesStatQueue[context_ptr->codedFramesStatQueueHeadIndex]);
                        queueEntryIndexTemp = context_ptr->codedFramesStatQueueHeadIndex;
                        // This is set to false, so the last frame would go inside the loop
                        end_of_sequence_flag = EB_FALSE;
                        frames_in_sw = 0;
                        context_ptr->totalBitActualPerSw = 0;

                        while (!end_of_sequence_flag &&
                            queueEntryIndexTemp < context_ptr->codedFramesStatQueueHeadIndex + context_ptr->rateAveragePeriodinFrames) {
                            frames_in_sw++;

                            queueEntryIndexTemp2 = (queueEntryIndexTemp > CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? queueEntryIndexTemp - CODED_FRAMES_STAT_QUEUE_MAX_DEPTH : queueEntryIndexTemp;

                            context_ptr->totalBitActualPerSw += context_ptr->codedFramesStatQueue[queueEntryIndexTemp2]->frameTotalBitActual;
                            end_of_sequence_flag = context_ptr->codedFramesStatQueue[queueEntryIndexTemp2]->end_of_sequence_flag;

                            queueEntryIndexTemp =
                                (queueEntryIndexTemp == CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? 0 : queueEntryIndexTemp + 1;

                        }
                        //

                        //if(frames_in_sw == context_ptr->rateAveragePeriodinFrames)
                        //    printf("POC:%d\t %.3f\n", queueEntryPtr->picture_number, (double)context_ptr->totalBitActualPerSw*(sequence_control_set_ptr->frame_rate>> RC_PRECISION)/(double)frames_in_sw/1000);
                        if (frames_in_sw == (uint32_t)sequence_control_set_ptr->intra_period_length + 1) {
                            context_ptr->maxBitActualPerSw = MAX(context_ptr->maxBitActualPerSw, context_ptr->totalBitActualPerSw*(sequence_control_set_ptr->frame_rate >> RC_PRECISION) / frames_in_sw / 1000);
                            if (queueEntryPtr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0) {
                                context_ptr->maxBitActualPerGop = MAX(context_ptr->maxBitActualPerGop, context_ptr->totalBitActualPerSw*(sequence_control_set_ptr->frame_rate >> RC_PRECISION) / frames_in_sw / 1000);
                                if (context_ptr->totalBitActualPerSw > sequence_control_set_ptr->static_config.maxBufferSize) {
                                    printf("\nPOC:%d\tOvershoot:%.0f%% \n",
                                        (int32_t)queueEntryPtr->picture_number,
                                        (double)((int64_t)context_ptr->totalBitActualPerSw * 100 / (int64_t)sequence_control_set_ptr->static_config.maxBufferSize - 100));
                                }
                            }
                        }
                        if (frames_in_sw == context_ptr->rateAveragePeriodinFrames - 1) {
                            printf("\n%d MAX\n", (int32_t)context_ptr->maxBitActualPerSw);
                            printf("\n%d GopMa\n", (int32_t)context_ptr->maxBitActualPerGop);
                        }

                        // Reset the Queue Entry
                        queueEntryPtr->picture_number += CODED_FRAMES_STAT_QUEUE_MAX_DEPTH;
                        queueEntryPtr->frameTotalBitActual = -1;

                        // Increment the Reorder Queue head Ptr
                        context_ptr->codedFramesStatQueueHeadIndex =
                            (context_ptr->codedFramesStatQueueHeadIndex == CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? 0 : context_ptr->codedFramesStatQueueHeadIndex + 1;

                        queueEntryPtr = (context_ptr->codedFramesStatQueue[context_ptr->codedFramesStatQueueHeadIndex]);

                    }
                }
            }
#endif
            totalNumberOfFbFrames++;

            // Release the SequenceControlSet
            EbReleaseObject(parentPictureControlSetPtr->sequence_control_set_wrapper_ptr);

            // Release the input buffer 
            EbReleaseObject(parentPictureControlSetPtr->input_picture_wrapper_ptr);

            // Release the ParentPictureControlSet
            EbReleaseObject(rateControlTasksPtr->pictureControlSetWrapperPtr);

            // Release Rate Control Tasks
            EbReleaseObject(rateControlTasksWrapperPtr);

            break;

        case RC_ENTROPY_CODING_ROW_FEEDBACK_RESULT:

            // Extract bits-per-lcu-row

            // Release Rate Control Tasks
            EbReleaseObject(rateControlTasksWrapperPtr);

            break;

        default:
            picture_control_set_ptr = (PictureControlSet_t*)rateControlTasksPtr->pictureControlSetWrapperPtr->objectPtr;
            sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
            //encode_context_ptr            = sequence_control_set_ptr->encode_context_ptr;
            //CHECK_REPORT_ERROR_NC(
            //             encode_context_ptr->app_callback_ptr,
            //             EB_ENC_RC_ERROR1);

            break;
        }
    }
    return EB_NULL;
}
