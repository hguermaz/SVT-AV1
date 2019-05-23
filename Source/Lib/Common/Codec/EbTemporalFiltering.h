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
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "EbPictureControlSet.h"
#include "EbMotionEstimationProcess.h"
#include "EbSequenceControlSet.h"
#include "EbDefinitions.h"
#include "EbBitstreamUnit.h"

#define ALTREF_MAX_NFRAMES 10
#if MOVE_TF


EbErrorType init_temporal_filtering(PictureParentControlSet **list_picture_control_set_ptr,
#if FIX_SHORT
    PictureParentControlSet *picture_control_set_ptr_central,
#endif
#if ME_CLEAN
	MotionEstimationContext_t *me_context_ptr,
#endif
	int32_t segment_index);

#else
EbErrorType init_temporal_filtering(PictureParentControlSet **list_picture_control_set_ptr);
#endif