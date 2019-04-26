/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPredictionStructure_h
#define EbPredictionStructure_h

#include "EbDefinitions.h"
#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif
    /************************************************
     * Defines
     ************************************************/
#define FLAT_PREDICTION_STRUCTURE_PERIOD                                    1
#define TWO_LEVEL_HIERARCHICAL_PREDICTION_STRUCTURE_PERIOD                  2
#define THREE_LEVEL_HIERARCHICAL_PREDICTION_STRUCTURE_PERIOD                4
#define FOUR_LEVEL_HIERARCHICAL_PREDICTION_STRUCTURE_PERIOD                 8
#define MAX_PREDICTION_STRUCTURE_PERIOD                                     64

     /************************************************
      * RPS defines
      ************************************************/
#define RPS_UNUSED                                                          ~0
#define MAX_NUM_OF_NEGATIVE_REF_PICS                                        5
#define MAX_NUM_OF_POSITIVE_REF_PICS                                        5
#define MAX_NUM_OF_REF_PICS_TOTAL                                           (MAX_NUM_OF_POSITIVE_REF_PICS + MAX_NUM_OF_NEGATIVE_REF_PICS)

      /************************************************
       * Reference List
       *
       *   reference_list - Contains the deltaPOCs of
       *    the pictures referenced by the current
       *    picture.
       ************************************************/
    typedef struct ReferenceList
    {
#if MRP_ME
        int32_t                              *reference_list;
#else
        int32_t                              reference_list;
#endif
        uint32_t                              reference_list_count;

    } ReferenceList;

    /************************************************
     * Dependent List
     *
     *   list_count - Contains count of how
     *     deep into list should be used
     *     depending on how many references are
     *     being used in the prediction structure.
     *
     *   list - Contains the deltaPOCs of
     *     pictures that reference the current picture.
     *     The dependent list pictures must be grouped
     *     by the referenceCount group in ascending
     *     order.  The grouping is not display order!
     ************************************************/
    typedef struct DependentList
    {
        int32_t                             *list;
        uint32_t                              list_count;

    } DependentList;

    /************************************************
     * Prediction Structure Config Entry
     *   Contains the basic reference lists and
     *   configurations for each Prediction Structure
     *   Config Entry.
     ************************************************/
    typedef struct PredictionStructureConfigEntry 
    {
        uint32_t temporal_layer_index;
        uint32_t decode_order;
#if MRP_ME
        int32_t                              ref_list0[REF_LIST_MAX_DEPTH];
        int32_t                              ref_list1[REF_LIST_MAX_DEPTH];
#else
        int32_t  ref_list0;
        int32_t  ref_list1;
#endif
    } PredictionStructureConfigEntry;

    /************************************************
     * Prediction Structure Config
     *   Contains a collection of basic control data
     *   for the basic prediction structure.
     ************************************************/
    typedef struct PredictionStructureConfig 
    {
        uint32_t                          entry_count;
        PredictionStructureConfigEntry   *entry_array;
    } PredictionStructureConfig;

    /************************************************
     * Prediction Structure Entry
     *   Contains the reference and dependent lists
     *   for a particular picture in the Prediction
     *   Structure.
     ************************************************/
    typedef struct PredictionStructureEntry 
    {
        ReferenceList                     ref_list0;
        ReferenceList                     ref_list1;
        DependentList                     dep_list0;
        DependentList                     dep_list1;
        uint32_t                           temporal_layer_index;
        uint32_t                           decode_order;
        EbBool                             is_referenced;

        // High-level RPS
        EbBool                             short_term_rps_in_sps_flag;
        uint32_t                              short_term_rps_in_sps_index;
        EbBool                             inter_rps_prediction_flag;
        EbBool                             long_term_rps_present_flag;
        uint32_t                              gop_position_least_significant_bits;

        // Predicted Short-Term RPS
        uint32_t                              delta_rps_index_minus1;
        uint32_t                              absolute_delta_rps_minus1;
        uint32_t                              delta_rps_sign;
        EbBool                             used_by_curr_pic_flag[MAX_NUM_OF_REF_PICS_TOTAL];
        EbBool                             used_by_future_pic_flag[MAX_NUM_OF_REF_PICS_TOTAL];
        // Non-Predicted Short-Term RPS
        uint32_t                              negative_ref_pics_total_count;
        uint32_t                              positive_ref_pics_total_count;
        uint32_t                              delta_negative_gop_pos_minus1[MAX_NUM_OF_NEGATIVE_REF_PICS];
        uint32_t                              delta_positive_gop_pos_minus1[MAX_NUM_OF_POSITIVE_REF_PICS];
        EbBool                             used_by_negative_curr_pic_flag[MAX_NUM_OF_NEGATIVE_REF_PICS];
        EbBool                             used_by_positive_curr_pic_flag[MAX_NUM_OF_POSITIVE_REF_PICS];

        // Long-Term RPS
        uint32_t                              long_term_ref_pics_total_count;
        uint32_t                              delta_gop_poslsb[MAX_NUM_OF_REF_PICS_TOTAL];
        EbBool                             delta_gop_pos_msb_present_flag[MAX_NUM_OF_REF_PICS_TOTAL];
        uint32_t                              delta_gop_pos_msb_minus1[MAX_NUM_OF_REF_PICS_TOTAL];
        EbBool                             used_by_lt_curr_pic_flag_array[MAX_NUM_OF_REF_PICS_TOTAL];
        // List Construction
        EbBool                             ref_pics_override_total_count_flag;
        int32_t                              ref_pics_list0_total_count_minus1;
        int32_t                              ref_pics_list1_total_count_minus1;
        EbBool                             lists_modification_present_flag;
        EbBool                             restricted_ref_pic_lists_flag;      // Same list enable flag (if set,
                                                                            //   it implies all slices of the
                                                                            //   same type in the same picture
                                                                            //   have identical lists)

        // List Modification
        // *Note - This should probably be moved to the slice header since its a dynamic control - JMJ Jan 2, 2013
        EbBool                             list0_modification_flag;
        EbBool                             list1_modification_flag;
        uint32_t                              list0_mod_index[MAX_NUM_OF_REF_PICS_TOTAL];
        uint32_t                              list1_mod_index[MAX_NUM_OF_REF_PICS_TOTAL];

        // Lists Combination (STUB)

    } PredictionStructureEntry;

    /************************************************
     * Prediction Structure
     *   Contains a collection of control and RPS
     *   data types for an entire Prediction Structure
     ************************************************/
    typedef struct PredictionStructure 
    {

        uint32_t                              pred_struct_entry_count;
        PredictionStructureEntry        **pred_struct_entry_ptr_array;
        EbPred                             pred_type;
        uint32_t                              temporal_layer_count;
        uint32_t                              pred_struct_period;
        uint32_t                              maximum_extent;

        // Section Indices
        uint32_t                              leading_pic_index;
        uint32_t                              init_pic_index;
        uint32_t                              steady_state_index;

        // RPS Related Entries
        EbBool                             restricted_ref_pic_lists_enable_flag;
        EbBool                             lists_modification_enable_flag;
        EbBool                             long_term_enable_flag;
        uint32_t                              default_ref_pics_list0_total_count_minus1;
        uint32_t                              default_ref_pics_list1_total_count_minus1;

    } PredictionStructure;

    /************************************************
     * Prediction Structure Group
     *   Contains the control structures for all
     *   supported prediction structures.
     ************************************************/
    typedef struct PredictionStructureGroup 
    {
        PredictionStructure             **prediction_structure_ptr_array;
        uint32_t                              prediction_structure_count;
    } PredictionStructureGroup;

    /************************************************
     * Declarations
     ************************************************/
    extern EbErrorType prediction_structure_group_ctor(
#if MRP_M1
        uint8_t          enc_mode,
#endif
        PredictionStructureGroup   **predictionStructureGroupDblPtr,
        uint32_t                        base_layer_switch_mode);

    extern PredictionStructure* get_prediction_structure(
        PredictionStructureGroup    *prediction_structure_group_ptr,
        EbPred                        pred_structure,
        uint32_t                         number_of_references,
        uint32_t                         levels_of_hierarchy);
#if NEW_RPS
    enum {
        LAST = 0,
        LAST2 = 1,
        LAST3 = 2,
        GOLD = 3,
        BWD = 4,
        ALT2 = 5,
        ALT = 6
    } REF_FRAME_MINUS1;
#endif
    typedef struct Av1RpsNode 
    {
        uint8_t refresh_frame_mask;
        uint8_t ref_dpb_index[7];//LAST-LAST2-LAST3-GOLDEN-BWD-ALT2-ALT
#if REF_ORDER
        uint64_t ref_poc_array[7]; //decoder based ref poc array //LAST-LAST2-LAST3-GOLDEN-BWD-ALT2-ALT
#endif
    } Av1RpsNode;

#if MRP_PRED_STRUCTURE
    /**********************************************************
    * Macros
    **********************************************************/
#define PRED_STRUCT_INDEX(hierarchicalLevelCount, predType, refCount) (((hierarchicalLevelCount) * EB_PRED_TOTAL_COUNT + (predType)) * REF_LIST_MAX_DEPTH + (refCount))

    //#define DEP_INDEX(predIdx, entryIndex, entryTotalCount) ((((predIdx) - ((int32_t) entryIndex)) % (entryTotalCount)))
#define DEP_INDEX(predIdx, entryIndex, entryTotalCount) (((((int32_t) entryIndex) - (predIdx)) % (entryTotalCount)))

    /**********************************************************
    * Instructions for how to create a Predicion Structure
    *
    * Overview:
    *   The prediction structure consists of a collection
    *   of Prediction Structure Entires, which themselves
    *   consist of reference and dependent lists.  The
    *   reference lists are exactly like those found in
    *   the standard and can be clipped in order to reduce
    *   the number of references.
    *
    *   Dependent lists are the corollary to reference lists,
    *   the describe how a particular picture is referenced.
    *   Dependent lists can also be clipped at predefined
    *   junctions (i.e. the list_count array) in order
    *   to reduce the number of references.  Note that the
    *   dependent deltaPOCs must be grouped together in order
    *   of ascending referencePicture in order for the Dependent
    *   List clip to work properly.
    *
    *   All control and RPS information is derived from
    *   these lists.  The lists for a structure are defined
    *   for both P & B-picture variants.  In the case of
    *   P-pictures, only Lists 0 are used.
    *
    *   Negative deltaPOCs are for backward-referencing pictures
    *   in display order and positive deltaPOCs are for
    *   forward-referencing pictures.
    *
    *   Please note that there is no assigned coding order,
    *   the PictureManager will start pictures as soon as
    *   their references become available.
    *
    *   Any prediction structure is possible; however, we are
    *     restricting usage to the following controls:
    *     # Hierarchical Levels
    *     # Number of References
    *     # B-pictures enabled
    *     # Intra Refresh Period
    *
    *  To Get Low Delay P, only use List 0
    *  To Get Low Delay B, replace List 1 with List 0
    *  To Get Random Access, use the preduction structure as is
    **********************************************************/

    /************************************************
    * Flat
    *
    *  I-B-B-B-B-B-B-B-B
    *
    * Display & Coding Order:
    *  0 1 2 3 4 5 6 7 8
    *
    ************************************************/
    static PredictionStructureConfigEntry flat_pred_struct[] = {

        {
            0,               // GOP Index 0 - Temporal Layer
            0,               // GOP Index 0 - Decode Order
    { 1, 2, 3, 4 },    // GOP Index 0 - Ref List 0
    { 1, 2, 3, 4 }     // GOP Index 0 - Ref List 1
        }
    };

    /************************************************
    * Random Access - Two-Level Hierarchical
    *
    *    b   b   b   b      Temporal Layer 1
    *   / \ / \ / \ / \
    *  I---B---B---B---B    Temporal Layer 0
    *
    * Display Order:
    *  0 1 2 3 4 5 6 7 8
    *
    * Coding Order:
    *  0 2 1 4 3 6 5 8 7
    ************************************************/
    static PredictionStructureConfigEntry two_level_hierarchical_pred_struct[] = {

        {
            0,                // GOP Index 0 - Temporal Layer
            0,                // GOP Index 0 - Decode Order
    { 2, 4, 6, 8 },     // GOP Index 0 - Ref List 0
    { 2, 4, 6, 8 }      // GOP Index 0 - Ref List 1
        },

    {
        1,                // GOP Index 1 - Temporal Layer
        1,                // GOP Index 1 - Decode Order
    { 1, 3, 5, 7 },    // GOP Index 1 - Ref List 0
    { -1, 1, 3, 5 }     // GOP Index 1 - Ref List 1
    }
    };

    /************************************************
    * Three-Level Hierarchical
    *
    *      b   b       b   b       b   b        Temporal Layer 2
    *     / \ / \     / \ / \     / \ / \
    *    /   B   \   /   B   \   /   B   \      Temporal Layer 1
    *   /   / \   \ /   / \   \ /   / \   \
    *  I-----------B-----------B-----------B    Temporal Layer 0
    *
    * Display Order:
    *  0   1 2 3   4   5 6 7   8   9 1 1   1
    *                                0 1   2
    *
    * Coding Order:
    *  0   3 2 4   1   7 6 8   5   1 1 1   9
    *                              1 0 2
    ************************************************/
    static PredictionStructureConfigEntry three_level_hierarchical_pred_struct[] = {

        {
            0,                  // GOP Index 0 - Temporal Layer
            0,                  // GOP Index 0 - Decode Order
    { 4, 8, 12, 16 },     // GOP Index 0 - Ref List 0
    { 4, 8, 12, 16 }      // GOP Index 0 - Ref List 1
        },

    {
        2,                  // GOP Index 1 - Temporal Layer
        2,                  // GOP Index 1 - Decode Order
    { 1,  3, 5, 0 },     // GOP Index 1 - Ref List 0
    { -1, -3, 1, 3 }      // GOP Index 1 - Ref List 1
    },

    {
        1,                  // GOP Index 2 - Temporal Layer
        1,                  // GOP Index 2 - Decode Order
    { 2, 6, 10, 14 },    // GOP Index 2 - Ref List 0
    { -2, 2,  6, 10 }     // GOP Index 2 - Ref List 1
    },

    {
        2,                  // GOP Index 3 - Temporal Layer
        3,                  // GOP Index 3 - Decode Order
    { 1, 3, 5, 7 },      // GOP Index 3 - Ref List 0
    { -1, 1, 3, 5 }       // GOP Index 3 - Ref List 1
    }
    };

    /************************************************************************************************************
    * Four-Level Hierarchical
    *
    *
    *          b     b           b     b               b     b           b     b           Temporal Layer 3
    *         / \   / \         / \   / \             / \   / \         / \   / \
    *        /   \ /   \       /   \ /   \           /   \ /   \       /   \ /   \
    *       /     B     \     /     B     \         /     B     \     /     B     \        Temporal Layer 2
    *      /     / \     \   /     / \     \       /     / \     \   /     / \     \
    *     /     /   \     \ /     /   \     \     /     /   \     \ /     /   \     \
    *    /     /     ------B------     \     \   /     /     ------B------     \     \     Temporal Layer 1
    *   /     /           / \           \     \ /     /           / \           \     \
    *  I---------------------------------------B---------------------------------------B   Temporal Layer 0
    *
    * Display Order:
    *  0       1  2  3     4     5  6  7       8       9  1  1     1     1  1  1       1
    *                                                     0  1     2     3  4  5       6
    *
    * Coding Order:
    *  0       5  3  6     2     7  4  8       1       1  1  1     1     1  1  1       9
    *                                                  3  1  4     0     5  2  6
    *
    ***********************************************************************************************************/
#if RPS_4L
    static PredictionStructureConfigEntry four_level_hierarchical_pred_struct[] = {
#else
    static PredictionStructureConfigEntry four_level_hierarchical_pred_struct[] = {
#endif

        {
            0,                  // GOP Index 0 - Temporal Layer
            0,                  // GOP Index 0 - Decode Order
#if RPS_4L
    { 8, 0,  0, 0 },     // GOP Index 0 - Ref List 0
    { 8, 0,  0, 0 }      // GOP Index 0 - Ref List 1
#else
    {8, 16,  0, 0},     // GOP Index 0 - Ref List 0
    { 8, 16,  0, 0 }      // GOP Index 0 - Ref List 1
#endif
        },

    {
        3,                  // GOP Index 1 - Temporal Layer
        3,                  // GOP Index 1 - Decode Order
#if RPS_4L
    { 1,  3,  5,  9 },   // GOP Index 1 - Ref List 0
    { -1, -3, -7,  0 }    // GOP Index 1 - Ref List 1
#else
    { 1,  3,  0,  0},   // GOP Index 1 - Ref List 0
    { -1, -3, -7,  1 }    // GOP Index 1 - Ref List 1
#endif
    },

    {
        2,                  // GOP Index 2 - Temporal Layer
        2,                  // GOP Index 2 - Decode Order
#if RPS_4L
    { 2,  4,  6,  10 },   // GOP Index 2 - Ref List 0
    { -2, -6,  0,   0 }    // GOP Index 2 - Ref List 1
#else
    {2,   6, 10,  0},   // GOP Index 2 - Ref List 0
    { -2, -6,  2,  6 }    // GOP Index 2 - Ref List 1
#endif
    },

    {
        3,                  // GOP Index 3 - Temporal Layer
        4,                  // GOP Index 3 - Decode Order
#if RPS_4L
    { 1,   3, 5,  7 },    //    GOP Index 3 - Ref List 0  
    { -1,  -5, 0,  0 }     //     GOP Index 3 - Ref List 1
#else
    { 1,  3,  5, 0},    // GOP Index 3 - Ref List 0
    { -1, -5,  1, 3 }     // GOP Index 3 - Ref List 1
#endif
    },

    {
        1,                  // GOP Index 4 - Temporal Layer
        1,                  // GOP Index 4 - Decode Order
#if RPS_4L
    { 4, 8, 12,  0 },    // GOP Index 4 - Ref List 0   
    { -4,   0, 0,  0 }   // GOP Index 4 - Ref List 1
#else
    { 4, 12,  0,  0},   // GOP Index 4 - Ref List 0
    { -4,  4, 12,  0 }    // GOP Index 4 - Ref List 1
#endif
    },

    {
        3,                  // GOP Index 5 - Temporal Layer
        6,                  // GOP Index 5 - Decode Order
#if RPS_4L
    { 1,   3, 5,  9 },    // GOP Index 5 - Ref List 0
    { -1,  -3, 0,  0 }     // GOP Index 5 - Ref List 1
#else
    { 1,  3, 5, 0},     // GOP Index 5 - Ref List 0
    { -1, -3, 1, 3 }      // GOP Index 5 - Ref List 1
#endif
    },

    {
        2,                  // GOP Index 6 - Temporal Layer
        5,                  // GOP Index 6 - Decode Order
#if RPS_4L
    { 2,  4, 6,  10 },   // GOP Index 6 - Ref List 0
    { -2,  0, 0,  0 }    // GOP Index 6 - Ref List 1
#else
    { 2, 6, 10,  0},    // GOP Index 6 - Ref List 0
    { -2, 2,  6, 10 }     // GOP Index 6 - Ref List 1
#endif
    },

    {
        3,                  // GOP Index 7 - Temporal Layer
        7,                  // GOP Index 7 - Decode Order
#if RPS_4L
    { 1,  3, 5,  7 },     //  GOP Index 7 - Ref List 0
    { -1,  0, 0,  0 }      // GOP Index 7 - Ref List 1
#else
    { 1, 3, 5, 7},      // GOP Index 7 - Ref List 0
    { -1, 1, 3, 5 }       // GOP Index 7 - Ref List 1
#endif
    }
    };

    /***********************************************************************************************************
    * Five-Level Level Hierarchical
    *
    *           b     b           b     b               b     b           b     b              Temporal Layer 4
    *          / \   / \         / \   / \             / \   / \         / \   / \
    *         /   \ /   \       /   \ /   \           /   \ /   \       /   \ /   \
    *        /     B     \     /     B     \         /     B     \     /     B     \           Temporal Layer 3
    *       /     / \     \   /     / \     \       /     / \     \   /     / \     \
    *      /     /   \     \ /     /   \     \     /     /   \     \ /     /   \     \
    *     /     /     ------B------     \     \   /     /     ------B------     \     \        Temporal Layer 2
    *    /     /           / \           \     \ /     /           / \           \     \
    *   /     /           /   \-----------------B------------------   \           \     \      Temporal Layer 1
    *  /     /           /                     / \                     \           \     \
    * I-----------------------------------------------------------------------------------B    Temporal Layer 0
    *
    * Display Order:
    *  0        1  2  3     4     5  6  7       8       9  1  1     1     1  1  1         1
    *                                                      0  1     2     3  4  5         6
    *
    * Coding Order:
    *  0        9  5  1     3     1  6  1       2       1  7  1     4     1  8  1         1
    *                 0           1     2               3     4           5     6
    *
    ***********************************************************************************************************/
#if REF_ORDER
    static PredictionStructureConfigEntry five_level_hierarchical_pred_struct[] = {
#else
    static PredictionStructureConfigEntry five_level_hierarchical_pred_struct[] = {
#endif

        {
            0,                  // GOP Index 0 - Temporal Layer
            0,                  // GOP Index 0 - Decode Order
#if MRP_BASE
    { 16, 48, 0, 0 },      // GOP Index 0 - Ref List 0
    { 16, 32, 0, 0 }       // GOP Index 0 - Ref List 1     //we need keep 16 as first entry in List1, this will ensure that for poc=16 there is a valid ref frame in List1.
#else
    {16, 0, 0, 0},      // GOP Index 0 - Ref List 0
    { 16, 0, 0, 0 }       // GOP Index 0 - Ref List 1
#endif
        },

    {
        4,                  // GOP Index 1 - Temporal Layer
        4,                  // GOP Index 1 - Decode Order

#if MRP_5L_STRUCT  
    { 1, 5, 9, 17 },  // GOP Index 1 - Ref List 0
    { -1, -3, -7, 0 }   // GOP Index 1 - Ref List 1
#else
    { 1,  0,  0,   0},  // GOP Index 1 - Ref List 0
    { -1, -3, -7, -15 }   // GOP Index 1 - Ref List 1
#endif
    },

    {
        3,                  // GOP Index 2 - Temporal Layer
        3,                  // GOP Index 2 - Decode Order

#if MRP_5L_STRUCT  
    { 2, 4, 6, 10 },        // GOP Index 2 - Ref List 0
    { -2, -6, -14,  0 }   // GOP Index 2 - Ref List 1
#else
    { 2,  0,   0,  0},  // GOP Index 2 - Ref List 0
    { -2, -6, -14,  0 }   // GOP Index 2 - Ref List 1
#endif
    },

    {
        4,                  // GOP Index 3 - Temporal Layer
        5,                  // GOP Index 3 - Decode Order
#if MRP_5L_STRUCT  
    { 1, 3, 7, 11 },     // GOP Index 3 - Ref List 0
    { -1, -5, -13, 0 }   // GOP Index 3 - Ref List 1
#else
    { 1,  3,   0, 0},   // GOP Index 3 - Ref List 0
    { -1, -5, -13, 1 }    // GOP Index 3 - Ref List 1
#endif
    },

    {
        2,                  // GOP Index 4 - Temporal Layer
        2,                  // GOP Index 4 - Decode Order
#if MRP_5L_STRUCT  
    { 4,   8,  12,  20 },  // GOP Index 4 - Ref List 0
    { -4, -12,  0,  0 }     // GOP Index 4 - Ref List 1
#else
    { 4,   0,  0,  0},  // GOP Index 4 - Ref List 0
    { -4, -12,  4,  0 }   // GOP Index 4 - Ref List 1
#endif
    },

    {
        4,                  // GOP Index 5 - Temporal Layer
        7,                  // GOP Index 5 - Decode Order
#if MRP_5L_STRUCT  
    { 1, 5, 9, 13 },      // GOP Index 5 - Ref List 0
    { -1, -3, -11, 0 }   // GOP Index 5 - Ref List 1
#else
    { 1,  5,   0, 0},   // GOP Index 5 - Ref List 0
    { -1, -3, -11, 1 }    // GOP Index 5 - Ref List 1
#endif
    },

    {
        3,                  // GOP Index 6 - Temporal Layer
        6,                  // GOP Index 6 - Decode Order
#if MRP_5L_STRUCT  
    { 2, 4, 6, 10 },        // GOP Index 6 - Ref List 0
    { -2, -10,  0,  0 }    // GOP Index 6 - Ref List 1
#else
    { 2,   6,  0,  0},  // GOP Index 6 - Ref List 0
    { -2, -10,  2,  6 }   // GOP Index 6 - Ref List 1
#endif
    },

    {
        4,                  // GOP Index 7 - Temporal Layer
        8,                  // GOP Index 7 - Decode Order
#if MRP_5L_STRUCT  
    { 1, 3, 7, 11 },    // GOP Index 7 - Ref List 0
    { -1, -9, 0, 0 }   // GOP Index 7 - Ref List 1
#else
    { 1,  3, 7, 0},     // GOP Index 7 - Ref List 0
    { -1, -9, 1, 3 }      // GOP Index 7 - Ref List 1
#endif
    },

    {
        1,                  // GOP Index 8 - Temporal Layer
        1,                  // GOP Index 8 - Decode Order
#if MRP_5L_STRUCT  
    { 8, 16, 24, 0 },   // GOP Index 8 - Ref List 0
    { -8, 0, 0, 0 }      // GOP Index 8 - Ref List 1
#else
    { 8,  0,  0,  0},   // GOP Index 8 - Ref List 0
    { -8,  8,  0,  0 }    // GOP Index 8 - Ref List 1
#endif
    },

    {
        4,                  // GOP Index 9 - Temporal Layer
        11,                 // GOP Index 9 - Decode Order
#if MRP_5L_STRUCT  
    { 1, 5, 9, 17 },     // GOP Index 9 - Ref List 0
    { -1, -3, -7, 0 }   // GOP Index 9 - Ref List 1
#else

    { 1,  9,  0,  0},   // GOP Index 9 - Ref List 0
    { -1, -3, -7,  1 }    // GOP Index 9 - Ref List 1
#endif
    },

    {
        3,                  // GOP Index 10 - Temporal Layer
        10,                 // GOP Index 10 - Decode Order
#if MRP_5L_STRUCT  
    { 2, 4, 6, 10 },       // GOP Index 10 - Ref List 0
    { -2, -6,  0,  0 }    // GOP Index 10 - Ref List 1
#else
    { 2, 10,  0,  0},   // GOP Index 10 - Ref List 0
    { -2, -6,  2,  0 }    // GOP Index 10 - Ref List 1
#endif
    },

    {
        4,                  // GOP Index 11 - Temporal Layer
        12,                 // GOP Index 11 - Decode Order
#if MRP_5L_STRUCT  
    { 1, 3, 7, 11 },    // GOP Index 11 - Ref List 0
    { -1, -5, 0, 0 }   // GOP Index 11 - Ref List 1
#else
    { 1,  3, 11, 0},    // GOP Index 11 - Ref List 0
    { -1, -5,  1, 3 }     // GOP Index 11 - Ref List 1
#endif
    },

    {
        2,                  // GOP Index 12 - Temporal Layer
        9,                  // GOP Index 12 - Decode Order
#if MRP_5L_STRUCT    
    { 4, 8, 12, 0 },       // GOP Index 12 - Ref List 0
    { -4,  0, 0,  0 }     // GOP Index 12 - Ref List 1
#else
    { 4, 12,  0,  0},   // GOP Index 12 - Ref List 0
    { -4,  4, 12,  0 }    // GOP Index 12 - Ref List 1
#endif
    },

    {
        4,                  // GOP Index 13 - Temporal Layer
        14,                 // GOP Index 13 - Decode Order
#if MRP_5L_STRUCT  
    { 1, 5, 9, 13 },  // GOP Index 13 - Ref List 0
    { -1, -3, 0, 0 }   // GOP Index 13 - Ref List 1
#else
    { 1,  3, 13, 0},    // GOP Index 13 - Ref List 0
    { -1, -3,  1, 3 }     // GOP Index 13 - Ref List 1
#endif
    },

    {
        3,                  // GOP Index 14 - Temporal Layer
        13,                 // GOP Index 14 - Decode Order
#if MRP_5L_STRUCT  
    { 2, 4, 6, 14 },    // GOP Index 14 - Ref List 0
    { -2, 0,  0, 0 }   // GOP Index 14 - Ref List 1

#else
    { 2, 6, 10, 14},    // GOP Index 14 - Ref List 0
    { -2, 2,  6, 10 }     // GOP Index 14 - Ref List 1
#endif
    },

    {
        4,                  // GOP Index 15 - Temporal Layer
        15,                 // GOP Index 15 - Decode Order
#if MRP_5L_STRUCT  
    { 1, 3, 7, 11 },  // GOP Index 15 - Ref List 0
    { -1, 0, 0, 0 }   // GOP Index 15 - Ref List 1
#else
    { 1, 3, 7, 0},      // GOP Index 15 - Ref List 0
    { -1, 1, 3, 7 }       // GOP Index 15 - Ref List 1
#endif
    }
    };

    /**********************************************************************************************************************************************************************************************************************
    * Six-Level Level Hierarchical
    *
    *
    *              b     b           b     b               b     b           b     b                   b     b           b     b               b     b           b     b               Temporal Layer 5
    *             / \   / \         / \   / \             / \   / \         / \   / \                 / \   / \         / \   / \             / \   / \         / \   / \
    *            /   \ /   \       /   \ /   \           /   \ /   \       /   \ /   \               /   \ /   \       /   \ /   \           /   \ /   \       /   \ /   \
    *           /     B     \     /     B     \         /     B     \     /     B     \             /     B     \     /     B     \         /     B     \     /     B     \            Temporal Layer 4
    *          /     / \     \   /     / \     \       /     / \     \   /     / \     \           /     / \     \   /     / \     \       /     / \     \   /     / \     \
    *         /     /   \     \ /     /   \     \     /     /   \     \ /     /   \     \         /     /   \     \ /     /   \     \     /     /   \     \ /     /   \     \
    *        /     /     ------B------     \     \   /     /     ------B------     \     \       /     /     ------B------     \     \   /     /     ------B------     \     \         Temporal Layer 3
    *       /     /           / \           \     \ /     /           / \           \     \     /     /           / \           \     \ /     /           / \           \     \
    *      /     /           /   \-----------------B------------------   \           \     \   /     /           /   \-----------------B------------------   \           \     \       Temporal Layer 2
    *     /     /           /                     / \                     \           \     \ /     /           /                     / \                     \           \     \
    *    /     /           /                     /   \---------------------------------------B---------------------------------------/   \                     \           \     \     Temporal Layer 1
    *   /     /           /                     /                                           / \                                           \                     \           \     \
    *  I---------------------------------------------------------------------------------------------------------------------------------------------------------------------------B   Temporal Layer 0
    *
    * Display Order:
    *  0           1  2  3     4     5  6  7       8       9  1  1     1     1  1  1         1         1  1  1     2     2  2  2       2       2  2  2     2     2  3  3           3
    *                                                         0  1     2     3  4  5         6         7  8  9     0     1  2  3       4       5  6  7     8     9  0  1           2
    *
    * Coding Order:
    *  0           1  9  1     5     1  1  2       3       2  1  2     6     2  1  2         2         2  1  2     7     2  1  2       4       2  1  3     8     3  1  3           1
    *              7     8           9  0  0               1  1  2           3  2  4                   5  3  6           7  4  8               9  5  0           1  6  2
    *
    **********************************************************************************************************************************************************************************************************************/
    static PredictionStructureConfigEntry six_level_hierarchical_pred_struct[] = {

        {
            0,                  // GOP Index 0 - Temporal Layer
            0,                  // GOP Index 0 - Decode Order
    { 32,  0, 0, 0 },     // GOP Index 0 - Ref List 0
    { 32,  0, 0, 0 }      // GOP Index 0 - Ref List 1
        },

    {
        5,                  // GOP Index 1 - Temporal Layer
        5,                  // GOP Index 1 - Decode Order
    { 1,  0,  0,   0 },  // GOP Index 1 - Ref List 0
    { -1, -3, -7, -15 }   // GOP Index 1 - Ref List 1
    },

    {
        4,                  // GOP Index 2 - Temporal Layer
        4,                  // GOP Index 2 - Decode Order
    { 2,   0,   0,   0 }, // GOP Index 2 - Ref List 0
    { -2, -6, -14, -30 }  // GOP Index 2 - Ref List 1
    },

    {
        5,                  // GOP Index 3 - Temporal Layer
        6,                  // GOP Index 3 - Decode Order
    { 1,  3,   0,   0 }, // GOP Index 3 - Ref List 0
    { -1, -5, -13,   0 }  // GOP Index 3 - Ref List 1
    },

    {
        3,                  // GOP Index 4 - Temporal Layer
        3,                  // GOP Index 4 - Decode Order
    { 4,   0,   0,  0 }, // GOP Index 4 - Ref List 0
    { -4, -12, -28,  4 }  // GOP Index 4 - Ref List 1
    },

    {
        5,                  // GOP Index 5 - Temporal Layer
        8,                  // GOP Index 5 - Decode Order
    { 1,  5,   0, 0 },   // GOP Index 5 - Ref List 0
    { -1, -3, -11, 0 }    // GOP Index 5 - Ref List 1
    },

    {
        4,                  // GOP Index 6 - Temporal Layer
        7,                  // GOP Index 6 - Decode Order
    { 2,   6,   0, 0 },  // GOP Index 6 - Ref List 0
    { -2, -10, -26, 2 }   // GOP Index 6 - Ref List 1
    },

    {
        5,                  // GOP Index 7 - Temporal Layer
        9,                  // GOP Index 7 - Decode Order
    { 1,  3,   0, 0 },   // GOP Index 7 - Ref List 0
    { -1, -9, -25, 1 }    // GOP Index 7 - Ref List 1
    },

    {
        2,                  // GOP Index 8 - Temporal Layer
        2,                  // GOP Index 8 - Decode Order
    { 8,   0, 0, 0 },    // GOP Index 8 - Ref List 0
    { -8, -24, 8, 0 }     // GOP Index 8 - Ref List 1
    },

    {
        5,                  // GOP Index 9 - Temporal Layer
        12,                 // GOP Index 9 - Decode Order
    { 1,  9,  0, 0 },    // GOP Index 9 - Ref List 0
    { -1, -3, -7, 0 }     // GOP Index 9 - Ref List 1
    },

    {
        4,                  // GOP Index 10 - Temporal Layer
        11,                 // GOP Index 10 - Decode Order
    { 2,  10,   0,  0 },  // GOP Index 10 - Ref List 0
    { -2, -6, -22,  2 }   // GOP Index 10 - Ref List 1
    },

    {
        5,                  // GOP Index 11 - Temporal Layer
        13,                 // GOP Index 11 - Decode Order
    { 1,  3,   0, 0 },   // GOP Index 11 - Ref List 0
    { -1, -5, -21, 1 }    // GOP Index 11 - Ref List 1
    },

    {
        3,                  // GOP Index 12 - Temporal Layer
        10,                 // GOP Index 12 - Decode Order
    { 4,  12,  0,  0 },  // GOP Index 12 - Ref List 0
    { -4, -20,  4, 12 }   // GOP Index 12 - Ref List 1
    },

    {
        5,                  // GOP Index 13 - Temporal Layer
        15,                 // GOP Index 13 - Decode Order
    { 1,  5,   0, 0 },  // GOP Index 13 - Ref List 0
    { -1, -3, -19, 1 }   // GOP Index 13 - Ref List 1
    },

    {
        4,                  // GOP Index 14 - Temporal Layer
        14,                 // GOP Index 14 - Decode Order
    { 2,   6, 14,  0 },  // GOP Index 14 - Ref List 0
    { -2, -18,  2,  6 }   // GOP Index 14 - Ref List 1
    },

    {
        5,                  // GOP Index 15 - Temporal Layer
        16,                 // GOP Index 15 - Decode Order
    { 1,   3, 7,  0 },   // GOP Index 15 - Ref List 0
    { -1, -17, 1,  3 }    // GOP Index 15 - Ref List 1
    },

    {
        1,                  // GOP Index 16 - Temporal Layer
        1,                  // GOP Index 16 - Decode Order
    { 16,  0, 0, 0 },    // GOP Index 16 - Ref List 0
    { -16, 16, 0, 0 }     // GOP Index 16 - Ref List 1
    },

    {
        5,                  // GOP Index 17 - Temporal Layer
        20,                 // GOP Index 17 - Decode Order
    { 1, 17,  0,  0 },   // GOP Index 17 - Ref List 0
    { -1, -3, -7,  0 }    // GOP Index 17 - Ref List 1
    },

    {
        4,                  // GOP Index 18 - Temporal Layer
        19,                 // GOP Index 18 - Decode Order
    { 2, 18,   0,  0 },  // GOP Index 18 - Ref List 0
    { -2, -6, -14,  2 }   // GOP Index 18 - Ref List 1
    },

    {
        5,                  // GOP Index 19 - Temporal Layer
        21,                 // GOP Index 19 - Decode Order
    { 1,  3,   0, 0 },   // GOP Index 19 - Ref List 0
    { -1, -5, -13, 1 }    // GOP Index 19 - Ref List 1
    },

    {
        3,                  // GOP Index 20 - Temporal Layer
        18,                 // GOP Index 20 - Decode Order
    { 4,  20, 0,  0 },   // GOP Index 20 - Ref List 0
    { -4, -12, 4, 20 }    // GOP Index 20 - Ref List 1
    },

    {
        5,                  // GOP Index 21 - Temporal Layer
        23,                 // GOP Index 21 - Decode Order
    { 1,  5,   0, 0 },   // GOP Index 21 - Ref List 0
    { -1, -3, -11, 1 }    // GOP Index 21 - Ref List 1
    },

    {
        4,                  // GOP Index 22 - Temporal Layer
        22,                 // GOP Index 22 - Decode Order
    { 2,   6, 22, 0 },   // GOP Index 22 - Ref List 0
    { -2, -10,  2, 6 }    // GOP Index 22 - Ref List 1
    },

    {
        5,                  // GOP Index 23 - Temporal Layer
        24,                 // GOP Index 23 - Decode Order
    { 1,  3, 7,  0 },    // GOP Index 23 - Ref List 0
    { -1, -9, 1,  3 }     // GOP Index 23 - Ref List 1
    },

    {
        2,                  // GOP Index 24 - Temporal Layer
        17,                 // GOP Index 24 - Decode Order
    { 8, 24,  0, 0 },    // GOP Index 24 - Ref List 0
    { -8,  8, 24, 0 }     // GOP Index 24 - Ref List 1
    },

    {
        5,                  // GOP Index 25 - Temporal Layer
        27,                 // GOP Index 25 - Decode Order
    { 1,  9,  0, 0 },    // GOP Index 25 - Ref List 0
    { -1, -3, -7, 1 }     // GOP Index 25 - Ref List 1
    },

    {
        4,                  // GOP Index 26 - Temporal Layer
        26,                 // GOP Index 26 - Decode Order
    { 2, 10, 26,  0 },   // GOP Index 26 - Ref List 0
    { -2, -6,  2, 10 }    // GOP Index 26 - Ref List 1
    },

    {
        5,                  // GOP Index 27 - Temporal Layer
        28,                 // GOP Index 27 - Decode Order
    { 1,  3, 11,  0 },   // GOP Index 27 - Ref List 0
    { -1, -5,  1,  3 }    // GOP Index 27 - Ref List 1
    },

    {
        3,                  // GOP Index 28 - Temporal Layer
        25,                 // GOP Index 28 - Decode Order
    { 4, 12, 28,  0 },   // GOP Index 28 - Ref List 0
    { -4,  4, 12, 28 }    // GOP Index 28 - Ref List 1
    },

    {
        5,                  // GOP Index 29 - Temporal Layer
        30,                 // GOP Index 29 - Decode Order
    { 1,  5, 13,  0 },   // GOP Index 29 - Ref List 0
    { -1, -3,  1,  5 }    // GOP Index 29 - Ref List 1
    },

    {
        4,                  // GOP Index 30 - Temporal Layer
        29,                 // GOP Index 30 - Decode Order
    { 2, 6, 14,  0 },    // GOP Index 30 - Ref List 0
    { -2, 2,  6, 14 }     // GOP Index 30 - Ref List 1
    },

    {
        5,                  // GOP Index 31 - Temporal Layer
        31,                 // GOP Index 31 - Decode Order
    { 1, 3, 5, 15 },     // GOP Index 31 - Ref List 0
    { -1, 1, 3,  5 }      // GOP Index 31 - Ref List 1
    }
};
#endif

#ifdef __cplusplus
}
#endif
#endif // EbPredictionStructure_h