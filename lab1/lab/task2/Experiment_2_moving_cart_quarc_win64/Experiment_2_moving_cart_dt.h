/*
 * Experiment_2_moving_cart_dt.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "Experiment_2_moving_cart".
 *
 * Model version              : 1.41
 * Simulink Coder version : 8.14 (R2018a) 06-Feb-2018
 * C source code generated on : Wed Dec 15 14:52:48 2021
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "ext_types.h"

/* data type size table */
static uint_T rtDataTypeSizes[] = {
  sizeof(real_T),
  sizeof(real32_T),
  sizeof(int8_T),
  sizeof(uint8_T),
  sizeof(int16_T),
  sizeof(uint16_T),
  sizeof(int32_T),
  sizeof(uint32_T),
  sizeof(boolean_T),
  sizeof(fcn_call_T),
  sizeof(int_T),
  sizeof(pointer_T),
  sizeof(action_T),
  2*sizeof(uint32_T),
  sizeof(t_card),
  sizeof(emxArray_real_T_3x3_Experimen_T),
  sizeof(emxArray_real_T_3_Experiment__T),
  sizeof(emxArray_real_T_2x2_Experimen_T),
  sizeof(emxArray_real_T_2_Experiment__T)
};

/* data type name table */
static const char_T * rtDataTypeNames[] = {
  "real_T",
  "real32_T",
  "int8_T",
  "uint8_T",
  "int16_T",
  "uint16_T",
  "int32_T",
  "uint32_T",
  "boolean_T",
  "fcn_call_T",
  "int_T",
  "pointer_T",
  "action_T",
  "timer_uint32_pair_T",
  "t_card",
  "emxArray_real_T_3x3_Experimen_T",
  "emxArray_real_T_3_Experiment__T",
  "emxArray_real_T_2x2_Experimen_T",
  "emxArray_real_T_2_Experiment__T"
};

/* data type transitions for block I/O structure */
static DataTypeTransition rtBTransitions[] = {
  { (char_T *)(&Experiment_2_moving_cart_B.Gain1), 0, 0, 16 }
  ,

  { (char_T *)(&Experiment_2_moving_cart_DW.Phi), 15, 0, 1 },

  { (char_T *)(&Experiment_2_moving_cart_DW.Phi_b), 17, 0, 1 },

  { (char_T *)(&Experiment_2_moving_cart_DW.bD), 16, 0, 3 },

  { (char_T *)(&Experiment_2_moving_cart_DW.bD_e), 18, 0, 3 },

  { (char_T *)(&Experiment_2_moving_cart_DW.HILInitialize_AIMinimums[0]), 0, 0,
    64 },

  { (char_T *)(&Experiment_2_moving_cart_DW.HILInitialize_Card), 14, 0, 1 },

  { (char_T *)(&Experiment_2_moving_cart_DW.TransportDelay_RWORK.modelTStart), 0,
    0, 1 },

  { (char_T *)(&Experiment_2_moving_cart_DW.HILReadEncoder_PWORK), 11, 0, 21 },

  { (char_T *)(&Experiment_2_moving_cart_DW.HILInitialize_ClockModes[0]), 6, 0,
    45 },

  { (char_T *)(&Experiment_2_moving_cart_DW.HILInitialize_POSortedChans[0]), 7,
    0, 8 },

  { (char_T *)(&Experiment_2_moving_cart_DW.TransportDelay_IWORK.Tail), 10, 0, 1
  },

  { (char_T *)(&Experiment_2_moving_cart_DW.zp_not_empty), 8, 0, 2 }
};

/* data type transition table for block I/O structure */
static DataTypeTransitionTable rtBTransTable = {
  13U,
  rtBTransitions
};

/* data type transitions for Parameters structure */
static DataTypeTransition rtPTransitions[] = {
  { (char_T *)(&Experiment_2_moving_cart_P.gain_cart), 0, 0, 8 },

  { (char_T *)(&Experiment_2_moving_cart_P.HILReadEncoder_channels[0]), 7, 0, 3
  },

  { (char_T *)(&Experiment_2_moving_cart_P.Gain2_Gain), 0, 0, 33 },

  { (char_T *)(&Experiment_2_moving_cart_P.HILInitialize_CKChannels[0]), 6, 0, 9
  },

  { (char_T *)(&Experiment_2_moving_cart_P.HILInitialize_AIChannels[0]), 7, 0,
    33 },

  { (char_T *)(&Experiment_2_moving_cart_P.HILInitialize_Active), 8, 0, 37 },

  { (char_T *)(&Experiment_2_moving_cart_P.ManualSwitch_CurrentSetting), 3, 0, 1
  }
};

/* data type transition table for Parameters structure */
static DataTypeTransitionTable rtPTransTable = {
  7U,
  rtPTransitions
};

/* [EOF] Experiment_2_moving_cart_dt.h */
