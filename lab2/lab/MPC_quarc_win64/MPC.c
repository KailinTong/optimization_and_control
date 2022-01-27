/*
 * MPC.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "MPC".
 *
 * Model version              : 1.175
 * Simulink Coder version : 8.14 (R2018a) 06-Feb-2018
 * C source code generated on : Thu Jan 27 12:26:03 2022
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Generic->Unspecified (assume 32-bit Generic)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "MPC.h"
#include "MPC_private.h"
#include "MPC_dt.h"

/* Block signals (default storage) */
B_MPC_T MPC_B;

/* Block states (default storage) */
DW_MPC_T MPC_DW;

/* Real-time model */
RT_MODEL_MPC_T MPC_M_;
RT_MODEL_MPC_T *const MPC_M = &MPC_M_;

/* Model output function */
void MPC_output(void)
{
  /* local block i/o variables */
  real_T rtb_HILReadAnalog1_o3;
  real_T rtb_HILReadAnalog1_o4;
  real_T rtb_Saturation1;
  real_T rtb_Saturation;
  int32_T iter;
  real_T A_iq_tr[445];
  real_T delta_u_hat_k[4];
  real_T lambda[94];
  real_T J;
  real_T rtb_c_delta[5];
  real_T rtb_b_iq[89];
  real_T rtb_gxn[40];
  real_T tmp[20];
  real_T tmp_0[20];
  real_T tmp_1[80];
  real_T tmp_2[20];
  int32_T i;
  int32_T i_0;
  real_T rtb_Add4_idx_0;
  real_T rtb_Add4_idx_1;

  /* Sum: '<Root>/Add2' incorporates:
   *  Constant: '<Root>/Constant3'
   */
  J = MPC_P.h_02 - MPC_P.xe2;

  /* S-Function (hil_read_analog_block): '<S1>/HIL Read Analog1' */

  /* S-Function Block: MPC/3Tank Modell 1/HIL Read Analog1 (hil_read_analog_block) */
  {
    t_error result = hil_read_analog(MPC_DW.HILInitialize_Card,
      MPC_P.HILReadAnalog1_channels, 4, &MPC_DW.HILReadAnalog1_Buffer[0]);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(MPC_M, _rt_error_message);
    }

    rtb_Saturation1 = MPC_DW.HILReadAnalog1_Buffer[0];
    rtb_Saturation = MPC_DW.HILReadAnalog1_Buffer[1];
    rtb_HILReadAnalog1_o3 = MPC_DW.HILReadAnalog1_Buffer[2];
    rtb_HILReadAnalog1_o4 = MPC_DW.HILReadAnalog1_Buffer[3];
  }

  /* Sum: '<S1>/Sum' incorporates:
   *  Constant: '<S1>/Constant1'
   *  Gain: '<S1>/Gain'
   */
  MPC_B.Sum = MPC_P.gainh1 * rtb_Saturation1 + MPC_P.offseth1;

  /* Sum: '<S1>/Sum1' incorporates:
   *  Constant: '<S1>/Constant2'
   *  Gain: '<S1>/Gain1'
   */
  MPC_B.Sum1 = MPC_P.gainh2 * rtb_Saturation + MPC_P.offseth2;

  /* Sum: '<Root>/Add4' incorporates:
   *  Constant: '<Root>/Constant4'
   *  Constant: '<Root>/Constant7'
   *  Sum: '<Root>/Sum6'
   */
  rtb_Add4_idx_0 = (MPC_P.Constant7_Value[0] + MPC_B.Sum) - MPC_P.xe[0];
  rtb_Add4_idx_1 = (MPC_P.Constant7_Value[1] + MPC_B.Sum1) - MPC_P.xe[1];
  for (iter = 0; iter < 40; iter++) {
    /* Sum: '<Root>/Sum' incorporates:
     *  Delay: '<Root>/Delay'
     *  Gain: '<Root>/Gain'
     *  Gain: '<Root>/Gain3'
     */
    rtb_gxn[iter] = (MPC_P.Fx[iter + 40] * rtb_Add4_idx_1 + MPC_P.Fx[iter] *
                     rtb_Add4_idx_0) + (MPC_P.Gx[iter + 40] *
      MPC_DW.Delay_DSTATE[1] + MPC_P.Gx[iter] * MPC_DW.Delay_DSTATE[0]);
  }

  /* MATLAB Function 'MATLAB Function': '<S2>:1' */
  /* '<S2>:1:2' */
  /* '<S2>:1:3' */
  /* '<S2>:1:4' */
  for (iter = 0; iter < 20; iter++) {
    /* Gain: '<Root>/Gain2' incorporates:
     *  Sum: '<Root>/Sum1'
     */
    tmp[iter] = MPC_P.Fy[iter + 20] * rtb_Add4_idx_1 + MPC_P.Fy[iter] *
      rtb_Add4_idx_0;

    /* Gain: '<Root>/Gain1' incorporates:
     *  Delay: '<Root>/Delay'
     *  Sum: '<Root>/Sum1'
     */
    tmp_0[iter] = MPC_P.Gy[iter + 20] * MPC_DW.Delay_DSTATE[1] + MPC_P.Gy[iter] *
      MPC_DW.Delay_DSTATE[0];
  }

  /* MATLAB Function: '<Root>/MATLAB Function' incorporates:
   *  Constant: '<Root>/Constant6'
   *  Delay: '<Root>/Delay'
   *  Sum: '<Root>/Add2'
   *  Sum: '<Root>/Sum1'
   */
  for (iter = 0; iter < 4; iter++) {
    for (i = 0; i < 20; i++) {
      tmp_1[iter + (i << 2)] = 0.0;
      for (i_0 = 0; i_0 < 20; i_0++) {
        tmp_1[iter + (i << 2)] += MPC_P.Hy[20 * iter + i_0] * 2.0 * MPC_P.Qy[20 *
          i + i_0];
      }
    }
  }

  for (iter = 0; iter < 20; iter++) {
    tmp_2[iter] = (tmp[iter] + tmp_0[iter]) - (MPC_P.Constant6_Value[iter] + J);
  }

  rtb_c_delta[4] = MPC_P.roh;

  /* '<S2>:1:5' */
  for (iter = 0; iter < 4; iter++) {
    delta_u_hat_k[iter] = 0.0;
    for (i = 0; i < 20; i++) {
      delta_u_hat_k[iter] += tmp_1[(i << 2) + iter] * tmp_2[i];
    }

    rtb_c_delta[iter] = delta_u_hat_k[iter];
    J = MPC_P.L[iter + 4] * MPC_DW.Delay_DSTATE[1] + MPC_P.L[iter] *
      MPC_DW.Delay_DSTATE[0];
    rtb_b_iq[iter] = J - MPC_P.u_min[iter];
    rtb_b_iq[iter + 4] = MPC_P.u_max[iter] - J;
  }

  for (iter = 0; iter < 40; iter++) {
    rtb_b_iq[iter + 8] = rtb_gxn[iter] - MPC_P.x_min[iter];
    rtb_b_iq[iter + 48] = MPC_P.x_max[iter] - rtb_gxn[iter];
  }

  rtb_b_iq[88] = MPC_P.roh;

  /* End of MATLAB Function: '<Root>/MATLAB Function' */

  /* MATLAB Function: '<Root>/qpOASES-native' incorporates:
   *  Constant: '<Root>/Constant1'
   *  Constant: '<Root>/Constant2'
   */
  /* MATLAB Function 'qpOASES-native': '<S3>:1' */
  /* '<S3>:1:22' */
  /* '<S3>:1:5' */
  for (iter = 0; iter < 89; iter++) {
    for (i = 0; i < 5; i++) {
      A_iq_tr[i + 5 * iter] = MPC_P.A_iq[89 * i + iter];
    }
  }

  /* '<S3>:1:7' */
  /* '<S3>:1:8' */
  /* '<S3>:1:9' */
  /* '<S3>:1:10' */
  iter = 180;
  if (!MPC_DW.initialized_not_empty) {
    /* '<S3>:1:12' */
    qpOASES_Options_init(&MPC_DW.options, 2);
    QProblem_setup(5, 89, 6);
    MPC_DW.initialized_not_empty = true;
  }

  /* '<S3>:1:20' */
  QProblem_init(MPC_P.W_delta, rtb_c_delta, A_iq_tr, NULL, NULL, NULL, rtb_b_iq,
                &iter, NULL, &MPC_DW.options, delta_u_hat_k, lambda, &J,
                &MPC_B.error);

  /* '<S3>:1:26' */
  MPC_B.iter = iter;

  /* DiscreteTransferFcn: '<Root>/Discrete Transfer Fcn' incorporates:
   *  MATLAB Function: '<Root>/qpOASES-native'
   */
  MPC_DW.DiscreteTransferFcn_tmp[0U] = 0.0;
  MPC_DW.DiscreteTransferFcn_tmp[0] = (delta_u_hat_k[0] -
    MPC_P.DiscreteTransferFcn_DenCoef[1] * MPC_DW.DiscreteTransferFcn_states[0])
    / MPC_P.DiscreteTransferFcn_DenCoef[0];
  MPC_B.u1[0] = MPC_P.DiscreteTransferFcn_NumCoef[0] *
    MPC_DW.DiscreteTransferFcn_tmp[0] + MPC_P.DiscreteTransferFcn_NumCoef[1] *
    MPC_DW.DiscreteTransferFcn_states[0];

  /* Sum: '<Root>/Add' incorporates:
   *  Constant: '<Root>/Constant5'
   */
  MPC_B.Add[0] = MPC_B.u1[0] + MPC_P.ue[0];

  /* DiscreteTransferFcn: '<Root>/Discrete Transfer Fcn' incorporates:
   *  MATLAB Function: '<Root>/qpOASES-native'
   */
  MPC_DW.DiscreteTransferFcn_tmp[1] = (delta_u_hat_k[1] -
    MPC_P.DiscreteTransferFcn_DenCoef[1] * MPC_DW.DiscreteTransferFcn_states[1])
    / MPC_P.DiscreteTransferFcn_DenCoef[0];
  MPC_B.u1[1] = MPC_P.DiscreteTransferFcn_NumCoef[0] *
    MPC_DW.DiscreteTransferFcn_tmp[1] + MPC_P.DiscreteTransferFcn_NumCoef[1] *
    MPC_DW.DiscreteTransferFcn_states[1];

  /* Sum: '<Root>/Add' incorporates:
   *  Constant: '<Root>/Constant5'
   */
  MPC_B.Add[1] = MPC_B.u1[1] + MPC_P.ue[1];

  /* Stop: '<S1>/Stop Simulation' incorporates:
   *  Constant: '<S5>/Constant'
   *  RelationalOperator: '<S5>/Compare'
   */
  if (MPC_B.Sum >= MPC_P.CompareToConstant2_const) {
    rtmSetStopRequested(MPC_M, 1);
  }

  /* End of Stop: '<S1>/Stop Simulation' */

  /* Stop: '<S1>/Stop Simulation1' incorporates:
   *  Constant: '<S4>/Constant'
   *  RelationalOperator: '<S4>/Compare'
   */
  if (MPC_B.Sum1 >= MPC_P.CompareToConstant1_const) {
    rtmSetStopRequested(MPC_M, 1);
  }

  /* End of Stop: '<S1>/Stop Simulation1' */

  /* Saturate: '<S1>/Saturation' */
  if (MPC_B.Add[0] > MPC_P.Saturation_UpperSat) {
    rtb_Saturation = MPC_P.Saturation_UpperSat;
  } else if (MPC_B.Add[0] < MPC_P.Saturation_LowerSat) {
    rtb_Saturation = MPC_P.Saturation_LowerSat;
  } else {
    rtb_Saturation = MPC_B.Add[0];
  }

  /* End of Saturate: '<S1>/Saturation' */

  /* Saturate: '<S1>/Saturation1' incorporates:
   *  Constant: '<Root>/Constant'
   */
  if (MPC_P.Constant_Value > MPC_P.Saturation1_UpperSat) {
    rtb_Saturation1 = MPC_P.Saturation1_UpperSat;
  } else if (MPC_P.Constant_Value < MPC_P.Saturation1_LowerSat) {
    rtb_Saturation1 = MPC_P.Saturation1_LowerSat;
  } else {
    rtb_Saturation1 = MPC_P.Constant_Value;
  }

  /* End of Saturate: '<S1>/Saturation1' */

  /* S-Function (hil_write_analog_block): '<S1>/HIL Write Analog2' */

  /* S-Function Block: MPC/3Tank Modell 1/HIL Write Analog2 (hil_write_analog_block) */
  {
    t_error result;
    MPC_DW.HILWriteAnalog2_Buffer[0] = rtb_Saturation;
    MPC_DW.HILWriteAnalog2_Buffer[1] = rtb_Saturation1;
    result = hil_write_analog(MPC_DW.HILInitialize_Card,
      MPC_P.HILWriteAnalog2_channels, 2, &MPC_DW.HILWriteAnalog2_Buffer[0]);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(MPC_M, _rt_error_message);
    }
  }
}

/* Model update function */
void MPC_update(void)
{
  /* Update for Delay: '<Root>/Delay' */
  MPC_DW.Delay_DSTATE[0] = MPC_B.u1[0];

  /* Update for DiscreteTransferFcn: '<Root>/Discrete Transfer Fcn' */
  MPC_DW.DiscreteTransferFcn_states[0] = MPC_DW.DiscreteTransferFcn_tmp[0];

  /* Update for Delay: '<Root>/Delay' */
  MPC_DW.Delay_DSTATE[1] = MPC_B.u1[1];

  /* Update for DiscreteTransferFcn: '<Root>/Discrete Transfer Fcn' */
  MPC_DW.DiscreteTransferFcn_states[1] = MPC_DW.DiscreteTransferFcn_tmp[1];

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++MPC_M->Timing.clockTick0)) {
    ++MPC_M->Timing.clockTickH0;
  }

  MPC_M->Timing.t[0] = MPC_M->Timing.clockTick0 * MPC_M->Timing.stepSize0 +
    MPC_M->Timing.clockTickH0 * MPC_M->Timing.stepSize0 * 4294967296.0;
}

/* Model initialize function */
void MPC_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: MPC/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &MPC_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(MPC_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(MPC_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(MPC_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(MPC_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(MPC_M, _rt_error_message);
      return;
    }

    if ((MPC_P.HILInitialize_AIPStart && !is_switching) ||
        (MPC_P.HILInitialize_AIPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &MPC_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = (MPC_P.HILInitialize_AILow);
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &MPC_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = MPC_P.HILInitialize_AIHigh;
        }
      }

      result = hil_set_analog_input_ranges(MPC_DW.HILInitialize_Card,
        MPC_P.HILInitialize_AIChannels, 8U,
        &MPC_DW.HILInitialize_AIMinimums[0], &MPC_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(MPC_M, _rt_error_message);
        return;
      }
    }

    if ((MPC_P.HILInitialize_AOPStart && !is_switching) ||
        (MPC_P.HILInitialize_AOPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &MPC_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = (MPC_P.HILInitialize_AOLow);
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &MPC_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = MPC_P.HILInitialize_AOHigh;
        }
      }

      result = hil_set_analog_output_ranges(MPC_DW.HILInitialize_Card,
        MPC_P.HILInitialize_AOChannels, 8U,
        &MPC_DW.HILInitialize_AOMinimums[0], &MPC_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(MPC_M, _rt_error_message);
        return;
      }
    }

    if ((MPC_P.HILInitialize_AOStart && !is_switching) ||
        (MPC_P.HILInitialize_AOEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &MPC_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = MPC_P.HILInitialize_AOInitial;
        }
      }

      result = hil_write_analog(MPC_DW.HILInitialize_Card,
        MPC_P.HILInitialize_AOChannels, 8U, &MPC_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(MPC_M, _rt_error_message);
        return;
      }
    }

    if (MPC_P.HILInitialize_AOReset) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &MPC_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = MPC_P.HILInitialize_AOWatchdog;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (MPC_DW.HILInitialize_Card, MPC_P.HILInitialize_AOChannels, 8U,
         &MPC_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(MPC_M, _rt_error_message);
        return;
      }
    }

    if ((MPC_P.HILInitialize_EIPStart && !is_switching) ||
        (MPC_P.HILInitialize_EIPEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes = &MPC_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = MPC_P.HILInitialize_EIQuadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode(MPC_DW.HILInitialize_Card,
        MPC_P.HILInitialize_EIChannels, 8U, (t_encoder_quadrature_mode *)
        &MPC_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(MPC_M, _rt_error_message);
        return;
      }
    }

    if ((MPC_P.HILInitialize_EIStart && !is_switching) ||
        (MPC_P.HILInitialize_EIEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts = &MPC_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] = MPC_P.HILInitialize_EIInitial;
        }
      }

      result = hil_set_encoder_counts(MPC_DW.HILInitialize_Card,
        MPC_P.HILInitialize_EIChannels, 8U,
        &MPC_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(MPC_M, _rt_error_message);
        return;
      }
    }

    if ((MPC_P.HILInitialize_POPStart && !is_switching) ||
        (MPC_P.HILInitialize_POPEnter && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues = &MPC_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = MPC_P.HILInitialize_POModes;
        }
      }

      result = hil_set_pwm_mode(MPC_DW.HILInitialize_Card,
        MPC_P.HILInitialize_POChannels, 8U, (t_pwm_mode *)
        &MPC_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(MPC_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_POChannels =
          MPC_P.HILInitialize_POChannels;
        int32_T *dw_POModeValues = &MPC_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE ||
              dw_POModeValues[i1] == PWM_RAW_MODE) {
            MPC_DW.HILInitialize_POSortedChans[num_duty_cycle_modes] =
              (p_HILInitialize_POChannels[i1]);
            MPC_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes] =
              MPC_P.HILInitialize_POFrequency;
            num_duty_cycle_modes++;
          } else {
            MPC_DW.HILInitialize_POSortedChans[7U - num_frequency_modes] =
              (p_HILInitialize_POChannels[i1]);
            MPC_DW.HILInitialize_POSortedFreqs[7U - num_frequency_modes] =
              MPC_P.HILInitialize_POFrequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(MPC_DW.HILInitialize_Card,
          &MPC_DW.HILInitialize_POSortedChans[0], num_duty_cycle_modes,
          &MPC_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(MPC_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(MPC_DW.HILInitialize_Card,
          &MPC_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &MPC_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(MPC_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues = &MPC_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = MPC_P.HILInitialize_POConfiguration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues = &MPC_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = MPC_P.HILInitialize_POAlignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals = &MPC_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = MPC_P.HILInitialize_POPolarity;
        }
      }

      result = hil_set_pwm_configuration(MPC_DW.HILInitialize_Card,
        MPC_P.HILInitialize_POChannels, 8U,
        (t_pwm_configuration *) &MPC_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &MPC_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &MPC_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(MPC_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs = &MPC_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = MPC_P.HILInitialize_POLeading;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &MPC_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = MPC_P.HILInitialize_POTrailing;
        }
      }

      result = hil_set_pwm_deadband(MPC_DW.HILInitialize_Card,
        MPC_P.HILInitialize_POChannels, 8U,
        &MPC_DW.HILInitialize_POSortedFreqs[0], &MPC_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(MPC_M, _rt_error_message);
        return;
      }
    }

    if ((MPC_P.HILInitialize_POStart && !is_switching) ||
        (MPC_P.HILInitialize_POEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &MPC_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = MPC_P.HILInitialize_POInitial;
        }
      }

      result = hil_write_pwm(MPC_DW.HILInitialize_Card,
        MPC_P.HILInitialize_POChannels, 8U, &MPC_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(MPC_M, _rt_error_message);
        return;
      }
    }

    if (MPC_P.HILInitialize_POReset) {
      {
        int_T i1;
        real_T *dw_POValues = &MPC_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = MPC_P.HILInitialize_POWatchdog;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state(MPC_DW.HILInitialize_Card,
        MPC_P.HILInitialize_POChannels, 8U, &MPC_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(MPC_M, _rt_error_message);
        return;
      }
    }
  }

  /* InitializeConditions for Delay: '<Root>/Delay' */
  MPC_DW.Delay_DSTATE[0] = MPC_P.Delay_InitialCondition[0];

  /* InitializeConditions for DiscreteTransferFcn: '<Root>/Discrete Transfer Fcn' */
  MPC_DW.DiscreteTransferFcn_states[0] = MPC_P.DiscreteTransferFcn_InitialStat[0];

  /* InitializeConditions for Delay: '<Root>/Delay' */
  MPC_DW.Delay_DSTATE[1] = MPC_P.Delay_InitialCondition[1];

  /* InitializeConditions for DiscreteTransferFcn: '<Root>/Discrete Transfer Fcn' */
  MPC_DW.DiscreteTransferFcn_states[1] = MPC_P.DiscreteTransferFcn_InitialStat[1];

  /* SystemInitialize for MATLAB Function: '<Root>/qpOASES-native' */
  MPC_DW.initialized_not_empty = false;
}

/* Model terminate function */
void MPC_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: MPC/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(MPC_DW.HILInitialize_Card);
    hil_monitor_stop_all(MPC_DW.HILInitialize_Card);
    is_switching = false;
    if ((MPC_P.HILInitialize_AOTerminate && !is_switching) ||
        (MPC_P.HILInitialize_AOExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &MPC_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = MPC_P.HILInitialize_AOFinal;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((MPC_P.HILInitialize_POTerminate && !is_switching) ||
        (MPC_P.HILInitialize_POExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &MPC_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = MPC_P.HILInitialize_POFinal;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(MPC_DW.HILInitialize_Card
                         , MPC_P.HILInitialize_AOChannels,
                         num_final_analog_outputs
                         , MPC_P.HILInitialize_POChannels, num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &MPC_DW.HILInitialize_AOVoltages[0]
                         , &MPC_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(MPC_DW.HILInitialize_Card,
            MPC_P.HILInitialize_AOChannels, num_final_analog_outputs,
            &MPC_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(MPC_DW.HILInitialize_Card,
            MPC_P.HILInitialize_POChannels, num_final_pwm_outputs,
            &MPC_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(MPC_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(MPC_DW.HILInitialize_Card);
    hil_monitor_delete_all(MPC_DW.HILInitialize_Card);
    hil_close(MPC_DW.HILInitialize_Card);
    MPC_DW.HILInitialize_Card = NULL;
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/
void MdlOutputs(int_T tid)
{
  MPC_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  MPC_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  MPC_initialize();
}

void MdlTerminate(void)
{
  MPC_terminate();
}

/* Registration function */
RT_MODEL_MPC_T *MPC(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)MPC_M, 0,
                sizeof(RT_MODEL_MPC_T));

  /* Initialize timing info */
  {
    int_T *mdlTsMap = MPC_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    MPC_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    MPC_M->Timing.sampleTimes = (&MPC_M->Timing.sampleTimesArray[0]);
    MPC_M->Timing.offsetTimes = (&MPC_M->Timing.offsetTimesArray[0]);

    /* task periods */
    MPC_M->Timing.sampleTimes[0] = (0.5);

    /* task offsets */
    MPC_M->Timing.offsetTimes[0] = (0.0);
  }

  rtmSetTPtr(MPC_M, &MPC_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = MPC_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    MPC_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(MPC_M, -1);
  MPC_M->Timing.stepSize0 = 0.5;

  /* External mode info */
  MPC_M->Sizes.checksums[0] = (3842294566U);
  MPC_M->Sizes.checksums[1] = (3972289132U);
  MPC_M->Sizes.checksums[2] = (2765226405U);
  MPC_M->Sizes.checksums[3] = (1246815104U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[3];
    MPC_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = &rtAlwaysEnabled;
    systemRan[2] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(MPC_M->extModeInfo,
      &MPC_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(MPC_M->extModeInfo, MPC_M->Sizes.checksums);
    rteiSetTPtr(MPC_M->extModeInfo, rtmGetTPtr(MPC_M));
  }

  MPC_M->solverInfoPtr = (&MPC_M->solverInfo);
  MPC_M->Timing.stepSize = (0.5);
  rtsiSetFixedStepSize(&MPC_M->solverInfo, 0.5);
  rtsiSetSolverMode(&MPC_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  MPC_M->blockIO = ((void *) &MPC_B);
  (void) memset(((void *) &MPC_B), 0,
                sizeof(B_MPC_T));

  /* parameters */
  MPC_M->defaultParam = ((real_T *)&MPC_P);

  /* states (dwork) */
  MPC_M->dwork = ((void *) &MPC_DW);
  (void) memset((void *)&MPC_DW, 0,
                sizeof(DW_MPC_T));

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    MPC_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  MPC_M->Sizes.numContStates = (0);    /* Number of continuous states */
  MPC_M->Sizes.numY = (0);             /* Number of model outputs */
  MPC_M->Sizes.numU = (0);             /* Number of model inputs */
  MPC_M->Sizes.sysDirFeedThru = (0);   /* The model is not direct feedthrough */
  MPC_M->Sizes.numSampTimes = (1);     /* Number of sample times */
  MPC_M->Sizes.numBlocks = (45);       /* Number of blocks */
  MPC_M->Sizes.numBlockIO = (7);       /* Number of block outputs */
  MPC_M->Sizes.numBlockPrms = (1437);  /* Sum of parameter "widths" */
  return MPC_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
