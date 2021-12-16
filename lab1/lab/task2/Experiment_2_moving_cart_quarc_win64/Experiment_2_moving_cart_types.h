/*
 * Experiment_2_moving_cart_types.h
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

#ifndef RTW_HEADER_Experiment_2_moving_cart_types_h_
#define RTW_HEADER_Experiment_2_moving_cart_types_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "zero_crossing_types.h"

/* Custom Type definition for MATLAB Function: '<S5>/Differentiator' */
#ifndef struct_emxArray_real_T_3x3
#define struct_emxArray_real_T_3x3

struct emxArray_real_T_3x3
{
  real_T data[9];
  int32_T size[2];
};

#endif                                 /*struct_emxArray_real_T_3x3*/

#ifndef typedef_emxArray_real_T_3x3_Experimen_T
#define typedef_emxArray_real_T_3x3_Experimen_T

typedef struct emxArray_real_T_3x3 emxArray_real_T_3x3_Experimen_T;

#endif                                 /*typedef_emxArray_real_T_3x3_Experimen_T*/

/* Custom Type definition for MATLAB Function: '<S6>/Differentiator' */
#ifndef struct_emxArray_real_T_2x2
#define struct_emxArray_real_T_2x2

struct emxArray_real_T_2x2
{
  real_T data[4];
  int32_T size[2];
};

#endif                                 /*struct_emxArray_real_T_2x2*/

#ifndef typedef_emxArray_real_T_2x2_Experimen_T
#define typedef_emxArray_real_T_2x2_Experimen_T

typedef struct emxArray_real_T_2x2 emxArray_real_T_2x2_Experimen_T;

#endif                                 /*typedef_emxArray_real_T_2x2_Experimen_T*/

#ifndef struct_emxArray_real_T_2x3
#define struct_emxArray_real_T_2x3

struct emxArray_real_T_2x3
{
  real_T data[6];
  int32_T size[2];
};

#endif                                 /*struct_emxArray_real_T_2x3*/

#ifndef typedef_emxArray_real_T_2x3_Experimen_T
#define typedef_emxArray_real_T_2x3_Experimen_T

typedef struct emxArray_real_T_2x3 emxArray_real_T_2x3_Experimen_T;

#endif                                 /*typedef_emxArray_real_T_2x3_Experimen_T*/

#ifndef struct_sdJmLrH9dQMRvcKhJSSLdyB_tag
#define struct_sdJmLrH9dQMRvcKhJSSLdyB_tag

struct sdJmLrH9dQMRvcKhJSSLdyB_tag
{
  emxArray_real_T_2x3_Experimen_T f1;
};

#endif                                 /*struct_sdJmLrH9dQMRvcKhJSSLdyB_tag*/

#ifndef typedef_b_cell_wrap_1_Experiment_2_mo_T
#define typedef_b_cell_wrap_1_Experiment_2_mo_T

typedef struct sdJmLrH9dQMRvcKhJSSLdyB_tag b_cell_wrap_1_Experiment_2_mo_T;

#endif                                 /*typedef_b_cell_wrap_1_Experiment_2_mo_T*/

#ifndef struct_emxArray_real_T_3x4
#define struct_emxArray_real_T_3x4

struct emxArray_real_T_3x4
{
  real_T data[12];
  int32_T size[2];
};

#endif                                 /*struct_emxArray_real_T_3x4*/

#ifndef typedef_emxArray_real_T_3x4_Experimen_T
#define typedef_emxArray_real_T_3x4_Experimen_T

typedef struct emxArray_real_T_3x4 emxArray_real_T_3x4_Experimen_T;

#endif                                 /*typedef_emxArray_real_T_3x4_Experimen_T*/

#ifndef struct_sgOmh9KUKPVuCREeTuTuOQD_tag
#define struct_sgOmh9KUKPVuCREeTuTuOQD_tag

struct sgOmh9KUKPVuCREeTuTuOQD_tag
{
  emxArray_real_T_3x4_Experimen_T f1;
};

#endif                                 /*struct_sgOmh9KUKPVuCREeTuTuOQD_tag*/

#ifndef typedef_b_cell_wrap_1_Experiment_2__n_T
#define typedef_b_cell_wrap_1_Experiment_2__n_T

typedef struct sgOmh9KUKPVuCREeTuTuOQD_tag b_cell_wrap_1_Experiment_2__n_T;

#endif                                 /*typedef_b_cell_wrap_1_Experiment_2__n_T*/

/* Custom Type definition for MATLAB Function: '<S5>/Differentiator' */
#ifndef struct_emxArray_real_T_3
#define struct_emxArray_real_T_3

struct emxArray_real_T_3
{
  real_T data[3];
  int32_T size;
};

#endif                                 /*struct_emxArray_real_T_3*/

#ifndef typedef_emxArray_real_T_3_Experiment__T
#define typedef_emxArray_real_T_3_Experiment__T

typedef struct emxArray_real_T_3 emxArray_real_T_3_Experiment__T;

#endif                                 /*typedef_emxArray_real_T_3_Experiment__T*/

/* Custom Type definition for MATLAB Function: '<S6>/Differentiator' */
#ifndef struct_emxArray_real_T_2
#define struct_emxArray_real_T_2

struct emxArray_real_T_2
{
  real_T data[2];
  int32_T size;
};

#endif                                 /*struct_emxArray_real_T_2*/

#ifndef typedef_emxArray_real_T_2_Experiment__T
#define typedef_emxArray_real_T_2_Experiment__T

typedef struct emxArray_real_T_2 emxArray_real_T_2_Experiment__T;

#endif                                 /*typedef_emxArray_real_T_2_Experiment__T*/

/* Parameters (default storage) */
typedef struct P_Experiment_2_moving_cart_T_ P_Experiment_2_moving_cart_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_Experiment_2_moving_c_T RT_MODEL_Experiment_2_moving__T;

#endif                                 /* RTW_HEADER_Experiment_2_moving_cart_types_h_ */
