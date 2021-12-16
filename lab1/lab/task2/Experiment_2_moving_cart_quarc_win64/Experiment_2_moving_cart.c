/*
 * Experiment_2_moving_cart.c
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

#include "Experiment_2_moving_cart.h"
#include "Experiment_2_moving_cart_private.h"
#include "Experiment_2_moving_cart_dt.h"

/* Named constants for MATLAB Function: '<S5>/Differentiator' */
#define Experiment_2_moving_cart_r     (15.0)

/* Named constants for MATLAB Function: '<S6>/Differentiator' */
#define Experiment_2_moving_cart_r_d   (75.0)

/* Block signals (default storage) */
B_Experiment_2_moving_cart_T Experiment_2_moving_cart_B;

/* Block states (default storage) */
DW_Experiment_2_moving_cart_T Experiment_2_moving_cart_DW;

/* Real-time model */
RT_MODEL_Experiment_2_moving__T Experiment_2_moving_cart_M_;
RT_MODEL_Experiment_2_moving__T *const Experiment_2_moving_cart_M =
  &Experiment_2_moving_cart_M_;

/* Forward declaration for local functions */
static real_T Experiment_2_moving_cart_norm_j(const real_T x_data[], const
  int32_T x_size[2]);
static void Expe_eml_signed_integer_colon_o(int32_T b, int32_T y_data[], int32_T
  y_size[2]);
static real_T Experiment_2_moving_cart_xnrm2(int32_T b_n, const real_T x_data[],
  int32_T ix0);
static void Experiment_2_moving_ca_xgeqp3_l(real_T A_data[], int32_T A_size[2],
  real_T tau_data[], int32_T *tau_size, int32_T jpvt_data[], int32_T jpvt_size[2]);
static void Experiment_2_moving_c_qrsolve_g(const real_T A_data[], const int32_T
  A_size[2], const real_T B_data[], const int32_T B_size[2], real_T Y_data[],
  int32_T Y_size[2]);
static void Experiment_2_moving_c_xzgetrf_b(int32_T b_m, int32_T b_n, real_T
  A_data[], int32_T A_size[2], int32_T lda, int32_T ipiv_data[], int32_T
  ipiv_size[2], int32_T *info);
static void Experiment_2_moving_car_xtrsm_b(const real_T A_data[], real_T
  B_data[], int32_T B_size[2]);
static void Experiment_2_moving__mldivide_g(const real_T A_data[], const int32_T
  A_size[2], const real_T B_data[], const int32_T B_size[2], real_T Y_data[],
  int32_T Y_size[2]);
static void Expe_PadeApproximantOfDegree_kh(const real_T A_data[], const int32_T
  A_size[2], real_T F_data[], int32_T F_size[2]);
static void Exper_PadeApproximantOfDegree_k(const real_T A_data[], const int32_T
  A_size[2], uint8_T b_m, real_T F_data[], int32_T F_size[2]);
static void Experiment_2_moving__mrdivide_f(const real_T A_data[], const int32_T
  A_size[2], const real_T B_data[], const int32_T B_size[2], real_T y_data[],
  int32_T y_size[2]);
static real_T Experim_disc_eigenvalues_URED_d(real_T x0, real_T b_mu, real_T tau);
static void Experime_ackerman_precomputed_g(real_T b_Ts, real_T z, real_T
  lambda_data[], int32_T *lambda_size);
static void Experiment_2_moving_cart_step_n(real_T x0, const real_T z_data[],
  real_T u, const real_T Phi_data[], const real_T bD_data[], const real_T
  lambda_data[], real_T z_new_data[], int32_T *z_new_size);
static real_T Experiment_2_moving_cart_norm(const real_T x_data[], const int32_T
  x_size[2]);
static void Experi_eml_signed_integer_colon(int32_T b, int32_T y_data[], int32_T
  y_size[2]);
static void Experiment_2_moving_cart_xgeqp3(real_T A_data[], int32_T A_size[2],
  real_T tau_data[], int32_T *tau_size, int32_T jpvt_data[], int32_T jpvt_size[2]);
static void Experiment_2_moving_car_qrsolve(const real_T A_data[], const int32_T
  A_size[2], const real_T B_data[], const int32_T B_size[2], real_T Y_data[],
  int32_T Y_size[2]);
static void Experiment_2_moving_car_xzgetrf(int32_T b_m, int32_T b_n, real_T
  A_data[], int32_T A_size[2], int32_T lda, int32_T ipiv_data[], int32_T
  ipiv_size[2], int32_T *info);
static void Experiment_2_moving_cart_xtrsm(const real_T A_data[], real_T B_data[],
  int32_T B_size[2]);
static void Experiment_2_moving_ca_mldivide(const real_T A_data[], const int32_T
  A_size[2], const real_T B_data[], const int32_T B_size[2], real_T Y_data[],
  int32_T Y_size[2]);
static void Exper_PadeApproximantOfDegree_e(const real_T A_data[], const int32_T
  A_size[2], real_T F_data[], int32_T F_size[2]);
static void Experim_PadeApproximantOfDegree(const real_T A_data[], const int32_T
  A_size[2], uint8_T b_m, real_T F_data[], int32_T F_size[2]);
static void Experiment_2_moving_ca_mrdivide(const real_T A_data[], const int32_T
  A_size[2], const real_T B_data[], const int32_T B_size[2], real_T y_data[],
  int32_T y_size[2]);
static real_T Experimen_disc_eigenvalues_URED(real_T x0, real_T b_mu, real_T tau);
static void Experiment_ackerman_precomputed(real_T b_Ts, real_T z, real_T
  lambda_data[], int32_T *lambda_size);
static void Experiment_2_moving_cart_step_p(real_T x0, const real_T z_data[],
  real_T u, const real_T Phi_data[], const real_T bD_data[], const real_T
  lambda_data[], real_T z_new_data[], int32_T *z_new_size);

/*
 * Time delay interpolation routine
 *
 * The linear interpolation is performed using the formula:
 *
 *          (t2 - tMinusDelay)         (tMinusDelay - t1)
 * u(t)  =  ----------------- * u1  +  ------------------- * u2
 *              (t2 - t1)                  (t2 - t1)
 */
real_T rt_TDelayInterpolate(
  real_T tMinusDelay,                  /* tMinusDelay = currentSimTime - delay */
  real_T tStart,
  real_T *tBuf,
  real_T *uBuf,
  int_T bufSz,
  int_T *lastIdx,
  int_T oldestIdx,
  int_T newIdx,
  real_T initOutput,
  boolean_T discrete,
  boolean_T minorStepAndTAtLastMajorOutput)
{
  int_T i;
  real_T yout, t1, t2, u1, u2;

  /*
   * If there is only one data point in the buffer, this data point must be
   * the t= 0 and tMinusDelay > t0, it ask for something unknown. The best
   * guess if initial output as well
   */
  if ((newIdx == 0) && (oldestIdx ==0 ) && (tMinusDelay > tStart))
    return initOutput;

  /*
   * If tMinusDelay is less than zero, should output initial value
   */
  if (tMinusDelay <= tStart)
    return initOutput;

  /* For fixed buffer extrapolation:
   * if tMinusDelay is small than the time at oldestIdx, if discrete, output
   * tailptr value,  else use tailptr and tailptr+1 value to extrapolate
   * It is also for fixed buffer. Note: The same condition can happen for transport delay block where
   * use tStart and and t[tail] other than using t[tail] and t[tail+1].
   * See below
   */
  if ((tMinusDelay <= tBuf[oldestIdx] ) ) {
    if (discrete) {
      return(uBuf[oldestIdx]);
    } else {
      int_T tempIdx= oldestIdx + 1;
      if (oldestIdx == bufSz-1)
        tempIdx = 0;
      t1= tBuf[oldestIdx];
      t2= tBuf[tempIdx];
      u1= uBuf[oldestIdx];
      u2= uBuf[tempIdx];
      if (t2 == t1) {
        if (tMinusDelay >= t2) {
          yout = u2;
        } else {
          yout = u1;
        }
      } else {
        real_T f1 = (t2-tMinusDelay) / (t2-t1);
        real_T f2 = 1.0 - f1;

        /*
         * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
         */
        yout = f1*u1 + f2*u2;
      }

      return yout;
    }
  }

  /*
   * When block does not have direct feedthrough, we use the table of
   * values to extrapolate off the end of the table for delays that are less
   * than 0 (less then step size).  This is not completely accurate.  The
   * chain of events is as follows for a given time t.  Major output - look
   * in table.  Update - add entry to table.  Now, if we call the output at
   * time t again, there is a new entry in the table. For very small delays,
   * this means that we will have a different answer from the previous call
   * to the output fcn at the same time t.  The following code prevents this
   * from happening.
   */
  if (minorStepAndTAtLastMajorOutput) {
    /* pretend that the new entry has not been added to table */
    if (newIdx != 0) {
      if (*lastIdx == newIdx) {
        (*lastIdx)--;
      }

      newIdx--;
    } else {
      if (*lastIdx == newIdx) {
        *lastIdx = bufSz-1;
      }

      newIdx = bufSz - 1;
    }
  }

  i = *lastIdx;
  if (tBuf[i] < tMinusDelay) {
    /* Look forward starting at last index */
    while (tBuf[i] < tMinusDelay) {
      /* May occur if the delay is less than step-size - extrapolate */
      if (i == newIdx)
        break;
      i = ( i < (bufSz-1) ) ? (i+1) : 0;/* move through buffer */
    }
  } else {
    /*
     * Look backwards starting at last index which can happen when the
     * delay time increases.
     */
    while (tBuf[i] >= tMinusDelay) {
      /*
       * Due to the entry condition at top of function, we
       * should never hit the end.
       */
      i = (i > 0) ? i-1 : (bufSz-1);   /* move through buffer */
    }

    i = ( i < (bufSz-1) ) ? (i+1) : 0;
  }

  *lastIdx = i;
  if (discrete) {
    /*
     * tempEps = 128 * eps;
     * localEps = max(tempEps, tempEps*fabs(tBuf[i]))/2;
     */
    double tempEps = (DBL_EPSILON) * 128.0;
    double localEps = tempEps * fabs(tBuf[i]);
    if (tempEps > localEps) {
      localEps = tempEps;
    }

    localEps = localEps / 2.0;
    if (tMinusDelay >= (tBuf[i] - localEps)) {
      yout = uBuf[i];
    } else {
      if (i == 0) {
        yout = uBuf[bufSz-1];
      } else {
        yout = uBuf[i-1];
      }
    }
  } else {
    if (i == 0) {
      t1 = tBuf[bufSz-1];
      u1 = uBuf[bufSz-1];
    } else {
      t1 = tBuf[i-1];
      u1 = uBuf[i-1];
    }

    t2 = tBuf[i];
    u2 = uBuf[i];
    if (t2 == t1) {
      if (tMinusDelay >= t2) {
        yout = u2;
      } else {
        yout = u1;
      }
    } else {
      real_T f1 = (t2-tMinusDelay) / (t2-t1);
      real_T f2 = 1.0 - f1;

      /*
       * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
       */
      yout = f1*u1 + f2*u2;
    }
  }

  return(yout);
}

/* Function for MATLAB Function: '<S6>/Differentiator' */
static real_T Experiment_2_moving_cart_norm_j(const real_T x_data[], const
  int32_T x_size[2])
{
  real_T y;
  real_T s;
  int32_T j;
  int32_T i;
  boolean_T exitg1;
  if (x_size[0] == 0) {
    y = 0.0;
  } else if (x_size[0] == 1) {
    y = (fabs(x_data[0]) + fabs(x_data[1])) + fabs(x_data[2]);
  } else {
    y = 0.0;
    j = 0;
    exitg1 = false;
    while ((!exitg1) && (j <= 2)) {
      s = 0.0;
      for (i = 0; i < x_size[0]; i++) {
        s += fabs(x_data[x_size[0] * j + i]);
      }

      if (rtIsNaN(s)) {
        y = (rtNaN);
        exitg1 = true;
      } else {
        if (s > y) {
          y = s;
        }

        j++;
      }
    }
  }

  return y;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  real_T tmp_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/* Function for MATLAB Function: '<S6>/Differentiator' */
static void Expe_eml_signed_integer_colon_o(int32_T b, int32_T y_data[], int32_T
  y_size[2])
{
  int32_T yk;
  int32_T k;
  y_size[0] = 1;
  y_size[1] = b;
  y_data[0] = 1;
  yk = 1;
  for (k = 2; k <= b; k++) {
    yk++;
    y_data[k - 1] = yk;
  }
}

/* Function for MATLAB Function: '<S5>/Differentiator' */
static real_T Experiment_2_moving_cart_xnrm2(int32_T b_n, const real_T x_data[],
  int32_T ix0)
{
  real_T y;
  real_T scale;
  int32_T kend;
  real_T absxk;
  real_T t;
  int32_T k;
  y = 0.0;
  if (!(b_n < 1)) {
    if (b_n == 1) {
      y = fabs(x_data[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = (ix0 + b_n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = fabs(x_data[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  return y;
}

real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T a;
  a = fabs(u0);
  y = fabs(u1);
  if (a < y) {
    a /= y;
    y *= sqrt(a * a + 1.0);
  } else if (a > y) {
    y /= a;
    y = sqrt(y * y + 1.0) * a;
  } else {
    if (!rtIsNaN(y)) {
      y = a * 1.4142135623730951;
    }
  }

  return y;
}

/* Function for MATLAB Function: '<S6>/Differentiator' */
static void Experiment_2_moving_ca_xgeqp3_l(real_T A_data[], int32_T A_size[2],
  real_T tau_data[], int32_T *tau_size, int32_T jpvt_data[], int32_T jpvt_size[2])
{
  int32_T b_m;
  int32_T b_n;
  int32_T mn;
  real_T work_data[3];
  real_T vn1_data[3];
  real_T vn2_data[3];
  int32_T k;
  int32_T i_i;
  int32_T nmi;
  int32_T mmi;
  int32_T ix;
  real_T smax;
  real_T s;
  int32_T b_ix;
  int32_T iy;
  int32_T b_i;
  int32_T c_ix;
  int32_T l;
  int32_T b_ia;
  int32_T d_ix;
  int32_T exitg1;
  boolean_T exitg2;
  b_m = A_size[0];
  b_n = A_size[1];
  if (A_size[0] < A_size[1]) {
    mn = A_size[0];
  } else {
    mn = A_size[1];
  }

  *tau_size = (int8_T)mn;
  Expe_eml_signed_integer_colon_o(A_size[1], jpvt_data, jpvt_size);
  if (A_size[0] != 0) {
    k = (int8_T)A_size[1];
    if (0 <= k - 1) {
      memset(&work_data[0], 0, k * sizeof(real_T));
    }

    k = 1;
    for (mmi = 0; mmi < b_n; mmi++) {
      vn1_data[mmi] = Experiment_2_moving_cart_xnrm2(b_m, A_data, k);
      vn2_data[mmi] = vn1_data[mmi];
      k += b_m;
    }

    for (k = 0; k < mn; k++) {
      i_i = k * b_m + k;
      nmi = b_n - k;
      mmi = b_m - k;
      if (nmi < 1) {
        b_i = 0;
      } else {
        b_i = 1;
        if (nmi > 1) {
          ix = k;
          smax = fabs(vn1_data[k]);
          for (b_ix = 2; b_ix <= nmi; b_ix++) {
            ix++;
            s = fabs(vn1_data[ix]);
            if (s > smax) {
              b_i = b_ix;
              smax = s;
            }
          }
        }
      }

      ix = (k + b_i) - 1;
      if (ix + 1 != k + 1) {
        b_ix = b_m * ix;
        iy = b_m * k;
        for (b_i = 1; b_i <= b_m; b_i++) {
          smax = A_data[b_ix];
          A_data[b_ix] = A_data[iy];
          A_data[iy] = smax;
          b_ix++;
          iy++;
        }

        b_ix = jpvt_data[ix];
        jpvt_data[ix] = jpvt_data[k];
        jpvt_data[k] = b_ix;
        vn1_data[ix] = vn1_data[k];
        vn2_data[ix] = vn2_data[k];
      }

      if (k + 1 < b_m) {
        smax = A_data[i_i];
        tau_data[k] = 0.0;
        if (!(mmi <= 0)) {
          s = Experiment_2_moving_cart_xnrm2(mmi - 1, A_data, i_i + 2);
          if (s != 0.0) {
            s = rt_hypotd_snf(A_data[i_i], s);
            if (A_data[i_i] >= 0.0) {
              s = -s;
            }

            if (fabs(s) < 1.0020841800044864E-292) {
              ix = 0;
              b_i = i_i + mmi;
              do {
                ix++;
                for (iy = i_i + 1; iy < b_i; iy++) {
                  A_data[iy] *= 9.9792015476736E+291;
                }

                s *= 9.9792015476736E+291;
                smax *= 9.9792015476736E+291;
              } while (!(fabs(s) >= 1.0020841800044864E-292));

              s = rt_hypotd_snf(smax, Experiment_2_moving_cart_xnrm2(mmi - 1,
                A_data, i_i + 2));
              if (smax >= 0.0) {
                s = -s;
              }

              tau_data[k] = (s - smax) / s;
              smax = 1.0 / (smax - s);
              b_i = i_i + mmi;
              for (iy = i_i + 1; iy < b_i; iy++) {
                A_data[iy] *= smax;
              }

              for (b_i = 1; b_i <= ix; b_i++) {
                s *= 1.0020841800044864E-292;
              }

              smax = s;
            } else {
              tau_data[k] = (s - A_data[i_i]) / s;
              smax = 1.0 / (A_data[i_i] - s);
              ix = i_i + mmi;
              for (b_ix = i_i + 1; b_ix < ix; b_ix++) {
                A_data[b_ix] *= smax;
              }

              smax = s;
            }
          }
        }

        A_data[i_i] = smax;
      } else {
        tau_data[k] = 0.0;
      }

      if (k + 1 < b_n) {
        smax = A_data[i_i];
        A_data[i_i] = 1.0;
        iy = ((k + 1) * b_m + k) + 1;
        if (tau_data[k] != 0.0) {
          ix = mmi;
          b_i = (i_i + mmi) - 1;
          while ((ix > 0) && (A_data[b_i] == 0.0)) {
            ix--;
            b_i--;
          }

          nmi--;
          exitg2 = false;
          while ((!exitg2) && (nmi > 0)) {
            b_i = (nmi - 1) * b_m + iy;
            b_ix = b_i;
            do {
              exitg1 = 0;
              if (b_ix <= (b_i + ix) - 1) {
                if (A_data[b_ix - 1] != 0.0) {
                  exitg1 = 1;
                } else {
                  b_ix++;
                }
              } else {
                nmi--;
                exitg1 = 2;
              }
            } while (exitg1 == 0);

            if (exitg1 == 1) {
              exitg2 = true;
            }
          }
        } else {
          ix = 0;
          nmi = 0;
        }

        if (ix > 0) {
          if (nmi != 0) {
            for (b_i = 1; b_i <= nmi; b_i++) {
              work_data[b_i - 1] = 0.0;
            }

            b_i = 0;
            b_ix = (nmi - 1) * b_m + iy;
            d_ix = iy;
            while ((b_m > 0) && (d_ix <= b_ix)) {
              c_ix = i_i;
              s = 0.0;
              l = (d_ix + ix) - 1;
              for (b_ia = d_ix; b_ia <= l; b_ia++) {
                s += A_data[b_ia - 1] * A_data[c_ix];
                c_ix++;
              }

              work_data[b_i] += s;
              b_i++;
              d_ix += b_m;
            }
          }

          if (!(-tau_data[k] == 0.0)) {
            iy--;
            b_i = 0;
            for (b_ix = 1; b_ix <= nmi; b_ix++) {
              if (work_data[b_i] != 0.0) {
                s = work_data[b_i] * -tau_data[k];
                d_ix = i_i;
                c_ix = ix + iy;
                for (l = iy; l < c_ix; l++) {
                  A_data[l] += A_data[d_ix] * s;
                  d_ix++;
                }
              }

              b_i++;
              iy += b_m;
            }
          }
        }

        A_data[i_i] = smax;
      }

      for (i_i = k + 1; i_i < b_n; i_i++) {
        if (vn1_data[i_i] != 0.0) {
          smax = fabs(A_data[A_size[0] * i_i + k]) / vn1_data[i_i];
          smax = 1.0 - smax * smax;
          if (smax < 0.0) {
            smax = 0.0;
          }

          s = vn1_data[i_i] / vn2_data[i_i];
          s = s * s * smax;
          if (s <= 1.4901161193847656E-8) {
            if (k + 1 < b_m) {
              vn1_data[i_i] = Experiment_2_moving_cart_xnrm2(mmi - 1, A_data,
                (b_m * i_i + k) + 2);
              vn2_data[i_i] = vn1_data[i_i];
            } else {
              vn1_data[i_i] = 0.0;
              vn2_data[i_i] = 0.0;
            }
          } else {
            vn1_data[i_i] *= sqrt(smax);
          }
        }
      }
    }
  }
}

/* Function for MATLAB Function: '<S6>/Differentiator' */
static void Experiment_2_moving_c_qrsolve_g(const real_T A_data[], const int32_T
  A_size[2], const real_T B_data[], const int32_T B_size[2], real_T Y_data[],
  int32_T Y_size[2])
{
  real_T b_A_data[9];
  real_T tau_data[3];
  int32_T jpvt_data[3];
  int32_T rankR;
  int32_T minmn;
  int32_T maxmn;
  real_T b_B_data[9];
  int32_T b_i;
  real_T wj;
  int32_T b_j;
  int32_T b_A_size[2];
  int32_T jpvt_size[2];
  int32_T b_B_size_idx_0;
  int8_T b_idx_0;
  int8_T b_idx_1;
  int32_T Y_data_tmp;
  int32_T wj_tmp;
  b_A_size[0] = A_size[0];
  b_A_size[1] = A_size[1];
  minmn = A_size[0] * A_size[1] - 1;
  if (0 <= minmn) {
    memcpy(&b_A_data[0], &A_data[0], (minmn + 1) * sizeof(real_T));
  }

  Experiment_2_moving_ca_xgeqp3_l(b_A_data, b_A_size, tau_data, &minmn,
    jpvt_data, jpvt_size);
  rankR = 0;
  if (b_A_size[0] < b_A_size[1]) {
    minmn = b_A_size[0];
    maxmn = b_A_size[1];
  } else {
    minmn = b_A_size[1];
    maxmn = b_A_size[0];
  }

  if (minmn > 0) {
    while ((rankR < minmn) && (!(fabs(b_A_data[b_A_size[0] * rankR + rankR]) <=
             (real_T)maxmn * fabs(b_A_data[0]) * 2.2204460492503131E-16))) {
      rankR++;
    }
  }

  b_idx_0 = (int8_T)b_A_size[1];
  b_idx_1 = (int8_T)B_size[1];
  Y_size[0] = b_idx_0;
  Y_size[1] = b_idx_1;
  minmn = b_idx_0 * b_idx_1 - 1;
  if (0 <= minmn) {
    memset(&Y_data[0], 0, (minmn + 1) * sizeof(real_T));
  }

  b_B_size_idx_0 = B_size[0];
  minmn = B_size[0] * B_size[1] - 1;
  if (0 <= minmn) {
    memcpy(&b_B_data[0], &B_data[0], (minmn + 1) * sizeof(real_T));
  }

  minmn = b_A_size[0];
  if (b_A_size[0] < b_A_size[1]) {
    maxmn = b_A_size[0];
  } else {
    maxmn = b_A_size[1];
  }

  for (b_j = 0; b_j < maxmn; b_j++) {
    if (tau_data[b_j] != 0.0) {
      for (Y_data_tmp = 0; Y_data_tmp < B_size[1]; Y_data_tmp++) {
        wj_tmp = b_B_size_idx_0 * Y_data_tmp;
        wj = b_B_data[wj_tmp + b_j];
        for (b_i = b_j + 1; b_i < minmn; b_i++) {
          wj += b_A_data[b_A_size[0] * b_j + b_i] * b_B_data[wj_tmp + b_i];
        }

        wj *= tau_data[b_j];
        if (wj != 0.0) {
          b_B_data[b_j + wj_tmp] = b_B_data[b_B_size_idx_0 * Y_data_tmp + b_j] -
            wj;
          for (b_i = b_j + 1; b_i < minmn; b_i++) {
            b_B_data[b_i + wj_tmp] -= b_A_data[b_A_size[0] * b_j + b_i] * wj;
          }
        }
      }
    }
  }

  for (minmn = 0; minmn < B_size[1]; minmn++) {
    for (maxmn = 0; maxmn < rankR; maxmn++) {
      Y_data[(jpvt_data[maxmn] + b_idx_0 * minmn) - 1] = b_B_data[b_B_size_idx_0
        * minmn + maxmn];
    }

    for (maxmn = rankR - 1; maxmn + 1 > 0; maxmn--) {
      b_j = b_idx_0 * minmn;
      Y_data_tmp = b_A_size[0] * maxmn;
      Y_data[(jpvt_data[maxmn] + b_j) - 1] /= b_A_data[Y_data_tmp + maxmn];
      for (b_i = 0; b_i < maxmn; b_i++) {
        Y_data[(jpvt_data[b_i] + b_j) - 1] -= Y_data[(b_idx_0 * minmn +
          jpvt_data[maxmn]) - 1] * b_A_data[Y_data_tmp + b_i];
      }
    }
  }
}

/* Function for MATLAB Function: '<S6>/Differentiator' */
static void Experiment_2_moving_c_xzgetrf_b(int32_T b_m, int32_T b_n, real_T
  A_data[], int32_T A_size[2], int32_T lda, int32_T ipiv_data[], int32_T
  ipiv_size[2], int32_T *info)
{
  int32_T mmj;
  int32_T j;
  int32_T b_c;
  int32_T ix;
  real_T smax;
  real_T s;
  int32_T iy;
  int32_T jA;
  int32_T c_ix;
  int32_T b_j;
  int32_T e;
  int32_T ijA;
  int32_T b_m_0;
  if (b_m < b_n) {
    b_m_0 = b_m;
  } else {
    b_m_0 = b_n;
  }

  Expe_eml_signed_integer_colon_o(b_m_0, ipiv_data, ipiv_size);
  *info = 0;
  b_m_0 = b_m - 1;
  if (!(b_m_0 < b_n)) {
    b_m_0 = b_n;
  }

  for (j = 1; j <= b_m_0; j++) {
    mmj = b_m - j;
    b_c = (j - 1) * (lda + 1);
    jA = 0;
    if (mmj + 1 > 1) {
      ix = b_c;
      smax = fabs(A_data[b_c]);
      for (iy = 1; iy < mmj + 1; iy++) {
        ix++;
        s = fabs(A_data[ix]);
        if (s > smax) {
          jA = iy;
          smax = s;
        }
      }
    }

    if (A_data[b_c + jA] != 0.0) {
      if (jA != 0) {
        iy = j + jA;
        ipiv_data[j - 1] = iy;
        ix = j - 1;
        iy--;
        for (jA = 1; jA <= b_n; jA++) {
          smax = A_data[ix];
          A_data[ix] = A_data[iy];
          A_data[iy] = smax;
          ix += lda;
          iy += lda;
        }
      }

      iy = (b_c + mmj) + 1;
      for (jA = b_c + 1; jA < iy; jA++) {
        A_data[jA] /= A_data[b_c];
      }
    } else {
      *info = j;
    }

    iy = b_n - j;
    ix = b_c + lda;
    jA = ix + 1;
    for (b_j = 1; b_j <= iy; b_j++) {
      smax = A_data[ix];
      if (A_data[ix] != 0.0) {
        c_ix = b_c + 1;
        e = mmj + jA;
        for (ijA = jA; ijA < e; ijA++) {
          A_data[ijA] += A_data[c_ix] * -smax;
          c_ix++;
        }
      }

      ix += lda;
      jA += lda;
    }
  }

  if ((*info == 0) && (b_m <= b_n) && (!(A_data[((b_m - 1) * A_size[0] + b_m) -
        1] != 0.0))) {
    *info = b_m;
  }
}

/* Function for MATLAB Function: '<S6>/Differentiator' */
static void Experiment_2_moving_car_xtrsm_b(const real_T A_data[], real_T
  B_data[], int32_T B_size[2])
{
  int32_T jBcol;
  int32_T j;
  int32_T i;
  int32_T tmp;
  if (B_size[0] != 0) {
    for (j = 0; j < 3; j++) {
      jBcol = 3 * j;
      if (B_data[2 + jBcol] != 0.0) {
        B_data[2 + jBcol] /= A_data[8];
        for (i = 0; i < 2; i++) {
          tmp = i + jBcol;
          B_data[tmp] -= B_data[2 + jBcol] * A_data[i + 6];
        }
      }

      if (B_data[1 + jBcol] != 0.0) {
        B_data[1 + jBcol] /= A_data[4];
        for (i = 0; i < 1; i++) {
          tmp = i + jBcol;
          B_data[tmp] -= B_data[1 + jBcol] * A_data[i + 3];
        }
      }

      if (B_data[jBcol] != 0.0) {
        B_data[jBcol] /= A_data[0];
      }
    }
  }
}

/* Function for MATLAB Function: '<S6>/Differentiator' */
static void Experiment_2_moving__mldivide_g(const real_T A_data[], const int32_T
  A_size[2], const real_T B_data[], const int32_T B_size[2], real_T Y_data[],
  int32_T Y_size[2])
{
  real_T temp;
  int32_T ip;
  real_T b_A_data[9];
  int32_T ipiv_data[3];
  int32_T loop_ub;
  int32_T b_A_size[2];
  int32_T ipiv_size[2];
  int32_T temp_tmp;
  int32_T Y_data_tmp;
  if ((A_size[0] == 0) || (B_size[0] == 0)) {
    Y_size[0] = 3;
    Y_size[1] = 3;
    memset(&Y_data[0], 0, 9U * sizeof(real_T));
  } else if (A_size[0] == 3) {
    b_A_size[0] = 3;
    b_A_size[1] = 3;
    loop_ub = 3 * A_size[1] - 1;
    if (0 <= loop_ub) {
      memcpy(&b_A_data[0], &A_data[0], (loop_ub + 1) * sizeof(real_T));
    }

    Experiment_2_moving_c_xzgetrf_b(3, 3, b_A_data, b_A_size, 3, ipiv_data,
      ipiv_size, &loop_ub);
    Y_size[0] = B_size[0];
    Y_size[1] = 3;
    loop_ub = B_size[0] * B_size[1] - 1;
    if (0 <= loop_ub) {
      memcpy(&Y_data[0], &B_data[0], (loop_ub + 1) * sizeof(real_T));
    }

    for (loop_ub = 0; loop_ub < 2; loop_ub++) {
      if (loop_ub + 1 != ipiv_data[loop_ub]) {
        ip = ipiv_data[loop_ub] - 1;
        temp = Y_data[loop_ub];
        Y_data[loop_ub] = Y_data[ip];
        Y_data[ip] = temp;
        temp_tmp = loop_ub + Y_size[0];
        temp = Y_data[temp_tmp];
        Y_data_tmp = ip + Y_size[0];
        Y_data[temp_tmp] = Y_data[Y_data_tmp];
        Y_data[Y_data_tmp] = temp;
        temp = Y_data[(Y_size[0] << 1) + loop_ub];
        Y_data[loop_ub + (Y_size[0] << 1)] = Y_data[(Y_size[0] << 1) + ip];
        Y_data[ip + (Y_size[0] << 1)] = temp;
      }
    }

    for (loop_ub = 0; loop_ub < 3; loop_ub++) {
      ip = 3 * loop_ub;
      if (Y_data[ip] != 0.0) {
        for (temp_tmp = 1; temp_tmp < 3; temp_tmp++) {
          Y_data_tmp = temp_tmp + ip;
          Y_data[Y_data_tmp] -= Y_data[ip] * b_A_data[temp_tmp];
        }
      }

      if (Y_data[1 + ip] != 0.0) {
        for (temp_tmp = 2; temp_tmp < 3; temp_tmp++) {
          Y_data_tmp = temp_tmp + ip;
          Y_data[Y_data_tmp] -= Y_data[1 + ip] * b_A_data[temp_tmp + 3];
        }
      }
    }

    Experiment_2_moving_car_xtrsm_b(b_A_data, Y_data, Y_size);
  } else {
    Experiment_2_moving_c_qrsolve_g(A_data, A_size, B_data, B_size, b_A_data,
      b_A_size);
    Y_size[0] = b_A_size[0];
    Y_size[1] = b_A_size[1];
    loop_ub = b_A_size[0] * b_A_size[1] - 1;
    if (0 <= loop_ub) {
      memcpy(&Y_data[0], &b_A_data[0], (loop_ub + 1) * sizeof(real_T));
    }
  }
}

/* Function for MATLAB Function: '<S6>/Differentiator' */
static void Expe_PadeApproximantOfDegree_kh(const real_T A_data[], const int32_T
  A_size[2], real_T F_data[], int32_T F_size[2])
{
  int32_T b_n;
  real_T U_data[9];
  real_T V_data[9];
  real_T y_data[9];
  int32_T b_m;
  int32_T coffset;
  int32_T boffset;
  int32_T aoffset;
  int32_T j;
  int32_T i;
  int32_T b_i;
  real_T b_y_data[9];
  real_T c_y_data[9];
  real_T d_y_data[9];
  real_T b_b_data[9];
  real_T e_y_data[9];
  int32_T U_size[2];
  int32_T V_size[2];
  int32_T y_size_idx_0;
  int32_T e_y_size_idx_0;
  int32_T y_data_tmp;
  int32_T d_y_data_tmp;
  b_n = A_size[0] - 1;
  if (A_size[0] == 1) {
    y_size_idx_0 = 1;
    for (j = 0; j < 1; j++) {
      for (coffset = 0; coffset < 3; coffset++) {
        y_data[coffset] = A_data[2 + coffset] * A_data[2] + (A_data[1 + coffset]
          * A_data[1] + A_data[0] * A_data[coffset]);
      }
    }
  } else {
    b_m = A_size[0];
    y_size_idx_0 = A_size[0];
    for (j = 0; j < 3; j++) {
      coffset = j * b_m;
      boffset = j * 3;
      for (i = 1; i <= b_m; i++) {
        y_data[(coffset + i) - 1] = 0.0;
      }

      for (i = 0; i < 3; i++) {
        if (A_data[boffset + i] != 0.0) {
          aoffset = i * b_m;
          for (b_i = 0; b_i < b_m; b_i++) {
            y_data_tmp = coffset + b_i;
            y_data[y_data_tmp] += A_data[boffset + i] * A_data[aoffset + b_i];
          }
        }
      }
    }
  }

  if (y_size_idx_0 == 1) {
    for (j = 0; j < 1; j++) {
      for (coffset = 0; coffset < 3; coffset++) {
        b_y_data[coffset] = 0.0;
        for (boffset = 0; boffset < 3; boffset++) {
          b_y_data[coffset] += y_data[boffset + coffset] * y_data[boffset];
        }
      }
    }

    y_size_idx_0 = 1;
    for (j = 0; j < 1; j++) {
      for (coffset = 0; coffset < 3; coffset++) {
        c_y_data[coffset] = 0.0;
        for (boffset = 0; boffset < 3; boffset++) {
          c_y_data[coffset] += y_data[boffset + coffset] * b_y_data[boffset];
        }
      }
    }
  } else {
    for (j = 0; j < 3; j++) {
      coffset = j * y_size_idx_0;
      boffset = j * 3;
      for (i = 1; i <= y_size_idx_0; i++) {
        b_y_data[(coffset + i) - 1] = 0.0;
      }

      for (i = 0; i < 3; i++) {
        if (y_data[boffset + i] != 0.0) {
          aoffset = i * y_size_idx_0;
          for (b_i = 0; b_i < y_size_idx_0; b_i++) {
            y_data_tmp = coffset + b_i;
            b_y_data[y_data_tmp] += y_data[boffset + i] * y_data[aoffset + b_i];
          }
        }
      }
    }

    for (j = 0; j < 3; j++) {
      coffset = j * y_size_idx_0;
      boffset = j * 3;
      for (i = 1; i <= y_size_idx_0; i++) {
        c_y_data[(coffset + i) - 1] = 0.0;
      }

      for (i = 0; i < 3; i++) {
        if (y_data[boffset + i] != 0.0) {
          aoffset = i * y_size_idx_0;
          for (b_i = 0; b_i < y_size_idx_0; b_i++) {
            y_data_tmp = coffset + b_i;
            c_y_data[y_data_tmp] += y_data[boffset + i] * b_y_data[aoffset + b_i];
          }
        }
      }
    }
  }

  y_data_tmp = y_size_idx_0 * 3 - 1;
  for (j = 0; j <= y_data_tmp; j++) {
    U_data[j] = (3.352212864E+10 * c_y_data[j] + 1.05594705216E+13 * b_y_data[j])
      + 1.1873537964288E+15 * y_data[j];
  }

  for (i = 0; i <= b_n; i++) {
    e_y_size_idx_0 = y_size_idx_0 * i;
    U_data[i + e_y_size_idx_0] += 3.238237626624E+16;
  }

  for (j = 0; j <= y_data_tmp; j++) {
    b_b_data[j] = (16380.0 * b_y_data[j] + c_y_data[j]) + 4.08408E+7 * y_data[j];
  }

  if (y_size_idx_0 == 1) {
    for (j = 0; j < y_size_idx_0; j++) {
      for (coffset = 0; coffset < 3; coffset++) {
        d_y_data_tmp = y_size_idx_0 * coffset;
        i = j + d_y_data_tmp;
        d_y_data[i] = 0.0;
        for (boffset = 0; boffset < 3; boffset++) {
          d_y_data[i] = c_y_data[y_size_idx_0 * boffset + j] * b_b_data[boffset
            + coffset] + d_y_data[d_y_data_tmp + j];
        }
      }
    }
  } else {
    for (j = 0; j < 3; j++) {
      coffset = j * y_size_idx_0;
      boffset = j * 3;
      for (i = 1; i <= y_size_idx_0; i++) {
        d_y_data[(coffset + i) - 1] = 0.0;
      }

      for (i = 0; i < 3; i++) {
        if (b_b_data[boffset + i] != 0.0) {
          aoffset = i * y_size_idx_0;
          for (b_i = 0; b_i < y_size_idx_0; b_i++) {
            d_y_data_tmp = coffset + b_i;
            d_y_data[d_y_data_tmp] += b_b_data[boffset + i] * c_y_data[aoffset +
              b_i];
          }
        }
      }
    }
  }

  boffset = y_size_idx_0 * 3 - 1;
  for (j = 0; j <= boffset; j++) {
    d_y_data[j] += U_data[j];
  }

  if (y_size_idx_0 == 1) {
    e_y_size_idx_0 = A_size[0];
    boffset = A_size[0];
    for (j = 0; j < boffset; j++) {
      for (coffset = 0; coffset < 3; coffset++) {
        d_y_data_tmp = e_y_size_idx_0 * coffset;
        i = j + d_y_data_tmp;
        e_y_data[i] = 0.0;
        e_y_data[i] = e_y_data[d_y_data_tmp + j] + A_data[j] * d_y_data[coffset];
        e_y_data[i] = e_y_data[e_y_size_idx_0 * coffset + j] + A_data[j +
          A_size[0]] * d_y_data[1 + coffset];
        e_y_data[i] = A_data[(A_size[0] << 1) + j] * d_y_data[2 + coffset] +
          e_y_data[e_y_size_idx_0 * coffset + j];
      }
    }
  } else {
    b_m = A_size[0];
    e_y_size_idx_0 = A_size[0];
    for (j = 0; j < 3; j++) {
      coffset = j * b_m;
      boffset = j * 3;
      for (i = 1; i <= b_m; i++) {
        e_y_data[(coffset + i) - 1] = 0.0;
      }

      for (i = 0; i < 3; i++) {
        if (d_y_data[boffset + i] != 0.0) {
          aoffset = i * b_m;
          for (b_i = 0; b_i < b_m; b_i++) {
            d_y_data_tmp = coffset + b_i;
            e_y_data[d_y_data_tmp] += d_y_data[boffset + i] * A_data[aoffset +
              b_i];
          }
        }
      }
    }
  }

  U_size[0] = e_y_size_idx_0;
  U_size[1] = 3;
  boffset = e_y_size_idx_0 * 3 - 1;
  if (0 <= boffset) {
    memcpy(&U_data[0], &e_y_data[0], (boffset + 1) * sizeof(real_T));
  }

  for (j = 0; j <= y_data_tmp; j++) {
    b_b_data[j] = (182.0 * c_y_data[j] + 960960.0 * b_y_data[j]) + 1.32324192E+9
      * y_data[j];
  }

  if (y_size_idx_0 == 1) {
    for (j = 0; j < y_size_idx_0; j++) {
      for (coffset = 0; coffset < 3; coffset++) {
        d_y_data_tmp = y_size_idx_0 * coffset;
        i = j + d_y_data_tmp;
        d_y_data[i] = 0.0;
        for (boffset = 0; boffset < 3; boffset++) {
          d_y_data[i] = c_y_data[y_size_idx_0 * boffset + j] * b_b_data[boffset
            + coffset] + d_y_data[d_y_data_tmp + j];
        }
      }
    }
  } else {
    for (j = 0; j < 3; j++) {
      coffset = j * y_size_idx_0;
      boffset = j * 3;
      for (i = 1; i <= y_size_idx_0; i++) {
        d_y_data[(coffset + i) - 1] = 0.0;
      }

      for (i = 0; i < 3; i++) {
        if (b_b_data[boffset + i] != 0.0) {
          aoffset = i * y_size_idx_0;
          for (b_i = 0; b_i < y_size_idx_0; b_i++) {
            d_y_data_tmp = coffset + b_i;
            d_y_data[d_y_data_tmp] += b_b_data[boffset + i] * c_y_data[aoffset +
              b_i];
          }
        }
      }
    }
  }

  V_size[0] = y_size_idx_0;
  V_size[1] = 3;
  boffset = y_size_idx_0 * 3 - 1;
  for (j = 0; j <= boffset; j++) {
    V_data[j] = ((6.704425728E+11 * c_y_data[j] + d_y_data[j]) +
                 1.29060195264E+14 * b_y_data[j]) + 7.7717703038976E+15 *
      y_data[j];
  }

  for (i = 0; i <= b_n; i++) {
    y_data_tmp = y_size_idx_0 * i;
    V_data[i + y_data_tmp] += 6.476475253248E+16;
  }

  y_size_idx_0 = e_y_size_idx_0 * 3;
  for (b_m = 0; b_m < y_size_idx_0; b_m++) {
    V_data[b_m] -= U_data[b_m];
    U_data[b_m] *= 2.0;
  }

  Experiment_2_moving__mldivide_g(V_data, V_size, U_data, U_size, F_data, F_size);
  for (i = 0; i <= b_n; i++) {
    y_size_idx_0 = F_size[0] * i;
    F_data[i + y_size_idx_0]++;
  }
}

/* Function for MATLAB Function: '<S6>/Differentiator' */
static void Exper_PadeApproximantOfDegree_k(const real_T A_data[], const int32_T
  A_size[2], uint8_T b_m, real_T F_data[], int32_T F_size[2])
{
  int32_T b_n;
  real_T U_data[9];
  real_T V_data[9];
  real_T b_d;
  real_T y_data[9];
  int32_T c_m;
  int32_T coffset;
  int32_T boffset;
  int32_T aoffset;
  int32_T j;
  int32_T i;
  int32_T b_i;
  real_T b_y_data[9];
  real_T c_y_data[9];
  real_T d_y_data[9];
  real_T b_b_data[9];
  int32_T U_size[2];
  int32_T V_size[2];
  int32_T y_size_idx_0;
  int32_T b_y_size_idx_0;
  int32_T y_data_tmp;
  int32_T b_y_data_tmp;
  int32_T c_y_data_tmp;
  b_n = A_size[0] - 1;
  if (A_size[0] == 1) {
    y_size_idx_0 = 1;
    for (c_m = 0; c_m < 1; c_m++) {
      for (j = 0; j < 3; j++) {
        y_data[j] = A_data[2 + j] * A_data[2] + (A_data[1 + j] * A_data[1] +
          A_data[0] * A_data[j]);
      }
    }
  } else {
    c_m = A_size[0];
    y_size_idx_0 = A_size[0];
    for (j = 0; j < 3; j++) {
      coffset = j * c_m;
      boffset = j * 3;
      for (i = 1; i <= c_m; i++) {
        y_data[(coffset + i) - 1] = 0.0;
      }

      for (i = 0; i < 3; i++) {
        if (A_data[boffset + i] != 0.0) {
          aoffset = i * c_m;
          for (b_i = 0; b_i < c_m; b_i++) {
            y_data_tmp = coffset + b_i;
            y_data[y_data_tmp] += A_data[boffset + i] * A_data[aoffset + b_i];
          }
        }
      }
    }
  }

  if (b_m == 3) {
    y_data_tmp = y_size_idx_0 * 3 - 1;
    if (0 <= y_data_tmp) {
      memcpy(&U_data[0], &y_data[0], (y_data_tmp + 1) * sizeof(real_T));
    }

    for (i = 0; i <= b_n; i++) {
      c_m = y_size_idx_0 * i;
      U_data[i + c_m] += 60.0;
    }

    if (y_size_idx_0 == 1) {
      b_y_size_idx_0 = A_size[0];
      coffset = A_size[0];
      for (c_m = 0; c_m < coffset; c_m++) {
        for (j = 0; j < 3; j++) {
          b_y_data_tmp = b_y_size_idx_0 * j;
          boffset = c_m + b_y_data_tmp;
          b_y_data[boffset] = 0.0;
          b_y_data[boffset] = b_y_data[b_y_data_tmp + c_m] + A_data[c_m] *
            U_data[j];
          b_y_data[boffset] = b_y_data[b_y_size_idx_0 * j + c_m] + A_data[c_m +
            A_size[0]] * U_data[1 + j];
          b_y_data[boffset] = A_data[(A_size[0] << 1) + c_m] * U_data[2 + j] +
            b_y_data[b_y_size_idx_0 * j + c_m];
        }
      }
    } else {
      c_m = A_size[0];
      b_y_size_idx_0 = A_size[0];
      for (j = 0; j < 3; j++) {
        coffset = j * c_m;
        boffset = j * 3;
        for (i = 1; i <= c_m; i++) {
          b_y_data[(coffset + i) - 1] = 0.0;
        }

        for (i = 0; i < 3; i++) {
          if (U_data[boffset + i] != 0.0) {
            aoffset = i * c_m;
            for (b_i = 0; b_i < c_m; b_i++) {
              b_y_data_tmp = coffset + b_i;
              b_y_data[b_y_data_tmp] += U_data[boffset + i] * A_data[aoffset +
                b_i];
            }
          }
        }
      }
    }

    U_size[0] = b_y_size_idx_0;
    U_size[1] = 3;
    coffset = b_y_size_idx_0 * 3 - 1;
    if (0 <= coffset) {
      memcpy(&U_data[0], &b_y_data[0], (coffset + 1) * sizeof(real_T));
    }

    V_size[0] = y_size_idx_0;
    V_size[1] = 3;
    for (c_m = 0; c_m <= y_data_tmp; c_m++) {
      V_data[c_m] = 12.0 * y_data[c_m];
    }

    b_d = 120.0;
  } else {
    if (y_size_idx_0 == 1) {
      b_y_size_idx_0 = 1;
      for (c_m = 0; c_m < 1; c_m++) {
        for (j = 0; j < 3; j++) {
          b_y_data[j] = 0.0;
          for (coffset = 0; coffset < 3; coffset++) {
            b_y_data[j] += y_data[coffset + j] * y_data[coffset];
          }
        }
      }
    } else {
      b_y_size_idx_0 = y_size_idx_0;
      for (j = 0; j < 3; j++) {
        coffset = j * y_size_idx_0;
        boffset = j * 3;
        for (i = 1; i <= y_size_idx_0; i++) {
          b_y_data[(coffset + i) - 1] = 0.0;
        }

        for (i = 0; i < 3; i++) {
          if (y_data[boffset + i] != 0.0) {
            aoffset = i * y_size_idx_0;
            for (b_i = 0; b_i < y_size_idx_0; b_i++) {
              b_y_data_tmp = coffset + b_i;
              b_y_data[b_y_data_tmp] += y_data[boffset + i] * y_data[aoffset +
                b_i];
            }
          }
        }
      }
    }

    if (b_m == 5) {
      y_data_tmp = b_y_size_idx_0 * 3 - 1;
      for (c_m = 0; c_m <= y_data_tmp; c_m++) {
        U_data[c_m] = 420.0 * y_data[c_m] + b_y_data[c_m];
      }

      for (i = 0; i <= b_n; i++) {
        c_m = b_y_size_idx_0 * i;
        U_data[i + c_m] += 15120.0;
      }

      if (b_y_size_idx_0 == 1) {
        b_y_data_tmp = A_size[0];
        coffset = A_size[0];
        for (c_m = 0; c_m < coffset; c_m++) {
          for (j = 0; j < 3; j++) {
            c_y_data_tmp = b_y_data_tmp * j;
            boffset = c_m + c_y_data_tmp;
            c_y_data[boffset] = 0.0;
            c_y_data[boffset] = c_y_data[c_y_data_tmp + c_m] + A_data[c_m] *
              U_data[j];
            c_y_data[boffset] = c_y_data[b_y_data_tmp * j + c_m] + A_data[c_m +
              A_size[0]] * U_data[1 + j];
            c_y_data[boffset] = A_data[(A_size[0] << 1) + c_m] * U_data[2 + j] +
              c_y_data[b_y_data_tmp * j + c_m];
          }
        }
      } else {
        c_m = A_size[0];
        b_y_data_tmp = A_size[0];
        for (j = 0; j < 3; j++) {
          coffset = j * c_m;
          boffset = j * 3;
          for (i = 1; i <= c_m; i++) {
            c_y_data[(coffset + i) - 1] = 0.0;
          }

          for (i = 0; i < 3; i++) {
            if (U_data[boffset + i] != 0.0) {
              aoffset = i * c_m;
              for (b_i = 0; b_i < c_m; b_i++) {
                c_y_data_tmp = coffset + b_i;
                c_y_data[c_y_data_tmp] += U_data[boffset + i] * A_data[aoffset +
                  b_i];
              }
            }
          }
        }
      }

      U_size[0] = b_y_data_tmp;
      U_size[1] = 3;
      coffset = b_y_data_tmp * 3 - 1;
      if (0 <= coffset) {
        memcpy(&U_data[0], &c_y_data[0], (coffset + 1) * sizeof(real_T));
      }

      V_size[0] = b_y_size_idx_0;
      V_size[1] = 3;
      for (c_m = 0; c_m <= y_data_tmp; c_m++) {
        V_data[c_m] = 30.0 * b_y_data[c_m] + 3360.0 * y_data[c_m];
      }

      b_d = 30240.0;
    } else {
      if (y_size_idx_0 == 1) {
        b_y_data_tmp = b_y_size_idx_0;
        for (c_m = 0; c_m < b_y_size_idx_0; c_m++) {
          for (j = 0; j < 3; j++) {
            c_y_data_tmp = b_y_size_idx_0 * j;
            boffset = c_m + c_y_data_tmp;
            c_y_data[boffset] = 0.0;
            for (coffset = 0; coffset < 3; coffset++) {
              c_y_data[boffset] = b_y_data[b_y_size_idx_0 * coffset + c_m] *
                y_data[coffset + j] + c_y_data[c_y_data_tmp + c_m];
            }
          }
        }
      } else {
        b_y_data_tmp = b_y_size_idx_0;
        for (j = 0; j < 3; j++) {
          coffset = j * b_y_size_idx_0;
          boffset = j * 3;
          for (i = 1; i <= b_y_size_idx_0; i++) {
            c_y_data[(coffset + i) - 1] = 0.0;
          }

          for (i = 0; i < 3; i++) {
            if (y_data[boffset + i] != 0.0) {
              aoffset = i * b_y_size_idx_0;
              for (b_i = 0; b_i < b_y_size_idx_0; b_i++) {
                c_y_data_tmp = coffset + b_i;
                c_y_data[c_y_data_tmp] += y_data[boffset + i] * b_y_data[aoffset
                  + b_i];
              }
            }
          }
        }
      }

      if (b_m == 7) {
        coffset = b_y_data_tmp * 3 - 1;
        for (c_m = 0; c_m <= coffset; c_m++) {
          U_data[c_m] = (1512.0 * b_y_data[c_m] + c_y_data[c_m]) + 277200.0 *
            y_data[c_m];
        }

        for (i = 0; i <= b_n; i++) {
          c_m = b_y_data_tmp * i;
          U_data[i + c_m] += 8.64864E+6;
        }

        if (b_y_data_tmp == 1) {
          y_size_idx_0 = A_size[0];
          coffset = A_size[0];
          for (c_m = 0; c_m < coffset; c_m++) {
            for (j = 0; j < 3; j++) {
              b_y_size_idx_0 = y_size_idx_0 * j;
              boffset = c_m + b_y_size_idx_0;
              d_y_data[boffset] = 0.0;
              d_y_data[boffset] = d_y_data[b_y_size_idx_0 + c_m] + A_data[c_m] *
                U_data[j];
              d_y_data[boffset] = d_y_data[y_size_idx_0 * j + c_m] + A_data[c_m
                + A_size[0]] * U_data[1 + j];
              d_y_data[boffset] = A_data[(A_size[0] << 1) + c_m] * U_data[2 + j]
                + d_y_data[y_size_idx_0 * j + c_m];
            }
          }
        } else {
          c_m = A_size[0];
          y_size_idx_0 = A_size[0];
          for (j = 0; j < 3; j++) {
            coffset = j * c_m;
            boffset = j * 3;
            for (i = 1; i <= c_m; i++) {
              d_y_data[(coffset + i) - 1] = 0.0;
            }

            for (i = 0; i < 3; i++) {
              if (U_data[boffset + i] != 0.0) {
                aoffset = i * c_m;
                for (b_i = 0; b_i < c_m; b_i++) {
                  b_y_size_idx_0 = coffset + b_i;
                  d_y_data[b_y_size_idx_0] += U_data[boffset + i] *
                    A_data[aoffset + b_i];
                }
              }
            }
          }
        }

        U_size[0] = y_size_idx_0;
        U_size[1] = 3;
        coffset = y_size_idx_0 * 3 - 1;
        if (0 <= coffset) {
          memcpy(&U_data[0], &d_y_data[0], (coffset + 1) * sizeof(real_T));
        }

        V_size[0] = b_y_data_tmp;
        V_size[1] = 3;
        coffset = b_y_data_tmp * 3 - 1;
        for (c_m = 0; c_m <= coffset; c_m++) {
          V_data[c_m] = (56.0 * c_y_data[c_m] + 25200.0 * b_y_data[c_m]) +
            1.99584E+6 * y_data[c_m];
        }

        b_d = 1.729728E+7;
      } else if (b_m == 9) {
        if (y_size_idx_0 == 1) {
          y_size_idx_0 = b_y_data_tmp;
          for (c_m = 0; c_m < b_y_data_tmp; c_m++) {
            for (j = 0; j < 3; j++) {
              b_y_size_idx_0 = b_y_data_tmp * j;
              boffset = c_m + b_y_size_idx_0;
              d_y_data[boffset] = 0.0;
              for (coffset = 0; coffset < 3; coffset++) {
                d_y_data[boffset] = c_y_data[b_y_data_tmp * coffset + c_m] *
                  y_data[coffset + j] + d_y_data[b_y_size_idx_0 + c_m];
              }
            }
          }
        } else {
          y_size_idx_0 = b_y_data_tmp;
          for (j = 0; j < 3; j++) {
            coffset = j * b_y_data_tmp;
            boffset = j * 3;
            for (i = 1; i <= b_y_data_tmp; i++) {
              d_y_data[(coffset + i) - 1] = 0.0;
            }

            for (i = 0; i < 3; i++) {
              if (y_data[boffset + i] != 0.0) {
                aoffset = i * b_y_data_tmp;
                for (b_i = 0; b_i < b_y_data_tmp; b_i++) {
                  b_y_size_idx_0 = coffset + b_i;
                  d_y_data[b_y_size_idx_0] += y_data[boffset + i] *
                    c_y_data[aoffset + b_i];
                }
              }
            }
          }
        }

        y_data_tmp = y_size_idx_0 * 3 - 1;
        for (c_m = 0; c_m <= y_data_tmp; c_m++) {
          U_data[c_m] = ((3960.0 * c_y_data[c_m] + d_y_data[c_m]) + 2.16216E+6 *
                         b_y_data[c_m]) + 3.027024E+8 * y_data[c_m];
        }

        for (i = 0; i <= b_n; i++) {
          c_m = y_size_idx_0 * i;
          U_data[i + c_m] += 8.8216128E+9;
        }

        if (y_size_idx_0 == 1) {
          b_y_size_idx_0 = A_size[0];
          coffset = A_size[0];
          for (c_m = 0; c_m < coffset; c_m++) {
            for (j = 0; j < 3; j++) {
              c_y_data_tmp = b_y_size_idx_0 * j;
              boffset = c_m + c_y_data_tmp;
              b_b_data[boffset] = 0.0;
              b_b_data[boffset] = b_b_data[c_y_data_tmp + c_m] + A_data[c_m] *
                U_data[j];
              b_b_data[boffset] = b_b_data[b_y_size_idx_0 * j + c_m] +
                A_data[c_m + A_size[0]] * U_data[1 + j];
              b_b_data[boffset] = A_data[(A_size[0] << 1) + c_m] * U_data[2 + j]
                + b_b_data[b_y_size_idx_0 * j + c_m];
            }
          }
        } else {
          c_m = A_size[0];
          b_y_size_idx_0 = A_size[0];
          for (j = 0; j < 3; j++) {
            coffset = j * c_m;
            boffset = j * 3;
            for (i = 1; i <= c_m; i++) {
              b_b_data[(coffset + i) - 1] = 0.0;
            }

            for (i = 0; i < 3; i++) {
              if (U_data[boffset + i] != 0.0) {
                aoffset = i * c_m;
                for (b_i = 0; b_i < c_m; b_i++) {
                  c_y_data_tmp = coffset + b_i;
                  b_b_data[c_y_data_tmp] += U_data[boffset + i] * A_data[aoffset
                    + b_i];
                }
              }
            }
          }
        }

        U_size[0] = b_y_size_idx_0;
        U_size[1] = 3;
        coffset = b_y_size_idx_0 * 3 - 1;
        if (0 <= coffset) {
          memcpy(&U_data[0], &b_b_data[0], (coffset + 1) * sizeof(real_T));
        }

        V_size[0] = y_size_idx_0;
        V_size[1] = 3;
        for (c_m = 0; c_m <= y_data_tmp; c_m++) {
          V_data[c_m] = ((90.0 * d_y_data[c_m] + 110880.0 * c_y_data[c_m]) +
                         3.027024E+7 * b_y_data[c_m]) + 2.0756736E+9 *
            y_data[c_m];
        }

        b_d = 1.76432256E+10;
      } else {
        y_data_tmp = b_y_data_tmp * 3 - 1;
        for (c_m = 0; c_m <= y_data_tmp; c_m++) {
          U_data[c_m] = (3.352212864E+10 * c_y_data[c_m] + 1.05594705216E+13 *
                         b_y_data[c_m]) + 1.1873537964288E+15 * y_data[c_m];
        }

        for (i = 0; i <= b_n; i++) {
          c_m = b_y_data_tmp * i;
          U_data[i + c_m] += 3.238237626624E+16;
        }

        for (c_m = 0; c_m <= y_data_tmp; c_m++) {
          b_b_data[c_m] = (16380.0 * b_y_data[c_m] + c_y_data[c_m]) + 4.08408E+7
            * y_data[c_m];
        }

        if (b_y_data_tmp == 1) {
          y_size_idx_0 = b_y_data_tmp;
          for (c_m = 0; c_m < b_y_data_tmp; c_m++) {
            for (j = 0; j < 3; j++) {
              b_y_size_idx_0 = b_y_data_tmp * j;
              boffset = c_m + b_y_size_idx_0;
              d_y_data[boffset] = 0.0;
              for (coffset = 0; coffset < 3; coffset++) {
                d_y_data[boffset] = c_y_data[b_y_data_tmp * coffset + c_m] *
                  b_b_data[coffset + j] + d_y_data[b_y_size_idx_0 + c_m];
              }
            }
          }
        } else {
          y_size_idx_0 = b_y_data_tmp;
          for (j = 0; j < 3; j++) {
            coffset = j * b_y_data_tmp;
            boffset = j * 3;
            for (i = 1; i <= b_y_data_tmp; i++) {
              d_y_data[(coffset + i) - 1] = 0.0;
            }

            for (i = 0; i < 3; i++) {
              if (b_b_data[boffset + i] != 0.0) {
                aoffset = i * b_y_data_tmp;
                for (b_i = 0; b_i < b_y_data_tmp; b_i++) {
                  b_y_size_idx_0 = coffset + b_i;
                  d_y_data[b_y_size_idx_0] += b_b_data[boffset + i] *
                    c_y_data[aoffset + b_i];
                }
              }
            }
          }
        }

        coffset = y_size_idx_0 * 3 - 1;
        for (c_m = 0; c_m <= coffset; c_m++) {
          d_y_data[c_m] += U_data[c_m];
        }

        if (y_size_idx_0 == 1) {
          b_y_size_idx_0 = A_size[0];
          coffset = A_size[0];
          for (c_m = 0; c_m < coffset; c_m++) {
            for (j = 0; j < 3; j++) {
              c_y_data_tmp = b_y_size_idx_0 * j;
              boffset = c_m + c_y_data_tmp;
              b_b_data[boffset] = 0.0;
              b_b_data[boffset] = b_b_data[c_y_data_tmp + c_m] + A_data[c_m] *
                d_y_data[j];
              b_b_data[boffset] = b_b_data[b_y_size_idx_0 * j + c_m] +
                A_data[c_m + A_size[0]] * d_y_data[1 + j];
              b_b_data[boffset] = A_data[(A_size[0] << 1) + c_m] * d_y_data[2 +
                j] + b_b_data[b_y_size_idx_0 * j + c_m];
            }
          }
        } else {
          c_m = A_size[0];
          b_y_size_idx_0 = A_size[0];
          for (j = 0; j < 3; j++) {
            coffset = j * c_m;
            boffset = j * 3;
            for (i = 1; i <= c_m; i++) {
              b_b_data[(coffset + i) - 1] = 0.0;
            }

            for (i = 0; i < 3; i++) {
              if (d_y_data[boffset + i] != 0.0) {
                aoffset = i * c_m;
                for (b_i = 0; b_i < c_m; b_i++) {
                  c_y_data_tmp = coffset + b_i;
                  b_b_data[c_y_data_tmp] += d_y_data[boffset + i] *
                    A_data[aoffset + b_i];
                }
              }
            }
          }
        }

        U_size[0] = b_y_size_idx_0;
        U_size[1] = 3;
        coffset = b_y_size_idx_0 * 3 - 1;
        if (0 <= coffset) {
          memcpy(&U_data[0], &b_b_data[0], (coffset + 1) * sizeof(real_T));
        }

        for (c_m = 0; c_m <= y_data_tmp; c_m++) {
          b_b_data[c_m] = (182.0 * c_y_data[c_m] + 960960.0 * b_y_data[c_m]) +
            1.32324192E+9 * y_data[c_m];
        }

        if (b_y_data_tmp == 1) {
          for (c_m = 0; c_m < b_y_data_tmp; c_m++) {
            for (j = 0; j < 3; j++) {
              b_y_size_idx_0 = b_y_data_tmp * j;
              boffset = c_m + b_y_size_idx_0;
              d_y_data[boffset] = 0.0;
              for (coffset = 0; coffset < 3; coffset++) {
                d_y_data[boffset] = c_y_data[b_y_data_tmp * coffset + c_m] *
                  b_b_data[coffset + j] + d_y_data[b_y_size_idx_0 + c_m];
              }
            }
          }
        } else {
          for (j = 0; j < 3; j++) {
            coffset = j * b_y_data_tmp;
            boffset = j * 3;
            for (i = 1; i <= b_y_data_tmp; i++) {
              d_y_data[(coffset + i) - 1] = 0.0;
            }

            for (i = 0; i < 3; i++) {
              if (b_b_data[boffset + i] != 0.0) {
                aoffset = i * b_y_data_tmp;
                for (b_i = 0; b_i < b_y_data_tmp; b_i++) {
                  b_y_size_idx_0 = coffset + b_i;
                  d_y_data[b_y_size_idx_0] += b_b_data[boffset + i] *
                    c_y_data[aoffset + b_i];
                }
              }
            }
          }
        }

        V_size[0] = y_size_idx_0;
        V_size[1] = 3;
        coffset = y_size_idx_0 * 3 - 1;
        for (c_m = 0; c_m <= coffset; c_m++) {
          V_data[c_m] = ((6.704425728E+11 * c_y_data[c_m] + d_y_data[c_m]) +
                         1.29060195264E+14 * b_y_data[c_m]) +
            7.7717703038976E+15 * y_data[c_m];
        }

        b_d = 6.476475253248E+16;
      }
    }
  }

  for (i = 0; i <= b_n; i++) {
    y_data_tmp = V_size[0] * i;
    V_data[i + y_data_tmp] += b_d;
  }

  y_data_tmp = U_size[0] * 3;
  for (c_m = 0; c_m < y_data_tmp; c_m++) {
    V_data[c_m] -= U_data[c_m];
    U_data[c_m] *= 2.0;
  }

  Experiment_2_moving__mldivide_g(V_data, V_size, U_data, U_size, F_data, F_size);
  for (i = 0; i <= b_n; i++) {
    y_data_tmp = F_size[0] * i;
    F_data[i + y_data_tmp]++;
  }
}

/* Function for MATLAB Function: '<S6>/Differentiator' */
static void Experiment_2_moving__mrdivide_f(const real_T A_data[], const int32_T
  A_size[2], const real_T B_data[], const int32_T B_size[2], real_T y_data[],
  int32_T y_size[2])
{
  int32_T b_n;
  int32_T nb;
  real_T temp;
  int32_T xi;
  real_T b_A_data[9];
  int32_T ipiv_data[3];
  int32_T j;
  int32_T b_jAcol;
  int32_T b_jBcol;
  real_T B_data_0[4];
  real_T A_data_0[4];
  int32_T loop_ub;
  int32_T b_A_size[2];
  int32_T B_size_0[2];
  int32_T A_size_0[2];
  int32_T X_data_tmp;
  if (B_size[0] == B_size[1]) {
    b_n = B_size[1];
    b_A_size[0] = B_size[0];
    b_A_size[1] = B_size[1];
    nb = B_size[0] * B_size[1] - 1;
    if (0 <= nb) {
      memcpy(&b_A_data[0], &B_data[0], (nb + 1) * sizeof(real_T));
    }

    Experiment_2_moving_c_xzgetrf_b(B_size[1], B_size[1], b_A_data, b_A_size,
      B_size[1], ipiv_data, B_size_0, &nb);
    nb = A_size[0];
    loop_ub = A_size[0] * A_size[1] - 1;
    if (0 <= loop_ub) {
      memcpy(&B_data_0[0], &A_data[0], (loop_ub + 1) * sizeof(real_T));
    }

    for (j = 0; j < b_n; j++) {
      b_jBcol = nb * j;
      b_jAcol = b_n * j;
      xi = 1;
      while (xi <= j) {
        if (b_A_data[b_jAcol] != 0.0) {
          for (xi = 0; xi < nb; xi++) {
            X_data_tmp = xi + b_jBcol;
            B_data_0[X_data_tmp] -= b_A_data[b_jAcol] * B_data_0[xi];
          }
        }

        xi = 2;
      }

      temp = 1.0 / b_A_data[j + b_jAcol];
      for (b_jAcol = 0; b_jAcol < nb; b_jAcol++) {
        X_data_tmp = b_jAcol + b_jBcol;
        B_data_0[X_data_tmp] *= temp;
      }
    }

    for (j = B_size[1]; j > 0; j--) {
      b_jBcol = (j - 1) * nb;
      b_jAcol = (j - 1) * b_n + 1;
      xi = j + 1;
      while (xi <= b_n) {
        if (b_A_data[b_jAcol] != 0.0) {
          for (xi = 0; xi < nb; xi++) {
            X_data_tmp = xi + b_jBcol;
            B_data_0[X_data_tmp] -= B_data_0[xi + nb] * b_A_data[b_jAcol];
          }
        }

        xi = 3;
      }
    }

    b_n = B_size[1] - 1;
    while (b_n > 0) {
      if (ipiv_data[0] != 1) {
        b_n = ipiv_data[0] - 1;
        for (xi = 0; xi < nb; xi++) {
          temp = B_data_0[xi];
          X_data_tmp = A_size[0] * b_n;
          B_data_0[xi] = B_data_0[X_data_tmp + xi];
          B_data_0[xi + X_data_tmp] = temp;
        }
      }

      b_n = 0;
    }

    y_size[0] = A_size[0];
    y_size[1] = A_size[1];
    if (0 <= loop_ub) {
      memcpy(&y_data[0], &B_data_0[0], (loop_ub + 1) * sizeof(real_T));
    }
  } else {
    B_size_0[0] = B_size[1];
    B_size_0[1] = B_size[0];
    nb = B_size[0];
    for (b_n = 0; b_n < nb; b_n++) {
      loop_ub = B_size[1];
      for (j = 0; j < loop_ub; j++) {
        B_data_0[j + B_size_0[0] * b_n] = B_data[B_size[0] * j + b_n];
      }
    }

    A_size_0[0] = A_size[1];
    A_size_0[1] = A_size[0];
    nb = A_size[0];
    for (b_n = 0; b_n < nb; b_n++) {
      loop_ub = A_size[1];
      for (j = 0; j < loop_ub; j++) {
        A_data_0[j + A_size_0[0] * b_n] = A_data[A_size[0] * j + b_n];
      }
    }

    Experiment_2_moving_c_qrsolve_g(B_data_0, B_size_0, A_data_0, A_size_0,
      b_A_data, b_A_size);
    y_size[0] = b_A_size[1];
    y_size[1] = b_A_size[0];
    nb = b_A_size[0];
    for (b_n = 0; b_n < nb; b_n++) {
      loop_ub = b_A_size[1];
      for (j = 0; j < loop_ub; j++) {
        y_data[j + y_size[0] * b_n] = b_A_data[b_A_size[0] * j + b_n];
      }
    }
  }
}

/* Function for MATLAB Function: '<S6>/Differentiator' */
static real_T Experim_disc_eigenvalues_URED_d(real_T x0, real_T b_mu, real_T tau)
{
  real_T z;
  real_T c;
  c = rt_powd_snf(fabs(x0), 0.5);
  z = c / ((tau * 75.0 * b_mu * fabs(x0) + c) + tau * 75.0);
  return z;
}

/* Function for MATLAB Function: '<S6>/Differentiator' */
static void Experime_ackerman_precomputed_g(real_T b_Ts, real_T z, real_T
  lambda_data[], int32_T *lambda_size)
{
  *lambda_size = 2;
  lambda_data[0] = (2.0 - z) - z;
  lambda_data[1] = (z - 1.0) * (z - 1.0) / b_Ts;
}

/* Function for MATLAB Function: '<S6>/Differentiator' */
static void Experiment_2_moving_cart_step_n(real_T x0, const real_T z_data[],
  real_T u, const real_T Phi_data[], const real_T bD_data[], const real_T
  lambda_data[], real_T z_new_data[], int32_T *z_new_size)
{
  real_T y_data[2];
  int32_T aoffset;
  int32_T i;
  int32_T b_i;
  for (i = 1; i < 3; i++) {
    y_data[i - 1] = 0.0;
  }

  for (i = 0; i < 2; i++) {
    if (z_data[i] != 0.0) {
      aoffset = i << 1;
      for (b_i = 0; b_i < 2; b_i++) {
        y_data[b_i] += Phi_data[aoffset + b_i] * z_data[i];
      }
    }
  }

  *z_new_size = 2;
  for (i = 0; i < 2; i++) {
    z_new_data[i] = (bD_data[i] * u + y_data[i]) + lambda_data[i] * x0;
  }
}

/* Function for MATLAB Function: '<S5>/Differentiator' */
static real_T Experiment_2_moving_cart_norm(const real_T x_data[], const int32_T
  x_size[2])
{
  real_T y;
  real_T s;
  int32_T j;
  int32_T i;
  boolean_T exitg1;
  if (x_size[0] == 0) {
    y = 0.0;
  } else if (x_size[0] == 1) {
    y = ((fabs(x_data[0]) + fabs(x_data[1])) + fabs(x_data[2])) + fabs(x_data[3]);
  } else {
    y = 0.0;
    j = 0;
    exitg1 = false;
    while ((!exitg1) && (j <= 3)) {
      s = 0.0;
      for (i = 0; i < x_size[0]; i++) {
        s += fabs(x_data[x_size[0] * j + i]);
      }

      if (rtIsNaN(s)) {
        y = (rtNaN);
        exitg1 = true;
      } else {
        if (s > y) {
          y = s;
        }

        j++;
      }
    }
  }

  return y;
}

/* Function for MATLAB Function: '<S5>/Differentiator' */
static void Experi_eml_signed_integer_colon(int32_T b, int32_T y_data[], int32_T
  y_size[2])
{
  int32_T yk;
  int32_T k;
  y_size[0] = 1;
  y_size[1] = b;
  y_data[0] = 1;
  yk = 1;
  for (k = 2; k <= b; k++) {
    yk++;
    y_data[k - 1] = yk;
  }
}

/* Function for MATLAB Function: '<S5>/Differentiator' */
static void Experiment_2_moving_cart_xgeqp3(real_T A_data[], int32_T A_size[2],
  real_T tau_data[], int32_T *tau_size, int32_T jpvt_data[], int32_T jpvt_size[2])
{
  int32_T b_m;
  int32_T b_n;
  int32_T mn;
  real_T work_data[4];
  real_T vn1_data[4];
  real_T vn2_data[4];
  int32_T k;
  int32_T i_i;
  int32_T nmi;
  int32_T mmi;
  int32_T ix;
  real_T smax;
  real_T s;
  int32_T b_ix;
  int32_T iy;
  int32_T b_i;
  int32_T c_ix;
  int32_T l;
  int32_T b_ia;
  int32_T d_ix;
  int32_T exitg1;
  boolean_T exitg2;
  b_m = A_size[0];
  b_n = A_size[1];
  if (A_size[0] < A_size[1]) {
    mn = A_size[0];
  } else {
    mn = A_size[1];
  }

  *tau_size = (int8_T)mn;
  Experi_eml_signed_integer_colon(A_size[1], jpvt_data, jpvt_size);
  if (A_size[0] != 0) {
    k = (int8_T)A_size[1];
    if (0 <= k - 1) {
      memset(&work_data[0], 0, k * sizeof(real_T));
    }

    k = 1;
    for (mmi = 0; mmi < b_n; mmi++) {
      vn1_data[mmi] = Experiment_2_moving_cart_xnrm2(b_m, A_data, k);
      vn2_data[mmi] = vn1_data[mmi];
      k += b_m;
    }

    for (k = 0; k < mn; k++) {
      i_i = k * b_m + k;
      nmi = b_n - k;
      mmi = b_m - k;
      if (nmi < 1) {
        b_i = 0;
      } else {
        b_i = 1;
        if (nmi > 1) {
          ix = k;
          smax = fabs(vn1_data[k]);
          for (b_ix = 2; b_ix <= nmi; b_ix++) {
            ix++;
            s = fabs(vn1_data[ix]);
            if (s > smax) {
              b_i = b_ix;
              smax = s;
            }
          }
        }
      }

      ix = (k + b_i) - 1;
      if (ix + 1 != k + 1) {
        b_ix = b_m * ix;
        iy = b_m * k;
        for (b_i = 1; b_i <= b_m; b_i++) {
          smax = A_data[b_ix];
          A_data[b_ix] = A_data[iy];
          A_data[iy] = smax;
          b_ix++;
          iy++;
        }

        b_ix = jpvt_data[ix];
        jpvt_data[ix] = jpvt_data[k];
        jpvt_data[k] = b_ix;
        vn1_data[ix] = vn1_data[k];
        vn2_data[ix] = vn2_data[k];
      }

      if (k + 1 < b_m) {
        smax = A_data[i_i];
        tau_data[k] = 0.0;
        if (!(mmi <= 0)) {
          s = Experiment_2_moving_cart_xnrm2(mmi - 1, A_data, i_i + 2);
          if (s != 0.0) {
            s = rt_hypotd_snf(A_data[i_i], s);
            if (A_data[i_i] >= 0.0) {
              s = -s;
            }

            if (fabs(s) < 1.0020841800044864E-292) {
              ix = 0;
              b_i = i_i + mmi;
              do {
                ix++;
                for (iy = i_i + 1; iy < b_i; iy++) {
                  A_data[iy] *= 9.9792015476736E+291;
                }

                s *= 9.9792015476736E+291;
                smax *= 9.9792015476736E+291;
              } while (!(fabs(s) >= 1.0020841800044864E-292));

              s = rt_hypotd_snf(smax, Experiment_2_moving_cart_xnrm2(mmi - 1,
                A_data, i_i + 2));
              if (smax >= 0.0) {
                s = -s;
              }

              tau_data[k] = (s - smax) / s;
              smax = 1.0 / (smax - s);
              b_i = i_i + mmi;
              for (iy = i_i + 1; iy < b_i; iy++) {
                A_data[iy] *= smax;
              }

              for (b_i = 1; b_i <= ix; b_i++) {
                s *= 1.0020841800044864E-292;
              }

              smax = s;
            } else {
              tau_data[k] = (s - A_data[i_i]) / s;
              smax = 1.0 / (A_data[i_i] - s);
              ix = i_i + mmi;
              for (b_ix = i_i + 1; b_ix < ix; b_ix++) {
                A_data[b_ix] *= smax;
              }

              smax = s;
            }
          }
        }

        A_data[i_i] = smax;
      } else {
        tau_data[k] = 0.0;
      }

      if (k + 1 < b_n) {
        smax = A_data[i_i];
        A_data[i_i] = 1.0;
        iy = ((k + 1) * b_m + k) + 1;
        if (tau_data[k] != 0.0) {
          ix = mmi;
          b_i = (i_i + mmi) - 1;
          while ((ix > 0) && (A_data[b_i] == 0.0)) {
            ix--;
            b_i--;
          }

          nmi--;
          exitg2 = false;
          while ((!exitg2) && (nmi > 0)) {
            b_i = (nmi - 1) * b_m + iy;
            b_ix = b_i;
            do {
              exitg1 = 0;
              if (b_ix <= (b_i + ix) - 1) {
                if (A_data[b_ix - 1] != 0.0) {
                  exitg1 = 1;
                } else {
                  b_ix++;
                }
              } else {
                nmi--;
                exitg1 = 2;
              }
            } while (exitg1 == 0);

            if (exitg1 == 1) {
              exitg2 = true;
            }
          }
        } else {
          ix = 0;
          nmi = 0;
        }

        if (ix > 0) {
          if (nmi != 0) {
            for (b_i = 1; b_i <= nmi; b_i++) {
              work_data[b_i - 1] = 0.0;
            }

            b_i = 0;
            b_ix = (nmi - 1) * b_m + iy;
            d_ix = iy;
            while ((b_m > 0) && (d_ix <= b_ix)) {
              c_ix = i_i;
              s = 0.0;
              l = (d_ix + ix) - 1;
              for (b_ia = d_ix; b_ia <= l; b_ia++) {
                s += A_data[b_ia - 1] * A_data[c_ix];
                c_ix++;
              }

              work_data[b_i] += s;
              b_i++;
              d_ix += b_m;
            }
          }

          if (!(-tau_data[k] == 0.0)) {
            iy--;
            b_i = 0;
            for (b_ix = 1; b_ix <= nmi; b_ix++) {
              if (work_data[b_i] != 0.0) {
                s = work_data[b_i] * -tau_data[k];
                d_ix = i_i;
                c_ix = ix + iy;
                for (l = iy; l < c_ix; l++) {
                  A_data[l] += A_data[d_ix] * s;
                  d_ix++;
                }
              }

              b_i++;
              iy += b_m;
            }
          }
        }

        A_data[i_i] = smax;
      }

      for (i_i = k + 1; i_i < b_n; i_i++) {
        if (vn1_data[i_i] != 0.0) {
          smax = fabs(A_data[A_size[0] * i_i + k]) / vn1_data[i_i];
          smax = 1.0 - smax * smax;
          if (smax < 0.0) {
            smax = 0.0;
          }

          s = vn1_data[i_i] / vn2_data[i_i];
          s = s * s * smax;
          if (s <= 1.4901161193847656E-8) {
            if (k + 1 < b_m) {
              vn1_data[i_i] = Experiment_2_moving_cart_xnrm2(mmi - 1, A_data,
                (b_m * i_i + k) + 2);
              vn2_data[i_i] = vn1_data[i_i];
            } else {
              vn1_data[i_i] = 0.0;
              vn2_data[i_i] = 0.0;
            }
          } else {
            vn1_data[i_i] *= sqrt(smax);
          }
        }
      }
    }
  }
}

/* Function for MATLAB Function: '<S5>/Differentiator' */
static void Experiment_2_moving_car_qrsolve(const real_T A_data[], const int32_T
  A_size[2], const real_T B_data[], const int32_T B_size[2], real_T Y_data[],
  int32_T Y_size[2])
{
  real_T b_A_data[16];
  real_T tau_data[4];
  int32_T jpvt_data[4];
  int32_T rankR;
  int32_T minmn;
  int32_T maxmn;
  real_T b_B_data[16];
  int32_T b_i;
  real_T wj;
  int32_T b_j;
  int32_T b_A_size[2];
  int32_T jpvt_size[2];
  int32_T b_B_size_idx_0;
  int8_T b_idx_0;
  int8_T b_idx_1;
  int32_T Y_data_tmp;
  int32_T wj_tmp;
  b_A_size[0] = A_size[0];
  b_A_size[1] = A_size[1];
  minmn = A_size[0] * A_size[1] - 1;
  if (0 <= minmn) {
    memcpy(&b_A_data[0], &A_data[0], (minmn + 1) * sizeof(real_T));
  }

  Experiment_2_moving_cart_xgeqp3(b_A_data, b_A_size, tau_data, &minmn,
    jpvt_data, jpvt_size);
  rankR = 0;
  if (b_A_size[0] < b_A_size[1]) {
    minmn = b_A_size[0];
    maxmn = b_A_size[1];
  } else {
    minmn = b_A_size[1];
    maxmn = b_A_size[0];
  }

  if (minmn > 0) {
    while ((rankR < minmn) && (!(fabs(b_A_data[b_A_size[0] * rankR + rankR]) <=
             (real_T)maxmn * fabs(b_A_data[0]) * 2.2204460492503131E-16))) {
      rankR++;
    }
  }

  b_idx_0 = (int8_T)b_A_size[1];
  b_idx_1 = (int8_T)B_size[1];
  Y_size[0] = b_idx_0;
  Y_size[1] = b_idx_1;
  minmn = b_idx_0 * b_idx_1 - 1;
  if (0 <= minmn) {
    memset(&Y_data[0], 0, (minmn + 1) * sizeof(real_T));
  }

  b_B_size_idx_0 = B_size[0];
  minmn = B_size[0] * B_size[1] - 1;
  if (0 <= minmn) {
    memcpy(&b_B_data[0], &B_data[0], (minmn + 1) * sizeof(real_T));
  }

  minmn = b_A_size[0];
  if (b_A_size[0] < b_A_size[1]) {
    maxmn = b_A_size[0];
  } else {
    maxmn = b_A_size[1];
  }

  for (b_j = 0; b_j < maxmn; b_j++) {
    if (tau_data[b_j] != 0.0) {
      for (Y_data_tmp = 0; Y_data_tmp < B_size[1]; Y_data_tmp++) {
        wj_tmp = b_B_size_idx_0 * Y_data_tmp;
        wj = b_B_data[wj_tmp + b_j];
        for (b_i = b_j + 1; b_i < minmn; b_i++) {
          wj += b_A_data[b_A_size[0] * b_j + b_i] * b_B_data[wj_tmp + b_i];
        }

        wj *= tau_data[b_j];
        if (wj != 0.0) {
          b_B_data[b_j + wj_tmp] = b_B_data[b_B_size_idx_0 * Y_data_tmp + b_j] -
            wj;
          for (b_i = b_j + 1; b_i < minmn; b_i++) {
            b_B_data[b_i + wj_tmp] -= b_A_data[b_A_size[0] * b_j + b_i] * wj;
          }
        }
      }
    }
  }

  for (minmn = 0; minmn < B_size[1]; minmn++) {
    for (maxmn = 0; maxmn < rankR; maxmn++) {
      Y_data[(jpvt_data[maxmn] + b_idx_0 * minmn) - 1] = b_B_data[b_B_size_idx_0
        * minmn + maxmn];
    }

    for (maxmn = rankR - 1; maxmn + 1 > 0; maxmn--) {
      b_j = b_idx_0 * minmn;
      Y_data_tmp = b_A_size[0] * maxmn;
      Y_data[(jpvt_data[maxmn] + b_j) - 1] /= b_A_data[Y_data_tmp + maxmn];
      for (b_i = 0; b_i < maxmn; b_i++) {
        Y_data[(jpvt_data[b_i] + b_j) - 1] -= Y_data[(b_idx_0 * minmn +
          jpvt_data[maxmn]) - 1] * b_A_data[Y_data_tmp + b_i];
      }
    }
  }
}

/* Function for MATLAB Function: '<S5>/Differentiator' */
static void Experiment_2_moving_car_xzgetrf(int32_T b_m, int32_T b_n, real_T
  A_data[], int32_T A_size[2], int32_T lda, int32_T ipiv_data[], int32_T
  ipiv_size[2], int32_T *info)
{
  int32_T mmj;
  int32_T j;
  int32_T b_c;
  int32_T ix;
  real_T smax;
  real_T s;
  int32_T iy;
  int32_T jA;
  int32_T c_ix;
  int32_T b_j;
  int32_T e;
  int32_T ijA;
  int32_T b_m_0;
  if (b_m < b_n) {
    b_m_0 = b_m;
  } else {
    b_m_0 = b_n;
  }

  Experi_eml_signed_integer_colon(b_m_0, ipiv_data, ipiv_size);
  *info = 0;
  b_m_0 = b_m - 1;
  if (!(b_m_0 < b_n)) {
    b_m_0 = b_n;
  }

  for (j = 0; j < b_m_0; j++) {
    mmj = b_m - j;
    b_c = (lda + 1) * j;
    if (mmj < 1) {
      jA = -1;
    } else {
      jA = 0;
      if (mmj > 1) {
        ix = b_c;
        smax = fabs(A_data[b_c]);
        for (iy = 1; iy < mmj; iy++) {
          ix++;
          s = fabs(A_data[ix]);
          if (s > smax) {
            jA = iy;
            smax = s;
          }
        }
      }
    }

    if (A_data[b_c + jA] != 0.0) {
      if (jA != 0) {
        iy = j + jA;
        ipiv_data[j] = iy + 1;
        ix = j;
        for (jA = 1; jA <= b_n; jA++) {
          smax = A_data[ix];
          A_data[ix] = A_data[iy];
          A_data[iy] = smax;
          ix += lda;
          iy += lda;
        }
      }

      iy = b_c + mmj;
      for (jA = b_c + 1; jA < iy; jA++) {
        A_data[jA] /= A_data[b_c];
      }
    } else {
      *info = j + 1;
    }

    iy = (b_n - j) - 1;
    ix = b_c + lda;
    jA = ix + 1;
    for (b_j = 1; b_j <= iy; b_j++) {
      smax = A_data[ix];
      if (A_data[ix] != 0.0) {
        c_ix = b_c + 1;
        e = mmj + jA;
        for (ijA = jA; ijA < e - 1; ijA++) {
          A_data[ijA] += A_data[c_ix] * -smax;
          c_ix++;
        }
      }

      ix += lda;
      jA += lda;
    }
  }

  if ((*info == 0) && (b_m <= b_n) && (!(A_data[((b_m - 1) * A_size[0] + b_m) -
        1] != 0.0))) {
    *info = b_m;
  }
}

/* Function for MATLAB Function: '<S5>/Differentiator' */
static void Experiment_2_moving_cart_xtrsm(const real_T A_data[], real_T B_data[],
  int32_T B_size[2])
{
  int32_T jBcol;
  int32_T j;
  int32_T i;
  if (B_size[0] != 0) {
    for (j = 0; j < 4; j++) {
      jBcol = j << 2;
      if (B_data[3 + jBcol] != 0.0) {
        B_data[3 + jBcol] /= A_data[15];
        for (i = 0; i < 3; i++) {
          B_data[i + jBcol] -= B_data[3 + jBcol] * A_data[i + 12];
        }
      }

      if (B_data[2 + jBcol] != 0.0) {
        B_data[2 + jBcol] /= A_data[10];
        for (i = 0; i < 2; i++) {
          B_data[i + jBcol] -= B_data[2 + jBcol] * A_data[i + 8];
        }
      }

      if (B_data[1 + jBcol] != 0.0) {
        B_data[1 + jBcol] /= A_data[5];
        for (i = 0; i < 1; i++) {
          B_data[i + jBcol] -= B_data[1 + jBcol] * A_data[i + 4];
        }
      }

      if (B_data[jBcol] != 0.0) {
        B_data[jBcol] /= A_data[0];
      }
    }
  }
}

/* Function for MATLAB Function: '<S5>/Differentiator' */
static void Experiment_2_moving_ca_mldivide(const real_T A_data[], const int32_T
  A_size[2], const real_T B_data[], const int32_T B_size[2], real_T Y_data[],
  int32_T Y_size[2])
{
  real_T temp;
  int32_T ip;
  real_T b_A_data[16];
  int32_T ipiv_data[4];
  int32_T loop_ub;
  int32_T b_A_size[2];
  int32_T ipiv_size[2];
  int32_T temp_tmp;
  int32_T Y_data_tmp;
  if ((A_size[0] == 0) || (B_size[0] == 0)) {
    Y_size[0] = 4;
    Y_size[1] = 4;
    memset(&Y_data[0], 0, sizeof(real_T) << 4U);
  } else if (A_size[0] == 4) {
    b_A_size[0] = 4;
    b_A_size[1] = 4;
    loop_ub = (A_size[1] << 2) - 1;
    if (0 <= loop_ub) {
      memcpy(&b_A_data[0], &A_data[0], (loop_ub + 1) * sizeof(real_T));
    }

    Experiment_2_moving_car_xzgetrf(4, 4, b_A_data, b_A_size, 4, ipiv_data,
      ipiv_size, &loop_ub);
    Y_size[0] = B_size[0];
    Y_size[1] = 4;
    loop_ub = B_size[0] * B_size[1] - 1;
    if (0 <= loop_ub) {
      memcpy(&Y_data[0], &B_data[0], (loop_ub + 1) * sizeof(real_T));
    }

    for (loop_ub = 0; loop_ub < 3; loop_ub++) {
      if (loop_ub + 1 != ipiv_data[loop_ub]) {
        ip = ipiv_data[loop_ub] - 1;
        temp = Y_data[loop_ub];
        Y_data[loop_ub] = Y_data[ip];
        Y_data[ip] = temp;
        temp_tmp = loop_ub + Y_size[0];
        temp = Y_data[temp_tmp];
        Y_data_tmp = ip + Y_size[0];
        Y_data[temp_tmp] = Y_data[Y_data_tmp];
        Y_data[Y_data_tmp] = temp;
        temp = Y_data[(Y_size[0] << 1) + loop_ub];
        Y_data[loop_ub + (Y_size[0] << 1)] = Y_data[(Y_size[0] << 1) + ip];
        Y_data[ip + (Y_size[0] << 1)] = temp;
        temp = Y_data[Y_size[0] * 3 + loop_ub];
        Y_data[loop_ub + Y_size[0] * 3] = Y_data[Y_size[0] * 3 + ip];
        Y_data[ip + Y_size[0] * 3] = temp;
      }
    }

    for (loop_ub = 0; loop_ub < 4; loop_ub++) {
      ip = loop_ub << 2;
      if (Y_data[ip] != 0.0) {
        for (temp_tmp = 1; temp_tmp < 4; temp_tmp++) {
          Y_data[temp_tmp + ip] -= Y_data[ip] * b_A_data[temp_tmp];
        }
      }

      if (Y_data[1 + ip] != 0.0) {
        for (temp_tmp = 2; temp_tmp < 4; temp_tmp++) {
          Y_data[temp_tmp + ip] -= Y_data[1 + ip] * b_A_data[temp_tmp + 4];
        }
      }

      if (Y_data[2 + ip] != 0.0) {
        for (temp_tmp = 3; temp_tmp < 4; temp_tmp++) {
          Y_data[temp_tmp + ip] -= Y_data[2 + ip] * b_A_data[temp_tmp + 8];
        }
      }
    }

    Experiment_2_moving_cart_xtrsm(b_A_data, Y_data, Y_size);
  } else {
    Experiment_2_moving_car_qrsolve(A_data, A_size, B_data, B_size, b_A_data,
      b_A_size);
    Y_size[0] = b_A_size[0];
    Y_size[1] = b_A_size[1];
    loop_ub = b_A_size[0] * b_A_size[1] - 1;
    if (0 <= loop_ub) {
      memcpy(&Y_data[0], &b_A_data[0], (loop_ub + 1) * sizeof(real_T));
    }
  }
}

/* Function for MATLAB Function: '<S5>/Differentiator' */
static void Exper_PadeApproximantOfDegree_e(const real_T A_data[], const int32_T
  A_size[2], real_T F_data[], int32_T F_size[2])
{
  int32_T b_n;
  real_T U_data[16];
  real_T V_data[16];
  real_T y_data[16];
  int32_T b_m;
  int32_T coffset;
  int32_T boffset;
  int32_T aoffset;
  int32_T j;
  int32_T i;
  int32_T b_i;
  real_T b_y_data[16];
  real_T c_y_data[16];
  real_T d_y_data[16];
  real_T b_b_data[16];
  real_T e_y_data[16];
  int32_T U_size[2];
  int32_T V_size[2];
  int32_T y_size_idx_0;
  int32_T e_y_size_idx_0;
  real_T y_data_0;
  int32_T y_data_tmp;
  int32_T d_y_data_tmp;
  b_n = A_size[0] - 1;
  if (A_size[0] == 1) {
    y_size_idx_0 = 1;
    for (j = 0; j < 1; j++) {
      for (coffset = 0; coffset < 4; coffset++) {
        y_data_0 = A_data[3 + coffset] * A_data[3] + (A_data[2 + coffset] *
          A_data[2] + (A_data[1 + coffset] * A_data[1] + A_data[0] *
                       A_data[coffset]));
        y_data[coffset] = y_data_0;
      }
    }
  } else {
    b_m = A_size[0];
    y_size_idx_0 = A_size[0];
    for (j = 0; j < 4; j++) {
      coffset = j * b_m;
      boffset = j << 2;
      for (i = 1; i <= b_m; i++) {
        y_data[(coffset + i) - 1] = 0.0;
      }

      for (i = 0; i < 4; i++) {
        if (A_data[boffset + i] != 0.0) {
          aoffset = i * b_m;
          for (b_i = 0; b_i < b_m; b_i++) {
            y_data_tmp = coffset + b_i;
            y_data[y_data_tmp] += A_data[boffset + i] * A_data[aoffset + b_i];
          }
        }
      }
    }
  }

  if (y_size_idx_0 == 1) {
    for (j = 0; j < 1; j++) {
      for (coffset = 0; coffset < 4; coffset++) {
        b_y_data[coffset] = 0.0;
        for (boffset = 0; boffset < 4; boffset++) {
          b_y_data[coffset] += y_data[boffset + coffset] * y_data[boffset];
        }
      }
    }

    y_size_idx_0 = 1;
    for (j = 0; j < 1; j++) {
      for (coffset = 0; coffset < 4; coffset++) {
        c_y_data[coffset] = 0.0;
        for (boffset = 0; boffset < 4; boffset++) {
          c_y_data[coffset] += y_data[boffset + coffset] * b_y_data[boffset];
        }
      }
    }
  } else {
    for (j = 0; j < 4; j++) {
      coffset = j * y_size_idx_0;
      boffset = j << 2;
      for (i = 1; i <= y_size_idx_0; i++) {
        b_y_data[(coffset + i) - 1] = 0.0;
      }

      for (i = 0; i < 4; i++) {
        if (y_data[boffset + i] != 0.0) {
          aoffset = i * y_size_idx_0;
          for (b_i = 0; b_i < y_size_idx_0; b_i++) {
            y_data_tmp = coffset + b_i;
            b_y_data[y_data_tmp] += y_data[boffset + i] * y_data[aoffset + b_i];
          }
        }
      }
    }

    for (j = 0; j < 4; j++) {
      coffset = j * y_size_idx_0;
      boffset = j << 2;
      for (i = 1; i <= y_size_idx_0; i++) {
        c_y_data[(coffset + i) - 1] = 0.0;
      }

      for (i = 0; i < 4; i++) {
        if (y_data[boffset + i] != 0.0) {
          aoffset = i * y_size_idx_0;
          for (b_i = 0; b_i < y_size_idx_0; b_i++) {
            y_data_tmp = coffset + b_i;
            c_y_data[y_data_tmp] += y_data[boffset + i] * b_y_data[aoffset + b_i];
          }
        }
      }
    }
  }

  y_data_tmp = (y_size_idx_0 << 2) - 1;
  for (j = 0; j <= y_data_tmp; j++) {
    U_data[j] = (3.352212864E+10 * c_y_data[j] + 1.05594705216E+13 * b_y_data[j])
      + 1.1873537964288E+15 * y_data[j];
  }

  for (i = 0; i <= b_n; i++) {
    e_y_size_idx_0 = y_size_idx_0 * i;
    U_data[i + e_y_size_idx_0] += 3.238237626624E+16;
  }

  for (j = 0; j <= y_data_tmp; j++) {
    b_b_data[j] = (16380.0 * b_y_data[j] + c_y_data[j]) + 4.08408E+7 * y_data[j];
  }

  if (y_size_idx_0 == 1) {
    for (j = 0; j < y_size_idx_0; j++) {
      for (coffset = 0; coffset < 4; coffset++) {
        d_y_data_tmp = y_size_idx_0 * coffset;
        i = j + d_y_data_tmp;
        d_y_data[i] = 0.0;
        for (boffset = 0; boffset < 4; boffset++) {
          d_y_data[i] = c_y_data[y_size_idx_0 * boffset + j] * b_b_data[boffset
            + coffset] + d_y_data[d_y_data_tmp + j];
        }
      }
    }
  } else {
    for (j = 0; j < 4; j++) {
      coffset = j * y_size_idx_0;
      boffset = j << 2;
      for (i = 1; i <= y_size_idx_0; i++) {
        d_y_data[(coffset + i) - 1] = 0.0;
      }

      for (i = 0; i < 4; i++) {
        if (b_b_data[boffset + i] != 0.0) {
          aoffset = i * y_size_idx_0;
          for (b_i = 0; b_i < y_size_idx_0; b_i++) {
            d_y_data_tmp = coffset + b_i;
            d_y_data[d_y_data_tmp] += b_b_data[boffset + i] * c_y_data[aoffset +
              b_i];
          }
        }
      }
    }
  }

  boffset = (y_size_idx_0 << 2) - 1;
  for (j = 0; j <= boffset; j++) {
    d_y_data[j] += U_data[j];
  }

  if (y_size_idx_0 == 1) {
    e_y_size_idx_0 = A_size[0];
    boffset = A_size[0];
    for (j = 0; j < boffset; j++) {
      for (coffset = 0; coffset < 4; coffset++) {
        d_y_data_tmp = e_y_size_idx_0 * coffset;
        i = j + d_y_data_tmp;
        e_y_data[i] = 0.0;
        e_y_data[i] = e_y_data[d_y_data_tmp + j] + A_data[j] * d_y_data[coffset];
        e_y_data[i] = e_y_data[e_y_size_idx_0 * coffset + j] + A_data[j +
          A_size[0]] * d_y_data[1 + coffset];
        e_y_data[i] = A_data[(A_size[0] << 1) + j] * d_y_data[2 + coffset] +
          e_y_data[e_y_size_idx_0 * coffset + j];
        e_y_data[i] = A_data[A_size[0] * 3 + j] * d_y_data[3 + coffset] +
          e_y_data[e_y_size_idx_0 * coffset + j];
      }
    }
  } else {
    b_m = A_size[0];
    e_y_size_idx_0 = A_size[0];
    for (j = 0; j < 4; j++) {
      coffset = j * b_m;
      boffset = j << 2;
      for (i = 1; i <= b_m; i++) {
        e_y_data[(coffset + i) - 1] = 0.0;
      }

      for (i = 0; i < 4; i++) {
        if (d_y_data[boffset + i] != 0.0) {
          aoffset = i * b_m;
          for (b_i = 0; b_i < b_m; b_i++) {
            d_y_data_tmp = coffset + b_i;
            e_y_data[d_y_data_tmp] += d_y_data[boffset + i] * A_data[aoffset +
              b_i];
          }
        }
      }
    }
  }

  U_size[0] = e_y_size_idx_0;
  U_size[1] = 4;
  boffset = (e_y_size_idx_0 << 2) - 1;
  if (0 <= boffset) {
    memcpy(&U_data[0], &e_y_data[0], (boffset + 1) * sizeof(real_T));
  }

  for (j = 0; j <= y_data_tmp; j++) {
    b_b_data[j] = (182.0 * c_y_data[j] + 960960.0 * b_y_data[j]) + 1.32324192E+9
      * y_data[j];
  }

  if (y_size_idx_0 == 1) {
    for (j = 0; j < y_size_idx_0; j++) {
      for (coffset = 0; coffset < 4; coffset++) {
        d_y_data_tmp = y_size_idx_0 * coffset;
        i = j + d_y_data_tmp;
        d_y_data[i] = 0.0;
        for (boffset = 0; boffset < 4; boffset++) {
          d_y_data[i] = c_y_data[y_size_idx_0 * boffset + j] * b_b_data[boffset
            + coffset] + d_y_data[d_y_data_tmp + j];
        }
      }
    }
  } else {
    for (j = 0; j < 4; j++) {
      coffset = j * y_size_idx_0;
      boffset = j << 2;
      for (i = 1; i <= y_size_idx_0; i++) {
        d_y_data[(coffset + i) - 1] = 0.0;
      }

      for (i = 0; i < 4; i++) {
        if (b_b_data[boffset + i] != 0.0) {
          aoffset = i * y_size_idx_0;
          for (b_i = 0; b_i < y_size_idx_0; b_i++) {
            d_y_data_tmp = coffset + b_i;
            d_y_data[d_y_data_tmp] += b_b_data[boffset + i] * c_y_data[aoffset +
              b_i];
          }
        }
      }
    }
  }

  V_size[0] = y_size_idx_0;
  V_size[1] = 4;
  boffset = (y_size_idx_0 << 2) - 1;
  for (j = 0; j <= boffset; j++) {
    V_data[j] = ((6.704425728E+11 * c_y_data[j] + d_y_data[j]) +
                 1.29060195264E+14 * b_y_data[j]) + 7.7717703038976E+15 *
      y_data[j];
  }

  for (i = 0; i <= b_n; i++) {
    y_data_tmp = y_size_idx_0 * i;
    V_data[i + y_data_tmp] += 6.476475253248E+16;
  }

  y_size_idx_0 = e_y_size_idx_0 << 2;
  for (b_m = 0; b_m < y_size_idx_0; b_m++) {
    V_data[b_m] -= U_data[b_m];
    U_data[b_m] *= 2.0;
  }

  Experiment_2_moving_ca_mldivide(V_data, V_size, U_data, U_size, F_data, F_size);
  for (i = 0; i <= b_n; i++) {
    y_size_idx_0 = F_size[0] * i;
    F_data[i + y_size_idx_0]++;
  }
}

/* Function for MATLAB Function: '<S5>/Differentiator' */
static void Experim_PadeApproximantOfDegree(const real_T A_data[], const int32_T
  A_size[2], uint8_T b_m, real_T F_data[], int32_T F_size[2])
{
  int32_T b_n;
  real_T U_data[16];
  real_T V_data[16];
  real_T b_d;
  real_T y_data[16];
  int32_T c_m;
  int32_T coffset;
  int32_T boffset;
  int32_T aoffset;
  int32_T j;
  int32_T i;
  int32_T b_i;
  real_T b_y_data[16];
  real_T c_y_data[16];
  real_T d_y_data[16];
  real_T b_b_data[16];
  int32_T U_size[2];
  int32_T V_size[2];
  int32_T y_size_idx_0;
  int32_T b_y_size_idx_0;
  int32_T y_data_tmp;
  int32_T b_y_data_tmp;
  int32_T c_y_data_tmp;
  b_n = A_size[0] - 1;
  if (A_size[0] == 1) {
    y_size_idx_0 = 1;
    for (c_m = 0; c_m < 1; c_m++) {
      for (j = 0; j < 4; j++) {
        b_d = A_data[3 + j] * A_data[3] + (A_data[2 + j] * A_data[2] + (A_data[1
          + j] * A_data[1] + A_data[0] * A_data[j]));
        y_data[j] = b_d;
      }
    }
  } else {
    c_m = A_size[0];
    y_size_idx_0 = A_size[0];
    for (j = 0; j < 4; j++) {
      coffset = j * c_m;
      boffset = j << 2;
      for (i = 1; i <= c_m; i++) {
        y_data[(coffset + i) - 1] = 0.0;
      }

      for (i = 0; i < 4; i++) {
        if (A_data[boffset + i] != 0.0) {
          aoffset = i * c_m;
          for (b_i = 0; b_i < c_m; b_i++) {
            y_data_tmp = coffset + b_i;
            y_data[y_data_tmp] += A_data[boffset + i] * A_data[aoffset + b_i];
          }
        }
      }
    }
  }

  if (b_m == 3) {
    y_data_tmp = (y_size_idx_0 << 2) - 1;
    if (0 <= y_data_tmp) {
      memcpy(&U_data[0], &y_data[0], (y_data_tmp + 1) * sizeof(real_T));
    }

    for (i = 0; i <= b_n; i++) {
      c_m = y_size_idx_0 * i;
      U_data[i + c_m] += 60.0;
    }

    if (y_size_idx_0 == 1) {
      b_y_size_idx_0 = A_size[0];
      coffset = A_size[0];
      for (c_m = 0; c_m < coffset; c_m++) {
        for (j = 0; j < 4; j++) {
          b_y_data_tmp = b_y_size_idx_0 * j;
          boffset = c_m + b_y_data_tmp;
          b_y_data[boffset] = 0.0;
          b_y_data[boffset] = b_y_data[b_y_data_tmp + c_m] + A_data[c_m] *
            U_data[j];
          b_y_data[boffset] = b_y_data[b_y_size_idx_0 * j + c_m] + A_data[c_m +
            A_size[0]] * U_data[1 + j];
          b_y_data[boffset] = A_data[(A_size[0] << 1) + c_m] * U_data[2 + j] +
            b_y_data[b_y_size_idx_0 * j + c_m];
          b_y_data[boffset] = A_data[A_size[0] * 3 + c_m] * U_data[3 + j] +
            b_y_data[b_y_size_idx_0 * j + c_m];
        }
      }
    } else {
      c_m = A_size[0];
      b_y_size_idx_0 = A_size[0];
      for (j = 0; j < 4; j++) {
        coffset = j * c_m;
        boffset = j << 2;
        for (i = 1; i <= c_m; i++) {
          b_y_data[(coffset + i) - 1] = 0.0;
        }

        for (i = 0; i < 4; i++) {
          if (U_data[boffset + i] != 0.0) {
            aoffset = i * c_m;
            for (b_i = 0; b_i < c_m; b_i++) {
              b_y_data_tmp = coffset + b_i;
              b_y_data[b_y_data_tmp] += U_data[boffset + i] * A_data[aoffset +
                b_i];
            }
          }
        }
      }
    }

    U_size[0] = b_y_size_idx_0;
    U_size[1] = 4;
    coffset = (b_y_size_idx_0 << 2) - 1;
    if (0 <= coffset) {
      memcpy(&U_data[0], &b_y_data[0], (coffset + 1) * sizeof(real_T));
    }

    V_size[0] = y_size_idx_0;
    V_size[1] = 4;
    for (c_m = 0; c_m <= y_data_tmp; c_m++) {
      V_data[c_m] = 12.0 * y_data[c_m];
    }

    b_d = 120.0;
  } else {
    if (y_size_idx_0 == 1) {
      b_y_size_idx_0 = 1;
      for (c_m = 0; c_m < 1; c_m++) {
        for (j = 0; j < 4; j++) {
          b_y_data[j] = 0.0;
          for (coffset = 0; coffset < 4; coffset++) {
            b_y_data[j] += y_data[coffset + j] * y_data[coffset];
          }
        }
      }
    } else {
      b_y_size_idx_0 = y_size_idx_0;
      for (j = 0; j < 4; j++) {
        coffset = j * y_size_idx_0;
        boffset = j << 2;
        for (i = 1; i <= y_size_idx_0; i++) {
          b_y_data[(coffset + i) - 1] = 0.0;
        }

        for (i = 0; i < 4; i++) {
          if (y_data[boffset + i] != 0.0) {
            aoffset = i * y_size_idx_0;
            for (b_i = 0; b_i < y_size_idx_0; b_i++) {
              b_y_data_tmp = coffset + b_i;
              b_y_data[b_y_data_tmp] += y_data[boffset + i] * y_data[aoffset +
                b_i];
            }
          }
        }
      }
    }

    if (b_m == 5) {
      y_data_tmp = (b_y_size_idx_0 << 2) - 1;
      for (c_m = 0; c_m <= y_data_tmp; c_m++) {
        U_data[c_m] = 420.0 * y_data[c_m] + b_y_data[c_m];
      }

      for (i = 0; i <= b_n; i++) {
        c_m = b_y_size_idx_0 * i;
        U_data[i + c_m] += 15120.0;
      }

      if (b_y_size_idx_0 == 1) {
        b_y_data_tmp = A_size[0];
        coffset = A_size[0];
        for (c_m = 0; c_m < coffset; c_m++) {
          for (j = 0; j < 4; j++) {
            c_y_data_tmp = b_y_data_tmp * j;
            boffset = c_m + c_y_data_tmp;
            c_y_data[boffset] = 0.0;
            c_y_data[boffset] = c_y_data[c_y_data_tmp + c_m] + A_data[c_m] *
              U_data[j];
            c_y_data[boffset] = c_y_data[b_y_data_tmp * j + c_m] + A_data[c_m +
              A_size[0]] * U_data[1 + j];
            c_y_data[boffset] = A_data[(A_size[0] << 1) + c_m] * U_data[2 + j] +
              c_y_data[b_y_data_tmp * j + c_m];
            c_y_data[boffset] = A_data[A_size[0] * 3 + c_m] * U_data[3 + j] +
              c_y_data[b_y_data_tmp * j + c_m];
          }
        }
      } else {
        c_m = A_size[0];
        b_y_data_tmp = A_size[0];
        for (j = 0; j < 4; j++) {
          coffset = j * c_m;
          boffset = j << 2;
          for (i = 1; i <= c_m; i++) {
            c_y_data[(coffset + i) - 1] = 0.0;
          }

          for (i = 0; i < 4; i++) {
            if (U_data[boffset + i] != 0.0) {
              aoffset = i * c_m;
              for (b_i = 0; b_i < c_m; b_i++) {
                c_y_data_tmp = coffset + b_i;
                c_y_data[c_y_data_tmp] += U_data[boffset + i] * A_data[aoffset +
                  b_i];
              }
            }
          }
        }
      }

      U_size[0] = b_y_data_tmp;
      U_size[1] = 4;
      coffset = (b_y_data_tmp << 2) - 1;
      if (0 <= coffset) {
        memcpy(&U_data[0], &c_y_data[0], (coffset + 1) * sizeof(real_T));
      }

      V_size[0] = b_y_size_idx_0;
      V_size[1] = 4;
      for (c_m = 0; c_m <= y_data_tmp; c_m++) {
        V_data[c_m] = 30.0 * b_y_data[c_m] + 3360.0 * y_data[c_m];
      }

      b_d = 30240.0;
    } else {
      if (y_size_idx_0 == 1) {
        b_y_data_tmp = b_y_size_idx_0;
        for (c_m = 0; c_m < b_y_size_idx_0; c_m++) {
          for (j = 0; j < 4; j++) {
            c_y_data_tmp = b_y_size_idx_0 * j;
            boffset = c_m + c_y_data_tmp;
            c_y_data[boffset] = 0.0;
            for (coffset = 0; coffset < 4; coffset++) {
              c_y_data[boffset] = b_y_data[b_y_size_idx_0 * coffset + c_m] *
                y_data[coffset + j] + c_y_data[c_y_data_tmp + c_m];
            }
          }
        }
      } else {
        b_y_data_tmp = b_y_size_idx_0;
        for (j = 0; j < 4; j++) {
          coffset = j * b_y_size_idx_0;
          boffset = j << 2;
          for (i = 1; i <= b_y_size_idx_0; i++) {
            c_y_data[(coffset + i) - 1] = 0.0;
          }

          for (i = 0; i < 4; i++) {
            if (y_data[boffset + i] != 0.0) {
              aoffset = i * b_y_size_idx_0;
              for (b_i = 0; b_i < b_y_size_idx_0; b_i++) {
                c_y_data_tmp = coffset + b_i;
                c_y_data[c_y_data_tmp] += y_data[boffset + i] * b_y_data[aoffset
                  + b_i];
              }
            }
          }
        }
      }

      if (b_m == 7) {
        coffset = (b_y_data_tmp << 2) - 1;
        for (c_m = 0; c_m <= coffset; c_m++) {
          U_data[c_m] = (1512.0 * b_y_data[c_m] + c_y_data[c_m]) + 277200.0 *
            y_data[c_m];
        }

        for (i = 0; i <= b_n; i++) {
          c_m = b_y_data_tmp * i;
          U_data[i + c_m] += 8.64864E+6;
        }

        if (b_y_data_tmp == 1) {
          y_size_idx_0 = A_size[0];
          coffset = A_size[0];
          for (c_m = 0; c_m < coffset; c_m++) {
            for (j = 0; j < 4; j++) {
              b_y_size_idx_0 = y_size_idx_0 * j;
              boffset = c_m + b_y_size_idx_0;
              d_y_data[boffset] = 0.0;
              d_y_data[boffset] = d_y_data[b_y_size_idx_0 + c_m] + A_data[c_m] *
                U_data[j];
              d_y_data[boffset] = d_y_data[y_size_idx_0 * j + c_m] + A_data[c_m
                + A_size[0]] * U_data[1 + j];
              d_y_data[boffset] = A_data[(A_size[0] << 1) + c_m] * U_data[2 + j]
                + d_y_data[y_size_idx_0 * j + c_m];
              d_y_data[boffset] = A_data[A_size[0] * 3 + c_m] * U_data[3 + j] +
                d_y_data[y_size_idx_0 * j + c_m];
            }
          }
        } else {
          c_m = A_size[0];
          y_size_idx_0 = A_size[0];
          for (j = 0; j < 4; j++) {
            coffset = j * c_m;
            boffset = j << 2;
            for (i = 1; i <= c_m; i++) {
              d_y_data[(coffset + i) - 1] = 0.0;
            }

            for (i = 0; i < 4; i++) {
              if (U_data[boffset + i] != 0.0) {
                aoffset = i * c_m;
                for (b_i = 0; b_i < c_m; b_i++) {
                  b_y_size_idx_0 = coffset + b_i;
                  d_y_data[b_y_size_idx_0] += U_data[boffset + i] *
                    A_data[aoffset + b_i];
                }
              }
            }
          }
        }

        U_size[0] = y_size_idx_0;
        U_size[1] = 4;
        coffset = (y_size_idx_0 << 2) - 1;
        if (0 <= coffset) {
          memcpy(&U_data[0], &d_y_data[0], (coffset + 1) * sizeof(real_T));
        }

        V_size[0] = b_y_data_tmp;
        V_size[1] = 4;
        coffset = (b_y_data_tmp << 2) - 1;
        for (c_m = 0; c_m <= coffset; c_m++) {
          V_data[c_m] = (56.0 * c_y_data[c_m] + 25200.0 * b_y_data[c_m]) +
            1.99584E+6 * y_data[c_m];
        }

        b_d = 1.729728E+7;
      } else if (b_m == 9) {
        if (y_size_idx_0 == 1) {
          y_size_idx_0 = b_y_data_tmp;
          for (c_m = 0; c_m < b_y_data_tmp; c_m++) {
            for (j = 0; j < 4; j++) {
              b_y_size_idx_0 = b_y_data_tmp * j;
              boffset = c_m + b_y_size_idx_0;
              d_y_data[boffset] = 0.0;
              for (coffset = 0; coffset < 4; coffset++) {
                d_y_data[boffset] = c_y_data[b_y_data_tmp * coffset + c_m] *
                  y_data[coffset + j] + d_y_data[b_y_size_idx_0 + c_m];
              }
            }
          }
        } else {
          y_size_idx_0 = b_y_data_tmp;
          for (j = 0; j < 4; j++) {
            coffset = j * b_y_data_tmp;
            boffset = j << 2;
            for (i = 1; i <= b_y_data_tmp; i++) {
              d_y_data[(coffset + i) - 1] = 0.0;
            }

            for (i = 0; i < 4; i++) {
              if (y_data[boffset + i] != 0.0) {
                aoffset = i * b_y_data_tmp;
                for (b_i = 0; b_i < b_y_data_tmp; b_i++) {
                  b_y_size_idx_0 = coffset + b_i;
                  d_y_data[b_y_size_idx_0] += y_data[boffset + i] *
                    c_y_data[aoffset + b_i];
                }
              }
            }
          }
        }

        y_data_tmp = (y_size_idx_0 << 2) - 1;
        for (c_m = 0; c_m <= y_data_tmp; c_m++) {
          U_data[c_m] = ((3960.0 * c_y_data[c_m] + d_y_data[c_m]) + 2.16216E+6 *
                         b_y_data[c_m]) + 3.027024E+8 * y_data[c_m];
        }

        for (i = 0; i <= b_n; i++) {
          c_m = y_size_idx_0 * i;
          U_data[i + c_m] += 8.8216128E+9;
        }

        if (y_size_idx_0 == 1) {
          b_y_size_idx_0 = A_size[0];
          coffset = A_size[0];
          for (c_m = 0; c_m < coffset; c_m++) {
            for (j = 0; j < 4; j++) {
              c_y_data_tmp = b_y_size_idx_0 * j;
              boffset = c_m + c_y_data_tmp;
              b_b_data[boffset] = 0.0;
              b_b_data[boffset] = b_b_data[c_y_data_tmp + c_m] + A_data[c_m] *
                U_data[j];
              b_b_data[boffset] = b_b_data[b_y_size_idx_0 * j + c_m] +
                A_data[c_m + A_size[0]] * U_data[1 + j];
              b_b_data[boffset] = A_data[(A_size[0] << 1) + c_m] * U_data[2 + j]
                + b_b_data[b_y_size_idx_0 * j + c_m];
              b_b_data[boffset] = A_data[A_size[0] * 3 + c_m] * U_data[3 + j] +
                b_b_data[b_y_size_idx_0 * j + c_m];
            }
          }
        } else {
          c_m = A_size[0];
          b_y_size_idx_0 = A_size[0];
          for (j = 0; j < 4; j++) {
            coffset = j * c_m;
            boffset = j << 2;
            for (i = 1; i <= c_m; i++) {
              b_b_data[(coffset + i) - 1] = 0.0;
            }

            for (i = 0; i < 4; i++) {
              if (U_data[boffset + i] != 0.0) {
                aoffset = i * c_m;
                for (b_i = 0; b_i < c_m; b_i++) {
                  c_y_data_tmp = coffset + b_i;
                  b_b_data[c_y_data_tmp] += U_data[boffset + i] * A_data[aoffset
                    + b_i];
                }
              }
            }
          }
        }

        U_size[0] = b_y_size_idx_0;
        U_size[1] = 4;
        coffset = (b_y_size_idx_0 << 2) - 1;
        if (0 <= coffset) {
          memcpy(&U_data[0], &b_b_data[0], (coffset + 1) * sizeof(real_T));
        }

        V_size[0] = y_size_idx_0;
        V_size[1] = 4;
        for (c_m = 0; c_m <= y_data_tmp; c_m++) {
          V_data[c_m] = ((90.0 * d_y_data[c_m] + 110880.0 * c_y_data[c_m]) +
                         3.027024E+7 * b_y_data[c_m]) + 2.0756736E+9 *
            y_data[c_m];
        }

        b_d = 1.76432256E+10;
      } else {
        y_data_tmp = (b_y_data_tmp << 2) - 1;
        for (c_m = 0; c_m <= y_data_tmp; c_m++) {
          U_data[c_m] = (3.352212864E+10 * c_y_data[c_m] + 1.05594705216E+13 *
                         b_y_data[c_m]) + 1.1873537964288E+15 * y_data[c_m];
        }

        for (i = 0; i <= b_n; i++) {
          c_m = b_y_data_tmp * i;
          U_data[i + c_m] += 3.238237626624E+16;
        }

        for (c_m = 0; c_m <= y_data_tmp; c_m++) {
          b_b_data[c_m] = (16380.0 * b_y_data[c_m] + c_y_data[c_m]) + 4.08408E+7
            * y_data[c_m];
        }

        if (b_y_data_tmp == 1) {
          y_size_idx_0 = b_y_data_tmp;
          for (c_m = 0; c_m < b_y_data_tmp; c_m++) {
            for (j = 0; j < 4; j++) {
              b_y_size_idx_0 = b_y_data_tmp * j;
              boffset = c_m + b_y_size_idx_0;
              d_y_data[boffset] = 0.0;
              for (coffset = 0; coffset < 4; coffset++) {
                d_y_data[boffset] = c_y_data[b_y_data_tmp * coffset + c_m] *
                  b_b_data[coffset + j] + d_y_data[b_y_size_idx_0 + c_m];
              }
            }
          }
        } else {
          y_size_idx_0 = b_y_data_tmp;
          for (j = 0; j < 4; j++) {
            coffset = j * b_y_data_tmp;
            boffset = j << 2;
            for (i = 1; i <= b_y_data_tmp; i++) {
              d_y_data[(coffset + i) - 1] = 0.0;
            }

            for (i = 0; i < 4; i++) {
              if (b_b_data[boffset + i] != 0.0) {
                aoffset = i * b_y_data_tmp;
                for (b_i = 0; b_i < b_y_data_tmp; b_i++) {
                  b_y_size_idx_0 = coffset + b_i;
                  d_y_data[b_y_size_idx_0] += b_b_data[boffset + i] *
                    c_y_data[aoffset + b_i];
                }
              }
            }
          }
        }

        coffset = (y_size_idx_0 << 2) - 1;
        for (c_m = 0; c_m <= coffset; c_m++) {
          d_y_data[c_m] += U_data[c_m];
        }

        if (y_size_idx_0 == 1) {
          b_y_size_idx_0 = A_size[0];
          coffset = A_size[0];
          for (c_m = 0; c_m < coffset; c_m++) {
            for (j = 0; j < 4; j++) {
              c_y_data_tmp = b_y_size_idx_0 * j;
              boffset = c_m + c_y_data_tmp;
              b_b_data[boffset] = 0.0;
              b_b_data[boffset] = b_b_data[c_y_data_tmp + c_m] + A_data[c_m] *
                d_y_data[j];
              b_b_data[boffset] = b_b_data[b_y_size_idx_0 * j + c_m] +
                A_data[c_m + A_size[0]] * d_y_data[1 + j];
              b_b_data[boffset] = A_data[(A_size[0] << 1) + c_m] * d_y_data[2 +
                j] + b_b_data[b_y_size_idx_0 * j + c_m];
              b_b_data[boffset] = A_data[A_size[0] * 3 + c_m] * d_y_data[3 + j]
                + b_b_data[b_y_size_idx_0 * j + c_m];
            }
          }
        } else {
          c_m = A_size[0];
          b_y_size_idx_0 = A_size[0];
          for (j = 0; j < 4; j++) {
            coffset = j * c_m;
            boffset = j << 2;
            for (i = 1; i <= c_m; i++) {
              b_b_data[(coffset + i) - 1] = 0.0;
            }

            for (i = 0; i < 4; i++) {
              if (d_y_data[boffset + i] != 0.0) {
                aoffset = i * c_m;
                for (b_i = 0; b_i < c_m; b_i++) {
                  c_y_data_tmp = coffset + b_i;
                  b_b_data[c_y_data_tmp] += d_y_data[boffset + i] *
                    A_data[aoffset + b_i];
                }
              }
            }
          }
        }

        U_size[0] = b_y_size_idx_0;
        U_size[1] = 4;
        coffset = (b_y_size_idx_0 << 2) - 1;
        if (0 <= coffset) {
          memcpy(&U_data[0], &b_b_data[0], (coffset + 1) * sizeof(real_T));
        }

        for (c_m = 0; c_m <= y_data_tmp; c_m++) {
          b_b_data[c_m] = (182.0 * c_y_data[c_m] + 960960.0 * b_y_data[c_m]) +
            1.32324192E+9 * y_data[c_m];
        }

        if (b_y_data_tmp == 1) {
          for (c_m = 0; c_m < b_y_data_tmp; c_m++) {
            for (j = 0; j < 4; j++) {
              b_y_size_idx_0 = b_y_data_tmp * j;
              boffset = c_m + b_y_size_idx_0;
              d_y_data[boffset] = 0.0;
              for (coffset = 0; coffset < 4; coffset++) {
                d_y_data[boffset] = c_y_data[b_y_data_tmp * coffset + c_m] *
                  b_b_data[coffset + j] + d_y_data[b_y_size_idx_0 + c_m];
              }
            }
          }
        } else {
          for (j = 0; j < 4; j++) {
            coffset = j * b_y_data_tmp;
            boffset = j << 2;
            for (i = 1; i <= b_y_data_tmp; i++) {
              d_y_data[(coffset + i) - 1] = 0.0;
            }

            for (i = 0; i < 4; i++) {
              if (b_b_data[boffset + i] != 0.0) {
                aoffset = i * b_y_data_tmp;
                for (b_i = 0; b_i < b_y_data_tmp; b_i++) {
                  b_y_size_idx_0 = coffset + b_i;
                  d_y_data[b_y_size_idx_0] += b_b_data[boffset + i] *
                    c_y_data[aoffset + b_i];
                }
              }
            }
          }
        }

        V_size[0] = y_size_idx_0;
        V_size[1] = 4;
        coffset = (y_size_idx_0 << 2) - 1;
        for (c_m = 0; c_m <= coffset; c_m++) {
          V_data[c_m] = ((6.704425728E+11 * c_y_data[c_m] + d_y_data[c_m]) +
                         1.29060195264E+14 * b_y_data[c_m]) +
            7.7717703038976E+15 * y_data[c_m];
        }

        b_d = 6.476475253248E+16;
      }
    }
  }

  for (i = 0; i <= b_n; i++) {
    y_data_tmp = V_size[0] * i;
    V_data[i + y_data_tmp] += b_d;
  }

  y_data_tmp = U_size[0] << 2;
  for (c_m = 0; c_m < y_data_tmp; c_m++) {
    V_data[c_m] -= U_data[c_m];
    U_data[c_m] *= 2.0;
  }

  Experiment_2_moving_ca_mldivide(V_data, V_size, U_data, U_size, F_data, F_size);
  for (i = 0; i <= b_n; i++) {
    y_data_tmp = F_size[0] * i;
    F_data[i + y_data_tmp]++;
  }
}

/* Function for MATLAB Function: '<S5>/Differentiator' */
static void Experiment_2_moving_ca_mrdivide(const real_T A_data[], const int32_T
  A_size[2], const real_T B_data[], const int32_T B_size[2], real_T y_data[],
  int32_T y_size[2])
{
  int32_T b_n;
  int32_T jp;
  int32_T nb;
  real_T temp;
  real_T b_A_data[16];
  int32_T ipiv_data[4];
  int32_T b_jAcol;
  int32_T b_jBcol;
  int32_T b_kBcol;
  int32_T b_k;
  int32_T c_i;
  real_T B_data_0[9];
  real_T A_data_0[9];
  int32_T loop_ub;
  int32_T b_A_size[2];
  int32_T B_size_0[2];
  int32_T A_size_0[2];
  int32_T loop_ub_tmp;
  int32_T X_data_tmp;
  if (B_size[0] == B_size[1]) {
    b_n = B_size[1];
    b_A_size[0] = B_size[0];
    b_A_size[1] = B_size[1];
    nb = B_size[0] * B_size[1] - 1;
    if (0 <= nb) {
      memcpy(&b_A_data[0], &B_data[0], (nb + 1) * sizeof(real_T));
    }

    Experiment_2_moving_car_xzgetrf(B_size[1], B_size[1], b_A_data, b_A_size,
      B_size[1], ipiv_data, B_size_0, &nb);
    nb = A_size[0];
    loop_ub = A_size[0];
    loop_ub_tmp = A_size[0] * A_size[1] - 1;
    if (0 <= loop_ub_tmp) {
      memcpy(&B_data_0[0], &A_data[0], (loop_ub_tmp + 1) * sizeof(real_T));
    }

    for (jp = 0; jp < b_n; jp++) {
      b_jBcol = nb * jp;
      b_jAcol = b_n * jp;
      for (b_k = 0; b_k < jp; b_k++) {
        b_kBcol = nb * b_k;
        if (b_A_data[b_k + b_jAcol] != 0.0) {
          for (c_i = 0; c_i < nb; c_i++) {
            X_data_tmp = c_i + b_jBcol;
            B_data_0[X_data_tmp] -= b_A_data[b_k + b_jAcol] * B_data_0[c_i +
              b_kBcol];
          }
        }
      }

      temp = 1.0 / b_A_data[jp + b_jAcol];
      for (b_jAcol = 0; b_jAcol < nb; b_jAcol++) {
        X_data_tmp = b_jAcol + b_jBcol;
        B_data_0[X_data_tmp] *= temp;
      }
    }

    for (jp = B_size[1]; jp > 0; jp--) {
      b_jBcol = (jp - 1) * nb;
      b_jAcol = (jp - 1) * b_n;
      for (b_k = jp; b_k < b_n; b_k++) {
        b_kBcol = nb * b_k;
        if (b_A_data[b_k + b_jAcol] != 0.0) {
          for (c_i = 0; c_i < nb; c_i++) {
            X_data_tmp = c_i + b_jBcol;
            B_data_0[X_data_tmp] -= b_A_data[b_k + b_jAcol] * B_data_0[c_i +
              b_kBcol];
          }
        }
      }
    }

    for (b_n = B_size[1] - 2; b_n + 1 > 0; b_n--) {
      if (b_n + 1 != ipiv_data[b_n]) {
        jp = ipiv_data[b_n] - 1;
        for (b_jAcol = 0; b_jAcol < nb; b_jAcol++) {
          b_jBcol = loop_ub * b_n;
          temp = B_data_0[b_jBcol + b_jAcol];
          X_data_tmp = loop_ub * jp;
          B_data_0[b_jAcol + b_jBcol] = B_data_0[X_data_tmp + b_jAcol];
          B_data_0[b_jAcol + X_data_tmp] = temp;
        }
      }
    }

    y_size[0] = A_size[0];
    y_size[1] = A_size[1];
    if (0 <= loop_ub_tmp) {
      memcpy(&y_data[0], &B_data_0[0], (loop_ub_tmp + 1) * sizeof(real_T));
    }
  } else {
    B_size_0[0] = B_size[1];
    B_size_0[1] = B_size[0];
    nb = B_size[0];
    for (b_n = 0; b_n < nb; b_n++) {
      loop_ub = B_size[1];
      for (loop_ub_tmp = 0; loop_ub_tmp < loop_ub; loop_ub_tmp++) {
        B_data_0[loop_ub_tmp + B_size_0[0] * b_n] = B_data[B_size[0] *
          loop_ub_tmp + b_n];
      }
    }

    A_size_0[0] = A_size[1];
    A_size_0[1] = A_size[0];
    nb = A_size[0];
    for (b_n = 0; b_n < nb; b_n++) {
      loop_ub = A_size[1];
      for (loop_ub_tmp = 0; loop_ub_tmp < loop_ub; loop_ub_tmp++) {
        A_data_0[loop_ub_tmp + A_size_0[0] * b_n] = A_data[A_size[0] *
          loop_ub_tmp + b_n];
      }
    }

    Experiment_2_moving_car_qrsolve(B_data_0, B_size_0, A_data_0, A_size_0,
      b_A_data, b_A_size);
    y_size[0] = b_A_size[1];
    y_size[1] = b_A_size[0];
    nb = b_A_size[0];
    for (b_n = 0; b_n < nb; b_n++) {
      loop_ub = b_A_size[1];
      for (loop_ub_tmp = 0; loop_ub_tmp < loop_ub; loop_ub_tmp++) {
        y_data[loop_ub_tmp + y_size[0] * b_n] = b_A_data[b_A_size[0] *
          loop_ub_tmp + b_n];
      }
    }
  }
}

/* Function for MATLAB Function: '<S5>/Differentiator' */
static real_T Experimen_disc_eigenvalues_URED(real_T x0, real_T b_mu, real_T tau)
{
  real_T z;
  real_T c;
  c = rt_powd_snf(fabs(x0), 0.33333333333333331);
  z = c / ((tau * 15.0 * b_mu * fabs(x0) + c) + tau * 15.0);
  return z;
}

/* Function for MATLAB Function: '<S5>/Differentiator' */
static void Experiment_ackerman_precomputed(real_T b_Ts, real_T z, real_T
  lambda_data[], int32_T *lambda_size)
{
  *lambda_size = 3;
  lambda_data[0] = (3.0 - z) - 2.0 * z;
  lambda_data[1] = (((3.0 * z + z) + z * z) - 5.0) * (z - 1.0) / (2.0 * b_Ts);
  lambda_data[2] = -((z - 1.0) * (z - 1.0) * (z - 1.0)) / (b_Ts * b_Ts);
}

/* Function for MATLAB Function: '<S5>/Differentiator' */
static void Experiment_2_moving_cart_step_p(real_T x0, const real_T z_data[],
  real_T u, const real_T Phi_data[], const real_T bD_data[], const real_T
  lambda_data[], real_T z_new_data[], int32_T *z_new_size)
{
  real_T y_data[3];
  int32_T aoffset;
  int32_T i;
  int32_T b_i;
  for (i = 1; i < 4; i++) {
    y_data[i - 1] = 0.0;
  }

  for (i = 0; i < 3; i++) {
    if (z_data[i] != 0.0) {
      aoffset = i * 3;
      for (b_i = 0; b_i < 3; b_i++) {
        y_data[b_i] += Phi_data[aoffset + b_i] * z_data[i];
      }
    }
  }

  *z_new_size = 3;
  for (i = 0; i < 3; i++) {
    z_new_data[i] = (bD_data[i] * u + y_data[i]) + lambda_data[i] * x0;
  }
}

/* Model output function */
void Experiment_2_moving_cart_output(void)
{
  /* local block i/o variables */
  real_T rtb_HILReadEncoder_o1;
  real_T rtb_HILReadEncoder_o2;
  real_T P_data[9];
  real_T So_inv_data[6];
  b_cell_wrap_1_Experiment_2_mo_T reshapes[2];
  real_T y_data[9];
  real_T normA;
  int32_T j;
  real_T c_y_data[9];
  int32_T b_m;
  int32_T coffset;
  int32_T boffset;
  int32_T aoffset;
  int32_T b_j;
  int32_T c_i;
  real_T d_y_data[4];
  int8_T e_y_data[2];
  int8_T f_y_data[4];
  real_T g_y_data[4];
  real_T b_n_data[2];
  real_T z1_data[2];
  int8_T b_b_data[2];
  static const real_T c[170] = { 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0,
    40320.0, 362880.0, 3.6288E+6, 3.99168E+7, 4.790016E+8, 6.2270208E+9,
    8.71782912E+10, 1.307674368E+12, 2.0922789888E+13, 3.55687428096E+14,
    6.402373705728E+15, 1.21645100408832E+17, 2.43290200817664E+18,
    5.109094217170944E+19, 1.1240007277776077E+21, 2.5852016738884978E+22,
    6.2044840173323941E+23, 1.5511210043330986E+25, 4.0329146112660565E+26,
    1.0888869450418352E+28, 3.0488834461171384E+29, 8.8417619937397008E+30,
    2.6525285981219103E+32, 8.2228386541779224E+33, 2.6313083693369352E+35,
    8.6833176188118859E+36, 2.9523279903960412E+38, 1.0333147966386144E+40,
    3.7199332678990118E+41, 1.3763753091226343E+43, 5.23022617466601E+44,
    2.0397882081197442E+46, 8.1591528324789768E+47, 3.3452526613163803E+49,
    1.4050061177528798E+51, 6.0415263063373834E+52, 2.6582715747884485E+54,
    1.1962222086548019E+56, 5.5026221598120885E+57, 2.5862324151116818E+59,
    1.2413915592536073E+61, 6.0828186403426752E+62, 3.0414093201713376E+64,
    1.5511187532873822E+66, 8.0658175170943877E+67, 4.2748832840600255E+69,
    2.3084369733924138E+71, 1.2696403353658276E+73, 7.1099858780486348E+74,
    4.0526919504877221E+76, 2.3505613312828789E+78, 1.3868311854568986E+80,
    8.3209871127413916E+81, 5.0758021387722484E+83, 3.1469973260387939E+85,
    1.98260831540444E+87, 1.2688693218588417E+89, 8.2476505920824715E+90,
    5.4434493907744307E+92, 3.6471110918188683E+94, 2.4800355424368305E+96,
    1.711224524281413E+98, 1.197857166996989E+100, 8.5047858856786218E+101,
    6.1234458376886077E+103, 4.4701154615126834E+105, 3.3078854415193856E+107,
    2.4809140811395391E+109, 1.8854947016660498E+111, 1.4518309202828584E+113,
    1.1324281178206295E+115, 8.9461821307829729E+116, 7.1569457046263779E+118,
    5.7971260207473655E+120, 4.75364333701284E+122, 3.9455239697206569E+124,
    3.314240134565352E+126, 2.8171041143805494E+128, 2.4227095383672724E+130,
    2.1077572983795269E+132, 1.8548264225739836E+134, 1.6507955160908452E+136,
    1.4857159644817607E+138, 1.3520015276784023E+140, 1.24384140546413E+142,
    1.1567725070816409E+144, 1.0873661566567424E+146, 1.0329978488239052E+148,
    9.916779348709491E+149, 9.6192759682482062E+151, 9.426890448883242E+153,
    9.33262154439441E+155, 9.33262154439441E+157, 9.4259477598383536E+159,
    9.6144667150351211E+161, 9.9029007164861754E+163, 1.0299016745145622E+166,
    1.0813967582402903E+168, 1.1462805637347078E+170, 1.2265202031961373E+172,
    1.3246418194518284E+174, 1.4438595832024928E+176, 1.5882455415227421E+178,
    1.7629525510902437E+180, 1.9745068572210728E+182, 2.2311927486598123E+184,
    2.5435597334721862E+186, 2.9250936934930141E+188, 3.3931086844518965E+190,
    3.969937160808719E+192, 4.6845258497542883E+194, 5.5745857612076033E+196,
    6.6895029134491239E+198, 8.09429852527344E+200, 9.8750442008335976E+202,
    1.2146304367025325E+205, 1.5061417415111404E+207, 1.8826771768889254E+209,
    2.3721732428800459E+211, 3.0126600184576582E+213, 3.8562048236258025E+215,
    4.9745042224772855E+217, 6.4668554892204716E+219, 8.4715806908788174E+221,
    1.1182486511960039E+224, 1.4872707060906852E+226, 1.9929427461615181E+228,
    2.6904727073180495E+230, 3.6590428819525472E+232, 5.01288874827499E+234,
    6.9177864726194859E+236, 9.6157231969410859E+238, 1.346201247571752E+241,
    1.89814375907617E+243, 2.6953641378881614E+245, 3.8543707171800706E+247,
    5.5502938327393013E+249, 8.0479260574719866E+251, 1.17499720439091E+254,
    1.7272458904546376E+256, 2.5563239178728637E+258, 3.8089226376305671E+260,
    5.7133839564458505E+262, 8.6272097742332346E+264, 1.3113358856834518E+267,
    2.0063439050956811E+269, 3.0897696138473489E+271, 4.7891429014633912E+273,
    7.47106292628289E+275, 1.1729568794264138E+278, 1.8532718694937338E+280,
    2.9467022724950369E+282, 4.714723635992059E+284, 7.5907050539472148E+286,
    1.2296942187394488E+289, 2.0044015765453015E+291, 3.2872185855342945E+293,
    5.423910666131586E+295, 9.0036917057784329E+297, 1.5036165148649983E+300,
    2.5260757449731969E+302, 4.2690680090047027E+304, 7.257415615307994E+306 };

  static const uint8_T g[5] = { 3U, 5U, 7U, 9U, 13U };

  static const real_T theta[5] = { 0.01495585217958292, 0.253939833006323,
    0.95041789961629319, 2.097847961257068, 5.3719203511481517 };

  int8_T A_data[9];
  real_T P_data_0[16];
  real_T So_inv_data_0[12];
  b_cell_wrap_1_Experiment_2__n_T reshapes_0[2];
  real_T y_data_0[16];
  int32_T e;
  real_T c_y_data_0[16];
  int8_T e_y_data_0[3];
  real_T f_y_data_0[9];
  real_T b_n_data_0[3];
  int32_T b_coffset;
  real_T z1_data_0[3];
  int8_T b_data[3];
  static const real_T c_0[170] = { 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0,
    40320.0, 362880.0, 3.6288E+6, 3.99168E+7, 4.790016E+8, 6.2270208E+9,
    8.71782912E+10, 1.307674368E+12, 2.0922789888E+13, 3.55687428096E+14,
    6.402373705728E+15, 1.21645100408832E+17, 2.43290200817664E+18,
    5.109094217170944E+19, 1.1240007277776077E+21, 2.5852016738884978E+22,
    6.2044840173323941E+23, 1.5511210043330986E+25, 4.0329146112660565E+26,
    1.0888869450418352E+28, 3.0488834461171384E+29, 8.8417619937397008E+30,
    2.6525285981219103E+32, 8.2228386541779224E+33, 2.6313083693369352E+35,
    8.6833176188118859E+36, 2.9523279903960412E+38, 1.0333147966386144E+40,
    3.7199332678990118E+41, 1.3763753091226343E+43, 5.23022617466601E+44,
    2.0397882081197442E+46, 8.1591528324789768E+47, 3.3452526613163803E+49,
    1.4050061177528798E+51, 6.0415263063373834E+52, 2.6582715747884485E+54,
    1.1962222086548019E+56, 5.5026221598120885E+57, 2.5862324151116818E+59,
    1.2413915592536073E+61, 6.0828186403426752E+62, 3.0414093201713376E+64,
    1.5511187532873822E+66, 8.0658175170943877E+67, 4.2748832840600255E+69,
    2.3084369733924138E+71, 1.2696403353658276E+73, 7.1099858780486348E+74,
    4.0526919504877221E+76, 2.3505613312828789E+78, 1.3868311854568986E+80,
    8.3209871127413916E+81, 5.0758021387722484E+83, 3.1469973260387939E+85,
    1.98260831540444E+87, 1.2688693218588417E+89, 8.2476505920824715E+90,
    5.4434493907744307E+92, 3.6471110918188683E+94, 2.4800355424368305E+96,
    1.711224524281413E+98, 1.197857166996989E+100, 8.5047858856786218E+101,
    6.1234458376886077E+103, 4.4701154615126834E+105, 3.3078854415193856E+107,
    2.4809140811395391E+109, 1.8854947016660498E+111, 1.4518309202828584E+113,
    1.1324281178206295E+115, 8.9461821307829729E+116, 7.1569457046263779E+118,
    5.7971260207473655E+120, 4.75364333701284E+122, 3.9455239697206569E+124,
    3.314240134565352E+126, 2.8171041143805494E+128, 2.4227095383672724E+130,
    2.1077572983795269E+132, 1.8548264225739836E+134, 1.6507955160908452E+136,
    1.4857159644817607E+138, 1.3520015276784023E+140, 1.24384140546413E+142,
    1.1567725070816409E+144, 1.0873661566567424E+146, 1.0329978488239052E+148,
    9.916779348709491E+149, 9.6192759682482062E+151, 9.426890448883242E+153,
    9.33262154439441E+155, 9.33262154439441E+157, 9.4259477598383536E+159,
    9.6144667150351211E+161, 9.9029007164861754E+163, 1.0299016745145622E+166,
    1.0813967582402903E+168, 1.1462805637347078E+170, 1.2265202031961373E+172,
    1.3246418194518284E+174, 1.4438595832024928E+176, 1.5882455415227421E+178,
    1.7629525510902437E+180, 1.9745068572210728E+182, 2.2311927486598123E+184,
    2.5435597334721862E+186, 2.9250936934930141E+188, 3.3931086844518965E+190,
    3.969937160808719E+192, 4.6845258497542883E+194, 5.5745857612076033E+196,
    6.6895029134491239E+198, 8.09429852527344E+200, 9.8750442008335976E+202,
    1.2146304367025325E+205, 1.5061417415111404E+207, 1.8826771768889254E+209,
    2.3721732428800459E+211, 3.0126600184576582E+213, 3.8562048236258025E+215,
    4.9745042224772855E+217, 6.4668554892204716E+219, 8.4715806908788174E+221,
    1.1182486511960039E+224, 1.4872707060906852E+226, 1.9929427461615181E+228,
    2.6904727073180495E+230, 3.6590428819525472E+232, 5.01288874827499E+234,
    6.9177864726194859E+236, 9.6157231969410859E+238, 1.346201247571752E+241,
    1.89814375907617E+243, 2.6953641378881614E+245, 3.8543707171800706E+247,
    5.5502938327393013E+249, 8.0479260574719866E+251, 1.17499720439091E+254,
    1.7272458904546376E+256, 2.5563239178728637E+258, 3.8089226376305671E+260,
    5.7133839564458505E+262, 8.6272097742332346E+264, 1.3113358856834518E+267,
    2.0063439050956811E+269, 3.0897696138473489E+271, 4.7891429014633912E+273,
    7.47106292628289E+275, 1.1729568794264138E+278, 1.8532718694937338E+280,
    2.9467022724950369E+282, 4.714723635992059E+284, 7.5907050539472148E+286,
    1.2296942187394488E+289, 2.0044015765453015E+291, 3.2872185855342945E+293,
    5.423910666131586E+295, 9.0036917057784329E+297, 1.5036165148649983E+300,
    2.5260757449731969E+302, 4.2690680090047027E+304, 7.257415615307994E+306 };

  static const uint8_T g_0[5] = { 3U, 5U, 7U, 9U, 13U };

  static const real_T theta_0[5] = { 0.01495585217958292, 0.253939833006323,
    0.95041789961629319, 2.097847961257068, 5.3719203511481517 };

  real_T rtb_Clock1;
  real_T y_data_1[16];
  real_T tmp_data[11];
  real_T y_data_2[9];
  int8_T A_data_0[6];
  int32_T P_size[2];
  int32_T y_size[2];
  int32_T c_y_size[2];
  int32_T d_y_size[2];
  int32_T P_size_0[2];
  int32_T y_size_0[2];
  int32_T c_y_size_0[2];
  int32_T y_size_1[2];
  int32_T y_size_2[2];
  real_T tmp;
  b_cell_wrap_1_Experiment_2__n_T reshapes_1;
  int32_T c_y_data_tmp;
  boolean_T exitg1;

  /* S-Function (hil_read_encoder_block): '<S3>/HIL Read Encoder' */

  /* S-Function Block: Experiment_2_moving_cart/Pendulum/Quanser IO/HIL Read Encoder (hil_read_encoder_block) */
  {
    t_error result = hil_read_encoder
      (Experiment_2_moving_cart_DW.HILInitialize_Card,
       Experiment_2_moving_cart_P.HILReadEncoder_channels, 2,
       &Experiment_2_moving_cart_DW.HILReadEncoder_Buffer[0]);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
    } else {
      rtb_HILReadEncoder_o1 = Experiment_2_moving_cart_DW.HILReadEncoder_Buffer
        [0];
      rtb_HILReadEncoder_o2 = Experiment_2_moving_cart_DW.HILReadEncoder_Buffer
        [1];
    }
  }

  /* Gain: '<S3>/Gain1' */
  Experiment_2_moving_cart_B.Gain1 = Experiment_2_moving_cart_P.gain_cart *
    rtb_HILReadEncoder_o2;

  /* Sum: '<S4>/Sum' incorporates:
   *  Constant: '<S4>/Constant'
   *  Gain: '<S3>/Gain'
   */
  Experiment_2_moving_cart_B.Sum = Experiment_2_moving_cart_P.gain_pendulum *
    rtb_HILReadEncoder_o1 + Experiment_2_moving_cart_P.Constant_Value;

  /* MATLAB Function: '<S4>/Embedded MATLAB Function' */
  Experiment_2_moving_cart_B.Phi = Experiment_2_moving_cart_B.Sum;
  while (Experiment_2_moving_cart_B.Phi > 6.2831853071795862) {
    Experiment_2_moving_cart_B.Phi -= 6.2831853071795862;
  }

  while (Experiment_2_moving_cart_B.Phi < 0.0) {
    Experiment_2_moving_cart_B.Phi += 6.2831853071795862;
  }

  /* End of MATLAB Function: '<S4>/Embedded MATLAB Function' */

  /* Outputs for Atomic SubSystem: '<S4>/Differentiator1' */
  /* MATLAB Function: '<S6>/Differentiator' incorporates:
   *  SampleTimeMath: '<S6>/Weighted Sample Time'
   *
   * About '<S6>/Weighted Sample Time':
   *  y = K where K = ( w * Ts )
   */
  if (!Experiment_2_moving_cart_DW.zp_not_empty) {
    Experiment_2_moving_cart_DW.zp_p.size = 2;
    Experiment_2_moving_cart_DW.zp_p.data[0] = 0.0;
    Experiment_2_moving_cart_DW.zp_p.data[1] = 0.0;
    Experiment_2_moving_cart_DW.zp_not_empty = true;
    A_data_0[0] = 0;
    A_data_0[1] = 0;
    A_data_0[2] = 1;
    A_data_0[3] = 0;
    A_data_0[4] = 0;
    A_data_0[5] = 0;
    for (b_m = 0; b_m < 3; b_m++) {
      for (coffset = 0; coffset < 2; coffset++) {
        reshapes[0].f1.data[coffset + (b_m << 1)] = A_data_0[(b_m << 1) +
          coffset];
      }
    }

    y_size[0] = 3;
    y_size[1] = 3;
    for (b_m = 0; b_m < 3; b_m++) {
      for (coffset = 0; coffset < 2; coffset++) {
        y_data[coffset + 3 * b_m] = reshapes[0].f1.data[(b_m << 1) + coffset] *
          Experiment_2_moving_cart_P.WeightedSampleTime_WtEt_b;
      }
    }

    for (b_m = 0; b_m < 3; b_m++) {
      y_data[2 + 3 * b_m] = Experiment_2_moving_cart_P.WeightedSampleTime_WtEt_b
        * 0.0;
    }

    P_size[0] = 3;
    normA = Experiment_2_moving_cart_norm_j(y_data, y_size);
    if (normA <= 5.3719203511481517) {
      b_coffset = 0;
      exitg1 = false;
      while ((!exitg1) && (b_coffset < 5)) {
        if (normA <= theta[b_coffset]) {
          Exper_PadeApproximantOfDegree_k(y_data, y_size, g[b_coffset], P_data,
            P_size);
          exitg1 = true;
        } else {
          b_coffset++;
        }
      }
    } else {
      rtb_Clock1 = normA / 5.3719203511481517;
      if ((!rtIsInf(rtb_Clock1)) && (!rtIsNaN(rtb_Clock1))) {
        rtb_Clock1 = frexp(rtb_Clock1, &b_coffset);
      } else {
        b_coffset = 0;
      }

      normA = b_coffset;
      if (rtb_Clock1 == 0.5) {
        normA = (real_T)b_coffset - 1.0;
      }

      rtb_Clock1 = rt_powd_snf(2.0, normA);
      y_size_2[0] = 3;
      y_size_2[1] = 3;
      for (b_m = 0; b_m < 9; b_m++) {
        y_data_2[b_m] = y_data[b_m] / rtb_Clock1;
      }

      Expe_PadeApproximantOfDegree_kh(y_data_2, y_size_2, P_data, P_size);
      for (j = 0; j < (int32_T)normA; j++) {
        if (P_size[0] == 1) {
          c_y_size[0] = 1;
          for (b_m = 0; b_m < 1; b_m++) {
            for (coffset = 0; coffset < 3; coffset++) {
              c_y_data[coffset] = P_data[2 + coffset] * P_data[2] + (P_data[1 +
                coffset] * P_data[1] + P_data[0] * P_data[coffset]);
            }
          }
        } else {
          b_m = P_size[0];
          c_y_size[0] = P_size[0];
          for (b_j = 0; b_j < 3; b_j++) {
            coffset = b_j * b_m - 1;
            boffset = b_j * 3;
            for (b_coffset = 1; b_coffset <= b_m; b_coffset++) {
              c_y_data[coffset + b_coffset] = 0.0;
            }

            for (b_coffset = 0; b_coffset < 3; b_coffset++) {
              if (P_data[boffset + b_coffset] != 0.0) {
                aoffset = b_coffset * b_m;
                for (c_i = 1; c_i <= b_m; c_i++) {
                  c_y_data_tmp = coffset + c_i;
                  c_y_data[c_y_data_tmp] += P_data[(aoffset + c_i) - 1] *
                    P_data[boffset + b_coffset];
                }
              }
            }
          }
        }

        P_size[0] = c_y_size[0];
        boffset = c_y_size[0] * 3 - 1;
        if (0 <= boffset) {
          memcpy(&P_data[0], &c_y_data[0], (boffset + 1) * sizeof(real_T));
        }
      }
    }

    Experiment_2_moving_cart_DW.Phi_b.size[0] = 2;
    Experiment_2_moving_cart_DW.Phi_b.size[1] = 2;
    Experiment_2_moving_cart_DW.bD_e.size = 2;
    Experiment_2_moving_cart_DW.bD_e.data[0] = P_data[P_size[0] << 1];
    Experiment_2_moving_cart_DW.Phi_b.data[0] = P_data[0];
    d_y_data[0] = 0.0;
    f_y_data[0] = 0;
    Experiment_2_moving_cart_DW.Phi_b.data[1] = P_data[1];
    d_y_data[2] = 0.0;
    f_y_data[2] = 1;
    Experiment_2_moving_cart_DW.bD_e.data[1] = P_data[(P_size[0] << 1) + 1];
    Experiment_2_moving_cart_DW.Phi_b.data[Experiment_2_moving_cart_DW.Phi_b.size
      [0]] = P_data[P_size[0]];
    d_y_data[1] = 1.0;
    f_y_data[1] = 0;
    Experiment_2_moving_cart_DW.Phi_b.data[1 +
      Experiment_2_moving_cart_DW.Phi_b.size[0]] = P_data[1 + P_size[0]];
    d_y_data[3] = 1.0;
    f_y_data[3] = 1;
    for (b_coffset = 0; b_coffset < 4; b_coffset++) {
      g_y_data[b_coffset] = rt_powd_snf(d_y_data[b_coffset], (real_T)
        f_y_data[b_coffset]);
    }

    b_n_data[0] = 0.0;
    b_n_data[1] = 1.0;
    for (b_coffset = 0; b_coffset < 2; b_coffset++) {
      normA = b_n_data[b_coffset];
      if ((!(b_n_data[b_coffset] != b_n_data[b_coffset])) && (!rtIsInf
           (b_n_data[b_coffset]))) {
        if (b_n_data[b_coffset] > 170.0) {
          normA = (rtInf);
        } else if (b_n_data[b_coffset] < 1.0) {
          normA = 1.0;
        } else {
          normA = c[(int32_T)b_n_data[b_coffset] - 1];
        }
      }

      b_n_data[b_coffset] = normA;
    }

    for (b_m = 0; b_m < 2; b_m++) {
      b_n_data[b_m] = 1.0 / b_n_data[b_m];
    }

    for (b_m = 0; b_m < 4; b_m++) {
      d_y_data[b_m] = 0.0;
    }

    for (j = 0; j < 2; j++) {
      d_y_data[j + (j << 1)] = b_n_data[j];
    }

    c_y_size[0] = 2;
    c_y_size[1] = 2;
    for (b_j = 0; b_j < 2; b_j++) {
      b_coffset = (b_j << 1) - 1;
      b_m = b_j << 1;
      c_y_data[b_coffset + 1] = 0.0;
      c_y_data[b_coffset + 2] = 0.0;
      for (coffset = 0; coffset < 2; coffset++) {
        if (d_y_data[b_m + coffset] != 0.0) {
          boffset = coffset << 1;
          c_y_data[b_coffset + 1] += d_y_data[b_m + coffset] * g_y_data[boffset];
          c_y_data[b_coffset + 2] += d_y_data[b_m + coffset] * g_y_data[boffset
            + 1];
        }
      }
    }

    e_y_data[0] = 0;
    e_y_data[1] = 1;
    for (b_m = 0; b_m < 2; b_m++) {
      z1_data[b_m] = b_n_data[b_m];
    }

    for (b_coffset = 0; b_coffset < 2; b_coffset++) {
      z1_data[b_coffset] = rt_powd_snf
        (Experiment_2_moving_cart_P.WeightedSampleTime_WtEt_b, (real_T)
         e_y_data[b_coffset]);
    }

    for (b_m = 0; b_m < 2; b_m++) {
      z1_data[b_m] = 1.0 / z1_data[b_m];
    }

    d_y_size[0] = 2;
    d_y_size[1] = 2;
    for (b_m = 0; b_m < 4; b_m++) {
      d_y_data[b_m] = 0.0;
    }

    for (b_j = 0; b_j < 2; b_j++) {
      d_y_data[b_j + (b_j << 1)] = z1_data[b_j];
    }

    Experiment_2_moving__mrdivide_f(d_y_data, d_y_size, c_y_data, c_y_size,
      So_inv_data, y_size);
    b_b_data[0] = 0;
    b_b_data[1] = 1;
    if (y_size[1] == 1) {
      Experiment_2_moving_cart_DW.So_inv_lastCol_e.size = y_size[0];
      boffset = y_size[0];
      for (b_m = 0; b_m < boffset; b_m++) {
        Experiment_2_moving_cart_DW.So_inv_lastCol_e.data[b_m] = 0.0;
        for (coffset = 0; coffset < 1; coffset++) {
          Experiment_2_moving_cart_DW.So_inv_lastCol_e.data[b_m] +=
            So_inv_data[b_m] * 0.0;
        }
      }
    } else {
      b_m = y_size[0];
      Experiment_2_moving_cart_DW.So_inv_lastCol_e.size = y_size[0];
      for (b_coffset = 1; b_coffset <= b_m; b_coffset++) {
        Experiment_2_moving_cart_DW.So_inv_lastCol_e.data[b_coffset - 1] = 0.0;
      }

      for (b_j = 0; b_j < y_size[1]; b_j++) {
        if (b_b_data[b_j] != 0) {
          coffset = b_j * b_m;
          for (boffset = 0; boffset < b_m; boffset++) {
            Experiment_2_moving_cart_DW.So_inv_lastCol_e.data[boffset] +=
              So_inv_data[coffset + boffset] * (real_T)b_b_data[b_j];
          }
        }
      }
    }
  }

  normA = Experiment_2_moving_cart_B.Gain1 -
    Experiment_2_moving_cart_DW.zp_p.data[0];
  Experiment_2_moving_cart_B.z = Experiment_2_moving_cart_DW.zp_p.data[1];

  /* SampleTimeMath: '<S6>/Weighted Sample Time' incorporates:
   *  MATLAB Function: '<S6>/Differentiator'
   *
   * About '<S6>/Weighted Sample Time':
   *  y = K where K = ( w * Ts )
   */
  Experime_ackerman_precomputed_g
    (Experiment_2_moving_cart_P.WeightedSampleTime_WtEt_b,
     Experim_disc_eigenvalues_URED_d(normA,
      Experiment_2_moving_cart_P.Differentiator1_mu,
      Experiment_2_moving_cart_P.WeightedSampleTime_WtEt_b), tmp_data, &b_m);

  /* MATLAB Function: '<S6>/Differentiator' */
  b_n_data[0] = Experiment_2_moving_cart_DW.zp_p.data[0];
  b_n_data[1] = Experiment_2_moving_cart_DW.zp_p.data[1];
  Experiment_2_moving_cart_step_n(normA, b_n_data,
    Experiment_2_moving_cart_B.Gain1, Experiment_2_moving_cart_DW.Phi_b.data,
    Experiment_2_moving_cart_DW.bD_e.data, tmp_data,
    Experiment_2_moving_cart_DW.zp_p.data,
    &Experiment_2_moving_cart_DW.zp_p.size);
  Experiment_2_moving_cart_DW.zp_not_empty = true;
  Experiment_2_moving_cart_B.x0 = normA;

  /* End of Outputs for SubSystem: '<S4>/Differentiator1' */

  /* Outputs for Atomic SubSystem: '<S4>/Differentiator' */
  /* MATLAB Function: '<S5>/Differentiator' incorporates:
   *  SampleTimeMath: '<S5>/Weighted Sample Time'
   *
   * About '<S5>/Weighted Sample Time':
   *  y = K where K = ( w * Ts )
   */
  if (!Experiment_2_moving_cart_DW.zp_not_empty_n) {
    Experiment_2_moving_cart_DW.zp.size = 3;
    Experiment_2_moving_cart_DW.zp.data[0] = 0.0;
    Experiment_2_moving_cart_DW.zp.data[1] = 0.0;
    Experiment_2_moving_cart_DW.zp.data[2] = 0.0;
    Experiment_2_moving_cart_DW.zp_not_empty_n = true;
    for (b_m = 0; b_m < 9; b_m++) {
      A_data[b_m] = 0;
    }

    A_data[3] = 1;
    A_data[7] = 1;
    for (b_m = 0; b_m < 3; b_m++) {
      y_data[3 * b_m] = A_data[3 * b_m];
      j = 3 * b_m + 1;
      b_j = 1 + 3 * b_m;
      y_data[b_j] = A_data[j];
      coffset = 3 * b_m + 2;
      boffset = 2 + 3 * b_m;
      y_data[boffset] = A_data[coffset];
      So_inv_data_0[3 * b_m] = y_data[3 * b_m];
      So_inv_data_0[b_j] = y_data[j];
      So_inv_data_0[boffset] = y_data[coffset];
    }

    So_inv_data_0[9] = 0.0;
    So_inv_data_0[10] = 0.0;
    So_inv_data_0[11] = 0.0;
    y_size_0[0] = 4;
    y_size_0[1] = 4;
    for (b_m = 0; b_m < 3; b_m++) {
      reshapes_1 = reshapes_0[0];
      reshapes_1.f1.data[b_m] = So_inv_data_0[b_m];
      reshapes_1.f1.data[b_m + 3] = So_inv_data_0[b_m + 3];
      reshapes_1.f1.data[b_m + 6] = So_inv_data_0[b_m + 6];
      reshapes_1.f1.data[b_m + 9] = So_inv_data_0[b_m + 9];
      y_data_0[b_m] = Experiment_2_moving_cart_P.WeightedSampleTime_WtEt *
        reshapes_1.f1.data[b_m];
      y_data_0[b_m + 4] = reshapes_1.f1.data[b_m + 3] *
        Experiment_2_moving_cart_P.WeightedSampleTime_WtEt;
      y_data_0[b_m + 8] = reshapes_1.f1.data[b_m + 6] *
        Experiment_2_moving_cart_P.WeightedSampleTime_WtEt;
      y_data_0[b_m + 12] = reshapes_1.f1.data[b_m + 9] *
        Experiment_2_moving_cart_P.WeightedSampleTime_WtEt;
      reshapes_0[0] = reshapes_1;
    }

    y_data_0[3] = Experiment_2_moving_cart_P.WeightedSampleTime_WtEt * 0.0;
    y_data_0[7] = Experiment_2_moving_cart_P.WeightedSampleTime_WtEt * 0.0;
    y_data_0[11] = Experiment_2_moving_cart_P.WeightedSampleTime_WtEt * 0.0;
    y_data_0[15] = Experiment_2_moving_cart_P.WeightedSampleTime_WtEt * 0.0;
    P_size_0[0] = 4;
    normA = Experiment_2_moving_cart_norm(y_data_0, y_size_0);
    if (normA <= 5.3719203511481517) {
      b_coffset = 0;
      exitg1 = false;
      while ((!exitg1) && (b_coffset < 5)) {
        if (normA <= theta_0[b_coffset]) {
          Experim_PadeApproximantOfDegree(y_data_0, y_size_0, g_0[b_coffset],
            P_data_0, P_size_0);
          exitg1 = true;
        } else {
          b_coffset++;
        }
      }
    } else {
      rtb_Clock1 = normA / 5.3719203511481517;
      if ((!rtIsInf(rtb_Clock1)) && (!rtIsNaN(rtb_Clock1))) {
        rtb_Clock1 = frexp(rtb_Clock1, &e);
      } else {
        e = 0;
      }

      normA = e;
      if (rtb_Clock1 == 0.5) {
        normA = (real_T)e - 1.0;
      }

      rtb_Clock1 = rt_powd_snf(2.0, normA);
      y_size_1[0] = 4;
      y_size_1[1] = 4;
      for (b_m = 0; b_m < 16; b_m++) {
        y_data_1[b_m] = y_data_0[b_m] / rtb_Clock1;
      }

      Exper_PadeApproximantOfDegree_e(y_data_1, y_size_1, P_data_0, P_size_0);
      for (b_j = 0; b_j < (int32_T)normA; b_j++) {
        if (P_size_0[0] == 1) {
          c_y_size_0[0] = 1;
          for (b_m = 0; b_m < 1; b_m++) {
            for (coffset = 0; coffset < 4; coffset++) {
              rtb_Clock1 = P_data_0[3 + coffset] * P_data_0[3] + (P_data_0[2 +
                coffset] * P_data_0[2] + (P_data_0[1 + coffset] * P_data_0[1] +
                P_data_0[0] * P_data_0[coffset]));
              c_y_data_0[coffset] = rtb_Clock1;
            }
          }
        } else {
          b_m = P_size_0[0];
          c_y_size_0[0] = P_size_0[0];
          for (j = 0; j < 4; j++) {
            coffset = j * b_m - 1;
            boffset = j << 2;
            for (b_coffset = 1; b_coffset <= b_m; b_coffset++) {
              c_y_data_0[coffset + b_coffset] = 0.0;
            }

            for (b_coffset = 0; b_coffset < 4; b_coffset++) {
              if (P_data_0[boffset + b_coffset] != 0.0) {
                aoffset = b_coffset * b_m;
                for (c_i = 1; c_i <= b_m; c_i++) {
                  c_y_data_tmp = coffset + c_i;
                  c_y_data_0[c_y_data_tmp] += P_data_0[(aoffset + c_i) - 1] *
                    P_data_0[boffset + b_coffset];
                }
              }
            }
          }
        }

        P_size_0[0] = c_y_size_0[0];
        boffset = (c_y_size_0[0] << 2) - 1;
        if (0 <= boffset) {
          memcpy(&P_data_0[0], &c_y_data_0[0], (boffset + 1) * sizeof(real_T));
        }
      }
    }

    Experiment_2_moving_cart_DW.Phi.size[0] = 3;
    Experiment_2_moving_cart_DW.Phi.size[1] = 3;
    Experiment_2_moving_cart_DW.bD.size = 3;
    for (b_m = 0; b_m < 3; b_m++) {
      Experiment_2_moving_cart_DW.bD.data[b_m] = P_data_0[P_size_0[0] * 3 + b_m];
      coffset = P_size_0[0] * b_m;
      e = Experiment_2_moving_cart_DW.Phi.size[0] * b_m;
      Experiment_2_moving_cart_DW.Phi.data[e] = P_data_0[coffset];
      c_y_data[b_m] = b_m;
      f_y_data_0[b_m] = 0.0;
      Experiment_2_moving_cart_DW.Phi.data[1 + e] = P_data_0[coffset + 1];
      c_y_data[b_m + 3] = b_m;
      f_y_data_0[b_m + 3] = 1.0;
      Experiment_2_moving_cart_DW.Phi.data[2 + e] = P_data_0[coffset + 2];
      c_y_data[b_m + 6] = b_m;
      f_y_data_0[b_m + 6] = 2.0;
    }

    for (b_coffset = 0; b_coffset < 9; b_coffset++) {
      y_data[b_coffset] = rt_powd_snf(c_y_data[b_coffset], f_y_data_0[b_coffset]);
    }

    b_n_data_0[0] = 0.0;
    b_n_data_0[1] = 1.0;
    b_n_data_0[2] = 2.0;
    for (b_coffset = 0; b_coffset < 3; b_coffset++) {
      normA = b_n_data_0[b_coffset];
      if ((!(b_n_data_0[b_coffset] != b_n_data_0[b_coffset])) && (!rtIsInf
           (b_n_data_0[b_coffset]))) {
        if (b_n_data_0[b_coffset] > 170.0) {
          normA = (rtInf);
        } else if (b_n_data_0[b_coffset] < 1.0) {
          normA = 1.0;
        } else {
          normA = c_0[(int32_T)b_n_data_0[b_coffset] - 1];
        }
      }

      b_n_data_0[b_coffset] = normA;
    }

    for (b_m = 0; b_m < 3; b_m++) {
      b_n_data_0[b_m] = 1.0 / b_n_data_0[b_m];
    }

    memset(&c_y_data[0], 0, 9U * sizeof(real_T));
    for (b_j = 0; b_j < 3; b_j++) {
      c_y_data[b_j + 3 * b_j] = b_n_data_0[b_j];
    }

    c_y_size_0[0] = 3;
    c_y_size_0[1] = 3;
    for (b_j = 0; b_j < 3; b_j++) {
      b_coffset = b_j * 3 - 1;
      b_m = b_j * 3;
      c_y_data_0[b_coffset + 1] = 0.0;
      c_y_data_0[b_coffset + 2] = 0.0;
      c_y_data_0[b_coffset + 3] = 0.0;
      for (coffset = 0; coffset < 3; coffset++) {
        if (c_y_data[b_m + coffset] != 0.0) {
          boffset = coffset * 3;
          c_y_data_0[b_coffset + 1] += c_y_data[b_m + coffset] * y_data[boffset];
          c_y_data_0[b_coffset + 2] += c_y_data[b_m + coffset] * y_data[boffset
            + 1];
          c_y_data_0[b_coffset + 3] += c_y_data[b_m + coffset] * y_data[boffset
            + 2];
        }
      }
    }

    e_y_data_0[0] = 0;
    e_y_data_0[1] = 1;
    e_y_data_0[2] = 2;
    for (b_m = 0; b_m < 3; b_m++) {
      z1_data_0[b_m] = b_n_data_0[b_m];
    }

    for (b_coffset = 0; b_coffset < 3; b_coffset++) {
      z1_data_0[b_coffset] = rt_powd_snf
        (Experiment_2_moving_cart_P.WeightedSampleTime_WtEt, (real_T)
         e_y_data_0[b_coffset]);
    }

    for (b_m = 0; b_m < 3; b_m++) {
      z1_data_0[b_m] = 1.0 / z1_data_0[b_m];
    }

    c_y_size[0] = 3;
    c_y_size[1] = 3;
    memset(&c_y_data[0], 0, 9U * sizeof(real_T));
    for (e = 0; e < 3; e++) {
      c_y_data[e + 3 * e] = z1_data_0[e];
    }

    Experiment_2_moving_ca_mrdivide(c_y_data, c_y_size, c_y_data_0, c_y_size_0,
      So_inv_data_0, y_size);
    b_data[0] = 0;
    b_data[1] = 0;
    b_data[2] = 1;
    if (y_size[1] == 1) {
      Experiment_2_moving_cart_DW.So_inv_lastCol.size = y_size[0];
      boffset = y_size[0];
      for (b_m = 0; b_m < boffset; b_m++) {
        Experiment_2_moving_cart_DW.So_inv_lastCol.data[b_m] = 0.0;
        for (coffset = 0; coffset < 1; coffset++) {
          Experiment_2_moving_cart_DW.So_inv_lastCol.data[b_m] +=
            So_inv_data_0[b_m] * 0.0;
        }
      }
    } else {
      b_m = y_size[0];
      Experiment_2_moving_cart_DW.So_inv_lastCol.size = y_size[0];
      for (b_coffset = 1; b_coffset <= b_m; b_coffset++) {
        Experiment_2_moving_cart_DW.So_inv_lastCol.data[b_coffset - 1] = 0.0;
      }

      for (b_j = 0; b_j < y_size[1]; b_j++) {
        if (b_data[b_j] != 0) {
          coffset = b_j * b_m;
          for (boffset = 0; boffset < b_m; boffset++) {
            Experiment_2_moving_cart_DW.So_inv_lastCol.data[boffset] +=
              So_inv_data_0[coffset + boffset] * (real_T)b_data[b_j];
          }
        }
      }
    }
  }

  normA = Experiment_2_moving_cart_B.Sum - Experiment_2_moving_cart_DW.zp.data[0];
  Experiment_2_moving_cart_B.z_h[0] = Experiment_2_moving_cart_DW.zp.data[1];
  Experiment_2_moving_cart_B.z_h[1] = Experiment_2_moving_cart_DW.zp.data[2];

  /* SampleTimeMath: '<S5>/Weighted Sample Time' incorporates:
   *  MATLAB Function: '<S5>/Differentiator'
   *
   * About '<S5>/Weighted Sample Time':
   *  y = K where K = ( w * Ts )
   */
  Experiment_ackerman_precomputed
    (Experiment_2_moving_cart_P.WeightedSampleTime_WtEt,
     Experimen_disc_eigenvalues_URED(normA,
      Experiment_2_moving_cart_P.Differentiator_mu,
      Experiment_2_moving_cart_P.WeightedSampleTime_WtEt), tmp_data, &b_m);

  /* MATLAB Function: '<S5>/Differentiator' */
  b_n_data_0[0] = Experiment_2_moving_cart_DW.zp.data[0];
  b_n_data_0[1] = Experiment_2_moving_cart_DW.zp.data[1];
  b_n_data_0[2] = Experiment_2_moving_cart_DW.zp.data[2];
  Experiment_2_moving_cart_step_p(normA, b_n_data_0,
    Experiment_2_moving_cart_B.Sum, Experiment_2_moving_cart_DW.Phi.data,
    Experiment_2_moving_cart_DW.bD.data, tmp_data,
    Experiment_2_moving_cart_DW.zp.data, &Experiment_2_moving_cart_DW.zp.size);
  Experiment_2_moving_cart_DW.zp_not_empty_n = true;
  Experiment_2_moving_cart_B.x0_n = normA;

  /* End of Outputs for SubSystem: '<S4>/Differentiator' */
  /* TransportDelay: '<Root>/Transport Delay' */
  {
    real_T **uBuffer = (real_T**)
      &Experiment_2_moving_cart_DW.TransportDelay_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &Experiment_2_moving_cart_DW.TransportDelay_PWORK.TUbufferPtrs[1];
    real_T simTime = Experiment_2_moving_cart_M->Timing.t[0];
    real_T tMinusDelay = simTime -
      Experiment_2_moving_cart_P.TransportDelay_Delay;
    Experiment_2_moving_cart_B.u = rt_TDelayInterpolate(
      tMinusDelay,
      0.0,
      *tBuffer,
      *uBuffer,
      Experiment_2_moving_cart_DW.TransportDelay_IWORK.CircularBufSize,
      &Experiment_2_moving_cart_DW.TransportDelay_IWORK.Last,
      Experiment_2_moving_cart_DW.TransportDelay_IWORK.Tail,
      Experiment_2_moving_cart_DW.TransportDelay_IWORK.Head,
      Experiment_2_moving_cart_P.TransportDelay_InitOutput,
      0,
      0);
  }

  /* SignalConversion: '<Root>/TmpSignal ConversionAtTo Workspace1Inport1' */
  Experiment_2_moving_cart_B.TmpSignalConversionAtToWorkspac[0] =
    Experiment_2_moving_cart_B.Gain1;
  Experiment_2_moving_cart_B.TmpSignalConversionAtToWorkspac[1] =
    Experiment_2_moving_cart_B.Phi;
  Experiment_2_moving_cart_B.TmpSignalConversionAtToWorkspac[2] =
    Experiment_2_moving_cart_B.z;
  Experiment_2_moving_cart_B.TmpSignalConversionAtToWorkspac[3] =
    Experiment_2_moving_cart_B.z_h[0];

  /* Step: '<Root>/Step2' incorporates:
   *  Step: '<Root>/Step'
   *  Step: '<Root>/Step1'
   */
  normA = Experiment_2_moving_cart_M->Timing.t[0];

  /* Clock: '<S1>/Clock1' */
  rtb_Clock1 = Experiment_2_moving_cart_M->Timing.t[0];

  /* Gain: '<S1>/Gain' incorporates:
   *  Constant: '<S1>/deltaFreq'
   *  Constant: '<S1>/targetTime'
   *  Product: '<S1>/Product'
   */
  Experiment_2_moving_cart_B.Gain = (Experiment_2_moving_cart_P.ChirpSignal_f2 -
    Experiment_2_moving_cart_P.ChirpSignal_f1) * 6.2831853071795862 /
    Experiment_2_moving_cart_P.ChirpSignal_T
    * Experiment_2_moving_cart_P.Gain_Gain;

  /* ManualSwitch: '<Root>/Manual Switch' incorporates:
   *  Constant: '<S1>/initialFreq'
   *  Gain: '<Root>/Gain2'
   *  Product: '<S1>/Product1'
   *  Product: '<S1>/Product2'
   *  Step: '<Root>/Step'
   *  Sum: '<Root>/Add'
   *  Sum: '<S1>/Sum'
   *  Trigonometry: '<S1>/Output'
   */
  if (Experiment_2_moving_cart_P.ManualSwitch_CurrentSetting == 1) {
    Experiment_2_moving_cart_B.iA = sin((rtb_Clock1 *
      Experiment_2_moving_cart_B.Gain + 6.2831853071795862 *
      Experiment_2_moving_cart_P.ChirpSignal_f1) * rtb_Clock1) *
      Experiment_2_moving_cart_P.Gain2_Gain;
  } else {
    if (normA < Experiment_2_moving_cart_P.Step_Time) {
      /* Step: '<Root>/Step' */
      rtb_Clock1 = Experiment_2_moving_cart_P.Step_Y0;
    } else {
      /* Step: '<Root>/Step' */
      rtb_Clock1 = Experiment_2_moving_cart_P.Step_YFinal;
    }

    /* Step: '<Root>/Step1' */
    if (normA < Experiment_2_moving_cart_P.Step1_Time) {
      tmp = Experiment_2_moving_cart_P.Step1_Y0;
    } else {
      tmp = Experiment_2_moving_cart_P.Step1_YFinal;
    }

    /* Step: '<Root>/Step2' */
    if (normA < Experiment_2_moving_cart_P.Step2_Time) {
      normA = Experiment_2_moving_cart_P.Step2_Y0;
    } else {
      normA = Experiment_2_moving_cart_P.Step2_YFinal;
    }

    Experiment_2_moving_cart_B.iA = (rtb_Clock1 + tmp) + normA;
  }

  /* End of ManualSwitch: '<Root>/Manual Switch' */

  /* Saturate: '<S2>/Saturation' */
  if (Experiment_2_moving_cart_B.u > Experiment_2_moving_cart_P.imax) {
    Experiment_2_moving_cart_B.Saturation = Experiment_2_moving_cart_P.imax;
  } else if (Experiment_2_moving_cart_B.u < -Experiment_2_moving_cart_P.imax) {
    Experiment_2_moving_cart_B.Saturation = -Experiment_2_moving_cart_P.imax;
  } else {
    Experiment_2_moving_cart_B.Saturation = Experiment_2_moving_cart_B.u;
  }

  /* End of Saturate: '<S2>/Saturation' */

  /* S-Function (hil_write_analog_block): '<S3>/HIL Write Analog' */

  /* S-Function Block: Experiment_2_moving_cart/Pendulum/Quanser IO/HIL Write Analog (hil_write_analog_block) */
  {
    t_error result;
    result = hil_write_analog(Experiment_2_moving_cart_DW.HILInitialize_Card,
      &Experiment_2_moving_cart_P.HILWriteAnalog_channels, 1,
      &Experiment_2_moving_cart_B.Saturation);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
    }
  }
}

/* Model update function */
void Experiment_2_moving_cart_update(void)
{
  /* Update for TransportDelay: '<Root>/Transport Delay' */
  {
    real_T **uBuffer = (real_T**)
      &Experiment_2_moving_cart_DW.TransportDelay_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &Experiment_2_moving_cart_DW.TransportDelay_PWORK.TUbufferPtrs[1];
    real_T simTime = Experiment_2_moving_cart_M->Timing.t[0];
    boolean_T bufferisfull = false;
    Experiment_2_moving_cart_DW.TransportDelay_IWORK.Head =
      ((Experiment_2_moving_cart_DW.TransportDelay_IWORK.Head <
        (Experiment_2_moving_cart_DW.TransportDelay_IWORK.CircularBufSize-1)) ?
       (Experiment_2_moving_cart_DW.TransportDelay_IWORK.Head+1) : 0);
    if (Experiment_2_moving_cart_DW.TransportDelay_IWORK.Head ==
        Experiment_2_moving_cart_DW.TransportDelay_IWORK.Tail) {
      bufferisfull = true;
      Experiment_2_moving_cart_DW.TransportDelay_IWORK.Tail =
        ((Experiment_2_moving_cart_DW.TransportDelay_IWORK.Tail <
          (Experiment_2_moving_cart_DW.TransportDelay_IWORK.CircularBufSize-1)) ?
         (Experiment_2_moving_cart_DW.TransportDelay_IWORK.Tail+1) : 0);
    }

    (*tBuffer)[Experiment_2_moving_cart_DW.TransportDelay_IWORK.Head] = simTime;
    (*uBuffer)[Experiment_2_moving_cart_DW.TransportDelay_IWORK.Head] =
      Experiment_2_moving_cart_B.iA;
    if (bufferisfull) {
      rtsiSetBlockStateForSolverChangedAtMajorStep
        (&Experiment_2_moving_cart_M->solverInfo, true);
      rtsiSetContTimeOutputInconsistentWithStateAtMajorStep
        (&Experiment_2_moving_cart_M->solverInfo, true);
    }
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++Experiment_2_moving_cart_M->Timing.clockTick0)) {
    ++Experiment_2_moving_cart_M->Timing.clockTickH0;
  }

  Experiment_2_moving_cart_M->Timing.t[0] =
    Experiment_2_moving_cart_M->Timing.clockTick0 *
    Experiment_2_moving_cart_M->Timing.stepSize0 +
    Experiment_2_moving_cart_M->Timing.clockTickH0 *
    Experiment_2_moving_cart_M->Timing.stepSize0 * 4294967296.0;

  {
    /* Update absolute timer for sample time: [0.01s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++Experiment_2_moving_cart_M->Timing.clockTick1)) {
      ++Experiment_2_moving_cart_M->Timing.clockTickH1;
    }

    Experiment_2_moving_cart_M->Timing.t[1] =
      Experiment_2_moving_cart_M->Timing.clockTick1 *
      Experiment_2_moving_cart_M->Timing.stepSize1 +
      Experiment_2_moving_cart_M->Timing.clockTickH1 *
      Experiment_2_moving_cart_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Model initialize function */
void Experiment_2_moving_cart_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: Experiment_2_moving_cart/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0",
                      &Experiment_2_moving_cart_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options
      (Experiment_2_moving_cart_DW.HILInitialize_Card,
       "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(Experiment_2_moving_cart_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
      return;
    }

    if ((Experiment_2_moving_cart_P.HILInitialize_AIPStart && !is_switching) ||
        (Experiment_2_moving_cart_P.HILInitialize_AIPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums =
          &Experiment_2_moving_cart_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = (Experiment_2_moving_cart_P.HILInitialize_AILow);
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums =
          &Experiment_2_moving_cart_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = Experiment_2_moving_cart_P.HILInitialize_AIHigh;
        }
      }

      result = hil_set_analog_input_ranges
        (Experiment_2_moving_cart_DW.HILInitialize_Card,
         Experiment_2_moving_cart_P.HILInitialize_AIChannels, 8U,
         &Experiment_2_moving_cart_DW.HILInitialize_AIMinimums[0],
         &Experiment_2_moving_cart_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
        return;
      }
    }

    if ((Experiment_2_moving_cart_P.HILInitialize_AOPStart && !is_switching) ||
        (Experiment_2_moving_cart_P.HILInitialize_AOPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums =
          &Experiment_2_moving_cart_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = (Experiment_2_moving_cart_P.HILInitialize_AOLow);
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums =
          &Experiment_2_moving_cart_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = Experiment_2_moving_cart_P.HILInitialize_AOHigh;
        }
      }

      result = hil_set_analog_output_ranges
        (Experiment_2_moving_cart_DW.HILInitialize_Card,
         Experiment_2_moving_cart_P.HILInitialize_AOChannels, 8U,
         &Experiment_2_moving_cart_DW.HILInitialize_AOMinimums[0],
         &Experiment_2_moving_cart_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
        return;
      }
    }

    if ((Experiment_2_moving_cart_P.HILInitialize_AOStart && !is_switching) ||
        (Experiment_2_moving_cart_P.HILInitialize_AOEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages =
          &Experiment_2_moving_cart_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = Experiment_2_moving_cart_P.HILInitialize_AOInitial;
        }
      }

      result = hil_write_analog(Experiment_2_moving_cart_DW.HILInitialize_Card,
        Experiment_2_moving_cart_P.HILInitialize_AOChannels, 8U,
        &Experiment_2_moving_cart_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
        return;
      }
    }

    if (Experiment_2_moving_cart_P.HILInitialize_AOReset) {
      {
        int_T i1;
        real_T *dw_AOVoltages =
          &Experiment_2_moving_cart_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            Experiment_2_moving_cart_P.HILInitialize_AOWatchdog;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (Experiment_2_moving_cart_DW.HILInitialize_Card,
         Experiment_2_moving_cart_P.HILInitialize_AOChannels, 8U,
         &Experiment_2_moving_cart_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
        return;
      }
    }

    if ((Experiment_2_moving_cart_P.HILInitialize_EIPStart && !is_switching) ||
        (Experiment_2_moving_cart_P.HILInitialize_EIPEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &Experiment_2_moving_cart_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] =
            Experiment_2_moving_cart_P.HILInitialize_EIQuadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (Experiment_2_moving_cart_DW.HILInitialize_Card,
         Experiment_2_moving_cart_P.HILInitialize_EIChannels, 8U,
         (t_encoder_quadrature_mode *)
         &Experiment_2_moving_cart_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
        return;
      }
    }

    if ((Experiment_2_moving_cart_P.HILInitialize_EIStart && !is_switching) ||
        (Experiment_2_moving_cart_P.HILInitialize_EIEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &Experiment_2_moving_cart_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] =
            Experiment_2_moving_cart_P.HILInitialize_EIInitial;
        }
      }

      result = hil_set_encoder_counts
        (Experiment_2_moving_cart_DW.HILInitialize_Card,
         Experiment_2_moving_cart_P.HILInitialize_EIChannels, 8U,
         &Experiment_2_moving_cart_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
        return;
      }
    }

    if ((Experiment_2_moving_cart_P.HILInitialize_POPStart && !is_switching) ||
        (Experiment_2_moving_cart_P.HILInitialize_POPEnter && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &Experiment_2_moving_cart_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = Experiment_2_moving_cart_P.HILInitialize_POModes;
        }
      }

      result = hil_set_pwm_mode(Experiment_2_moving_cart_DW.HILInitialize_Card,
        Experiment_2_moving_cart_P.HILInitialize_POChannels, 8U, (t_pwm_mode *)
        &Experiment_2_moving_cart_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_POChannels =
          Experiment_2_moving_cart_P.HILInitialize_POChannels;
        int32_T *dw_POModeValues =
          &Experiment_2_moving_cart_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE ||
              dw_POModeValues[i1] == PWM_RAW_MODE) {
            Experiment_2_moving_cart_DW.HILInitialize_POSortedChans[num_duty_cycle_modes]
              = (p_HILInitialize_POChannels[i1]);
            Experiment_2_moving_cart_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]
              = Experiment_2_moving_cart_P.HILInitialize_POFrequency;
            num_duty_cycle_modes++;
          } else {
            Experiment_2_moving_cart_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = (p_HILInitialize_POChannels[i1]);
            Experiment_2_moving_cart_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] =
              Experiment_2_moving_cart_P.HILInitialize_POFrequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency
          (Experiment_2_moving_cart_DW.HILInitialize_Card,
           &Experiment_2_moving_cart_DW.HILInitialize_POSortedChans[0],
           num_duty_cycle_modes,
           &Experiment_2_moving_cart_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle
          (Experiment_2_moving_cart_DW.HILInitialize_Card,
           &Experiment_2_moving_cart_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
           num_frequency_modes,
           &Experiment_2_moving_cart_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &Experiment_2_moving_cart_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] =
            Experiment_2_moving_cart_P.HILInitialize_POConfiguration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &Experiment_2_moving_cart_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] =
            Experiment_2_moving_cart_P.HILInitialize_POAlignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &Experiment_2_moving_cart_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] =
            Experiment_2_moving_cart_P.HILInitialize_POPolarity;
        }
      }

      result = hil_set_pwm_configuration
        (Experiment_2_moving_cart_DW.HILInitialize_Card,
         Experiment_2_moving_cart_P.HILInitialize_POChannels, 8U,
         (t_pwm_configuration *)
         &Experiment_2_moving_cart_DW.HILInitialize_POModeValues[0],
         (t_pwm_alignment *)
         &Experiment_2_moving_cart_DW.HILInitialize_POAlignValues[0],
         (t_pwm_polarity *)
         &Experiment_2_moving_cart_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &Experiment_2_moving_cart_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] =
            Experiment_2_moving_cart_P.HILInitialize_POLeading;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues =
          &Experiment_2_moving_cart_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = Experiment_2_moving_cart_P.HILInitialize_POTrailing;
        }
      }

      result = hil_set_pwm_deadband
        (Experiment_2_moving_cart_DW.HILInitialize_Card,
         Experiment_2_moving_cart_P.HILInitialize_POChannels, 8U,
         &Experiment_2_moving_cart_DW.HILInitialize_POSortedFreqs[0],
         &Experiment_2_moving_cart_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
        return;
      }
    }

    if ((Experiment_2_moving_cart_P.HILInitialize_POStart && !is_switching) ||
        (Experiment_2_moving_cart_P.HILInitialize_POEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues =
          &Experiment_2_moving_cart_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = Experiment_2_moving_cart_P.HILInitialize_POInitial;
        }
      }

      result = hil_write_pwm(Experiment_2_moving_cart_DW.HILInitialize_Card,
        Experiment_2_moving_cart_P.HILInitialize_POChannels, 8U,
        &Experiment_2_moving_cart_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
        return;
      }
    }

    if (Experiment_2_moving_cart_P.HILInitialize_POReset) {
      {
        int_T i1;
        real_T *dw_POValues =
          &Experiment_2_moving_cart_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = Experiment_2_moving_cart_P.HILInitialize_POWatchdog;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (Experiment_2_moving_cart_DW.HILInitialize_Card,
         Experiment_2_moving_cart_P.HILInitialize_POChannels, 8U,
         &Experiment_2_moving_cart_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for TransportDelay: '<Root>/Transport Delay' */
  {
    real_T *pBuffer =
      &Experiment_2_moving_cart_DW.TransportDelay_RWORK.TUbufferArea[0];
    Experiment_2_moving_cart_DW.TransportDelay_IWORK.Tail = 0;
    Experiment_2_moving_cart_DW.TransportDelay_IWORK.Head = 0;
    Experiment_2_moving_cart_DW.TransportDelay_IWORK.Last = 0;
    Experiment_2_moving_cart_DW.TransportDelay_IWORK.CircularBufSize = 1024;
    pBuffer[0] = Experiment_2_moving_cart_P.TransportDelay_InitOutput;
    pBuffer[1024] = Experiment_2_moving_cart_M->Timing.t[0];
    Experiment_2_moving_cart_DW.TransportDelay_PWORK.TUbufferPtrs[0] = (void *)
      &pBuffer[0];
    Experiment_2_moving_cart_DW.TransportDelay_PWORK.TUbufferPtrs[1] = (void *)
      &pBuffer[1024];
  }

  /* SystemInitialize for Atomic SubSystem: '<S4>/Differentiator1' */
  /* SystemInitialize for MATLAB Function: '<S6>/Differentiator' */
  Experiment_2_moving_cart_DW.zp_not_empty = false;

  /* End of SystemInitialize for SubSystem: '<S4>/Differentiator1' */
  Experiment_2_moving_cart_DW.Phi_b.size[1] = 0;
  Experiment_2_moving_cart_DW.bD_e.size = 0;
  Experiment_2_moving_cart_DW.So_inv_lastCol_e.size = 0;
  Experiment_2_moving_cart_DW.zp_p.size = 0;

  /* SystemInitialize for Atomic SubSystem: '<S4>/Differentiator' */
  /* SystemInitialize for MATLAB Function: '<S5>/Differentiator' */
  Experiment_2_moving_cart_DW.zp_not_empty_n = false;

  /* End of SystemInitialize for SubSystem: '<S4>/Differentiator' */
  Experiment_2_moving_cart_DW.Phi.size[1] = 0;
  Experiment_2_moving_cart_DW.bD.size = 0;
  Experiment_2_moving_cart_DW.So_inv_lastCol.size = 0;
  Experiment_2_moving_cart_DW.zp.size = 0;
}

/* Model terminate function */
void Experiment_2_moving_cart_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: Experiment_2_moving_cart/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(Experiment_2_moving_cart_DW.HILInitialize_Card);
    hil_monitor_stop_all(Experiment_2_moving_cart_DW.HILInitialize_Card);
    is_switching = false;
    if ((Experiment_2_moving_cart_P.HILInitialize_AOTerminate && !is_switching) ||
        (Experiment_2_moving_cart_P.HILInitialize_AOExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages =
          &Experiment_2_moving_cart_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = Experiment_2_moving_cart_P.HILInitialize_AOFinal;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((Experiment_2_moving_cart_P.HILInitialize_POTerminate && !is_switching) ||
        (Experiment_2_moving_cart_P.HILInitialize_POExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues =
          &Experiment_2_moving_cart_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = Experiment_2_moving_cart_P.HILInitialize_POFinal;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(Experiment_2_moving_cart_DW.HILInitialize_Card
                         , Experiment_2_moving_cart_P.HILInitialize_AOChannels,
                         num_final_analog_outputs
                         , Experiment_2_moving_cart_P.HILInitialize_POChannels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         ,
                         &Experiment_2_moving_cart_DW.HILInitialize_AOVoltages[0]
                         , &Experiment_2_moving_cart_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog
            (Experiment_2_moving_cart_DW.HILInitialize_Card,
             Experiment_2_moving_cart_P.HILInitialize_AOChannels,
             num_final_analog_outputs,
             &Experiment_2_moving_cart_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm
            (Experiment_2_moving_cart_DW.HILInitialize_Card,
             Experiment_2_moving_cart_P.HILInitialize_POChannels,
             num_final_pwm_outputs,
             &Experiment_2_moving_cart_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(Experiment_2_moving_cart_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(Experiment_2_moving_cart_DW.HILInitialize_Card);
    hil_monitor_delete_all(Experiment_2_moving_cart_DW.HILInitialize_Card);
    hil_close(Experiment_2_moving_cart_DW.HILInitialize_Card);
    Experiment_2_moving_cart_DW.HILInitialize_Card = NULL;
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/
void MdlOutputs(int_T tid)
{
  Experiment_2_moving_cart_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  Experiment_2_moving_cart_update();
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
  Experiment_2_moving_cart_initialize();
}

void MdlTerminate(void)
{
  Experiment_2_moving_cart_terminate();
}

/* Registration function */
RT_MODEL_Experiment_2_moving__T *Experiment_2_moving_cart(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)Experiment_2_moving_cart_M, 0,
                sizeof(RT_MODEL_Experiment_2_moving__T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&Experiment_2_moving_cart_M->solverInfo,
                          &Experiment_2_moving_cart_M->Timing.simTimeStep);
    rtsiSetTPtr(&Experiment_2_moving_cart_M->solverInfo, &rtmGetTPtr
                (Experiment_2_moving_cart_M));
    rtsiSetStepSizePtr(&Experiment_2_moving_cart_M->solverInfo,
                       &Experiment_2_moving_cart_M->Timing.stepSize0);
    rtsiSetErrorStatusPtr(&Experiment_2_moving_cart_M->solverInfo,
                          (&rtmGetErrorStatus(Experiment_2_moving_cart_M)));
    rtsiSetRTModelPtr(&Experiment_2_moving_cart_M->solverInfo,
                      Experiment_2_moving_cart_M);
  }

  rtsiSetSimTimeStep(&Experiment_2_moving_cart_M->solverInfo, MAJOR_TIME_STEP);
  rtsiSetSolverName(&Experiment_2_moving_cart_M->solverInfo,"FixedStepDiscrete");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = Experiment_2_moving_cart_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    Experiment_2_moving_cart_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    Experiment_2_moving_cart_M->Timing.sampleTimes =
      (&Experiment_2_moving_cart_M->Timing.sampleTimesArray[0]);
    Experiment_2_moving_cart_M->Timing.offsetTimes =
      (&Experiment_2_moving_cart_M->Timing.offsetTimesArray[0]);

    /* task periods */
    Experiment_2_moving_cart_M->Timing.sampleTimes[0] = (0.0);
    Experiment_2_moving_cart_M->Timing.sampleTimes[1] = (0.01);

    /* task offsets */
    Experiment_2_moving_cart_M->Timing.offsetTimes[0] = (0.0);
    Experiment_2_moving_cart_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(Experiment_2_moving_cart_M,
             &Experiment_2_moving_cart_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = Experiment_2_moving_cart_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    Experiment_2_moving_cart_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(Experiment_2_moving_cart_M, 11.0);
  Experiment_2_moving_cart_M->Timing.stepSize0 = 0.01;
  Experiment_2_moving_cart_M->Timing.stepSize1 = 0.01;

  /* External mode info */
  Experiment_2_moving_cart_M->Sizes.checksums[0] = (1163847484U);
  Experiment_2_moving_cart_M->Sizes.checksums[1] = (2834647576U);
  Experiment_2_moving_cart_M->Sizes.checksums[2] = (2657793900U);
  Experiment_2_moving_cart_M->Sizes.checksums[3] = (2784945U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[8];
    Experiment_2_moving_cart_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = &rtAlwaysEnabled;
    systemRan[2] = &rtAlwaysEnabled;
    systemRan[3] = &rtAlwaysEnabled;
    systemRan[4] = &rtAlwaysEnabled;
    systemRan[5] = &rtAlwaysEnabled;
    systemRan[6] = &rtAlwaysEnabled;
    systemRan[7] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(Experiment_2_moving_cart_M->extModeInfo,
      &Experiment_2_moving_cart_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(Experiment_2_moving_cart_M->extModeInfo,
                        Experiment_2_moving_cart_M->Sizes.checksums);
    rteiSetTPtr(Experiment_2_moving_cart_M->extModeInfo, rtmGetTPtr
                (Experiment_2_moving_cart_M));
  }

  Experiment_2_moving_cart_M->solverInfoPtr =
    (&Experiment_2_moving_cart_M->solverInfo);
  Experiment_2_moving_cart_M->Timing.stepSize = (0.01);
  rtsiSetFixedStepSize(&Experiment_2_moving_cart_M->solverInfo, 0.01);
  rtsiSetSolverMode(&Experiment_2_moving_cart_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  Experiment_2_moving_cart_M->blockIO = ((void *) &Experiment_2_moving_cart_B);
  (void) memset(((void *) &Experiment_2_moving_cart_B), 0,
                sizeof(B_Experiment_2_moving_cart_T));

  /* parameters */
  Experiment_2_moving_cart_M->defaultParam = ((real_T *)
    &Experiment_2_moving_cart_P);

  /* states (dwork) */
  Experiment_2_moving_cart_M->dwork = ((void *) &Experiment_2_moving_cart_DW);
  (void) memset((void *)&Experiment_2_moving_cart_DW, 0,
                sizeof(DW_Experiment_2_moving_cart_T));

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    Experiment_2_moving_cart_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 19;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  Experiment_2_moving_cart_M->Sizes.numContStates = (0);/* Number of continuous states */
  Experiment_2_moving_cart_M->Sizes.numY = (0);/* Number of model outputs */
  Experiment_2_moving_cart_M->Sizes.numU = (0);/* Number of model inputs */
  Experiment_2_moving_cart_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  Experiment_2_moving_cart_M->Sizes.numSampTimes = (2);/* Number of sample times */
  Experiment_2_moving_cart_M->Sizes.numBlocks = (52);/* Number of blocks */
  Experiment_2_moving_cart_M->Sizes.numBlockIO = (12);/* Number of block outputs */
  Experiment_2_moving_cart_M->Sizes.numBlockPrms = (124);/* Sum of parameter "widths" */
  return Experiment_2_moving_cart_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
