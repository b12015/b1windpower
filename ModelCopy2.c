/*
 * File: ModelCopy2.c
 *
 * Code generated for Simulink model 'ModelCopy2'.
 *
 * Model version                  : 1.29
 * Simulink Coder version         : 8.7 (R2014b) 08-Sep-2014
 * C/C++ source code generated on : Wed Apr 15 14:56:58 2015
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: ARM Compatible->ARM 11
 * Code generation objectives:
 *    1. Execution efficiency
 *    2. RAM efficiency
 *    3. Traceability
 * Validation result: Passed (8), Warnings (5), Error (0)
 */

#include "ModelCopy2.h"
#include "ModelCopy2_private.h"

/* Block signals (auto storage) */
B_ModelCopy2_T ModelCopy2_B;

/* Continuous states */
X_ModelCopy2_T ModelCopy2_X;

/* Block states (auto storage) */
DW_ModelCopy2_T ModelCopy2_DW;

/* Real-time model */
RT_MODEL_ModelCopy2_T ModelCopy2_M_;
RT_MODEL_ModelCopy2_T *const ModelCopy2_M = &ModelCopy2_M_;

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

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si , real_T arg_VIND)
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = (ODE3_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 25;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  ModelCopy2_derivatives(arg_VIND);

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  ModelCopy2_custom(arg_VIND);
  ModelCopy2_derivatives(arg_VIND);

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  ModelCopy2_custom(arg_VIND);
  ModelCopy2_derivatives(arg_VIND);

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
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
        y = (rtNaN);
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

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  int32_T u0_0;
  int32_T u1_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      u0_0 = 1;
    } else {
      u0_0 = -1;
    }

    if (u1 > 0.0) {
      u1_0 = 1;
    } else {
      u1_0 = -1;
    }

    y = atan2(u0_0, u1_0);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

/* Model step function */
real_T ModelCopy2_custom(real_T arg_VIND)
{
  /* local block i/o variables */
  real_T rtb_Clock;
  real_T rtb_LogicalOperator;
  real_T rtb_integrator;
  real_T rtb_TransportDelay;
  real_T rtb_integrator_m;
  real_T rtb_TransportDelay_a;
  real_T rtb_integrator_c;
  real_T rtb_TransportDelay_e;
  real_T rtb_integrator_mf;
  real_T rtb_TransportDelay_p;
  real_T rtb_integrator_b;
  real_T rtb_TransportDelay_n;
  real_T rtb_integrator_mk;
  real_T rtb_TransportDelay_m;
  real_T rtb_Integrator[2];
  real_T rtb_Integrator_g[2];
  real_T rtb_Sum7_c;
  real_T rtb_Sum7_a;
  real_T rtb_Sum7_p;
  real_T rtb_Sum7_i;
  real_T rtb_Sum_co;
  boolean_T rtb_Compare;
  real_T rtb_Llr2;
  real_T rtb_ids;
  real_T rtb_Llr1_o;
  real_T rtb_iqs;
  real_T rtb_Product2_dm;
  boolean_T rtb_LogicalOperator_c;
  real_T rtb_Powerspeed;
  real_T rtb_ComplextoMagnitudeAngle_o2;
  real_T rtb_iqr;
  real_T rtb_idr;
  real_T rtb_Vdc;
  real_T rtb_Switch1;
  real_T rtb_Switch2;
  real_T rtb_Switch3;
  real_T rtb_Switch1_c;
  real_T rtb_Ton;
  real_T rtb_Switch3_a;
  real_T rtb_wwr;
  real_T rtb_ComplextoRealImag1_o2;
  real_T rtb_Product4;
  real_T rtb_Product5;
  real_T rtb_Product6;
  real_T rtb_Sum_b;
  real_T rtb_ComplextoRealImag_o2;
  real_T rtb_Sum7_l;
  real_T rtb_Switch5_idx_0;
  real_T rtb_Switch1_l_idx_0;
  real_T rtb_Switch1_l_idx_1;
  real_T rtb_Sum_l_idx_0;
  real_T rtb_Sum_l_idx_1;
  real_T rtb_Sum_l_idx_2;
  real_T re;
  real_T rtb_Sum1_re;
  real_T rtb_Sum1_im;
  real_T rtb_donotdeletethisgain_re;
  real_T rtb_donotdeletethisgain_im;
  real_T rtb_a23_re;
  real_T rtb_a23_im;
  real_T rtb_Kv1_idx_0_re;
  real_T rtb_Kv1_idx_0_im;
  real_T rtb_Kv1_idx_1_re;
  real_T rtb_Kv1_idx_1_im;
  real_T rtb_Kv1_idx_2_re;
  real_T rtb_Kv1_idx_2_im;
  real_T im;

  /* specified return value */
  real_T arg_PUT;
  if (rtmIsMajorTimeStep(ModelCopy2_M)) {
    /* set solver stop time */
    rtsiSetSolverStopTime(&ModelCopy2_M->solverInfo,
                          ((ModelCopy2_M->Timing.clockTick0+1)*
      ModelCopy2_M->Timing.stepSize0));
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(ModelCopy2_M)) {
    ModelCopy2_M->Timing.t[0] = rtsiGetT(&ModelCopy2_M->solverInfo);
  }

  /* RelationalOperator: '<S3>/Compare' incorporates:
   *  Constant: '<S3>/Constant'
   *  Inport: '<Root>/VIND'
   */
  rtb_Compare = (arg_VIND >= ModelCopy2_P.CompareToConstant_const);
  if (rtmIsMajorTimeStep(ModelCopy2_M)) {
    /* RelationalOperator: '<S43>/Relational Operator1' incorporates:
     *  Constant: '<S43>/Constant'
     *  Constant: '<S43>/Constant2'
     */
    ModelCopy2_B.RelationalOperator1 = (ModelCopy2_P.Bistable1_Qpriority >
      ModelCopy2_P.Constant_Value_j);

    /* Outputs for Enabled SubSystem: '<S43>/SET  Priority' incorporates:
     *  EnablePort: '<S45>/Enable'
     */
    if (rtmIsMajorTimeStep(ModelCopy2_M)) {
      if (ModelCopy2_B.RelationalOperator1) {
        if (!ModelCopy2_DW.SETPriority_MODE) {
          ModelCopy2_DW.SETPriority_MODE = true;
        }
      } else {
        if (ModelCopy2_DW.SETPriority_MODE) {
          ModelCopy2_DW.SETPriority_MODE = false;
        }
      }
    }

    /* End of Outputs for SubSystem: '<S43>/SET  Priority' */
  }

  /* Outputs for Enabled SubSystem: '<S43>/SET  Priority' incorporates:
   *  EnablePort: '<S45>/Enable'
   */
  if (ModelCopy2_DW.SETPriority_MODE) {
    if (rtmIsMajorTimeStep(ModelCopy2_M)) {
      /* Memory: '<S45>/IC=ic' */
      ModelCopy2_B.ICic_c = ModelCopy2_DW.ICic_PreviousInput_n;
    }

    /* Logic: '<S45>/Logical Operator' incorporates:
     *  Logic: '<S45>/Logical Operator2'
     */
    ModelCopy2_B.LogicalOperator_e = (rtb_Compare || ModelCopy2_B.ICic_c);
  }

  /* End of Outputs for SubSystem: '<S43>/SET  Priority' */

  /* Outputs for Enabled SubSystem: '<S43>/RESET Priority' incorporates:
   *  EnablePort: '<S44>/Enable'
   */
  if (rtmIsMajorTimeStep(ModelCopy2_M) && rtmIsMajorTimeStep(ModelCopy2_M)) {
    /* Logic: '<S43>/Logical Operator2' */
    if (!ModelCopy2_B.RelationalOperator1) {
      if (!ModelCopy2_DW.RESETPriority_MODE) {
        ModelCopy2_DW.RESETPriority_MODE = true;
      }
    } else {
      if (ModelCopy2_DW.RESETPriority_MODE) {
        ModelCopy2_DW.RESETPriority_MODE = false;
      }
    }

    /* End of Logic: '<S43>/Logical Operator2' */
  }

  if (ModelCopy2_DW.RESETPriority_MODE) {
    if (rtmIsMajorTimeStep(ModelCopy2_M)) {
      /* Memory: '<S44>/IC=ic' */
      ModelCopy2_B.ICic_cd = ModelCopy2_DW.ICic_PreviousInput_c;
    }

    /* Logic: '<S44>/Logical Operator2' incorporates:
     *  Logic: '<S44>/Logical Operator'
     */
    ModelCopy2_B.LogicalOperator2_i = (rtb_Compare || ModelCopy2_B.ICic_cd);
  }

  /* End of Outputs for SubSystem: '<S43>/RESET Priority' */

  /* Logic: '<S43>/Logical Operator1' */
  rtb_Compare = (ModelCopy2_B.LogicalOperator_e ||
                 ModelCopy2_B.LogicalOperator2_i);

  /* Switch: '<S41>/Switch1' incorporates:
   *  Constant: '<S41>/Constant'
   *  Constant: '<S41>/Constant1'
   *  Constant: '<S42>/Constant'
   *  Constant: '<S42>/Constant2'
   *  Integrator: '<S41>/phidr'
   *  Integrator: '<S41>/phiqr'
   *  Integrator: '<S42>/phids'
   *  Integrator: '<S42>/phiqs'
   *  Switch: '<S41>/Switch2'
   *  Switch: '<S42>/Switch1'
   *  Switch: '<S42>/Switch3'
   */
  if (rtb_Compare) {
    rtb_Switch1 = ModelCopy2_P.Constant_Value_e;
    rtb_Switch2 = ModelCopy2_P.Constant1_Value_k;
    rtb_Switch3 = ModelCopy2_P.Constant2_Value_h;
    rtb_Switch1_c = ModelCopy2_P.Constant_Value_o;
  } else {
    rtb_Switch1 = ModelCopy2_X.phiqr_CSTATE;
    rtb_Switch2 = ModelCopy2_X.phidr_CSTATE;
    rtb_Switch3 = ModelCopy2_X.phiqs_CSTATE;
    rtb_Switch1_c = ModelCopy2_X.phids_CSTATE;
  }

  /* End of Switch: '<S41>/Switch1' */

  /* Gain: '<S40>/1\Llr2' incorporates:
   *  Gain: '<S40>/1\Llr'
   *  Gain: '<S40>/1\Lls'
   *  Sum: '<S40>/Sum1'
   */
  rtb_Llr2 = (ModelCopy2_P.Llr_Gain * rtb_Switch2 + ModelCopy2_P.Lls_Gain *
              rtb_Switch1_c) * ModelCopy2_P.Llr2_Gain;

  /* Gain: '<S42>/1\Llr2' incorporates:
   *  Sum: '<S42>/Sum3'
   */
  rtb_ids = (rtb_Switch1_c - rtb_Llr2) * ModelCopy2_P.Llr2_Gain_b;

  /* Gain: '<S40>/1\Llr1' incorporates:
   *  Gain: '<S40>/1\Llr'
   *  Gain: '<S40>/1\Lls'
   *  Sum: '<S40>/Sum'
   */
  rtb_Llr1_o = (ModelCopy2_P.Llr_Gain * rtb_Switch1 + ModelCopy2_P.Lls_Gain *
                rtb_Switch3) * ModelCopy2_P.Llr1_Gain_c;

  /* Gain: '<S42>/1\Llr' incorporates:
   *  Sum: '<S42>/Sum4'
   */
  rtb_iqs = (rtb_Switch3 - rtb_Llr1_o) * ModelCopy2_P.Llr_Gain_i;

  /* Switch: '<S89>/Switch1' incorporates:
   *  Constant: '<S89>/Constant'
   *  Integrator: '<S89>/Integrator'
   */
  if (rtb_Compare) {
    rtb_Switch1_l_idx_0 = ModelCopy2_P.Constant_Value_k;
    rtb_Switch1_l_idx_1 = ModelCopy2_P.Constant_Value_k;
  } else {
    rtb_Switch1_l_idx_0 = ModelCopy2_X.Integrator_CSTATE[0];
    rtb_Switch1_l_idx_1 = ModelCopy2_X.Integrator_CSTATE[1];
  }

  /* End of Switch: '<S89>/Switch1' */

  /* Sum: '<S16>/Sum' incorporates:
   *  RealImagToComplex: '<S36>/Real-Imag to Complex'
   *  RealImagToComplex: '<S90>/Real-Imag to Complex1'
   */
  rtb_a23_re = rtb_ids + rtb_Switch1_l_idx_0;
  rtb_a23_im = rtb_iqs + rtb_Switch1_l_idx_1;

  /* Sum: '<S16>/Sum1' incorporates:
   *  Gain: '<S36>/a^2'
   *  Gain: '<S90>/a^2'
   *  RealImagToComplex: '<S36>/Real-Imag to Complex'
   *  RealImagToComplex: '<S90>/Real-Imag to Complex1'
   */
  rtb_Sum1_re = (ModelCopy2_P.a2_Gain.re * rtb_ids - ModelCopy2_P.a2_Gain.im *
                 rtb_iqs) + (ModelCopy2_P.a2_Gain_b.re * rtb_Switch1_l_idx_0 -
    ModelCopy2_P.a2_Gain_b.im * rtb_Switch1_l_idx_1);
  rtb_Sum1_im = (ModelCopy2_P.a2_Gain.re * rtb_iqs + ModelCopy2_P.a2_Gain.im *
                 rtb_ids) + (ModelCopy2_P.a2_Gain_b.re * rtb_Switch1_l_idx_1 +
    ModelCopy2_P.a2_Gain_b.im * rtb_Switch1_l_idx_0);

  /* Clock: '<S9>/Clock' */
  rtb_Clock = ModelCopy2_M->Timing.t[0];
  if (rtmIsMajorTimeStep(ModelCopy2_M)) {
    /* RelationalOperator: '<S6>/Relational Operator' incorporates:
     *  Constant: '<S6>/Constant'
     *  Constant: '<S6>/valp_nom3'
     */
    ModelCopy2_B.Amplitude = (ModelCopy2_P.ukV_VariationEntity ==
      ModelCopy2_P.Constant_Value_n);

    /* Logic: '<S6>/Logical Operator1' incorporates:
     *  Constant: '<S6>/Constant6'
     *  Constant: '<S6>/valp_nom5'
     *  RelationalOperator: '<S6>/Relational Operator3'
     */
    ModelCopy2_B.LogicalOperator1 = ((ModelCopy2_P.valp_nom5_Value ==
      ModelCopy2_P.Constant6_Value_j) && ModelCopy2_B.Amplitude);
  }

  /* Step: '<S10>/Step1' */
  if (ModelCopy2_M->Timing.t[0] < ModelCopy2_P.VariationSubSystem_Toff_Variati)
  {
    rtb_Product2_dm = ModelCopy2_P.Step1_Y0;
  } else {
    rtb_Product2_dm = ModelCopy2_P.Step1_YFinal;
  }

  /* End of Step: '<S10>/Step1' */
  if (rtmIsMajorTimeStep(ModelCopy2_M)) {
    /* RelationalOperator: '<S10>/Relational Operator1' incorporates:
     *  Constant: '<S10>/Constant3'
     *  Constant: '<S6>/valp_nom5'
     */
    ModelCopy2_B.RelationalOperator1_b = (ModelCopy2_P.valp_nom5_Value ==
      ModelCopy2_P.Constant3_Value_d);
  }

  /* Step: '<S10>/Step' */
  if (ModelCopy2_M->Timing.t[0] < ModelCopy2_P.VariationSubSystem_Ton_Variatio)
  {
    rtb_Ton = ModelCopy2_P.Step_Y0_l;
  } else {
    rtb_Ton = ModelCopy2_P.Step_YFinal_g;
  }

  /* End of Step: '<S10>/Step' */
  if (rtmIsMajorTimeStep(ModelCopy2_M)) {
    /* Gain: '<S10>/Gain1' incorporates:
     *  Constant: '<S10>/valp_nom9'
     */
    ModelCopy2_B.Gain1 = ModelCopy2_P.Gain1_Gain *
      ModelCopy2_P.VariationSubSystem_VariationFre;

    /* Memory: '<S10>/Memory' */
    ModelCopy2_B.Memory = ModelCopy2_DW.Memory_PreviousInput;
  }

  /* Switch: '<S10>/Switch2' */
  if (rtb_Product2_dm >= ModelCopy2_P.Switch2_Threshold) {
    /* MultiPortSwitch: '<S10>/Multiport Switch1' incorporates:
     *  Constant: '<S10>/Constant5'
     *  Constant: '<S10>/valp_nom6'
     *  Constant: '<S10>/valp_nom8'
     *  Constant: '<S6>/valp_nom5'
     *  Integrator: '<S10>/Integrator'
     *  Product: '<S10>/Product'
     *  Product: '<S10>/Product1'
     *  Product: '<S10>/Product2'
     *  Trigonometry: '<S10>/Trigonometric Function1'
     */
    switch ((int32_T)ModelCopy2_P.valp_nom5_Value) {
     case 1:
      ModelCopy2_B.Switch2 = ModelCopy2_P.VariationSubSystem_VariationSte *
        rtb_Ton;
      break;

     case 2:
      ModelCopy2_B.Switch2 = ModelCopy2_X.Integrator_CSTATE_j;
      break;

     case 3:
      ModelCopy2_B.Switch2 = sin(ModelCopy2_X.Integrator_CSTATE_j *
        ModelCopy2_B.Gain1) * ModelCopy2_P.VariationSubSystem_VariationMag;
      break;

     default:
      ModelCopy2_B.Switch2 = ModelCopy2_P.Constant5_Value_d;
      break;
    }

    /* End of MultiPortSwitch: '<S10>/Multiport Switch1' */
  } else {
    ModelCopy2_B.Switch2 = ModelCopy2_B.Memory;
  }

  /* End of Switch: '<S10>/Switch2' */

  /* Switch: '<S10>/Switch3' incorporates:
   *  Constant: '<S10>/Constant1'
   *  DataTypeConversion: '<S10>/Data Type  Conversion2'
   *  Logic: '<S10>/Logical Operator'
   *  Logic: '<S10>/Logical Operator1'
   */
  if ((!(rtb_Product2_dm != 0.0)) && ModelCopy2_B.RelationalOperator1_b) {
    rtb_Switch3_a = ModelCopy2_P.Constant1_Value_m;
  } else {
    rtb_Switch3_a = ModelCopy2_B.Switch2;
  }

  /* End of Switch: '<S10>/Switch3' */

  /* Switch: '<S6>/Switch2' incorporates:
   *  Constant: '<S6>/Constant1'
   */
  if (ModelCopy2_B.Amplitude) {
    rtb_Sum_l_idx_1 = rtb_Switch3_a;
  } else {
    rtb_Sum_l_idx_1 = ModelCopy2_P.Constant1_Value_n;
  }

  /* End of Switch: '<S6>/Switch2' */

  /* Sum: '<S6>/Sum3' incorporates:
   *  Constant: '<S6>/valp_nom2'
   */
  rtb_Product2_dm = rtb_Sum_l_idx_1 + ModelCopy2_P.valp_nom2_Value;

  /* Switch: '<S6>/Switch1' incorporates:
   *  Lookup: '<S9>/Look-Up Table'
   */
  if (ModelCopy2_B.LogicalOperator1) {
    rtb_Powerspeed = rt_Lookup(ModelCopy2_P.LookUpTable_XData, 7, rtb_Clock,
      ModelCopy2_P.LookUpTable_YData);
  } else {
    rtb_Powerspeed = rtb_Product2_dm;
  }

  /* End of Switch: '<S6>/Switch1' */

  /* Switch: '<S6>/Switch5' incorporates:
   *  Constant: '<S6>/SinglePhase'
   */
  rtb_Switch5_idx_0 = rtb_Powerspeed;
  if (!(ModelCopy2_P.SinglePhase_Value >= ModelCopy2_P.Switch5_Threshold)) {
    rtb_Product2_dm = rtb_Powerspeed;
  }

  /* End of Switch: '<S6>/Switch5' */
  if (rtmIsMajorTimeStep(ModelCopy2_M)) {
    /* Gain: '<S6>/Gain3' incorporates:
     *  Constant: '<S6>/valp_nom'
     */
    ModelCopy2_B.Gain3 = ModelCopy2_P.Gain3_Gain * ModelCopy2_P.valp_nom_Value;

    /* RelationalOperator: '<S6>/Relational Operator1' incorporates:
     *  Constant: '<S6>/Constant2'
     *  Constant: '<S6>/valp_nom3'
     */
    ModelCopy2_B.Phase = (ModelCopy2_P.ukV_VariationEntity ==
                          ModelCopy2_P.Constant2_Value_l);
  }

  /* Switch: '<S6>/Switch3' incorporates:
   *  Constant: '<S6>/Constant4'
   *  Gain: '<S6>/Gain4'
   */
  if (ModelCopy2_B.Phase) {
    rtb_Powerspeed = ModelCopy2_P.Gain4_Gain * rtb_Switch3_a;
  } else {
    rtb_Powerspeed = ModelCopy2_P.Constant4_Value_j;
  }

  /* End of Switch: '<S6>/Switch3' */
  if (rtmIsMajorTimeStep(ModelCopy2_M)) {
    /* DataTypeConversion: '<S6>/Data Type  Conversion2' incorporates:
     *  Constant: '<S6>/valp_nom7'
     */
    ModelCopy2_B.DataTypeConversion2 = (ModelCopy2_P.ukV_HarmonicGeneration !=
      0.0);
  }

  /* Step: '<S6>/Step' */
  if (ModelCopy2_M->Timing.t[0] < ModelCopy2_P.Step_Time) {
    rtb_Sum_l_idx_1 = ModelCopy2_P.Step_Y0_f;
  } else {
    rtb_Sum_l_idx_1 = ModelCopy2_P.Step_YFinal_ih;
  }

  /* End of Step: '<S6>/Step' */

  /* Step: '<S6>/Step1' */
  if (ModelCopy2_M->Timing.t[0] < ModelCopy2_P.Step1_Time) {
    rtb_Sum_l_idx_2 = ModelCopy2_P.Step1_Y0_l;
  } else {
    rtb_Sum_l_idx_2 = ModelCopy2_P.Step1_YFinal_p;
  }

  /* End of Step: '<S6>/Step1' */

  /* Logic: '<S6>/Logical Operator' incorporates:
   *  DataTypeConversion: '<S6>/Data Type  Conversion1'
   *  Sum: '<S6>/Sum4'
   */
  rtb_LogicalOperator_c = ((rtb_Sum_l_idx_1 + rtb_Sum_l_idx_2 != 0.0) &&
    ModelCopy2_B.DataTypeConversion2);
  if (rtmIsMajorTimeStep(ModelCopy2_M)) {
    /* Gain: '<S7>/Gain3' incorporates:
     *  Constant: '<S7>/Phase_Harmo'
     */
    rtb_Switch3_a = ModelCopy2_P.Gain3_Gain_c *
      ModelCopy2_P.SeqAGeneration_Phase_Harmo;

    /* MultiPortSwitch: '<S7>/Multiport Switch' incorporates:
     *  Constant: '<S7>/Negative-sequence'
     *  Constant: '<S7>/Phase_Harmo2'
     *  Constant: '<S7>/Positive-sequence'
     *  Constant: '<S7>/Zero-sequence'
     *  Constant: '<S7>/valp_nom2'
     *  Sum: '<S7>/Sum1'
     */
    switch ((int32_T)(ModelCopy2_P.SeqAGeneration_Seq_Harmo +
                      ModelCopy2_P.valp_nom2_Value_d)) {
     case 1:
      rtb_Sum_l_idx_0 = ModelCopy2_P.Zerosequence_Value[0];
      rtb_Sum_l_idx_1 = ModelCopy2_P.Zerosequence_Value[1];
      rtb_Sum_l_idx_2 = ModelCopy2_P.Zerosequence_Value[2];
      break;

     case 2:
      rtb_Sum_l_idx_0 = ModelCopy2_P.Positivesequence_Value[0];
      rtb_Sum_l_idx_1 = ModelCopy2_P.Positivesequence_Value[1];
      rtb_Sum_l_idx_2 = ModelCopy2_P.Positivesequence_Value[2];
      break;

     default:
      rtb_Sum_l_idx_0 = ModelCopy2_P.Negativesequence_Value[0];
      rtb_Sum_l_idx_1 = ModelCopy2_P.Negativesequence_Value[1];
      rtb_Sum_l_idx_2 = ModelCopy2_P.Negativesequence_Value[2];
      break;
    }

    /* End of MultiPortSwitch: '<S7>/Multiport Switch' */

    /* Sum: '<S7>/Sum' */
    rtb_Sum_l_idx_0 += rtb_Switch3_a;
    rtb_Sum_l_idx_1 += rtb_Switch3_a;
    rtb_Sum_l_idx_2 += rtb_Switch3_a;

    /* MagnitudeAngleToComplex: '<S7>/Magnitude-Angle to Complex' incorporates:
     *  Constant: '<S7>/Phase_Harmo1'
     */
    ModelCopy2_B.MagnitudeAngletoComplex[0].re =
      ModelCopy2_P.SeqAGeneration_Mag_Harmo * cos(rtb_Sum_l_idx_0);
    ModelCopy2_B.MagnitudeAngletoComplex[0].im =
      ModelCopy2_P.SeqAGeneration_Mag_Harmo * sin(rtb_Sum_l_idx_0);
    ModelCopy2_B.MagnitudeAngletoComplex[1].re =
      ModelCopy2_P.SeqAGeneration_Mag_Harmo * cos(rtb_Sum_l_idx_1);
    ModelCopy2_B.MagnitudeAngletoComplex[1].im =
      ModelCopy2_P.SeqAGeneration_Mag_Harmo * sin(rtb_Sum_l_idx_1);
    ModelCopy2_B.MagnitudeAngletoComplex[2].re =
      ModelCopy2_P.SeqAGeneration_Mag_Harmo * cos(rtb_Sum_l_idx_2);
    ModelCopy2_B.MagnitudeAngletoComplex[2].im =
      ModelCopy2_P.SeqAGeneration_Mag_Harmo * sin(rtb_Sum_l_idx_2);

    /* Gain: '<S8>/Gain3' incorporates:
     *  Constant: '<S8>/Phase_Harmo'
     */
    rtb_Switch3_a = ModelCopy2_P.Gain3_Gain_m *
      ModelCopy2_P.SeqBGeneration_Phase_Harmo;

    /* MultiPortSwitch: '<S8>/Multiport Switch' incorporates:
     *  Constant: '<S8>/Negative-sequence'
     *  Constant: '<S8>/Phase_Harmo2'
     *  Constant: '<S8>/Positive-sequence'
     *  Constant: '<S8>/Zero-sequence'
     *  Constant: '<S8>/valp_nom2'
     *  Sum: '<S8>/Sum1'
     */
    switch ((int32_T)(ModelCopy2_P.SeqBGeneration_Seq_Harmo +
                      ModelCopy2_P.valp_nom2_Value_j)) {
     case 1:
      rtb_Sum_l_idx_0 = ModelCopy2_P.Zerosequence_Value_m[0];
      rtb_Sum_l_idx_1 = ModelCopy2_P.Zerosequence_Value_m[1];
      rtb_Sum_l_idx_2 = ModelCopy2_P.Zerosequence_Value_m[2];
      break;

     case 2:
      rtb_Sum_l_idx_0 = ModelCopy2_P.Positivesequence_Value_h[0];
      rtb_Sum_l_idx_1 = ModelCopy2_P.Positivesequence_Value_h[1];
      rtb_Sum_l_idx_2 = ModelCopy2_P.Positivesequence_Value_h[2];
      break;

     default:
      rtb_Sum_l_idx_0 = ModelCopy2_P.Negativesequence_Value_a[0];
      rtb_Sum_l_idx_1 = ModelCopy2_P.Negativesequence_Value_a[1];
      rtb_Sum_l_idx_2 = ModelCopy2_P.Negativesequence_Value_a[2];
      break;
    }

    /* End of MultiPortSwitch: '<S8>/Multiport Switch' */

    /* Sum: '<S8>/Sum' */
    rtb_Sum_l_idx_0 += rtb_Switch3_a;
    rtb_Sum_l_idx_1 += rtb_Switch3_a;
    rtb_Sum_l_idx_2 += rtb_Switch3_a;

    /* MagnitudeAngleToComplex: '<S8>/Magnitude-Angle to Complex' incorporates:
     *  Constant: '<S8>/Phase_Harmo1'
     */
    ModelCopy2_B.MagnitudeAngletoComplex_l[0].re =
      ModelCopy2_P.SeqBGeneration_Mag_Harmo * cos(rtb_Sum_l_idx_0);
    ModelCopy2_B.MagnitudeAngletoComplex_l[0].im =
      ModelCopy2_P.SeqBGeneration_Mag_Harmo * sin(rtb_Sum_l_idx_0);
    ModelCopy2_B.MagnitudeAngletoComplex_l[1].re =
      ModelCopy2_P.SeqBGeneration_Mag_Harmo * cos(rtb_Sum_l_idx_1);
    ModelCopy2_B.MagnitudeAngletoComplex_l[1].im =
      ModelCopy2_P.SeqBGeneration_Mag_Harmo * sin(rtb_Sum_l_idx_1);
    ModelCopy2_B.MagnitudeAngletoComplex_l[2].re =
      ModelCopy2_P.SeqBGeneration_Mag_Harmo * cos(rtb_Sum_l_idx_2);
    ModelCopy2_B.MagnitudeAngletoComplex_l[2].im =
      ModelCopy2_P.SeqBGeneration_Mag_Harmo * sin(rtb_Sum_l_idx_2);
  }

  /* Sum: '<S6>/Sum2' incorporates:
   *  Constant: '<S6>/P1'
   */
  rtb_Sum_l_idx_1 = (ModelCopy2_B.Gain3 + ModelCopy2_P.P1_Value[0]) +
    rtb_Powerspeed;
  rtb_Sum_l_idx_0 = (ModelCopy2_B.Gain3 + ModelCopy2_P.P1_Value[1]) +
    rtb_Powerspeed;
  rtb_Sum_l_idx_2 = (ModelCopy2_B.Gain3 + ModelCopy2_P.P1_Value[2]) +
    rtb_Powerspeed;

  /* ComplexToRealImag: '<S96>/Complex to Real-Imag' incorporates:
   *  DataTypeConversion: '<S6>/Data Type  Conversion'
   *  Gain: '<S4>/pu->A '
   *  Gain: '<S4>/pu->A  '
   *  MagnitudeAngleToComplex: '<S6>/Magnitude-Angle to Complex'
   *  Product: '<S7>/Product1'
   *  Product: '<S8>/Product1'
   *  Sum: '<S6>/Sum5'
   */
  ModelCopy2_B.ComplextoRealImag_o1[0] = ModelCopy2_P.puA_Gain * rtb_a23_re;
  ModelCopy2_B.ComplextoRealImag_o1[1] = ModelCopy2_P.puA_Gain_o * rtb_Sum1_re;
  ModelCopy2_B.ComplextoRealImag_o1[2] = (rtb_Switch5_idx_0 * cos
    (rtb_Sum_l_idx_1) + (real_T)rtb_LogicalOperator_c *
    ModelCopy2_B.MagnitudeAngletoComplex[0].re) + (real_T)rtb_LogicalOperator_c *
    ModelCopy2_B.MagnitudeAngletoComplex_l[0].re;
  ModelCopy2_B.ComplextoRealImag_o1[3] = (rtb_Product2_dm * cos(rtb_Sum_l_idx_0)
    + (real_T)rtb_LogicalOperator_c * ModelCopy2_B.MagnitudeAngletoComplex[1].re)
    + (real_T)rtb_LogicalOperator_c * ModelCopy2_B.MagnitudeAngletoComplex_l[1].
    re;
  ModelCopy2_B.ComplextoRealImag_o1[4] = (rtb_Product2_dm * cos(rtb_Sum_l_idx_2)
    + (real_T)rtb_LogicalOperator_c * ModelCopy2_B.MagnitudeAngletoComplex[2].re)
    + (real_T)rtb_LogicalOperator_c * ModelCopy2_B.MagnitudeAngletoComplex_l[2].
    re;
  ModelCopy2_B.ComplextoRealImag_o2[0] = ModelCopy2_P.puA_Gain * rtb_a23_im;
  ModelCopy2_B.ComplextoRealImag_o2[1] = ModelCopy2_P.puA_Gain_o * rtb_Sum1_im;
  ModelCopy2_B.ComplextoRealImag_o2[2] = (rtb_Switch5_idx_0 * sin
    (rtb_Sum_l_idx_1) + (real_T)rtb_LogicalOperator_c *
    ModelCopy2_B.MagnitudeAngletoComplex[0].im) + (real_T)rtb_LogicalOperator_c *
    ModelCopy2_B.MagnitudeAngletoComplex_l[0].im;
  ModelCopy2_B.ComplextoRealImag_o2[3] = (rtb_Product2_dm * sin(rtb_Sum_l_idx_0)
    + (real_T)rtb_LogicalOperator_c * ModelCopy2_B.MagnitudeAngletoComplex[1].im)
    + (real_T)rtb_LogicalOperator_c * ModelCopy2_B.MagnitudeAngletoComplex_l[1].
    im;
  ModelCopy2_B.ComplextoRealImag_o2[4] = (rtb_Product2_dm * sin(rtb_Sum_l_idx_2)
    + (real_T)rtb_LogicalOperator_c * ModelCopy2_B.MagnitudeAngletoComplex[2].im)
    + (real_T)rtb_LogicalOperator_c * ModelCopy2_B.MagnitudeAngletoComplex_l[2].
    im;

  /* Level2 S-Function Block: '<S92>/State-Space' (sfun_psbdqc) */
  {
    SimStruct *rts = ModelCopy2_M->childSfunctions[0];
    sfcnOutputs(rts, 0);
  }

  /* Gain: '<S15>/Kv1' incorporates:
   *  Gain: '<S23>/do not delete this gain'
   *  Gain: '<S24>/do not delete this gain'
   *  Gain: '<S25>/do not delete this gain'
   *  RealImagToComplex: '<S95>/Real-Imag to Complex'
   */
  rtb_Kv1_idx_0_re = ModelCopy2_P.donotdeletethisgain_Gain *
    ModelCopy2_B.StateSpace[0] * ModelCopy2_P.Kv1_Gain;
  rtb_Kv1_idx_0_im = ModelCopy2_P.donotdeletethisgain_Gain *
    ModelCopy2_B.StateSpace[3] * ModelCopy2_P.Kv1_Gain;
  rtb_Kv1_idx_1_re = ModelCopy2_P.donotdeletethisgain_Gain_e *
    ModelCopy2_B.StateSpace[1] * ModelCopy2_P.Kv1_Gain;
  rtb_Kv1_idx_1_im = ModelCopy2_P.donotdeletethisgain_Gain_e *
    ModelCopy2_B.StateSpace[4] * ModelCopy2_P.Kv1_Gain;
  rtb_Kv1_idx_2_re = ModelCopy2_P.donotdeletethisgain_Gain_f *
    ModelCopy2_B.StateSpace[2] * ModelCopy2_P.Kv1_Gain;
  rtb_Kv1_idx_2_im = ModelCopy2_P.donotdeletethisgain_Gain_f *
    ModelCopy2_B.StateSpace[5] * ModelCopy2_P.Kv1_Gain;

  /* Gain: '<S33>/pu->V' */
  rtb_Sum_l_idx_1 = ModelCopy2_P.puV_Gain * rtb_Kv1_idx_0_re;
  rtb_Switch5_idx_0 = ModelCopy2_P.puV_Gain * rtb_Kv1_idx_0_im;

  /* Math: '<S88>/Math Function' incorporates:
   *  Gain: '<S33>/pu->A'
   */
  rtb_Product2_dm = ModelCopy2_P.puA_Gain_i * rtb_a23_re;
  im = -(ModelCopy2_P.puA_Gain_i * rtb_a23_im);

  /* Gain: '<S33>/pu->V' */
  re = ModelCopy2_P.puV_Gain * rtb_Kv1_idx_1_re;
  rtb_wwr = ModelCopy2_P.puV_Gain * rtb_Kv1_idx_1_im;

  /* Math: '<S88>/Math Function' incorporates:
   *  Gain: '<S33>/pu->A'
   */
  rtb_ComplextoRealImag_o2 = ModelCopy2_P.puA_Gain_i * rtb_Sum1_re;
  rtb_Product6 = -(ModelCopy2_P.puA_Gain_i * rtb_Sum1_im);

  /* Gain: '<S33>/pu->V' */
  rtb_donotdeletethisgain_re = ModelCopy2_P.puV_Gain * rtb_Kv1_idx_2_re;
  rtb_Sum_b = ModelCopy2_P.puV_Gain * rtb_Kv1_idx_2_im;

  /* Math: '<S88>/Math Function' incorporates:
   *  Gain: '<S33>/-1 '
   *  Gain: '<S33>/pu->A'
   *  Sum: '<S33>/Sum2'
   */
  rtb_Sum_l_idx_0 = (rtb_a23_re + rtb_Sum1_re) * ModelCopy2_P.u_Gain *
    ModelCopy2_P.puA_Gain_i;
  rtb_Sum_l_idx_2 = -((rtb_a23_im + rtb_Sum1_im) * ModelCopy2_P.u_Gain *
                      ModelCopy2_P.puA_Gain_i);

  /* Gain: '<S33>/-1' incorporates:
   *  Gain: '<S33>/-->pu'
   *  Gain: '<S88>/K'
   *  Product: '<S88>/Product2'
   *  Sum: '<S88>/Sum3'
   */
  ModelCopy2_B.Ppu = (((rtb_Sum_l_idx_1 * rtb_Product2_dm - rtb_Switch5_idx_0 *
                        im) + (re * rtb_ComplextoRealImag_o2 - rtb_wwr *
    rtb_Product6)) + (rtb_donotdeletethisgain_re * rtb_Sum_l_idx_0 - rtb_Sum_b *
                      rtb_Sum_l_idx_2)) * ModelCopy2_P.K_Gain *
    ModelCopy2_P.pu_Gain * ModelCopy2_P.u_Gain_f;

  /* Outport: '<Root>/PUT' incorporates:
   *  Gain: '<Root>/900k'
   */
  arg_PUT = ModelCopy2_P.u0k_Gain * ModelCopy2_B.Ppu;

  /* Gain: '<S33>/-2' incorporates:
   *  Gain: '<S33>/-->pu1'
   *  Gain: '<S88>/K'
   *  Product: '<S88>/Product2'
   *  Sum: '<S88>/Sum3'
   */
  ModelCopy2_B.Qpu = (((rtb_Sum_l_idx_1 * im + rtb_Switch5_idx_0 *
                        rtb_Product2_dm) + (re * rtb_Product6 + rtb_wwr *
    rtb_ComplextoRealImag_o2)) + (rtb_donotdeletethisgain_re * rtb_Sum_l_idx_2 +
    rtb_Sum_b * rtb_Sum_l_idx_0)) * ModelCopy2_P.K_Gain * ModelCopy2_P.pu1_Gain *
    ModelCopy2_P.u_Gain_b;

  /* Gain: '<Root>/900k_ ' */
  ModelCopy2_B.QUT = ModelCopy2_P.u0k__Gain * ModelCopy2_B.Qpu;

  /* Saturate: '<Root>/Saturation' incorporates:
   *  Inport: '<Root>/VIND'
   */
  if (arg_VIND > ModelCopy2_P.Saturation_UpperSat_c) {
    rtb_Powerspeed = ModelCopy2_P.Saturation_UpperSat_c;
  } else if (arg_VIND < ModelCopy2_P.Saturation_LowerSat_l) {
    rtb_Powerspeed = ModelCopy2_P.Saturation_LowerSat_l;
  } else {
    rtb_Powerspeed = arg_VIND;
  }

  /* End of Saturate: '<Root>/Saturation' */

  /* Gain: '<S19>/1//wind_base' */
  rtb_Product2_dm = 1.0 / ModelCopy2_P.WindTurbineDoublyFedInductio_of *
    rtb_Powerspeed;

  /* Fcn: '<S19>/wind_speed^3' */
  rtb_ComplextoRealImag1_o2 = rt_powd_snf(rtb_Product2_dm, 3.0);

  /* Integrator: '<S35>/Integrator' */
  rtb_Sum_l_idx_0 = ModelCopy2_X.Integrator_CSTATE_c;

  /* Saturate: '<S19>/Avoid division by zero ' */
  if (rtb_Product2_dm > ModelCopy2_P.Avoiddivisionbyzero_UpperSat_o) {
    rtb_Product2_dm = ModelCopy2_P.Avoiddivisionbyzero_UpperSat_o;
  } else {
    if (rtb_Product2_dm < ModelCopy2_P.Avoiddivisionbyzero_LowerSat_g) {
      rtb_Product2_dm = ModelCopy2_P.Avoiddivisionbyzero_LowerSat_g;
    }
  }

  /* Gain: '<S19>/lambda_nom' incorporates:
   *  Gain: '<S19>/pu->pu '
   *  Integrator: '<S35>/Integrator'
   *  Product: '<S19>/Product'
   *  Saturate: '<S19>/Avoid division by zero '
   */
  rtb_Product2_dm = 1.0 / ModelCopy2_P.WindTurbine_speed_nom *
    ModelCopy2_X.Integrator_CSTATE_c * (1.0 / rtb_Product2_dm) *
    ModelCopy2_P.lambda_nom_Gain;

  /* Saturate: '<S19>/Saturation1' */
  if (rtb_Product2_dm > ModelCopy2_P.Saturation1_UpperSat) {
    rtb_Product2_dm = ModelCopy2_P.Saturation1_UpperSat;
  } else {
    if (rtb_Product2_dm < ModelCopy2_P.Saturation1_LowerSat) {
      rtb_Product2_dm = ModelCopy2_P.Saturation1_LowerSat;
    }
  }

  /* End of Saturate: '<S19>/Saturation1' */

  /* Gain: '<S31>/pitch_gain' incorporates:
   *  Constant: '<S31>/Constant2'
   *  Integrator: '<S35>/Integrator'
   *  Sum: '<S31>/Sum'
   */
  rtb_Switch3_a = (ModelCopy2_X.Integrator_CSTATE_c -
                   ModelCopy2_P.Constant2_Value_m) *
    ModelCopy2_P.WindTurbineDoublyFedInduction_n;

  /* Saturate: '<S31>/0-pitch_max' */
  if (rtb_Switch3_a > ModelCopy2_P.WindTurbineDoublyFedInduction_m) {
    rtb_Switch3_a = ModelCopy2_P.WindTurbineDoublyFedInduction_m;
  } else {
    if (rtb_Switch3_a < ModelCopy2_P.pitch_max_LowerSat) {
      rtb_Switch3_a = ModelCopy2_P.pitch_max_LowerSat;
    }
  }

  /* End of Saturate: '<S31>/0-pitch_max' */

  /* RateLimiter: '<S31>/Rate Limiter   ' */
  if (ModelCopy2_DW.LastMajorTime == (rtInf)) {
    ModelCopy2_B.RateLimiter = rtb_Switch3_a;
  } else {
    im = ModelCopy2_M->Timing.t[0] - ModelCopy2_DW.LastMajorTime;
    re = im * ModelCopy2_P.WindTurbineDoublyFedInduction_o;
    rtb_Sum_l_idx_2 = rtb_Switch3_a - ModelCopy2_DW.PrevY;
    if (rtb_Sum_l_idx_2 > re) {
      ModelCopy2_B.RateLimiter = ModelCopy2_DW.PrevY + re;
    } else {
      im *= -ModelCopy2_P.WindTurbineDoublyFedInduction_o;
      if (rtb_Sum_l_idx_2 < im) {
        ModelCopy2_B.RateLimiter = ModelCopy2_DW.PrevY + im;
      } else {
        ModelCopy2_B.RateLimiter = rtb_Switch3_a;
      }
    }
  }

  /* End of RateLimiter: '<S31>/Rate Limiter   ' */

  /* Fcn: '<S91>/Fcn' */
  rtb_ComplextoMagnitudeAngle_o2 = 1.0 / (1.0 / (0.08 * ModelCopy2_B.RateLimiter
    + rtb_Product2_dm) - 0.035 / (rt_powd_snf(ModelCopy2_B.RateLimiter, 3.0) +
    1.0));

  /* Fcn: '<S91>/Fcn1' */
  rtb_Product2_dm = ((116.0 / rtb_ComplextoMagnitudeAngle_o2 - 0.4 *
                      ModelCopy2_B.RateLimiter) - 5.0) * 0.51763 * exp(-21.0 /
    rtb_ComplextoMagnitudeAngle_o2) + 0.006795 * rtb_Product2_dm;

  /* Switch: '<S4>/Switch' incorporates:
   *  Constant: '<S4>/ExternalTorque'
   *  Fcn: '<S19>/wind_speed^3'
   *  Gain: '<S19>/1//cp_nom'
   *  Gain: '<S19>/Gain'
   *  Gain: '<S19>/pu->pu'
   *  Product: '<S19>/Product '
   *  Product: '<S19>/Product2'
   *  Saturate: '<S19>/Avoid division by zero'
   */
  if (ModelCopy2_P.ExternalTorque_Value_p >= ModelCopy2_P.Switch_Threshold_c) {
    /* Saturate: '<S19>/Avoid division by zero' incorporates:
     *  Integrator: '<S35>/Integrator'
     */
    if (ModelCopy2_X.Integrator_CSTATE_c >
        ModelCopy2_P.Avoiddivisionbyzero_UpperSat) {
      rtb_Sum_l_idx_1 = ModelCopy2_P.Avoiddivisionbyzero_UpperSat;
    } else if (ModelCopy2_X.Integrator_CSTATE_c <
               ModelCopy2_P.Avoiddivisionbyzero_LowerSat) {
      rtb_Sum_l_idx_1 = ModelCopy2_P.Avoiddivisionbyzero_LowerSat;
    } else {
      rtb_Sum_l_idx_1 = ModelCopy2_X.Integrator_CSTATE_c;
    }

    rtb_Powerspeed = ModelCopy2_P.WindTurbineDoublyFedInduction_k *
      ModelCopy2_P.WindTurbineDoublyFedInductionGe /
      ModelCopy2_P.WindTurbine_Pelec_base * (ModelCopy2_P.cp_nom_Gain *
      rtb_Product2_dm * rtb_ComplextoRealImag1_o2) / rtb_Sum_l_idx_1 *
      ModelCopy2_P.Gain_Gain;
  }

  /* End of Switch: '<S4>/Switch' */

  /* Gain: '<S35>/1_2H' incorporates:
   *  Gain: '<S35>/F'
   *  Integrator: '<S35>/Integrator'
   *  Product: '<S35>/Product1'
   *  Product: '<S35>/Product2'
   *  Sum: '<S35>/Sum'
   *  Sum: '<S35>/Sum1'
   */
  ModelCopy2_B._2H = (((rtb_Switch1_c * rtb_iqs - rtb_Switch3 * rtb_ids) -
                       rtb_Powerspeed) - ModelCopy2_P.F_Gain *
                      ModelCopy2_X.Integrator_CSTATE_c) * ModelCopy2_P._2H_Gain;

  /* Gain: '<S35>/web' incorporates:
   *  Integrator: '<S35>/Integrator'
   */
  ModelCopy2_B.web = ModelCopy2_P.web_Gain_p * ModelCopy2_X.Integrator_CSTATE_c;

  /* Gain: '<S41>/1\Llr' incorporates:
   *  Sum: '<S41>/Sum3'
   */
  rtb_iqr = (rtb_Switch1 - rtb_Llr1_o) * ModelCopy2_P.Llr_Gain_if;

  /* Gain: '<S41>/1\Llr2' incorporates:
   *  Sum: '<S41>/Sum4'
   */
  rtb_idr = (rtb_Switch2 - rtb_Llr2) * ModelCopy2_P.Llr2_Gain_f;

  /* Sum: '<S29>/Sum1' incorporates:
   *  Constant: '<S29>/ws'
   *  Integrator: '<S35>/Integrator'
   */
  rtb_Sum1_im = ModelCopy2_P.ws_Value - ModelCopy2_X.Integrator_CSTATE_c;

  /* Integrator: '<S32>/Integrator' */
  rtb_Vdc = ModelCopy2_X.Integrator_CSTATE_cb;

  /* Sum: '<S47>/Sum3' incorporates:
   *  Gain: '<S47>/(a^2)//3'
   *  Gain: '<S47>/Gain'
   *  Gain: '<S47>/a//3'
   */
  rtb_donotdeletethisgain_re = ((ModelCopy2_P.a3_Gain.re * rtb_Kv1_idx_1_re -
    ModelCopy2_P.a3_Gain.im * rtb_Kv1_idx_1_im) + ModelCopy2_P.Gain_Gain_c *
    rtb_Kv1_idx_0_re) + (ModelCopy2_P.a23_Gain.re * rtb_Kv1_idx_2_re -
    ModelCopy2_P.a23_Gain.im * rtb_Kv1_idx_2_im);
  rtb_donotdeletethisgain_im = ((ModelCopy2_P.a3_Gain.re * rtb_Kv1_idx_1_im +
    ModelCopy2_P.a3_Gain.im * rtb_Kv1_idx_1_re) + ModelCopy2_P.Gain_Gain_c *
    rtb_Kv1_idx_0_im) + (ModelCopy2_P.a23_Gain.re * rtb_Kv1_idx_2_im +
    ModelCopy2_P.a23_Gain.im * rtb_Kv1_idx_2_re);

  /* ComplexToMagnitudeAngle: '<S47>/Complex to Magnitude-Angle' */
  ModelCopy2_B.ComplextoMagnitudeAngle_o1 = rt_hypotd_snf
    (rtb_donotdeletethisgain_re, rtb_donotdeletethisgain_im);
  rtb_ComplextoMagnitudeAngle_o2 = rt_atan2d_snf(rtb_donotdeletethisgain_im,
    rtb_donotdeletethisgain_re);

  /* Relay: '<S47>/Relay' */
  if (rtmIsMajorTimeStep(ModelCopy2_M)) {
    if (ModelCopy2_B.ComplextoMagnitudeAngle_o1 >= ModelCopy2_P.Relay_OnVal_o) {
      ModelCopy2_DW.Relay_Mode = true;
    } else {
      if (ModelCopy2_B.ComplextoMagnitudeAngle_o1 <= ModelCopy2_P.Relay_OffVal_c)
      {
        ModelCopy2_DW.Relay_Mode = false;
      }
    }
  }

  if (rtmIsMajorTimeStep(ModelCopy2_M)) {
    /* Memory: '<S47>/IC=ic' */
    ModelCopy2_B.ICic = ModelCopy2_DW.ICic_PreviousInput;
  }

  /* Relay: '<S47>/Relay' */
  if (ModelCopy2_DW.Relay_Mode) {
    rtb_Sum_l_idx_1 = ModelCopy2_P.Relay_YOn_b;
  } else {
    rtb_Sum_l_idx_1 = ModelCopy2_P.Relay_YOff_p;
  }

  /* Switch: '<S47>/20%' */
  if (rtb_Sum_l_idx_1 >= ModelCopy2_P.u_Threshold_e) {
    ModelCopy2_B.u = rtb_ComplextoMagnitudeAngle_o2;
  } else {
    ModelCopy2_B.u = ModelCopy2_B.ICic;
  }

  /* End of Switch: '<S47>/20%' */

  /* Lookup: '<S31>/Power(speed)' incorporates:
   *  Integrator: '<S35>/Integrator'
   */
  rtb_Switch3_a = rt_Lookup(ModelCopy2_P.Powerspeed_XData, 6,
    ModelCopy2_X.Integrator_CSTATE_c, ModelCopy2_P.Powerspeed_YData);

  /* Saturate: '<S31>/0-inf' */
  if (rtb_Switch3_a > ModelCopy2_P.inf_UpperSat_d) {
    rtb_Switch3_a = ModelCopy2_P.inf_UpperSat_d;
  } else {
    if (rtb_Switch3_a < ModelCopy2_P.inf_LowerSat_b) {
      rtb_Switch3_a = ModelCopy2_P.inf_LowerSat_b;
    }
  }

  /* End of Saturate: '<S31>/0-inf' */

  /* RateLimiter: '<S31>/Rate Limiter ' */
  if (ModelCopy2_DW.LastMajorTime_i == (rtInf)) {
    ModelCopy2_B.RateLimiter_b = rtb_Switch3_a;
  } else {
    im = ModelCopy2_M->Timing.t[0] - ModelCopy2_DW.LastMajorTime_i;
    re = im * ModelCopy2_P.RateLimiter_RisingLim_f;
    rtb_Sum_l_idx_2 = rtb_Switch3_a - ModelCopy2_DW.PrevY_k;
    if (rtb_Sum_l_idx_2 > re) {
      ModelCopy2_B.RateLimiter_b = ModelCopy2_DW.PrevY_k + re;
    } else {
      im *= ModelCopy2_P.RateLimiter_FallingLim_l;
      if (rtb_Sum_l_idx_2 < im) {
        ModelCopy2_B.RateLimiter_b = ModelCopy2_DW.PrevY_k + im;
      } else {
        ModelCopy2_B.RateLimiter_b = rtb_Switch3_a;
      }
    }
  }

  /* End of RateLimiter: '<S31>/Rate Limiter ' */

  /* Logic: '<S31>/Logical Operator' */
  rtb_LogicalOperator = !rtb_Compare;

  /* Outputs for Enabled SubSystem: '<S31>/wind_dfig_rotor' incorporates:
   *  EnablePort: '<S49>/Enable'
   */
  if (rtmIsMajorTimeStep(ModelCopy2_M) && rtmIsMajorTimeStep(ModelCopy2_M)) {
    if (rtb_LogicalOperator > 0.0) {
      if (!ModelCopy2_DW.wind_dfig_rotor_MODE) {
        ModelCopy2_DW.wind_dfig_rotor_MODE = true;
      }
    } else {
      if (ModelCopy2_DW.wind_dfig_rotor_MODE) {
        /* Disable for Enabled SubSystem: '<S49>/V Regulator' */
        if (ModelCopy2_DW.VRegulator_MODE) {
          ModelCopy2_DW.VRegulator_MODE = false;
        }

        /* End of Disable for SubSystem: '<S49>/V Regulator' */

        /* Disable for Enabled SubSystem: '<S49>/Q Regulator' */
        if (ModelCopy2_DW.QRegulator_MODE) {
          ModelCopy2_DW.QRegulator_MODE = false;
        }

        /* End of Disable for SubSystem: '<S49>/Q Regulator' */

        /* Disable for Enabled SubSystem: '<S64>/Subsystem' */
        if (ModelCopy2_DW.Subsystem_MODE) {
          ModelCopy2_DW.Subsystem_MODE = false;
        }

        /* End of Disable for SubSystem: '<S64>/Subsystem' */

        /* Disable for Enabled SubSystem: '<S64>/Subsystem ' */
        if (ModelCopy2_DW.Subsystem_MODE_o) {
          ModelCopy2_DW.Subsystem_MODE_o = false;
        }

        /* End of Disable for SubSystem: '<S64>/Subsystem ' */
        ModelCopy2_DW.wind_dfig_rotor_MODE = false;
      }
    }
  }

  if (ModelCopy2_DW.wind_dfig_rotor_MODE) {
    /* Sum: '<S60>/Sum1' incorporates:
     *  Constant: '<S60>/w=1pu'
     *  Integrator: '<S35>/Integrator'
     */
    rtb_wwr = ModelCopy2_P.w1pu_Value - ModelCopy2_X.Integrator_CSTATE_c;

    /* Sum: '<S61>/Sum' incorporates:
     *  ComplexToMagnitudeAngle: '<S61>/Complex to Magnitude-Angle'
     *  RealImagToComplex: '<S61>/Real-Imag to Complex'
     */
    rtb_Switch3_a = rt_atan2d_snf(rtb_iqs, rtb_ids) - ModelCopy2_B.u;

    /* MagnitudeAngleToComplex: '<S61>/Magnitude-Angle to Complex' incorporates:
     *  ComplexToMagnitudeAngle: '<S61>/Complex to Magnitude-Angle'
     *  RealImagToComplex: '<S61>/Real-Imag to Complex'
     */
    rtb_donotdeletethisgain_re = rt_hypotd_snf(rtb_ids, rtb_iqs) * cos
      (rtb_Switch3_a);
    rtb_donotdeletethisgain_im = rt_hypotd_snf(rtb_ids, rtb_iqs) * sin
      (rtb_Switch3_a);

    /* ComplexToRealImag: '<S61>/Complex to Real-Imag' */
    rtb_Switch3_a = rtb_donotdeletethisgain_re;
    rtb_ComplextoRealImag_o2 = rtb_donotdeletethisgain_im;

    /* ComplexToMagnitudeAngle: '<S62>/Complex to Magnitude-Angle1' */
    rtb_Product2_dm = rt_hypotd_snf(rtb_donotdeletethisgain_re,
      rtb_donotdeletethisgain_im);
    rtb_Llr1_o = rt_atan2d_snf(rtb_donotdeletethisgain_im,
      rtb_donotdeletethisgain_re);

    /* Sum: '<S61>/Sum1' incorporates:
     *  ComplexToMagnitudeAngle: '<S61>/Complex to Magnitude-Angle1'
     *  RealImagToComplex: '<S61>/Real-Imag to Complex1'
     */
    rtb_Sum_b = rt_atan2d_snf(rtb_iqr, rtb_idr) - ModelCopy2_B.u;

    /* MagnitudeAngleToComplex: '<S61>/Magnitude-Angle to Complex1' incorporates:
     *  ComplexToMagnitudeAngle: '<S61>/Complex to Magnitude-Angle1'
     *  RealImagToComplex: '<S61>/Real-Imag to Complex1'
     */
    rtb_donotdeletethisgain_re = rt_hypotd_snf(rtb_idr, rtb_iqr) * cos(rtb_Sum_b);
    rtb_donotdeletethisgain_im = rt_hypotd_snf(rtb_idr, rtb_iqr) * sin(rtb_Sum_b);

    /* ComplexToRealImag: '<S61>/Complex to Real-Imag1' */
    rtb_Sum_b = rtb_donotdeletethisgain_re;
    rtb_Sum7_l = rtb_donotdeletethisgain_im;

    /* Product: '<S63>/Product1' incorporates:
     *  ComplexToRealImag: '<S61>/Complex to Real-Imag1'
     *  Constant: '<S63>/Constant1'
     *  Sum: '<S63>/Sum1'
     */
    rtb_Powerspeed = (rtb_Switch3_a + rtb_donotdeletethisgain_re) *
      ModelCopy2_P.Constant1_Value_h;

    /* Product: '<S63>/Product5' incorporates:
     *  ComplexToRealImag: '<S61>/Complex to Real-Imag1'
     *  Constant: '<S63>/Constant1'
     *  Sum: '<S63>/Sum2'
     */
    rtb_Llr2 = (rtb_ComplextoRealImag_o2 + rtb_donotdeletethisgain_im) *
      ModelCopy2_P.Constant1_Value_h;

    /* Fcn: '<S70>/x->r' */
    rtb_Switch3_a = rt_hypotd_snf(rtb_Powerspeed, rtb_Llr2);

    /* Relay: '<S63>/Relay' */
    if (rtmIsMajorTimeStep(ModelCopy2_M)) {
      if (rtb_Switch3_a >= ModelCopy2_P.Relay_OnVal) {
        ModelCopy2_DW.Relay_Mode_p = true;
      } else {
        if (rtb_Switch3_a <= ModelCopy2_P.Relay_OffVal) {
          ModelCopy2_DW.Relay_Mode_p = false;
        }
      }
    }

    if (rtmIsMajorTimeStep(ModelCopy2_M)) {
      /* Memory: '<S63>/IC=ic' */
      ModelCopy2_B.ICic_a = ModelCopy2_DW.ICic_PreviousInput_b;
    }

    /* Relay: '<S63>/Relay' */
    if (ModelCopy2_DW.Relay_Mode_p) {
      rtb_Sum_l_idx_1 = ModelCopy2_P.Relay_YOn;
    } else {
      rtb_Sum_l_idx_1 = ModelCopy2_P.Relay_YOff;
    }

    /* Switch: '<S63>/20%' incorporates:
     *  Fcn: '<S70>/x->theta'
     */
    if (rtb_Sum_l_idx_1 >= ModelCopy2_P.u_Threshold) {
      ModelCopy2_B.u_g = rt_atan2d_snf(rtb_Llr2, rtb_Powerspeed);
    } else {
      ModelCopy2_B.u_g = ModelCopy2_B.ICic_a;
    }

    /* End of Switch: '<S63>/20%' */

    /* Sum: '<S62>/Sum1' */
    rtb_Llr1_o -= ModelCopy2_B.u_g;

    /* MagnitudeAngleToComplex: '<S62>/Magnitude-Angle to Complex1' */
    rtb_donotdeletethisgain_im = rtb_Product2_dm * sin(rtb_Llr1_o);

    /* ComplexToRealImag: '<S62>/Complex to Real-Imag1' incorporates:
     *  MagnitudeAngleToComplex: '<S62>/Magnitude-Angle to Complex1'
     */
    rtb_Llr2 = rtb_Product2_dm * cos(rtb_Llr1_o);
    rtb_ComplextoRealImag1_o2 = rtb_donotdeletethisgain_im;

    /* Product: '<S60>/Product1' incorporates:
     *  ComplexToRealImag: '<S62>/Complex to Real-Imag1'
     *  Constant: '<S60>/Lm1'
     */
    rtb_a23_re = rtb_wwr * ModelCopy2_P.Lm1_Value * rtb_donotdeletethisgain_im;

    /* RealImagToComplex: '<S62>/Real-Imag to Complex' */
    rtb_donotdeletethisgain_im = rtb_Sum7_l;

    /* Sum: '<S62>/Sum' incorporates:
     *  ComplexToMagnitudeAngle: '<S62>/Complex to Magnitude-Angle'
     *  ComplexToRealImag: '<S61>/Complex to Real-Imag1'
     *  RealImagToComplex: '<S62>/Real-Imag to Complex'
     */
    rtb_Sum7_l = rt_atan2d_snf(rtb_Sum7_l, rtb_donotdeletethisgain_re) -
      ModelCopy2_B.u_g;

    /* MagnitudeAngleToComplex: '<S62>/Magnitude-Angle to Complex' incorporates:
     *  ComplexToMagnitudeAngle: '<S62>/Complex to Magnitude-Angle'
     *  ComplexToRealImag: '<S61>/Complex to Real-Imag1'
     *  RealImagToComplex: '<S62>/Real-Imag to Complex'
     */
    rtb_donotdeletethisgain_re = rt_hypotd_snf(rtb_donotdeletethisgain_re,
      rtb_donotdeletethisgain_im) * cos(rtb_Sum7_l);
    rtb_donotdeletethisgain_im = rt_hypotd_snf(rtb_Sum_b,
      rtb_donotdeletethisgain_im) * sin(rtb_Sum7_l);

    /* ComplexToRealImag: '<S62>/Complex to Real-Imag' */
    ModelCopy2_B.ComplextoRealImag_o1_b = rtb_donotdeletethisgain_re;
    ModelCopy2_B.ComplextoRealImag_o2_p = rtb_donotdeletethisgain_im;

    /* Product: '<S60>/Product2' incorporates:
     *  Constant: '<S60>/Llr+Lm1'
     */
    rtb_Sum1_re = rtb_wwr * ModelCopy2_P.LlrLm1_Value *
      ModelCopy2_B.ComplextoRealImag_o2_p;

    /* Product: '<S60>/Product3' incorporates:
     *  Constant: '<S60>/Llr+Lm2'
     */
    rtb_a23_im = rtb_wwr * ModelCopy2_P.LlrLm2_Value *
      ModelCopy2_B.ComplextoRealImag_o1_b;

    /* Product: '<S60>/Product4' incorporates:
     *  Constant: '<S60>/Lm2'
     */
    rtb_Product4 = rtb_wwr * ModelCopy2_P.Lm2_Value * rtb_Llr2;

    /* Product: '<S60>/Product5' incorporates:
     *  Constant: '<S60>/Lm3'
     */
    rtb_Product5 = ModelCopy2_P.Lm3_Value * ModelCopy2_B.ComplextoRealImag_o2_p;

    /* Product: '<S60>/Product6' incorporates:
     *  Constant: '<S60>/Lm4'
     */
    rtb_Product6 = ModelCopy2_P.Lm4_Value * ModelCopy2_B.ComplextoRealImag_o1_b;

    /* Sum: '<S61>/Sum2' incorporates:
     *  ComplexToMagnitudeAngle: '<S61>/Complex to Magnitude-Angle2'
     *  RealImagToComplex: '<S61>/Real-Imag to Complex2'
     */
    rtb_Sum7_l = rt_atan2d_snf(rtb_Switch1_l_idx_1, rtb_Switch1_l_idx_0) -
      ModelCopy2_B.u;

    /* MagnitudeAngleToComplex: '<S61>/Magnitude-Angle to Complex2' incorporates:
     *  ComplexToMagnitudeAngle: '<S61>/Complex to Magnitude-Angle2'
     *  RealImagToComplex: '<S61>/Real-Imag to Complex2'
     */
    rtb_donotdeletethisgain_re = rt_hypotd_snf(rtb_Switch1_l_idx_0,
      rtb_Switch1_l_idx_1) * cos(rtb_Sum7_l);
    rtb_donotdeletethisgain_im = rt_hypotd_snf(rtb_Switch1_l_idx_0,
      rtb_Switch1_l_idx_1) * sin(rtb_Sum7_l);

    /* Outputs for Enabled SubSystem: '<S49>/V Regulator' incorporates:
     *  EnablePort: '<S67>/Enable'
     */
    if (rtmIsMajorTimeStep(ModelCopy2_M) && rtmIsMajorTimeStep(ModelCopy2_M)) {
      /* Constant: '<S49>/ExternalTorque1' */
      if (ModelCopy2_P.ExternalTorque1_Value > 0.0) {
        if (!ModelCopy2_DW.VRegulator_MODE) {
          ModelCopy2_DW.VRegulator_MODE = true;
        }
      } else {
        if (ModelCopy2_DW.VRegulator_MODE) {
          ModelCopy2_DW.VRegulator_MODE = false;
        }
      }
    }

    if (ModelCopy2_DW.VRegulator_MODE) {
      /* Gain: '<S67>/Droop' incorporates:
       *  ComplexToRealImag: '<S61>/Complex to Real-Imag2'
       *  Sum: '<S67>/Sum1'
       */
      ModelCopy2_B.Droop = (rtb_ComplextoRealImag_o2 +
                            rtb_donotdeletethisgain_im) *
        ModelCopy2_P.Droop_Gain;

      /* Integrator: '<S84>/integrator' */
      rtb_integrator = ModelCopy2_X.integrator_CSTATE;

      /* TransportDelay: '<S84>/Transport Delay' */
      {
        real_T **uBuffer = (real_T**)
          &ModelCopy2_DW.TransportDelay_PWORK.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)
          &ModelCopy2_DW.TransportDelay_PWORK.TUbufferPtrs[1];
        real_T simTime = ModelCopy2_M->Timing.t[0];
        real_T tMinusDelay = simTime - ModelCopy2_P.MeanValue_Period_d;
        rtb_TransportDelay = rt_TDelayInterpolate(
          tMinusDelay,
          0.0,
          *tBuffer,
          *uBuffer,
          ModelCopy2_DW.TransportDelay_IWORK.CircularBufSize,
          &ModelCopy2_DW.TransportDelay_IWORK.Last,
          ModelCopy2_DW.TransportDelay_IWORK.Tail,
          ModelCopy2_DW.TransportDelay_IWORK.Head,
          ModelCopy2_P.TransportDelay_InitOutput_i,
          0,
          0);
      }

      /* Step: '<S84>/Step' */
      if (ModelCopy2_M->Timing.t[0] < ModelCopy2_P.MeanValue_Period_d +
          2.2204460492503131E-16) {
        rtb_Sum7_p = ModelCopy2_P.Step_Y0_o;
      } else {
        rtb_Sum7_p = ModelCopy2_P.Step_YFinal_l;
      }

      /* End of Step: '<S84>/Step' */

      /* RateLimiter: '<S67>/Rate Limiter ' incorporates:
       *  Constant: '<S4>/Vref '
       */
      if (ModelCopy2_DW.LastMajorTime_kl == (rtInf)) {
        ModelCopy2_B.RateLimiter_b5 = ModelCopy2_P.Vref_Value;
      } else {
        im = ModelCopy2_M->Timing.t[0] - ModelCopy2_DW.LastMajorTime_kl;
        re = im * ModelCopy2_P.RateLimiter_RisingLim_a;
        rtb_Sum_l_idx_2 = ModelCopy2_P.Vref_Value - ModelCopy2_DW.PrevY_i;
        if (rtb_Sum_l_idx_2 > re) {
          ModelCopy2_B.RateLimiter_b5 = ModelCopy2_DW.PrevY_i + re;
        } else {
          im *= ModelCopy2_P.RateLimiter_FallingLim_o;
          if (rtb_Sum_l_idx_2 < im) {
            ModelCopy2_B.RateLimiter_b5 = ModelCopy2_DW.PrevY_i + im;
          } else {
            ModelCopy2_B.RateLimiter_b5 = ModelCopy2_P.Vref_Value;
          }
        }
      }

      /* End of RateLimiter: '<S67>/Rate Limiter ' */

      /* Switch: '<S84>/Switch' incorporates:
       *  Gain: '<S84>/Gain'
       *  Sum: '<S84>/Sum'
       */
      if (rtb_Sum7_p >= ModelCopy2_P.Switch_Threshold_a) {
        rtb_Sum7_p = 1.0 / ModelCopy2_P.MeanValue_Period_d * (rtb_integrator -
          rtb_TransportDelay);
      } else {
        rtb_Sum7_p = ModelCopy2_B.RateLimiter_b5;
      }

      /* End of Switch: '<S84>/Switch' */

      /* Integrator: '<S85>/integrator' */
      rtb_integrator_m = ModelCopy2_X.integrator_CSTATE_p;

      /* TransportDelay: '<S85>/Transport Delay' */
      {
        real_T **uBuffer = (real_T**)
          &ModelCopy2_DW.TransportDelay_PWORK_g.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)
          &ModelCopy2_DW.TransportDelay_PWORK_g.TUbufferPtrs[1];
        real_T simTime = ModelCopy2_M->Timing.t[0];
        real_T tMinusDelay = simTime - ModelCopy2_P.MeanValue1_Period;
        rtb_TransportDelay_a = rt_TDelayInterpolate(
          tMinusDelay,
          0.0,
          *tBuffer,
          *uBuffer,
          ModelCopy2_DW.TransportDelay_IWORK_g.CircularBufSize,
          &ModelCopy2_DW.TransportDelay_IWORK_g.Last,
          ModelCopy2_DW.TransportDelay_IWORK_g.Tail,
          ModelCopy2_DW.TransportDelay_IWORK_g.Head,
          ModelCopy2_P.TransportDelay_InitOutput_a,
          0,
          0);
      }

      /* Step: '<S85>/Step' */
      if (ModelCopy2_M->Timing.t[0] < ModelCopy2_P.MeanValue1_Period +
          2.2204460492503131E-16) {
        rtb_Sum_l_idx_1 = ModelCopy2_P.Step_Y0_a;
      } else {
        rtb_Sum_l_idx_1 = ModelCopy2_P.Step_YFinal_f;
      }

      /* End of Step: '<S85>/Step' */

      /* Switch: '<S85>/Switch' incorporates:
       *  Gain: '<S85>/Gain'
       *  Sum: '<S85>/Sum'
       */
      if (rtb_Sum_l_idx_1 >= ModelCopy2_P.Switch_Threshold_j) {
        rtb_Product2_dm = 1.0 / ModelCopy2_P.MeanValue1_Period *
          (rtb_integrator_m - rtb_TransportDelay_a);
      } else {
        rtb_Product2_dm = ModelCopy2_B.Droop;
      }

      /* End of Switch: '<S85>/Switch' */

      /* Sum: '<S67>/Sum' */
      rtb_Product2_dm = (ModelCopy2_B.RateLimiter_b5 - rtb_Product2_dm) -
        rtb_Sum7_p;

      /* Gain: '<S86>/Gain1' */
      rtb_Switch3_a = ModelCopy2_P.Subsystem1_Kp_e * rtb_Product2_dm;

      /* Integrator: '<S86>/Integrator'
       *
       * Regarding '<S86>/Integrator':
       *  Limited Integrator
       */
      if (ModelCopy2_DW.Integrator_IWORK.IcNeedsLoading) {
        ModelCopy2_X.Integrator_CSTATE_p = ModelCopy2_B.ComplextoRealImag_o1_b;
      }

      if (ModelCopy2_X.Integrator_CSTATE_p >= ModelCopy2_P.Integrator_UpperSat_i
          ) {
        ModelCopy2_X.Integrator_CSTATE_p = ModelCopy2_P.Integrator_UpperSat_i;
      } else if (ModelCopy2_X.Integrator_CSTATE_p <=
                 (ModelCopy2_P.Integrator_LowerSat_o) ) {
        ModelCopy2_X.Integrator_CSTATE_p = (ModelCopy2_P.Integrator_LowerSat_o);
      }

      rtb_Sum7_p = ModelCopy2_X.Integrator_CSTATE_p;

      /* Sum: '<S86>/Sum7' */
      rtb_Sum7_p += rtb_Switch3_a;

      /* Saturate: '<S86>/Saturation' */
      if (rtb_Sum7_p > ModelCopy2_P.Saturation_UpperSat_h) {
        rtb_Switch3_a = ModelCopy2_P.Saturation_UpperSat_h;
      } else if (rtb_Sum7_p < ModelCopy2_P.Saturation_LowerSat_k) {
        rtb_Switch3_a = ModelCopy2_P.Saturation_LowerSat_k;
      } else {
        rtb_Switch3_a = rtb_Sum7_p;
      }

      /* End of Saturate: '<S86>/Saturation' */

      /* RateLimiter: '<S67>/Rate Limiter' */
      if (ModelCopy2_DW.LastMajorTime_a == (rtInf)) {
        ModelCopy2_B.RateLimiter_g = rtb_Switch3_a;
      } else {
        im = ModelCopy2_M->Timing.t[0] - ModelCopy2_DW.LastMajorTime_a;
        re = im * ModelCopy2_P.RateLimiter_RisingLim_as;
        rtb_Sum_l_idx_2 = rtb_Switch3_a - ModelCopy2_DW.PrevY_hi;
        if (rtb_Sum_l_idx_2 > re) {
          ModelCopy2_B.RateLimiter_g = ModelCopy2_DW.PrevY_hi + re;
        } else {
          im *= ModelCopy2_P.RateLimiter_FallingLim_bd;
          if (rtb_Sum_l_idx_2 < im) {
            ModelCopy2_B.RateLimiter_g = ModelCopy2_DW.PrevY_hi + im;
          } else {
            ModelCopy2_B.RateLimiter_g = rtb_Switch3_a;
          }
        }
      }

      /* End of RateLimiter: '<S67>/Rate Limiter' */

      /* Gain: '<S86>/Gain' */
      ModelCopy2_B.Gain_i = ModelCopy2_P.Subsystem1_Ki_i * rtb_Product2_dm;
    }

    /* End of Outputs for SubSystem: '<S49>/V Regulator' */

    /* Outputs for Enabled SubSystem: '<S49>/Q Regulator' incorporates:
     *  EnablePort: '<S66>/Enable'
     */
    if (rtmIsMajorTimeStep(ModelCopy2_M) && rtmIsMajorTimeStep(ModelCopy2_M)) {
      /* Logic: '<S49>/Logical Operator' incorporates:
       *  Constant: '<S49>/ExternalTorque1'
       */
      if (!(ModelCopy2_P.ExternalTorque1_Value != 0.0)) {
        if (!ModelCopy2_DW.QRegulator_MODE) {
          ModelCopy2_DW.QRegulator_MODE = true;
        }
      } else {
        if (ModelCopy2_DW.QRegulator_MODE) {
          ModelCopy2_DW.QRegulator_MODE = false;
        }
      }

      /* End of Logic: '<S49>/Logical Operator' */
    }

    if (ModelCopy2_DW.QRegulator_MODE) {
      /* Integrator: '<S82>/integrator' */
      rtb_integrator_c = ModelCopy2_X.integrator_CSTATE_f;

      /* TransportDelay: '<S82>/Transport Delay' */
      {
        real_T **uBuffer = (real_T**)
          &ModelCopy2_DW.TransportDelay_PWORK_a.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)
          &ModelCopy2_DW.TransportDelay_PWORK_a.TUbufferPtrs[1];
        real_T simTime = ModelCopy2_M->Timing.t[0];
        real_T tMinusDelay = simTime - ModelCopy2_P.MeanValue_Period_i;
        rtb_TransportDelay_e = rt_TDelayInterpolate(
          tMinusDelay,
          0.0,
          *tBuffer,
          *uBuffer,
          ModelCopy2_DW.TransportDelay_IWORK_h.CircularBufSize,
          &ModelCopy2_DW.TransportDelay_IWORK_h.Last,
          ModelCopy2_DW.TransportDelay_IWORK_h.Tail,
          ModelCopy2_DW.TransportDelay_IWORK_h.Head,
          ModelCopy2_P.TransportDelay_InitOutput_d,
          0,
          0);
      }

      /* RateLimiter: '<S66>/Rate Limiter ' incorporates:
       *  Constant: '<S4>/Qref '
       */
      if (ModelCopy2_DW.LastMajorTime_m == (rtInf)) {
        ModelCopy2_B.RateLimiter_nz = ModelCopy2_P.Qref_Value;
      } else {
        im = ModelCopy2_M->Timing.t[0] - ModelCopy2_DW.LastMajorTime_m;
        re = im * ModelCopy2_P.RateLimiter_RisingLim_dv;
        rtb_Sum_l_idx_2 = ModelCopy2_P.Qref_Value - ModelCopy2_DW.PrevY_n;
        if (rtb_Sum_l_idx_2 > re) {
          ModelCopy2_B.RateLimiter_nz = ModelCopy2_DW.PrevY_n + re;
        } else {
          im *= ModelCopy2_P.RateLimiter_FallingLim_g;
          if (rtb_Sum_l_idx_2 < im) {
            ModelCopy2_B.RateLimiter_nz = ModelCopy2_DW.PrevY_n + im;
          } else {
            ModelCopy2_B.RateLimiter_nz = ModelCopy2_P.Qref_Value;
          }
        }
      }

      /* End of RateLimiter: '<S66>/Rate Limiter ' */

      /* Step: '<S82>/Step' */
      if (ModelCopy2_M->Timing.t[0] < ModelCopy2_P.MeanValue_Period_i +
          2.2204460492503131E-16) {
        rtb_Sum_l_idx_1 = ModelCopy2_P.Step_Y0_e;
      } else {
        rtb_Sum_l_idx_1 = ModelCopy2_P.Step_YFinal_i;
      }

      /* End of Step: '<S82>/Step' */

      /* Switch: '<S82>/Switch' incorporates:
       *  Gain: '<S82>/Gain'
       *  Sum: '<S82>/Sum'
       */
      if (rtb_Sum_l_idx_1 >= ModelCopy2_P.Switch_Threshold_p) {
        rtb_Product2_dm = 1.0 / ModelCopy2_P.MeanValue_Period_i *
          (rtb_integrator_c - rtb_TransportDelay_e);
      } else {
        rtb_Product2_dm = ModelCopy2_B.RateLimiter_nz;
      }

      /* End of Switch: '<S82>/Switch' */

      /* Sum: '<S66>/Sum' */
      rtb_Product2_dm = ModelCopy2_B.RateLimiter_nz - rtb_Product2_dm;

      /* Gain: '<S83>/Gain1' */
      rtb_Switch3_a = ModelCopy2_P.Subsystem1_Kp_a * rtb_Product2_dm;

      /* Integrator: '<S83>/Integrator'
       *
       * Regarding '<S83>/Integrator':
       *  Limited Integrator
       */
      if (ModelCopy2_DW.Integrator_IWORK_o.IcNeedsLoading) {
        ModelCopy2_X.Integrator_CSTATE_a = ModelCopy2_B.ComplextoRealImag_o1_b;
      }

      if (ModelCopy2_X.Integrator_CSTATE_a >=
          ModelCopy2_P.Integrator_UpperSat_gb ) {
        ModelCopy2_X.Integrator_CSTATE_a = ModelCopy2_P.Integrator_UpperSat_gb;
      } else if (ModelCopy2_X.Integrator_CSTATE_a <=
                 (ModelCopy2_P.Integrator_LowerSat_n) ) {
        ModelCopy2_X.Integrator_CSTATE_a = (ModelCopy2_P.Integrator_LowerSat_n);
      }

      rtb_Sum7_c = ModelCopy2_X.Integrator_CSTATE_a;

      /* Sum: '<S83>/Sum7' */
      rtb_Sum7_c += rtb_Switch3_a;

      /* Saturate: '<S83>/Saturation' */
      if (rtb_Sum7_c > ModelCopy2_P.Saturation_UpperSat_i) {
        rtb_Switch3_a = ModelCopy2_P.Saturation_UpperSat_i;
      } else if (rtb_Sum7_c < ModelCopy2_P.Saturation_LowerSat_m) {
        rtb_Switch3_a = ModelCopy2_P.Saturation_LowerSat_m;
      } else {
        rtb_Switch3_a = rtb_Sum7_c;
      }

      /* End of Saturate: '<S83>/Saturation' */

      /* RateLimiter: '<S66>/Rate Limiter' */
      if (ModelCopy2_DW.LastMajorTime_m5 == (rtInf)) {
        ModelCopy2_B.RateLimiter_j = rtb_Switch3_a;
      } else {
        im = ModelCopy2_M->Timing.t[0] - ModelCopy2_DW.LastMajorTime_m5;
        re = im * ModelCopy2_P.RateLimiter_RisingLim_i;
        rtb_Sum_l_idx_2 = rtb_Switch3_a - ModelCopy2_DW.PrevY_c;
        if (rtb_Sum_l_idx_2 > re) {
          ModelCopy2_B.RateLimiter_j = ModelCopy2_DW.PrevY_c + re;
        } else {
          im *= ModelCopy2_P.RateLimiter_FallingLim_e;
          if (rtb_Sum_l_idx_2 < im) {
            ModelCopy2_B.RateLimiter_j = ModelCopy2_DW.PrevY_c + im;
          } else {
            ModelCopy2_B.RateLimiter_j = rtb_Switch3_a;
          }
        }
      }

      /* End of RateLimiter: '<S66>/Rate Limiter' */

      /* Gain: '<S83>/Gain' */
      ModelCopy2_B.Gain_j = ModelCopy2_P.Subsystem1_Ki_c * rtb_Product2_dm;
    }

    /* End of Outputs for SubSystem: '<S49>/Q Regulator' */

    /* Switch: '<S49>/Switch' incorporates:
     *  Constant: '<S49>/ExternalTorque1'
     */
    if (ModelCopy2_P.ExternalTorque1_Value >= ModelCopy2_P.Switch_Threshold_nd)
    {
      rtb_Sum7_l = ModelCopy2_B.RateLimiter_g;
    } else {
      rtb_Sum7_l = ModelCopy2_B.RateLimiter_j;
    }

    /* End of Switch: '<S49>/Switch' */

    /* Outputs for Enabled SubSystem: '<S64>/Subsystem' incorporates:
     *  EnablePort: '<S71>/Enable'
     */
    if (rtmIsMajorTimeStep(ModelCopy2_M) && rtmIsMajorTimeStep(ModelCopy2_M)) {
      /* Constant: '<S64>/ExternalTorque' */
      if (ModelCopy2_P.ExternalTorque_Value > 0.0) {
        if (!ModelCopy2_DW.Subsystem_MODE) {
          ModelCopy2_DW.Subsystem_MODE = true;
        }
      } else {
        if (ModelCopy2_DW.Subsystem_MODE) {
          ModelCopy2_DW.Subsystem_MODE = false;
        }
      }
    }

    if (ModelCopy2_DW.Subsystem_MODE) {
      if (rtmIsMajorTimeStep(ModelCopy2_M)) {
        /* Product: '<S76>/Product' incorporates:
         *  Constant: '<S76>/Constant2'
         *  Constant: '<S76>/Constant3'
         *  Constant: '<S76>/Constant4'
         *  Constant: '<S76>/Constant5'
         *  Sum: '<S76>/Sum6'
         *  Sum: '<S76>/Sum7'
         */
        ModelCopy2_B.Product = (ModelCopy2_P.Constant2_Value_p -
          ModelCopy2_P.WindTurbineDoublyFedInduction_k) /
          (ModelCopy2_P.Constant4_Value_h - ModelCopy2_P.Constant5_Value);

        /* Product: '<S75>/Product2' incorporates:
         *  Constant: '<S75>/Constant6'
         *  Constant: '<S75>/Constant7'
         *  Constant: '<S75>/Constant8'
         *  Constant: '<S75>/Constant9'
         *  Sum: '<S75>/Sum8'
         *  Sum: '<S75>/Sum9'
         */
        ModelCopy2_B.Product2 = (ModelCopy2_P.Constant6_Value -
          ModelCopy2_P.Constant7_Value) / (ModelCopy2_P.Constant8_Value -
          ModelCopy2_P.Constant9_Value);
      }

      /* Switch: '<S71>/Switch' incorporates:
       *  Constant: '<S71>/Constant1'
       *  Product: '<S71>/Product3'
       *  Saturate: '<S71>/speed_A-speed_B'
       *  Sum: '<S71>/Sum5'
       */
      if (rtb_Sum_l_idx_0 >= ModelCopy2_P.Switch_Threshold_n) {
        /* Gain: '<S71>/Gain ' incorporates:
         *  Fcn: '<S71>/wm^3'
         *  Gain: '<S71>/pu->pu '
         */
        rtb_Product2_dm = rt_powd_snf(ModelCopy2_P.pupu_Gain * rtb_Sum_l_idx_0,
          3.0) * ModelCopy2_P.WindTurbineDoublyFedInduction_k;

        /* Saturate: '<S71>/0-power_C' */
        if (rtb_Product2_dm > ModelCopy2_P.WindTurbineDoublyFedInduction_k) {
          rtb_Product2_dm = ModelCopy2_P.WindTurbineDoublyFedInduction_k;
        } else {
          if (rtb_Product2_dm < ModelCopy2_P.power_C_LowerSat) {
            rtb_Product2_dm = ModelCopy2_P.power_C_LowerSat;
          }
        }

        /* End of Saturate: '<S71>/0-power_C' */
      } else {
        if (rtb_Sum_l_idx_0 > ModelCopy2_P.speed_Aspeed_B_UpperSat) {
          /* Saturate: '<S71>/speed_A-speed_B' */
          rtb_Sum_l_idx_2 = ModelCopy2_P.speed_Aspeed_B_UpperSat;
        } else if (rtb_Sum_l_idx_0 < ModelCopy2_P.speed_Aspeed_B_LowerSat) {
          /* Saturate: '<S71>/speed_A-speed_B' */
          rtb_Sum_l_idx_2 = ModelCopy2_P.speed_Aspeed_B_LowerSat;
        } else {
          /* Saturate: '<S71>/speed_A-speed_B' */
          rtb_Sum_l_idx_2 = rtb_Sum_l_idx_0;
        }

        rtb_Product2_dm = (rtb_Sum_l_idx_2 - ModelCopy2_P.Constant1_Value_b) *
          ModelCopy2_B.Product2;
      }

      /* End of Switch: '<S71>/Switch' */

      /* Sum: '<S71>/Sum2' incorporates:
       *  Constant: '<S71>/Constant'
       */
      rtb_Sum_l_idx_1 = rtb_Sum_l_idx_0 - ModelCopy2_P.Constant_Value;

      /* Saturate: '<S71>/0-inf' */
      if (rtb_Sum_l_idx_1 > ModelCopy2_P.inf_UpperSat) {
        rtb_Sum_l_idx_1 = ModelCopy2_P.inf_UpperSat;
      } else {
        if (rtb_Sum_l_idx_1 < ModelCopy2_P.inf_LowerSat) {
          rtb_Sum_l_idx_1 = ModelCopy2_P.inf_LowerSat;
        }
      }

      /* Sum: '<S71>/Sum3' incorporates:
       *  Product: '<S71>/Product1'
       *  Saturate: '<S71>/0-inf'
       */
      rtb_Switch3_a = ModelCopy2_B.Product * rtb_Sum_l_idx_1 + rtb_Product2_dm;

      /* RateLimiter: '<S71>/Rate Limiter ' */
      if (ModelCopy2_DW.LastMajorTime_i5 == (rtInf)) {
        ModelCopy2_B.RateLimiter_e2 = rtb_Switch3_a;
      } else {
        im = ModelCopy2_M->Timing.t[0] - ModelCopy2_DW.LastMajorTime_i5;
        re = im * ModelCopy2_P.RateLimiter_RisingLim_e;
        rtb_Sum_l_idx_2 = rtb_Switch3_a - ModelCopy2_DW.PrevY_p;
        if (rtb_Sum_l_idx_2 > re) {
          ModelCopy2_B.RateLimiter_e2 = ModelCopy2_DW.PrevY_p + re;
        } else {
          im *= ModelCopy2_P.RateLimiter_FallingLim_de;
          if (rtb_Sum_l_idx_2 < im) {
            ModelCopy2_B.RateLimiter_e2 = ModelCopy2_DW.PrevY_p + im;
          } else {
            ModelCopy2_B.RateLimiter_e2 = rtb_Switch3_a;
          }
        }
      }

      /* End of RateLimiter: '<S71>/Rate Limiter ' */

      /* Saturate: '<S71>/0-1' */
      if (ModelCopy2_B.RateLimiter_e2 > ModelCopy2_P.u_UpperSat_m) {
        rtb_Llr1_o = ModelCopy2_P.u_UpperSat_m;
      } else if (ModelCopy2_B.RateLimiter_e2 < ModelCopy2_P.u_LowerSat_h) {
        rtb_Llr1_o = ModelCopy2_P.u_LowerSat_h;
      } else {
        rtb_Llr1_o = ModelCopy2_B.RateLimiter_e2;
      }

      /* End of Saturate: '<S71>/0-1' */

      /* Gain: '<S71>/pu->pu  ' */
      rtb_Powerspeed = ModelCopy2_P.WindTurbineDoublyFedInductionGe / 870000.0;

      /* Integrator: '<S78>/integrator' */
      rtb_integrator_b = ModelCopy2_X.integrator_CSTATE_k;

      /* TransportDelay: '<S78>/Transport Delay' */
      {
        real_T **uBuffer = (real_T**)
          &ModelCopy2_DW.TransportDelay_PWORK_gm.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)
          &ModelCopy2_DW.TransportDelay_PWORK_gm.TUbufferPtrs[1];
        real_T simTime = ModelCopy2_M->Timing.t[0];
        real_T tMinusDelay = simTime - ModelCopy2_P.MeanValue_Period;
        rtb_TransportDelay_n = rt_TDelayInterpolate(
          tMinusDelay,
          0.0,
          *tBuffer,
          *uBuffer,
          ModelCopy2_DW.TransportDelay_IWORK_ao.CircularBufSize,
          &ModelCopy2_DW.TransportDelay_IWORK_ao.Last,
          ModelCopy2_DW.TransportDelay_IWORK_ao.Tail,
          ModelCopy2_DW.TransportDelay_IWORK_ao.Head,
          ModelCopy2_P.TransportDelay_InitOutput,
          0,
          0);
      }

      /* Fcn: '<S73>/id^2+iq^2' */
      rtb_Sum7_i = rt_powd_snf(rtb_Llr2, 2.0) + rt_powd_snf
        (rtb_ComplextoRealImag1_o2, 2.0);

      /* Product: '<S73>/Product3' incorporates:
       *  Constant: '<S73>/Rs'
       */
      rtb_Product2_dm = ModelCopy2_P.Rs_Value * rtb_Sum7_i;

      /* Fcn: '<S73>/id^2+iq^2 ' */
      rtb_Sum7_i = rt_powd_snf(ModelCopy2_B.ComplextoRealImag_o1_b, 2.0) +
        rt_powd_snf(ModelCopy2_B.ComplextoRealImag_o2_p, 2.0);

      /* Product: '<S73>/Product4' incorporates:
       *  Constant: '<S73>/Rr'
       */
      rtb_Switch3_a = ModelCopy2_P.Rr_Value * rtb_Sum7_i;

      /* Fcn: '<S73>/id^2+iq^2   ' incorporates:
       *  ComplexToRealImag: '<S61>/Complex to Real-Imag2'
       */
      rtb_Sum7_i = rt_powd_snf(rtb_donotdeletethisgain_re, 2.0) + rt_powd_snf
        (rtb_donotdeletethisgain_im, 2.0);

      /* Sum: '<S73>/Sum2' incorporates:
       *  Constant: '<S73>/Constant3'
       *  Gain: '<S73>/Friction Factor'
       *  Product: '<S73>/Product1'
       *  Product: '<S73>/Product2'
       *  Sum: '<S73>/Sum5'
       */
      ModelCopy2_B.Looses = (ModelCopy2_P.FrictionFactor_Gain * rtb_Sum_l_idx_0 *
        rtb_Sum_l_idx_0 + rtb_Product2_dm) + (ModelCopy2_P.Constant3_Value_p *
        rtb_Sum7_i + rtb_Switch3_a);

      /* Step: '<S78>/Step' */
      if (ModelCopy2_M->Timing.t[0] < ModelCopy2_P.MeanValue_Period +
          2.2204460492503131E-16) {
        rtb_Sum_l_idx_1 = ModelCopy2_P.Step_Y0;
      } else {
        rtb_Sum_l_idx_1 = ModelCopy2_P.Step_YFinal;
      }

      /* End of Step: '<S78>/Step' */

      /* Switch: '<S78>/Switch' incorporates:
       *  Gain: '<S78>/Gain'
       *  Sum: '<S78>/Sum'
       */
      if (rtb_Sum_l_idx_1 >= ModelCopy2_P.Switch_Threshold_f) {
        rtb_Product2_dm = 1.0 / ModelCopy2_P.MeanValue_Period *
          (rtb_integrator_b - rtb_TransportDelay_n);
      } else {
        rtb_Product2_dm = ModelCopy2_B.Looses;
      }

      /* End of Switch: '<S78>/Switch' */

      /* Sum: '<S71>/Sum4' incorporates:
       *  Gain: '<S71>/pu->pu  '
       */
      rtb_Product2_dm = rtb_Powerspeed * rtb_Llr1_o - rtb_Product2_dm;

      /* Saturate: '<S71>/0-inf ' */
      if (rtb_Product2_dm > ModelCopy2_P.inf_UpperSat_o) {
        rtb_Product2_dm = ModelCopy2_P.inf_UpperSat_o;
      } else {
        if (rtb_Product2_dm < ModelCopy2_P.inf_LowerSat_a) {
          rtb_Product2_dm = ModelCopy2_P.inf_LowerSat_a;
        }
      }

      /* End of Saturate: '<S71>/0-inf ' */

      /* Integrator: '<S74>/integrator' */
      rtb_integrator_mk = ModelCopy2_X.integrator_CSTATE_o;

      /* TransportDelay: '<S74>/Transport Delay' */
      {
        real_T **uBuffer = (real_T**)
          &ModelCopy2_DW.TransportDelay_PWORK_n.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)
          &ModelCopy2_DW.TransportDelay_PWORK_n.TUbufferPtrs[1];
        real_T simTime = ModelCopy2_M->Timing.t[0];
        real_T tMinusDelay = simTime - ModelCopy2_P.MeanValue_Period_o;
        rtb_TransportDelay_m = rt_TDelayInterpolate(
          tMinusDelay,
          0.0,
          *tBuffer,
          *uBuffer,
          ModelCopy2_DW.TransportDelay_IWORK_ak.CircularBufSize,
          &ModelCopy2_DW.TransportDelay_IWORK_ak.Last,
          ModelCopy2_DW.TransportDelay_IWORK_ak.Tail,
          ModelCopy2_DW.TransportDelay_IWORK_ak.Head,
          ModelCopy2_P.TransportDelay_InitOutput_e,
          0,
          0);
      }

      /* Step: '<S74>/Step' */
      if (ModelCopy2_M->Timing.t[0] < ModelCopy2_P.MeanValue_Period_o +
          2.2204460492503131E-16) {
        rtb_Sum7_i = ModelCopy2_P.Step_Y0_k;
      } else {
        rtb_Sum7_i = ModelCopy2_P.Step_YFinal_a;
      }

      /* End of Step: '<S74>/Step' */

      /* Switch: '<S74>/Switch' incorporates:
       *  Gain: '<S74>/Gain'
       *  Sum: '<S74>/Sum'
       */
      if (rtb_Sum7_i >= ModelCopy2_P.Switch_Threshold_k) {
        rtb_Sum7_i = 1.0 / ModelCopy2_P.MeanValue_Period_o * (rtb_integrator_mk
          - rtb_TransportDelay_m);
      } else {
        rtb_Sum7_i = rtb_Product2_dm;
      }

      /* End of Switch: '<S74>/Switch' */

      /* Sum: '<S71>/Sum1' */
      rtb_Product2_dm -= rtb_Sum7_i;

      /* Gain: '<S77>/Gain1' */
      rtb_Switch3_a = ModelCopy2_P.Subsystem2_Kp * rtb_Product2_dm;

      /* Integrator: '<S77>/Integrator'
       *
       * Regarding '<S77>/Integrator':
       *  Limited Integrator
       */
      if (ModelCopy2_DW.Integrator_IWORK_a.IcNeedsLoading) {
        ModelCopy2_X.Integrator_CSTATE_k = ModelCopy2_B.ComplextoRealImag_o2_p;
      }

      if (ModelCopy2_X.Integrator_CSTATE_k >= ModelCopy2_P.Integrator_UpperSat_a
          ) {
        ModelCopy2_X.Integrator_CSTATE_k = ModelCopy2_P.Integrator_UpperSat_a;
      } else if (ModelCopy2_X.Integrator_CSTATE_k <=
                 ModelCopy2_P.Integrator_LowerSat_b ) {
        ModelCopy2_X.Integrator_CSTATE_k = ModelCopy2_P.Integrator_LowerSat_b;
      }

      rtb_Sum7_i = ModelCopy2_X.Integrator_CSTATE_k;

      /* Sum: '<S77>/Sum7' */
      rtb_Sum7_i += rtb_Switch3_a;

      /* Saturate: '<S77>/Saturation' */
      if (rtb_Sum7_i > ModelCopy2_P.Saturation_UpperSat_b) {
        rtb_Switch3_a = ModelCopy2_P.Saturation_UpperSat_b;
      } else if (rtb_Sum7_i < ModelCopy2_P.Saturation_LowerSat_g) {
        rtb_Switch3_a = ModelCopy2_P.Saturation_LowerSat_g;
      } else {
        rtb_Switch3_a = rtb_Sum7_i;
      }

      /* End of Saturate: '<S77>/Saturation' */

      /* RateLimiter: '<S71>/Rate Limiter' */
      if (ModelCopy2_DW.LastMajorTime_l == (rtInf)) {
        ModelCopy2_B.RateLimiter_m = rtb_Switch3_a;
      } else {
        im = ModelCopy2_M->Timing.t[0] - ModelCopy2_DW.LastMajorTime_l;
        re = im * ModelCopy2_P.RateLimiter_RisingLim_n;
        rtb_Sum_l_idx_2 = rtb_Switch3_a - ModelCopy2_DW.PrevY_c2;
        if (rtb_Sum_l_idx_2 > re) {
          ModelCopy2_B.RateLimiter_m = ModelCopy2_DW.PrevY_c2 + re;
        } else {
          im *= ModelCopy2_P.RateLimiter_FallingLim_dw;
          if (rtb_Sum_l_idx_2 < im) {
            ModelCopy2_B.RateLimiter_m = ModelCopy2_DW.PrevY_c2 + im;
          } else {
            ModelCopy2_B.RateLimiter_m = rtb_Switch3_a;
          }
        }
      }

      /* End of RateLimiter: '<S71>/Rate Limiter' */

      /* Gain: '<S77>/Gain' */
      ModelCopy2_B.Gain_b = ModelCopy2_P.Subsystem2_Ki * rtb_Product2_dm;
    }

    /* End of Outputs for SubSystem: '<S64>/Subsystem' */

    /* Outputs for Enabled SubSystem: '<S64>/Subsystem ' incorporates:
     *  EnablePort: '<S72>/Enable'
     */
    if (rtmIsMajorTimeStep(ModelCopy2_M) && rtmIsMajorTimeStep(ModelCopy2_M)) {
      /* Logic: '<S64>/Logical Operator' incorporates:
       *  Constant: '<S64>/ExternalTorque'
       */
      if (!(ModelCopy2_P.ExternalTorque_Value != 0.0)) {
        if (!ModelCopy2_DW.Subsystem_MODE_o) {
          ModelCopy2_DW.Subsystem_MODE_o = true;
        }
      } else {
        if (ModelCopy2_DW.Subsystem_MODE_o) {
          ModelCopy2_DW.Subsystem_MODE_o = false;
        }
      }

      /* End of Logic: '<S64>/Logical Operator' */
    }

    if (ModelCopy2_DW.Subsystem_MODE_o) {
      /* Integrator: '<S79>/integrator' */
      rtb_integrator_mf = ModelCopy2_X.integrator_CSTATE_d;

      /* TransportDelay: '<S79>/Transport Delay' */
      {
        real_T **uBuffer = (real_T**)
          &ModelCopy2_DW.TransportDelay_PWORK_c.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)
          &ModelCopy2_DW.TransportDelay_PWORK_c.TUbufferPtrs[1];
        real_T simTime = ModelCopy2_M->Timing.t[0];
        real_T tMinusDelay = simTime - ModelCopy2_P.MeanValue_Period_h;
        rtb_TransportDelay_p = rt_TDelayInterpolate(
          tMinusDelay,
          0.0,
          *tBuffer,
          *uBuffer,
          ModelCopy2_DW.TransportDelay_IWORK_a.CircularBufSize,
          &ModelCopy2_DW.TransportDelay_IWORK_a.Last,
          ModelCopy2_DW.TransportDelay_IWORK_a.Tail,
          ModelCopy2_DW.TransportDelay_IWORK_a.Head,
          ModelCopy2_P.TransportDelay_InitOutput_c,
          0,
          0);
      }

      /* Step: '<S79>/Step' */
      if (ModelCopy2_M->Timing.t[0] < ModelCopy2_P.MeanValue_Period_h +
          2.2204460492503131E-16) {
        rtb_Sum_l_idx_1 = ModelCopy2_P.Step_Y0_i;
      } else {
        rtb_Sum_l_idx_1 = ModelCopy2_P.Step_YFinal_m;
      }

      /* End of Step: '<S79>/Step' */

      /* Switch: '<S79>/Switch' incorporates:
       *  Gain: '<S79>/Gain'
       *  Sum: '<S79>/Sum'
       */
      if (rtb_Sum_l_idx_1 >= ModelCopy2_P.Switch_Threshold_e) {
        rtb_Product2_dm = 1.0 / ModelCopy2_P.MeanValue_Period_h *
          (rtb_integrator_mf - rtb_TransportDelay_p);
      } else {
        rtb_Product2_dm = ModelCopy2_B.RateLimiter_b;
      }

      /* End of Switch: '<S79>/Switch' */

      /* Sum: '<S72>/Sum1' */
      rtb_Product2_dm = ModelCopy2_B.RateLimiter_b - rtb_Product2_dm;

      /* Gain: '<S80>/Gain1' */
      rtb_Switch3_a = ModelCopy2_P.Subsystem2_Kp_c * rtb_Product2_dm;

      /* Integrator: '<S80>/Integrator'
       *
       * Regarding '<S80>/Integrator':
       *  Limited Integrator
       */
      if (ModelCopy2_DW.Integrator_IWORK_l.IcNeedsLoading) {
        ModelCopy2_X.Integrator_CSTATE_e = ModelCopy2_B.ComplextoRealImag_o2_p;
      }

      if (ModelCopy2_X.Integrator_CSTATE_e >= ModelCopy2_P.Integrator_UpperSat_g
          ) {
        ModelCopy2_X.Integrator_CSTATE_e = ModelCopy2_P.Integrator_UpperSat_g;
      } else if (ModelCopy2_X.Integrator_CSTATE_e <=
                 ModelCopy2_P.Integrator_LowerSat_g ) {
        ModelCopy2_X.Integrator_CSTATE_e = ModelCopy2_P.Integrator_LowerSat_g;
      }

      rtb_Sum7_a = ModelCopy2_X.Integrator_CSTATE_e;

      /* Sum: '<S80>/Sum7' */
      rtb_Sum7_a += rtb_Switch3_a;

      /* Saturate: '<S80>/Saturation' */
      if (rtb_Sum7_a > ModelCopy2_P.Saturation_UpperSat_e) {
        rtb_Switch3_a = ModelCopy2_P.Saturation_UpperSat_e;
      } else if (rtb_Sum7_a < ModelCopy2_P.Saturation_LowerSat_j) {
        rtb_Switch3_a = ModelCopy2_P.Saturation_LowerSat_j;
      } else {
        rtb_Switch3_a = rtb_Sum7_a;
      }

      /* End of Saturate: '<S80>/Saturation' */

      /* RateLimiter: '<S72>/Rate Limiter' */
      if (ModelCopy2_DW.LastMajorTime_c == (rtInf)) {
        ModelCopy2_B.RateLimiter_h = rtb_Switch3_a;
      } else {
        im = ModelCopy2_M->Timing.t[0] - ModelCopy2_DW.LastMajorTime_c;
        re = im * ModelCopy2_P.RateLimiter_RisingLim_d;
        rtb_Sum_l_idx_2 = rtb_Switch3_a - ModelCopy2_DW.PrevY_hm;
        if (rtb_Sum_l_idx_2 > re) {
          ModelCopy2_B.RateLimiter_h = ModelCopy2_DW.PrevY_hm + re;
        } else {
          im *= ModelCopy2_P.RateLimiter_FallingLim_a;
          if (rtb_Sum_l_idx_2 < im) {
            ModelCopy2_B.RateLimiter_h = ModelCopy2_DW.PrevY_hm + im;
          } else {
            ModelCopy2_B.RateLimiter_h = rtb_Switch3_a;
          }
        }
      }

      /* End of RateLimiter: '<S72>/Rate Limiter' */

      /* Gain: '<S80>/Gain' */
      ModelCopy2_B.Gain_h = ModelCopy2_P.Subsystem2_Ki_d * rtb_Product2_dm;
    }

    /* End of Outputs for SubSystem: '<S64>/Subsystem ' */

    /* Switch: '<S64>/Switch' incorporates:
     *  Constant: '<S64>/ExternalTorque'
     */
    if (ModelCopy2_P.ExternalTorque_Value >= ModelCopy2_P.Switch_Threshold_i) {
      rtb_Switch3_a = ModelCopy2_B.RateLimiter_m;
    } else {
      rtb_Switch3_a = ModelCopy2_B.RateLimiter_h;
    }

    /* End of Switch: '<S64>/Switch' */

    /* Saturate: '<S64>/Saturation' */
    if (rtb_Switch3_a > ModelCopy2_P.Saturation_UpperSat_a) {
      rtb_Switch3_a = ModelCopy2_P.Saturation_UpperSat_a;
    } else {
      if (rtb_Switch3_a < ModelCopy2_P.Saturation_LowerSat_o) {
        rtb_Switch3_a = ModelCopy2_P.Saturation_LowerSat_o;
      }
    }

    /* End of Saturate: '<S64>/Saturation' */

    /* RateLimiter: '<S64>/Rate Limiter' */
    if (ModelCopy2_DW.LastMajorTime_k == (rtInf)) {
      ModelCopy2_B.RateLimiter_e = rtb_Switch3_a;
    } else {
      im = ModelCopy2_M->Timing.t[0] - ModelCopy2_DW.LastMajorTime_k;
      re = im * ModelCopy2_P.RateLimiter_RisingLim_o;
      rtb_Sum_l_idx_2 = rtb_Switch3_a - ModelCopy2_DW.PrevY_h;
      if (rtb_Sum_l_idx_2 > re) {
        ModelCopy2_B.RateLimiter_e = ModelCopy2_DW.PrevY_h + re;
      } else {
        im *= ModelCopy2_P.RateLimiter_FallingLim_es;
        if (rtb_Sum_l_idx_2 < im) {
          ModelCopy2_B.RateLimiter_e = ModelCopy2_DW.PrevY_h + im;
        } else {
          ModelCopy2_B.RateLimiter_e = rtb_Switch3_a;
        }
      }
    }

    /* End of RateLimiter: '<S64>/Rate Limiter' */

    /* Switch: '<S65>/Switch' incorporates:
     *  Fcn: '<S81>/x->r'
     *  Product: '<S65>/Product'
     *  Signum: '<S65>/Sign'
     */
    if (rt_hypotd_snf(rtb_Sum7_l, ModelCopy2_B.RateLimiter_e) >=
        ModelCopy2_P.Switch_Threshold_ew) {
      /* Sum: '<S65>/Sum' incorporates:
       *  Constant: '<S65>/Irotor_max^2'
       *  Math: '<S65>/Math Function'
       */
      rtb_Switch3_a = ModelCopy2_P.Irotor_max2_Value - rtb_Sum7_l * rtb_Sum7_l;

      /* Math: '<S65>/Math Function1'
       *
       * About '<S65>/Math Function1':
       *  Operator: sqrt
       */
      if (rtb_Switch3_a < 0.0) {
        rtb_Switch3_a = -sqrt(fabs(rtb_Switch3_a));
      } else {
        rtb_Switch3_a = sqrt(rtb_Switch3_a);
      }

      /* End of Math: '<S65>/Math Function1' */
      rtb_Product2_dm = rtb_Sum7_l;

      /* Signum: '<S65>/Sign' */
      if (ModelCopy2_B.RateLimiter_e < 0.0) {
        rtb_Sum_l_idx_1 = -1.0;
      } else if (ModelCopy2_B.RateLimiter_e > 0.0) {
        rtb_Sum_l_idx_1 = 1.0;
      } else if (ModelCopy2_B.RateLimiter_e == 0.0) {
        rtb_Sum_l_idx_1 = 0.0;
      } else {
        rtb_Sum_l_idx_1 = ModelCopy2_B.RateLimiter_e;
      }

      rtb_Sum_l_idx_2 = rtb_Sum_l_idx_1 * rtb_Switch3_a;
    } else {
      rtb_Product2_dm = rtb_Sum7_l;
      rtb_Sum_l_idx_2 = ModelCopy2_B.RateLimiter_e;
    }

    /* End of Switch: '<S65>/Switch' */

    /* RateLimiter: '<S65>/Rate Limiter' */
    if (ModelCopy2_DW.LastMajorTime_n == (rtInf)) {
      ModelCopy2_B.RateLimiter_n[0] = rtb_Product2_dm;
      ModelCopy2_B.RateLimiter_n[1] = rtb_Sum_l_idx_2;
    } else {
      im = ModelCopy2_M->Timing.t[0] - ModelCopy2_DW.LastMajorTime_n;
      re = im * ModelCopy2_P.RateLimiter_RisingLim_il;
      rtb_Sum_l_idx_1 = rtb_Product2_dm - ModelCopy2_DW.PrevY_g[0];
      rtb_Switch5_idx_0 = rtb_Sum_l_idx_2 - ModelCopy2_DW.PrevY_g[1];
      rtb_Switch3_a = im * ModelCopy2_P.RateLimiter_FallingLim_h;
      if (rtb_Sum_l_idx_1 > re) {
        ModelCopy2_B.RateLimiter_n[0] = ModelCopy2_DW.PrevY_g[0] + re;
      } else if (rtb_Sum_l_idx_1 < rtb_Switch3_a) {
        ModelCopy2_B.RateLimiter_n[0] = ModelCopy2_DW.PrevY_g[0] + rtb_Switch3_a;
      } else {
        ModelCopy2_B.RateLimiter_n[0] = rtb_Product2_dm;
      }

      if (rtb_Switch5_idx_0 > re) {
        ModelCopy2_B.RateLimiter_n[1] = ModelCopy2_DW.PrevY_g[1] + re;
      } else if (rtb_Switch5_idx_0 < rtb_Switch3_a) {
        ModelCopy2_B.RateLimiter_n[1] = ModelCopy2_DW.PrevY_g[1] + rtb_Switch3_a;
      } else {
        ModelCopy2_B.RateLimiter_n[1] = rtb_Sum_l_idx_2;
      }
    }

    /* End of RateLimiter: '<S65>/Rate Limiter' */

    /* Sum: '<S60>/Sum4' */
    rtb_Sum_b = ModelCopy2_B.RateLimiter_n[0] -
      ModelCopy2_B.ComplextoRealImag_o1_b;

    /* Sum: '<S60>/Sum3' */
    rtb_Sum7_l = ModelCopy2_B.RateLimiter_n[1] -
      ModelCopy2_B.ComplextoRealImag_o2_p;

    /* Gain: '<S69>/Gain' */
    ModelCopy2_B.Gain[0] = ModelCopy2_P.Subsystem1_Ki_e * rtb_Sum_b;
    ModelCopy2_B.Gain[1] = ModelCopy2_P.Subsystem1_Ki_e * rtb_Sum7_l;

    /* Gain: '<S69>/Gain1' */
    rtb_Product2_dm = ModelCopy2_P.Subsystem1_Kp_i * rtb_Sum_b;
    rtb_Sum_l_idx_2 = ModelCopy2_P.Subsystem1_Kp_i * rtb_Sum7_l;

    /* Integrator: '<S69>/Integrator'
     *
     * Regarding '<S69>/Integrator':
     *  Limited Integrator
     */
    if (ModelCopy2_X.Integrator_CSTATE_c4[0] >=
        ModelCopy2_P.Integrator_UpperSat_l ) {
      ModelCopy2_X.Integrator_CSTATE_c4[0] = ModelCopy2_P.Integrator_UpperSat_l;
    } else if (ModelCopy2_X.Integrator_CSTATE_c4[0] <=
               (ModelCopy2_P.Integrator_LowerSat_l) ) {
      ModelCopy2_X.Integrator_CSTATE_c4[0] = (ModelCopy2_P.Integrator_LowerSat_l);
    }

    if (ModelCopy2_X.Integrator_CSTATE_c4[1] >=
        ModelCopy2_P.Integrator_UpperSat_l ) {
      ModelCopy2_X.Integrator_CSTATE_c4[1] = ModelCopy2_P.Integrator_UpperSat_l;
    } else if (ModelCopy2_X.Integrator_CSTATE_c4[1] <=
               (ModelCopy2_P.Integrator_LowerSat_l) ) {
      ModelCopy2_X.Integrator_CSTATE_c4[1] = (ModelCopy2_P.Integrator_LowerSat_l);
    }

    rtb_Integrator[0] = ModelCopy2_X.Integrator_CSTATE_c4[0];
    rtb_Integrator[1] = ModelCopy2_X.Integrator_CSTATE_c4[1];

    /* Sum: '<S69>/Sum7' */
    rtb_Product2_dm += rtb_Integrator[0];
    rtb_Sum_l_idx_2 += rtb_Integrator[1];

    /* Saturate: '<S69>/Saturation' */
    if (rtb_Product2_dm > ModelCopy2_P.Saturation_UpperSat_k) {
      rtb_Product2_dm = ModelCopy2_P.Saturation_UpperSat_k;
    } else {
      if (rtb_Product2_dm < ModelCopy2_P.Saturation_LowerSat_p) {
        rtb_Product2_dm = ModelCopy2_P.Saturation_LowerSat_p;
      }
    }

    /* Sum: '<S60>/Sum5' incorporates:
     *  Sum: '<S60>/Sum2'
     */
    rtb_Sum_b = ((rtb_Product6 - rtb_Sum1_re) - rtb_a23_re) + rtb_Product2_dm;

    /* Saturate: '<S69>/Saturation' */
    if (rtb_Sum_l_idx_2 > ModelCopy2_P.Saturation_UpperSat_k) {
      rtb_Sum_l_idx_2 = ModelCopy2_P.Saturation_UpperSat_k;
    } else {
      if (rtb_Sum_l_idx_2 < ModelCopy2_P.Saturation_LowerSat_p) {
        rtb_Sum_l_idx_2 = ModelCopy2_P.Saturation_LowerSat_p;
      }
    }

    /* Sum: '<S60>/Sum7' incorporates:
     *  Sum: '<S60>/Sum6'
     */
    rtb_Sum7_l = ((rtb_a23_im + rtb_Product4) + rtb_Product5) + rtb_Sum_l_idx_2;

    /* Saturate: '<S68>/Avoid division by zero' */
    if (rtb_Vdc > ModelCopy2_P.Avoiddivisionbyzero_UpperSat_c) {
      rtb_Sum_l_idx_2 = ModelCopy2_P.Avoiddivisionbyzero_UpperSat_c;
    } else if (rtb_Vdc < ModelCopy2_P.Avoiddivisionbyzero_LowerSat_lp) {
      rtb_Sum_l_idx_2 = ModelCopy2_P.Avoiddivisionbyzero_LowerSat_lp;
    } else {
      rtb_Sum_l_idx_2 = rtb_Vdc;
    }

    /* Product: '<S68>/Product2' incorporates:
     *  Constant: '<S68>/K'
     *  Fcn: '<S87>/x->r'
     *  Saturate: '<S68>/Avoid division by zero'
     */
    rtb_Llr1_o = rt_hypotd_snf(rtb_Sum_b, rtb_Sum7_l) * ModelCopy2_P.K_Value_i /
      rtb_Sum_l_idx_2;

    /* Saturate: '<S68>/0-1' */
    if (rtb_Llr1_o > ModelCopy2_P.u_UpperSat_b) {
      rtb_Llr1_o = ModelCopy2_P.u_UpperSat_b;
    } else {
      if (rtb_Llr1_o < ModelCopy2_P.u_LowerSat_l) {
        rtb_Llr1_o = ModelCopy2_P.u_LowerSat_l;
      }
    }

    /* End of Saturate: '<S68>/0-1' */

    /* Fcn: '<S87>/x->theta' */
    rtb_Sum_b = rt_atan2d_snf(rtb_Sum7_l, rtb_Sum_b);

    /* Sum: '<S68>/Sum' incorporates:
     *  Sum: '<S49>/Sum'
     */
    rtb_Sum_b += ModelCopy2_B.u + ModelCopy2_B.u_g;

    /* MagnitudeAngleToComplex: '<S68>/Magnitude-Angle to Complex' */
    ModelCopy2_B.MagnitudeAngletoComplex_i.re = rtb_Llr1_o * cos(rtb_Sum_b);
    ModelCopy2_B.MagnitudeAngletoComplex_i.im = rtb_Llr1_o * sin(rtb_Sum_b);
  }

  /* End of Outputs for SubSystem: '<S31>/wind_dfig_rotor' */

  /* Gain: '<S46>/->pu ' incorporates:
   *  Gain: '<S46>/Gain1'
   *  Product: '<S46>/Product1'
   */
  rtb_donotdeletethisgain_re = rtb_Vdc *
    ModelCopy2_B.MagnitudeAngletoComplex_i.re * ModelCopy2_P.Gain1_Gain_h *
    ModelCopy2_P.pu_Gain_d;
  rtb_donotdeletethisgain_im = rtb_Vdc *
    ModelCopy2_B.MagnitudeAngletoComplex_i.im * ModelCopy2_P.Gain1_Gain_h *
    ModelCopy2_P.pu_Gain_d;

  /* ComplexToRealImag: '<S46>/Complex to Real-Imag1' */
  rtb_Product4 = rtb_donotdeletethisgain_re;
  rtb_Product5 = rtb_donotdeletethisgain_im;

  /* Switch: '<S41>/Switch' incorporates:
   *  ComplexToRealImag: '<S46>/Complex to Real-Imag1'
   *  Constant: '<S41>/Constant'
   *  Constant: '<S41>/Constant1'
   *  Gain: '<S41>/1\Llr1'
   *  Gain: '<S41>/1\Llr3'
   *  Gain: '<S41>/web'
   *  Gain: '<S41>/web1'
   *  Product: '<S41>/Product'
   *  Product: '<S41>/Product1'
   *  Sum: '<S41>/Sum1'
   *  Sum: '<S41>/Sum2'
   *  Switch: '<S41>/Switch3'
   */
  if (rtb_Compare) {
    ModelCopy2_B.Switch = ModelCopy2_P.Constant_Value_e;
    ModelCopy2_B.Switch3 = ModelCopy2_P.Constant1_Value_k;
  } else {
    ModelCopy2_B.Switch = ((rtb_donotdeletethisgain_im - rtb_Switch2 *
      rtb_Sum1_im) - ModelCopy2_P.Llr1_Gain * rtb_iqr) * ModelCopy2_P.web_Gain;
    ModelCopy2_B.Switch3 = ((rtb_Sum1_im * rtb_Switch1 +
      rtb_donotdeletethisgain_re) - ModelCopy2_P.Llr3_Gain * rtb_idr) *
      ModelCopy2_P.web1_Gain;
  }

  /* End of Switch: '<S41>/Switch' */

  /* Sum: '<S38>/Sum1' */
  rtb_Sum_l_idx_0 = rtb_Kv1_idx_1_re - rtb_Kv1_idx_2_re;
  rtb_Sum_l_idx_2 = rtb_Kv1_idx_1_im - rtb_Kv1_idx_2_im;

  /* Switch: '<S42>/Switch' incorporates:
   *  Constant: '<S37>/Constant7'
   *  Constant: '<S42>/Constant'
   *  Constant: '<S42>/Constant2'
   *  Gain: '<S39>/(a^2)//3'
   *  Gain: '<S39>/Gain'
   *  Gain: '<S42>/1\Llr1'
   *  Gain: '<S42>/1\Llr3'
   *  Gain: '<S42>/web'
   *  Gain: '<S42>/web1'
   *  Product: '<S42>/Product'
   *  Product: '<S42>/Product1'
   *  Sum: '<S38>/Sum'
   *  Sum: '<S39>/Sum3'
   *  Sum: '<S42>/Sum1'
   *  Sum: '<S42>/Sum2'
   *  Switch: '<S42>/Switch2'
   */
  if (rtb_Compare) {
    ModelCopy2_B.Switch_n = ModelCopy2_P.Constant_Value_o;
    ModelCopy2_B.Switch2_d = ModelCopy2_P.Constant2_Value_h;
  } else {
    ModelCopy2_B.Switch_n = ((((rtb_Kv1_idx_0_re - rtb_Kv1_idx_1_re) *
      ModelCopy2_P.Gain_Gain_b - (ModelCopy2_P.a23_Gain_a.re * rtb_Sum_l_idx_0 -
      ModelCopy2_P.a23_Gain_a.im * rtb_Sum_l_idx_2)) +
      ModelCopy2_P.Constant7_Value_d * rtb_Switch3) - ModelCopy2_P.Llr3_Gain_l *
      rtb_ids) * ModelCopy2_P.web1_Gain_h;
    ModelCopy2_B.Switch2_d = ((((rtb_Kv1_idx_0_im - rtb_Kv1_idx_1_im) *
      ModelCopy2_P.Gain_Gain_b - (ModelCopy2_P.a23_Gain_a.re * rtb_Sum_l_idx_2 +
      ModelCopy2_P.a23_Gain_a.im * rtb_Sum_l_idx_0)) - ModelCopy2_P.Llr1_Gain_o *
      rtb_iqs) - ModelCopy2_P.Constant7_Value_d * rtb_Switch1_c) *
      ModelCopy2_P.web_Gain_d;
  }

  /* End of Switch: '<S42>/Switch' */

  /* ComplexToRealImag: '<S47>/Complex to Real-Imag' incorporates:
   *  MagnitudeAngleToComplex: '<S47>/Magnitude-Angle to Complex'
   */
  rtb_ComplextoRealImag1_o2 = ModelCopy2_B.ComplextoMagnitudeAngle_o1 * cos
    (rtb_ComplextoMagnitudeAngle_o2);
  rtb_Sum_l_idx_0 = ModelCopy2_B.ComplextoMagnitudeAngle_o1 * sin
    (rtb_ComplextoMagnitudeAngle_o2);

  /* Outputs for Enabled SubSystem: '<S31>/wind_dfig_grid' incorporates:
   *  EnablePort: '<S48>/Enable'
   */
  if (rtmIsMajorTimeStep(ModelCopy2_M) && rtmIsMajorTimeStep(ModelCopy2_M)) {
    if (rtb_LogicalOperator > 0.0) {
      if (!ModelCopy2_DW.wind_dfig_grid_MODE) {
        ModelCopy2_DW.wind_dfig_grid_MODE = true;
      }
    } else {
      if (ModelCopy2_DW.wind_dfig_grid_MODE) {
        ModelCopy2_DW.wind_dfig_grid_MODE = false;
      }
    }
  }

  if (ModelCopy2_DW.wind_dfig_grid_MODE) {
    /* Sum: '<S52>/Sum' incorporates:
     *  ComplexToMagnitudeAngle: '<S52>/Complex to Magnitude-Angle'
     *  RealImagToComplex: '<S52>/Real-Imag to Complex'
     */
    rtb_Llr1_o = rt_atan2d_snf(rtb_Switch1_l_idx_1, rtb_Switch1_l_idx_0) -
      ModelCopy2_B.u;

    /* MagnitudeAngleToComplex: '<S52>/Magnitude-Angle to Complex' incorporates:
     *  ComplexToMagnitudeAngle: '<S52>/Complex to Magnitude-Angle'
     *  RealImagToComplex: '<S52>/Real-Imag to Complex'
     */
    rtb_donotdeletethisgain_re = rt_hypotd_snf(rtb_Switch1_l_idx_0,
      rtb_Switch1_l_idx_1) * cos(rtb_Llr1_o);
    rtb_donotdeletethisgain_im = rt_hypotd_snf(rtb_Switch1_l_idx_0,
      rtb_Switch1_l_idx_1) * sin(rtb_Llr1_o);

    /* Product: '<S50>/Product' incorporates:
     *  ComplexToRealImag: '<S52>/Complex to Real-Imag'
     *  Constant: '<S50>/Constant1'
     */
    rtb_Llr2 = ModelCopy2_P.Constant1_Value * rtb_donotdeletethisgain_im;

    /* Product: '<S50>/Product1' incorporates:
     *  ComplexToRealImag: '<S52>/Complex to Real-Imag'
     *  Constant: '<S50>/Constant2'
     */
    rtb_wwr = rtb_donotdeletethisgain_re * ModelCopy2_P.Constant2_Value_c;

    /* Product: '<S50>/Product2' incorporates:
     *  ComplexToRealImag: '<S52>/Complex to Real-Imag'
     *  Constant: '<S50>/Constant3'
     */
    rtb_Sum_b = ModelCopy2_P.Constant3_Value * rtb_donotdeletethisgain_re;

    /* Product: '<S50>/Product3' incorporates:
     *  ComplexToRealImag: '<S52>/Complex to Real-Imag'
     *  Constant: '<S50>/Constant4'
     */
    rtb_Sum7_l = rtb_donotdeletethisgain_im * ModelCopy2_P.Constant4_Value;

    /* Sum: '<S51>/Sum7' incorporates:
     *  Constant: '<S48>/Vdc_ref (V)'
     */
    rtb_ComplextoRealImag_o2 = ModelCopy2_P.Vdc_refV_Value - rtb_Vdc;

    /* Gain: '<S57>/Gain1' */
    rtb_Product2_dm = ModelCopy2_P.Subsystem1_Kp * rtb_ComplextoRealImag_o2;

    /* Integrator: '<S57>/Integrator'
     *
     * Regarding '<S57>/Integrator':
     *  Limited Integrator
     */
    if (ModelCopy2_X.Integrator_CSTATE_e1 >= ModelCopy2_P.Integrator_UpperSat )
    {
      ModelCopy2_X.Integrator_CSTATE_e1 = ModelCopy2_P.Integrator_UpperSat;
    } else if (ModelCopy2_X.Integrator_CSTATE_e1 <=
               (ModelCopy2_P.Integrator_LowerSat) ) {
      ModelCopy2_X.Integrator_CSTATE_e1 = (ModelCopy2_P.Integrator_LowerSat);
    }

    rtb_Sum_co = ModelCopy2_X.Integrator_CSTATE_e1;

    /* Sum: '<S57>/Sum7' */
    rtb_Sum_co += rtb_Product2_dm;

    /* Saturate: '<S57>/Saturation' */
    if (rtb_Sum_co > ModelCopy2_P.Saturation_UpperSat) {
      rtb_Switch3_a = ModelCopy2_P.Saturation_UpperSat;
    } else if (rtb_Sum_co < ModelCopy2_P.Saturation_LowerSat) {
      rtb_Switch3_a = ModelCopy2_P.Saturation_LowerSat;
    } else {
      rtb_Switch3_a = rtb_Sum_co;
    }

    /* End of Saturate: '<S57>/Saturation' */

    /* RateLimiter: '<S51>/Rate Limiter' */
    if (ModelCopy2_DW.LastMajorTime_h == (rtInf)) {
      ModelCopy2_B.RateLimiter_ek = rtb_Switch3_a;
    } else {
      im = ModelCopy2_M->Timing.t[0] - ModelCopy2_DW.LastMajorTime_h;
      re = im * ModelCopy2_P.RateLimiter_RisingLim;
      rtb_Sum_l_idx_2 = rtb_Switch3_a - ModelCopy2_DW.PrevY_l;
      if (rtb_Sum_l_idx_2 > re) {
        ModelCopy2_B.RateLimiter_ek = ModelCopy2_DW.PrevY_l + re;
      } else {
        im *= ModelCopy2_P.RateLimiter_FallingLim;
        if (rtb_Sum_l_idx_2 < im) {
          ModelCopy2_B.RateLimiter_ek = ModelCopy2_DW.PrevY_l + im;
        } else {
          ModelCopy2_B.RateLimiter_ek = rtb_Switch3_a;
        }
      }
    }

    /* End of RateLimiter: '<S51>/Rate Limiter' */

    /* RateLimiter: '<S48>/Rate Limiter' incorporates:
     *  Constant: '<S4>/Iq_ref '
     */
    if (ModelCopy2_DW.LastMajorTime_ml == (rtInf)) {
      ModelCopy2_B.RateLimiter_l = ModelCopy2_P.Iq_ref_Value;
    } else {
      im = ModelCopy2_M->Timing.t[0] - ModelCopy2_DW.LastMajorTime_ml;
      re = im * ModelCopy2_P.RateLimiter_RisingLim_h;
      rtb_Sum_l_idx_2 = ModelCopy2_P.Iq_ref_Value - ModelCopy2_DW.PrevY_b;
      if (rtb_Sum_l_idx_2 > re) {
        ModelCopy2_B.RateLimiter_l = ModelCopy2_DW.PrevY_b + re;
      } else {
        im *= ModelCopy2_P.RateLimiter_FallingLim_b;
        if (rtb_Sum_l_idx_2 < im) {
          ModelCopy2_B.RateLimiter_l = ModelCopy2_DW.PrevY_b + im;
        } else {
          ModelCopy2_B.RateLimiter_l = ModelCopy2_P.Iq_ref_Value;
        }
      }
    }

    /* End of RateLimiter: '<S48>/Rate Limiter' */

    /* Fcn: '<S58>/x->r' */
    rtb_Sum_co = rt_hypotd_snf(ModelCopy2_B.RateLimiter_ek,
      ModelCopy2_B.RateLimiter_l);

    /* Switch: '<S53>/Switch' incorporates:
     *  Product: '<S53>/Product'
     *  Signum: '<S53>/Sign'
     */
    if (rtb_Sum_co >= ModelCopy2_P.Switch_Threshold) {
      /* Sum: '<S53>/Sum' incorporates:
       *  Constant: '<S53>/Igrid_conv_max^2'
       *  Math: '<S53>/Math Function'
       */
      rtb_Switch3_a = ModelCopy2_P.Igrid_conv_max2_Value -
        ModelCopy2_B.RateLimiter_ek * ModelCopy2_B.RateLimiter_ek;

      /* Math: '<S53>/Math Function1'
       *
       * About '<S53>/Math Function1':
       *  Operator: sqrt
       */
      if (rtb_Switch3_a < 0.0) {
        rtb_Switch3_a = -sqrt(fabs(rtb_Switch3_a));
      } else {
        rtb_Switch3_a = sqrt(rtb_Switch3_a);
      }

      /* End of Math: '<S53>/Math Function1' */
      rtb_Product2_dm = ModelCopy2_B.RateLimiter_ek;

      /* Signum: '<S53>/Sign' */
      if (ModelCopy2_B.RateLimiter_l < 0.0) {
        rtb_Sum_l_idx_1 = -1.0;
      } else if (ModelCopy2_B.RateLimiter_l > 0.0) {
        rtb_Sum_l_idx_1 = 1.0;
      } else if (ModelCopy2_B.RateLimiter_l == 0.0) {
        rtb_Sum_l_idx_1 = 0.0;
      } else {
        rtb_Sum_l_idx_1 = ModelCopy2_B.RateLimiter_l;
      }

      rtb_Sum_l_idx_2 = rtb_Sum_l_idx_1 * rtb_Switch3_a;
    } else {
      rtb_Product2_dm = ModelCopy2_B.RateLimiter_ek;
      rtb_Sum_l_idx_2 = ModelCopy2_B.RateLimiter_l;
    }

    /* End of Switch: '<S53>/Switch' */

    /* RateLimiter: '<S53>/Rate Limiter' */
    if (ModelCopy2_DW.LastMajorTime_ck == (rtInf)) {
      ModelCopy2_B.RateLimiter_jx[0] = rtb_Product2_dm;
      ModelCopy2_B.RateLimiter_jx[1] = rtb_Sum_l_idx_2;
    } else {
      im = ModelCopy2_M->Timing.t[0] - ModelCopy2_DW.LastMajorTime_ck;
      re = im * ModelCopy2_P.RateLimiter_RisingLim_p;
      rtb_Sum_l_idx_1 = rtb_Product2_dm - ModelCopy2_DW.PrevY_d[0];
      rtb_Switch5_idx_0 = rtb_Sum_l_idx_2 - ModelCopy2_DW.PrevY_d[1];
      rtb_Switch3_a = im * ModelCopy2_P.RateLimiter_FallingLim_d;
      if (rtb_Sum_l_idx_1 > re) {
        ModelCopy2_B.RateLimiter_jx[0] = ModelCopy2_DW.PrevY_d[0] + re;
      } else if (rtb_Sum_l_idx_1 < rtb_Switch3_a) {
        ModelCopy2_B.RateLimiter_jx[0] = ModelCopy2_DW.PrevY_d[0] +
          rtb_Switch3_a;
      } else {
        ModelCopy2_B.RateLimiter_jx[0] = rtb_Product2_dm;
      }

      if (rtb_Switch5_idx_0 > re) {
        ModelCopy2_B.RateLimiter_jx[1] = ModelCopy2_DW.PrevY_d[1] + re;
      } else if (rtb_Switch5_idx_0 < rtb_Switch3_a) {
        ModelCopy2_B.RateLimiter_jx[1] = ModelCopy2_DW.PrevY_d[1] +
          rtb_Switch3_a;
      } else {
        ModelCopy2_B.RateLimiter_jx[1] = rtb_Sum_l_idx_2;
      }
    }

    /* End of RateLimiter: '<S53>/Rate Limiter' */

    /* Sum: '<S50>/Sum4' incorporates:
     *  ComplexToRealImag: '<S52>/Complex to Real-Imag'
     */
    rtb_Llr1_o = ModelCopy2_B.RateLimiter_jx[0] - rtb_donotdeletethisgain_re;

    /* Sum: '<S50>/Sum2' incorporates:
     *  ComplexToRealImag: '<S52>/Complex to Real-Imag'
     */
    rtb_Powerspeed = ModelCopy2_B.RateLimiter_jx[1] - rtb_donotdeletethisgain_im;

    /* Gain: '<S56>/Gain' */
    ModelCopy2_B.Gain_iy[0] = ModelCopy2_P.Subsystem3_Ki * rtb_Llr1_o;
    ModelCopy2_B.Gain_iy[1] = ModelCopy2_P.Subsystem3_Ki * rtb_Powerspeed;

    /* Gain: '<S56>/Gain1' */
    rtb_Product2_dm = ModelCopy2_P.Subsystem3_Kp * rtb_Llr1_o;
    rtb_Sum_l_idx_2 = ModelCopy2_P.Subsystem3_Kp * rtb_Powerspeed;

    /* Integrator: '<S56>/Integrator'
     *
     * Regarding '<S56>/Integrator':
     *  Limited Integrator
     */
    if (ModelCopy2_X.Integrator_CSTATE_g[0] >=
        ModelCopy2_P.Integrator_UpperSat_d ) {
      ModelCopy2_X.Integrator_CSTATE_g[0] = ModelCopy2_P.Integrator_UpperSat_d;
    } else if (ModelCopy2_X.Integrator_CSTATE_g[0] <=
               (ModelCopy2_P.Integrator_LowerSat_d) ) {
      ModelCopy2_X.Integrator_CSTATE_g[0] = (ModelCopy2_P.Integrator_LowerSat_d);
    }

    if (ModelCopy2_X.Integrator_CSTATE_g[1] >=
        ModelCopy2_P.Integrator_UpperSat_d ) {
      ModelCopy2_X.Integrator_CSTATE_g[1] = ModelCopy2_P.Integrator_UpperSat_d;
    } else if (ModelCopy2_X.Integrator_CSTATE_g[1] <=
               (ModelCopy2_P.Integrator_LowerSat_d) ) {
      ModelCopy2_X.Integrator_CSTATE_g[1] = (ModelCopy2_P.Integrator_LowerSat_d);
    }

    rtb_Integrator_g[0] = ModelCopy2_X.Integrator_CSTATE_g[0];
    rtb_Integrator_g[1] = ModelCopy2_X.Integrator_CSTATE_g[1];

    /* Sum: '<S56>/Sum7' */
    rtb_Product2_dm += rtb_Integrator_g[0];
    rtb_Sum_l_idx_2 += rtb_Integrator_g[1];

    /* Saturate: '<S56>/Saturation' */
    if (rtb_Product2_dm > ModelCopy2_P.Saturation_UpperSat_m) {
      rtb_Product2_dm = ModelCopy2_P.Saturation_UpperSat_m;
    } else {
      if (rtb_Product2_dm < ModelCopy2_P.Saturation_LowerSat_c) {
        rtb_Product2_dm = ModelCopy2_P.Saturation_LowerSat_c;
      }
    }

    /* ComplexToMagnitudeAngle: '<S55>/Complex to Magnitude-Angle' incorporates:
     *  RealImagToComplex: '<S55>/Real-Imag to Complex'
     */
    rtb_Sum_co = rt_hypotd_snf(rtb_ComplextoRealImag1_o2, rtb_Sum_l_idx_0);

    /* Sum: '<S55>/Sum' incorporates:
     *  ComplexToMagnitudeAngle: '<S55>/Complex to Magnitude-Angle'
     *  RealImagToComplex: '<S55>/Real-Imag to Complex'
     */
    rtb_Llr1_o = rt_atan2d_snf(rtb_Sum_l_idx_0, rtb_ComplextoRealImag1_o2) -
      ModelCopy2_B.u;

    /* MagnitudeAngleToComplex: '<S55>/Magnitude-Angle to Complex' */
    rtb_donotdeletethisgain_im = rtb_Sum_co * sin(rtb_Llr1_o);

    /* ComplexToRealImag: '<S55>/Complex to Real-Imag' incorporates:
     *  MagnitudeAngleToComplex: '<S55>/Magnitude-Angle to Complex'
     */
    rtb_Sum_co *= cos(rtb_Llr1_o);

    /* Sum: '<S50>/Sum1' incorporates:
     *  Sum: '<S50>/Sum3'
     */
    rtb_Sum_co = ((rtb_Sum_co + rtb_Llr2) - rtb_Sum_b) - rtb_Product2_dm;

    /* Saturate: '<S56>/Saturation' */
    if (rtb_Sum_l_idx_2 > ModelCopy2_P.Saturation_UpperSat_m) {
      rtb_Sum_l_idx_2 = ModelCopy2_P.Saturation_UpperSat_m;
    } else {
      if (rtb_Sum_l_idx_2 < ModelCopy2_P.Saturation_LowerSat_c) {
        rtb_Sum_l_idx_2 = ModelCopy2_P.Saturation_LowerSat_c;
      }
    }

    /* Sum: '<S50>/Sum6' incorporates:
     *  ComplexToRealImag: '<S55>/Complex to Real-Imag'
     *  Sum: '<S50>/Sum5'
     */
    rtb_Llr1_o = ((rtb_donotdeletethisgain_im - rtb_Sum7_l) - rtb_wwr) -
      rtb_Sum_l_idx_2;

    /* Gain: '<S57>/Gain' */
    ModelCopy2_B.Gain_k = ModelCopy2_P.Subsystem1_Ki * rtb_ComplextoRealImag_o2;

    /* Saturate: '<S54>/Avoid division by zero' */
    if (rtb_Vdc > ModelCopy2_P.Avoiddivisionbyzero_UpperSat_p) {
      rtb_Sum_l_idx_2 = ModelCopy2_P.Avoiddivisionbyzero_UpperSat_p;
    } else if (rtb_Vdc < ModelCopy2_P.Avoiddivisionbyzero_LowerSat_l) {
      rtb_Sum_l_idx_2 = ModelCopy2_P.Avoiddivisionbyzero_LowerSat_l;
    } else {
      rtb_Sum_l_idx_2 = rtb_Vdc;
    }

    /* Product: '<S54>/Product1' incorporates:
     *  Constant: '<S54>/K'
     *  Fcn: '<S59>/x->r'
     *  Saturate: '<S54>/Avoid division by zero'
     */
    rtb_ComplextoRealImag_o2 = rt_hypotd_snf(rtb_Sum_co, rtb_Llr1_o) *
      ModelCopy2_P.K_Value / rtb_Sum_l_idx_2;

    /* Saturate: '<S54>/0-1' */
    if (rtb_ComplextoRealImag_o2 > ModelCopy2_P.u_UpperSat) {
      rtb_ComplextoRealImag_o2 = ModelCopy2_P.u_UpperSat;
    } else {
      if (rtb_ComplextoRealImag_o2 < ModelCopy2_P.u_LowerSat) {
        rtb_ComplextoRealImag_o2 = ModelCopy2_P.u_LowerSat;
      }
    }

    /* End of Saturate: '<S54>/0-1' */

    /* Fcn: '<S59>/x->theta' */
    rtb_Sum_co = rt_atan2d_snf(rtb_Llr1_o, rtb_Sum_co);

    /* Sum: '<S54>/Sum' */
    rtb_Sum_co += ModelCopy2_B.u;

    /* MagnitudeAngleToComplex: '<S54>/Magnitude-Angle to Complex' */
    ModelCopy2_B.MagnitudeAngletoComplex_c.re = rtb_ComplextoRealImag_o2 * cos
      (rtb_Sum_co);
    ModelCopy2_B.MagnitudeAngletoComplex_c.im = rtb_ComplextoRealImag_o2 * sin
      (rtb_Sum_co);
  }

  /* End of Outputs for SubSystem: '<S31>/wind_dfig_grid' */

  /* Gain: '<S46>/->pu' incorporates:
   *  Gain: '<S46>/Gain'
   *  Product: '<S46>/Product2'
   */
  rtb_donotdeletethisgain_re = rtb_Vdc *
    ModelCopy2_B.MagnitudeAngletoComplex_c.re * ModelCopy2_P.Gain_Gain_m *
    ModelCopy2_P.pu_Gain_l;
  rtb_donotdeletethisgain_im = rtb_Vdc *
    ModelCopy2_B.MagnitudeAngletoComplex_c.im * ModelCopy2_P.Gain_Gain_m *
    ModelCopy2_P.pu_Gain_l;

  /* Gain: '<S32>/deg->rad1' incorporates:
   *  ComplexToRealImag: '<S46>/Complex to Real-Imag'
   *  Gain: '<S32>/pu->W'
   *  Product: '<S32>/Product2'
   *  Product: '<S34>/Product1'
   *  Product: '<S34>/Product3'
   *  Sum: '<S32>/Sum2'
   *  Sum: '<S34>/Sum1'
   *  Sum: '<S34>/Sum3'
   */
  ModelCopy2_B.degrad1 = ((rtb_Switch1_l_idx_0 * rtb_donotdeletethisgain_re +
    rtb_Switch1_l_idx_1 * rtb_donotdeletethisgain_im) - (rtb_Product4 * rtb_idr
    + rtb_Product5 * rtb_iqr)) * ModelCopy2_P.DCbusmodel_Pnom / rtb_Vdc * (1.0 /
    ModelCopy2_P.DCbusmodel_C);

  /* Switch: '<S89>/Switch' incorporates:
   *  ComplexToRealImag: '<S46>/Complex to Real-Imag'
   *  Constant: '<S89>/Constant'
   *  Gain: '<S89>/R_choke'
   *  Gain: '<S89>/R_choke1'
   *  Gain: '<S89>/R_choke2'
   *  Gain: '<S89>/R_choke3'
   *  Sum: '<S89>/Sum'
   */
  if (rtb_Compare) {
    ModelCopy2_B.Switch_nz[0] = ModelCopy2_P.Constant_Value_k;
    ModelCopy2_B.Switch_nz[1] = ModelCopy2_P.Constant_Value_k;
  } else {
    /* Gain: '<S89>/R_choke3' */
    rtb_Powerspeed = 376.99111843077515 /
      ModelCopy2_P.dqaxismodelofa3phaseseriesRLb_m;
    ModelCopy2_B.Switch_nz[0] = (((rtb_ComplextoRealImag1_o2 -
      rtb_donotdeletethisgain_re) - ModelCopy2_P.dqaxismodelofa3phaseseriesRLb_p
      * rtb_Switch1_l_idx_0) + ModelCopy2_P.dqaxismodelofa3phaseseriesRLb_m *
      rtb_Switch1_l_idx_1 * ModelCopy2_P.R_choke2_Gain[0]) * rtb_Powerspeed;
    ModelCopy2_B.Switch_nz[1] = (((rtb_Sum_l_idx_0 - rtb_donotdeletethisgain_im)
      - ModelCopy2_P.dqaxismodelofa3phaseseriesRLb_p * rtb_Switch1_l_idx_1) +
      ModelCopy2_P.dqaxismodelofa3phaseseriesRLb_m * rtb_Switch1_l_idx_0 *
      ModelCopy2_P.R_choke2_Gain[1]) * rtb_Powerspeed;
  }

  /* End of Switch: '<S89>/Switch' */
  if (rtmIsMajorTimeStep(ModelCopy2_M)) {
    /* Switch: '<S10>/Switch1' incorporates:
     *  Constant: '<S10>/Constant'
     *  Constant: '<S10>/Constant2'
     *  Constant: '<S10>/valp_nom7'
     *  Constant: '<S6>/valp_nom5'
     *  RelationalOperator: '<S10>/Relational Operator'
     */
    if (ModelCopy2_P.valp_nom5_Value == ModelCopy2_P.Constant_Value_c) {
      ModelCopy2_B.Switch1 = ModelCopy2_P.VariationSubSystem_VariationRat;
    } else {
      ModelCopy2_B.Switch1 = ModelCopy2_P.Constant2_Value;
    }

    /* End of Switch: '<S10>/Switch1' */
  }

  /* Switch: '<S10>/Switch' incorporates:
   *  Constant: '<S10>/Constant4'
   */
  if (rtb_Ton >= ModelCopy2_P.Switch_Threshold_l) {
    ModelCopy2_B.Switch_j = ModelCopy2_B.Switch1;
  } else {
    ModelCopy2_B.Switch_j = ModelCopy2_P.Constant4_Value_jk;
  }

  /* End of Switch: '<S10>/Switch' */
  if (rtmIsMajorTimeStep(ModelCopy2_M)) {
    /* Update for Enabled SubSystem: '<S43>/SET  Priority' incorporates:
     *  Update for EnablePort: '<S45>/Enable'
     */
    if (ModelCopy2_DW.SETPriority_MODE && rtmIsMajorTimeStep(ModelCopy2_M)) {
      /* Update for Memory: '<S45>/IC=ic' */
      ModelCopy2_DW.ICic_PreviousInput_n = ModelCopy2_B.LogicalOperator_e;
    }

    /* End of Update for SubSystem: '<S43>/SET  Priority' */

    /* Update for Enabled SubSystem: '<S43>/RESET Priority' incorporates:
     *  Update for EnablePort: '<S44>/Enable'
     */
    if (ModelCopy2_DW.RESETPriority_MODE && rtmIsMajorTimeStep(ModelCopy2_M)) {
      /* Update for Memory: '<S44>/IC=ic' */
      ModelCopy2_DW.ICic_PreviousInput_c = ModelCopy2_B.LogicalOperator2_i;
    }

    /* End of Update for SubSystem: '<S43>/RESET Priority' */
    if (rtmIsMajorTimeStep(ModelCopy2_M)) {
      /* Update for Memory: '<S10>/Memory' */
      ModelCopy2_DW.Memory_PreviousInput = ModelCopy2_B.Switch2;
    }

    /* Update for RateLimiter: '<S31>/Rate Limiter   ' */
    ModelCopy2_DW.PrevY = ModelCopy2_B.RateLimiter;
    ModelCopy2_DW.LastMajorTime = ModelCopy2_M->Timing.t[0];
    if (rtmIsMajorTimeStep(ModelCopy2_M)) {
      /* Update for Memory: '<S47>/IC=ic' */
      ModelCopy2_DW.ICic_PreviousInput = ModelCopy2_B.u;
    }

    /* Update for RateLimiter: '<S31>/Rate Limiter ' */
    ModelCopy2_DW.PrevY_k = ModelCopy2_B.RateLimiter_b;
    ModelCopy2_DW.LastMajorTime_i = ModelCopy2_M->Timing.t[0];

    /* Update for Enabled SubSystem: '<S31>/wind_dfig_rotor' incorporates:
     *  Update for EnablePort: '<S49>/Enable'
     */
    if (ModelCopy2_DW.wind_dfig_rotor_MODE) {
      if (rtmIsMajorTimeStep(ModelCopy2_M)) {
        /* Update for Memory: '<S63>/IC=ic' */
        ModelCopy2_DW.ICic_PreviousInput_b = ModelCopy2_B.u_g;
      }

      /* Update for Enabled SubSystem: '<S49>/V Regulator' incorporates:
       *  Update for EnablePort: '<S67>/Enable'
       */
      if (ModelCopy2_DW.VRegulator_MODE) {
        /* Update for TransportDelay: '<S84>/Transport Delay' */
        {
          real_T **uBuffer = (real_T**)
            &ModelCopy2_DW.TransportDelay_PWORK.TUbufferPtrs[0];
          real_T **tBuffer = (real_T**)
            &ModelCopy2_DW.TransportDelay_PWORK.TUbufferPtrs[1];
          real_T simTime = ModelCopy2_M->Timing.t[0];
          ModelCopy2_DW.TransportDelay_IWORK.Head =
            ((ModelCopy2_DW.TransportDelay_IWORK.Head <
              (ModelCopy2_DW.TransportDelay_IWORK.CircularBufSize-1)) ?
             (ModelCopy2_DW.TransportDelay_IWORK.Head+1) : 0);
          if (ModelCopy2_DW.TransportDelay_IWORK.Head ==
              ModelCopy2_DW.TransportDelay_IWORK.Tail) {
            ModelCopy2_DW.TransportDelay_IWORK.Tail =
              ((ModelCopy2_DW.TransportDelay_IWORK.Tail <
                (ModelCopy2_DW.TransportDelay_IWORK.CircularBufSize-1)) ?
               (ModelCopy2_DW.TransportDelay_IWORK.Tail+1) : 0);
          }

          (*tBuffer)[ModelCopy2_DW.TransportDelay_IWORK.Head] = simTime;
          (*uBuffer)[ModelCopy2_DW.TransportDelay_IWORK.Head] = rtb_integrator;
        }

        /* Update for RateLimiter: '<S67>/Rate Limiter ' */
        ModelCopy2_DW.PrevY_i = ModelCopy2_B.RateLimiter_b5;
        ModelCopy2_DW.LastMajorTime_kl = ModelCopy2_M->Timing.t[0];

        /* Update for TransportDelay: '<S85>/Transport Delay' */
        {
          real_T **uBuffer = (real_T**)
            &ModelCopy2_DW.TransportDelay_PWORK_g.TUbufferPtrs[0];
          real_T **tBuffer = (real_T**)
            &ModelCopy2_DW.TransportDelay_PWORK_g.TUbufferPtrs[1];
          real_T simTime = ModelCopy2_M->Timing.t[0];
          ModelCopy2_DW.TransportDelay_IWORK_g.Head =
            ((ModelCopy2_DW.TransportDelay_IWORK_g.Head <
              (ModelCopy2_DW.TransportDelay_IWORK_g.CircularBufSize-1)) ?
             (ModelCopy2_DW.TransportDelay_IWORK_g.Head+1) : 0);
          if (ModelCopy2_DW.TransportDelay_IWORK_g.Head ==
              ModelCopy2_DW.TransportDelay_IWORK_g.Tail) {
            ModelCopy2_DW.TransportDelay_IWORK_g.Tail =
              ((ModelCopy2_DW.TransportDelay_IWORK_g.Tail <
                (ModelCopy2_DW.TransportDelay_IWORK_g.CircularBufSize-1)) ?
               (ModelCopy2_DW.TransportDelay_IWORK_g.Tail+1) : 0);
          }

          (*tBuffer)[ModelCopy2_DW.TransportDelay_IWORK_g.Head] = simTime;
          (*uBuffer)[ModelCopy2_DW.TransportDelay_IWORK_g.Head] =
            rtb_integrator_m;
        }

        /* Update for Integrator: '<S86>/Integrator' */
        ModelCopy2_DW.Integrator_IWORK.IcNeedsLoading = 0;

        /* Update for RateLimiter: '<S67>/Rate Limiter' */
        ModelCopy2_DW.PrevY_hi = ModelCopy2_B.RateLimiter_g;
        ModelCopy2_DW.LastMajorTime_a = ModelCopy2_M->Timing.t[0];
      }

      /* End of Update for SubSystem: '<S49>/V Regulator' */

      /* Update for Enabled SubSystem: '<S49>/Q Regulator' incorporates:
       *  Update for EnablePort: '<S66>/Enable'
       */
      if (ModelCopy2_DW.QRegulator_MODE) {
        /* Update for TransportDelay: '<S82>/Transport Delay' */
        {
          real_T **uBuffer = (real_T**)
            &ModelCopy2_DW.TransportDelay_PWORK_a.TUbufferPtrs[0];
          real_T **tBuffer = (real_T**)
            &ModelCopy2_DW.TransportDelay_PWORK_a.TUbufferPtrs[1];
          real_T simTime = ModelCopy2_M->Timing.t[0];
          ModelCopy2_DW.TransportDelay_IWORK_h.Head =
            ((ModelCopy2_DW.TransportDelay_IWORK_h.Head <
              (ModelCopy2_DW.TransportDelay_IWORK_h.CircularBufSize-1)) ?
             (ModelCopy2_DW.TransportDelay_IWORK_h.Head+1) : 0);
          if (ModelCopy2_DW.TransportDelay_IWORK_h.Head ==
              ModelCopy2_DW.TransportDelay_IWORK_h.Tail) {
            ModelCopy2_DW.TransportDelay_IWORK_h.Tail =
              ((ModelCopy2_DW.TransportDelay_IWORK_h.Tail <
                (ModelCopy2_DW.TransportDelay_IWORK_h.CircularBufSize-1)) ?
               (ModelCopy2_DW.TransportDelay_IWORK_h.Tail+1) : 0);
          }

          (*tBuffer)[ModelCopy2_DW.TransportDelay_IWORK_h.Head] = simTime;
          (*uBuffer)[ModelCopy2_DW.TransportDelay_IWORK_h.Head] =
            rtb_integrator_c;
        }

        /* Update for RateLimiter: '<S66>/Rate Limiter ' */
        ModelCopy2_DW.PrevY_n = ModelCopy2_B.RateLimiter_nz;
        ModelCopy2_DW.LastMajorTime_m = ModelCopy2_M->Timing.t[0];

        /* Update for Integrator: '<S83>/Integrator' */
        ModelCopy2_DW.Integrator_IWORK_o.IcNeedsLoading = 0;

        /* Update for RateLimiter: '<S66>/Rate Limiter' */
        ModelCopy2_DW.PrevY_c = ModelCopy2_B.RateLimiter_j;
        ModelCopy2_DW.LastMajorTime_m5 = ModelCopy2_M->Timing.t[0];
      }

      /* End of Update for SubSystem: '<S49>/Q Regulator' */

      /* Update for Enabled SubSystem: '<S64>/Subsystem' incorporates:
       *  Update for EnablePort: '<S71>/Enable'
       */
      if (ModelCopy2_DW.Subsystem_MODE) {
        /* Update for RateLimiter: '<S71>/Rate Limiter ' */
        ModelCopy2_DW.PrevY_p = ModelCopy2_B.RateLimiter_e2;
        ModelCopy2_DW.LastMajorTime_i5 = ModelCopy2_M->Timing.t[0];

        /* Update for TransportDelay: '<S78>/Transport Delay' */
        {
          real_T **uBuffer = (real_T**)
            &ModelCopy2_DW.TransportDelay_PWORK_gm.TUbufferPtrs[0];
          real_T **tBuffer = (real_T**)
            &ModelCopy2_DW.TransportDelay_PWORK_gm.TUbufferPtrs[1];
          real_T simTime = ModelCopy2_M->Timing.t[0];
          ModelCopy2_DW.TransportDelay_IWORK_ao.Head =
            ((ModelCopy2_DW.TransportDelay_IWORK_ao.Head <
              (ModelCopy2_DW.TransportDelay_IWORK_ao.CircularBufSize-1)) ?
             (ModelCopy2_DW.TransportDelay_IWORK_ao.Head+1) : 0);
          if (ModelCopy2_DW.TransportDelay_IWORK_ao.Head ==
              ModelCopy2_DW.TransportDelay_IWORK_ao.Tail) {
            ModelCopy2_DW.TransportDelay_IWORK_ao.Tail =
              ((ModelCopy2_DW.TransportDelay_IWORK_ao.Tail <
                (ModelCopy2_DW.TransportDelay_IWORK_ao.CircularBufSize-1)) ?
               (ModelCopy2_DW.TransportDelay_IWORK_ao.Tail+1) : 0);
          }

          (*tBuffer)[ModelCopy2_DW.TransportDelay_IWORK_ao.Head] = simTime;
          (*uBuffer)[ModelCopy2_DW.TransportDelay_IWORK_ao.Head] =
            rtb_integrator_b;
        }

        /* Update for TransportDelay: '<S74>/Transport Delay' */
        {
          real_T **uBuffer = (real_T**)
            &ModelCopy2_DW.TransportDelay_PWORK_n.TUbufferPtrs[0];
          real_T **tBuffer = (real_T**)
            &ModelCopy2_DW.TransportDelay_PWORK_n.TUbufferPtrs[1];
          real_T simTime = ModelCopy2_M->Timing.t[0];
          ModelCopy2_DW.TransportDelay_IWORK_ak.Head =
            ((ModelCopy2_DW.TransportDelay_IWORK_ak.Head <
              (ModelCopy2_DW.TransportDelay_IWORK_ak.CircularBufSize-1)) ?
             (ModelCopy2_DW.TransportDelay_IWORK_ak.Head+1) : 0);
          if (ModelCopy2_DW.TransportDelay_IWORK_ak.Head ==
              ModelCopy2_DW.TransportDelay_IWORK_ak.Tail) {
            ModelCopy2_DW.TransportDelay_IWORK_ak.Tail =
              ((ModelCopy2_DW.TransportDelay_IWORK_ak.Tail <
                (ModelCopy2_DW.TransportDelay_IWORK_ak.CircularBufSize-1)) ?
               (ModelCopy2_DW.TransportDelay_IWORK_ak.Tail+1) : 0);
          }

          (*tBuffer)[ModelCopy2_DW.TransportDelay_IWORK_ak.Head] = simTime;
          (*uBuffer)[ModelCopy2_DW.TransportDelay_IWORK_ak.Head] =
            rtb_integrator_mk;
        }

        /* Update for Integrator: '<S77>/Integrator' */
        ModelCopy2_DW.Integrator_IWORK_a.IcNeedsLoading = 0;

        /* Update for RateLimiter: '<S71>/Rate Limiter' */
        ModelCopy2_DW.PrevY_c2 = ModelCopy2_B.RateLimiter_m;
        ModelCopy2_DW.LastMajorTime_l = ModelCopy2_M->Timing.t[0];
      }

      /* End of Update for SubSystem: '<S64>/Subsystem' */

      /* Update for Enabled SubSystem: '<S64>/Subsystem ' incorporates:
       *  Update for EnablePort: '<S72>/Enable'
       */
      if (ModelCopy2_DW.Subsystem_MODE_o) {
        /* Update for TransportDelay: '<S79>/Transport Delay' */
        {
          real_T **uBuffer = (real_T**)
            &ModelCopy2_DW.TransportDelay_PWORK_c.TUbufferPtrs[0];
          real_T **tBuffer = (real_T**)
            &ModelCopy2_DW.TransportDelay_PWORK_c.TUbufferPtrs[1];
          real_T simTime = ModelCopy2_M->Timing.t[0];
          ModelCopy2_DW.TransportDelay_IWORK_a.Head =
            ((ModelCopy2_DW.TransportDelay_IWORK_a.Head <
              (ModelCopy2_DW.TransportDelay_IWORK_a.CircularBufSize-1)) ?
             (ModelCopy2_DW.TransportDelay_IWORK_a.Head+1) : 0);
          if (ModelCopy2_DW.TransportDelay_IWORK_a.Head ==
              ModelCopy2_DW.TransportDelay_IWORK_a.Tail) {
            ModelCopy2_DW.TransportDelay_IWORK_a.Tail =
              ((ModelCopy2_DW.TransportDelay_IWORK_a.Tail <
                (ModelCopy2_DW.TransportDelay_IWORK_a.CircularBufSize-1)) ?
               (ModelCopy2_DW.TransportDelay_IWORK_a.Tail+1) : 0);
          }

          (*tBuffer)[ModelCopy2_DW.TransportDelay_IWORK_a.Head] = simTime;
          (*uBuffer)[ModelCopy2_DW.TransportDelay_IWORK_a.Head] =
            rtb_integrator_mf;
        }

        /* Update for Integrator: '<S80>/Integrator' */
        ModelCopy2_DW.Integrator_IWORK_l.IcNeedsLoading = 0;

        /* Update for RateLimiter: '<S72>/Rate Limiter' */
        ModelCopy2_DW.PrevY_hm = ModelCopy2_B.RateLimiter_h;
        ModelCopy2_DW.LastMajorTime_c = ModelCopy2_M->Timing.t[0];
      }

      /* End of Update for SubSystem: '<S64>/Subsystem ' */

      /* Update for RateLimiter: '<S64>/Rate Limiter' */
      ModelCopy2_DW.PrevY_h = ModelCopy2_B.RateLimiter_e;
      ModelCopy2_DW.LastMajorTime_k = ModelCopy2_M->Timing.t[0];

      /* Update for RateLimiter: '<S65>/Rate Limiter' */
      ModelCopy2_DW.PrevY_g[0] = ModelCopy2_B.RateLimiter_n[0];
      ModelCopy2_DW.PrevY_g[1] = ModelCopy2_B.RateLimiter_n[1];
      ModelCopy2_DW.LastMajorTime_n = ModelCopy2_M->Timing.t[0];
    }

    /* End of Update for SubSystem: '<S31>/wind_dfig_rotor' */

    /* Update for Enabled SubSystem: '<S31>/wind_dfig_grid' incorporates:
     *  Update for EnablePort: '<S48>/Enable'
     */
    if (ModelCopy2_DW.wind_dfig_grid_MODE) {
      /* Update for RateLimiter: '<S51>/Rate Limiter' */
      ModelCopy2_DW.PrevY_l = ModelCopy2_B.RateLimiter_ek;
      ModelCopy2_DW.LastMajorTime_h = ModelCopy2_M->Timing.t[0];

      /* Update for RateLimiter: '<S48>/Rate Limiter' */
      ModelCopy2_DW.PrevY_b = ModelCopy2_B.RateLimiter_l;
      ModelCopy2_DW.LastMajorTime_ml = ModelCopy2_M->Timing.t[0];

      /* Update for RateLimiter: '<S53>/Rate Limiter' */
      ModelCopy2_DW.PrevY_d[0] = ModelCopy2_B.RateLimiter_jx[0];
      ModelCopy2_DW.PrevY_d[1] = ModelCopy2_B.RateLimiter_jx[1];
      ModelCopy2_DW.LastMajorTime_ck = ModelCopy2_M->Timing.t[0];
    }

    /* End of Update for SubSystem: '<S31>/wind_dfig_grid' */
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(ModelCopy2_M)) {
    rt_ertODEUpdateContinuousStates(&ModelCopy2_M->solverInfo, arg_VIND);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     */
    ++ModelCopy2_M->Timing.clockTick0;
    ModelCopy2_M->Timing.t[0] = rtsiGetSolverStopTime(&ModelCopy2_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.001s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.001, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       */
      ModelCopy2_M->Timing.clockTick1++;
    }
  }                                    /* end MajorTimeStep */

  return arg_PUT;
}

/* Derivatives for root system: '<Root>' */
void ModelCopy2_derivatives(real_T arg_VIND)
{
  XDot_ModelCopy2_T *_rtXdot;
  _rtXdot = ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs);

  /* Derivatives for Integrator: '<S41>/phiqr' */
  _rtXdot->phiqr_CSTATE = ModelCopy2_B.Switch;

  /* Derivatives for Integrator: '<S41>/phidr' */
  _rtXdot->phidr_CSTATE = ModelCopy2_B.Switch3;

  /* Derivatives for Integrator: '<S42>/phiqs' */
  _rtXdot->phiqs_CSTATE = ModelCopy2_B.Switch2_d;

  /* Derivatives for Integrator: '<S42>/phids' */
  _rtXdot->phids_CSTATE = ModelCopy2_B.Switch_n;

  /* Derivatives for Integrator: '<S89>/Integrator' */
  _rtXdot->Integrator_CSTATE[0] = ModelCopy2_B.Switch_nz[0];
  _rtXdot->Integrator_CSTATE[1] = ModelCopy2_B.Switch_nz[1];

  /* Derivatives for Integrator: '<S10>/Integrator' */
  _rtXdot->Integrator_CSTATE_j = ModelCopy2_B.Switch_j;

  /* Derivatives for Integrator: '<S35>/Integrator' */
  _rtXdot->Integrator_CSTATE_c = ModelCopy2_B._2H;

  /* Derivatives for Integrator: '<S35>/Integrator1' */
  _rtXdot->Integrator1_CSTATE = ModelCopy2_B.web;

  /* Derivatives for Integrator: '<S32>/Integrator' */
  _rtXdot->Integrator_CSTATE_cb = ModelCopy2_B.degrad1;

  /* Derivatives for Enabled SubSystem: '<S31>/wind_dfig_rotor' */
  if (ModelCopy2_DW.wind_dfig_rotor_MODE) {
    /* Derivatives for Enabled SubSystem: '<S49>/V Regulator' */
    if (ModelCopy2_DW.VRegulator_MODE) {
      /* Derivatives for Integrator: '<S84>/integrator' */
      _rtXdot->integrator_CSTATE = ModelCopy2_B.ComplextoMagnitudeAngle_o1;

      /* Derivatives for Integrator: '<S85>/integrator' */
      _rtXdot->integrator_CSTATE_p = ModelCopy2_B.Droop;

      /* Derivatives for Integrator: '<S86>/Integrator' */
      {
        boolean_T lsat;
        boolean_T usat;
        lsat = ( ModelCopy2_X.Integrator_CSTATE_p <=
                (ModelCopy2_P.Integrator_LowerSat_o) );
        usat = ( ModelCopy2_X.Integrator_CSTATE_p >=
                ModelCopy2_P.Integrator_UpperSat_i );
        if ((!lsat && !usat) ||
            (lsat && (ModelCopy2_B.Gain_i > 0)) ||
            (usat && (ModelCopy2_B.Gain_i < 0)) ) {
          ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
            ->Integrator_CSTATE_p = ModelCopy2_B.Gain_i;
        } else {
          /* in saturation */
          ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
            ->Integrator_CSTATE_p = 0.0;
        }
      }
    } else {
      {
        real_T *dx;
        int_T i;
        dx = &(((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
               ->integrator_CSTATE);
        for (i=0; i < 3; i++) {
          dx[i] = 0.0;
        }
      }
    }

    /* End of Derivatives for SubSystem: '<S49>/V Regulator' */

    /* Derivatives for Enabled SubSystem: '<S49>/Q Regulator' */
    if (ModelCopy2_DW.QRegulator_MODE) {
      /* Derivatives for Integrator: '<S82>/integrator' */
      _rtXdot->integrator_CSTATE_f = ModelCopy2_B.Qpu;

      /* Derivatives for Integrator: '<S83>/Integrator' */
      {
        boolean_T lsat;
        boolean_T usat;
        lsat = ( ModelCopy2_X.Integrator_CSTATE_a <=
                (ModelCopy2_P.Integrator_LowerSat_n) );
        usat = ( ModelCopy2_X.Integrator_CSTATE_a >=
                ModelCopy2_P.Integrator_UpperSat_gb );
        if ((!lsat && !usat) ||
            (lsat && (ModelCopy2_B.Gain_j > 0)) ||
            (usat && (ModelCopy2_B.Gain_j < 0)) ) {
          ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
            ->Integrator_CSTATE_a = ModelCopy2_B.Gain_j;
        } else {
          /* in saturation */
          ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
            ->Integrator_CSTATE_a = 0.0;
        }
      }
    } else {
      {
        real_T *dx;
        int_T i;
        dx = &(((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
               ->integrator_CSTATE_f);
        for (i=0; i < 2; i++) {
          dx[i] = 0.0;
        }
      }
    }

    /* End of Derivatives for SubSystem: '<S49>/Q Regulator' */

    /* Derivatives for Enabled SubSystem: '<S64>/Subsystem' */
    if (ModelCopy2_DW.Subsystem_MODE) {
      /* Derivatives for Integrator: '<S78>/integrator' */
      _rtXdot->integrator_CSTATE_k = ModelCopy2_B.Looses;

      /* Derivatives for Integrator: '<S74>/integrator' */
      _rtXdot->integrator_CSTATE_o = ModelCopy2_B.Ppu;

      /* Derivatives for Integrator: '<S77>/Integrator' */
      {
        boolean_T lsat;
        boolean_T usat;
        lsat = ( ModelCopy2_X.Integrator_CSTATE_k <=
                ModelCopy2_P.Integrator_LowerSat_b );
        usat = ( ModelCopy2_X.Integrator_CSTATE_k >=
                ModelCopy2_P.Integrator_UpperSat_a );
        if ((!lsat && !usat) ||
            (lsat && (ModelCopy2_B.Gain_b > 0)) ||
            (usat && (ModelCopy2_B.Gain_b < 0)) ) {
          ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
            ->Integrator_CSTATE_k = ModelCopy2_B.Gain_b;
        } else {
          /* in saturation */
          ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
            ->Integrator_CSTATE_k = 0.0;
        }
      }
    } else {
      {
        real_T *dx;
        int_T i;
        dx = &(((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
               ->integrator_CSTATE_k);
        for (i=0; i < 3; i++) {
          dx[i] = 0.0;
        }
      }
    }

    /* End of Derivatives for SubSystem: '<S64>/Subsystem' */

    /* Derivatives for Enabled SubSystem: '<S64>/Subsystem ' */
    if (ModelCopy2_DW.Subsystem_MODE_o) {
      /* Derivatives for Integrator: '<S79>/integrator' */
      _rtXdot->integrator_CSTATE_d = ModelCopy2_B.Ppu;

      /* Derivatives for Integrator: '<S80>/Integrator' */
      {
        boolean_T lsat;
        boolean_T usat;
        lsat = ( ModelCopy2_X.Integrator_CSTATE_e <=
                ModelCopy2_P.Integrator_LowerSat_g );
        usat = ( ModelCopy2_X.Integrator_CSTATE_e >=
                ModelCopy2_P.Integrator_UpperSat_g );
        if ((!lsat && !usat) ||
            (lsat && (ModelCopy2_B.Gain_h > 0)) ||
            (usat && (ModelCopy2_B.Gain_h < 0)) ) {
          ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
            ->Integrator_CSTATE_e = ModelCopy2_B.Gain_h;
        } else {
          /* in saturation */
          ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
            ->Integrator_CSTATE_e = 0.0;
        }
      }
    } else {
      {
        real_T *dx;
        int_T i;
        dx = &(((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
               ->integrator_CSTATE_d);
        for (i=0; i < 2; i++) {
          dx[i] = 0.0;
        }
      }
    }

    /* End of Derivatives for SubSystem: '<S64>/Subsystem ' */
    /* Derivatives for Integrator: '<S69>/Integrator' */
    {
      boolean_T lsat;
      boolean_T usat;
      lsat = ( ModelCopy2_X.Integrator_CSTATE_c4[0] <=
              (ModelCopy2_P.Integrator_LowerSat_l) );
      usat = ( ModelCopy2_X.Integrator_CSTATE_c4[0] >=
              ModelCopy2_P.Integrator_UpperSat_l );
      if ((!lsat && !usat) ||
          (lsat && (ModelCopy2_B.Gain[0] > 0)) ||
          (usat && (ModelCopy2_B.Gain[0] < 0)) ) {
        ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
          ->Integrator_CSTATE_c4[0] = ModelCopy2_B.Gain[0];
      } else {
        /* in saturation */
        ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
          ->Integrator_CSTATE_c4[0] = 0.0;
      }

      lsat = ( ModelCopy2_X.Integrator_CSTATE_c4[1] <=
              (ModelCopy2_P.Integrator_LowerSat_l) );
      usat = ( ModelCopy2_X.Integrator_CSTATE_c4[1] >=
              ModelCopy2_P.Integrator_UpperSat_l );
      if ((!lsat && !usat) ||
          (lsat && (ModelCopy2_B.Gain[1] > 0)) ||
          (usat && (ModelCopy2_B.Gain[1] < 0)) ) {
        ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
          ->Integrator_CSTATE_c4[1] = ModelCopy2_B.Gain[1];
      } else {
        /* in saturation */
        ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
          ->Integrator_CSTATE_c4[1] = 0.0;
      }
    }
  } else {
    {
      real_T *dx;
      int_T i;
      dx = &(((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
             ->Integrator_CSTATE_c4[0]);
      for (i=0; i < 12; i++) {
        dx[i] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S31>/wind_dfig_rotor' */

  /* Derivatives for Enabled SubSystem: '<S31>/wind_dfig_grid' */
  if (ModelCopy2_DW.wind_dfig_grid_MODE) {
    /* Derivatives for Integrator: '<S57>/Integrator' */
    {
      boolean_T lsat;
      boolean_T usat;
      lsat = ( ModelCopy2_X.Integrator_CSTATE_e1 <=
              (ModelCopy2_P.Integrator_LowerSat) );
      usat = ( ModelCopy2_X.Integrator_CSTATE_e1 >=
              ModelCopy2_P.Integrator_UpperSat );
      if ((!lsat && !usat) ||
          (lsat && (ModelCopy2_B.Gain_k > 0)) ||
          (usat && (ModelCopy2_B.Gain_k < 0)) ) {
        ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
          ->Integrator_CSTATE_e1 = ModelCopy2_B.Gain_k;
      } else {
        /* in saturation */
        ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
          ->Integrator_CSTATE_e1 = 0.0;
      }
    }

    /* Derivatives for Integrator: '<S56>/Integrator' */
    {
      boolean_T lsat;
      boolean_T usat;
      lsat = ( ModelCopy2_X.Integrator_CSTATE_g[0] <=
              (ModelCopy2_P.Integrator_LowerSat_d) );
      usat = ( ModelCopy2_X.Integrator_CSTATE_g[0] >=
              ModelCopy2_P.Integrator_UpperSat_d );
      if ((!lsat && !usat) ||
          (lsat && (ModelCopy2_B.Gain_iy[0] > 0)) ||
          (usat && (ModelCopy2_B.Gain_iy[0] < 0)) ) {
        ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
          ->Integrator_CSTATE_g[0] = ModelCopy2_B.Gain_iy[0];
      } else {
        /* in saturation */
        ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
          ->Integrator_CSTATE_g[0] = 0.0;
      }

      lsat = ( ModelCopy2_X.Integrator_CSTATE_g[1] <=
              (ModelCopy2_P.Integrator_LowerSat_d) );
      usat = ( ModelCopy2_X.Integrator_CSTATE_g[1] >=
              ModelCopy2_P.Integrator_UpperSat_d );
      if ((!lsat && !usat) ||
          (lsat && (ModelCopy2_B.Gain_iy[1] > 0)) ||
          (usat && (ModelCopy2_B.Gain_iy[1] < 0)) ) {
        ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
          ->Integrator_CSTATE_g[1] = ModelCopy2_B.Gain_iy[1];
      } else {
        /* in saturation */
        ((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
          ->Integrator_CSTATE_g[1] = 0.0;
      }
    }
  } else {
    {
      real_T *dx;
      int_T i;
      dx = &(((XDot_ModelCopy2_T *) ModelCopy2_M->ModelData.derivs)
             ->Integrator_CSTATE_e1);
      for (i=0; i < 3; i++) {
        dx[i] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S31>/wind_dfig_grid' */
}

/* Model initialize function */
void ModelCopy2_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  ModelCopy2_P.Integrator_UpperSat_d = rtInf;
  ModelCopy2_P.Integrator_LowerSat_d = rtMinusInf;
  ModelCopy2_P.Saturation_UpperSat_m = rtInf;
  ModelCopy2_P.Saturation_LowerSat_c = rtMinusInf;
  ModelCopy2_P.inf_UpperSat = rtInf;
  ModelCopy2_P.inf_UpperSat_o = rtInf;
  ModelCopy2_P.Integrator_UpperSat_l = rtInf;
  ModelCopy2_P.Integrator_LowerSat_l = rtMinusInf;
  ModelCopy2_P.Saturation_UpperSat_k = rtInf;
  ModelCopy2_P.Saturation_LowerSat_p = rtMinusInf;
  ModelCopy2_P.Saturation1_UpperSat = rtInf;
  ModelCopy2_P.inf_UpperSat_d = rtInf;

  /* initialize real-time model */
  (void) memset((void *)ModelCopy2_M, 0,
                sizeof(RT_MODEL_ModelCopy2_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&ModelCopy2_M->solverInfo,
                          &ModelCopy2_M->Timing.simTimeStep);
    rtsiSetTPtr(&ModelCopy2_M->solverInfo, &rtmGetTPtr(ModelCopy2_M));
    rtsiSetStepSizePtr(&ModelCopy2_M->solverInfo,
                       &ModelCopy2_M->Timing.stepSize0);
    rtsiSetdXPtr(&ModelCopy2_M->solverInfo, &ModelCopy2_M->ModelData.derivs);
    rtsiSetContStatesPtr(&ModelCopy2_M->solverInfo, (real_T **)
                         &ModelCopy2_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&ModelCopy2_M->solverInfo,
      &ModelCopy2_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&ModelCopy2_M->solverInfo, (&rtmGetErrorStatus
      (ModelCopy2_M)));
    rtsiSetRTModelPtr(&ModelCopy2_M->solverInfo, ModelCopy2_M);
  }

  rtsiSetSimTimeStep(&ModelCopy2_M->solverInfo, MAJOR_TIME_STEP);
  ModelCopy2_M->ModelData.intgData.y = ModelCopy2_M->ModelData.odeY;
  ModelCopy2_M->ModelData.intgData.f[0] = ModelCopy2_M->ModelData.odeF[0];
  ModelCopy2_M->ModelData.intgData.f[1] = ModelCopy2_M->ModelData.odeF[1];
  ModelCopy2_M->ModelData.intgData.f[2] = ModelCopy2_M->ModelData.odeF[2];
  ModelCopy2_M->ModelData.contStates = ((X_ModelCopy2_T *) &ModelCopy2_X);
  rtsiSetSolverData(&ModelCopy2_M->solverInfo, (void *)
                    &ModelCopy2_M->ModelData.intgData);
  rtsiSetSolverName(&ModelCopy2_M->solverInfo,"ode3");
  ModelCopy2_M->solverInfoPtr = (&ModelCopy2_M->solverInfo);

  /* Initialize timing info */
  {
    int_T *mdlTsMap = ModelCopy2_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    ModelCopy2_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    ModelCopy2_M->Timing.sampleTimes = (&ModelCopy2_M->Timing.sampleTimesArray[0]);
    ModelCopy2_M->Timing.offsetTimes = (&ModelCopy2_M->Timing.offsetTimesArray[0]);

    /* task periods */
    ModelCopy2_M->Timing.sampleTimes[0] = (0.0);
    ModelCopy2_M->Timing.sampleTimes[1] = (0.001);

    /* task offsets */
    ModelCopy2_M->Timing.offsetTimes[0] = (0.0);
    ModelCopy2_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(ModelCopy2_M, &ModelCopy2_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = ModelCopy2_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    ModelCopy2_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(ModelCopy2_M, -1);
  ModelCopy2_M->Timing.stepSize0 = 0.001;
  rtmSetFirstInitCond(ModelCopy2_M, 1);
  ModelCopy2_M->solverInfoPtr = (&ModelCopy2_M->solverInfo);
  ModelCopy2_M->Timing.stepSize = (0.001);
  rtsiSetFixedStepSize(&ModelCopy2_M->solverInfo, 0.001);
  rtsiSetSolverMode(&ModelCopy2_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  (void) memset(((void *) &ModelCopy2_B), 0,
                sizeof(B_ModelCopy2_T));

  /* states (continuous) */
  {
    (void) memset((void *)&ModelCopy2_X, 0,
                  sizeof(X_ModelCopy2_T));
  }

  /* states (dwork) */
  (void) memset((void *)&ModelCopy2_DW, 0,
                sizeof(DW_ModelCopy2_T));

  /* child S-Function registration */
  {
    RTWSfcnInfo *sfcnInfo = &ModelCopy2_M->NonInlinedSFcns.sfcnInfo;
    ModelCopy2_M->sfcnInfo = (sfcnInfo);
    rtssSetErrorStatusPtr(sfcnInfo, (&rtmGetErrorStatus(ModelCopy2_M)));
    rtssSetNumRootSampTimesPtr(sfcnInfo, &ModelCopy2_M->Sizes.numSampTimes);
    ModelCopy2_M->NonInlinedSFcns.taskTimePtrs[0] = &(rtmGetTPtr(ModelCopy2_M)[0]);
    ModelCopy2_M->NonInlinedSFcns.taskTimePtrs[1] = &(rtmGetTPtr(ModelCopy2_M)[1]);
    rtssSetTPtrPtr(sfcnInfo,ModelCopy2_M->NonInlinedSFcns.taskTimePtrs);
    rtssSetTStartPtr(sfcnInfo, &rtmGetTStart(ModelCopy2_M));
    rtssSetTFinalPtr(sfcnInfo, &rtmGetTFinal(ModelCopy2_M));
    rtssSetTimeOfLastOutputPtr(sfcnInfo, &rtmGetTimeOfLastOutput(ModelCopy2_M));
    rtssSetStepSizePtr(sfcnInfo, &ModelCopy2_M->Timing.stepSize);
    rtssSetStopRequestedPtr(sfcnInfo, &rtmGetStopRequested(ModelCopy2_M));
    rtssSetDerivCacheNeedsResetPtr(sfcnInfo,
      &ModelCopy2_M->ModelData.derivCacheNeedsReset);
    rtssSetZCCacheNeedsResetPtr(sfcnInfo,
      &ModelCopy2_M->ModelData.zCCacheNeedsReset);
    rtssSetBlkStateChangePtr(sfcnInfo, &ModelCopy2_M->ModelData.blkStateChange);
    rtssSetSampleHitsPtr(sfcnInfo, &ModelCopy2_M->Timing.sampleHits);
    rtssSetPerTaskSampleHitsPtr(sfcnInfo,
      &ModelCopy2_M->Timing.perTaskSampleHits);
    rtssSetSimModePtr(sfcnInfo, &ModelCopy2_M->simMode);
    rtssSetSolverInfoPtr(sfcnInfo, &ModelCopy2_M->solverInfoPtr);
  }

  ModelCopy2_M->Sizes.numSFcns = (1);

  /* register each child */
  {
    (void) memset((void *)&ModelCopy2_M->NonInlinedSFcns.childSFunctions[0], 0,
                  1*sizeof(SimStruct));
    ModelCopy2_M->childSfunctions =
      (&ModelCopy2_M->NonInlinedSFcns.childSFunctionPtrs[0]);
    ModelCopy2_M->childSfunctions[0] =
      (&ModelCopy2_M->NonInlinedSFcns.childSFunctions[0]);

    /* Level2 S-Function Block: ModelCopy2/<S92>/State-Space (sfun_psbdqc) */
    {
      SimStruct *rts = ModelCopy2_M->childSfunctions[0];

      /* timing info */
      time_T *sfcnPeriod = ModelCopy2_M->NonInlinedSFcns.Sfcn0.sfcnPeriod;
      time_T *sfcnOffset = ModelCopy2_M->NonInlinedSFcns.Sfcn0.sfcnOffset;
      int_T *sfcnTsMap = ModelCopy2_M->NonInlinedSFcns.Sfcn0.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts, &ModelCopy2_M->NonInlinedSFcns.blkInfo2[0]);
      }

      ssSetRTWSfcnInfo(rts, ModelCopy2_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &ModelCopy2_M->NonInlinedSFcns.methods2[0]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &ModelCopy2_M->NonInlinedSFcns.methods3[0]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &ModelCopy2_M->NonInlinedSFcns.statesInfo2[0]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &ModelCopy2_M->NonInlinedSFcns.Sfcn0.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &ModelCopy2_M->NonInlinedSFcns.Sfcn0.UPtrs0;

          {
            int_T i1;
            const real_T *u0 = &ModelCopy2_B.ComplextoRealImag_o1[0];
            for (i1=0; i1 < 5; i1++) {
              sfcnUPtrs[i1] = &u0[i1];
            }

            u0 = &ModelCopy2_B.ComplextoRealImag_o2[0];
            for (i1=0; i1 < 5; i1++) {
              sfcnUPtrs[i1+ 5] = &u0[i1];
            }
          }

          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 10);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &ModelCopy2_M->NonInlinedSFcns.Sfcn0.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 6);
          ssSetOutputPortSignal(rts, 0, ((real_T *) ModelCopy2_B.StateSpace));
        }
      }

      /* path info */
      ssSetModelName(rts, "State-Space");
      ssSetPath(rts, "ModelCopy2/powergui/EquivalentModel1/State-Space");
      ssSetRTModel(rts,ModelCopy2_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &ModelCopy2_M->NonInlinedSFcns.Sfcn0.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)ModelCopy2_P.StateSpace_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)ModelCopy2_P.StateSpace_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)ModelCopy2_P.StateSpace_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)ModelCopy2_P.StateSpace_P4_Size);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *) &ModelCopy2_DW.StateSpace_RWORK[0]);
      ssSetIWork(rts, (int_T *) &ModelCopy2_DW.StateSpace_IWORK[0]);
      ssSetPWork(rts, (void **) &ModelCopy2_DW.StateSpace_PWORK[0]);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &ModelCopy2_M->NonInlinedSFcns.Sfcn0.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &ModelCopy2_M->NonInlinedSFcns.Sfcn0.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 6);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &ModelCopy2_DW.StateSpace_RWORK[0]);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 4);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &ModelCopy2_DW.StateSpace_IWORK[0]);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 15);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &ModelCopy2_DW.StateSpace_PWORK[0]);
      }

      /* registration */
      sfun_psbdqc(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.0);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }
  }

  /* InitializeConditions for Enabled SubSystem: '<S43>/SET  Priority' */
  /* InitializeConditions for Memory: '<S45>/IC=ic' */
  ModelCopy2_DW.ICic_PreviousInput_n = ModelCopy2_P.Bistable1_ic;

  /* End of InitializeConditions for SubSystem: '<S43>/SET  Priority' */

  /* Start for Enabled SubSystem: '<S43>/SET  Priority' */
  /* VirtualOutportStart for Outport: '<S45>/Q' */
  ModelCopy2_B.LogicalOperator_e = ModelCopy2_P.Q_Y0_f;

  /* End of Start for SubSystem: '<S43>/SET  Priority' */

  /* InitializeConditions for Enabled SubSystem: '<S43>/RESET Priority' */
  /* InitializeConditions for Memory: '<S44>/IC=ic' */
  ModelCopy2_DW.ICic_PreviousInput_c = ModelCopy2_P.Bistable1_ic;

  /* End of InitializeConditions for SubSystem: '<S43>/RESET Priority' */

  /* Start for Enabled SubSystem: '<S43>/RESET Priority' */
  /* VirtualOutportStart for Outport: '<S44>/Q' */
  ModelCopy2_B.LogicalOperator2_i = ModelCopy2_P.Q_Y0;

  /* End of Start for SubSystem: '<S43>/RESET Priority' */
  /* Level2 S-Function Block: '<S92>/State-Space' (sfun_psbdqc) */
  {
    SimStruct *rts = ModelCopy2_M->childSfunctions[0];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Start for Enabled SubSystem: '<S31>/wind_dfig_rotor' */

  /* Start for Enabled SubSystem: '<S49>/V Regulator' */

  /* Start for TransportDelay: '<S84>/Transport Delay' */
  {
    real_T *pBuffer = &ModelCopy2_DW.TransportDelay_RWORK.TUbufferArea[0];
    ModelCopy2_DW.TransportDelay_IWORK.Tail = 0;
    ModelCopy2_DW.TransportDelay_IWORK.Head = 0;
    ModelCopy2_DW.TransportDelay_IWORK.Last = 0;
    ModelCopy2_DW.TransportDelay_IWORK.CircularBufSize = 1024;
    pBuffer[0] = ModelCopy2_P.TransportDelay_InitOutput_i;
    pBuffer[1024] = ModelCopy2_M->Timing.t[0];
    ModelCopy2_DW.TransportDelay_PWORK.TUbufferPtrs[0] = (void *) &pBuffer[0];
    ModelCopy2_DW.TransportDelay_PWORK.TUbufferPtrs[1] = (void *) &pBuffer[1024];
  }

  /* Start for TransportDelay: '<S85>/Transport Delay' */
  {
    real_T *pBuffer = &ModelCopy2_DW.TransportDelay_RWORK_o.TUbufferArea[0];
    ModelCopy2_DW.TransportDelay_IWORK_g.Tail = 0;
    ModelCopy2_DW.TransportDelay_IWORK_g.Head = 0;
    ModelCopy2_DW.TransportDelay_IWORK_g.Last = 0;
    ModelCopy2_DW.TransportDelay_IWORK_g.CircularBufSize = 1024;
    pBuffer[0] = ModelCopy2_P.TransportDelay_InitOutput_a;
    pBuffer[1024] = ModelCopy2_M->Timing.t[0];
    ModelCopy2_DW.TransportDelay_PWORK_g.TUbufferPtrs[0] = (void *) &pBuffer[0];
    ModelCopy2_DW.TransportDelay_PWORK_g.TUbufferPtrs[1] = (void *) &pBuffer
      [1024];
  }

  /* End of Start for SubSystem: '<S49>/V Regulator' */

  /* InitializeConditions for Enabled SubSystem: '<S49>/V Regulator' */
  /* InitializeConditions for Integrator: '<S84>/integrator' */
  ModelCopy2_X.integrator_CSTATE = ModelCopy2_P.integrator_IC_a;

  /* InitializeConditions for RateLimiter: '<S67>/Rate Limiter ' */
  ModelCopy2_DW.LastMajorTime_kl = (rtInf);

  /* InitializeConditions for Integrator: '<S85>/integrator' */
  ModelCopy2_X.integrator_CSTATE_p = ModelCopy2_P.integrator_IC_d;

  /* InitializeConditions for Integrator: '<S86>/Integrator' */
  if (rtmIsFirstInitCond(ModelCopy2_M)) {
    ModelCopy2_X.Integrator_CSTATE_p = 0.0;
  }

  ModelCopy2_DW.Integrator_IWORK.IcNeedsLoading = 1;

  /* InitializeConditions for RateLimiter: '<S67>/Rate Limiter' */
  ModelCopy2_DW.LastMajorTime_a = (rtInf);

  /* End of InitializeConditions for SubSystem: '<S49>/V Regulator' */

  /* Start for Enabled SubSystem: '<S49>/V Regulator' */
  /* VirtualOutportStart for Outport: '<S67>/Idr+' */
  ModelCopy2_B.RateLimiter_g = ModelCopy2_P.Idr_Y0_m;

  /* End of Start for SubSystem: '<S49>/V Regulator' */

  /* Start for Enabled SubSystem: '<S49>/Q Regulator' */

  /* Start for TransportDelay: '<S82>/Transport Delay' */
  {
    real_T *pBuffer = &ModelCopy2_DW.TransportDelay_RWORK_e.TUbufferArea[0];
    ModelCopy2_DW.TransportDelay_IWORK_h.Tail = 0;
    ModelCopy2_DW.TransportDelay_IWORK_h.Head = 0;
    ModelCopy2_DW.TransportDelay_IWORK_h.Last = 0;
    ModelCopy2_DW.TransportDelay_IWORK_h.CircularBufSize = 1024;
    pBuffer[0] = ModelCopy2_P.TransportDelay_InitOutput_d;
    pBuffer[1024] = ModelCopy2_M->Timing.t[0];
    ModelCopy2_DW.TransportDelay_PWORK_a.TUbufferPtrs[0] = (void *) &pBuffer[0];
    ModelCopy2_DW.TransportDelay_PWORK_a.TUbufferPtrs[1] = (void *) &pBuffer
      [1024];
  }

  /* End of Start for SubSystem: '<S49>/Q Regulator' */

  /* InitializeConditions for Enabled SubSystem: '<S49>/Q Regulator' */
  /* InitializeConditions for Integrator: '<S82>/integrator' */
  ModelCopy2_X.integrator_CSTATE_f = ModelCopy2_P.integrator_IC_k;

  /* InitializeConditions for RateLimiter: '<S66>/Rate Limiter ' */
  ModelCopy2_DW.LastMajorTime_m = (rtInf);

  /* InitializeConditions for Integrator: '<S83>/Integrator' */
  if (rtmIsFirstInitCond(ModelCopy2_M)) {
    ModelCopy2_X.Integrator_CSTATE_a = 0.0;
  }

  ModelCopy2_DW.Integrator_IWORK_o.IcNeedsLoading = 1;

  /* InitializeConditions for RateLimiter: '<S66>/Rate Limiter' */
  ModelCopy2_DW.LastMajorTime_m5 = (rtInf);

  /* End of InitializeConditions for SubSystem: '<S49>/Q Regulator' */

  /* Start for Enabled SubSystem: '<S49>/Q Regulator' */
  /* VirtualOutportStart for Outport: '<S66>/Idr+' */
  ModelCopy2_B.RateLimiter_j = ModelCopy2_P.Idr_Y0;

  /* End of Start for SubSystem: '<S49>/Q Regulator' */

  /* Start for Enabled SubSystem: '<S64>/Subsystem' */

  /* Start for TransportDelay: '<S78>/Transport Delay' */
  {
    real_T *pBuffer = &ModelCopy2_DW.TransportDelay_RWORK_k.TUbufferArea[0];
    ModelCopy2_DW.TransportDelay_IWORK_ao.Tail = 0;
    ModelCopy2_DW.TransportDelay_IWORK_ao.Head = 0;
    ModelCopy2_DW.TransportDelay_IWORK_ao.Last = 0;
    ModelCopy2_DW.TransportDelay_IWORK_ao.CircularBufSize = 1024;
    pBuffer[0] = ModelCopy2_P.TransportDelay_InitOutput;
    pBuffer[1024] = ModelCopy2_M->Timing.t[0];
    ModelCopy2_DW.TransportDelay_PWORK_gm.TUbufferPtrs[0] = (void *) &pBuffer[0];
    ModelCopy2_DW.TransportDelay_PWORK_gm.TUbufferPtrs[1] = (void *) &pBuffer
      [1024];
  }

  /* Start for TransportDelay: '<S74>/Transport Delay' */
  {
    real_T *pBuffer = &ModelCopy2_DW.TransportDelay_RWORK_f.TUbufferArea[0];
    ModelCopy2_DW.TransportDelay_IWORK_ak.Tail = 0;
    ModelCopy2_DW.TransportDelay_IWORK_ak.Head = 0;
    ModelCopy2_DW.TransportDelay_IWORK_ak.Last = 0;
    ModelCopy2_DW.TransportDelay_IWORK_ak.CircularBufSize = 1024;
    pBuffer[0] = ModelCopy2_P.TransportDelay_InitOutput_e;
    pBuffer[1024] = ModelCopy2_M->Timing.t[0];
    ModelCopy2_DW.TransportDelay_PWORK_n.TUbufferPtrs[0] = (void *) &pBuffer[0];
    ModelCopy2_DW.TransportDelay_PWORK_n.TUbufferPtrs[1] = (void *) &pBuffer
      [1024];
  }

  /* End of Start for SubSystem: '<S64>/Subsystem' */

  /* InitializeConditions for Enabled SubSystem: '<S64>/Subsystem' */
  /* InitializeConditions for RateLimiter: '<S71>/Rate Limiter ' */
  ModelCopy2_DW.LastMajorTime_i5 = (rtInf);

  /* InitializeConditions for Integrator: '<S78>/integrator' */
  ModelCopy2_X.integrator_CSTATE_k = ModelCopy2_P.integrator_IC;

  /* InitializeConditions for Integrator: '<S74>/integrator' */
  ModelCopy2_X.integrator_CSTATE_o = ModelCopy2_P.integrator_IC_l;

  /* InitializeConditions for Integrator: '<S77>/Integrator' */
  if (rtmIsFirstInitCond(ModelCopy2_M)) {
    ModelCopy2_X.Integrator_CSTATE_k = 0.0;
  }

  ModelCopy2_DW.Integrator_IWORK_a.IcNeedsLoading = 1;

  /* InitializeConditions for RateLimiter: '<S71>/Rate Limiter' */
  ModelCopy2_DW.LastMajorTime_l = (rtInf);

  /* End of InitializeConditions for SubSystem: '<S64>/Subsystem' */

  /* Start for Enabled SubSystem: '<S64>/Subsystem' */
  /* VirtualOutportStart for Outport: '<S71>/Iqr+' */
  ModelCopy2_B.RateLimiter_m = ModelCopy2_P.Iqr_Y0;

  /* End of Start for SubSystem: '<S64>/Subsystem' */

  /* Start for Enabled SubSystem: '<S64>/Subsystem ' */

  /* Start for TransportDelay: '<S79>/Transport Delay' */
  {
    real_T *pBuffer = &ModelCopy2_DW.TransportDelay_RWORK_i.TUbufferArea[0];
    ModelCopy2_DW.TransportDelay_IWORK_a.Tail = 0;
    ModelCopy2_DW.TransportDelay_IWORK_a.Head = 0;
    ModelCopy2_DW.TransportDelay_IWORK_a.Last = 0;
    ModelCopy2_DW.TransportDelay_IWORK_a.CircularBufSize = 1024;
    pBuffer[0] = ModelCopy2_P.TransportDelay_InitOutput_c;
    pBuffer[1024] = ModelCopy2_M->Timing.t[0];
    ModelCopy2_DW.TransportDelay_PWORK_c.TUbufferPtrs[0] = (void *) &pBuffer[0];
    ModelCopy2_DW.TransportDelay_PWORK_c.TUbufferPtrs[1] = (void *) &pBuffer
      [1024];
  }

  /* End of Start for SubSystem: '<S64>/Subsystem ' */

  /* InitializeConditions for Enabled SubSystem: '<S64>/Subsystem ' */
  /* InitializeConditions for Integrator: '<S79>/integrator' */
  ModelCopy2_X.integrator_CSTATE_d = ModelCopy2_P.integrator_IC_o;

  /* InitializeConditions for Integrator: '<S80>/Integrator' */
  if (rtmIsFirstInitCond(ModelCopy2_M)) {
    ModelCopy2_X.Integrator_CSTATE_e = 0.0;
  }

  ModelCopy2_DW.Integrator_IWORK_l.IcNeedsLoading = 1;

  /* InitializeConditions for RateLimiter: '<S72>/Rate Limiter' */
  ModelCopy2_DW.LastMajorTime_c = (rtInf);

  /* End of InitializeConditions for SubSystem: '<S64>/Subsystem ' */

  /* Start for Enabled SubSystem: '<S64>/Subsystem ' */
  /* VirtualOutportStart for Outport: '<S72>/Iqr+' */
  ModelCopy2_B.RateLimiter_h = ModelCopy2_P.Iqr_Y0_b;

  /* End of Start for SubSystem: '<S64>/Subsystem ' */
  /* End of Start for SubSystem: '<S31>/wind_dfig_rotor' */

  /* InitializeConditions for Enabled SubSystem: '<S31>/wind_dfig_rotor' */
  /* InitializeConditions for Memory: '<S63>/IC=ic' */
  ModelCopy2_DW.ICic_PreviousInput_b = ModelCopy2_P.ICic_X0;

  /* InitializeConditions for RateLimiter: '<S64>/Rate Limiter' */
  ModelCopy2_DW.LastMajorTime_k = (rtInf);

  /* InitializeConditions for RateLimiter: '<S65>/Rate Limiter' */
  ModelCopy2_DW.LastMajorTime_n = (rtInf);

  /* InitializeConditions for Integrator: '<S69>/Integrator' */
  ModelCopy2_X.Integrator_CSTATE_c4[0] = ModelCopy2_P.Subsystem1_Init_d;
  ModelCopy2_X.Integrator_CSTATE_c4[1] = ModelCopy2_P.Subsystem1_Init_d;

  /* End of InitializeConditions for SubSystem: '<S31>/wind_dfig_rotor' */

  /* Start for Enabled SubSystem: '<S31>/wind_dfig_rotor' */
  /* VirtualOutportStart for Outport: '<S49>/Vdq_ctrl_rotor_conv' */
  ModelCopy2_B.MagnitudeAngletoComplex_i.re =
    ModelCopy2_P.Vdq_ctrl_rotor_conv_Y0;
  ModelCopy2_B.MagnitudeAngletoComplex_i.im = 0.0;

  /* End of Start for SubSystem: '<S31>/wind_dfig_rotor' */

  /* InitializeConditions for Enabled SubSystem: '<S31>/wind_dfig_grid' */
  /* InitializeConditions for Integrator: '<S57>/Integrator' */
  ModelCopy2_X.Integrator_CSTATE_e1 = ModelCopy2_P.Subsystem1_Init;

  /* InitializeConditions for RateLimiter: '<S51>/Rate Limiter' */
  ModelCopy2_DW.LastMajorTime_h = (rtInf);

  /* InitializeConditions for RateLimiter: '<S48>/Rate Limiter' */
  ModelCopy2_DW.LastMajorTime_ml = (rtInf);

  /* InitializeConditions for RateLimiter: '<S53>/Rate Limiter' */
  ModelCopy2_DW.LastMajorTime_ck = (rtInf);

  /* InitializeConditions for Integrator: '<S56>/Integrator' */
  ModelCopy2_X.Integrator_CSTATE_g[0] = ModelCopy2_P.Subsystem3_Init;
  ModelCopy2_X.Integrator_CSTATE_g[1] = ModelCopy2_P.Subsystem3_Init;

  /* End of InitializeConditions for SubSystem: '<S31>/wind_dfig_grid' */

  /* Start for Enabled SubSystem: '<S31>/wind_dfig_grid' */
  /* VirtualOutportStart for Outport: '<S48>/Vdq_ctrl_grid_conv' */
  ModelCopy2_B.MagnitudeAngletoComplex_c.re = ModelCopy2_P.Vdq_ctrl_grid_conv_Y0;
  ModelCopy2_B.MagnitudeAngletoComplex_c.im = 0.0;

  /* End of Start for SubSystem: '<S31>/wind_dfig_grid' */

  /* InitializeConditions for Integrator: '<S41>/phiqr' */
  ModelCopy2_X.phiqr_CSTATE = ModelCopy2_P.phiqr_IC;

  /* InitializeConditions for Integrator: '<S41>/phidr' */
  ModelCopy2_X.phidr_CSTATE = ModelCopy2_P.phidr_IC;

  /* InitializeConditions for Integrator: '<S42>/phiqs' */
  ModelCopy2_X.phiqs_CSTATE = ModelCopy2_P.phiqs_IC;

  /* InitializeConditions for Integrator: '<S42>/phids' */
  ModelCopy2_X.phids_CSTATE = ModelCopy2_P.phids_IC;

  /* InitializeConditions for Integrator: '<S89>/Integrator' */
  ModelCopy2_X.Integrator_CSTATE[0] =
    ModelCopy2_P.dqaxismodelofa3phaseseriesRLbra[0];
  ModelCopy2_X.Integrator_CSTATE[1] =
    ModelCopy2_P.dqaxismodelofa3phaseseriesRLbra[1];

  /* InitializeConditions for Integrator: '<S10>/Integrator' */
  ModelCopy2_X.Integrator_CSTATE_j = ModelCopy2_P.Integrator_IC;

  /* InitializeConditions for Memory: '<S10>/Memory' */
  ModelCopy2_DW.Memory_PreviousInput = ModelCopy2_P.Memory_X0;

  /* InitializeConditions for Integrator: '<S35>/Integrator' */
  ModelCopy2_X.Integrator_CSTATE_c = ModelCopy2_P.Integrator_IC_m;

  /* InitializeConditions for RateLimiter: '<S31>/Rate Limiter   ' */
  ModelCopy2_DW.LastMajorTime = (rtInf);

  /* InitializeConditions for Integrator: '<S35>/Integrator1' */
  ModelCopy2_X.Integrator1_CSTATE = ModelCopy2_P.Integrator1_IC;

  /* InitializeConditions for Integrator: '<S32>/Integrator' */
  ModelCopy2_X.Integrator_CSTATE_cb = ModelCopy2_P.DCbusmodel_Vdc_Init;

  /* InitializeConditions for Memory: '<S47>/IC=ic' */
  ModelCopy2_DW.ICic_PreviousInput = ModelCopy2_P.ICic_X0_f;

  /* InitializeConditions for RateLimiter: '<S31>/Rate Limiter ' */
  ModelCopy2_DW.LastMajorTime_i = (rtInf);

  /* set "at time zero" to false */
  if (rtmIsFirstInitCond(ModelCopy2_M)) {
    rtmSetFirstInitCond(ModelCopy2_M, 0);
  }
}

/* Model terminate function */
void ModelCopy2_terminate(void)
{
  /* Level2 S-Function Block: '<S92>/State-Space' (sfun_psbdqc) */
  {
    SimStruct *rts = ModelCopy2_M->childSfunctions[0];
    sfcnTerminate(rts);
  }
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
