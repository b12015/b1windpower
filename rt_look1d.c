/*
 * File: rt_look1d.c
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

#include "rt_look1d.h"

/* 1D lookup routine for data type of real_T. */
real_T rt_Lookup(const real_T *x, int_T xlen, real_T u, const real_T *y)
{
  int_T idx = rt_GetLookupIndex(x, xlen, u);
  real_T num = y[idx+1] - y[idx];
  real_T den = x[idx+1] - x[idx];

  /* Due to the way the binary search is implemented
     in rt_look.c (rt_GetLookupIndex), den cannot be
     0.  Equivalently, m cannot be inf or nan. */
  real_T m = num/den;
  return (y[idx] + (m * (u - x[idx])));
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
