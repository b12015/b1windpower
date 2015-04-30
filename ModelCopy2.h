/*
 * File: ModelCopy2.h
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

#ifndef RTW_HEADER_ModelCopy2_h_
#define RTW_HEADER_ModelCopy2_h_
#include <stddef.h>
#include <float.h>
#include <math.h>
#include <string.h>
#ifndef ModelCopy2_COMMON_INCLUDES_
# define ModelCopy2_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#endif                                 /* ModelCopy2_COMMON_INCLUDES_ */
#include "sfcn_bridge.h"
#include "ModelCopy2_types.h"
#include "rtGetInf.h"
#include "rt_nonfinite.h"
#include "rt_defines.h"
#include "rt_look.h"
#include "rt_look1d.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetFinalTime
# define rtmGetFinalTime(rtm)          ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetSampleHitArray
# define rtmGetSampleHitArray(rtm)     ((rtm)->Timing.sampleHitArray)
#endif

#ifndef rtmGetStepSize
# define rtmGetStepSize(rtm)           ((rtm)->Timing.stepSize)
#endif

#ifndef rtmGet_TimeOfLastOutput
# define rtmGet_TimeOfLastOutput(rtm)  ((rtm)->Timing.timeOfLastOutput)
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
# define rtmGetTFinal(rtm)             ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTStart
# define rtmGetTStart(rtm)             ((rtm)->Timing.tStart)
#endif

#ifndef rtmGetTimeOfLastOutput
# define rtmGetTimeOfLastOutput(rtm)   ((rtm)->Timing.timeOfLastOutput)
#endif

/* Block signals (auto storage) */
typedef struct {
  creal_T MagnitudeAngletoComplex[3];  /* '<S7>/Magnitude-Angle to Complex' */
  creal_T MagnitudeAngletoComplex_l[3];/* '<S8>/Magnitude-Angle to Complex' */
  creal_T MagnitudeAngletoComplex_i;   /* '<S68>/Magnitude-Angle to Complex' */
  creal_T MagnitudeAngletoComplex_c;   /* '<S54>/Magnitude-Angle to Complex' */
  real_T Gain1;                        /* '<S10>/Gain1' */
  real_T Memory;                       /* '<S10>/Memory' */
  real_T Switch2;                      /* '<S10>/Switch2' */
  real_T Gain3;                        /* '<S6>/Gain3' */
  real_T ComplextoRealImag_o1[5];      /* '<S96>/Complex to Real-Imag' */
  real_T ComplextoRealImag_o2[5];      /* '<S96>/Complex to Real-Imag' */
  real_T StateSpace[6];                /* '<S92>/State-Space' */
  real_T Ppu;                          /* '<S33>/-1' */
  real_T Qpu;                          /* '<S33>/-2' */
  real_T QUT;                          /* '<Root>/900k_ ' */
  real_T RateLimiter;                  /* '<S31>/Rate Limiter   ' */
  real_T _2H;                          /* '<S35>/1_2H' */
  real_T web;                          /* '<S35>/web' */
  real_T ComplextoMagnitudeAngle_o1;   /* '<S47>/Complex to Magnitude-Angle' */
  real_T ICic;                         /* '<S47>/IC=ic' */
  real_T u;                            /* '<S47>/20%' */
  real_T RateLimiter_b;                /* '<S31>/Rate Limiter ' */
  real_T Switch;                       /* '<S41>/Switch' */
  real_T Switch3;                      /* '<S41>/Switch3' */
  real_T Switch_n;                     /* '<S42>/Switch' */
  real_T Switch2_d;                    /* '<S42>/Switch2' */
  real_T degrad1;                      /* '<S32>/deg->rad1' */
  real_T Switch_nz[2];                 /* '<S89>/Switch' */
  real_T Switch1;                      /* '<S10>/Switch1' */
  real_T Switch_j;                     /* '<S10>/Switch' */
  real_T ICic_a;                       /* '<S63>/IC=ic' */
  real_T u_g;                          /* '<S63>/20%' */
  real_T ComplextoRealImag_o1_b;       /* '<S62>/Complex to Real-Imag' */
  real_T ComplextoRealImag_o2_p;       /* '<S62>/Complex to Real-Imag' */
  real_T RateLimiter_e;                /* '<S64>/Rate Limiter' */
  real_T RateLimiter_n[2];             /* '<S65>/Rate Limiter' */
  real_T Gain[2];                      /* '<S69>/Gain' */
  real_T Droop;                        /* '<S67>/Droop' */
  real_T RateLimiter_b5;               /* '<S67>/Rate Limiter ' */
  real_T RateLimiter_g;                /* '<S67>/Rate Limiter' */
  real_T Gain_i;                       /* '<S86>/Gain' */
  real_T RateLimiter_nz;               /* '<S66>/Rate Limiter ' */
  real_T RateLimiter_j;                /* '<S66>/Rate Limiter' */
  real_T Gain_j;                       /* '<S83>/Gain' */
  real_T RateLimiter_h;                /* '<S72>/Rate Limiter' */
  real_T Gain_h;                       /* '<S80>/Gain' */
  real_T Product;                      /* '<S76>/Product' */
  real_T Product2;                     /* '<S75>/Product2' */
  real_T RateLimiter_e2;               /* '<S71>/Rate Limiter ' */
  real_T Looses;                       /* '<S73>/Sum2' */
  real_T RateLimiter_m;                /* '<S71>/Rate Limiter' */
  real_T Gain_b;                       /* '<S77>/Gain' */
  real_T RateLimiter_ek;               /* '<S51>/Rate Limiter' */
  real_T RateLimiter_l;                /* '<S48>/Rate Limiter' */
  real_T RateLimiter_jx[2];            /* '<S53>/Rate Limiter' */
  real_T Gain_iy[2];                   /* '<S56>/Gain' */
  real_T Gain_k;                       /* '<S57>/Gain' */
  boolean_T RelationalOperator1;       /* '<S43>/Relational Operator1' */
  boolean_T Amplitude;                 /* '<S6>/Relational Operator' */
  boolean_T LogicalOperator1;          /* '<S6>/Logical Operator1' */
  boolean_T RelationalOperator1_b;     /* '<S10>/Relational Operator1' */
  boolean_T Phase;                     /* '<S6>/Relational Operator1' */
  boolean_T DataTypeConversion2;       /* '<S6>/Data Type  Conversion2' */
  boolean_T ICic_c;                    /* '<S45>/IC=ic' */
  boolean_T LogicalOperator_e;         /* '<S45>/Logical Operator' */
  boolean_T ICic_cd;                   /* '<S44>/IC=ic' */
  boolean_T LogicalOperator2_i;        /* '<S44>/Logical Operator2' */
} B_ModelCopy2_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T Memory_PreviousInput;         /* '<S10>/Memory' */
  real_T PrevY;                        /* '<S31>/Rate Limiter   ' */
  real_T LastMajorTime;                /* '<S31>/Rate Limiter   ' */
  real_T ICic_PreviousInput;           /* '<S47>/IC=ic' */
  real_T PrevY_k;                      /* '<S31>/Rate Limiter ' */
  real_T LastMajorTime_i;              /* '<S31>/Rate Limiter ' */
  real_T ICic_PreviousInput_b;         /* '<S63>/IC=ic' */
  real_T PrevY_h;                      /* '<S64>/Rate Limiter' */
  real_T LastMajorTime_k;              /* '<S64>/Rate Limiter' */
  real_T PrevY_g[2];                   /* '<S65>/Rate Limiter' */
  real_T LastMajorTime_n;              /* '<S65>/Rate Limiter' */
  real_T PrevY_i;                      /* '<S67>/Rate Limiter ' */
  real_T LastMajorTime_kl;             /* '<S67>/Rate Limiter ' */
  real_T PrevY_hi;                     /* '<S67>/Rate Limiter' */
  real_T LastMajorTime_a;              /* '<S67>/Rate Limiter' */
  real_T PrevY_n;                      /* '<S66>/Rate Limiter ' */
  real_T LastMajorTime_m;              /* '<S66>/Rate Limiter ' */
  real_T PrevY_c;                      /* '<S66>/Rate Limiter' */
  real_T LastMajorTime_m5;             /* '<S66>/Rate Limiter' */
  real_T PrevY_hm;                     /* '<S72>/Rate Limiter' */
  real_T LastMajorTime_c;              /* '<S72>/Rate Limiter' */
  real_T PrevY_p;                      /* '<S71>/Rate Limiter ' */
  real_T LastMajorTime_i5;             /* '<S71>/Rate Limiter ' */
  real_T PrevY_c2;                     /* '<S71>/Rate Limiter' */
  real_T LastMajorTime_l;              /* '<S71>/Rate Limiter' */
  real_T PrevY_l;                      /* '<S51>/Rate Limiter' */
  real_T LastMajorTime_h;              /* '<S51>/Rate Limiter' */
  real_T PrevY_b;                      /* '<S48>/Rate Limiter' */
  real_T LastMajorTime_ml;             /* '<S48>/Rate Limiter' */
  real_T PrevY_d[2];                   /* '<S53>/Rate Limiter' */
  real_T LastMajorTime_ck;             /* '<S53>/Rate Limiter' */
  real_T StateSpace_RWORK[6];          /* '<S92>/State-Space' */
  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay_RWORK;              /* '<S84>/Transport Delay' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay_RWORK_o;            /* '<S85>/Transport Delay' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay_RWORK_e;            /* '<S82>/Transport Delay' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay_RWORK_i;            /* '<S79>/Transport Delay' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay_RWORK_k;            /* '<S78>/Transport Delay' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay_RWORK_f;            /* '<S74>/Transport Delay' */

  void *StateSpace_PWORK[15];          /* '<S92>/State-Space' */
  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK;              /* '<S84>/Transport Delay' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK_g;            /* '<S85>/Transport Delay' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK_a;            /* '<S82>/Transport Delay' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK_c;            /* '<S79>/Transport Delay' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK_gm;           /* '<S78>/Transport Delay' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK_n;            /* '<S74>/Transport Delay' */

  int_T StateSpace_IWORK[4];           /* '<S92>/State-Space' */
  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK;              /* '<S84>/Transport Delay' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK_g;            /* '<S85>/Transport Delay' */

  struct {
    int_T IcNeedsLoading;
  } Integrator_IWORK;                  /* '<S86>/Integrator' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK_h;            /* '<S82>/Transport Delay' */

  struct {
    int_T IcNeedsLoading;
  } Integrator_IWORK_o;                /* '<S83>/Integrator' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK_a;            /* '<S79>/Transport Delay' */

  struct {
    int_T IcNeedsLoading;
  } Integrator_IWORK_l;                /* '<S80>/Integrator' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK_ao;           /* '<S78>/Transport Delay' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK_ak;           /* '<S74>/Transport Delay' */

  struct {
    int_T IcNeedsLoading;
  } Integrator_IWORK_a;                /* '<S77>/Integrator' */

  boolean_T Relay_Mode;                /* '<S47>/Relay' */
  boolean_T Relay_Mode_p;              /* '<S63>/Relay' */
  boolean_T ICic_PreviousInput_n;      /* '<S45>/IC=ic' */
  boolean_T ICic_PreviousInput_c;      /* '<S44>/IC=ic' */
  boolean_T SETPriority_MODE;          /* '<S43>/SET  Priority' */
  boolean_T RESETPriority_MODE;        /* '<S43>/RESET Priority' */
  boolean_T wind_dfig_rotor_MODE;      /* '<S31>/wind_dfig_rotor' */
  boolean_T wind_dfig_grid_MODE;       /* '<S31>/wind_dfig_grid' */
  boolean_T VRegulator_MODE;           /* '<S49>/V Regulator' */
  boolean_T QRegulator_MODE;           /* '<S49>/Q Regulator' */
  boolean_T Subsystem_MODE;            /* '<S64>/Subsystem' */
  boolean_T Subsystem_MODE_o;          /* '<S64>/Subsystem ' */
} DW_ModelCopy2_T;

/* Continuous states (auto storage) */
typedef struct {
  real_T phiqr_CSTATE;                 /* '<S41>/phiqr' */
  real_T phidr_CSTATE;                 /* '<S41>/phidr' */
  real_T phiqs_CSTATE;                 /* '<S42>/phiqs' */
  real_T phids_CSTATE;                 /* '<S42>/phids' */
  real_T Integrator_CSTATE[2];         /* '<S89>/Integrator' */
  real_T Integrator_CSTATE_j;          /* '<S10>/Integrator' */
  real_T Integrator_CSTATE_c;          /* '<S35>/Integrator' */
  real_T Integrator1_CSTATE;           /* '<S35>/Integrator1' */
  real_T Integrator_CSTATE_cb;         /* '<S32>/Integrator' */
  real_T Integrator_CSTATE_c4[2];      /* '<S69>/Integrator' */
  real_T integrator_CSTATE;            /* '<S84>/integrator' */
  real_T integrator_CSTATE_p;          /* '<S85>/integrator' */
  real_T Integrator_CSTATE_p;          /* '<S86>/Integrator' */
  real_T integrator_CSTATE_f;          /* '<S82>/integrator' */
  real_T Integrator_CSTATE_a;          /* '<S83>/Integrator' */
  real_T integrator_CSTATE_d;          /* '<S79>/integrator' */
  real_T Integrator_CSTATE_e;          /* '<S80>/Integrator' */
  real_T integrator_CSTATE_k;          /* '<S78>/integrator' */
  real_T integrator_CSTATE_o;          /* '<S74>/integrator' */
  real_T Integrator_CSTATE_k;          /* '<S77>/Integrator' */
  real_T Integrator_CSTATE_e1;         /* '<S57>/Integrator' */
  real_T Integrator_CSTATE_g[2];       /* '<S56>/Integrator' */
} X_ModelCopy2_T;

/* State derivatives (auto storage) */
typedef struct {
  real_T phiqr_CSTATE;                 /* '<S41>/phiqr' */
  real_T phidr_CSTATE;                 /* '<S41>/phidr' */
  real_T phiqs_CSTATE;                 /* '<S42>/phiqs' */
  real_T phids_CSTATE;                 /* '<S42>/phids' */
  real_T Integrator_CSTATE[2];         /* '<S89>/Integrator' */
  real_T Integrator_CSTATE_j;          /* '<S10>/Integrator' */
  real_T Integrator_CSTATE_c;          /* '<S35>/Integrator' */
  real_T Integrator1_CSTATE;           /* '<S35>/Integrator1' */
  real_T Integrator_CSTATE_cb;         /* '<S32>/Integrator' */
  real_T Integrator_CSTATE_c4[2];      /* '<S69>/Integrator' */
  real_T integrator_CSTATE;            /* '<S84>/integrator' */
  real_T integrator_CSTATE_p;          /* '<S85>/integrator' */
  real_T Integrator_CSTATE_p;          /* '<S86>/Integrator' */
  real_T integrator_CSTATE_f;          /* '<S82>/integrator' */
  real_T Integrator_CSTATE_a;          /* '<S83>/Integrator' */
  real_T integrator_CSTATE_d;          /* '<S79>/integrator' */
  real_T Integrator_CSTATE_e;          /* '<S80>/Integrator' */
  real_T integrator_CSTATE_k;          /* '<S78>/integrator' */
  real_T integrator_CSTATE_o;          /* '<S74>/integrator' */
  real_T Integrator_CSTATE_k;          /* '<S77>/Integrator' */
  real_T Integrator_CSTATE_e1;         /* '<S57>/Integrator' */
  real_T Integrator_CSTATE_g[2];       /* '<S56>/Integrator' */
} XDot_ModelCopy2_T;

/* State disabled  */
typedef struct {
  boolean_T phiqr_CSTATE;              /* '<S41>/phiqr' */
  boolean_T phidr_CSTATE;              /* '<S41>/phidr' */
  boolean_T phiqs_CSTATE;              /* '<S42>/phiqs' */
  boolean_T phids_CSTATE;              /* '<S42>/phids' */
  boolean_T Integrator_CSTATE[2];      /* '<S89>/Integrator' */
  boolean_T Integrator_CSTATE_j;       /* '<S10>/Integrator' */
  boolean_T Integrator_CSTATE_c;       /* '<S35>/Integrator' */
  boolean_T Integrator1_CSTATE;        /* '<S35>/Integrator1' */
  boolean_T Integrator_CSTATE_cb;      /* '<S32>/Integrator' */
  boolean_T Integrator_CSTATE_c4[2];   /* '<S69>/Integrator' */
  boolean_T integrator_CSTATE;         /* '<S84>/integrator' */
  boolean_T integrator_CSTATE_p;       /* '<S85>/integrator' */
  boolean_T Integrator_CSTATE_p;       /* '<S86>/Integrator' */
  boolean_T integrator_CSTATE_f;       /* '<S82>/integrator' */
  boolean_T Integrator_CSTATE_a;       /* '<S83>/Integrator' */
  boolean_T integrator_CSTATE_d;       /* '<S79>/integrator' */
  boolean_T Integrator_CSTATE_e;       /* '<S80>/Integrator' */
  boolean_T integrator_CSTATE_k;       /* '<S78>/integrator' */
  boolean_T integrator_CSTATE_o;       /* '<S74>/integrator' */
  boolean_T Integrator_CSTATE_k;       /* '<S77>/Integrator' */
  boolean_T Integrator_CSTATE_e1;      /* '<S57>/Integrator' */
  boolean_T Integrator_CSTATE_g[2];    /* '<S56>/Integrator' */
} XDis_ModelCopy2_T;

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
} ODE3_IntgData;

#endif

/* Parameters (auto storage) */
struct P_ModelCopy2_T_ {
  real_T DCbusmodel_C;                 /* Mask Parameter: DCbusmodel_C
                                        * Referenced by: '<S32>/deg->rad1'
                                        */
  real_T ukV_HarmonicGeneration;       /* Mask Parameter: ukV_HarmonicGeneration
                                        * Referenced by: '<S6>/valp_nom7'
                                        */
  real_T dqaxismodelofa3phaseseriesRLbra[2];/* Mask Parameter: dqaxismodelofa3phaseseriesRLbra
                                             * Referenced by: '<S89>/Integrator'
                                             */
  real_T Subsystem1_Init;              /* Mask Parameter: Subsystem1_Init
                                        * Referenced by: '<S57>/Integrator'
                                        */
  real_T Subsystem3_Init;              /* Mask Parameter: Subsystem3_Init
                                        * Referenced by: '<S56>/Integrator'
                                        */
  real_T Subsystem1_Init_d;            /* Mask Parameter: Subsystem1_Init_d
                                        * Referenced by: '<S69>/Integrator'
                                        */
  real_T Subsystem3_Ki;                /* Mask Parameter: Subsystem3_Ki
                                        * Referenced by: '<S56>/Gain'
                                        */
  real_T Subsystem1_Ki;                /* Mask Parameter: Subsystem1_Ki
                                        * Referenced by: '<S57>/Gain'
                                        */
  real_T Subsystem2_Ki;                /* Mask Parameter: Subsystem2_Ki
                                        * Referenced by: '<S77>/Gain'
                                        */
  real_T Subsystem2_Ki_d;              /* Mask Parameter: Subsystem2_Ki_d
                                        * Referenced by: '<S80>/Gain'
                                        */
  real_T Subsystem1_Ki_c;              /* Mask Parameter: Subsystem1_Ki_c
                                        * Referenced by: '<S83>/Gain'
                                        */
  real_T Subsystem1_Ki_i;              /* Mask Parameter: Subsystem1_Ki_i
                                        * Referenced by: '<S86>/Gain'
                                        */
  real_T Subsystem1_Ki_e;              /* Mask Parameter: Subsystem1_Ki_e
                                        * Referenced by: '<S69>/Gain'
                                        */
  real_T Subsystem1_Kp;                /* Mask Parameter: Subsystem1_Kp
                                        * Referenced by: '<S57>/Gain1'
                                        */
  real_T Subsystem3_Kp;                /* Mask Parameter: Subsystem3_Kp
                                        * Referenced by: '<S56>/Gain1'
                                        */
  real_T Subsystem2_Kp;                /* Mask Parameter: Subsystem2_Kp
                                        * Referenced by: '<S77>/Gain1'
                                        */
  real_T Subsystem2_Kp_c;              /* Mask Parameter: Subsystem2_Kp_c
                                        * Referenced by: '<S80>/Gain1'
                                        */
  real_T Subsystem1_Kp_a;              /* Mask Parameter: Subsystem1_Kp_a
                                        * Referenced by: '<S83>/Gain1'
                                        */
  real_T Subsystem1_Kp_e;              /* Mask Parameter: Subsystem1_Kp_e
                                        * Referenced by: '<S86>/Gain1'
                                        */
  real_T Subsystem1_Kp_i;              /* Mask Parameter: Subsystem1_Kp_i
                                        * Referenced by: '<S69>/Gain1'
                                        */
  real_T dqaxismodelofa3phaseseriesRLb_m;/* Mask Parameter: dqaxismodelofa3phaseseriesRLb_m
                                          * Referenced by:
                                          *   '<S89>/R_choke1'
                                          *   '<S89>/R_choke3'
                                          */
  real_T SeqAGeneration_Mag_Harmo;     /* Mask Parameter: SeqAGeneration_Mag_Harmo
                                        * Referenced by: '<S7>/Phase_Harmo1'
                                        */
  real_T SeqBGeneration_Mag_Harmo;     /* Mask Parameter: SeqBGeneration_Mag_Harmo
                                        * Referenced by: '<S8>/Phase_Harmo1'
                                        */
  real_T WindTurbine_Pelec_base;       /* Mask Parameter: WindTurbine_Pelec_base
                                        * Referenced by: '<S19>/pu->pu'
                                        */
  real_T MeanValue_Period;             /* Mask Parameter: MeanValue_Period
                                        * Referenced by:
                                        *   '<S78>/Gain'
                                        *   '<S78>/Step'
                                        *   '<S78>/Transport Delay'
                                        */
  real_T MeanValue_Period_o;           /* Mask Parameter: MeanValue_Period_o
                                        * Referenced by:
                                        *   '<S74>/Gain'
                                        *   '<S74>/Step'
                                        *   '<S74>/Transport Delay'
                                        */
  real_T MeanValue_Period_h;           /* Mask Parameter: MeanValue_Period_h
                                        * Referenced by:
                                        *   '<S79>/Gain'
                                        *   '<S79>/Step'
                                        *   '<S79>/Transport Delay'
                                        */
  real_T MeanValue_Period_i;           /* Mask Parameter: MeanValue_Period_i
                                        * Referenced by:
                                        *   '<S82>/Gain'
                                        *   '<S82>/Step'
                                        *   '<S82>/Transport Delay'
                                        */
  real_T MeanValue_Period_d;           /* Mask Parameter: MeanValue_Period_d
                                        * Referenced by:
                                        *   '<S84>/Gain'
                                        *   '<S84>/Step'
                                        *   '<S84>/Transport Delay'
                                        */
  real_T MeanValue1_Period;            /* Mask Parameter: MeanValue1_Period
                                        * Referenced by:
                                        *   '<S85>/Gain'
                                        *   '<S85>/Step'
                                        *   '<S85>/Transport Delay'
                                        */
  real_T SeqAGeneration_Phase_Harmo;   /* Mask Parameter: SeqAGeneration_Phase_Harmo
                                        * Referenced by: '<S7>/Phase_Harmo'
                                        */
  real_T SeqBGeneration_Phase_Harmo;   /* Mask Parameter: SeqBGeneration_Phase_Harmo
                                        * Referenced by: '<S8>/Phase_Harmo'
                                        */
  real_T WindTurbineDoublyFedInductionGe;/* Mask Parameter: WindTurbineDoublyFedInductionGe
                                          * Referenced by:
                                          *   '<S19>/pu->pu'
                                          *   '<S71>/pu->pu  '
                                          */
  real_T DCbusmodel_Pnom;              /* Mask Parameter: DCbusmodel_Pnom
                                        * Referenced by: '<S32>/pu->W'
                                        */
  real_T Bistable1_Qpriority;          /* Mask Parameter: Bistable1_Qpriority
                                        * Referenced by: '<S43>/Constant2'
                                        */
  real_T dqaxismodelofa3phaseseriesRLb_p;/* Mask Parameter: dqaxismodelofa3phaseseriesRLb_p
                                          * Referenced by: '<S89>/R_choke'
                                          */
  real_T SeqAGeneration_Seq_Harmo;     /* Mask Parameter: SeqAGeneration_Seq_Harmo
                                        * Referenced by: '<S7>/Phase_Harmo2'
                                        */
  real_T SeqBGeneration_Seq_Harmo;     /* Mask Parameter: SeqBGeneration_Seq_Harmo
                                        * Referenced by: '<S8>/Phase_Harmo2'
                                        */
  real_T VariationSubSystem_Toff_Variati;/* Mask Parameter: VariationSubSystem_Toff_Variati
                                          * Referenced by: '<S10>/Step1'
                                          */
  real_T VariationSubSystem_Ton_Variatio;/* Mask Parameter: VariationSubSystem_Ton_Variatio
                                          * Referenced by: '<S10>/Step'
                                          */
  real_T ukV_VariationEntity;          /* Mask Parameter: ukV_VariationEntity
                                        * Referenced by: '<S6>/valp_nom3'
                                        */
  real_T VariationSubSystem_VariationFre;/* Mask Parameter: VariationSubSystem_VariationFre
                                          * Referenced by: '<S10>/valp_nom9'
                                          */
  real_T VariationSubSystem_VariationMag;/* Mask Parameter: VariationSubSystem_VariationMag
                                          * Referenced by: '<S10>/valp_nom8'
                                          */
  real_T VariationSubSystem_VariationRat;/* Mask Parameter: VariationSubSystem_VariationRat
                                          * Referenced by: '<S10>/valp_nom7'
                                          */
  real_T VariationSubSystem_VariationSte;/* Mask Parameter: VariationSubSystem_VariationSte
                                          * Referenced by: '<S10>/valp_nom6'
                                          */
  real_T DCbusmodel_Vdc_Init;          /* Mask Parameter: DCbusmodel_Vdc_Init
                                        * Referenced by: '<S32>/Integrator'
                                        */
  real_T CompareToConstant_const;      /* Mask Parameter: CompareToConstant_const
                                        * Referenced by: '<S3>/Constant'
                                        */
  real_T WindTurbineDoublyFedInduction_n;/* Mask Parameter: WindTurbineDoublyFedInduction_n
                                          * Referenced by: '<S31>/pitch_gain'
                                          */
  real_T WindTurbineDoublyFedInduction_m;/* Mask Parameter: WindTurbineDoublyFedInduction_m
                                          * Referenced by: '<S31>/0-pitch_max'
                                          */
  real_T WindTurbineDoublyFedInduction_o;/* Mask Parameter: WindTurbineDoublyFedInduction_o
                                          * Referenced by: '<S31>/Rate Limiter   '
                                          */
  real_T WindTurbineDoublyFedInduction_k;/* Mask Parameter: WindTurbineDoublyFedInduction_k
                                          * Referenced by:
                                          *   '<S19>/pu->pu'
                                          *   '<S71>/Gain '
                                          *   '<S71>/0-power_C'
                                          *   '<S76>/Constant3'
                                          */
  real_T WindTurbine_speed_nom;        /* Mask Parameter: WindTurbine_speed_nom
                                        * Referenced by: '<S19>/pu->pu '
                                        */
  real_T WindTurbineDoublyFedInductio_of;/* Mask Parameter: WindTurbineDoublyFedInductio_of
                                          * Referenced by: '<S19>/1//wind_base'
                                          */
  boolean_T Bistable1_ic;              /* Mask Parameter: Bistable1_ic
                                        * Referenced by:
                                        *   '<S44>/IC=ic'
                                        *   '<S45>/IC=ic'
                                        */
  real_T Gain4_Gain;                   /* Expression: pi/180
                                        * Referenced by: '<S6>/Gain4'
                                        */
  real_T LookUpTable_XData[7];         /* Expression: tv
                                        * Referenced by: '<S9>/Look-Up Table'
                                        */
  real_T LookUpTable_YData[7];         /* Expression: opv
                                        * Referenced by: '<S9>/Look-Up Table'
                                        */
  real_T Negativesequence_Value[3];    /* Expression: [0 2*pi/3 -2*pi/3]
                                        * Referenced by: '<S7>/Negative-sequence'
                                        */
  real_T Positivesequence_Value[3];    /* Expression: [0 -2*pi/3 2*pi/3]
                                        * Referenced by: '<S7>/Positive-sequence'
                                        */
  real_T Zerosequence_Value[3];        /* Expression: [0 0 0]
                                        * Referenced by: '<S7>/Zero-sequence'
                                        */
  real_T Negativesequence_Value_a[3];  /* Expression: [0 2*pi/3 -2*pi/3]
                                        * Referenced by: '<S8>/Negative-sequence'
                                        */
  real_T Positivesequence_Value_h[3];  /* Expression: [0 -2*pi/3 2*pi/3]
                                        * Referenced by: '<S8>/Positive-sequence'
                                        */
  real_T Zerosequence_Value_m[3];      /* Expression: [0 0 0]
                                        * Referenced by: '<S8>/Zero-sequence'
                                        */
  real_T Constant2_Value;              /* Expression: 1
                                        * Referenced by: '<S10>/Constant2'
                                        */
  real_T Avoiddivisionbyzero_UpperSat; /* Expression: 1e6
                                        * Referenced by: '<S19>/Avoid division by zero'
                                        */
  real_T Avoiddivisionbyzero_LowerSat; /* Expression: 1e-6
                                        * Referenced by: '<S19>/Avoid division by zero'
                                        */
  real_T Gain_Gain;                    /* Expression: -1
                                        * Referenced by: '<S19>/Gain'
                                        */
  real_T Llr3_Gain;                    /* Expression: Rr
                                        * Referenced by: '<S41>/1\Llr3'
                                        */
  real_T web1_Gain;                    /* Expression: web
                                        * Referenced by: '<S41>/web1'
                                        */
  real_T Llr1_Gain;                    /* Expression: Rr
                                        * Referenced by: '<S41>/1\Llr1'
                                        */
  real_T web_Gain;                     /* Expression: web
                                        * Referenced by: '<S41>/web'
                                        */
  real_T Llr1_Gain_o;                  /* Expression: Rs
                                        * Referenced by: '<S42>/1\Llr1'
                                        */
  real_T web_Gain_d;                   /* Expression: web
                                        * Referenced by: '<S42>/web'
                                        */
  real_T Llr3_Gain_l;                  /* Expression: Rs
                                        * Referenced by: '<S42>/1\Llr3'
                                        */
  real_T web1_Gain_h;                  /* Expression: web
                                        * Referenced by: '<S42>/web1'
                                        */
  real_T Vdq_ctrl_grid_conv_Y0;        /* Computed Parameter: Vdq_ctrl_grid_conv_Y0
                                        * Referenced by: '<S48>/Vdq_ctrl_grid_conv'
                                        */
  real_T Constant1_Value;              /* Expression: L_RL
                                        * Referenced by: '<S50>/Constant1'
                                        */
  real_T Constant2_Value_c;            /* Expression: L_RL
                                        * Referenced by: '<S50>/Constant2'
                                        */
  real_T Constant3_Value;              /* Expression: R_RL
                                        * Referenced by: '<S50>/Constant3'
                                        */
  real_T Constant4_Value;              /* Expression: R_RL
                                        * Referenced by: '<S50>/Constant4'
                                        */
  real_T Vdc_refV_Value;               /* Expression: Vdc_nom
                                        * Referenced by: '<S48>/Vdc_ref (V)'
                                        */
  real_T Integrator_UpperSat;          /* Expression: Upper_Limit
                                        * Referenced by: '<S57>/Integrator'
                                        */
  real_T Integrator_LowerSat;          /* Expression: Lower_Limit
                                        * Referenced by: '<S57>/Integrator'
                                        */
  real_T Saturation_UpperSat;          /* Expression: Upper_Limit
                                        * Referenced by: '<S57>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: Lower_Limit
                                        * Referenced by: '<S57>/Saturation'
                                        */
  real_T RateLimiter_RisingLim;        /* Expression: current_slew_rate
                                        * Referenced by: '<S51>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim;       /* Expression: -current_slew_rate
                                        * Referenced by: '<S51>/Rate Limiter'
                                        */
  real_T RateLimiter_RisingLim_h;      /* Expression: current_slew_rate
                                        * Referenced by: '<S48>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim_b;     /* Expression: -current_slew_rate
                                        * Referenced by: '<S48>/Rate Limiter'
                                        */
  real_T Igrid_conv_max2_Value;        /* Expression: Imax_grid_conv*Imax_grid_conv
                                        * Referenced by: '<S53>/Igrid_conv_max^2'
                                        */
  real_T Switch_Threshold;             /* Expression: Imax_grid_conv
                                        * Referenced by: '<S53>/Switch'
                                        */
  real_T RateLimiter_RisingLim_p;      /* Expression: current_slew_rate
                                        * Referenced by: '<S53>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim_d;     /* Expression: -current_slew_rate
                                        * Referenced by: '<S53>/Rate Limiter'
                                        */
  real_T Integrator_UpperSat_d;        /* Expression: Upper_Limit
                                        * Referenced by: '<S56>/Integrator'
                                        */
  real_T Integrator_LowerSat_d;        /* Expression: Lower_Limit
                                        * Referenced by: '<S56>/Integrator'
                                        */
  real_T Saturation_UpperSat_m;        /* Expression: Upper_Limit
                                        * Referenced by: '<S56>/Saturation'
                                        */
  real_T Saturation_LowerSat_c;        /* Expression: Lower_Limit
                                        * Referenced by: '<S56>/Saturation'
                                        */
  real_T K_Value;                      /* Expression: Vnom*2*sqrt(2/3)
                                        * Referenced by: '<S54>/K'
                                        */
  real_T Avoiddivisionbyzero_UpperSat_p;/* Expression: 1e6
                                         * Referenced by: '<S54>/Avoid division by zero'
                                         */
  real_T Avoiddivisionbyzero_LowerSat_l;/* Expression: 1e-6
                                         * Referenced by: '<S54>/Avoid division by zero'
                                         */
  real_T u_UpperSat;                   /* Expression: 1
                                        * Referenced by: '<S54>/0-1'
                                        */
  real_T u_LowerSat;                   /* Expression: 0
                                        * Referenced by: '<S54>/0-1'
                                        */
  real_T speed_Aspeed_B_UpperSat;      /* Expression: speed_B
                                        * Referenced by: '<S71>/speed_A-speed_B'
                                        */
  real_T speed_Aspeed_B_LowerSat;      /* Expression: speed_A
                                        * Referenced by: '<S71>/speed_A-speed_B'
                                        */
  real_T pupu_Gain;                    /* Expression: 1/speed_C
                                        * Referenced by: '<S71>/pu->pu '
                                        */
  real_T power_C_LowerSat;             /* Expression: 0
                                        * Referenced by: '<S71>/0-power_C'
                                        */
  real_T Iqr_Y0;                       /* Computed Parameter: Iqr_Y0
                                        * Referenced by: '<S71>/Iqr+'
                                        */
  real_T Constant2_Value_p;            /* Expression: power_D
                                        * Referenced by: '<S76>/Constant2'
                                        */
  real_T Constant4_Value_h;            /* Expression: speed_D
                                        * Referenced by: '<S76>/Constant4'
                                        */
  real_T Constant5_Value;              /* Expression: speed_C
                                        * Referenced by: '<S76>/Constant5'
                                        */
  real_T Constant_Value;               /* Expression: speed_C
                                        * Referenced by: '<S71>/Constant'
                                        */
  real_T inf_UpperSat;                 /* Expression: inf
                                        * Referenced by: '<S71>/0-inf'
                                        */
  real_T inf_LowerSat;                 /* Expression: 0
                                        * Referenced by: '<S71>/0-inf'
                                        */
  real_T Constant1_Value_b;            /* Expression: speed_A
                                        * Referenced by: '<S71>/Constant1'
                                        */
  real_T Constant6_Value;              /* Expression: power_B
                                        * Referenced by: '<S75>/Constant6'
                                        */
  real_T Constant7_Value;              /* Expression: power_A
                                        * Referenced by: '<S75>/Constant7'
                                        */
  real_T Constant8_Value;              /* Expression: speed_B
                                        * Referenced by: '<S75>/Constant8'
                                        */
  real_T Constant9_Value;              /* Expression: speed_A
                                        * Referenced by: '<S75>/Constant9'
                                        */
  real_T Switch_Threshold_n;           /* Expression: speed_B
                                        * Referenced by: '<S71>/Switch'
                                        */
  real_T RateLimiter_RisingLim_e;      /* Expression: power_slew_rate
                                        * Referenced by: '<S71>/Rate Limiter '
                                        */
  real_T RateLimiter_FallingLim_de;    /* Expression: -power_slew_rate
                                        * Referenced by: '<S71>/Rate Limiter '
                                        */
  real_T u_UpperSat_m;                 /* Expression: 1
                                        * Referenced by: '<S71>/0-1'
                                        */
  real_T u_LowerSat_h;                 /* Expression: 0
                                        * Referenced by: '<S71>/0-1'
                                        */
  real_T integrator_IC;                /* Expression: 0
                                        * Referenced by: '<S78>/integrator'
                                        */
  real_T TransportDelay_InitOutput;    /* Expression: 0
                                        * Referenced by: '<S78>/Transport Delay'
                                        */
  real_T Step_Y0;                      /* Expression: 0
                                        * Referenced by: '<S78>/Step'
                                        */
  real_T Step_YFinal;                  /* Expression: 1
                                        * Referenced by: '<S78>/Step'
                                        */
  real_T FrictionFactor_Gain;          /* Expression: F
                                        * Referenced by: '<S73>/Friction Factor'
                                        */
  real_T Rs_Value;                     /* Expression: Rs
                                        * Referenced by: '<S73>/Rs'
                                        */
  real_T Rr_Value;                     /* Expression: Rr
                                        * Referenced by: '<S73>/Rr'
                                        */
  real_T Constant3_Value_p;            /* Expression: R_RL
                                        * Referenced by: '<S73>/Constant3'
                                        */
  real_T Switch_Threshold_f;           /* Expression: 0.5
                                        * Referenced by: '<S78>/Switch'
                                        */
  real_T inf_UpperSat_o;               /* Expression: inf
                                        * Referenced by: '<S71>/0-inf '
                                        */
  real_T inf_LowerSat_a;               /* Expression: 0
                                        * Referenced by: '<S71>/0-inf '
                                        */
  real_T integrator_IC_l;              /* Expression: 0
                                        * Referenced by: '<S74>/integrator'
                                        */
  real_T TransportDelay_InitOutput_e;  /* Expression: 0
                                        * Referenced by: '<S74>/Transport Delay'
                                        */
  real_T Step_Y0_k;                    /* Expression: 0
                                        * Referenced by: '<S74>/Step'
                                        */
  real_T Step_YFinal_a;                /* Expression: 1
                                        * Referenced by: '<S74>/Step'
                                        */
  real_T Switch_Threshold_k;           /* Expression: 0.5
                                        * Referenced by: '<S74>/Switch'
                                        */
  real_T Integrator_UpperSat_a;        /* Expression: Upper_Limit
                                        * Referenced by: '<S77>/Integrator'
                                        */
  real_T Integrator_LowerSat_b;        /* Expression: Lower_Limit
                                        * Referenced by: '<S77>/Integrator'
                                        */
  real_T Saturation_UpperSat_b;        /* Expression: Upper_Limit
                                        * Referenced by: '<S77>/Saturation'
                                        */
  real_T Saturation_LowerSat_g;        /* Expression: Lower_Limit
                                        * Referenced by: '<S77>/Saturation'
                                        */
  real_T RateLimiter_RisingLim_n;      /* Expression: current_slew_rate
                                        * Referenced by: '<S71>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim_dw;    /* Expression: -current_slew_rate
                                        * Referenced by: '<S71>/Rate Limiter'
                                        */
  real_T Iqr_Y0_b;                     /* Computed Parameter: Iqr_Y0_b
                                        * Referenced by: '<S72>/Iqr+'
                                        */
  real_T integrator_IC_o;              /* Expression: 0
                                        * Referenced by: '<S79>/integrator'
                                        */
  real_T TransportDelay_InitOutput_c;  /* Expression: 0
                                        * Referenced by: '<S79>/Transport Delay'
                                        */
  real_T Step_Y0_i;                    /* Expression: 0
                                        * Referenced by: '<S79>/Step'
                                        */
  real_T Step_YFinal_m;                /* Expression: 1
                                        * Referenced by: '<S79>/Step'
                                        */
  real_T Switch_Threshold_e;           /* Expression: 0.5
                                        * Referenced by: '<S79>/Switch'
                                        */
  real_T Integrator_UpperSat_g;        /* Expression: Upper_Limit
                                        * Referenced by: '<S80>/Integrator'
                                        */
  real_T Integrator_LowerSat_g;        /* Expression: Lower_Limit
                                        * Referenced by: '<S80>/Integrator'
                                        */
  real_T Saturation_UpperSat_e;        /* Expression: Upper_Limit
                                        * Referenced by: '<S80>/Saturation'
                                        */
  real_T Saturation_LowerSat_j;        /* Expression: Lower_Limit
                                        * Referenced by: '<S80>/Saturation'
                                        */
  real_T RateLimiter_RisingLim_d;      /* Expression: current_slew_rate
                                        * Referenced by: '<S72>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim_a;     /* Expression: -current_slew_rate
                                        * Referenced by: '<S72>/Rate Limiter'
                                        */
  real_T Idr_Y0;                       /* Computed Parameter: Idr_Y0
                                        * Referenced by: '<S66>/Idr+'
                                        */
  real_T integrator_IC_k;              /* Expression: 0
                                        * Referenced by: '<S82>/integrator'
                                        */
  real_T TransportDelay_InitOutput_d;  /* Expression: 0
                                        * Referenced by: '<S82>/Transport Delay'
                                        */
  real_T Step_Y0_e;                    /* Expression: 0
                                        * Referenced by: '<S82>/Step'
                                        */
  real_T Step_YFinal_i;                /* Expression: 1
                                        * Referenced by: '<S82>/Step'
                                        */
  real_T RateLimiter_RisingLim_dv;     /* Expression: Q_slew_rate
                                        * Referenced by: '<S66>/Rate Limiter '
                                        */
  real_T RateLimiter_FallingLim_g;     /* Expression: -Q_slew_rate
                                        * Referenced by: '<S66>/Rate Limiter '
                                        */
  real_T Switch_Threshold_p;           /* Expression: 0.5
                                        * Referenced by: '<S82>/Switch'
                                        */
  real_T Integrator_UpperSat_gb;       /* Expression: Upper_Limit
                                        * Referenced by: '<S83>/Integrator'
                                        */
  real_T Integrator_LowerSat_n;        /* Expression: Lower_Limit
                                        * Referenced by: '<S83>/Integrator'
                                        */
  real_T Saturation_UpperSat_i;        /* Expression: Upper_Limit
                                        * Referenced by: '<S83>/Saturation'
                                        */
  real_T Saturation_LowerSat_m;        /* Expression: Lower_Limit
                                        * Referenced by: '<S83>/Saturation'
                                        */
  real_T RateLimiter_RisingLim_i;      /* Expression: current_slew_rate
                                        * Referenced by: '<S66>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim_e;     /* Expression: -current_slew_rate
                                        * Referenced by: '<S66>/Rate Limiter'
                                        */
  real_T Idr_Y0_m;                     /* Computed Parameter: Idr_Y0_m
                                        * Referenced by: '<S67>/Idr+'
                                        */
  real_T Droop_Gain;                   /* Expression: Xs
                                        * Referenced by: '<S67>/Droop'
                                        */
  real_T integrator_IC_a;              /* Expression: 0
                                        * Referenced by: '<S84>/integrator'
                                        */
  real_T TransportDelay_InitOutput_i;  /* Expression: 0
                                        * Referenced by: '<S84>/Transport Delay'
                                        */
  real_T Step_Y0_o;                    /* Expression: 0
                                        * Referenced by: '<S84>/Step'
                                        */
  real_T Step_YFinal_l;                /* Expression: 1
                                        * Referenced by: '<S84>/Step'
                                        */
  real_T RateLimiter_RisingLim_a;      /* Expression: V_slew_rate
                                        * Referenced by: '<S67>/Rate Limiter '
                                        */
  real_T RateLimiter_FallingLim_o;     /* Expression: -V_slew_rate
                                        * Referenced by: '<S67>/Rate Limiter '
                                        */
  real_T Switch_Threshold_a;           /* Expression: 0.5
                                        * Referenced by: '<S84>/Switch'
                                        */
  real_T integrator_IC_d;              /* Expression: 0
                                        * Referenced by: '<S85>/integrator'
                                        */
  real_T TransportDelay_InitOutput_a;  /* Expression: 0
                                        * Referenced by: '<S85>/Transport Delay'
                                        */
  real_T Step_Y0_a;                    /* Expression: 0
                                        * Referenced by: '<S85>/Step'
                                        */
  real_T Step_YFinal_f;                /* Expression: 1
                                        * Referenced by: '<S85>/Step'
                                        */
  real_T Switch_Threshold_j;           /* Expression: 0.5
                                        * Referenced by: '<S85>/Switch'
                                        */
  real_T Integrator_UpperSat_i;        /* Expression: Upper_Limit
                                        * Referenced by: '<S86>/Integrator'
                                        */
  real_T Integrator_LowerSat_o;        /* Expression: Lower_Limit
                                        * Referenced by: '<S86>/Integrator'
                                        */
  real_T Saturation_UpperSat_h;        /* Expression: Upper_Limit
                                        * Referenced by: '<S86>/Saturation'
                                        */
  real_T Saturation_LowerSat_k;        /* Expression: Lower_Limit
                                        * Referenced by: '<S86>/Saturation'
                                        */
  real_T RateLimiter_RisingLim_as;     /* Expression: current_slew_rate
                                        * Referenced by: '<S67>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim_bd;    /* Expression: -current_slew_rate
                                        * Referenced by: '<S67>/Rate Limiter'
                                        */
  real_T Vdq_ctrl_rotor_conv_Y0;       /* Computed Parameter: Vdq_ctrl_rotor_conv_Y0
                                        * Referenced by: '<S49>/Vdq_ctrl_rotor_conv'
                                        */
  real_T LlrLm1_Value;                 /* Expression: Llr+Lm
                                        * Referenced by: '<S60>/Llr+Lm1'
                                        */
  real_T LlrLm2_Value;                 /* Expression: Llr+Lm
                                        * Referenced by: '<S60>/Llr+Lm2'
                                        */
  real_T Lm1_Value;                    /* Expression: Lm
                                        * Referenced by: '<S60>/Lm1'
                                        */
  real_T Lm2_Value;                    /* Expression: Lm
                                        * Referenced by: '<S60>/Lm2'
                                        */
  real_T Lm3_Value;                    /* Expression: Rr
                                        * Referenced by: '<S60>/Lm3'
                                        */
  real_T Lm4_Value;                    /* Expression: Rr
                                        * Referenced by: '<S60>/Lm4'
                                        */
  real_T w1pu_Value;                   /* Expression: 1
                                        * Referenced by: '<S60>/w=1pu'
                                        */
  real_T Constant1_Value_h;            /* Expression: Lm
                                        * Referenced by: '<S63>/Constant1'
                                        */
  real_T Relay_OnVal;                  /* Expression: 0.5
                                        * Referenced by: '<S63>/Relay'
                                        */
  real_T Relay_OffVal;                 /* Expression: 0.1
                                        * Referenced by: '<S63>/Relay'
                                        */
  real_T Relay_YOn;                    /* Expression: 1
                                        * Referenced by: '<S63>/Relay'
                                        */
  real_T Relay_YOff;                   /* Expression: 0
                                        * Referenced by: '<S63>/Relay'
                                        */
  real_T ICic_X0;                      /* Expression: -pi/2
                                        * Referenced by: '<S63>/IC=ic'
                                        */
  real_T u_Threshold;                  /* Expression: 0.5
                                        * Referenced by: '<S63>/20%'
                                        */
  real_T ExternalTorque1_Value;        /* Expression: Vcontrol
                                        * Referenced by: '<S49>/ExternalTorque1'
                                        */
  real_T Switch_Threshold_nd;          /* Expression: 0.5
                                        * Referenced by: '<S49>/Switch'
                                        */
  real_T ExternalTorque_Value;         /* Expression: Wind_On
                                        * Referenced by: '<S64>/ExternalTorque'
                                        */
  real_T Switch_Threshold_i;           /* Expression: 0.5
                                        * Referenced by: '<S64>/Switch'
                                        */
  real_T Saturation_UpperSat_a;        /* Expression: 1
                                        * Referenced by: '<S64>/Saturation'
                                        */
  real_T Saturation_LowerSat_o;        /* Expression: 0
                                        * Referenced by: '<S64>/Saturation'
                                        */
  real_T RateLimiter_RisingLim_o;      /* Expression: current_slew_rate
                                        * Referenced by: '<S64>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim_es;    /* Expression: -current_slew_rate
                                        * Referenced by: '<S64>/Rate Limiter'
                                        */
  real_T Irotor_max2_Value;            /* Expression: 1
                                        * Referenced by: '<S65>/Irotor_max^2'
                                        */
  real_T Switch_Threshold_ew;          /* Expression: 1
                                        * Referenced by: '<S65>/Switch'
                                        */
  real_T RateLimiter_RisingLim_il;     /* Expression: current_slew_rate
                                        * Referenced by: '<S65>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim_h;     /* Expression: -current_slew_rate
                                        * Referenced by: '<S65>/Rate Limiter'
                                        */
  real_T Integrator_UpperSat_l;        /* Expression: Upper_Limit
                                        * Referenced by: '<S69>/Integrator'
                                        */
  real_T Integrator_LowerSat_l;        /* Expression: Lower_Limit
                                        * Referenced by: '<S69>/Integrator'
                                        */
  real_T Saturation_UpperSat_k;        /* Expression: Upper_Limit
                                        * Referenced by: '<S69>/Saturation'
                                        */
  real_T Saturation_LowerSat_p;        /* Expression: Lower_Limit
                                        * Referenced by: '<S69>/Saturation'
                                        */
  real_T K_Value_i;                    /* Expression: Vnom*2*sqrt(2/3)
                                        * Referenced by: '<S68>/K'
                                        */
  real_T Avoiddivisionbyzero_UpperSat_c;/* Expression: 1e6
                                         * Referenced by: '<S68>/Avoid division by zero'
                                         */
  real_T Avoiddivisionbyzero_LowerSat_lp;/* Expression: 1e-6
                                          * Referenced by: '<S68>/Avoid division by zero'
                                          */
  real_T u_UpperSat_b;                 /* Expression: 1
                                        * Referenced by: '<S68>/0-1'
                                        */
  real_T u_LowerSat_l;                 /* Expression: 0
                                        * Referenced by: '<S68>/0-1'
                                        */
  real_T R_choke2_Gain[2];             /* Expression: [1 -1]
                                        * Referenced by: '<S89>/R_choke2'
                                        */
  real_T Constant_Value_e;             /* Expression: 0
                                        * Referenced by: '<S41>/Constant'
                                        */
  real_T Constant_Value_j;             /* Expression: 1
                                        * Referenced by: '<S43>/Constant'
                                        */
  real_T phiqr_IC;                     /* Expression: phirqo
                                        * Referenced by: '<S41>/phiqr'
                                        */
  real_T Constant1_Value_k;            /* Expression: 0
                                        * Referenced by: '<S41>/Constant1'
                                        */
  real_T phidr_IC;                     /* Expression: phirdo
                                        * Referenced by: '<S41>/phidr'
                                        */
  real_T Llr_Gain;                     /* Expression: 1/Llr
                                        * Referenced by: '<S40>/1\Llr'
                                        */
  real_T Constant2_Value_h;            /* Expression: 0
                                        * Referenced by: '<S42>/Constant2'
                                        */
  real_T phiqs_IC;                     /* Expression: phisqo
                                        * Referenced by: '<S42>/phiqs'
                                        */
  real_T Constant_Value_o;             /* Expression: 0
                                        * Referenced by: '<S42>/Constant'
                                        */
  real_T phids_IC;                     /* Expression: phisdo
                                        * Referenced by: '<S42>/phids'
                                        */
  real_T Lls_Gain;                     /* Expression: 1/Lls
                                        * Referenced by: '<S40>/1\Lls'
                                        */
  real_T Llr2_Gain;                    /* Expression: Lad
                                        * Referenced by: '<S40>/1\Llr2'
                                        */
  real_T Llr2_Gain_b;                  /* Expression: 1/Lls
                                        * Referenced by: '<S42>/1\Llr2'
                                        */
  real_T Llr1_Gain_c;                  /* Expression: Laq
                                        * Referenced by: '<S40>/1\Llr1'
                                        */
  real_T Llr_Gain_i;                   /* Expression: 1/Lls
                                        * Referenced by: '<S42>/1\Llr'
                                        */
  real_T Constant_Value_k;             /* Expression: 0
                                        * Referenced by: '<S89>/Constant'
                                        */
  real_T puA_Gain;                     /* Expression: Pnom/sqrt(3)/Vnom*sqrt(2)
                                        * Referenced by: '<S4>/pu->A '
                                        */
  real_T puA_Gain_o;                   /* Expression: Pnom/sqrt(3)/Vnom*sqrt(2)
                                        * Referenced by: '<S4>/pu->A  '
                                        */
  real_T valp_nom5_Value;              /* Expression: VariationType
                                        * Referenced by: '<S6>/valp_nom5'
                                        */
  real_T Constant6_Value_j;            /* Expression: 4
                                        * Referenced by: '<S6>/Constant6'
                                        */
  real_T Constant_Value_n;             /* Expression: 2
                                        * Referenced by: '<S6>/Constant'
                                        */
  real_T Constant1_Value_m;            /* Expression: 0
                                        * Referenced by: '<S10>/Constant1'
                                        */
  real_T Step1_Y0;                     /* Expression: 1
                                        * Referenced by: '<S10>/Step1'
                                        */
  real_T Step1_YFinal;                 /* Expression: 0
                                        * Referenced by: '<S10>/Step1'
                                        */
  real_T Constant3_Value_d;            /* Expression: 1
                                        * Referenced by: '<S10>/Constant3'
                                        */
  real_T Step_Y0_l;                    /* Expression: 0
                                        * Referenced by: '<S10>/Step'
                                        */
  real_T Step_YFinal_g;                /* Expression: 1
                                        * Referenced by: '<S10>/Step'
                                        */
  real_T Integrator_IC;                /* Expression: 0
                                        * Referenced by: '<S10>/Integrator'
                                        */
  real_T Gain1_Gain;                   /* Expression: 2*pi
                                        * Referenced by: '<S10>/Gain1'
                                        */
  real_T Constant5_Value_d;            /* Expression: 0
                                        * Referenced by: '<S10>/Constant5'
                                        */
  real_T Memory_X0;                    /* Expression: 0
                                        * Referenced by: '<S10>/Memory'
                                        */
  real_T Switch2_Threshold;            /* Expression: 0.5
                                        * Referenced by: '<S10>/Switch2'
                                        */
  real_T Constant1_Value_n;            /* Expression: 0
                                        * Referenced by: '<S6>/Constant1'
                                        */
  real_T valp_nom2_Value;              /* Expression: MagnitudeVps
                                        * Referenced by: '<S6>/valp_nom2'
                                        */
  real_T SinglePhase_Value;            /* Expression: VariationPhaseA
                                        * Referenced by: '<S6>/SinglePhase'
                                        */
  real_T Switch5_Threshold;            /* Expression: 0.5
                                        * Referenced by: '<S6>/Switch5'
                                        */
  real_T valp_nom_Value;               /* Expression: PhaseVps
                                        * Referenced by: '<S6>/valp_nom'
                                        */
  real_T Gain3_Gain;                   /* Expression: pi/180
                                        * Referenced by: '<S6>/Gain3'
                                        */
  real_T P1_Value[3];                  /* Expression: [0  -2*pi/3  2*pi/3]
                                        * Referenced by: '<S6>/P1'
                                        */
  real_T Constant2_Value_l;            /* Expression: 3
                                        * Referenced by: '<S6>/Constant2'
                                        */
  real_T Constant4_Value_j;            /* Expression: 0
                                        * Referenced by: '<S6>/Constant4'
                                        */
  real_T Step_Time;                    /* Expression: Ton_Harmo
                                        * Referenced by: '<S6>/Step'
                                        */
  real_T Step_Y0_f;                    /* Expression: 0
                                        * Referenced by: '<S6>/Step'
                                        */
  real_T Step_YFinal_ih;               /* Expression: 1
                                        * Referenced by: '<S6>/Step'
                                        */
  real_T Step1_Time;                   /* Expression: Toff_Harmo
                                        * Referenced by: '<S6>/Step1'
                                        */
  real_T Step1_Y0_l;                   /* Expression: 0
                                        * Referenced by: '<S6>/Step1'
                                        */
  real_T Step1_YFinal_p;               /* Expression: -1
                                        * Referenced by: '<S6>/Step1'
                                        */
  real_T Gain3_Gain_c;                 /* Expression: pi/180
                                        * Referenced by: '<S7>/Gain3'
                                        */
  real_T valp_nom2_Value_d;            /* Expression: 1
                                        * Referenced by: '<S7>/valp_nom2'
                                        */
  real_T Gain3_Gain_m;                 /* Expression: pi/180
                                        * Referenced by: '<S8>/Gain3'
                                        */
  real_T valp_nom2_Value_j;            /* Expression: 1
                                        * Referenced by: '<S8>/valp_nom2'
                                        */
  real_T StateSpace_P1_Size[2];        /* Computed Parameter: StateSpace_P1_Size
                                        * Referenced by: '<S92>/State-Space'
                                        */
  real_T StateSpace_P1[15];            /* Expression: real(H)
                                        * Referenced by: '<S92>/State-Space'
                                        */
  real_T StateSpace_P2_Size[2];        /* Computed Parameter: StateSpace_P2_Size
                                        * Referenced by: '<S92>/State-Space'
                                        */
  real_T StateSpace_P2[15];            /* Expression: imag(H)
                                        * Referenced by: '<S92>/State-Space'
                                        */
  real_T StateSpace_P3_Size[2];        /* Computed Parameter: StateSpace_P3_Size
                                        * Referenced by: '<S92>/State-Space'
                                        */
  real_T StateSpace_P4_Size[2];        /* Computed Parameter: StateSpace_P4_Size
                                        * Referenced by: '<S92>/State-Space'
                                        */
  real_T StateSpace_P4[5];             /* Expression: InputsNonZero
                                        * Referenced by: '<S92>/State-Space'
                                        */
  real_T donotdeletethisgain_Gain;     /* Expression: 1
                                        * Referenced by: '<S23>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_e;   /* Expression: 1
                                        * Referenced by: '<S24>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_f;   /* Expression: 1
                                        * Referenced by: '<S25>/do not delete this gain'
                                        */
  real_T Kv1_Gain;                     /* Expression: Kv
                                        * Referenced by: '<S15>/Kv1'
                                        */
  real_T puV_Gain;                     /* Expression: Vnom*sqrt(2)/sqrt(3)
                                        * Referenced by: '<S33>/pu->V'
                                        */
  real_T u_Gain;                       /* Expression: -1
                                        * Referenced by: '<S33>/-1 '
                                        */
  real_T puA_Gain_i;                   /* Expression: Pnom/sqrt(3)/Vnom*sqrt(2)
                                        * Referenced by: '<S33>/pu->A'
                                        */
  real_T K_Gain;                       /* Expression: 1/2
                                        * Referenced by: '<S88>/K'
                                        */
  real_T pu_Gain;                      /* Expression: 1/Pnom
                                        * Referenced by: '<S33>/-->pu'
                                        */
  real_T u_Gain_f;                     /* Expression: -1
                                        * Referenced by: '<S33>/-1'
                                        */
  real_T u0k_Gain;                     /* Expression: 9e5
                                        * Referenced by: '<Root>/900k'
                                        */
  real_T pu1_Gain;                     /* Expression: 1/Pnom
                                        * Referenced by: '<S33>/-->pu1'
                                        */
  real_T u_Gain_b;                     /* Expression: -1
                                        * Referenced by: '<S33>/-2'
                                        */
  real_T u0k__Gain;                    /* Expression: 9e5
                                        * Referenced by: '<Root>/900k_ '
                                        */
  real_T Saturation_UpperSat_c;        /* Expression: 16
                                        * Referenced by: '<Root>/Saturation'
                                        */
  real_T Saturation_LowerSat_l;        /* Expression: 0
                                        * Referenced by: '<Root>/Saturation'
                                        */
  real_T ExternalTorque_Value_p;       /* Expression: Wind_On
                                        * Referenced by: '<S4>/ExternalTorque'
                                        */
  real_T Avoiddivisionbyzero_UpperSat_o;/* Expression: 1e6
                                         * Referenced by: '<S19>/Avoid division by zero '
                                         */
  real_T Avoiddivisionbyzero_LowerSat_g;/* Expression: 1e-6
                                         * Referenced by: '<S19>/Avoid division by zero '
                                         */
  real_T Integrator_IC_m;              /* Expression: wmo
                                        * Referenced by: '<S35>/Integrator'
                                        */
  real_T lambda_nom_Gain;              /* Expression: lambda_nom
                                        * Referenced by: '<S19>/lambda_nom'
                                        */
  real_T Saturation1_UpperSat;         /* Expression: inf
                                        * Referenced by: '<S19>/Saturation1'
                                        */
  real_T Saturation1_LowerSat;         /* Expression: 1e-6
                                        * Referenced by: '<S19>/Saturation1'
                                        */
  real_T Constant2_Value_m;            /* Expression: speed_D
                                        * Referenced by: '<S31>/Constant2'
                                        */
  real_T pitch_max_LowerSat;           /* Expression: 0
                                        * Referenced by: '<S31>/0-pitch_max'
                                        */
  real_T cp_nom_Gain;                  /* Expression: 1/cp_nom
                                        * Referenced by: '<S19>/1//cp_nom'
                                        */
  real_T Switch_Threshold_c;           /* Expression: 0.5
                                        * Referenced by: '<S4>/Switch'
                                        */
  real_T F_Gain;                       /* Expression: F
                                        * Referenced by: '<S35>/F'
                                        */
  real_T _2H_Gain;                     /* Expression: 1/(2*H)
                                        * Referenced by: '<S35>/1_2H'
                                        */
  real_T Integrator1_IC;               /* Expression: tho
                                        * Referenced by: '<S35>/Integrator1'
                                        */
  real_T web_Gain_p;                   /* Expression: web
                                        * Referenced by: '<S35>/web'
                                        */
  real_T Constant7_Value_d;            /* Expression: 1
                                        * Referenced by: '<S37>/Constant7'
                                        */
  real_T Llr_Gain_if;                  /* Expression: 1/Llr
                                        * Referenced by: '<S41>/1\Llr'
                                        */
  real_T Llr2_Gain_f;                  /* Expression: 1/Llr
                                        * Referenced by: '<S41>/1\Llr2'
                                        */
  real_T ws_Value;                     /* Expression: 1
                                        * Referenced by: '<S29>/ws'
                                        */
  real_T Gain_Gain_c;                  /* Expression: 1/3
                                        * Referenced by: '<S47>/Gain'
                                        */
  real_T Relay_OnVal_o;                /* Expression: 0.5
                                        * Referenced by: '<S47>/Relay'
                                        */
  real_T Relay_OffVal_c;               /* Expression: 0.1
                                        * Referenced by: '<S47>/Relay'
                                        */
  real_T Relay_YOn_b;                  /* Expression: 1
                                        * Referenced by: '<S47>/Relay'
                                        */
  real_T Relay_YOff_p;                 /* Expression: 0
                                        * Referenced by: '<S47>/Relay'
                                        */
  real_T ICic_X0_f;                    /* Expression: -pi/2
                                        * Referenced by: '<S47>/IC=ic'
                                        */
  real_T u_Threshold_e;                /* Expression: 0.5
                                        * Referenced by: '<S47>/20%'
                                        */
  real_T Vref_Value;                   /* Expression: Vref
                                        * Referenced by: '<S4>/Vref '
                                        */
  real_T Qref_Value;                   /* Expression: Qref
                                        * Referenced by: '<S4>/Qref '
                                        */
  real_T Powerspeed_XData[6];          /* Expression: speed_power(:,1)
                                        * Referenced by: '<S31>/Power(speed)'
                                        */
  real_T Powerspeed_YData[6];          /* Expression: speed_power(:,2)
                                        * Referenced by: '<S31>/Power(speed)'
                                        */
  real_T inf_UpperSat_d;               /* Expression: inf
                                        * Referenced by: '<S31>/0-inf'
                                        */
  real_T inf_LowerSat_b;               /* Expression: 0
                                        * Referenced by: '<S31>/0-inf'
                                        */
  real_T RateLimiter_RisingLim_f;      /* Expression: power_slew_rate
                                        * Referenced by: '<S31>/Rate Limiter '
                                        */
  real_T RateLimiter_FallingLim_l;     /* Expression: -power_slew_rate
                                        * Referenced by: '<S31>/Rate Limiter '
                                        */
  real_T Gain1_Gain_h;                 /* Expression: 1/2
                                        * Referenced by: '<S46>/Gain1'
                                        */
  real_T pu_Gain_d;                    /* Expression: 1/(Vnom*sqrt(2)/sqrt(3))
                                        * Referenced by: '<S46>/->pu '
                                        */
  real_T Gain_Gain_b;                  /* Expression: 1/3
                                        * Referenced by: '<S39>/Gain'
                                        */
  real_T Iq_ref_Value;                 /* Expression: Iq_ref
                                        * Referenced by: '<S4>/Iq_ref '
                                        */
  real_T Gain_Gain_m;                  /* Expression: 1/2
                                        * Referenced by: '<S46>/Gain'
                                        */
  real_T pu_Gain_l;                    /* Expression: 1/(Vnom*sqrt(2)/sqrt(3))
                                        * Referenced by: '<S46>/->pu'
                                        */
  real_T Constant_Value_c;             /* Expression: 2
                                        * Referenced by: '<S10>/Constant'
                                        */
  real_T Constant4_Value_jk;           /* Expression: 0
                                        * Referenced by: '<S10>/Constant4'
                                        */
  real_T Switch_Threshold_l;           /* Expression: 0.5
                                        * Referenced by: '<S10>/Switch'
                                        */
  creal_T a2_Gain;                     /* Expression: exp(-j*2*pi/3)
                                        * Referenced by: '<S36>/a^2'
                                        */
  creal_T a2_Gain_b;                   /* Expression: exp(-j*2*pi/3)
                                        * Referenced by: '<S90>/a^2'
                                        */
  creal_T a3_Gain;                     /* Expression: 1/3*exp(j*2*pi/3)
                                        * Referenced by: '<S47>/a//3'
                                        */
  creal_T a23_Gain;                    /* Expression: 1/3*exp(-j*2*pi/3)
                                        * Referenced by: '<S47>/(a^2)//3'
                                        */
  creal_T a23_Gain_a;                  /* Expression: 1/3*exp(-j*2*pi/3)
                                        * Referenced by: '<S39>/(a^2)//3'
                                        */
  boolean_T Q_Y0;                      /* Computed Parameter: Q_Y0
                                        * Referenced by: '<S44>/Q'
                                        */
  boolean_T Q_Y0_b;                    /* Computed Parameter: Q_Y0_b
                                        * Referenced by: '<S44>/!Q'
                                        */
  boolean_T Q_Y0_f;                    /* Computed Parameter: Q_Y0_f
                                        * Referenced by: '<S45>/Q'
                                        */
  boolean_T Q_Y0_l;                    /* Computed Parameter: Q_Y0_l
                                        * Referenced by: '<S45>/!Q'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_ModelCopy2_T {
  struct SimStruct_tag * *childSfunctions;
  const char_T *errorStatus;
  SS_SimMode simMode;
  RTWSolverInfo solverInfo;
  RTWSolverInfo *solverInfoPtr;
  void *sfcnInfo;

  /*
   * NonInlinedSFcns:
   * The following substructure contains information regarding
   * non-inlined s-functions used in the model.
   */
  struct {
    RTWSfcnInfo sfcnInfo;
    time_T *taskTimePtrs[2];
    SimStruct childSFunctions[1];
    SimStruct *childSFunctionPtrs[1];
    struct _ssBlkInfo2 blkInfo2[1];
    struct _ssSFcnModelMethods2 methods2[1];
    struct _ssSFcnModelMethods3 methods3[1];
    struct _ssStatesInfo2 statesInfo2[1];
    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[10];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn0;
  } NonInlinedSFcns;

  /*
   * ModelData:
   * The following substructure contains information regarding
   * the data used in the model.
   */
  struct {
    X_ModelCopy2_T *contStates;
    real_T *derivs;
    boolean_T *contStateDisabled;
    boolean_T zCCacheNeedsReset;
    boolean_T derivCacheNeedsReset;
    boolean_T blkStateChange;
    real_T odeY[25];
    real_T odeF[3][25];
    ODE3_IntgData intgData;
  } ModelData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    uint32_T options;
    int_T numContStates;
    int_T numU;
    int_T numY;
    int_T numSampTimes;
    int_T numBlocks;
    int_T numBlockIO;
    int_T numBlockPrms;
    int_T numDwork;
    int_T numSFcnPrms;
    int_T numSFcns;
    int_T numIports;
    int_T numOports;
    int_T numNonSampZCs;
    int_T sysDirFeedThru;
    int_T rtwGenSfcn;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T stepSize;
    uint32_T clockTick0;
    time_T stepSize0;
    uint32_T clockTick1;
    boolean_T firstInitCondFlag;
    time_T tStart;
    time_T tFinal;
    time_T timeOfLastOutput;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *sampleTimes;
    time_T *offsetTimes;
    int_T *sampleTimeTaskIDPtr;
    int_T *sampleHits;
    int_T *perTaskSampleHits;
    time_T *t;
    time_T sampleTimesArray[2];
    time_T offsetTimesArray[2];
    int_T sampleTimeTaskIDArray[2];
    int_T sampleHitArray[2];
    int_T perTaskSampleHitsArray[4];
    time_T tArray[2];
  } Timing;
};

/* Block parameters (auto storage) */
extern P_ModelCopy2_T ModelCopy2_P;

/* Block signals (auto storage) */
extern B_ModelCopy2_T ModelCopy2_B;

/* Continuous states (auto storage) */
extern X_ModelCopy2_T ModelCopy2_X;

/* Block states (auto storage) */
extern DW_ModelCopy2_T ModelCopy2_DW;

/* Model entry point functions */
extern void ModelCopy2_initialize(void);
extern void ModelCopy2_terminate(void);

/* Customized model step function */
extern real_T ModelCopy2_custom(real_T arg_VIND);

/* Real-time Model object */
extern RT_MODEL_ModelCopy2_T *const ModelCopy2_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'ModelCopy2'
 * '<S1>'   : 'ModelCopy2/25 kV'
 * '<S2>'   : 'ModelCopy2/25 kV// 575 V 4 MVA'
 * '<S3>'   : 'ModelCopy2/Compare To Constant'
 * '<S4>'   : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)'
 * '<S5>'   : 'ModelCopy2/powergui'
 * '<S6>'   : 'ModelCopy2/25 kV/Model'
 * '<S7>'   : 'ModelCopy2/25 kV/Model/Seq A Generation'
 * '<S8>'   : 'ModelCopy2/25 kV/Model/Seq B Generation'
 * '<S9>'   : 'ModelCopy2/25 kV/Model/Timer'
 * '<S10>'  : 'ModelCopy2/25 kV/Model/Variation SubSystem'
 * '<S11>'  : 'ModelCopy2/25 kV// 575 V 4 MVA/Model'
 * '<S12>'  : 'ModelCopy2/25 kV// 575 V 4 MVA/Model/Linear'
 * '<S13>'  : 'ModelCopy2/25 kV// 575 V 4 MVA/Model/Linear1'
 * '<S14>'  : 'ModelCopy2/25 kV// 575 V 4 MVA/Model/Linear2'
 * '<S15>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/B1'
 * '<S16>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters'
 * '<S17>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Ia3'
 * '<S18>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Ia4'
 * '<S19>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Wind Turbine'
 * '<S20>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/B1/Mode I'
 * '<S21>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/B1/Mode V'
 * '<S22>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/B1/Model'
 * '<S23>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/B1/Model/U A:'
 * '<S24>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/B1/Model/U B:'
 * '<S25>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/B1/Model/U C:'
 * '<S26>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/B1/Model/U A:/Model'
 * '<S27>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/B1/Model/U B:/Model'
 * '<S28>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/B1/Model/U C:/Model'
 * '<S29>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Asynchronous machine Positive sequence phasor model Inputs and outputs are in per unit (pu)'
 * '<S30>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Bistable1'
 * '<S31>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control'
 * '<S32>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/DC bus model'
 * '<S33>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Data acquisition'
 * '<S34>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Grid-side converter currents & Converters power'
 * '<S35>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Asynchronous machine Positive sequence phasor model Inputs and outputs are in per unit (pu)/Mechanical'
 * '<S36>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Asynchronous machine Positive sequence phasor model Inputs and outputs are in per unit (pu)/Phasors Ia Ib stator'
 * '<S37>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Asynchronous machine Positive sequence phasor model Inputs and outputs are in per unit (pu)/Positive Sequence Model'
 * '<S38>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Asynchronous machine Positive sequence phasor model Inputs and outputs are in per unit (pu)/Vab Vbc'
 * '<S39>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Asynchronous machine Positive sequence phasor model Inputs and outputs are in per unit (pu)/abc2dq'
 * '<S40>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Asynchronous machine Positive sequence phasor model Inputs and outputs are in per unit (pu)/Positive Sequence Model/Mutual fluxes'
 * '<S41>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Asynchronous machine Positive sequence phasor model Inputs and outputs are in per unit (pu)/Positive Sequence Model/Rotor'
 * '<S42>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Asynchronous machine Positive sequence phasor model Inputs and outputs are in per unit (pu)/Positive Sequence Model/Stator'
 * '<S43>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Bistable1/Model'
 * '<S44>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Bistable1/Model/RESET Priority'
 * '<S45>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Bistable1/Model/SET  Priority'
 * '<S46>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/Converters voltages'
 * '<S47>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/abc2dq & Positive sequence Voltage phase angle'
 * '<S48>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_grid'
 * '<S49>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor'
 * '<S50>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_grid/Current regulator'
 * '<S51>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_grid/DC bus voltage Regulator'
 * '<S52>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_grid/Idq B1 voltage Reference frame'
 * '<S53>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_grid/Idq references'
 * '<S54>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_grid/Subsystem'
 * '<S55>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_grid/Vdq B1 voltage Reference frame'
 * '<S56>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_grid/Current regulator/Subsystem3'
 * '<S57>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_grid/DC bus voltage Regulator/Subsystem1'
 * '<S58>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_grid/Idq references/Cartesian to Polar'
 * '<S59>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_grid/Subsystem/Cartesian to Polar1'
 * '<S60>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Current Regulator'
 * '<S61>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Idq B1 voltage Reference frame'
 * '<S62>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Idq Mutual flux Reference frame'
 * '<S63>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Mutual flux'
 * '<S64>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Power Regulator'
 * '<S65>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Priority Idr'
 * '<S66>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Q Regulator'
 * '<S67>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/V Regulator'
 * '<S68>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/dq --> abc'
 * '<S69>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Current Regulator/Subsystem1'
 * '<S70>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Mutual flux/Cartesian to Polar'
 * '<S71>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Power Regulator/Subsystem'
 * '<S72>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Power Regulator/Subsystem '
 * '<S73>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Power Regulator/Subsystem/Looses'
 * '<S74>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Power Regulator/Subsystem/Mean Value'
 * '<S75>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Power Regulator/Subsystem/Slope_B_A'
 * '<S76>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Power Regulator/Subsystem/Slope_D_C'
 * '<S77>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Power Regulator/Subsystem/Subsystem2'
 * '<S78>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Power Regulator/Subsystem/Looses/Mean Value'
 * '<S79>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Power Regulator/Subsystem /Mean Value'
 * '<S80>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Power Regulator/Subsystem /Subsystem2'
 * '<S81>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Priority Idr/Cartesian to Polar'
 * '<S82>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Q Regulator/Mean Value'
 * '<S83>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/Q Regulator/Subsystem1'
 * '<S84>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/V Regulator/Mean Value'
 * '<S85>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/V Regulator/Mean Value1'
 * '<S86>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/V Regulator/Subsystem1'
 * '<S87>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Control/wind_dfig_rotor/dq --> abc/Cartesian to Polar'
 * '<S88>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Data acquisition/Power (3ph, Phasor)'
 * '<S89>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Grid-side converter currents & Converters power/d-q axis model of a 3-phase series R-L branch'
 * '<S90>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Generator & Converters/Grid-side converter currents & Converters power/dq2ab'
 * '<S91>'  : 'ModelCopy2/Wind Turbine Doubly-Fed Induction Generator (Phasor Type)/Wind Turbine/cp(lambda,beta)'
 * '<S92>'  : 'ModelCopy2/powergui/EquivalentModel1'
 * '<S93>'  : 'ModelCopy2/powergui/EquivalentModel1/Sources'
 * '<S94>'  : 'ModelCopy2/powergui/EquivalentModel1/Yout'
 * '<S95>'  : 'ModelCopy2/powergui/EquivalentModel1/conversion'
 * '<S96>'  : 'ModelCopy2/powergui/EquivalentModel1/conversion '
 */
#endif                                 /* RTW_HEADER_ModelCopy2_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
