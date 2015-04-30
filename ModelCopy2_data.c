/*
 * File: ModelCopy2_data.c
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

/* Block parameters (auto storage) */
P_ModelCopy2_T ModelCopy2_P = {
  0.01,                                /* Mask Parameter: DCbusmodel_C
                                        * Referenced by: '<S32>/deg->rad1'
                                        */
  0.0,                                 /* Mask Parameter: ukV_HarmonicGeneration
                                        * Referenced by: '<S6>/valp_nom7'
                                        */

  /*  Mask Parameter: dqaxismodelofa3phaseseriesRLbra
   * Referenced by: '<S89>/Integrator'
   */
  { 0.0, 0.0 },
  0.0,                                 /* Mask Parameter: Subsystem1_Init
                                        * Referenced by: '<S57>/Integrator'
                                        */
  0.0,                                 /* Mask Parameter: Subsystem3_Init
                                        * Referenced by: '<S56>/Integrator'
                                        */
  0.0,                                 /* Mask Parameter: Subsystem1_Init_d
                                        * Referenced by: '<S69>/Integrator'
                                        */
  100.0,                               /* Mask Parameter: Subsystem3_Ki
                                        * Referenced by: '<S56>/Gain'
                                        */
  0.05,                                /* Mask Parameter: Subsystem1_Ki
                                        * Referenced by: '<S57>/Gain'
                                        */
  100.0,                               /* Mask Parameter: Subsystem2_Ki
                                        * Referenced by: '<S77>/Gain'
                                        */
  100.0,                               /* Mask Parameter: Subsystem2_Ki_d
                                        * Referenced by: '<S80>/Gain'
                                        */
  5.0,                                 /* Mask Parameter: Subsystem1_Ki_c
                                        * Referenced by: '<S83>/Gain'
                                        */
  300.0,                               /* Mask Parameter: Subsystem1_Ki_i
                                        * Referenced by: '<S86>/Gain'
                                        */
  8.0,                                 /* Mask Parameter: Subsystem1_Ki_e
                                        * Referenced by: '<S69>/Gain'
                                        */
  0.002,                               /* Mask Parameter: Subsystem1_Kp
                                        * Referenced by: '<S57>/Gain1'
                                        */
  1.0,                                 /* Mask Parameter: Subsystem3_Kp
                                        * Referenced by: '<S56>/Gain1'
                                        */
  2.0,                                 /* Mask Parameter: Subsystem2_Kp
                                        * Referenced by: '<S77>/Gain1'
                                        */
  2.0,                                 /* Mask Parameter: Subsystem2_Kp_c
                                        * Referenced by: '<S80>/Gain1'
                                        */
  0.05,                                /* Mask Parameter: Subsystem1_Kp_a
                                        * Referenced by: '<S83>/Gain1'
                                        */
  1.25,                                /* Mask Parameter: Subsystem1_Kp_e
                                        * Referenced by: '<S86>/Gain1'
                                        */
  0.3,                                 /* Mask Parameter: Subsystem1_Kp_i
                                        * Referenced by: '<S69>/Gain1'
                                        */
  0.15,                                /* Mask Parameter: dqaxismodelofa3phaseseriesRLb_m
                                        * Referenced by:
                                        *   '<S89>/R_choke1'
                                        *   '<S89>/R_choke3'
                                        */
  -2041.2414523193152,                 /* Mask Parameter: SeqAGeneration_Mag_Harmo
                                        * Referenced by: '<S7>/Phase_Harmo1'
                                        */
  0.0,                                 /* Mask Parameter: SeqBGeneration_Mag_Harmo
                                        * Referenced by: '<S8>/Phase_Harmo1'
                                        */
  870000.0,                            /* Mask Parameter: WindTurbine_Pelec_base
                                        * Referenced by: '<S19>/pu->pu'
                                        */
  0.016666666666666666,                /* Mask Parameter: MeanValue_Period
                                        * Referenced by:
                                        *   '<S78>/Gain'
                                        *   '<S78>/Step'
                                        *   '<S78>/Transport Delay'
                                        */
  0.016666666666666666,                /* Mask Parameter: MeanValue_Period_o
                                        * Referenced by:
                                        *   '<S74>/Gain'
                                        *   '<S74>/Step'
                                        *   '<S74>/Transport Delay'
                                        */
  0.016666666666666666,                /* Mask Parameter: MeanValue_Period_h
                                        * Referenced by:
                                        *   '<S79>/Gain'
                                        *   '<S79>/Step'
                                        *   '<S79>/Transport Delay'
                                        */
  0.016666666666666666,                /* Mask Parameter: MeanValue_Period_i
                                        * Referenced by:
                                        *   '<S82>/Gain'
                                        *   '<S82>/Step'
                                        *   '<S82>/Transport Delay'
                                        */
  0.016666666666666666,                /* Mask Parameter: MeanValue_Period_d
                                        * Referenced by:
                                        *   '<S84>/Gain'
                                        *   '<S84>/Step'
                                        *   '<S84>/Transport Delay'
                                        */
  0.016666666666666666,                /* Mask Parameter: MeanValue1_Period
                                        * Referenced by:
                                        *   '<S85>/Gain'
                                        *   '<S85>/Step'
                                        *   '<S85>/Transport Delay'
                                        */
  0.0,                                 /* Mask Parameter: SeqAGeneration_Phase_Harmo
                                        * Referenced by: '<S7>/Phase_Harmo'
                                        */
  35.0,                                /* Mask Parameter: SeqBGeneration_Phase_Harmo
                                        * Referenced by: '<S8>/Phase_Harmo'
                                        */
  900000.0,                            /* Mask Parameter: WindTurbineDoublyFedInductionGe
                                        * Referenced by:
                                        *   '<S19>/pu->pu'
                                        *   '<S71>/pu->pu  '
                                        */
  870000.0,                            /* Mask Parameter: DCbusmodel_Pnom
                                        * Referenced by: '<S32>/pu->W'
                                        */
  2.0,                                 /* Mask Parameter: Bistable1_Qpriority
                                        * Referenced by: '<S43>/Constant2'
                                        */
  0.0015,                              /* Mask Parameter: dqaxismodelofa3phaseseriesRLb_p
                                        * Referenced by: '<S89>/R_choke'
                                        */
  1.0,                                 /* Mask Parameter: SeqAGeneration_Seq_Harmo
                                        * Referenced by: '<S7>/Phase_Harmo2'
                                        */
  2.0,                                 /* Mask Parameter: SeqBGeneration_Seq_Harmo
                                        * Referenced by: '<S8>/Phase_Harmo2'
                                        */
  5.5,                                 /* Mask Parameter: VariationSubSystem_Toff_Variati
                                        * Referenced by: '<S10>/Step1'
                                        */
  5.0,                                 /* Mask Parameter: VariationSubSystem_Ton_Variatio
                                        * Referenced by: '<S10>/Step'
                                        */
  1.0,                                 /* Mask Parameter: ukV_VariationEntity
                                        * Referenced by: '<S6>/valp_nom3'
                                        */
  2.0,                                 /* Mask Parameter: VariationSubSystem_VariationFre
                                        * Referenced by: '<S10>/valp_nom9'
                                        */
  0.3,                                 /* Mask Parameter: VariationSubSystem_VariationMag
                                        * Referenced by: '<S10>/valp_nom8'
                                        */
  10.0,                                /* Mask Parameter: VariationSubSystem_VariationRat
                                        * Referenced by: '<S10>/valp_nom7'
                                        */
  -0.15,                               /* Mask Parameter: VariationSubSystem_VariationSte
                                        * Referenced by: '<S10>/valp_nom6'
                                        */
  1200.0,                              /* Mask Parameter: DCbusmodel_Vdc_Init
                                        * Referenced by: '<S32>/Integrator'
                                        */
  30.0,                                /* Mask Parameter: CompareToConstant_const
                                        * Referenced by: '<S3>/Constant'
                                        */
  500.0,                               /* Mask Parameter: WindTurbineDoublyFedInduction_n
                                        * Referenced by: '<S31>/pitch_gain'
                                        */
  90.0,                                /* Mask Parameter: WindTurbineDoublyFedInduction_m
                                        * Referenced by: '<S31>/0-pitch_max'
                                        */
  2.0,                                 /* Mask Parameter: WindTurbineDoublyFedInduction_o
                                        * Referenced by: '<S31>/Rate Limiter   '
                                        */
  0.78021978021978022,                 /* Mask Parameter: WindTurbineDoublyFedInduction_k
                                        * Referenced by:
                                        *   '<S19>/pu->pu'
                                        *   '<S71>/Gain '
                                        *   '<S71>/0-power_C'
                                        *   '<S76>/Constant3'
                                        */
  1.0,                                 /* Mask Parameter: WindTurbine_speed_nom
                                        * Referenced by: '<S19>/pu->pu '
                                        */
  12.0,                                /* Mask Parameter: WindTurbineDoublyFedInductio_of
                                        * Referenced by: '<S19>/1//wind_base'
                                        */
  0,                                   /* Mask Parameter: Bistable1_ic
                                        * Referenced by:
                                        *   '<S44>/IC=ic'
                                        *   '<S45>/IC=ic'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S6>/Gain4'
                                        */

  /*  Expression: tv
   * Referenced by: '<S9>/Look-Up Table'
   */
  { -1.0, 0.0, 10.0, 10.0, 10.1, 10.1, 11.1 },

  /*  Expression: opv
   * Referenced by: '<S9>/Look-Up Table'
   */
  { 1.0, 1.0, 1.0, 0.5, 0.5, 1.0, 1.0 },

  /*  Expression: [0 2*pi/3 -2*pi/3]
   * Referenced by: '<S7>/Negative-sequence'
   */
  { 0.0, 2.0943951023931953, -2.0943951023931953 },

  /*  Expression: [0 -2*pi/3 2*pi/3]
   * Referenced by: '<S7>/Positive-sequence'
   */
  { 0.0, -2.0943951023931953, 2.0943951023931953 },

  /*  Expression: [0 0 0]
   * Referenced by: '<S7>/Zero-sequence'
   */
  { 0.0, 0.0, 0.0 },

  /*  Expression: [0 2*pi/3 -2*pi/3]
   * Referenced by: '<S8>/Negative-sequence'
   */
  { 0.0, 2.0943951023931953, -2.0943951023931953 },

  /*  Expression: [0 -2*pi/3 2*pi/3]
   * Referenced by: '<S8>/Positive-sequence'
   */
  { 0.0, -2.0943951023931953, 2.0943951023931953 },

  /*  Expression: [0 0 0]
   * Referenced by: '<S8>/Zero-sequence'
   */
  { 0.0, 0.0, 0.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S10>/Constant2'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S19>/Avoid division by zero'
                                        */
  1.0E-6,                              /* Expression: 1e-6
                                        * Referenced by: '<S19>/Avoid division by zero'
                                        */
  -1.0,                                /* Expression: -1
                                        * Referenced by: '<S19>/Gain'
                                        */
  0.005,                               /* Expression: Rr
                                        * Referenced by: '<S41>/1\Llr3'
                                        */
  376.99111843077515,                  /* Expression: web
                                        * Referenced by: '<S41>/web1'
                                        */
  0.005,                               /* Expression: Rr
                                        * Referenced by: '<S41>/1\Llr1'
                                        */
  376.99111843077515,                  /* Expression: web
                                        * Referenced by: '<S41>/web'
                                        */
  0.00706,                             /* Expression: Rs
                                        * Referenced by: '<S42>/1\Llr1'
                                        */
  376.99111843077515,                  /* Expression: web
                                        * Referenced by: '<S42>/web'
                                        */
  0.00706,                             /* Expression: Rs
                                        * Referenced by: '<S42>/1\Llr3'
                                        */
  376.99111843077515,                  /* Expression: web
                                        * Referenced by: '<S42>/web1'
                                        */
  0.0,                                 /* Computed Parameter: Vdq_ctrl_grid_conv_Y0
                                        * Referenced by: '<S48>/Vdq_ctrl_grid_conv'
                                        */
  0.15,                                /* Expression: L_RL
                                        * Referenced by: '<S50>/Constant1'
                                        */
  0.15,                                /* Expression: L_RL
                                        * Referenced by: '<S50>/Constant2'
                                        */
  0.0015,                              /* Expression: R_RL
                                        * Referenced by: '<S50>/Constant3'
                                        */
  0.0015,                              /* Expression: R_RL
                                        * Referenced by: '<S50>/Constant4'
                                        */
  1200.0,                              /* Expression: Vdc_nom
                                        * Referenced by: '<S48>/Vdc_ref (V)'
                                        */
  1.0,                                 /* Expression: Upper_Limit
                                        * Referenced by: '<S57>/Integrator'
                                        */
  -1.0,                                /* Expression: Lower_Limit
                                        * Referenced by: '<S57>/Integrator'
                                        */
  1.0,                                 /* Expression: Upper_Limit
                                        * Referenced by: '<S57>/Saturation'
                                        */
  -1.0,                                /* Expression: Lower_Limit
                                        * Referenced by: '<S57>/Saturation'
                                        */
  200.0,                               /* Expression: current_slew_rate
                                        * Referenced by: '<S51>/Rate Limiter'
                                        */
  -200.0,                              /* Expression: -current_slew_rate
                                        * Referenced by: '<S51>/Rate Limiter'
                                        */
  200.0,                               /* Expression: current_slew_rate
                                        * Referenced by: '<S48>/Rate Limiter'
                                        */
  -200.0,                              /* Expression: -current_slew_rate
                                        * Referenced by: '<S48>/Rate Limiter'
                                        */
  1.0,                                 /* Expression: Imax_grid_conv*Imax_grid_conv
                                        * Referenced by: '<S53>/Igrid_conv_max^2'
                                        */
  1.0,                                 /* Expression: Imax_grid_conv
                                        * Referenced by: '<S53>/Switch'
                                        */
  200.0,                               /* Expression: current_slew_rate
                                        * Referenced by: '<S53>/Rate Limiter'
                                        */
  -200.0,                              /* Expression: -current_slew_rate
                                        * Referenced by: '<S53>/Rate Limiter'
                                        */
  0.0,                                 /* Expression: Upper_Limit
                                        * Referenced by: '<S56>/Integrator'
                                        */
  0.0,                                 /* Expression: Lower_Limit
                                        * Referenced by: '<S56>/Integrator'
                                        */
  0.0,                                 /* Expression: Upper_Limit
                                        * Referenced by: '<S56>/Saturation'
                                        */
  0.0,                                 /* Expression: Lower_Limit
                                        * Referenced by: '<S56>/Saturation'
                                        */
  938.971068066885,                    /* Expression: Vnom*2*sqrt(2/3)
                                        * Referenced by: '<S54>/K'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S54>/Avoid division by zero'
                                        */
  1.0E-6,                              /* Expression: 1e-6
                                        * Referenced by: '<S54>/Avoid division by zero'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S54>/0-1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S54>/0-1'
                                        */
  0.71,                                /* Expression: speed_B
                                        * Referenced by: '<S71>/speed_A-speed_B'
                                        */
  0.3,                                 /* Expression: speed_A
                                        * Referenced by: '<S71>/speed_A-speed_B'
                                        */
  1.0,                                 /* Expression: 1/speed_C
                                        * Referenced by: '<S71>/pu->pu '
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S71>/0-power_C'
                                        */
  0.0,                                 /* Computed Parameter: Iqr_Y0
                                        * Referenced by: '<S71>/Iqr+'
                                        */
  1.0,                                 /* Expression: power_D
                                        * Referenced by: '<S76>/Constant2'
                                        */
  1.5,                                 /* Expression: speed_D
                                        * Referenced by: '<S76>/Constant4'
                                        */
  1.0,                                 /* Expression: speed_C
                                        * Referenced by: '<S76>/Constant5'
                                        */
  1.0,                                 /* Expression: speed_C
                                        * Referenced by: '<S71>/Constant'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S71>/0-inf'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S71>/0-inf'
                                        */
  0.3,                                 /* Expression: speed_A
                                        * Referenced by: '<S71>/Constant1'
                                        */
  0.27924924175824173,                 /* Expression: power_B
                                        * Referenced by: '<S75>/Constant6'
                                        */
  0.0,                                 /* Expression: power_A
                                        * Referenced by: '<S75>/Constant7'
                                        */
  0.71,                                /* Expression: speed_B
                                        * Referenced by: '<S75>/Constant8'
                                        */
  0.3,                                 /* Expression: speed_A
                                        * Referenced by: '<S75>/Constant9'
                                        */
  0.71,                                /* Expression: speed_B
                                        * Referenced by: '<S71>/Switch'
                                        */
  1.0,                                 /* Expression: power_slew_rate
                                        * Referenced by: '<S71>/Rate Limiter '
                                        */
  -1.0,                                /* Expression: -power_slew_rate
                                        * Referenced by: '<S71>/Rate Limiter '
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S71>/0-1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S71>/0-1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S78>/integrator'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S78>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S78>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S78>/Step'
                                        */
  0.01,                                /* Expression: F
                                        * Referenced by: '<S73>/Friction Factor'
                                        */
  0.00706,                             /* Expression: Rs
                                        * Referenced by: '<S73>/Rs'
                                        */
  0.005,                               /* Expression: Rr
                                        * Referenced by: '<S73>/Rr'
                                        */
  0.0015,                              /* Expression: R_RL
                                        * Referenced by: '<S73>/Constant3'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S78>/Switch'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S71>/0-inf '
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S71>/0-inf '
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S74>/integrator'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S74>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S74>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S74>/Step'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S74>/Switch'
                                        */
  1.0,                                 /* Expression: Upper_Limit
                                        * Referenced by: '<S77>/Integrator'
                                        */
  0.0,                                 /* Expression: Lower_Limit
                                        * Referenced by: '<S77>/Integrator'
                                        */
  1.0,                                 /* Expression: Upper_Limit
                                        * Referenced by: '<S77>/Saturation'
                                        */
  0.0,                                 /* Expression: Lower_Limit
                                        * Referenced by: '<S77>/Saturation'
                                        */
  200.0,                               /* Expression: current_slew_rate
                                        * Referenced by: '<S71>/Rate Limiter'
                                        */
  -200.0,                              /* Expression: -current_slew_rate
                                        * Referenced by: '<S71>/Rate Limiter'
                                        */
  0.0,                                 /* Computed Parameter: Iqr_Y0_b
                                        * Referenced by: '<S72>/Iqr+'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S79>/integrator'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S79>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S79>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S79>/Step'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S79>/Switch'
                                        */
  1.0,                                 /* Expression: Upper_Limit
                                        * Referenced by: '<S80>/Integrator'
                                        */
  0.0,                                 /* Expression: Lower_Limit
                                        * Referenced by: '<S80>/Integrator'
                                        */
  1.0,                                 /* Expression: Upper_Limit
                                        * Referenced by: '<S80>/Saturation'
                                        */
  0.0,                                 /* Expression: Lower_Limit
                                        * Referenced by: '<S80>/Saturation'
                                        */
  200.0,                               /* Expression: current_slew_rate
                                        * Referenced by: '<S72>/Rate Limiter'
                                        */
  -200.0,                              /* Expression: -current_slew_rate
                                        * Referenced by: '<S72>/Rate Limiter'
                                        */
  0.0,                                 /* Computed Parameter: Idr_Y0
                                        * Referenced by: '<S66>/Idr+'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S82>/integrator'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S82>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S82>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S82>/Step'
                                        */
  100.0,                               /* Expression: Q_slew_rate
                                        * Referenced by: '<S66>/Rate Limiter '
                                        */
  -100.0,                              /* Expression: -Q_slew_rate
                                        * Referenced by: '<S66>/Rate Limiter '
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S82>/Switch'
                                        */
  1.0,                                 /* Expression: Upper_Limit
                                        * Referenced by: '<S83>/Integrator'
                                        */
  -1.0,                                /* Expression: Lower_Limit
                                        * Referenced by: '<S83>/Integrator'
                                        */
  1.0,                                 /* Expression: Upper_Limit
                                        * Referenced by: '<S83>/Saturation'
                                        */
  -1.0,                                /* Expression: Lower_Limit
                                        * Referenced by: '<S83>/Saturation'
                                        */
  200.0,                               /* Expression: current_slew_rate
                                        * Referenced by: '<S66>/Rate Limiter'
                                        */
  -200.0,                              /* Expression: -current_slew_rate
                                        * Referenced by: '<S66>/Rate Limiter'
                                        */
  0.0,                                 /* Computed Parameter: Idr_Y0_m
                                        * Referenced by: '<S67>/Idr+'
                                        */
  0.02,                                /* Expression: Xs
                                        * Referenced by: '<S67>/Droop'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S84>/integrator'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S84>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S84>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S84>/Step'
                                        */
  100.0,                               /* Expression: V_slew_rate
                                        * Referenced by: '<S67>/Rate Limiter '
                                        */
  -100.0,                              /* Expression: -V_slew_rate
                                        * Referenced by: '<S67>/Rate Limiter '
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S84>/Switch'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S85>/integrator'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S85>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S85>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S85>/Step'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S85>/Switch'
                                        */
  1.0,                                 /* Expression: Upper_Limit
                                        * Referenced by: '<S86>/Integrator'
                                        */
  -1.0,                                /* Expression: Lower_Limit
                                        * Referenced by: '<S86>/Integrator'
                                        */
  1.0,                                 /* Expression: Upper_Limit
                                        * Referenced by: '<S86>/Saturation'
                                        */
  -1.0,                                /* Expression: Lower_Limit
                                        * Referenced by: '<S86>/Saturation'
                                        */
  200.0,                               /* Expression: current_slew_rate
                                        * Referenced by: '<S67>/Rate Limiter'
                                        */
  -200.0,                              /* Expression: -current_slew_rate
                                        * Referenced by: '<S67>/Rate Limiter'
                                        */
  0.0,                                 /* Computed Parameter: Vdq_ctrl_rotor_conv_Y0
                                        * Referenced by: '<S49>/Vdq_ctrl_rotor_conv'
                                        */
  3.056,                               /* Expression: Llr+Lm
                                        * Referenced by: '<S60>/Llr+Lm1'
                                        */
  3.056,                               /* Expression: Llr+Lm
                                        * Referenced by: '<S60>/Llr+Lm2'
                                        */
  2.9,                                 /* Expression: Lm
                                        * Referenced by: '<S60>/Lm1'
                                        */
  2.9,                                 /* Expression: Lm
                                        * Referenced by: '<S60>/Lm2'
                                        */
  0.005,                               /* Expression: Rr
                                        * Referenced by: '<S60>/Lm3'
                                        */
  0.005,                               /* Expression: Rr
                                        * Referenced by: '<S60>/Lm4'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S60>/w=1pu'
                                        */
  2.9,                                 /* Expression: Lm
                                        * Referenced by: '<S63>/Constant1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S63>/Relay'
                                        */
  0.1,                                 /* Expression: 0.1
                                        * Referenced by: '<S63>/Relay'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S63>/Relay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S63>/Relay'
                                        */
  -1.5707963267948966,                 /* Expression: -pi/2
                                        * Referenced by: '<S63>/IC=ic'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S63>/20%'
                                        */
  1.0,                                 /* Expression: Vcontrol
                                        * Referenced by: '<S49>/ExternalTorque1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S49>/Switch'
                                        */
  1.0,                                 /* Expression: Wind_On
                                        * Referenced by: '<S64>/ExternalTorque'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S64>/Switch'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S64>/Saturation'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S64>/Saturation'
                                        */
  200.0,                               /* Expression: current_slew_rate
                                        * Referenced by: '<S64>/Rate Limiter'
                                        */
  -200.0,                              /* Expression: -current_slew_rate
                                        * Referenced by: '<S64>/Rate Limiter'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S65>/Irotor_max^2'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S65>/Switch'
                                        */
  200.0,                               /* Expression: current_slew_rate
                                        * Referenced by: '<S65>/Rate Limiter'
                                        */
  -200.0,                              /* Expression: -current_slew_rate
                                        * Referenced by: '<S65>/Rate Limiter'
                                        */
  0.0,                                 /* Expression: Upper_Limit
                                        * Referenced by: '<S69>/Integrator'
                                        */
  0.0,                                 /* Expression: Lower_Limit
                                        * Referenced by: '<S69>/Integrator'
                                        */
  0.0,                                 /* Expression: Upper_Limit
                                        * Referenced by: '<S69>/Saturation'
                                        */
  0.0,                                 /* Expression: Lower_Limit
                                        * Referenced by: '<S69>/Saturation'
                                        */
  938.971068066885,                    /* Expression: Vnom*2*sqrt(2/3)
                                        * Referenced by: '<S68>/K'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S68>/Avoid division by zero'
                                        */
  1.0E-6,                              /* Expression: 1e-6
                                        * Referenced by: '<S68>/Avoid division by zero'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S68>/0-1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S68>/0-1'
                                        */

  /*  Expression: [1 -1]
   * Referenced by: '<S89>/R_choke2'
   */
  { 1.0, -1.0 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S41>/Constant'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S43>/Constant'
                                        */
  0.0,                                 /* Expression: phirqo
                                        * Referenced by: '<S41>/phiqr'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S41>/Constant1'
                                        */
  0.0,                                 /* Expression: phirdo
                                        * Referenced by: '<S41>/phidr'
                                        */
  6.4102564102564106,                  /* Expression: 1/Llr
                                        * Referenced by: '<S40>/1\Llr'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S42>/Constant2'
                                        */
  0.0,                                 /* Expression: phisqo
                                        * Referenced by: '<S42>/phiqs'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S42>/Constant'
                                        */
  0.0,                                 /* Expression: phisdo
                                        * Referenced by: '<S42>/phids'
                                        */
  5.8479532163742682,                  /* Expression: 1/Lls
                                        * Referenced by: '<S40>/1\Lls'
                                        */
  0.079345953131154,                   /* Expression: Lad
                                        * Referenced by: '<S40>/1\Llr2'
                                        */
  5.8479532163742682,                  /* Expression: 1/Lls
                                        * Referenced by: '<S42>/1\Llr2'
                                        */
  0.079345953131154,                   /* Expression: Laq
                                        * Referenced by: '<S40>/1\Llr1'
                                        */
  5.8479532163742682,                  /* Expression: 1/Lls
                                        * Referenced by: '<S42>/1\Llr'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S89>/Constant'
                                        */
  1235.3948267949945,                  /* Expression: Pnom/sqrt(3)/Vnom*sqrt(2)
                                        * Referenced by: '<S4>/pu->A '
                                        */
  1235.3948267949945,                  /* Expression: Pnom/sqrt(3)/Vnom*sqrt(2)
                                        * Referenced by: '<S4>/pu->A  '
                                        */
  4.0,                                 /* Expression: VariationType
                                        * Referenced by: '<S6>/valp_nom5'
                                        */
  4.0,                                 /* Expression: 4
                                        * Referenced by: '<S6>/Constant6'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S6>/Constant'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S10>/Constant1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S10>/Step1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S10>/Step1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S10>/Constant3'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S10>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S10>/Step'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S10>/Integrator'
                                        */
  6.2831853071795862,                  /* Expression: 2*pi
                                        * Referenced by: '<S10>/Gain1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S10>/Constant5'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S10>/Memory'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S10>/Switch2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S6>/Constant1'
                                        */
  20412.414523193151,                  /* Expression: MagnitudeVps
                                        * Referenced by: '<S6>/valp_nom2'
                                        */
  0.0,                                 /* Expression: VariationPhaseA
                                        * Referenced by: '<S6>/SinglePhase'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S6>/Switch5'
                                        */
  0.0,                                 /* Expression: PhaseVps
                                        * Referenced by: '<S6>/valp_nom'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S6>/Gain3'
                                        */

  /*  Expression: [0  -2*pi/3  2*pi/3]
   * Referenced by: '<S6>/P1'
   */
  { 0.0, -2.0943951023931953, 2.0943951023931953 },
  3.0,                                 /* Expression: 3
                                        * Referenced by: '<S6>/Constant2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S6>/Constant4'
                                        */
  5.0,                                 /* Expression: Ton_Harmo
                                        * Referenced by: '<S6>/Step'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S6>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S6>/Step'
                                        */
  5.5,                                 /* Expression: Toff_Harmo
                                        * Referenced by: '<S6>/Step1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S6>/Step1'
                                        */
  -1.0,                                /* Expression: -1
                                        * Referenced by: '<S6>/Step1'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S7>/Gain3'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S7>/valp_nom2'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S8>/Gain3'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S8>/valp_nom2'
                                        */

  /*  Computed Parameter: StateSpace_P1_Size
   * Referenced by: '<S92>/State-Space'
   */
  { 3.0, 5.0 },

  /*  Expression: real(H)
   * Referenced by: '<S92>/State-Space'
   */
  { -0.00013925266077793777, -3.8978355590933711E-15, 0.00013925266076906105,
    -1.7952178914433514E-14, -0.00013925266076280233, 0.00013925266077662875,
    0.023077012715567207, -3.85255531681811E-5, -3.85255531681811E-5,
    -3.8525553168181238E-5, 0.023077012715567207, -3.8525553168182546E-5,
    -3.8525553168182722E-5, -3.8525553168182722E-5, 0.023077012715567207 },

  /*  Computed Parameter: StateSpace_P2_Size
   * Referenced by: '<S92>/State-Space'
   */
  { 3.0, 5.0 },

  /*  Expression: imag(H)
   * Referenced by: '<S92>/State-Space'
   */
  { -0.0041535500965102154, 0.0, 0.0041535500965102189, 3.0357660829594124E-18,
    -0.00415355009651022, 0.0041535500965102189, -3.7370821722957008E-6,
    1.2935430040715094E-6, 1.2935430041552181E-6, 1.2935430039878005E-6,
    -3.7370821721492102E-6, 1.2935430043854175E-6, 1.2935430044376364E-6,
    1.2935430041446556E-6, -3.7370821724318266E-6 },

  /*  Computed Parameter: StateSpace_P3_Size
   * Referenced by: '<S92>/State-Space'
   */
  { 0.0, 0.0 },

  /*  Computed Parameter: StateSpace_P4_Size
   * Referenced by: '<S92>/State-Space'
   */
  { 1.0, 5.0 },

  /*  Expression: InputsNonZero
   * Referenced by: '<S92>/State-Space'
   */
  { 1.0, 2.0, 3.0, 4.0, 5.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S23>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S24>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S25>/do not delete this gain'
                                        */
  0.0021299910806810243,               /* Expression: Kv
                                        * Referenced by: '<S15>/Kv1'
                                        */
  469.48553403344255,                  /* Expression: Vnom*sqrt(2)/sqrt(3)
                                        * Referenced by: '<S33>/pu->V'
                                        */
  -1.0,                                /* Expression: -1
                                        * Referenced by: '<S33>/-1 '
                                        */
  1235.3948267949945,                  /* Expression: Pnom/sqrt(3)/Vnom*sqrt(2)
                                        * Referenced by: '<S33>/pu->A'
                                        */
  0.5,                                 /* Expression: 1/2
                                        * Referenced by: '<S88>/K'
                                        */
  1.1494252873563219E-6,               /* Expression: 1/Pnom
                                        * Referenced by: '<S33>/-->pu'
                                        */
  -1.0,                                /* Expression: -1
                                        * Referenced by: '<S33>/-1'
                                        */
  900000.0,                            /* Expression: 9e5
                                        * Referenced by: '<Root>/900k'
                                        */
  1.1494252873563219E-6,               /* Expression: 1/Pnom
                                        * Referenced by: '<S33>/-->pu1'
                                        */
  -1.0,                                /* Expression: -1
                                        * Referenced by: '<S33>/-2'
                                        */
  900000.0,                            /* Expression: 9e5
                                        * Referenced by: '<Root>/900k_ '
                                        */
  16.0,                                /* Expression: 16
                                        * Referenced by: '<Root>/Saturation'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<Root>/Saturation'
                                        */
  1.0,                                 /* Expression: Wind_On
                                        * Referenced by: '<S4>/ExternalTorque'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S19>/Avoid division by zero '
                                        */
  1.0E-6,                              /* Expression: 1e-6
                                        * Referenced by: '<S19>/Avoid division by zero '
                                        */
  0.8,                                 /* Expression: wmo
                                        * Referenced by: '<S35>/Integrator'
                                        */
  8.1,                                 /* Expression: lambda_nom
                                        * Referenced by: '<S19>/lambda_nom'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S19>/Saturation1'
                                        */
  1.0E-6,                              /* Expression: 1e-6
                                        * Referenced by: '<S19>/Saturation1'
                                        */
  1.5,                                 /* Expression: speed_D
                                        * Referenced by: '<S31>/Constant2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S31>/0-pitch_max'
                                        */
  2.0833333333333335,                  /* Expression: 1/cp_nom
                                        * Referenced by: '<S19>/1//cp_nom'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S4>/Switch'
                                        */
  0.01,                                /* Expression: F
                                        * Referenced by: '<S35>/F'
                                        */
  0.0992063492063492,                  /* Expression: 1/(2*H)
                                        * Referenced by: '<S35>/1_2H'
                                        */
  0.0,                                 /* Expression: tho
                                        * Referenced by: '<S35>/Integrator1'
                                        */
  376.99111843077515,                  /* Expression: web
                                        * Referenced by: '<S35>/web'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S37>/Constant7'
                                        */
  6.4102564102564106,                  /* Expression: 1/Llr
                                        * Referenced by: '<S41>/1\Llr'
                                        */
  6.4102564102564106,                  /* Expression: 1/Llr
                                        * Referenced by: '<S41>/1\Llr2'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S29>/ws'
                                        */
  0.33333333333333331,                 /* Expression: 1/3
                                        * Referenced by: '<S47>/Gain'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S47>/Relay'
                                        */
  0.1,                                 /* Expression: 0.1
                                        * Referenced by: '<S47>/Relay'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S47>/Relay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S47>/Relay'
                                        */
  -1.5707963267948966,                 /* Expression: -pi/2
                                        * Referenced by: '<S47>/IC=ic'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S47>/20%'
                                        */
  1.0,                                 /* Expression: Vref
                                        * Referenced by: '<S4>/Vref '
                                        */
  0.0,                                 /* Expression: Qref
                                        * Referenced by: '<S4>/Qref '
                                        */

  /*  Expression: speed_power(:,1)
   * Referenced by: '<S31>/Power(speed)'
   */
  { 0.0, 0.2494, 0.6197, 0.96, 1.2, 1.5 },

  /*  Expression: speed_power(:,2)
   * Referenced by: '<S31>/Power(speed)'
   */
  { 0.0, 0.055555555555555552, 0.78888888888888886, 1.0, 1.0, 1.0 },
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S31>/0-inf'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S31>/0-inf'
                                        */
  1.0,                                 /* Expression: power_slew_rate
                                        * Referenced by: '<S31>/Rate Limiter '
                                        */
  -1.0,                                /* Expression: -power_slew_rate
                                        * Referenced by: '<S31>/Rate Limiter '
                                        */
  0.5,                                 /* Expression: 1/2
                                        * Referenced by: '<S46>/Gain1'
                                        */
  0.0021299910806810243,               /* Expression: 1/(Vnom*sqrt(2)/sqrt(3))
                                        * Referenced by: '<S46>/->pu '
                                        */
  0.33333333333333331,                 /* Expression: 1/3
                                        * Referenced by: '<S39>/Gain'
                                        */
  0.0,                                 /* Expression: Iq_ref
                                        * Referenced by: '<S4>/Iq_ref '
                                        */
  0.5,                                 /* Expression: 1/2
                                        * Referenced by: '<S46>/Gain'
                                        */
  0.0021299910806810243,               /* Expression: 1/(Vnom*sqrt(2)/sqrt(3))
                                        * Referenced by: '<S46>/->pu'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S10>/Constant'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S10>/Constant4'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S10>/Switch'
                                        */

  { -0.49999999999999978, -0.86602540378443871 },/* Expression: exp(-j*2*pi/3)
                                                  * Referenced by: '<S36>/a^2'
                                                  */

  { -0.49999999999999978, -0.86602540378443871 },/* Expression: exp(-j*2*pi/3)
                                                  * Referenced by: '<S90>/a^2'
                                                  */

  { -0.16666666666666657, 0.28867513459481287 },/* Expression: 1/3*exp(j*2*pi/3)
                                                 * Referenced by: '<S47>/a//3'
                                                 */

  { -0.16666666666666657, -0.28867513459481287 },/* Expression: 1/3*exp(-j*2*pi/3)
                                                  * Referenced by: '<S47>/(a^2)//3'
                                                  */

  { -0.16666666666666657, -0.28867513459481287 },/* Expression: 1/3*exp(-j*2*pi/3)
                                                  * Referenced by: '<S39>/(a^2)//3'
                                                  */
  0,                                   /* Computed Parameter: Q_Y0
                                        * Referenced by: '<S44>/Q'
                                        */
  0,                                   /* Computed Parameter: Q_Y0_b
                                        * Referenced by: '<S44>/!Q'
                                        */
  0,                                   /* Computed Parameter: Q_Y0_f
                                        * Referenced by: '<S45>/Q'
                                        */
  0                                    /* Computed Parameter: Q_Y0_l
                                        * Referenced by: '<S45>/!Q'
                                        */
};

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
