#ifndef _COMPLIANT_STABILIZER_DEFINES_H_
#define _COMPLIANT_STABILIZER_DEFINES_H_

/**
  * If Fz < Fzmin the leg is considered not in contact with the ground
  **/
#define DEFAULT_Fzmin 10

/**
  * Proportional gain for the CoP error
  **/
#define DEFAULT_Kx 0.1
#define DEFAULT_Ky 0.2

/**
  * Derivative gain for the CoP velocity
  **/
#define DEFAULT_Cx -0.005
#define DEFAULT_Cy -0.01

/**
  * Max and Min allowed CoM correction
  **/
#define DEFAULT_MaxLimsx 0.3
#define DEFAULT_MaxLimsy 0.15
#define DEFAULT_MaxLimsz 0.001
#define DEFAULT_MinLimsx -0.2
#define DEFAULT_MinLimsy -0.15
#define DEFAULT_MinLimsz -0.1

/**
  * Lenght and frequency of the low pass filter
  **/
#define DEFAULT_smaples2ODE 50
#define DEFAULT_freq 10

#define DEFAULT_ANKLE_X_OFFSET 0.0

#endif
