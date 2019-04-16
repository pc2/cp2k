/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2019  CP2K developers group                         *
 *****************************************************************************/

#ifndef FFT_FPGA_H
#define FFT_FPGA_H
/******************************************************************************
 *  Authors: Arjun Ramaswami
 *****************************************************************************/
#if defined ( __PW_FPGA )

typedef struct {
  double x;
  double y;
} double2;

typedef struct {
  float x;
  float y;
} float2;

// Initialize FPGA
extern "C" int pw_fpga_initialize_();

// Finalize FPGA
extern "C" void pw_fpga_final_();

// Single precision FFT3d procedure
extern "C" void pw_fpga_fft3d_sp_(uint32_t direction, uint32_t *N, float2 *din);

// Check fpga bitstream present in directory
extern "C" bool pw_fpga_check_bitstream_(uint32_t *N);

#endif
#endif
