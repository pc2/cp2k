/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2019  CP2K developers group                         *
 *****************************************************************************/

/******************************************************************************
 *  Author: Arjun Ramaswami
 *****************************************************************************/

#if defined ( __PW_FPGA )

// global dependencies
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <unistd.h>

// common dependencies
#include "../../common/aocl/aocl_utils.h"
#include "CL/opencl.h"
using namespace aocl_utils;

// local dependencies
#include "fft_config.h"

// host variables
static cl_platform_id platform = NULL;
static cl_device_id device = NULL;
static cl_context context = NULL;
static cl_program program = NULL;

static cl_command_queue queue1 = NULL, queue2 = NULL, queue3 = NULL;
static cl_command_queue queue4 = NULL, queue5 = NULL, queue6 = NULL;

static cl_kernel fft_kernel = NULL, fft_kernel_2 = NULL;
static cl_kernel transpose_kernel_2 = NULL, fetch_kernel = NULL, transpose_kernel = NULL;

// Device memory buffers
cl_mem d_inData, d_outData, d_tmp;

// Global Variables
static cl_int status = 0;
static uint32_t flag = 0;
static uint32_t fft_size[3] = {0,0,0};
bool fft_size_changed = true;
std::string binary_file = "";

cmplx *h_outData, *h_inData;
cmplx *h_verify_tmp, *h_verify;

// Function prototypes
bool init();
void cleanup();
static void init_program(const uint32_t *N);
static void queue_setup();
static void queue_cleanup();
static void fftfpga_run_3d(const bool fsign, const uint32_t *N, cmplx *c_in);
static bool select_binary(const uint32_t *N);

// --- CODE -------------------------------------------------------------------

extern "C" int pw_fpga_initialize_(){
   status = init();
   return status;
}

extern "C" void pw_fpga_final_(){
   cleanup();
}

/******************************************************************************
 * \brief  check whether FFT3d can be computed on the FPGA or not. This depends 
 *         on the availability of bitstreams whose sizes are for now listed here 
 * \author Arjun Ramaswami
 * \param  N - integer pointer to the size of the FFT3d
 * \retval true if no board binary found in the location
 *****************************************************************************/
extern "C" bool pw_fpga_check_bitstream_(uint32_t *N){
    
    if( (N[0] == 16 && N[1] == 16 && N[2] == 16) ||
        (N[0] == 32 && N[1] == 32 && N[2] == 32) ||
        (N[0] == 64 && N[1] == 64 && N[2] == 64)  ){

        if( fft_size[0] == N[0] && fft_size[1] == N[1] && fft_size[2] == N[2] ){
            fft_size_changed = false;
        }
        else{
            fft_size[0] = N[0];
            fft_size[1] = N[1];
            fft_size[2] = N[2];
            fft_size_changed = true;
        }
        return 1;
    }
    else{
        return 0;
    }
}

/******************************************************************************
 * \brief   compute an in-place single precision complex 3D-FFT on the FPGA
 * \author  Arjun Ramaswami
 * \param   direction : direction - 1/forward, otherwise/backward FFT3d
 * \param   N   : integer pointer to size of FFT3d  
 * \param   din : complex input/output single precision data pointer 
 *****************************************************************************/
extern "C" void pw_fpga_fft3d_sp_(uint32_t direction, uint32_t *N, cmplx *din) {

  queue_setup();

  printf("In SP FPGA C function\n");
  // If fft size changes, need to rebuild program using another binary
  if(fft_size_changed == true){
    status = select_binary(N);
    checkError(status, "Failed to select binary as no relevant FFT3d binaries found in the directory!");

    init_program(N);
  }

  // setup device specific constructs 
  if(direction == 1){
    fftfpga_run_3d(0, N, din);
  }
  else{
    fftfpga_run_3d(1, N, din);
  }

  queue_cleanup();

}

extern "C" void pw_fpga_fft3d_dp_(uint32_t direction, uint32_t *N, cmplx *din) {

  queue_setup();

  printf("In DP FPGA C function\n");
  // If fft size changes, need to rebuild program using another binary
  if(fft_size_changed == true){
    status = select_binary(N);
    checkError(status, "Failed to select binary as no relevant FFT3d binaries found in the directory!");

    init_program(N);
  }

  // setup device specific constructs 
  if(direction == 1){
    fftfpga_run_3d(0, N, din);
  }
  else{
    fftfpga_run_3d(1, N, din);
  }

  queue_cleanup();

}
/******************************************************************************
 * \brief   Execute a single precision complex FFT3d
 * \author  Arjun Ramaswami
 * \param   inverse : boolean
 * \param   N       : integer pointer to size of FFT3d  
 * \param   din     : complex input/output single precision data pointer 
 *****************************************************************************/
static void fftfpga_run_3d(bool inverse, const uint32_t *N, cmplx *c_in) {

  int inverse_int = inverse;

  h_inData = (cmplx *)alignedMalloc(sizeof(cmplx) * N[0] * N[1] * N[2]);
  h_outData = (cmplx *)alignedMalloc(sizeof(cmplx) * N[0] * N[1] * N[2]);

  memcpy(h_inData, c_in, sizeof(cmplx) * N[0] * N[1] * N[2]);

  // Copy data from host to device
  status = clEnqueueWriteBuffer(queue6, d_inData, CL_TRUE, 0, sizeof(cmplx) * N[0] * N[1] * N[2], h_inData, 0, NULL, NULL);
  checkError(status, "Failed to copy data to device");

  status = clFinish(queue6);
  checkError(status, "failed to finish");
  // Can't pass bool to device, so convert it to int

  status = clSetKernelArg(fetch_kernel, 0, sizeof(cl_mem), (void *)&d_inData);
  checkError(status, "Failed to set kernel arg 0");
  status = clSetKernelArg(fft_kernel, 0, sizeof(cl_int), (void*)&inverse_int);
  checkError(status, "Failed to set kernel arg 1");
  status = clSetKernelArg(transpose_kernel, 0, sizeof(cl_mem), (void *)&d_outData);
  checkError(status, "Failed to set kernel arg 2");
  status = clSetKernelArg(fft_kernel_2, 0, sizeof(cl_int), (void*)&inverse_int);
  checkError(status, "Failed to set kernel arg 3");

  status = clEnqueueTask(queue1, fetch_kernel, 0, NULL, NULL);
  checkError(status, "Failed to launch fetch kernel");

  // Launch the fft kernel - we launch a single work item hence enqueue a task
  status = clEnqueueTask(queue2, fft_kernel, 0, NULL, NULL);
  checkError(status, "Failed to launch fft kernel");

  status = clEnqueueTask(queue3, transpose_kernel, 0, NULL, NULL);
  checkError(status, "Failed to launch transpose kernel");

  status = clEnqueueTask(queue4, fft_kernel_2, 0, NULL, NULL);
  checkError(status, "Failed to launch second fft kernel");

  status = clEnqueueTask(queue5, transpose_kernel_2, 0, NULL, NULL);
  checkError(status, "Failed to launch second transpose kernel");

  // Wait for all command queues to complete pending events
  status = clFinish(queue1);
  checkError(status, "failed to finish");
  status = clFinish(queue2);
  checkError(status, "failed to finish");
  status = clFinish(queue3);
  checkError(status, "failed to finish");
  status = clFinish(queue4);
  checkError(status, "failed to finish");
  status = clFinish(queue5);
  checkError(status, "failed to finish");

  // Copy results from device to host
  status = clEnqueueReadBuffer(queue3, d_outData, CL_TRUE, 0, sizeof(cmplx) * N[0] * N[1] * N[2], h_outData, 0, NULL, NULL);
  checkError(status, "Failed to read data from device");

  memcpy(c_in, h_outData, sizeof(cmplx) * N[0] * N[1] * N[2] );

  if (h_outData)
	alignedFree(h_outData);
  if (h_inData)
	alignedFree(h_inData);

}

bool getCwdDir(char * path, int len) {
  // Get path of executable.
  ssize_t n = readlink("/proc/self/exe", path, len -1 );
  //ssize_t n = readlink("/proc/self/exe", path, sizeof(path)/sizeof(path[0]) - 1);
  
  if(n == -1) {
    return false;
  }
  path[n] = 0;

  // Find the last '\' or '/' and terminate the path there; it is now
  // the directory containing the executable.
  size_t i;
  for(i = strlen(path) - 1; i > 0 && path[i] != '/' && path[i] != '\\'; --i);
  path[i] = '\0';

  //memcpy(cwd, path, n);
  return true;
}
/******************************************************************************
 * \brief  select an FPGA board binary from the default location
 * \author Arjun Ramaswami
 * \param  N - integer pointer to the size of the FFT3d
 * \retval true if no board binary found in the location
 *****************************************************************************/
static bool select_binary(const uint32_t *N){
    char cwd[300];
    if(!getCwdDir(cwd, sizeof(cwd)/sizeof(cwd[0]) )){
        printf("Problem finding executable\n");
    }

#ifdef __PW_FPGA_SP
    //printf("Selecting Single Precision binaries - (%d, %d, %d)\n", N[0], N[1], N[2]);
    switch(N[0]){
        case 16 :
          strcat(cwd, "/../../fpgabitstream/fft3d/synthesis/syn16/fft3d");
          printf("Selecting Single Precision binaries - (%d, %d, %d) in path %s\n", N[0], N[1], N[2], cwd);
          binary_file = getBoardBinaryFile(path, device);
        break;

        case 32 :
          strcat(cwd, "/../../fpgabitstream/fft3d/synthesis/syn32/fft3d");
          printf("Selecting Single Precision binaries - (%d, %d, %d) in path %s\n", N[0], N[1], N[2], cwd);
          binary_file = getBoardBinaryFile(path, device);
          break;
     
        case 64 :
          strcat(cwd, "/../../fpgabitstream/fft3d/synthesis/syn64/fft3d");
          printf("Selecting Single Precision binaries - (%d, %d, %d) in path %s\n", N[0], N[1], N[2], cwd);
          binary_file = getBoardBinaryFile(path, device);
          break;
        
        default:
          return true;
    }
#else
    //printf("Selecting emulation Double Precision binaries  - (%d, %d, %d)\n", N[0], N[1], N[2] );
    switch(N[0]){
        case 16 :
          strcat(cwd, "/../../fpgabitstream/fft3d/emulation_dp/emu16/fft3d");
          printf("Selecting Double Precision binaries - (%d, %d, %d) in path %s\n", N[0], N[1], N[2], cwd);
          binary_file = getBoardBinaryFile(cwd, device);
        break;

        case 32 :
          strcat(cwd, "/../../fpgabitstream/fft3d/emulation_dp/emu32/fft3d");
          printf("Selecting Double Precision binaries - (%d, %d, %d) in path %s\n", N[0], N[1], N[2], cwd);
          binary_file = getBoardBinaryFile(cwd, device);
          break;
     
        case 64 :
          strcat(cwd, "/../../fpgabitstream/fft3d/emulation_dp/emu64/fft3d");
          printf("Selecting Double Precision binaries - (%d, %d, %d) in path %s\n", N[0], N[1], N[2], cwd);
          binary_file = getBoardBinaryFile(cwd, device);
          break;
        
        default:
          return true;
    }
#endif
    
    return false;
}

/******************************************************************************
 * \brief   Initialize the OpenCL FPGA environment
 * \author  Arjun Ramaswami
 * \retval  true if error in initialization
 *****************************************************************************/
bool init() {
  cl_int status;

  // Get the OpenCL platform.
  platform = findPlatform("Intel(R) FPGA");
  if(platform == NULL) {
    printf("ERROR: Unable to find Intel(R) FPGA OpenCL platform\n");
    return true;
  }

  // Query the available OpenCL devices.
  scoped_array<cl_device_id> devices;
  cl_uint num_devices;

  devices.reset(getDevices(platform, CL_DEVICE_TYPE_ALL, &num_devices));

  // use the first device.
  device = devices[0];

  // Create the context.
  context = clCreateContext(NULL, 1, &device, &oclContextCallback, NULL, &status);
  checkError(status, "Failed to create context");

  return false;
}

/******************************************************************************
 * \brief   Initialize the program and its kernels
 * \author  Arjun Ramaswami
 * \param   N - integer pointer to the size of the FFT
 *****************************************************************************/
static void init_program(const uint32_t *N) {

  // Create the program.
  program = createProgramFromBinary(context, binary_file.c_str(), &device, 1);

  // Build the program that was just created.
  status = clBuildProgram(program, 0, NULL, "", NULL, NULL);
  checkError(status, "Failed to build program");

  // Create the kernel - name passed in here must match kernel name in the
  // original CL file, that was compiled into an AOCX file using the AOC tool
  fft_kernel = clCreateKernel(program, "fft3da", &status);
  checkError(status, "Failed to create fft3da kernel");
  fft_kernel_2 = clCreateKernel(program, "fft3db", &status);
  checkError(status, "Failed to create fft3db kernel");
  fetch_kernel = clCreateKernel(program, "fetch", &status);
  checkError(status, "Failed to create fetch kernel");
  transpose_kernel = clCreateKernel(program, "transpose", &status);
  checkError(status, "Failed to create transpose kernel");
  transpose_kernel_2 = clCreateKernel(program, "transpose3d", &status);
  checkError(status, "Failed to create transpose3d kernel");

  d_inData = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cmplx) * N[0] * N[1] * N[2], NULL, &status);
  checkError(status, "Failed to allocate input device buffer\n");
  d_outData = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cmplx) * N[0] * N[1] * N[2], NULL, &status);
  checkError(status, "Failed to allocate output device buffer\n");
}

/******************************************************************************
 * \brief   Create a command queue for each kernel
 * \author  Arjun Ramaswami
 *****************************************************************************/
static void queue_setup(){
  cl_int status;
  // Create one command queue for each kernel.
  queue1 = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
  checkError(status, "Failed to create command queue1");
  queue2 = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
  checkError(status, "Failed to create command queue2");
  queue3 = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
  checkError(status, "Failed to create command queue3");
  queue4 = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
  checkError(status, "Failed to create command queue4");
  queue5 = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
  checkError(status, "Failed to create command queue5");
  queue6 = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
  checkError(status, "Failed to create command queue6");

}

/******************************************************************************
 * \brief   Free resources allocated during initialization
 * \author  Arjun Ramaswami
 *****************************************************************************/
void cleanup(){

  if(context)
    clReleaseContext(context);

  if(program) 
    clReleaseProgram(program);

  if(fft_kernel) 
    clReleaseKernel(fft_kernel);  
  if(fft_kernel_2) 
    clReleaseKernel(fft_kernel_2);  
  if(fetch_kernel) 
    clReleaseKernel(fetch_kernel);  
  if(transpose_kernel) 
    clReleaseKernel(transpose_kernel);  
  if(transpose_kernel_2) 
    clReleaseKernel(transpose_kernel_2);  

  if (d_inData)
	clReleaseMemObject(d_inData);
  if (d_outData) 
	clReleaseMemObject(d_outData);
}

/******************************************************************************
 * \brief   Release all command queues
 * \author  Arjun Ramaswami
 *****************************************************************************/
void queue_cleanup() {
  if(queue1) 
    clReleaseCommandQueue(queue1);
  if(queue2) 
    clReleaseCommandQueue(queue2);
  if(queue3) 
    clReleaseCommandQueue(queue3);
  if(queue4) 
    clReleaseCommandQueue(queue4);
  if(queue5) 
    clReleaseCommandQueue(queue5);
  if(queue6) 
    clReleaseCommandQueue(queue6);
}

#endif
