# FPGA

## Steps to generate bitstream for the given FFT3d kernel code 

### Set the size of the FFT3d 
 - The *fft3d_config.h* file contains the parameter *LOGN* to specify the log
    of the size of the FFT.
 - For example, LOGN of 6 represents a 64^3 FFT3d, considering 2^6 = 64
   signifies a dimension of the 3d FFT.
 - Not every kernel size may fit the FPGA device for the given kernel design.

### Synthesize kernels for a particular board 
	aoc -fp-relaxed --board $(BOARDNAME) -g -v fft3d.cl -o fft3d.aocx

### Copy the bitstream generated to particular path
 - The bitstream generated can be found as the *.aocx* file. 
 - Copy this file to the path *~/cp2k/bitstream/fft3d/synthesis/syn??/*. 
    - The sizes supported are 16^3, 32^3, 64^3.
    - Therefore, a 16^3 FFT3d file should be copied to syn16 folder.
 - If required size is different from the default options, 
    - modify the switch case in the *fft_fpga.cpp* to include the required size
      and the path to the location of the bitstream.


## NOTE
 - Kernel code has only been tested with single precision floating point
   numbers.
 - Kernel code works for uniform FFT3d sizes.
