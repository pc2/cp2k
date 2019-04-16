!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \note
!> This module contains routines necessary to operate on plane waves on INTEL
!> FPGAs using OpenCL. It depends at execution time on the board support
!> packages of the specific FPGA
!> \author Arjun Ramaswami
! **************************************************************************************************

MODULE pw_fpga
   USE ISO_C_BINDING,                   ONLY: C_BOOL,&
                                              C_INT,&
                                              C_LOC,&
                                              C_PTR
   USE kinds,                           ONLY: dp,&
                                              sp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pw_fpga_init, pw_fpga_finalize
   PUBLIC :: pw_fpga_check_bitstream
   PUBLIC :: pw_fpga_r3dc1d_3d, pw_fpga_c1dr3d_3d

   INTERFACE
! **************************************************************************************************
!> \brief Initialize FPGA
!> \retval status if the routine failed or not
! **************************************************************************************************
      FUNCTION pw_fpga_initialize() RESULT(stat) &
         BIND(C, name="pw_fpga_initialize_")
         IMPORT
         INTEGER(KIND=C_INT)                    :: stat
      END FUNCTION pw_fpga_initialize

! **************************************************************************************************
!> \brief Destroy FPGA
! **************************************************************************************************
      SUBROUTINE pw_fpga_final() &
         BIND(C, name="pw_fpga_final_")
      END SUBROUTINE pw_fpga_final

   END INTERFACE

   INTERFACE
! **************************************************************************************************
!> \brief Check whether an fpga bitstream for the given FFT3d size is present
!> \param npts ...
!> \return res ...
! **************************************************************************************************
      FUNCTION pw_fpga_check_bitstream(n) RESULT(res) &
         BIND(C, name="pw_fpga_check_bitstream_")
         IMPORT
         INTEGER(KIND=C_INT), DIMENSION(*), &
            INTENT(IN)                            :: n
         LOGICAL(KIND=C_BOOL)                     :: res
      END FUNCTION pw_fpga_check_bitstream

   END INTERFACE

   INTERFACE
! **************************************************************************************************
!> \brief FFT3d using FPGA
!> \param dir - direction of FFT3d ...
!> \param npts - dimensions of FFT3d
!> \param single precision c_in...
! **************************************************************************************************
      SUBROUTINE pw_fpga_fft3d_sp(dir, n, c_in_sp) &
         BIND(C, name="pw_fpga_fft3d_sp_")
         IMPORT
         INTEGER(KIND=C_INT), VALUE, INTENT(IN)   :: dir
         INTEGER(KIND=C_INT), DIMENSION(*), &
            INTENT(IN)                            :: n
         TYPE(C_PTR), VALUE                       :: c_in_sp
      END SUBROUTINE pw_fpga_fft3d_sp
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Allocates resources on the fpga device
! **************************************************************************************************
   SUBROUTINE pw_fpga_init()
#if defined (__PW_FPGA)
      INTEGER :: stat
      stat = pw_fpga_initialize()
      IF (stat /= 0) &
         CPABORT("pw_fpga_init: failed")
#endif
   END SUBROUTINE pw_fpga_init

! **************************************************************************************************
!> \brief Releases resources on the fpga device
! **************************************************************************************************
   SUBROUTINE pw_fpga_finalize()
#if defined (__PW_FPGA)
      CALL pw_fpga_final()
#endif
   END SUBROUTINE pw_fpga_finalize

! **************************************************************************************************
!> \brief perform an in-place fft3d on the FPGA
!> \param n ...
!> \param c_out  ...
! **************************************************************************************************
   SUBROUTINE pw_fpga_r3dc1d_3d(n, c_out)
      COMPLEX(KIND=dp), DIMENSION(:, :, :), &
         INTENT(INOUT)                               :: c_out
      INTEGER, DIMENSION(:), INTENT(IN)                 :: n

#if ! defined (__PW_FPGA)
      MARK_USED(c_out)
      MARK_USED(n)
#else
      COMPLEX, DIMENSION(:, :, :), POINTER              :: c_in_sp
      INTEGER                                           :: handle3

      CHARACTER(len=*), PARAMETER :: routineX = 'fw_fft_fpga_r3dc1d'

      ALLOCATE (c_in_sp(n(1), n(2), n(3)))
      ! pointer to single precision complex array
      c_in_sp(:, :, :) = CMPLX(c_out(:, :, :), KIND=sp)

      CALL timeset(routineX, handle3)
      CALL pw_fpga_fft3d_sp(+1, n, c_loc(c_in_sp))
      CALL timestop(handle3)

      ! typecast sp back to dp
      c_out(:, :, :) = CMPLX(c_in_sp(:, :, :), KIND=dp)

      DEALLOCATE (c_in_sp)
#endif
   END SUBROUTINE

! **************************************************************************************************
!> \brief perform an in-place inverse fft3d on the FPGA
!> \param n ...
!> \param c_out  ...
! **************************************************************************************************
   SUBROUTINE pw_fpga_c1dr3d_3d(n, c_out)
      COMPLEX(KIND=dp), DIMENSION(:, :, :), &
         INTENT(INOUT)                               :: c_out
      INTEGER, DIMENSION(:), INTENT(IN)                 :: n

#if ! defined (__PW_FPGA)
      MARK_USED(c_out)
      MARK_USED(n)

#else
      COMPLEX, DIMENSION(:, :, :), POINTER              :: c_in_sp
      INTEGER                                           :: handle3

      CHARACTER(len=*), PARAMETER :: routineX = 'bw_fft_fpga_c1dr3d'

      ALLOCATE (c_in_sp(n(1), n(2), n(3)))
      ! pointer to single precision complex array
      c_in_sp(:, :, :) = CMPLX(c_out(:, :, :), KIND=sp)

      CALL timeset(routineX, handle3)
      CALL pw_fpga_fft3d_sp(-1, n, c_loc(c_in_sp))
      CALL timestop(handle3)

      ! typecast sp back to dp
      c_out(:, :, :) = CMPLX(c_in_sp(:, :, :), KIND=dp)

      DEALLOCATE (c_in_sp)
#endif
   END SUBROUTINE

END MODULE pw_fpga

