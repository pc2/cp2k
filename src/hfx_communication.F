!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for data exchange between MPI processes
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
MODULE hfx_communication
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_block_p,&
                                              dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop,&
                                              dbcsr_iterator_type,&
                                              dbcsr_get_stored_coordinates,&
                                              dbcsr_put_block,&
                                              dbcsr_p_type,&
                                              dbcsr_type,&
                                              dbcsr_finalize
   USE hfx_contract_block,              ONLY: contract_block
   USE hfx_screening_methods,           ONLY: update_pmax_mat
   USE hfx_types,                       ONLY: hfx_2D_map,&
                                              hfx_basis_type,&
                                              hfx_type, &
                                              hfx_libheric_basis_type
   USE kinds,                           ONLY: dp,&
                                              int_8,&
                                              int_1
   USE libheric_interface,              ONLY: atom_quartet_type
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_request_type,&
                                              mp_waitall
   USE orbital_pointers,                ONLY: nso
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: get_full_density, &
             distribute_ks_matrix, &
             scale_and_add_fock_to_ks_matrix, &
             get_atomic_block_maps, &
             hfx_libheric_get_pmatrix_from_dbcsr_matrix, &
             hfx_libheric_initialize_kmatrix, &
             hfx_libheric_reconstruct_dbcsr_kmatrix, &
             hfx_libheric_free_sparse_matrix, &
             hfx_libheric_pmatrix_verify, &
             hfx_libheric_sparse_matrix, &
             hfx_libheric_get_column_index_map, &
             hfx_libheric_digestion_state, &
             hfx_libheric_init_digestion_state, &
             hfx_libheric_digestion_state_set_angular_momentum, &
             hfx_libheric_digestion_state_flush, &
             hfx_libheric_digest_atom_quartet, &
             hfx_libheric_free_digestion_state
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_communication'

!***
   TYPE hfx_libheric_block_row_info
      INTEGER, DIMENSION(:), ALLOCATABLE              :: key
      INTEGER(int_8), DIMENSION(:, :, :), ALLOCATABLE :: offset
      INTEGER, DIMENSION(:), ALLOCATABLE              :: origin
   END TYPE

   TYPE hfx_libheric_sparse_matrix
      REAL(dp), DIMENSION(:), POINTER                         :: data
      TYPE(hfx_libheric_block_row_info), DIMENSION(:), &
         ALLOCATABLE                                          :: rows
      INTEGER, DIMENSION(:, :, :, :, :, :, :, :), POINTER     :: shellset_offset
      INTEGER, DIMENSION(:, :, :, :), ALLOCATABLE             :: screening_info_size, atom_ll_block_size
      INTEGER(int_8)                                          :: size 
      REAL(dp), DIMENSION(:), POINTER                         :: zero_block
      INTEGER                                                 :: invalid_col = 0 !! Do not modify; should be considered constant !!
   END TYPE

   TYPE hfx_libheric_queue
      INTEGER, DIMENSION(:), ALLOCATABLE :: queue
      INTEGER                            :: read_ptr, write_ptr, count, size
   END TYPE

   TYPE hfx_libheric_kcache_entry
      INTEGER(int_8) :: hash 
      INTEGER        :: iatom, katom_col, la, lc, offset, atom_block_size
      INTEGER        :: left_block, right_block 
   END TYPE

   TYPE hfx_libheric_sparse_matrix_block
      INTEGER                                 :: iatom, katom, ikind, kkind, la, lc, atom_block_size, n_ac
      INTEGER, DIMENSION(:), POINTER          :: cols
      INTEGER, DIMENSION(:, :, :, :), POINTER :: shellset_offset 
      INTEGER(int_8)                          :: hash, hash_postfix
      LOGICAL                                 :: required, present
      REAL(dp), DIMENSION(:), POINTER         :: data 
   END TYPE

   TYPE hfx_libheric_digestion_state
      REAL(dp), DIMENSION(:), POINTER                            :: data
      TYPE(hfx_libheric_kcache_entry), DIMENSION(:), ALLOCATABLE :: kcache
      INTEGER                                                    :: first_block, last_block
      INTEGER                                                    :: null_block = -1 !! Do not modify; should be considered constant !!
      TYPE(hfx_libheric_queue)                                   :: free_kblocks
      !! p(1) == p_ac, p(2) == p_ad, p(3) == p_bc, p(4) == p_bd -- same for k(1:4)
      TYPE(hfx_libheric_sparse_matrix_block)                     :: p(4), k(4)
      !! Only necessary because contract_block zeros out k buffers before contraction
      REAL(dp), DIMENSION(:), POINTER                            :: k_ac, k_ad, k_bc, k_bd
      INTEGER                                                    :: k_tmp_offset
      INTEGER                                                    :: cache_block_size
      INTEGER, DIMENSION(:), POINTER                             :: pa_cols, pb_cols, ka_cols, kb_cols
      INTEGER                                                    :: prev_iatom, prev_jatom
      INTEGER                                                    :: la, lb, lc, ld
      INTEGER                                                    :: na, nb, nc, nd
      INTEGER                                                    :: lock_id
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief - Collects full density matrix from all CPUs
!> \param para_env ...
!> \param full_density The full Density matrix
!> \param rho Distributed density
!> \param number_of_p_entries Maximal buffer size
!> \param block_offset ...
!> \param kind_of ...
!> \param basis_parameter ...
!> \param get_max_vals_spin ...
!> \param rho_beta ...
!> \param antisymmetric ...
!> \par History
!>      11.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      - Communication with left/right node only
!>        added a mpi_sync before and after the ring of isendrecv. This *speed up* the
!>        communication, and might protect against idle neighbors flooding a busy node
!>        with messages [Joost]
! **************************************************************************************************
   SUBROUTINE get_full_density(para_env, full_density, rho, number_of_p_entries, &
                               block_offset, kind_of, basis_parameter, &
                               get_max_vals_spin, rho_beta, antisymmetric)

      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(dp), DIMENSION(:)                             :: full_density
      TYPE(dbcsr_type), POINTER                          :: rho
      INTEGER, INTENT(IN)                                :: number_of_p_entries
      INTEGER, DIMENSION(:), POINTER                     :: block_offset
      INTEGER                                            :: kind_of(*)
      TYPE(hfx_basis_type), DIMENSION(:), POINTER        :: basis_parameter
      LOGICAL, INTENT(IN)                                :: get_max_vals_spin
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: rho_beta
      LOGICAL, INTENT(IN)                                :: antisymmetric

      INTEGER :: block_size, data_from, dest, i, iatom, icpu, ikind, iset, jatom, jkind, jset, &
         mepos, ncpu, nseta, nsetb, pa, pa1, pb, pb1, source, source_cpu
      INTEGER, DIMENSION(:), POINTER                     :: nsgfa, nsgfb
      LOGICAL                                            :: found
      REAL(dp)                                           :: symmfac
      REAL(dp), DIMENSION(:), POINTER                    :: recbuffer, sendbuffer, swapbuffer
      REAL(dp), DIMENSION(:, :), POINTER                 :: sparse_block, sparse_block_beta
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_request_type), DIMENSION(2)                :: req

      full_density = 0.0_dp
      ALLOCATE (sendbuffer(number_of_p_entries))
      ALLOCATE (recbuffer(number_of_p_entries))
      sendbuffer = 0.0_dp
      recbuffer = 0.0_dp

      i = 1
      CALL dbcsr_iterator_start(iter, rho, shared=.FALSE.)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block)
         ! the resulting vector will be only the upper triangle.
         ! in case of antisymmetry take care to change signs if a lower block gets copied
         symmfac = 1.0_dp
         IF (antisymmetric .AND. iatom > jatom) symmfac = -1.0_dp
         ikind = kind_of(iatom)
         nseta = basis_parameter(ikind)%nset
         nsgfa => basis_parameter(ikind)%nsgf
         jkind = kind_of(jatom)
         nsetb = basis_parameter(jkind)%nset
         nsgfb => basis_parameter(jkind)%nsgf
         IF (get_max_vals_spin) THEN
            CALL dbcsr_get_block_p(rho_beta, &
                                   row=iatom, col=jatom, BLOCK=sparse_block_beta, found=found)
            pa = 0
            DO iset = 1, nseta
               pb = 0
               DO jset = 1, nsetb
                  DO pa1 = pa + 1, pa + nsgfa(iset)
                     DO pb1 = pb + 1, pb + nsgfb(jset)
                        sendbuffer(i) = MAX(ABS(sparse_block(pa1, pb1)), ABS(sparse_block_beta(pa1, pb1)))
                        i = i + 1
                     END DO
                  END DO
                  pb = pb + nsgfb(jset)
               END DO
               pa = pa + nsgfa(iset)
            END DO
         ELSE
            pa = 0
            DO iset = 1, nseta
               pb = 0
               DO jset = 1, nsetb
                  DO pa1 = pa + 1, pa + nsgfa(iset)
                     DO pb1 = pb + 1, pb + nsgfb(jset)
                        sendbuffer(i) = sparse_block(pa1, pb1)*symmfac
                        i = i + 1
                     END DO
                  END DO
                  pb = pb + nsgfb(jset)
               END DO
               pa = pa + nsgfa(iset)
            END DO
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      ! sync before/after a ring of isendrecv
      CALL para_env%sync()
      ncpu = para_env%num_pe
      mepos = para_env%mepos
      dest = MODULO(mepos + 1, ncpu)
      source = MODULO(mepos - 1, ncpu)
      DO icpu = 0, ncpu - 1
         IF (icpu /= ncpu - 1) THEN
            CALL para_env%isendrecv(sendbuffer, dest, recbuffer, source, &
                                    req(1), req(2), 13)
         END IF
         data_from = MODULO(mepos - icpu, ncpu)
         source_cpu = MODULO(data_from, ncpu) + 1
         block_size = block_offset(source_cpu + 1) - block_offset(source_cpu)
         full_density(block_offset(source_cpu):block_offset(source_cpu) + block_size - 1) = sendbuffer(1:block_size)

         IF (icpu /= ncpu - 1) THEN
            CALL mp_waitall(req)
         END IF
         swapbuffer => sendbuffer
         sendbuffer => recbuffer
         recbuffer => swapbuffer
      END DO
      DEALLOCATE (sendbuffer, recbuffer)
      ! sync before/after a ring of isendrecv
      CALL para_env%sync()

   END SUBROUTINE get_full_density

! **************************************************************************************************
!> \brief - Distributes the local full Kohn-Sham matrix to all CPUS
!> \param para_env ...
!> \param full_ks The full Kohn-Sham matrix
!> \param ks_matrix Distributed Kohn-Sham matrix
!> \param number_of_p_entries Maximal buffer size
!> \param block_offset ...
!> \param kind_of ...
!> \param basis_parameter ...
!> \param off_diag_fac ...
!> \param diag_fac ...
!> \par History
!>      11.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      - Communication with left/right node only
! **************************************************************************************************
   SUBROUTINE distribute_ks_matrix(para_env, full_ks, ks_matrix, number_of_p_entries, &
                                   block_offset, kind_of, basis_parameter, &
                                   off_diag_fac, diag_fac)

      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(dp), DIMENSION(:)                             :: full_ks
      TYPE(dbcsr_type), POINTER                          :: ks_matrix
      INTEGER, INTENT(IN)                                :: number_of_p_entries
      INTEGER, DIMENSION(:), POINTER                     :: block_offset
      INTEGER                                            :: kind_of(*)
      TYPE(hfx_basis_type), DIMENSION(:), POINTER        :: basis_parameter
      REAL(dp), INTENT(IN), OPTIONAL                     :: off_diag_fac, diag_fac

      INTEGER :: block_size, data_to, dest, dest_cpu, i, iatom, icpu, ikind, iset, jatom, jkind, &
         jset, mepos, ncpu, nseta, nsetb, pa, pa1, pb, pb1, source
      INTEGER, DIMENSION(:), POINTER                     :: nsgfa, nsgfb
      REAL(dp)                                           :: my_fac, myd_fac
      REAL(dp), DIMENSION(:), POINTER                    :: recbuffer, sendbuffer, swapbuffer
      REAL(dp), DIMENSION(:, :), POINTER                 :: sparse_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_request_type), DIMENSION(2)                :: req

      my_fac = 1.0_dp; myd_fac = 1.0_dp
      IF (PRESENT(off_diag_fac)) my_fac = off_diag_fac
      IF (PRESENT(diag_fac)) myd_fac = diag_fac

      ALLOCATE (sendbuffer(number_of_p_entries))
      sendbuffer = 0.0_dp
      ALLOCATE (recbuffer(number_of_p_entries))
      recbuffer = 0.0_dp

      ncpu = para_env%num_pe
      mepos = para_env%mepos
      dest = MODULO(mepos + 1, ncpu)
      source = MODULO(mepos - 1, ncpu)

      ! sync before/after a ring of isendrecv
      CALL para_env%sync()
      DO icpu = 1, ncpu
         i = 1
         data_to = mepos - icpu
         dest_cpu = MODULO(data_to, ncpu) + 1
         block_size = block_offset(dest_cpu + 1) - block_offset(dest_cpu)
       sendbuffer(1:block_size) = sendbuffer(1:block_size) + full_ks(block_offset(dest_cpu):block_offset(dest_cpu) + block_size - 1)
         IF (icpu == ncpu) EXIT
         CALL para_env%isendrecv(sendbuffer, dest, recbuffer, source, &
                                 req(1), req(2), 13)

         CALL mp_waitall(req)
         swapbuffer => sendbuffer
         sendbuffer => recbuffer
         recbuffer => swapbuffer
      END DO
      ! sync before/after a ring of isendrecv
      CALL para_env%sync()

      i = 1
      CALL dbcsr_iterator_start(iter, ks_matrix, shared=.FALSE.)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block)

         ikind = kind_of(iatom)
         nseta = basis_parameter(ikind)%nset
         nsgfa => basis_parameter(ikind)%nsgf
         jkind = kind_of(jatom)
         nsetb = basis_parameter(jkind)%nset
         nsgfb => basis_parameter(jkind)%nsgf
         pa = 0
         DO iset = 1, nseta
            pb = 0
            DO jset = 1, nsetb
               DO pa1 = pa + 1, pa + nsgfa(iset)
                  DO pb1 = pb + 1, pb + nsgfb(jset)
                     IF (iatom == jatom .AND. pa1 == pb1) THEN
                        sparse_block(pa1, pb1) = sendbuffer(i)*myd_fac + sparse_block(pa1, pb1)
                     ELSE
                        sparse_block(pa1, pb1) = sendbuffer(i)*my_fac + sparse_block(pa1, pb1)
                     END IF
                     i = i + 1
                  END DO
               END DO
               pb = pb + nsgfb(jset)
            END DO
            pa = pa + nsgfa(iset)
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      DEALLOCATE (sendbuffer, recbuffer)

   END SUBROUTINE distribute_ks_matrix

! **************************************************************************************************
!> \brief - Distributes the local full Kohn-Sham matrix to all CPUS. Is called in
!>        case of adiabatic rescaling. This is just a refactored version of
!>        distribute_ks_matrix
!> \param para_env ...
!> \param qs_env ...
!> \param ks_matrix Distributed Kohn-Sham matrix
!> \param irep ...
!> \param scaling_factor ...
!> \par History
!>      11.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      - Communication with left/right node only
! **************************************************************************************************
   SUBROUTINE scale_and_add_fock_to_ks_matrix(para_env, qs_env, ks_matrix, irep, &
                                              scaling_factor)

      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      INTEGER, INTENT(IN)                                :: irep
      REAL(dp), INTENT(IN)                               :: scaling_factor

      INTEGER                                            :: iatom, ikind, img, natom, nimages, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of, last_sgf_global
      REAL(dp), DIMENSION(:, :), POINTER                 :: full_ks
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_basis_type), DIMENSION(:), POINTER        :: basis_parameter
      TYPE(hfx_type), POINTER                            :: actual_x_data
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

!! All shared data is saved in i_thread = 1!

      NULLIFY (dft_control)
      actual_x_data => qs_env%x_data(irep, 1)
      basis_parameter => actual_x_data%basis_parameter

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, &
                      dft_control=dft_control)

      nspins = dft_control%nspins
      nimages = dft_control%nimages
      CPASSERT(nimages == 1)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)

      natom = SIZE(particle_set, 1)
      ALLOCATE (last_sgf_global(0:natom))
      last_sgf_global(0) = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         last_sgf_global(iatom) = last_sgf_global(iatom - 1) + basis_parameter(ikind)%nsgf_total
      END DO
      full_ks => actual_x_data%full_ks_alpha
      IF (scaling_factor /= 1.0_dp) THEN
         full_ks = full_ks*scaling_factor
      END IF
      DO img = 1, nimages
         CALL distribute_ks_matrix(para_env, full_ks(:, img), ks_matrix(1, img)%matrix, actual_x_data%number_of_p_entries, &
                                   actual_x_data%block_offset, kind_of, basis_parameter, &
                                   off_diag_fac=0.5_dp)
      END DO
      DEALLOCATE (actual_x_data%full_ks_alpha)

      IF (nspins == 2) THEN
         full_ks => actual_x_data%full_ks_beta
         IF (scaling_factor /= 1.0_dp) THEN
            full_ks = full_ks*scaling_factor
         END IF
         DO img = 1, nimages
            CALL distribute_ks_matrix(para_env, full_ks(:, img), ks_matrix(2, img)%matrix, actual_x_data%number_of_p_entries, &
                                      actual_x_data%block_offset, kind_of, basis_parameter, &
                                      off_diag_fac=0.5_dp)
         END DO
         DEALLOCATE (actual_x_data%full_ks_beta)
      END IF

      DEALLOCATE (last_sgf_global)

   END SUBROUTINE scale_and_add_fock_to_ks_matrix

! **************************************************************************************************
!> \brief Given a 2d index pair, this function returns a 1d index pair for
!>        a symmetric upper triangle NxN matrix
!>        The compiler should inline this function, therefore it appears in
!>        several modules
!> \param i 2d index
!> \param j 2d index
!> \param N matrix size
!> \return ...
!> \par History
!>      03.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   PURE FUNCTION get_1D_idx(i, j, N)
      INTEGER, INTENT(IN)                                :: i, j
      INTEGER(int_8), INTENT(IN)                         :: N
      INTEGER(int_8)                                     :: get_1D_idx

      INTEGER(int_8)                                     :: min_ij

      min_ij = MIN(i, j)
      get_1D_idx = min_ij*N + MAX(i, j) - (min_ij - 1)*min_ij/2 - N

   END FUNCTION get_1D_idx

! **************************************************************************************************
!> \brief create a several maps array that reflects the ks matrix sparsity
!> \param matrix ...
!> \param basis_parameter ...
!> \param kind_of ...
!> \param is_assoc_atomic_block ...
!> \param number_of_p_entries ...
!> \param para_env ...
!> \param atomic_block_offset ...
!> \param set_offset ...
!> \param block_offset ...
!> \param map_atoms_to_cpus ...
!> \param nkind ...
!> \par History
!>      11.2007 refactored [Joost VandeVondele]
!>      07.2009 add new maps
!> \author Manuel Guidon
!> \notes
!>      is_assoc_atomic_block returns the mpi rank + 1 for associated blocks,
!>      zero for unassiated blocks
! **************************************************************************************************
   SUBROUTINE get_atomic_block_maps(matrix, basis_parameter, kind_of, &
                                    is_assoc_atomic_block, number_of_p_entries, &
                                    para_env, atomic_block_offset, set_offset, &
                                    block_offset, map_atoms_to_cpus, nkind)

      TYPE(dbcsr_type), POINTER                          :: matrix
      TYPE(hfx_basis_type), DIMENSION(:)                 :: basis_parameter
      INTEGER, DIMENSION(:)                              :: kind_of
      INTEGER, DIMENSION(:, :), INTENT(OUT)              :: is_assoc_atomic_block
      INTEGER, INTENT(OUT)                               :: number_of_p_entries
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, DIMENSION(:, :), POINTER                  :: atomic_block_offset
      INTEGER, DIMENSION(:, :, :, :), POINTER            :: set_offset
      INTEGER, DIMENSION(:), POINTER                     :: block_offset
      TYPE(hfx_2D_map), DIMENSION(:), POINTER            :: map_atoms_to_cpus
      INTEGER                                            :: nkind

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_atomic_block_maps'

      INTEGER :: handle, iatom, ibuf, icpu, ikind, ilist, iset, itask, jatom, jkind, jset, natom, &
         ncpu, nseta, nsetb, number_of_p_blocks, offset, tmp(2)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: buffer_in, buffer_out, counter, rcount, &
                                                            rdispl
      INTEGER, DIMENSION(:), POINTER                     :: iatom_list, jatom_list, nsgfa, nsgfb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sparse_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      is_assoc_atomic_block = 0
      number_of_p_entries = 0
      number_of_p_blocks = 0

      !
      ! count number_of_p_blocks and number_of_p_entries
      !
      CALL dbcsr_iterator_start(iter, matrix, shared=.FALSE.)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block)
         ikind = kind_of(iatom)
         jkind = kind_of(jatom)
         number_of_p_blocks = number_of_p_blocks + 1
         number_of_p_entries = number_of_p_entries + &
                               basis_parameter(ikind)%nsgf_total*basis_parameter(jkind)%nsgf_total
      END DO
      CALL dbcsr_iterator_stop(iter)

      tmp = [number_of_p_entries, number_of_p_blocks]
      CALL para_env%max(tmp)
      number_of_p_entries = tmp(1)
      number_of_p_blocks = tmp(2)
      !
      ! send this info around, so we can construct is_assoc_atomic_block
      ! pack all data in buffers and use allgatherv
      !
      ALLOCATE (buffer_in(3*number_of_p_blocks))
      ALLOCATE (buffer_out(3*number_of_p_blocks*para_env%num_pe))
      ALLOCATE (rcount(para_env%num_pe), rdispl(para_env%num_pe))

      buffer_in = 0
      ibuf = 0

      CALL dbcsr_iterator_start(iter, matrix, shared=.FALSE.)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block)

         buffer_in(ibuf + 1) = iatom
         buffer_in(ibuf + 2) = jatom
         buffer_in(ibuf + 3) = para_env%mepos + 1
         ibuf = ibuf + 3
      END DO
      CALL dbcsr_iterator_stop(iter)

      rcount = SIZE(buffer_in)
      rdispl(1) = 0
      DO icpu = 2, para_env%num_pe
         rdispl(icpu) = rdispl(icpu - 1) + rcount(icpu - 1)
      END DO
      CALL para_env%allgatherv(buffer_in, buffer_out, rcount, rdispl)

      DO ibuf = 0, number_of_p_blocks*para_env%num_pe*3 - 3, 3
         itask = buffer_out(ibuf + 3)
         ! buffer_out can be 0 if buffer_in contained less elements than the max number of atom pairs
         ! is_assoc_atomic_block is a map for atom pairs to a processor (assumes symmetry, i,j on the ame as j,i)
         IF (itask /= 0) THEN
            iatom = buffer_out(ibuf + 1)
            jatom = buffer_out(ibuf + 2)
            is_assoc_atomic_block(iatom, jatom) = itask
            is_assoc_atomic_block(jatom, iatom) = itask
         END IF
      END DO

      IF (ASSOCIATED(map_atoms_to_cpus)) THEN
         DO icpu = 1, para_env%num_pe
            DEALLOCATE (map_atoms_to_cpus(icpu)%iatom_list)
            DEALLOCATE (map_atoms_to_cpus(icpu)%jatom_list)
         END DO
         DEALLOCATE (map_atoms_to_cpus)
      END IF

      natom = SIZE(is_assoc_atomic_block, 1)
      ALLOCATE (map_atoms_to_cpus(para_env%num_pe))
      ALLOCATE (counter(para_env%num_pe))
      counter = 0

      DO iatom = 1, natom
         DO jatom = iatom, natom
            icpu = is_assoc_atomic_block(jatom, iatom)
            IF (icpu > 0) counter(icpu) = counter(icpu) + 1
         END DO
      END DO
      DO icpu = 1, para_env%num_pe
         ALLOCATE (map_atoms_to_cpus(icpu)%iatom_list(counter(icpu)))
         ALLOCATE (map_atoms_to_cpus(icpu)%jatom_list(counter(icpu)))
      END DO
      counter = 0
      DO iatom = 1, natom
         DO jatom = iatom, natom
            icpu = is_assoc_atomic_block(jatom, iatom)
            IF (icpu > 0) THEN
               counter(icpu) = counter(icpu) + 1
               map_atoms_to_cpus(icpu)%jatom_list(counter(icpu)) = jatom
               map_atoms_to_cpus(icpu)%iatom_list(counter(icpu)) = iatom
            END IF
         END DO
      END DO

      DEALLOCATE (counter)

      ncpu = para_env%num_pe
      offset = 1
      atomic_block_offset = 0
      block_offset = 0
      DO icpu = 1, ncpu
         iatom_list => map_atoms_to_cpus(icpu)%iatom_list
         jatom_list => map_atoms_to_cpus(icpu)%jatom_list
         block_offset(icpu) = offset
         DO ilist = 1, SIZE(iatom_list)
            iatom = iatom_list(ilist)
            ikind = kind_of(iatom)
            jatom = jatom_list(ilist)
            jkind = kind_of(jatom)
            atomic_block_offset(iatom, jatom) = offset
            atomic_block_offset(jatom, iatom) = offset
            offset = offset + basis_parameter(ikind)%nsgf_total*basis_parameter(jkind)%nsgf_total
         END DO
      END DO
      block_offset(ncpu + 1) = offset
      set_offset = 0
      DO ikind = 1, nkind
         nseta = basis_parameter(ikind)%nset
         nsgfa => basis_parameter(ikind)%nsgf
         DO jkind = 1, nkind
            nsetb = basis_parameter(jkind)%nset
            nsgfb => basis_parameter(jkind)%nsgf
            offset = 1
            DO iset = 1, nseta
               DO jset = 1, nsetb
                  set_offset(jset, iset, jkind, ikind) = offset
                  offset = offset + nsgfa(iset)*nsgfb(jset)
               END DO
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE get_atomic_block_maps
   
! **************************************************************************************************
!> \brief initialize index map into an atom block (the shellset_offset), valid also if screening data
!>        is stored along side the matrix values (relevant for pmatrix screening)
!> \param max_shl - maximum number of shells in any shellset
!> \param max_set - maximum number of shellsets in any basis set
!> \param nkind - number of unique basis sets 
!> \param l_max - maximum angular momentum
!> \param basis_parameter_l - all hfx_libheric basis sets
!> \param shellset_offset - index map from (ll-pair, kk-pair, set-set-pair, shl-shl-pair) 
!>                          to linear index into atom block
!> \param add_screening_info - if .TRUE., the linear index is adjusted according to number of 
!>                             screening values stored in front of the actual data 
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_initialize_shellset_offset(max_shl, max_set, nkind, l_max, basis_parameter_l, &
                                                      shellset_offset, add_screening_info)
      INTEGER, INTENT(IN)                       :: max_shl, max_set, nkind, l_max
      TYPE(hfx_libheric_basis_type), & 
          DIMENSION(nkind, 0:l_max), INTENT(IN)     :: basis_parameter_l
      INTEGER, DIMENSION(:, :, :, :, :, :, :, :), &
          POINTER, INTENT(OUT)                      :: shellset_offset
      LOGICAL, INTENT(IN)                           :: add_screening_info

      INTEGER                                       :: la, ikind, iset, ishl, &
                                                       lb, jkind, jset, jshl, &
                                                       shellset_off

      ALLOCATE(shellset_offset(max_shl, max_shl, max_set, max_set, nkind, nkind, 0:l_max, 0:l_max), source=-1)
      
      DO la = 0, l_max
      DO lb = 0, l_max
         DO ikind = 1, nkind
         DO jkind = 1, nkind
            ASSOCIATE(bas_la => basis_parameter_l(ikind, la), bas_lb => basis_parameter_l(jkind, lb))
            !! screening info is stored at the beginning of a block 
            shellset_off = MERGE(1 + bas_la%nset * bas_lb%nset, 0, add_screening_info)
            DO iset = 1, bas_la%nset
            DO jset = 1, bas_lb%nset
               DO ishl = 1, bas_la%nshell(iset)
               DO jshl = 1, bas_lb%nshell(jset)
                  shellset_offset(jshl, ishl, jset, iset, jkind, ikind, lb, la) = shellset_off
                  shellset_off = shellset_off + nso(la) * nso(lb)
               END DO
               END DO
            END DO
            END DO
            END ASSOCIATE
         END DO
         END DO
      END DO
      END DO

   END SUBROUTINE hfx_libheric_initialize_shellset_offset 

! **************************************************************************************************
!> \brief gather distributed dbcsr pmatrix and store in sparse format 
!> \param dbcsr_matrix - local part of distributed dbcsr pmatrix
!> \param para_env - MPI environment
!> \param kind_of - map of atom index to atom kind
!> \param natom - number of atoms in the simulated system
!> \param nkind - number of unique basis sets 
!> \param l_max - maximum angular momentum
!> \param max_set - maximum number of shellsets in any basis set
!> \param max_shl - maximum number of shells in any shellset
!> \param basis_parameter - all hfx basis sets
!> \param basis_parameter_l - all hfx_libheric basis sets
!> \param p_matrix - sparse pmatrix replicated on all MPI-ranks
!> \param add_screening_info - if .TRUE., add additional screening values in front of pmatrix blocks
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_get_pmatrix_from_dbcsr_matrix(dbcsr_matrix, para_env, kind_of, &
                                                         natom, nkind, l_max, max_set, max_shl, &
                                                         basis_parameter, basis_parameter_l, &
                                                         p_matrix, add_screening_info)
      TYPE(dbcsr_type), POINTER                          :: dbcsr_matrix
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of
      INTEGER, INTENT(IN)                                :: natom, nkind, l_max, max_set, max_shl
      TYPE(hfx_basis_type), & 
          DIMENSION(nkind), INTENT(IN)                   :: basis_parameter
      TYPE(hfx_libheric_basis_type), & 
          DIMENSION(nkind, 0:l_max), INTENT(IN)          :: basis_parameter_l
      TYPE(hfx_libheric_sparse_matrix), INTENT(OUT)      :: p_matrix
      LOGICAL, INTENT(IN)                                :: add_screening_info
      
      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_libheric_get_pmatrix_from_dbcsr_matrix'

      INTEGER                                            :: handle 
      INTEGER                                            :: iatom, la, lb, ikind, jkind, screening_info_size, &
                                                            global_nset(nkind), atom_block_size(nkind, nkind)
      INTEGER(int_8)                                     :: max_local_matrix_blocks, max_local_matrix_entries, &
                                                            num_matrix_entries, offset       
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_row_blocks
      INTEGER(int_8), DIMENSION(natom, 0:l_max, 0:l_max) :: row_offset_ll, num_row_ll_entries


      CALL timeset(routineN, handle)
      
      !! global_nset(ikind) is not equivalent to basis_parameter(ikind)%nset because there might be sets
      !! with different angular momentum for different shells.
        

      global_nset = 0
      DO la = 0, l_max
         DO ikind = 1, nkind
            global_nset(ikind) = global_nset(ikind) + basis_parameter_l(ikind, la)%nset
         END DO
      END DO

      DO ikind = 1, nkind
         DO jkind = 1, nkind
            !! screening info: (l_max + 1) ** 2 pmax_atom entries and global_nset(ikind) * global_nset(jkind) pmax_set entries per atom pair
            screening_info_size = MERGE((l_max + 1) ** 2 + global_nset(ikind) * global_nset(jkind), 0, add_screening_info)
            atom_block_size(jkind, ikind) = basis_parameter(ikind)%nsgf_total * basis_parameter(jkind)%nsgf_total & 
                                          + screening_info_size
         END DO
      END DO

      ALLOCATE(p_matrix%screening_info_size(nkind, nkind, 0:l_max, 0:l_max))
      ALLOCATE(p_matrix%atom_ll_block_size(nkind, nkind, 0:l_max, 0:l_max))

      DO la = 0, l_max
         DO lb = 0, l_max
            DO ikind = 1, nkind
                DO jkind = 1, nkind
                ASSOCIATE(bas_lai => basis_parameter_l(ikind, la),  bas_lbj => basis_parameter_l(jkind, lb))
                    p_matrix%atom_ll_block_size(jkind, ikind, lb, la) = bas_lai%nsgf_total * bas_lbj%nsgf_total
                    !! screening info: 1 pmax_atom entry + bas_lai%nset * bas_lbj%nset pmax_set entries per atom and angular momentum pair 
                    p_matrix%screening_info_size(jkind, ikind, lb, la) = &
                        MERGE(1 + bas_lai%nset * bas_lbj%nset, 0, add_screening_info)
                END ASSOCIATE 
                END DO
            END DO
         END DO
      END DO


      ALLOCATE(p_matrix%rows(natom))
      ALLOCATE(num_row_blocks(natom))
      CALL hfx_libheric_collect_dbcsr_matrix_properties(dbcsr_matrix, para_env, p_matrix, kind_of, natom, nkind, l_max, &
                                                        atom_block_size, max_local_matrix_blocks, max_local_matrix_entries, &
                                                        num_matrix_entries, num_row_ll_entries, num_row_blocks) 
      
      
      !!              actual matrix      + mapping for invalid blocks, i.e. guaranteed screening + zero block 
      p_matrix%size = num_matrix_entries + MERGE(1 + max_set ** 2, 0, add_screening_info) + MAXVAL(p_matrix%atom_ll_block_size)
      ALLOCATE(p_matrix%data(p_matrix%size))
      DO iatom = 1, natom
         ALLOCATE(p_matrix%rows(iatom)%key(num_row_blocks(iatom)))
         ALLOCATE(p_matrix%rows(iatom)%offset(0:num_row_blocks(iatom), 0:l_max, 0:l_max))
         !! Set all invalid indices to invalid block (meaningless if no screening info is added)
         p_matrix%rows(iatom)%offset(0, 0:l_max, 0:l_max) = num_matrix_entries + 1
         ALLOCATE(p_matrix%rows(iatom)%origin(num_row_blocks(iatom)))
      END DO

      !! Ensure that all invalid blocks are screened out by p-screening
      IF(add_screening_info) p_matrix%data(num_matrix_entries + 1:num_matrix_entries + 1 + max_set ** 2) = -1000.0_dp
      !! Add zero block
      p_matrix%zero_block => p_matrix%data(num_matrix_entries + MERGE(1 + max_set ** 2, 0, add_screening_info) + 1:p_matrix%size)
      p_matrix%zero_block(:) = 0.0_dp
      
      offset = 1
      DO la = 0, l_max
      DO lb = 0, l_max 
         DO iatom = 1, natom
            row_offset_ll(iatom, lb, la) = offset
            offset = offset + num_row_ll_entries(iatom, lb, la)
         END DO
      END DO
      END DO

      CALL hfx_libheric_gather_pmatrix(dbcsr_matrix, para_env, kind_of, &
                                       max_local_matrix_blocks, max_local_matrix_entries, &
                                       natom, nkind, l_max, basis_parameter, basis_parameter_l, & 
                                       atom_block_size, row_offset_ll, p_matrix, add_screening_info)
      CALL hfx_libheric_initialize_shellset_offset(max_shl, max_set, nkind, l_max, basis_parameter_l, &
                                                   p_matrix%shellset_offset, add_screening_info)

      CALL timestop(handle)

   END SUBROUTINE hfx_libheric_get_pmatrix_from_dbcsr_matrix
  
! **************************************************************************************************
!> \brief collect information about the size and dimension of the locally stored part of the 
!>        dbcsr matrix
!> \param dbcsr_matrix - local part of distributed dbcsr pmatrix
!> \param para_env - MPI environment
!> \param hfx_matrix - sparse matrix to be replicated on all ranks
!> \param kind_of - map of atom index to atom kind
!> \param natom - number of atoms in the simulated system
!> \param nkind - number of unique basis sets 
!> \param l_max - maximum angular momentum
!> \param atom_block_size - holds all possible number of values in a dbcsr_block + screening values
!>                          in that block, accessible via two atom kind indices 
!> \param max_local_matrix_blocks - maximum number of matrix blocks in any local part of 
!>                                  dbcsr matrix (output) 
!> \param max_local_matrix_entries - maximum number of matrix values + number of screening values
!>                                   in any local part of dbcsr matrix (output) 
!> \param num_matrix_entries - number of matrix values in global dbcsr matrix 
!>                             + number of screening values (output)
!> \param num_row_ll_entries - number of matrix values per ll-row in global dbcsr matrix
!>                             + number of screening values in that ll-row (output)
!> \param num_row_blocks - number of blocks per row in global dbcsr matrix (output)
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_collect_dbcsr_matrix_properties(dbcsr_matrix, para_env, hfx_matrix, kind_of, natom, nkind, l_max, &
                                                           atom_block_size, max_local_matrix_blocks, max_local_matrix_entries, &
                                                           num_matrix_entries, num_row_ll_entries, num_row_blocks)    
      TYPE(dbcsr_type), POINTER                    :: dbcsr_matrix
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(hfx_libheric_sparse_matrix), INTENT(IN) :: hfx_matrix
      INTEGER, DIMENSION(:), INTENT(IN)            :: kind_of
      INTEGER, INTENT(IN)                          :: natom, nkind, l_max
      INTEGER, DIMENSION(nkind, nkind), INTENT(IN) :: atom_block_size
      INTEGER(int_8), INTENT(OUT)                  :: max_local_matrix_blocks, &
                                                      max_local_matrix_entries, &
                                                      num_matrix_entries
      INTEGER(int_8), &
         DIMENSION(natom, 0:l_max, 0:l_max), &
         INTENT(OUT)                               :: num_row_ll_entries
      INTEGER, DIMENSION(natom), INTENT(OUT)       :: num_row_blocks

      INTEGER                                      :: iatom, jatom, ikind, jkind, &
                                                      la, lb, block_multiplicity
      TYPE(dbcsr_iterator_type)                    :: iter
      
      max_local_matrix_blocks = 0
      num_matrix_entries = 0
      num_row_ll_entries = 0
      num_row_blocks = 0

      
      CALL dbcsr_iterator_start(iter, dbcsr_matrix, shared=.FALSE.)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, iatom, jatom)
         ikind = kind_of(iatom)
         jkind = kind_of(jatom)

         !! Take care not to double count the diagonal
         block_multiplicity = MERGE(1, 2, iatom .EQ. jatom)
         max_local_matrix_blocks = max_local_matrix_blocks + block_multiplicity
         num_matrix_entries = num_matrix_entries + block_multiplicity * atom_block_size(jkind, ikind) 

         DO la = 0, l_max
         DO lb = 0, l_max
            num_row_ll_entries(iatom, lb, la) = num_row_ll_entries(iatom, lb, la) &
                                              + hfx_matrix%atom_ll_block_size(jkind, ikind, lb, la) & 
                                              + hfx_matrix%screening_info_size(jkind, ikind, lb, la)
            IF(iatom .NE. jatom) num_row_ll_entries(jatom, la, lb) = num_row_ll_entries(jatom, la, lb) &
                                                                   + hfx_matrix%atom_ll_block_size(jkind, ikind, lb, la) & 
                                                                   + hfx_matrix%screening_info_size(jkind, ikind, lb, la)
         END DO 
         END DO
         
         num_row_blocks(iatom) = num_row_blocks(iatom) + 1
         IF(iatom .NE. jatom) num_row_blocks(jatom) = num_row_blocks(jatom) + 1
      END DO
      CALL dbcsr_iterator_stop(iter)
      CPASSERT(SUM(num_row_ll_entries) .EQ. num_matrix_entries)
      
      max_local_matrix_entries = num_matrix_entries
      CALL para_env%sum(num_matrix_entries)
      CALL para_env%sum(num_row_ll_entries)
      CALL para_env%sum(num_row_blocks)
      CALL para_env%max(max_local_matrix_blocks)
      CALL para_env%max(max_local_matrix_entries)
      CPASSERT(SUM(num_row_ll_entries) .EQ. num_matrix_entries)

   END SUBROUTINE 
   
! **************************************************************************************************
!> \brief gather distributed dbcsr pmatrix and store in sparse format 
!> \param dbcsr_matrix - local part distributed dbcsr pmatrix
!> \param para_env - MPI environment
!> \param kind_of - map of atom index to atom kind
!> \param max_local_matrix_blocks - maximum number of matrix blocks in any local part of 
!>                                  dbcsr matrix 
!> \param max_local_matrix_entries - maximum number of matrix values + number of screening values
!>                                   in any local part of dbcsr matrix  
!> \param natom - number of atoms in the simulated system
!> \param nkind - number of unique basis sets 
!> \param l_max - maximum angular momentum
!> \param basis_parameter - all hfx basis sets
!> \param basis_parameter_l - all hfx_libheric basis sets
!> \param atom_block_size - holds all possible number of values in a dbcsr_block + screening values
!>                          in that block, accessible via two atom kind indices 
!> \param row_offset_ll - start offset of each ll-row into p_matrix
!> \param p_matrix - sparse pmatrix replicated on all MPI-ranks
!> \param add_screening_info - if .TRUE., add additional screening values in front of pmatrix blocks
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_gather_pmatrix(dbcsr_matrix, para_env, kind_of, &
                                          max_local_matrix_blocks, max_local_matrix_entries, &
                                          natom, nkind, l_max, basis_parameter, basis_parameter_l, &
                                          atom_block_size, row_offset_ll, p_matrix, add_screening_info)

      TYPE(dbcsr_type), POINTER                          :: dbcsr_matrix
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, DIMENSION(:)                              :: kind_of
      INTEGER(int_8), INTENT(IN)                         :: max_local_matrix_blocks, &
                                                            max_local_matrix_entries
      INTEGER, INTENT(IN)                                :: natom, nkind, l_max
      TYPE(hfx_basis_type), & 
          DIMENSION(nkind), INTENT(IN)                   :: basis_parameter
      TYPE(hfx_libheric_basis_type), & 
          DIMENSION(nkind, 0:l_max), INTENT(IN)          :: basis_parameter_l
      INTEGER, DIMENSION(nkind, nkind), INTENT(IN)       :: atom_block_size
      INTEGER(int_8), DIMENSION(natom, 0:l_max, 0:l_max), &
          INTENT(INOUT)                                  :: row_offset_ll
      TYPE(hfx_libheric_sparse_matrix), INTENT(INOUT)    :: p_matrix
      LOGICAL, INTENT(IN)                                :: add_screening_info

      INTEGER                                            :: iatom, jatom, ikind, jkind, la, lb, &
                                                            jatom_idx, num_ij_atoms, ij_atom_idx
      INTEGER                                            :: icpu, ncpu, mepos, dest, source, id_ctr, &
                                                            source_cpu, block_size, max_dbcsr_block_dim_size
      INTEGER(int_8)                                     :: i, j
      INTEGER, DIMENSION(:), POINTER                     :: id_sendbuffer, id_recvbuffer, id_swapbuffer
      REAL(dp), DIMENSION(:, :), POINTER                 :: sparse_block, sparse_block_T
      REAL(dp), DIMENSION(:), POINTER                    :: data_sendbuffer, data_recvbuffer, data_swapbuffer
      INTEGER, DIMENSION(natom)                          :: elem_in_block_row
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_request_type), DIMENSION(2)                :: data_req, id_req

      
      max_dbcsr_block_dim_size = MAXVAL(basis_parameter(:)%nsgf_total)

      ALLOCATE (data_sendbuffer(max_local_matrix_entries), source=-1.0_dp)
      ALLOCATE (data_recvbuffer(max_local_matrix_entries), source=-1.0_dp)
      ALLOCATE (id_sendbuffer(0:2 * max_local_matrix_blocks), source=-1)
      ALLOCATE (id_recvbuffer(0:2 * max_local_matrix_blocks), source=-1)
      ALLOCATE (sparse_block_T(max_dbcsr_block_dim_size, max_dbcsr_block_dim_size))
        
      CALL dbcsr_iterator_start(iter, dbcsr_matrix, shared=.FALSE.)
      id_ctr = 0
      i = 1
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block)
         sparse_block_T(1:SIZE(sparse_block, 2), 1:SIZE(sparse_block, 1)) = TRANSPOSE(sparse_block)
         ikind = kind_of(iatom)
         jkind = kind_of(jatom)
         
         ASSOCIATE(bas_i => basis_parameter(ikind), bas_j => basis_parameter(jkind))
         id_sendbuffer(2 * id_ctr + 1) = iatom
         id_sendbuffer(2 * id_ctr + 2) = jatom
         id_ctr = id_ctr + 1
         IF(iatom .NE. jatom) THEN 
            id_sendbuffer(2 * id_ctr + 1) = jatom
            id_sendbuffer(2 * id_ctr + 2) = iatom
            id_ctr = id_ctr + 1
            j = i + atom_block_size(jkind, ikind) 
         END IF
         DO la = 0, l_max
         DO lb = 0, l_max
            ASSOCIATE(bas_lai => basis_parameter_l(ikind, la), bas_laj => basis_parameter_l(jkind, la), &
                      bas_lbi => basis_parameter_l(ikind, lb), bas_lbj => basis_parameter_l(jkind, lb))
            block_size = p_matrix%atom_ll_block_size(jkind, ikind, lb, la) &
                       + p_matrix%screening_info_size(jkind, ikind, lb, la) 
            CALL hfx_libheric_extract_ll_blocks_from_dbcsr_block(la, lb, block_size, data_sendbuffer(i:i+block_size-1), &
                                                                 bas_lai, bas_lbj, bas_i, bas_j, sparse_block, &
                                                                 add_screening_info)
            i = i + block_size
            IF(iatom .NE. jatom) THEN
               block_size = p_matrix%atom_ll_block_size(ikind, jkind, lb, la) &
                          + p_matrix%screening_info_size(ikind, jkind, lb, la) 
               CALL hfx_libheric_extract_ll_blocks_from_dbcsr_block(la, lb, block_size, data_sendbuffer(j:j+block_size-1), &
                                                                    bas_laj, bas_lbi, bas_j, bas_i, sparse_block_T, &
                                                                    add_screening_info)
               j = j + block_size
            END IF
            END ASSOCIATE
         END DO
         END DO
         IF(iatom .NE. jatom) i = j
         END ASSOCIATE
      END DO
      CALL dbcsr_iterator_stop(iter)
      id_sendbuffer(0) = id_ctr 
      
      ! sync before/after a ring of isendrecv
      CALL para_env%sync()
      ncpu = para_env%num_pe
      mepos = para_env%mepos
      dest = MODULO(mepos + 1, ncpu)
      source = MODULO(mepos - 1, ncpu)
      elem_in_block_row = 0
      DO icpu = 0, ncpu - 1
         IF (icpu /= ncpu - 1) THEN
            CALL para_env%isendrecv(data_sendbuffer, dest, data_recvbuffer, source, &
                                    data_req(1), data_req(2), 13)
            CALL para_env%isendrecv(id_sendbuffer, dest, id_recvbuffer, source, &
                                    id_req(1), id_req(2), 14)
         END IF
         source_cpu = MODULO(mepos - icpu, ncpu) + 1

         num_ij_atoms = id_sendbuffer(0)
         i = 1
         DO ij_atom_idx = 0, num_ij_atoms - 1
             iatom = id_sendbuffer(2 * ij_atom_idx + 1)
             jatom = id_sendbuffer(2 * ij_atom_idx + 2)
             ikind = kind_of(iatom)
             jkind = kind_of(jatom)

             ASSOCIATE(hfx_row => p_matrix%rows(iatom))
             jatom_idx = elem_in_block_row(iatom) + 1
             hfx_row%key(jatom_idx) = jatom
             hfx_row%origin(jatom_idx) = source_cpu

             DO la = 0, l_max
             DO lb = 0, l_max
                block_size = p_matrix%atom_ll_block_size(jkind, ikind, lb, la) &
                           + p_matrix%screening_info_size(jkind, ikind, lb, la) 
                j = row_offset_ll(iatom, lb, la)
                hfx_row%offset(jatom_idx, lb, la) = j
                p_matrix%data(j:j + block_size - 1) = data_sendbuffer(i:i + block_size - 1)
                row_offset_ll(iatom, lb, la) = j + block_size 
                i = i + block_size
             END DO
             END DO
             END ASSOCIATE
             elem_in_block_row(iatom) = elem_in_block_row(iatom) + 1
         END DO

         IF (icpu /= ncpu - 1) THEN
            CALL mp_waitall(data_req)
            CALL mp_waitall(id_req)
         END IF
         
         data_swapbuffer => data_sendbuffer
         data_sendbuffer => data_recvbuffer
         data_recvbuffer => data_swapbuffer

         data_recvbuffer = -1.0_dp
         
         id_swapbuffer => id_sendbuffer
         id_sendbuffer => id_recvbuffer
         id_recvbuffer => id_swapbuffer

         id_recvbuffer = -1
      END DO
      DEALLOCATE (id_sendbuffer, id_recvbuffer)
      DEALLOCATE (data_sendbuffer, data_recvbuffer)
      DEALLOCATE (sparse_block_T)
      ! sync before/after a ring of isendrecv
      CALL para_env%sync()
   END SUBROUTINE hfx_libheric_gather_pmatrix
   
! **************************************************************************************************
!> \brief extract all ll-blocks belonging to angular momemntum pair (la, lb) from dbcsr block 
!> \param la,lb - requested angular momentum pair
!> \param block_size - number of extracted values
!> \param lin_block - linearized block holding the (la, lb) blocks (output)
!> \param bas_la, bas_lb - hfx_libheric basis set belonging to the to the atoms in the dbcsr block
!> \param bas_a, bas_b - hfx basis set belonging to the to the atoms in the dbcsr block
!> \param dbcsr_block - ...
!> \param add_screening_info - if .TRUE., add additional screening values in front of blocks
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_extract_ll_blocks_from_dbcsr_block(la, lb, block_size, lin_block, bas_la, bas_lb, &
                                                              bas_a, bas_b, dbcsr_block, add_screening_info)
      INTEGER                                        :: la, lb, block_size 
      REAL(dp), DIMENSION(block_size), INTENT(INOUT) :: lin_block
      TYPE(hfx_libheric_basis_type), INTENT(IN)      :: bas_la, bas_lb
      TYPE(hfx_basis_type), INTENT(IN)               :: bas_a, bas_b
      REAL(dp), DIMENSION(:, :), POINTER             :: dbcsr_block
      LOGICAL, INTENT(IN)                            :: add_screening_info

      INTEGER                                        :: iset, jset, ishl, jshl, pa, pb, pa1, pb1, i_set, i_val
      REAL(dp)                                       :: dbcsr_entry, max_atom_entry, max_set_entry

      !! Layout of a linearized atom block (assuming screening info is collected)
      !! lin_block(1) -> log10 of maxval in block
      !! lin_block(2:1 + nseta * nsetb) -> log10 of maxval in set-set block
      !! lin_block(1 + nseta * nsetb + 1:end) -> linearized values of dbscr_block
      i_val = MERGE(1 + 1 + bas_la%nset * bas_lb%nset, 1, add_screening_info)
      i_set = 2
      max_atom_entry = -HUGE(0.0_dp)
      DO iset = 1, bas_la%nset
         DO jset = 1, bas_lb%nset
            max_set_entry = -HUGE(0.0_dp)
            DO ishl = 1, bas_la%nshell(iset)
               DO jshl = 1, bas_lb%nshell(jset)
                  ! Copy shellset into lin_block
                  ! Look up position of shellset in global basis set
                  ! If a set has multiple shells, then the values of a shellset block are not stored contiguously
                  ! Note that the order of pa and pb are "column major" here, even though every surrounding loop is row major
                  pa = bas_a%first_sgf(bas_la%global_shl(ishl, iset), bas_la%global_set(iset))
                  pb = bas_b%first_sgf(bas_lb%global_shl(jshl, jset), bas_lb%global_set(jset))
                  DO pb1 = pb, pb + nso(lb) - 1
                    DO pa1 = pa, pa + nso(la) - 1
                        dbcsr_entry = dbcsr_block(pa1, pb1)
                        max_atom_entry = MAX(max_atom_entry, ABS(dbcsr_entry))
                        max_set_entry = MAX(max_set_entry, ABS(dbcsr_entry))
                        lin_block(i_val) = dbcsr_block(pa1, pb1) 
                        i_val = i_val + 1
                     END DO
                  END DO
               END DO
            END DO
            IF(add_screening_info) THEN
               lin_block(i_set) = MERGE(LOG10(max_set_entry), -1000.0_dp, max_set_entry .GT. 0.0_dp)
               i_set = i_set + 1
            END IF
         END DO
      END DO
      IF(add_screening_info) lin_block(1) = MERGE(LOG10(max_atom_entry), -1000.0_dp, max_atom_entry .GT. 0.0_dp)
   END SUBROUTINE hfx_libheric_extract_ll_blocks_from_dbcsr_block

! **************************************************************************************************
!> \brief Compare new pmatrix format to "old" format; check that all matrix entries and all screening
!>        values match
!> \param x_data - data related HFX
!> \param para_env - MPI environment
!> \param rho - local part of dbcsr p_matrix 
!> \param ks - local part of dbcsr k_matrix
!> \param kind_of - map of atom index to atom kind
!> \param natom - number of atoms in the simulated system
!> \param nkind - number of unique basis sets 
!> \param l_max - maximum angular momentum
!> \param max_nsgf_total - maximum dbcsr block side-length
!> \param basis_parameter_l - all hfx_libheric basis sets
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_pmatrix_verify(x_data, para_env, rho, ks, kind_of, natom, nkind, l_max,  &
                                          max_nsgf_total, basis_parameter_l, p_matrix)
      TYPE(hfx_type), POINTER, INTENT(IN)          :: x_data
      TYPE(mp_para_env_type), POINTER              :: para_env
      TYPE(dbcsr_p_type), INTENT(INOUT)            :: rho, ks
      INTEGER, DIMENSION(:), INTENT(IN)            :: kind_of
      INTEGER, INTENT(IN)                          :: natom, nkind, l_max, max_nsgf_total
      TYPE(hfx_libheric_basis_type), &
         DIMENSION(nkind, 0:l_max), INTENT(IN)     :: basis_parameter_l 
      TYPE(hfx_libheric_sparse_matrix), INTENT(IN) :: p_matrix

      LOGICAL                             :: do_p_screening
      INTEGER                             :: la, iatom, ikind, iset, ishl, isetg, ishlg, &
                                             lb, jatom, jkind, jset, jshl, jsetg, jshlg, &
                                             prev_iatom, iatom_ref, jatom_ref, &
                                             ikind_ref, jkind_ref, la_ref, lb_ref, &
                                             pa, pb, i, pmax_atom_off, block_multiplicity, &
                                             num_pmax_exact_matches
      INTEGER(int_8)                      :: den_off, ref_off, ref_idx, kind_kind_idx
      REAL(dp), DIMENSION(:, :), POINTER  :: rho_ref
      REAL(dp), DIMENSION(:), ALLOCATABLE :: ref_subblock
      INTEGER, DIMENSION(natom)           :: p_cols
      LOGICAL                             :: assert_flag, ref_not_transposed
      INTEGER, DIMENSION(:), ALLOCATABLE  :: ref_checked
      REAL(dp)                            :: pmax_set_den, pmax_set_ref
      CHARACTER(len=2048)                 :: kind_str, ll_str, set_str, shl_str
      CHARACTER(len=:), ALLOCATABLE       :: error_str, tmp_str

      ALLOCATE(CHARACTER(2 * 32 * max_nsgf_total ** 2 + 2048) :: error_str) 
      ALLOCATE(CHARACTER(2 * 32 * max_nsgf_total + 2048) :: tmp_str)

      
      do_p_screening = x_data%screening_parameter%do_initial_p_screening
      ALLOCATE(ref_subblock(nso(l_max) ** 2))
      CALL get_atomic_block_maps(ks%matrix, x_data%basis_parameter, kind_of, &
                                 x_data%is_assoc_atomic_block, &
                                 x_data%number_of_p_entries, &
                                 para_env, &
                                 x_data%atomic_block_offset, &
                                 x_data%set_offset, &
                                 x_data%block_offset, &
                                 x_data%map_atoms_to_cpus, &
                                 nkind)


      ALLOCATE(rho_ref(x_data%block_offset(para_env%num_pe + 1) - 1, 1))
      CALL get_full_density(para_env, rho_ref(:, 1), rho%matrix, x_data%number_of_p_entries, &
                            x_data%block_offset, &
                            kind_of, x_data%basis_parameter, get_max_vals_spin=.FALSE., antisymmetric=.FALSE.)
      IF(do_p_screening) &
         CALL update_pmax_mat(x_data%initial_p, &
                              x_data%map_atom_to_kind_atom, &
                              x_data%set_offset, &
                              x_data%atomic_block_offset, &
                              x_data%pmax_atom, &
                              rho_ref, NULL(), &
                              natom, kind_of, x_data%basis_parameter, &
                              nkind, x_data%is_assoc_atomic_block)
      ALLOCATE(ref_checked(SIZE(rho_ref)))
      ref_checked = 0
      pmax_atom_off = MERGE(1, 0, do_p_screening)

      prev_iatom = -1
      DO iatom = 1, natom
         ASSOCIATE(p_row => p_matrix%rows(iatom))
         CALL hfx_libheric_get_column_index_map(p_matrix, iatom, prev_iatom, p_cols)
         DO jatom = 1, natom
            ref_not_transposed = iatom .LE. jatom
            ikind = kind_of(iatom)
            jkind = kind_of(jatom)

            iatom_ref = MERGE(iatom, jatom, ref_not_transposed)
            jatom_ref = MERGE(jatom, iatom, ref_not_transposed)
            ikind_ref = MERGE(ikind, jkind, ref_not_transposed)
            jkind_ref = MERGE(jkind, ikind, ref_not_transposed)
             
            WRITE(kind_str, "(3(A,I5),4(A,I2))") "rank = ", para_env%mepos, &
                                                 ", iatom = ", iatom, ", jatom = ", jatom, &
                                                 ", ikind = ", ikind, ", jkind = ", jkind, &
                                                 ", ikind_ref = ", ikind_ref, ", jkind_ref = ", jkind_ref
            num_pmax_exact_matches = 0
            DO la = 0, l_max
            DO lb = 0, l_max
               WRITE(ll_str, "(A,3X,2(A,I1))") TRIM(kind_str) // NEW_LINE('A'), "la = ", la, ", lb = ", lb
               ASSOCIATE(bas_a => x_data%basis_parameter(ikind_ref), bas_b => x_data%basis_parameter(jkind_ref), &
                         bas_la => basis_parameter_l(ikind, la), bas_lb => basis_parameter_l(jkind, lb), &
                         ref_block_ij_off => x_data%atomic_block_offset(iatom, jatom), &
                         den_block_ij_off => p_row%offset(p_cols(jatom), lb, la), &
                         set_off => p_matrix%shellset_offset(:, :, :, :, jkind, ikind, lb, la), &
                         block_size => nso(la) * nso(lb))
                  IF(do_p_screening) THEN
                     assert_flag = x_data%pmax_atom(jatom, iatom) .GE. p_matrix%data(den_block_ij_off)
                     IF(.NOT. assert_flag) THEN
                         WRITE(error_str, "(A,3X,A,2(A,ES28.20))") TRIM(ll_str) // NEW_LINE('A'), "pmax_atom mismatch: ", &
                                                                    "ref_pmax = ", x_data%pmax_atom(jatom, iatom), &
                                                                    ", den_pmax = ", p_matrix%data(den_block_ij_off)
                         CALL cp_abort(__LOCATION__, TRIM(error_str))
                     END IF
                     IF(x_data%pmax_atom(jatom, iatom) .EQ. p_matrix%data(den_block_ij_off)) &
                         num_pmax_exact_matches = num_pmax_exact_matches + 1
                  END IF
                  DO iset = 1, bas_la%nset
                  DO jset = 1, bas_lb%nset
                     isetg = MERGE(bas_la%global_set(iset), bas_lb%global_set(jset), ref_not_transposed)
                     jsetg = MERGE(bas_lb%global_set(jset), bas_la%global_set(iset), ref_not_transposed)
                     la_ref = MERGE(la, lb, ref_not_transposed) 
                     lb_ref = MERGE(lb, la, ref_not_transposed) 

                     WRITE(set_str, "(A,6X,4(A,I2),A)") TRIM(ll_str) // NEW_LINE('A'), "iset = ", iset, "(", isetg, &
                                                                                    "), jset = ", jset, "(", jsetg, ")"
                     assert_flag = bas_a%lmin(isetg) .LE. la_ref .AND. la_ref .LE. bas_a%lmax(isetg)
                     CPASSERT(assert_flag) 
                     assert_flag = bas_b%lmin(jsetg) .LE. lb_ref .AND. lb_ref .LE. bas_b%lmax(jsetg)
                     CPASSERT(assert_flag) 

                     ref_off = ref_block_ij_off + x_data%set_offset(jsetg, isetg, jkind_ref, ikind_ref) - 1
                     den_off = den_block_ij_off + set_off(1, 1, jset, iset)

                     IF(do_p_screening) THEN
                        pmax_set_den = p_matrix%data(den_block_ij_off + (iset - 1) * bas_lb%nset + jset)
                        kind_kind_idx = get_1D_idx(ikind_ref, jkind_ref, INT(nkind, kind=int_8))
                        IF (ikind_ref .GT. jkind_ref) THEN
                           pmax_set_ref = x_data%initial_p(kind_kind_idx)%p_kind(jsetg, isetg, &
                                                                                 x_data%map_atom_to_kind_atom(jatom_ref), &
                                                                                 x_data%map_atom_to_kind_atom(iatom_ref))
                        ELSE
                           pmax_set_ref = x_data%initial_p(kind_kind_idx)%p_kind(isetg, jsetg, &
                                                                                 x_data%map_atom_to_kind_atom(iatom_ref), &
                                                                                 x_data%map_atom_to_kind_atom(jatom_ref))
                        END IF
                        !! Check for sets with shell on different angular momentum
                        IF(bas_a%lmin(isetg) .EQ. bas_a%lmax(isetg) .AND. &
                           bas_b%lmin(jsetg) .EQ. bas_b%lmax(jsetg)) THEN
                           assert_flag = pmax_set_ref .EQ. pmax_set_den
                        ELSE
                           assert_flag = pmax_set_ref .GE. pmax_set_den
                        END IF
                        IF(.NOT. assert_flag) THEN
                            WRITE(error_str, "(A,6X,3(A,ES28.20))") TRIM(set_str) // NEW_LINE('A'), "pmax_set mismatch: ", &
                                                                      "ref_set_pmax = ", pmax_set_ref, &
                                                                      ", den_set_pmax = ", pmax_set_den, &
                                                                      ", diff = ", ABS(pmax_set_ref - pmax_set_den)
                            CALL cp_abort(__LOCATION__, TRIM(error_str))
                        END IF
                     END IF 

                     DO ishl = 1, bas_la%nshell(iset)
                     DO jshl = 1, bas_lb%nshell(jset)
                        ishlg = MERGE(bas_la%global_shl(ishl, iset), bas_lb%global_shl(jshl, jset), ref_not_transposed)
                        jshlg = MERGE(bas_lb%global_shl(jshl, jset), bas_la%global_shl(ishl, iset), ref_not_transposed)


                        WRITE(shl_str, "(A,9X,4(A,I1),A,3(A,I10))") TRIM(set_str) // NEW_LINE('A'), & 
                                                                     "ishl = ", ishl, "(", ishlg, &
                                                                     "), jshl = ", jshl, "(", jshlg, ")", &
                                                                     " -- den_block_ij_of = ", den_block_ij_off, &
                                                                     ", set_off = ", set_off(jshl, ishl, jset, iset), &
                                                                     ", den_off = ", den_off

                        ASSOCIATE(sgfa => bas_a%first_sgf(ishlg, isetg) - bas_a%first_sgf(1, isetg), &
                                  sgfb => bas_b%first_sgf(jshlg, jsetg) - bas_b%first_sgf(1, jsetg), &
                                  nsgfa => bas_a%nsgf(isetg), nsgfb => bas_b%nsgf(jsetg), &
                                  den_block => p_matrix%data(den_off:den_off + block_size - 1))
                        
                        i = 1
                        ref_subblock = -123.456
                        !! A bit confusingly, the individual set-blocks in full_density are stored in row-major
                        !! While the order of the sets is column major (see L211-227)
                        IF (.NOT. ref_not_transposed) THEN
                           DO pa = 0, nso(la_ref) - 1
                              DO pb = 0, nso(lb_ref) - 1
                                 ref_idx = ref_off + (sgfa + pa) * nsgfb + sgfb + pb
                                 ref_subblock(i) = rho_ref(ref_idx, 1) 
                                 ref_checked(ref_idx) = ref_checked(ref_idx) + 1
                                 i = i + 1
                              END DO
                           END DO
                        ELSE
                           DO pb = 0, nso(lb_ref) - 1
                              DO pa = 0, nso(la_ref) - 1
                                 ref_idx = ref_off + (sgfa + pa) * nsgfb + sgfb + pb
                                 ref_subblock(i) = rho_ref(ref_idx, 1) 
                                 ref_checked(ref_idx) = ref_checked(ref_idx) + 1
                                 i = i + 1
                              END DO
                           END DO
                        END IF
                         
                        assert_flag = ALL(ref_subblock(1:block_size) .EQ. den_block(1:block_size))
                        IF(.NOT. assert_flag) THEN
                           WRITE(error_str, "(A,9X,3(A,ES28.20))") TRIM(shl_str) // NEW_LINE('A'), "density block mismatch: " 
                           WRITE(tmp_str, "(9X,A)") "Density:"
                           error_str = TRIM(error_str) // NEW_LINE('A') // TRIM(tmp_str) // NEW_LINE('A')
                           DO pa = 0, nso(la) - 1
                              WRITE(tmp_str, "(9X,*(ES28.20,2X))") den_block(1 + pa * nso(lb):(pa + 1) * nso(lb))
                              error_str = TRIM(error_str) // TRIM(tmp_str) // NEW_LINE('A')
                           END DO
                           
                           WRITE(tmp_str, "(9X,A)") "Reference Density:"
                           error_str = TRIM(error_str) // NEW_LINE('A') // TRIM(tmp_str) // NEW_LINE('A')
                           DO pa = 0, nso(la) - 1
                              WRITE(tmp_str, "(9X,*(ES28.20,2X))") ref_subblock(1 + pa * nso(lb):(pa + 1) * nso(lb))
                              error_str = TRIM(error_str) // TRIM(tmp_str) // NEW_LINE('A')
                           END DO
                           CALL cp_abort(__LOCATION__, TRIM(error_str))
                        END IF

                        END ASSOCIATE
                        CPASSERT(bas_a%nl(ishlg, isetg) .EQ. la_ref) 
                        CPASSERT(bas_b%nl(jshlg, jsetg) .EQ. lb_ref)
                        den_off = den_off + block_size
                     END DO
                     END DO
                  END DO
                  END DO
               END ASSOCIATE
            END DO
            END DO
            IF(do_p_screening) THEN
               assert_flag = num_pmax_exact_matches .GT. 0
               IF(.NOT. assert_flag) THEN
                   WRITE(error_str, "(A,A)") &
                       TRIM(kind_str) // NEW_LINE('A'), "pmax_atom mismatch: None of the den_pmax values are equal to ref_pmax"
                   CALL cp_abort(__LOCATION__, TRIM(error_str))
               END IF
            END IF
         END DO
         END ASSOCIATE
      END DO
      DO iatom = 1, natom
      DO jatom = iatom, natom
         ikind = kind_of(iatom)
         jkind = kind_of(jatom)
         WRITE(kind_str, "(3(A,I5),2(A,I2))") "rank = ", para_env%mepos, &
                                              ", iatom = ", iatom, ", jatom = ", jatom, &
                                              ", ikind = ", ikind, ", jkind = ", jkind
         ASSOCIATE(bas_a => x_data%basis_parameter(ikind_ref), bas_b => x_data%basis_parameter(jkind_ref), &
                   ref_off => x_data%atomic_block_offset(iatom, jatom))
         block_multiplicity = MERGE(1, 2, iatom .EQ. jatom)
         assert_flag = ALL(ref_checked(ref_off:ref_off + bas_a%nsgf_total * bas_b%nsgf_total - 1) .EQ. block_multiplicity)
         IF(.NOT. assert_flag) THEN
             WRITE(error_str, "(A,A)") &
                 TRIM(kind_str) // NEW_LINE('A'), "Not all entries in ref block have been checked"
             CALL cp_abort(__LOCATION__, TRIM(error_str))
         END IF
         END ASSOCIATE
      END DO 
      END DO
      DEALLOCATE(rho_ref)
      CALL para_env%sync()
      
   END SUBROUTINE hfx_libheric_pmatrix_verify
   
! **************************************************************************************************
!> \brief initialize matrix for hfx contribution to kohn-sham matrix; matrix is initialized with 0s
!>        and is identical on each MPI rank
!> \param p_matrix - p_matrix replicated on all MPI ranks
!> \param kind_of - map of atom index to atom kind
!> \param natom - number of atoms in the simulated system
!> \param nkind - number of unique basis sets 
!> \param l_max - maximum angular momentum
!> \param max_set - maximum number of shellsets in any basis set
!> \param max_shl - maximum number of shells in any shellset
!> \param basis_parameter_l - all hfx_libheric basis sets
!> \param k_matrix - k_matrix replicated on all MPI-ranks
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_initialize_kmatrix(p_matrix, kind_of, natom, nkind, l_max, max_set, max_shl, &
                                              basis_parameter_l, k_matrix)
      TYPE(hfx_libheric_sparse_matrix), INTENT(IN)       :: p_matrix
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of
      INTEGER, INTENT(IN)                                :: natom, nkind, l_max, max_set, max_shl
      TYPE(hfx_libheric_basis_type), & 
          DIMENSION(nkind, 0:l_max), INTENT(IN)          :: basis_parameter_l
      TYPE(hfx_libheric_sparse_matrix), INTENT(OUT)      :: k_matrix
      
      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_libheric_initialize_kmatrix'

      INTEGER                                            :: handle 
      INTEGER                                            :: iatom, jatom, la, lb
      INTEGER(int_8)                                     :: offset


      CALL timeset(routineN, handle)

      ALLOCATE(k_matrix%screening_info_size(nkind, nkind, 0:l_max, 0:l_max))
      ALLOCATE(k_matrix%atom_ll_block_size(nkind, nkind, 0:l_max, 0:l_max))

      k_matrix%screening_info_size = 0
      k_matrix%atom_ll_block_size(:, :, :, :) = p_matrix%atom_ll_block_size

      ALLOCATE(k_matrix%rows(natom))

      offset = 1
      DO iatom = 1, natom
         ALLOCATE(k_matrix%rows(iatom)%key(natom))
         ALLOCATE(k_matrix%rows(iatom)%offset(0:natom,0:l_max,0:l_max))
         ALLOCATE(k_matrix%rows(iatom)%origin(natom))
         DO la = 0, l_max
            DO lb = 0, l_max
               DO jatom = 1, natom
                  ASSOCIATE(ikind => kind_of(iatom), jkind => kind_of(jatom))
                  k_matrix%rows(iatom)%key(jatom) = jatom
                  k_matrix%rows(iatom)%offset(jatom, lb, la) = offset
                  k_matrix%rows(iatom)%origin = -1
                  offset = offset + k_matrix%atom_ll_block_size(jkind, ikind, lb, la)
                  END ASSOCIATE
               END DO
            END DO
         END DO
      END DO
      k_matrix%size = offset + MAXVAL(k_matrix%atom_ll_block_size) - 1
      ALLOCATE(k_matrix%data(k_matrix%size), source=0.0_dp)
      DO iatom = 1, natom
         !! Map invalid indices to zero block 
         k_matrix%rows(iatom)%offset(0, 0:l_max, 0:l_max) = offset
      END DO
      k_matrix%zero_block => k_matrix%data(offset:k_matrix%size)

      CALL hfx_libheric_initialize_shellset_offset(max_shl, max_set, nkind, l_max, basis_parameter_l, &
                                                   k_matrix%shellset_offset, .FALSE.)

      CALL timestop(handle)

   END SUBROUTINE hfx_libheric_initialize_kmatrix

! **************************************************************************************************
!> \brief gather replicated k_matrix and add blocks belonging to the caller's MPI rank back into
!>        local part of the dbcsr kohn-sham matrix
!> \param k_matrix - k_matrix replicated on all MPI ranks
!> \param para_env - MPI environment
!> \param kind_of - map of atom index to atom kind
!> \param natom - number of atoms in the simulated system
!> \param nkind - number of unique basis sets 
!> \param l_max - maximum angular momentum
!> \param basis_parameter - all hfx basis sets
!> \param basis_parameter_l - all hfx_libheric basis sets
!> \param eps_filter - could be used to filter the resulting dbcsr matrix, so that it is sparse, but
!                      unused right now
!> \param dbcsr_matrix - local part of distributed dbcsr kmatrix
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_reconstruct_dbcsr_kmatrix(k_matrix, para_env, kind_of, natom, nkind, l_max, &
                                                     basis_parameter, basis_parameter_l, eps_filter, dbcsr_matrix)
      TYPE(hfx_libheric_sparse_matrix), INTENT(IN)       :: k_matrix
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, DIMENSION(:)                              :: kind_of
      INTEGER, INTENT(IN)                                :: natom, nkind, l_max 
      TYPE(hfx_basis_type), & 
          DIMENSION(nkind), INTENT(IN)                   :: basis_parameter
      TYPE(hfx_libheric_basis_type), & 
          DIMENSION(nkind, 0:l_max), INTENT(IN)          :: basis_parameter_l
      REAL(dp), INTENT(IN)                               :: eps_filter
      TYPE(dbcsr_type), POINTER, INTENT(INOUT)           :: dbcsr_matrix
      
      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_libheric_reconstruct_dbcsr_kmatrix'
    
      INTEGER                                            :: max_dbcsr_block_size, node_id, &
                                                            la, iatom, ikind, &
                                                            lb, jatom, jkind, &
                                                            jkey_idx, pa, pb, handle
      LOGICAL                                            :: need_transpose, is_diag_block
      REAL(dp), DIMENSION(:), ALLOCATABLE, TARGET        :: lin_dbcsr_blocks
      REAL(dp), DIMENSION(:), POINTER                    :: lin_dbcsr_block, lin_dbcsr_block_T 
      REAL(dp), DIMENSION(:, :), POINTER                 :: dbcsr_block, dbcsr_block_T

    
      MARK_USED(eps_filter)
      CALL timeset(routineN, handle)
      !! For now assumption is that k_matrix is not sparse and replicated on each rank with the same offsets
      CALL para_env%sum(k_matrix%data)

      max_dbcsr_block_size = MAXVAL([((basis_parameter(ikind)%nsgf_total * basis_parameter(jkind)%nsgf_total, &
                                       jkind=1,nkind), ikind=1,nkind)])
      ALLOCATE(lin_dbcsr_blocks(2 * max_dbcsr_block_size))
      lin_dbcsr_block => lin_dbcsr_blocks(1:max_dbcsr_block_size)
      lin_dbcsr_block_T => lin_dbcsr_blocks(max_dbcsr_block_size + 1:2 * max_dbcsr_block_size)

      DO iatom = 1, natom
         ikind = kind_of(iatom)
         ASSOCIATE(hfx_row => k_matrix%rows(iatom))
         DO jkey_idx = 1, SIZE(hfx_row%key)
            jatom = hfx_row%key(jkey_idx)
            jkind = kind_of(jatom)
            CALL dbcsr_get_stored_coordinates(dbcsr_matrix, iatom, jatom, node_id)
            is_diag_block = iatom .EQ. jatom
            ASSOCIATE(bas_i => basis_parameter(ikind), bas_j => basis_parameter(jkind))
            IF(node_id .EQ. para_env%mepos) THEN
               lin_dbcsr_block = 0.0_dp
               dbcsr_block(1:bas_i%nsgf_total,1:bas_j%nsgf_total) => lin_dbcsr_block(1:bas_i%nsgf_total * bas_j%nsgf_total)
               dbcsr_block_T(1:bas_j%nsgf_total,1:bas_i%nsgf_total) => lin_dbcsr_block_T(1:bas_i%nsgf_total * bas_j%nsgf_total)
               DO la = 0, l_max
               DO lb = 0, l_max
                  ASSOCIATE(offset => hfx_row%offset(jkey_idx, lb, la), & 
                            block_size => k_matrix%atom_ll_block_size(jkind, ikind, lb, la), &
                            bas_lai => basis_parameter_l(ikind, la), bas_lbj => basis_parameter_l(jkind, lb))
                  CALL hfx_libheric_reassemble_dbcsr_block_from_ll_block(la, lb, k_matrix%data(offset:offset + block_size - 1), &
                                                                         bas_lai, bas_lbj, bas_i, bas_j, dbcsr_block, is_diag_block)
                  END ASSOCIATE
               END DO
               END DO
               need_transpose = iatom .GT. jatom
               IF(.NOT. need_transpose) THEN
                   IF(is_diag_block) THEN
                      DO pb = 1, bas_j%nsgf_total
                          DO pa = pb + 1, bas_i%nsgf_total
                             dbcsr_block(pa, pb) = dbcsr_block(pa, pb) + dbcsr_block(pb, pa)
                             dbcsr_block(pb, pa) = dbcsr_block(pa, pb) 
                          END DO
                      END DO
                   END IF
                   CALL dbcsr_put_block(dbcsr_matrix, iatom, jatom, dbcsr_block, summation=.TRUE.)
               ELSE
                   dbcsr_block_T(:, :) = TRANSPOSE(dbcsr_block)
                   CALL dbcsr_put_block(dbcsr_matrix, jatom, iatom, dbcsr_block_T, summation=.TRUE.)
               END IF
            END IF
            END ASSOCIATE
         END DO
         END ASSOCIATE
      END DO
      CALL dbcsr_finalize(dbcsr_matrix)
      !! dbcsr_filter reshuffles the matrix structure
      !CALL dbcsr_filter(dbcsr_matrix, eps_filter)
      CALL timestop(handle)
   END SUBROUTINE hfx_libheric_reconstruct_dbcsr_kmatrix
   
! **************************************************************************************************
!> \brief insert (la, lb)-blocks into kohn-sham dbcsr matrix block 
!> \param la,lb - requested angular momentum pair
!> \param lin_block - linearized block holding the (la, lb) blocks
!> \param bas_la, bas_lb - hfx_libheric basis set belonging to the to the atoms in the dbcsr block
!> \param bas_a, bas_b - hfx basis set belonging to the to the atoms in the dbcsr block
!> \param dbcsr_block - ...
!> \param is_diag_block - special treatment for diagonal block for a symmetric matrix
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_reassemble_dbcsr_block_from_ll_block(la, lb, lin_ll_block, bas_la, bas_lb, &
                                                                bas_a, bas_b, dbcsr_block, is_diag_block)
      INTEGER, INTENT(IN)                            :: la, lb 
      REAL(dp), DIMENSION(:), INTENT(IN)             :: lin_ll_block
      TYPE(hfx_libheric_basis_type), INTENT(IN)      :: bas_la, bas_lb
      TYPE(hfx_basis_type), INTENT(IN)               :: bas_a, bas_b
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)       :: dbcsr_block
      LOGICAL, INTENT(IN)                            :: is_diag_block

      INTEGER                                        :: iset, jset, ishl, jshl, pa, pb, pa1, pb1, i
      REAL(dp)                                       :: diag_fac

      i = 1
      DO iset = 1, bas_la%nset
         DO jset = 1, bas_lb%nset
            DO ishl = 1, bas_la%nshell(iset)
               DO jshl = 1, bas_lb%nshell(jset)
                  ! Copy shellset into dbcsr_block
                  ! Look up position of shellset in global basis set
                  ! If a set has multiple shells, then the values of a shellset block are not stored contiguously
                  pa = bas_a%first_sgf(bas_la%global_shl(ishl, iset), bas_la%global_set(iset))
                  pb = bas_b%first_sgf(bas_lb%global_shl(jshl, jset), bas_lb%global_set(jset))
                  DO pb1 = pb, pb + nso(lb) - 1
                     DO pa1 = pa, pa + nso(la) - 1
                        diag_fac = MERGE(1.0_dp, 0.5_dp, is_diag_block .AND. pa1 .EQ. pb1)
                        dbcsr_block(pa1, pb1) = lin_ll_block(i) * diag_fac
                        i = i + 1
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE hfx_libheric_reassemble_dbcsr_block_from_ll_block

! **************************************************************************************************
!> \brief deallocate hfx_libheric_sparse_matrix 
!> \param hfx_matrix - ...
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_free_sparse_matrix(hfx_matrix)
      TYPE(hfx_libheric_sparse_matrix) :: hfx_matrix

      INTEGER :: i

      IF(ALLOCATED(hfx_matrix%rows)) THEN
         DO i = 1, SIZE(hfx_matrix%rows)
            IF(ALLOCATED(hfx_matrix%rows(i)%key))DEALLOCATE(hfx_matrix%rows(i)%key)
            IF(ALLOCATED(hfx_matrix%rows(i)%offset))DEALLOCATE(hfx_matrix%rows(i)%offset)
            IF(ALLOCATED(hfx_matrix%rows(i)%origin))DEALLOCATE(hfx_matrix%rows(i)%origin)
         END DO
         DEALLOCATE(hfx_matrix%rows)
      END IF
      IF(ASSOCIATED(hfx_matrix%zero_block)) hfx_matrix%zero_block => NULL()
      IF(ASSOCIATED(hfx_matrix%data)) DEALLOCATE(hfx_matrix%data)
      IF(ASSOCIATED(hfx_matrix%shellset_offset)) DEALLOCATE(hfx_matrix%shellset_offset)
      IF(ALLOCATED(hfx_matrix%screening_info_size)) DEALLOCATE(hfx_matrix%screening_info_size)
      IF(ALLOCATED(hfx_matrix%atom_ll_block_size)) DEALLOCATE(hfx_matrix%atom_ll_block_size)

   END SUBROUTINE hfx_libheric_free_sparse_matrix


! **************************************************************************************************
!> \brief update column index map into hfx_libheric sparse matrix 
!> \param hfx_matrix - hfx_libheric sparse matrix for which to create the column index map
!> \param row - row in hfx_matrix for which to create the column index map
!> \param prev_row - row that was previously mapped, put -1 if no row was previously mapped
!> \param col_map - colum index map for requested row (output)
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_get_column_index_map(hfx_matrix, row, prev_row, col_map)
      TYPE(hfx_libheric_sparse_matrix), INTENT(IN)  :: hfx_matrix
      INTEGER, INTENT(IN)                           :: row, prev_row
      INTEGER, DIMENSION(:), INTENT(INOUT)          :: col_map

      INTEGER :: col

      IF(row .EQ. prev_row)RETURN
      IF(prev_row .EQ. -1) THEN
         col_map = hfx_matrix%invalid_col
      ELSE
         !! Invalidate the previously set column mapping
         DO col = 1, SIZE(hfx_matrix%rows(prev_row)%key)
             col_map(hfx_matrix%rows(prev_row)%key(col)) = hfx_matrix%invalid_col
         END DO
      END IF

      DO col = 1, SIZE(hfx_matrix%rows(row)%key)
          col_map(hfx_matrix%rows(row)%key(col)) = col
      END DO

   END SUBROUTINE hfx_libheric_get_column_index_map

   SUBROUTINE hfx_libheric_queue_init(queue, num_spots)
      TYPE(hfx_libheric_queue), INTENT(INOUT) :: queue 
      INTEGER, INTENT(IN)            :: num_spots
      ALLOCATE(queue%queue(num_spots))
      queue%read_ptr = 1
      queue%write_ptr = 1
      queue%count = 0
      queue%size = num_spots
   END SUBROUTINE hfx_libheric_queue_init
   
   FUNCTION hfx_libheric_queue_empty(queue)
      LOGICAL                     :: hfx_libheric_queue_empty
      TYPE(hfx_libheric_queue), INTENT(IN) :: queue
      hfx_libheric_queue_empty = queue%count .EQ. 0
   END FUNCTION hfx_libheric_queue_empty
   
   FUNCTION hfx_libheric_queue_full(queue)
      LOGICAL                     :: hfx_libheric_queue_full
      TYPE(hfx_libheric_queue), INTENT(IN) :: queue
      hfx_libheric_queue_full = queue%count .EQ. queue%size
   END FUNCTION hfx_libheric_queue_full

   FUNCTION hfx_libheric_queue_pop(queue, pop_val)
      LOGICAL                        :: hfx_libheric_queue_pop
      TYPE(hfx_libheric_queue), INTENT(INOUT) :: queue
      INTEGER, INTENT(OUT)           :: pop_val
      hfx_libheric_queue_pop = .NOT. hfx_libheric_queue_empty(queue)
      IF(hfx_libheric_queue_pop) THEN
         pop_val = queue%queue(queue%read_ptr)
         queue%read_ptr = MODULO(queue%read_ptr, queue%size) + 1
         queue%count = queue%count - 1
      END IF
   END FUNCTION hfx_libheric_queue_pop
   
   FUNCTION hfx_libheric_queue_push(queue, push_val)
      LOGICAL                        :: hfx_libheric_queue_push
      TYPE(hfx_libheric_queue), INTENT(INOUT) :: queue
      INTEGER, INTENT(IN)            :: push_val
      hfx_libheric_queue_push = .NOT. hfx_libheric_queue_full(queue)
      IF(hfx_libheric_queue_push) THEN
         queue%queue(queue%write_ptr) = push_val
         queue%write_ptr = MODULO(queue%write_ptr, queue%size) + 1
         queue%count = queue%count + 1
      END IF
   END FUNCTION hfx_libheric_queue_push
   
   SUBROUTINE hfx_libheric_queue_release(queue)
      TYPE(hfx_libheric_queue), INTENT(INOUT) :: queue
      IF(ALLOCATED(queue%queue))DEALLOCATE(queue%queue)
   END SUBROUTINE 

! **************************************************************************************************
!> \brief initialize hfx_libheric_digestion_state by allocating buffers, however all initial values,
!>        pointers, indices are invalid. Call hfx_libheric_digestion_state_set_angular_momentum to
!>        bring digestion state into valid configuration
!> \param digestion_state - ... (output)
!> \param k_matrix - k_matrix replicated on all MPI ranks
!> \param l_max - maximum angular momentum
!> \param thread_id - OMP-thread-id of thread that is supposed to use the digestion state
!> \param tile_size - used to determine size of kcache; meaningless right now, but should be > 1
!> \param col_map - colum index map for requested row (output)
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_init_digestion_state(digestion_state, k_matrix, l_max, thread_id, tile_size)
      TYPE(hfx_libheric_digestion_state), INTENT(INOUT) :: digestion_state
      TYPE(hfx_libheric_sparse_matrix), INTENT(IN)      :: k_matrix  
      INTEGER, INTENT(IN)                               :: l_max, thread_id, tile_size

      INTEGER                                  :: num_entries, max_atom_block_size, i

      ! 2 * tile_size for ac, bc block, 2 * 4 for ad, bd blocks that are switched out regularly
      num_entries = 2 * tile_size + 2 * 4
      CALL hfx_libheric_queue_init(digestion_state%free_kblocks, num_entries)
      DO i = 1, num_entries
         IF(.NOT. hfx_libheric_queue_push(digestion_state%free_kblocks, i)) &
             CPABORT("Push into full queue!")
      END DO
      
      max_atom_block_size = MAXVAL(k_matrix%atom_ll_block_size)
      digestion_state%cache_block_size = max_atom_block_size 
      ! Also allocate space for a zero block and for temporary k buffers
      ALLOCATE(digestion_state%data(num_entries * max_atom_block_size + 4 * (nso(l_max) ** 4)))
      digestion_state%k_tmp_offset = num_entries * max_atom_block_size
      digestion_state%k_ac => NULL()
      digestion_state%k_ad => NULL()
      digestion_state%k_bc => NULL()
      digestion_state%k_bd => NULL()

      ALLOCATE(digestion_state%kcache(num_entries))
      DO i = 1, num_entries
         digestion_state%kcache(i)%hash = -1
         digestion_state%kcache(i)%iatom = -1
         digestion_state%kcache(i)%katom_col = -1
         digestion_state%kcache(i)%la = -1
         digestion_state%kcache(i)%lc = -1
         digestion_state%kcache(i)%atom_block_size = -1
         digestion_state%kcache(i)%left_block = digestion_state%null_block
         digestion_state%kcache(i)%right_block = digestion_state%null_block
      END DO
      digestion_state%first_block = digestion_state%null_block
      digestion_state%last_block = digestion_state%null_block

      !! SIZE(k_matrix%rows, 1) == natom
      ALLOCATE(digestion_state%pa_cols(SIZE(k_matrix%rows, 1)))
      ALLOCATE(digestion_state%pb_cols(SIZE(k_matrix%rows, 1)))
      ALLOCATE(digestion_state%ka_cols(SIZE(k_matrix%rows, 1)))
      ALLOCATE(digestion_state%kb_cols(SIZE(k_matrix%rows, 1)))
      
      digestion_state%p(1)%cols => digestion_state%pa_cols(:)
      digestion_state%p(2)%cols => digestion_state%pa_cols(:)
      digestion_state%p(3)%cols => digestion_state%pb_cols(:)
      digestion_state%p(4)%cols => digestion_state%pb_cols(:)
      
      digestion_state%k(1)%cols => digestion_state%ka_cols(:)
      digestion_state%k(2)%cols => digestion_state%ka_cols(:)
      digestion_state%k(3)%cols => digestion_state%kb_cols(:)
      digestion_state%k(4)%cols => digestion_state%kb_cols(:)
     
      digestion_state%na = -1
      digestion_state%nb = -1
      digestion_state%nc = -1
      digestion_state%nd = -1
      
      digestion_state%lock_id = -thread_id

   END SUBROUTINE hfx_libheric_init_digestion_state

! **************************************************************************************************
!> \brief initialize digestion state for quartet class {lalb|lcld}, only call if digestion state
!>        is either invalid or emtpy, i.e. after initialization or after full flush
!> \param digestion_state - ... (output)
!> \param la,lb,lc,ld - quartet_class
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_digestion_state_set_angular_momentum(digestion_state, la, lb, lc, ld)
      TYPE(hfx_libheric_digestion_state), INTENT(INOUT) :: digestion_state
      INTEGER, INTENT(IN)                               :: la, lb, lc, ld

      INTEGER                                           :: las(4), lcs(4), i, off

      las = (/ la, la, lb, lb /)
      lcs = (/ lc, ld, lc, ld /)

      CPASSERT(digestion_state%free_kblocks%count .EQ. digestion_state%free_kblocks%size)

      DO i = 1, 4
         digestion_state%p(i)%la = las(i)
         digestion_state%p(i)%lc = lcs(i)
         digestion_state%p(i)%iatom = -1
         digestion_state%p(i)%katom = -1
         digestion_state%p(i)%atom_block_size = -1
         digestion_state%p(i)%n_ac = nso(las(i)) * nso(lcs(i))
         digestion_state%p(i)%data => NULL()
         
         digestion_state%k(i)%la = las(i)
         digestion_state%k(i)%lc = lcs(i)
         digestion_state%k(i)%iatom = -1
         digestion_state%k(i)%katom = -1
         digestion_state%k(i)%atom_block_size = -1
         digestion_state%k(i)%n_ac = nso(las(i)) * nso(lcs(i))
         digestion_state%k(i)%data => NULL()

         !! block hash is encoded as | --- blk --- | --- lc --- | --- la --- | --- katom --- | --- iatom --- |
         !!                          | <- 2-bit -> | <- 4-bit ->| <- 4-bit ->| <-  24-bit -> | <-  24-bit -> |
         !! blk indicates whether a block is ac, ad, bc, or bd, i.e. we treat them seperately even if they technically alias
         digestion_state%k(i)%hash_postfix = IOR(INT(digestion_state%k(i)%la, int_8), ISHFT(INT(digestion_state%k(i)%lc, int_8), 4))
         digestion_state%k(i)%hash_postfix = IOR(digestion_state%k(i)%hash_postfix, ISHFT(INT(i - 1, int_8), 8))
         digestion_state%k(i)%hash_postfix = ISHFT(digestion_state%k(i)%hash_postfix, 48)

      END DO

      digestion_state%prev_iatom = -1
      digestion_state%prev_jatom = -1
      
      digestion_state%la = la 
      digestion_state%lb = lb 
      digestion_state%lc = lc 
      digestion_state%ld = ld 
      
      digestion_state%na = nso(la) 
      digestion_state%nb = nso(lb) 
      digestion_state%nc = nso(lc) 
      digestion_state%nd = nso(ld)

      off = digestion_state%k_tmp_offset
      digestion_state%k_ac => digestion_state%data(off + 1:off + digestion_state%k(1)%n_ac)
      off = off + digestion_state%k(1)%n_ac
      digestion_state%k_ad => digestion_state%data(off + 1:off + digestion_state%k(2)%n_ac)
      off = off + digestion_state%k(2)%n_ac
      digestion_state%k_bc => digestion_state%data(off + 1:off + digestion_state%k(3)%n_ac)
      off = off + digestion_state%k(3)%n_ac
      digestion_state%k_bd => digestion_state%data(off + 1:off + digestion_state%k(4)%n_ac)
      
   END SUBROUTINE hfx_libheric_digestion_state_set_angular_momentum

! **************************************************************************************************
!> \brief flush requested amount of kcache blocks back into k_matrix
!> \param ds - digestion_state
!> \param k_matrix - k_matrix replicated on all MPI ranks
!> \param flush_count - number of blocks to flush (for full flush use a value >= ds%free_blocks%size)
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_digestion_state_flush(ds, k_matrix, flush_count)
      TYPE(hfx_libheric_digestion_state), INTENT(INOUT) :: ds
      TYPE(hfx_libheric_sparse_matrix), INTENT(INOUT)   :: k_matrix  
      INTEGER, INTENT(IN)                               :: flush_count

      INTEGER                                 :: flushes_remaining, num_flagged, curr_block, i, num_in_cache, unlocked_list_size
      INTEGER(int_8)                          :: eviction_list(2, flush_count), unlocked_list(flush_count)
      LOGICAL                                 :: offset_locked

      num_in_cache = ds%free_kblocks%size - ds%free_kblocks%count
      flushes_remaining = MIN(flush_count, num_in_cache)
      DO WHILE(flushes_remaining .GT. 0)
         !! Iterate through cache in reverse for LRU eviction
         curr_block = ds%last_block
         num_flagged = 0
         !$OMP CRITICAL
         DO WHILE(curr_block .NE. ds%null_block .AND. num_flagged .LT. flushes_remaining)
            IF(num_flagged .GE. flushes_remaining) EXIT
            ASSOCIATE(blk => ds%kcache(curr_block))
            !! an offset < 1 indicates that the requested block is currently locked
            !! Each thread locks with its negative thread_id
            IF(k_matrix%rows(blk%iatom)%offset(blk%katom_col, blk%lc, blk%la) .GT. 0) THEN
               num_flagged = num_flagged + 1
               eviction_list(1, num_flagged) = k_matrix%rows(blk%iatom)%offset(blk%katom_col, blk%lc, blk%la)
               eviction_list(2, num_flagged) = curr_block
            END IF
            END ASSOCIATE
            curr_block = ds%kcache(curr_block)%left_block
         END DO

         !! Lock offset only *after* collecting all eviction candidates to avoid writing wrong offsets into the eviction list
         !! , e.g. if to blocks point to the same offset, the loop above would overwrite the correct offset for the second block
         DO i = 1, num_flagged
            ASSOCIATE(blk => ds%kcache(eviction_list(2, i)))
               k_matrix%rows(blk%iatom)%offset(blk%katom_col, blk%lc, blk%la) = ds%lock_id
            END ASSOCIATE
         END DO
         !$OMP END CRITICAL
         DO i = 1, num_flagged
            ASSOCIATE(k_offset => eviction_list(1, i), blk => ds%kcache(eviction_list(2, i)))
            k_matrix%data(k_offset:k_offset + blk%atom_block_size - 1) = k_matrix%data(k_offset:k_offset + blk%atom_block_size - 1)&
                                                                       + ds%data(blk%offset:blk%offset + blk%atom_block_size - 1)
            END ASSOCIATE
         END DO

         unlocked_list_size = 0
         !! Unlock offset only *after* writing back all data, because a thread might have locked an offset multiple times
         !! Also no offset should be unlocked twice, hence the "unlocked_list" to keep track of already unlocked offsets
         DO i = 1, num_flagged
            ASSOCIATE(k_offset => eviction_list(1, i), blk_id => eviction_list(2, i), blk => ds%kcache(eviction_list(2, i)))
            offset_locked = .NOT. ANY(unlocked_list(1:unlocked_list_size) .EQ. k_offset)

            IF(offset_locked) THEN
               !! Unlock by writing back correct offset
               k_matrix%rows(blk%iatom)%offset(blk%katom_col, blk%lc, blk%la) = k_offset
               unlocked_list_size = unlocked_list_size + 1
               unlocked_list(unlocked_list_size) = k_offset
            END IF
            
            !! Remove block from cache and invalidate
            blk%hash = -1
            blk%iatom = -1
            blk%katom_col = -1
            blk%la = -1
            blk%lc = -1
            blk%offset = -1
            blk%atom_block_size = -1
            IF(blk_id .EQ. ds%last_block .AND. blk_id .EQ. ds%first_block) THEN
               ds%first_block = ds%null_block
               ds%last_block = ds%null_block
            ELSE IF(blk_id .EQ. ds%last_block) THEN
               ds%last_block = blk%left_block
               ds%kcache(blk%left_block)%right_block = ds%null_block
            ELSE IF(eviction_list(2, i) .EQ. ds%first_block) THEN
               ds%first_block = blk%right_block
               ds%kcache(blk%right_block)%left_block = ds%null_block
            ELSE
               ds%kcache(blk%left_block)%right_block = blk%right_block  
               ds%kcache(blk%right_block)%left_block = blk%left_block  
            END IF
            blk%left_block = ds%null_block
            blk%right_block = ds%null_block
            IF(.NOT. hfx_libheric_queue_push(ds%free_kblocks, INT(blk_id))) &
                CPABORT("Push into full queue")

            flushes_remaining = flushes_remaining - 1
            END ASSOCIATE
         END DO
      END DO

   END SUBROUTINE hfx_libheric_digestion_state_flush

! **************************************************************************************************
!> \brief fetch p_matrix blocks needed for digestion of atom quartet (iatom,jatom|katom,latom)
!> \param ds - digestion_state
!> \param p_matrix - p_matrix replicated on all MPI ranks
!> \param kind_of - map of atom index to atom kind
!> \param iatom,jatom,katom,latom - requested atom quartet
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_digestion_state_get_pblocks(ds, p_matrix, kind_of, iatom, jatom, katom, latom)
      TYPE(hfx_libheric_digestion_state), INTENT(INOUT) :: ds
      TYPE(hfx_libheric_sparse_matrix), INTENT(IN)      :: p_matrix  
      INTEGER, DIMENSION(:), INTENT(IN)                 :: kind_of
      INTEGER, INTENT(IN)                               :: iatom, jatom, katom, latom
      
      INTEGER                                  :: iatoms(4), katoms(4), i
      INTEGER(int_8)                           :: offset, screening_size
      
      iatoms = (/ iatom, iatom, jatom, jatom /)
      katoms = (/ katom, latom, katom, latom /)

      DO i = 1, 4
         ASSOCIATE(p => ds%p(i))
         p%iatom = iatoms(i)
         p%katom = katoms(i)
         p%ikind = kind_of(p%iatom)
         p%kkind = kind_of(p%katom)
         p%required = p%cols(p%katom) .NE. p_matrix%invalid_col
         p%present = .TRUE.
         p%atom_block_size = p_matrix%atom_ll_block_size(p%kkind, p%ikind, p%lc, p%la)
         p%shellset_offset => p_matrix%shellset_offset(:, :, :, :, p%kkind, p%ikind, p%lc, p%la)
         screening_size = p_matrix%screening_info_size(p%kkind, p%ikind, p%lc, p%la)

         IF(p%required) THEN
             offset = p_matrix%rows(p%iatom)%offset(p%cols(p%katom), p%lc, p%la)
             p%data => p_matrix%data(offset:offset + p%atom_block_size + screening_size - 1)
         ELSE
             p%data => p_matrix%zero_block(1:p%atom_block_size + screening_size)
         END IF
         END ASSOCIATE
      END DO
      
   END SUBROUTINE hfx_libheric_digestion_state_get_pblocks
  
! **************************************************************************************************
!> \brief fetch k_matrix blocks needed for digestion of atom quartet (iatom,jatom|katom,latom); does
!>        not work properly without calling hfx_libheric_digestion_state_get_pblocks
!> \param ds - digestion_state
!> \param k_matrix - k_matrix replicated on all MPI ranks
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_digestion_state_get_kblocks(ds, k_matrix)
      TYPE(hfx_libheric_digestion_state), INTENT(INOUT) :: ds
      TYPE(hfx_libheric_sparse_matrix), INTENT(INOUT)   :: k_matrix  
      
      INTEGER                                  :: i, num_requested, num_hits, num_missing, &
                                                  curr_block, next_block, new_block, offset
      num_requested = 0
      DO i = 1, 4
         !! p_other is the p_block that the k_block will be matched with, e.g. k(2) == k_ad, and p(3) == p_bc
         ASSOCIATE(k => ds%k(i), p => ds%p(i), p_other => ds%p(4 - (i - 1)))
         k%iatom = p%iatom
         k%katom = p%katom
         k%ikind = p%ikind
         k%kkind = p%kkind
         k%required = p_other%required 
         k%atom_block_size = p%atom_block_size
         k%shellset_offset => k_matrix%shellset_offset(:, :, :, :, k%kkind, k%ikind, k%lc, k%la)
          
         IF(k%required) THEN
            num_requested = num_requested + 1
            k%present = .FALSE.
         ELSE
            k%data => k_matrix%zero_block(1:k%atom_block_size)
            k%present = .TRUE.
         END IF
         k%hash = IOR(IOR(INT(k%iatom, int_8), ISHFT(INT(k%katom, int_8), 24 )), k%hash_postfix)
         END ASSOCIATE
      END DO

      num_hits = 0
      !! Search for required k_blocks in kcache
      curr_block = ds%first_block
      DO WHILE(curr_block .NE. ds%null_block .AND. num_hits .LT. num_requested) 
         next_block = ds%kcache(curr_block)%right_block
         DO i = 1, 4
            ASSOCIATE(k => ds%k(i))
            IF(.NOT. k%present .AND. k%hash .EQ. ds%kcache(curr_block)%hash) THEN
               !! Cache hit
               num_hits = num_hits + 1
               offset = (curr_block - 1) * ds%cache_block_size
               k%data => ds%data(offset + 1:offset + k%atom_block_size)

               !! If curr_block is already the first block nothing has to be done
               IF(curr_block .NE. ds%first_block) THEN
                  !! Connect left and right block
                  IF(curr_block .NE. ds%last_block) THEN
                     ds%kcache(ds%kcache(curr_block)%left_block)%right_block = ds%kcache(curr_block)%right_block
                     ds%kcache(ds%kcache(curr_block)%right_block)%left_block = ds%kcache(curr_block)%left_block
                  ELSE !! If the block is the last block, the last pointer and the new last block need to be updated
                     ds%kcache(ds%kcache(curr_block)%left_block)%right_block = ds%null_block
                     ds%last_block = ds%kcache(curr_block)%left_block
                  END IF
                  !! Append to the front
                  ds%kcache(ds%first_block)%left_block = curr_block 
                  ds%kcache(curr_block)%left_block = ds%null_block
                  ds%kcache(curr_block)%right_block = ds%first_block
                  ds%first_block = curr_block
                  
               END IF
               k%present = .TRUE.
               EXIT
            END IF
            END ASSOCIATE
         END DO
         curr_block = next_block
      END DO

      num_missing = num_requested - num_hits
      IF(ds%free_kblocks%count .LT. num_missing) THEN
         CALL hfx_libheric_digestion_state_flush(ds, k_matrix, num_missing)
      END IF

      DO i = 1, 4
         ASSOCIATE(k => ds%k(i))
         IF(.NOT. k%present) THEN
            !! Sanity check
            IF(.NOT. hfx_libheric_queue_pop(ds%free_kblocks, new_block)) &
                CPABORT("Empty free spot queue!")
            ds%kcache(new_block)%hash = k%hash
            ds%kcache(new_block)%iatom = k%iatom
            ds%kcache(new_block)%katom_col = k%cols(k%katom)
            ds%kcache(new_block)%la = k%la
            ds%kcache(new_block)%lc = k%lc
            ds%kcache(new_block)%offset = (new_block - 1) * ds%cache_block_size + 1
            ds%kcache(new_block)%atom_block_size = k%atom_block_size
            ds%kcache(new_block)%left_block = ds%null_block
            IF(ds%first_block .EQ. ds%null_block) THEN
                ds%last_block = new_block
                ds%kcache(new_block)%right_block = ds%null_block
            ELSE
                ds%kcache(new_block)%right_block = ds%first_block
                ds%kcache(ds%first_block)%left_block = new_block
            END IF
            ds%first_block = new_block

            k%present = .TRUE.
            k%data => ds%data(ds%kcache(new_block)%offset:ds%kcache(new_block)%offset + k%atom_block_size - 1)
            k%data = 0.0_dp
         END IF
         END ASSOCIATE
      END DO
      
   END SUBROUTINE hfx_libheric_digestion_state_get_kblocks

! **************************************************************************************************
!> \brief deallocate digestion_state 
!> \param digestion_state - ...
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_free_digestion_state(digestion_state)
      TYPE(hfx_libheric_digestion_state), INTENT(INOUT) :: digestion_state

      INTEGER :: i
        
      DO i = 1, 4
         digestion_state%p(i)%data => NULL()
         digestion_state%p(i)%cols => NULL()
         digestion_state%p(i)%shellset_offset => NULL()
         
         digestion_state%k(i)%data => NULL()
         digestion_state%k(i)%cols => NULL()
         digestion_state%k(i)%shellset_offset => NULL()
      END DO
      digestion_state%k_ac => NULL()
      digestion_state%k_ad => NULL()
      digestion_state%k_bc => NULL()
      digestion_state%k_bd => NULL()

      IF(ALLOCATED(digestion_state%kcache)) DEALLOCATE(digestion_state%kcache)
      IF(ASSOCIATED(digestion_state%pa_cols)) DEALLOCATE(digestion_state%pa_cols)
      IF(ASSOCIATED(digestion_state%pb_cols)) DEALLOCATE(digestion_state%pb_cols)
      IF(ASSOCIATED(digestion_state%ka_cols)) DEALLOCATE(digestion_state%ka_cols)
      IF(ASSOCIATED(digestion_state%kb_cols)) DEALLOCATE(digestion_state%kb_cols)
      IF(ASSOCIATED(digestion_state%data)) DEALLOCATE(digestion_state%data)
      CALL hfx_libheric_queue_release(digestion_state%free_kblocks)
   END SUBROUTINE hfx_libheric_free_digestion_state

! **************************************************************************************************
!> \brief digest (contract) libheric atom_quartet into k_matrix 
!> \param ds - digestion state
!> \param atom_quartet - computed by libheric; contains multiple shellset quartets of a specific 
!>                       quartet class
!> \param nkind - number of unique basis sets 
!> \param l_max - maximum angular momentum
!> \param p_matrix - p_matrix replicated on all MPI ranks
!> \param k_matrix - k_matrix replicated on all MPI ranks
!> \param basis_parameter_l - all hfx_libheric basis sets
!> \param kind_of - map of atom index to atom kind
!> \param fac - scaling factor (hfx contribution fraction)
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_digest_atom_quartet(ds, atom_quartet, nkind, l_max, p_matrix, k_matrix, &
                                               basis_parameter_l, kind_of, fac)
      TYPE(hfx_libheric_digestion_state), INTENT(INOUT)  :: ds
      TYPE(atom_quartet_type), INTENT(IN)                :: atom_quartet
      INTEGER, INTENT(IN)                                :: nkind, l_max
      TYPE(hfx_libheric_sparse_matrix), INTENT(IN)       :: p_matrix
      TYPE(hfx_libheric_sparse_matrix), INTENT(INOUT)    :: k_matrix
      TYPE(hfx_libheric_basis_type), & 
          DIMENSION(nkind, 0:l_max), INTENT(IN)          :: basis_parameter_l
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of
      REAL(dp), INTENT(IN)                               :: fac

      INTEGER                                            :: iset, jset, kset, lset, &
                                                            ishl, jshl, kshl, lshl, &
                                                            shlset
      INTEGER(int_1)                                     :: shellset_id(8)
      REAL(dp)                                           :: symm_fac

      ASSOCIATE(iatom => atom_quartet%quartet_id(1) + 1, jatom => atom_quartet%quartet_id(2) + 1, &
                katom => atom_quartet%quartet_id(3) + 1, latom => atom_quartet%quartet_id(4) + 1)
      
      !! Must also check equivalence of angular momentum! 
      symm_fac = 0.5_dp
      IF (iatom == jatom) symm_fac = symm_fac*2.0_dp
      IF (katom == latom) symm_fac = symm_fac*2.0_dp
      IF (iatom == katom .AND. jatom == latom .AND. &
          iatom /= jatom .AND. katom /= latom) symm_fac = symm_fac*2.0_dp
      IF (iatom == katom .AND. iatom == jatom .AND. katom == latom) symm_fac = symm_fac*2.0_dp
      symm_fac = fac/symm_fac
      
      CALL hfx_libheric_get_column_index_map(p_matrix, iatom, ds%prev_iatom, ds%pa_cols)
      CALL hfx_libheric_get_column_index_map(p_matrix, jatom, ds%prev_jatom, ds%pb_cols)
      CALL hfx_libheric_get_column_index_map(k_matrix, iatom, ds%prev_iatom, ds%ka_cols)
      CALL hfx_libheric_get_column_index_map(k_matrix, jatom, ds%prev_jatom, ds%kb_cols)

      ds%prev_iatom = iatom
      ds%prev_jatom = jatom

      CALL hfx_libheric_digestion_state_get_pblocks(ds, p_matrix, kind_of, iatom, jatom, katom, latom)
      CALL hfx_libheric_digestion_state_get_kblocks(ds, k_matrix)

      DO shlset = 1, SIZE(atom_quartet%shellset_ids, 2)
         shellset_id = atom_quartet%shellset_ids(:, shlset)
         iset = shellset_id(1) + 1
         jset = shellset_id(2) + 1  
         kset = shellset_id(3) + 1
         lset = shellset_id(4) + 1
         ishl = shellset_id(5) + 1 
         jshl = shellset_id(6) + 1
         kshl = shellset_id(7) + 1
         lshl = shellset_id(8) + 1
         ASSOCIATE(isetg => basis_parameter_l(kind_of(iatom), ds%la)%global_set(iset), &
                   jsetg => basis_parameter_l(kind_of(jatom), ds%lb)%global_set(jset), &
                   ksetg => basis_parameter_l(kind_of(katom), ds%lc)%global_set(kset), &
                   lsetg => basis_parameter_l(kind_of(latom), ds%ld)%global_set(lset))
        
         !! shelset offset is 0-indexed
         ASSOCIATE(pac_off => ds%p(1)%shellset_offset(kshl, ishl, kset, iset), &
                   pad_off => ds%p(2)%shellset_offset(lshl, ishl, lset, iset), &
                   pbc_off => ds%p(3)%shellset_offset(kshl, jshl, kset, jset), &
                   pbd_off => ds%p(4)%shellset_offset(lshl, jshl, lset, jset), &
                   kac_off => ds%k(1)%shellset_offset(kshl, ishl, kset, iset), &
                   kad_off => ds%k(2)%shellset_offset(lshl, ishl, lset, iset), &
                   kbc_off => ds%k(3)%shellset_offset(kshl, jshl, kset, jset), &
                   kbd_off => ds%k(4)%shellset_offset(lshl, jshl, lset, jset))

         CALL contract_block(ds%na, ds%nb, ds%nc, ds%nd, &
                             ds%k_bd, ds%k_bc, ds%k_ad, ds%k_ac, &
                             ds%p(4)%data(pbd_off + 1:pbd_off + ds%p(4)%n_ac), ds%p(3)%data(pbc_off + 1:pbc_off + ds%p(3)%n_ac), &
                             ds%p(2)%data(pad_off + 1:pad_off + ds%p(2)%n_ac), ds%p(1)%data(pac_off + 1:pac_off + ds%p(1)%n_ac), &
                             atom_quartet%integrals(:, shlset), symm_fac)
         !! Only necessary because contract_block zeros out the k buffers every time
         !! Otherwise ds%k(i) could be passed directly to contract_block
         ds%k(1)%data(kac_off + 1:kac_off + ds%k(1)%n_ac) = ds%k(1)%data(kac_off + 1:kac_off + ds%k(1)%n_ac) &
                                                          + ds%k_ac(1:ds%k(1)%n_ac)
         ds%k(2)%data(kad_off + 1:kad_off + ds%k(2)%n_ac) = ds%k(2)%data(kad_off + 1:kad_off + ds%k(2)%n_ac) &
                                                          + ds%k_ad(1:ds%k(2)%n_ac)
         ds%k(3)%data(kbc_off + 1:kbc_off + ds%k(3)%n_ac) = ds%k(3)%data(kbc_off + 1:kbc_off + ds%k(3)%n_ac) &
                                                          + ds%k_bc(1:ds%k(3)%n_ac)
         ds%k(4)%data(kbd_off + 1:kbd_off + ds%k(4)%n_ac) = ds%k(4)%data(kbd_off + 1:kbd_off + ds%k(4)%n_ac) &
                                                          + ds%k_bd(1:ds%k(4)%n_ac)
         
         END ASSOCIATE
         END ASSOCIATE
      END DO
      END ASSOCIATE
   END SUBROUTINE hfx_libheric_digest_atom_quartet

END MODULE hfx_communication
