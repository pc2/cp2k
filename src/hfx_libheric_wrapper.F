MODULE hfx_libheric_wrapper
   USE hfx_types, ONLY: hfx_libheric_basis_type
#if defined(__LIBHERIC)
   USE libheric_interface, ONLY: hfx_libheric_atom_quartet => atom_quartet_type, &
                                 hfx_libheric_acc_kind => acc_kind, &
                                 hfx_libheric_acc_device => acc_device, &
                                 hfx_libheric_initialize => libheric_initialize, &
                                 hfx_libheric_release => libheric_release, &
                                 hfx_libheric_query_accelerator_devices => query_accelerator_devices, &
                                 hfx_libheric_allocate_accelerator_contexts => create_accelerator_contexts, &
                                 hfx_libheric_free_accelerator_contexts => free_accelerator_contexts, &
                                 hfx_libheric_is_functional => accelerator_context_is_functional, &
                                 hfx_libheric_set_potential_operator => accelerator_context_set_potential_operator, &
                                 hfx_libheric_set_quartet_class => accelerator_context_set_quartet_class, &
                                 hfx_libheric_set_atom_ids => accelerator_context_set_atom_ids, &
                                 hfx_libheric_set_set_ids => accelerator_context_set_set_ids, &
                                 hfx_libheric_set_shell_ids => accelerator_context_set_shell_ids, &
                                 hfx_libheric_set_B_offset => accelerator_context_set_B_offset, &
                                 hfx_libheric_set_D_offset => accelerator_context_set_D_offset, &
                                 hfx_libheric_set_exponent_ids => accelerator_context_set_exponent_ids, &
                                 hfx_libheric_set_Q_offset => accelerator_context_set_Q_offset, &
                                 hfx_libheric_submit_peri => accelerator_context_submit_peri, &
                                 hfx_libheric_finish_quartet_class => accelerator_context_finish_quartet_class, &
                                 hfx_libheric_get_next_staged_atom_quartet => accelerator_context_get_next_staged_atom_quartet, &
                                 hfx_libheric_prefetch_cached_atom_quartet => accelerator_context_prefetch_cached_atom_quartet, &
                                 hfx_libheric_update_positions => accelerator_context_update_positions, &
                                 hfx_libheric_update_cell => accelerator_context_update_cell, &
                                 hfx_libheric_optimize_screening_function => optimize_screening_function, &
                                 LIBHERIC_MAX_AM
                                 
   USE, INTRINSIC :: iso_c_binding, ONLY: C_BOOL, hfx_libheric_accelerator_context => C_PTR
#else
   USE input_constants,                 ONLY: &
        do_potential_TShPSC, do_potential_coulomb, do_potential_gaussian, do_potential_long, &
        do_potential_mix_cl, do_potential_mix_cl_trunc, do_potential_mix_lg, do_potential_short, &
        do_potential_truncated
!$ USE OMP_LIB, ONLY: omp_get_max_threads
#endif
   USE particle_types, ONLY: particle_type
   USE kinds, ONLY: dp, &
                    int_1, &
                    int_8, &
                    default_string_length
    
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: hfx_libheric_accelerator_context, &
             hfx_libheric_atom_quartet, &
             hfx_libheric_acc_kind, &
             hfx_libheric_acc_device, &
             hfx_libheric_initialize, &
             hfx_libheric_release, &
             hfx_libheric_query_accelerator_devices, &
             hfx_libheric_create_accelerator_contexts, &
             hfx_libheric_free_accelerator_contexts, &
             hfx_libheric_is_functional, &
             hfx_libheric_set_potential_operator, &
             hfx_libheric_set_quartet_class, &
             hfx_libheric_set_atom_ids, &
             hfx_libheric_set_set_ids, &
             hfx_libheric_set_shell_ids, &
             hfx_libheric_set_B_offset, &
             hfx_libheric_set_D_offset, &
             hfx_libheric_set_exponent_ids, &
             hfx_libheric_set_Q_offset, &
             hfx_libheric_submit_peri, &
             hfx_libheric_finish_quartet_class, &
             hfx_libheric_get_next_staged_atom_quartet, &
             hfx_libheric_prefetch_cached_atom_quartet, &
             hfx_libheric_update_positions, &
             hfx_libheric_update_cell, &
             hfx_libheric_optimize_screening_function

#if !defined (__LIBHERIC)

   TYPE hfx_libheric_accelerator_context
      TYPE(hfx_libheric_basis_type), DIMENSION(:, :), POINTER :: basis_parameter_l => NULL()
      REAL(dp), DIMENSION(:, :), ALLOCATABLE                  :: atom_positions
      INTEGER, DIMENSION(:), ALLOCATABLE                      :: atom_kinds
      REAL(dp), DIMENSION(3, 3)                               :: cell_mat

      LOGICAL                                                 :: functional = .FALSE.
      INTEGER                                                 :: l_max = -1, op_id = -1
      REAL(dp)                                                :: cutoff_radius = -1.0_dp, omega = -1.0_dp
   END TYPE hfx_libheric_accelerator_context

   TYPE hfx_libheric_acc_kind
      CHARACTER(len=default_string_length) :: name = ""
      CHARACTER(len=default_string_length) :: type = ""
      CHARACTER(len=default_string_length) :: library = ""
   END TYPE hfx_libheric_acc_kind

   TYPE hfx_libheric_acc_device
      CHARACTER(len=default_string_length) :: name = ""
      CHARACTER(len=default_string_length) :: bdf = ""
      CHARACTER(len=default_string_length) :: cpu_name = ""
      INTEGER                              :: num_cores = 0
      TYPE(hfx_libheric_acc_kind)          :: kind
   END TYPE hfx_libheric_acc_device
   
   TYPE hfx_libheric_atom_quartet
       INTEGER                                      :: quartet_id(4), n_shellset_quartets
       REAL(dp), DIMENSION(:, :), ALLOCATABLE       :: integrals
       INTEGER(int_1), DIMENSION(:, :), ALLOCATABLE :: shellset_ids
       REAL(dp), DIMENSION(:), ALLOCATABLE          :: pmax_set_vals
   END TYPE hfx_libheric_atom_quartet

   TYPE hfx_libheric_quartet_staging
       TYPE(hfx_libheric_atom_quartet), DIMENSION(:), ALLOCATABLE :: data
       INTEGER                                                    :: read_idx, write_idx, i_staging, n_staging  
   END TYPE

#endif
CONTAINS
   SUBROUTINE hfx_libheric_create_accelerator_contexts(max_set, nkind, l_max, basis_parameter_l, &
                                                       cell_mat, acc_dev, particle_set, & 
                                                       mem_available, compute_raw_cartesian, ctx_handles) 
      INTEGER, INTENT(IN)                                 :: max_set, nkind, l_max
      TYPE(hfx_libheric_basis_type), &
         DIMENSION(nkind, 0:l_max), TARGET, INTENT(IN)    :: basis_parameter_l
      REAL(dp), DIMENSION(3, 3), INTENT(IN)               :: cell_mat
      TYPE(hfx_libheric_acc_device), INTENT(IN)           :: acc_dev
      TYPE(particle_type), DIMENSION(:), INTENT(IN)       :: particle_set
      INTEGER(int_8), INTENT(IN)                          :: mem_available
      LOGICAL, INTENT(IN)                                 :: compute_raw_cartesian
      TYPE(hfx_libheric_accelerator_context), &
         DIMENSION(:), ALLOCATABLE, INTENT(OUT)           :: ctx_handles
#if defined(__LIBHERIC)
      INTEGER                                  :: i, ikind, iset, &
                                                  ishell, l, coeff_ctr, n_coeff, & 
                                                  iexp, exp_ctr, n_exp, natom
      LOGICAL(C_BOOL)                          :: internal_compute_raw_cartesian
      INTEGER, DIMENSION(:), ALLOCATABLE       :: atom_kinds
      INTEGER, DIMENSION(:, :), ALLOCATABLE    :: kind_nset
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: set_nexp, set_nshl
      REAL(dp), DIMENSION(:, :), ALLOCATABLE   :: atom_positions
      REAL(dp), DIMENSION(:), ALLOCATABLE      :: exponents, coefficients
      
      CPASSERT(l_max .LE. LIBHERIC_MAX_AM)
      natom = SIZE(particle_set)
      ALLOCATE(atom_positions(3, natom))
      ALLOCATE(atom_kinds(natom))
      DO i = 1, natom
          atom_positions(:, i) = particle_set(i)%r(:)
          atom_kinds(i) = particle_set(i)%atomic_kind%kind_number - 1 ! Make 0-indexed
          CPASSERT(particle_set(i)%atom_index .EQ. i)
      END DO
      
      ALLOCATE(kind_nset(nkind, 0:l_max))
      n_exp = 0
      n_coeff = 0
      DO l = 0, l_max
         DO ikind = 1, nkind
            ASSOCIATE(bp => basis_parameter_l(ikind, l))
            kind_nset(ikind, l) = bp%nset
            DO iset = 1, bp%nset
                n_exp = n_exp + bp%npgf(iset)
                n_coeff = n_coeff + bp%nshell(iset) * bp%npgf(iset)
            END DO
            END ASSOCIATE
         END DO
      END DO

      ALLOCATE(exponents(n_exp))
      ALLOCATE(coefficients(n_coeff))
      ALLOCATE(set_nexp(max_set, nkind, 0:l_max))
      ALLOCATE(set_nshl(max_set, nkind, 0:l_max))

      exp_ctr = 0
      coeff_ctr = 0
      DO l = 0, l_max
         DO ikind = 1, nkind
            ASSOCIATE(bp => basis_parameter_l(ikind, l))
            DO iset = 1, bp%nset
               set_nexp(iset, ikind, l) = bp%npgf(iset)
               set_nshl(iset, ikind, l) = bp%nshell(iset)
               DO iexp = 1, bp%npgf(iset)
                  exp_ctr = exp_ctr + 1
                  exponents(exp_ctr) = bp%zet(iexp, iset)
               END DO
               DO ishell = 1, bp%nshell(iset)
                  DO iexp = 1, bp%npgf(iset)
                     coeff_ctr = coeff_ctr + 1
                     coefficients(coeff_ctr) = bp%gcc(iexp, ishell, iset)
                  END DO
               END DO
            END DO
            END ASSOCIATE
         END DO
      END DO
      ALLOCATE(ctx_handles(acc_dev%num_cores))
      internal_compute_raw_cartesian = LOGICAL(compute_raw_cartesian, KIND=C_BOOL)
      CALL hfx_libheric_allocate_accelerator_contexts(ctx_handles(1), acc_dev, SIZE(ctx_handles), &
                                                      mem_available, &
                                                      internal_compute_raw_cartesian, &
                                                      cell_mat, atom_positions, atom_kinds, &
                                                      exponents, coefficients, &
                                                      kind_nset, set_nexp, set_nshl)
#else
      INTEGER :: i, j, natom
      natom = SIZE(particle_set)
      ALLOCATE(ctx_handles(acc_dev%num_cores))
      DO i = 1, SIZE(ctx_handles)
         ctx_handles(i)%basis_parameter_l => basis_parameter_l
         ALLOCATE(ctx_handles(i)%atom_positions(3, natom)) 
         ALLOCATE(ctx_handles(i)%atom_kinds(natom))
         DO j = 1, natom
             ctx_handles(i)%atom_positions(:, j) = particle_set(i)%r(:)
             ctx_handles(i)%atom_kinds(j) = particle_set(i)%atomic_kind%kind_number - 1 ! Make 0-indexed
             CPASSERT(particle_set(j)%atom_index .EQ. j)
         END DO
         ctx_handles(i)%cell_mat = cell_mat
         ctx_handles(i)%l_max = l_max
         ctx_handles(i)%functional = .FALSE.
      END DO
#endif
   END SUBROUTINE hfx_libheric_create_accelerator_contexts
#if !defined(__LIBHERIC)

   SUBROUTINE hfx_libheric_quartet_staging_reset(staging, la, lb, lc, ld, max_shellset_quartets)
       TYPE(hfx_libheric_quartet_staging), INTENT(INOUT) :: staging
       INTEGER, INTENT(IN)                               :: la, lb, lc, ld


       INTEGER :: i
       staging%n_staging = 2
       staging%i_staging = 0
       staging%read_ptr = 0
       staging%write_ptr = 0
       CALL hfx_libheric_quartet_staging_free(staging)
       ALLOCATE(staging%data(n_staging))
       DO i = 1, staging%n_staging
          ALLOCATE(staging%data(i)%integrals(nco(la) * nco(lb) * nco(lc) * nco(ld), max_shellset_quartets), SOURCE=0.0_dp)
          ALLOCATE(staging%data(i)%shellset_ids(8, max_shellset_quartets))
          ALLOCATE(staging%data(i)%pmax_set_vals(max_shellset_quartets))
       END DO
   END SUBROUTINE hfx_libheric_quartet_staging_reset 

   LOGICAL FUNCTION hfx_libheric_quartet_staging_full(staging)
       TYPE(hfx_libheric_quartet_staging), INTENT(IN) :: staging
       RETURN staging%i_staging .EQ. staging%n_staging
   END FUNCTION hfx_libheric_quartet_staging_full
   
   LOGICAL FUNCTION hfx_libheric_quartet_staging_empty(staging)
       TYPE(hfx_libheric_quartet_staging), INTENT(IN) :: staging
       RETURN staging%i_staging .EQ. 0
   END FUNCTION hfx_libheric_quartet_staging_empty

   INTEGER FUNCTION hfx_libheric_quartet_staging_reserve(staging, A, B, C, D)
      TYPE(hfx_libheric_quartet_staging), INTENT(INOUT) :: staging
      INTEGER, INTENT(IN)                               :: A, B, C, D

      IF(hfx_libheric_quartet_staging_full(staging))RETURN 0
      staging%write_ptr = MODULO(staging%write_ptr + 1, staging%n_staging)
      staging%i_staging = staging%i_staging + 1

      ASSOCIATE(curr_data => staging%data(staging%write_ptr + 1))
      IF(curr_data%n_shellset_quartets > 0) &
         curr_data%integrals(1:curr_data%n_shellset_quartets * staging%eris_allocated_per_shellset) = 0.0_dp
      
      curr_data%quartet_id(:) = [A, B, C, D]
      curr_data%n_shellset_quartets = 0
      curr_data%prims_requested = 0
      curr_data%prims_computed = 0
      END ASSOCIATE

      RETURN staging%write_ptr + 1

   END FUNCTION hfx_libheric_quartet_staging_reserve
   
   INTEGER FUNCTION hfx_libheric_quartet_staging_get_next_atom_quartet(staging)
      TYPE(hfx_libheric_quartet_staging), INTENT(INOUT) :: staging
      
      INTEGER :: new_read_idx
      IF(hfx_libheric_quartet_staging_empty(staging))RETURN 0

      new_read_ptr = MODULO(staging%read_ptr + 1, staging%n_staging)
      ASSOCIATE(candidate => staging%data(staging%read_ptr + 1))
      IF(candidate%prims_requested .NE. candidate%prims_computed)RETURN 0
      staging%read_ptr = new_read_ptr
      RETURN new_read_ptr + 1
      END ASSOCIATE
   END FUNCTION hfx_libheric_quartet_staging_reserve
   
   INTEGER FUNCTION hfx_libheric_quartet_staging_reserve(staging, A, B, C, D)
      IF(hfx_libheric_quartet_staging_full(staging))RETURN -1
      staging%write_ptr = MODULO(staging%write_ptr + 1, staging%n_staging)

      ASSOCIATE(curr_data => staging%data(staging%write_ptr))
      IF(curr_data%n_shellset_quartets > 0) &
         curr_data%integrals(1:curr_data%n_shellset_quartets * staging%eris_allocated_per_shellset) = 0.0_dp
      
      curr_data%quartet_id(:) = [A, B, C, D]
      curr_data%n_shellset_quartets = 0
      curr_data%prims_requested = 0
      curr_data%prims_computed = 0
      END ASSOCIATE

      RETURN staging%write_ptr

   END FUNCTION hfx_libheric_quartet_staging_reserve



   SUBROUTINE hfx_libheric_quartet_staging_free(staging)

      INTEGER :: i
      DO i = 1, staging%n_staging
         IF(ALLOCATED(staging%data(i)%integrals))DEALLOCATE(staging%data(i)%integrals)
         IF(ALLOCATED(staging%data(i)%shellset_ids))DEALLOCATE(staging%data(i)%shellset_ids)
         IF(ALLOCATED(staging%data(i)%pmax_set_vals))DEALLOCATE(staging%data(i)%pmax_set_vals)
      END DO
      IF(ALLOCATED(staging%data))DEALLOCATE(staging%data(n_staging))
   END SUBROUTINE hfx_libheric_quartet_staging_free

   LOGICAL FUNCTION hfx_libheric_initialize(mpi_rank) RESULT(initialized)
      INTEGER, INTENT(IN) :: mpi_rank
      initialized = .TRUE.
   END FUNCTION hfx_libheric_initialize

   SUBROUTINE hfx_libheric_release()
   END SUBROUTINE hfx_libheric_release
   
   SUBROUTINE hfx_libheric_query_accelerator_devices(kinds, devices)
       TYPE(hfx_libheric_acc_kind), DIMENSION(:), INTENT(IN)                 :: kinds
       TYPE(hfx_libheric_acc_device), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: devices

       ALLOCATE(devices(1))
       devices(1)%num_cores = omp_get_max_threads()

   END SUBROUTINE hfx_libheric_query_accelerator_devices

   SUBROUTINE hfx_libheric_free_accelerator_contexts(ctx_handles, num_contexts) 
      TYPE(hfx_libheric_accelerator_context), &
         DIMENSION(:), INTENT(INOUT)            :: ctx_handles
      INTEGER, INTENT(IN)                       :: num_contexts

      INTEGER :: i, j
      DO i = 1, SIZE(ctx_handles)
         ctx_handles(i)%basis_parameter_l => NULL()
         DEALLOCATE(ctx_handles(i)%atom_positions) 
         DEALLOCATE(ctx_handles(i)%atom_kinds)
      END DO
   END SUBROUTINE hfx_libheric_free_accelerator_contexts

   LOGICAL FUNCTION hfx_libheric_is_functional(ctx) RESULT(ctx_functional) 
      TYPE(hfx_libheric_accelerator_context) :: ctx
      ctx_functional = ctx%functional
   END FUNCTION hfx_libheric_is_functional

   SUBROUTINE hfx_libheric_set_potential_operator(ctx, op_id, op_params)
      TYPE(hfx_libheric_accelerator_context) :: ctx
      INTEGER                                :: op_id
      REAL(dp), DIMENSION(:)                 :: op_params

      ctx%op_id = op_id
      ctx%omega = -1.0_dp
      ctx%cutoff_radius = -1.0_dp
      SELECT CASE (op_id)
      CASE (do_potential_truncated, do_potential_TShPSC)
         ctx%cutoff_radius = op_params(1)
      CASE (do_potential_short, do_potential_long, do_potential_gaussian, do_potential_mix_cl, do_potential_mix_lg)
         ctx%omega = op_params(1)
      CASE (do_potential_mix_cl_trunc)
         ctx%cutoff_radius = op_params(1)
         ctx%omega = op_params(2)
      CASE DEFAULT
         CPABORT("Unsupported operator for radial r12 fundamentals")
      END SELECT

   END SUBROUTINE hfx_libheric_set_potential_operator

   SUBROUTINE hfx_libheric_set_quartet_class(ctx, la, lb, lc, ld) 
      TYPE(hfx_libheric_accelerator_context) :: ctx
      INTEGER                                :: la
      INTEGER                                :: lb
      INTEGER                                :: lc
      INTEGER                                :: ld

      INTEGER                                :: max_shellset_quartets

      ctx%functional = .FALSE.

      IF(la < 0 .OR. la > ctx%l_max .OR. &
         lb < 0 .OR. lb > ctx%l_max .OR. &
         lc < 0 .OR. lc > ctx%l_max .OR. &
         ld < 0 .OR. ld > ctx%l_max) CPABORT("Requested quartet class out of bounds")
      IF(ctx%quartet_class_finalizing .AND. .NOT. hfx_libheric_staging_empty(ctx%staging)) &
         CPABORT("Requested new quartet class, but staging is not empty")
      ctx%quartet_class = [la, lb, lc, ld]
      ctx%total_am = la + lb + lc + ld
      max_shellset_quartets = MAXVAL(basis_parameter(:, la)%nshellsets) & 
                            * MAXVAL(basis_parameter(:, lb)%nshellsets) & 
                            * MAXVAL(basis_parameter(:, lc)%nshellsets) & 
                            * MAXVAL(basis_parameter(:, ld)%nshellsets) 
      ctx%cart_size = nco(la) * nco(lb) * nco(lc) * nco(ld)
      ctx%sph_size  = nso(la) * nso(lb) * nso(lc) * nso(ld)
      
      CALL hfx_libheric_staging_reset(ctx%staging, ctx%quartet_class, max_shellset_quartets)

      !! Libint/cpu_acc set_quartet_class missing!

      ctx%functional = .TRUE.
      ctx%quartet_class_finalizing = .FALSE.


   END SUBROUTINE hfx_libheric_set_quartet_class

   SUBROUTINE hfx_libheric_set_atom_ids(ctx, atom_a, atom_b, atom_c, atom_d)
      TYPE(hfx_libheric_accelerator_context) :: ctx
      INTEGER                                :: atom_a, atom_b, atom_c, atom_d

      ctx%A(:) = ctx%atom_positions(:, atom_a + 1)
      ctx%B(:) = ctx%atom_positions(:, atom_b + 1)
      ctx%C(:) = ctx%atom_positions(:, atom_c + 1)
      ctx%D(:) = ctx%atom_positions(:, atom_d + 1)

      ctx%bas_a => ctx%basis_parameter_l(ctx%atom_kinds(atom_a + 1), ctx%la)
      ctx%bas_b => ctx%basis_parameter_l(ctx%atom_kinds(atom_b + 1), ctx%lb)
      ctx%bas_c => ctx%basis_parameter_l(ctx%atom_kinds(atom_c + 1), ctx%lc)
      ctx%bas_d => ctx%basis_parameter_l(ctx%atom_kinds(atom_d + 1), ctx%ld)

   END SUBROUTINE hfx_libheric_set_atom_ids

   SUBROUTINE hfx_libheric_set_B_offset(ctx, off_x, off_y, off_z) 
      TYPE(hfx_libheric_accelerator_context) :: ctx
      INTEGER                                :: off_x
      INTEGER                                :: off_y
      INTEGER                                :: off_z
   END SUBROUTINE hfx_libheric_set_B_offset

   SUBROUTINE hfx_libheric_set_D_offset(ctx, off_x, off_y, off_z) 
      TYPE(hfx_libheric_accelerator_context) :: ctx
      INTEGER                                :: off_x
      INTEGER                                :: off_y
      INTEGER                                :: off_z
   END SUBROUTINE hfx_libheric_set_D_offset

   SUBROUTINE hfx_libheric_set_Q_offset(ctx, off_x, off_y, off_z) 
      TYPE(hfx_libheric_accelerator_context) :: ctx
      INTEGER                                :: off_x
      INTEGER                                :: off_y
      INTEGER                                :: off_z
   END SUBROUTINE hfx_libheric_set_Q_offset

   SUBROUTINE hfx_libheric_set_set_ids(ctx, set_a, set_b, set_c, set_d, pmax)
      TYPE(hfx_libheric_accelerator_context) :: ctx
      INTEGER                                :: set_a, set_b, set_c, set_d
      REAL(dp)                               :: pmax 
   END SUBROUTINE hfx_libheric_set_set_ids

   SUBROUTINE hfx_libheric_set_shell_ids(ctx, shl_a, shl_b, shl_c, shl_d)
      TYPE(hfx_libheric_accelerator_context) :: ctx
      INTEGER                                :: shl_a, shl_b, shl_c, shl_d
   END SUBROUTINE hfx_libheric_set_shell_ids

   SUBROUTINE hfx_libheric_set_exponent_ids(ctx, exp_a, exp_b, exp_c, exp_d)
      TYPE(hfx_libheric_accelerator_context) :: ctx
      INTEGER                                :: exp_a, exp_b, exp_c, exp_d
   END SUBROUTINE hfx_libheric_set_exponent_ids

   SUBROUTINE hfx_libheric_submit_peri(ctx) 
      TYPE(hfx_libheric_accelerator_context) :: ctx
   END SUBROUTINE hfx_libheric_submit_peri

   SUBROUTINE hfx_libheric_finish_quartet_class(ctx) 
      TYPE(hfx_libheric_accelerator_context) :: ctx
   END SUBROUTINE hfx_libheric_finish_quartet_class

   LOGICAL FUNCTION hfx_libheric_prefetch_cached_atom_quartet(ctx, A, B, C, D) RESULT(quartet_available)
      TYPE(hfx_libheric_accelerator_context) :: ctx
      INTEGER                                :: A, B, C, D
      quartet_available = .FALSE.
   END FUNCTION hfx_libheric_prefetch_cached_atom_quartet
    
   LOGICAL FUNCTION hfx_libheric_get_next_staged_atom_quartet(ctx, atom_quartet) RESULT(quartet_available)
       TYPE(hfx_libheric_accelerator_context), INTENT(INOUT) :: ctx
       TYPE(hfx_libheric_atom_quartet), INTENT(INOUT)        :: atom_quartet
       quartet_available = .FALSE.
   END FUNCTION hfx_libheric_get_next_staged_atom_quartet 

   SUBROUTINE hfx_libheric_update_positions(ctx, raw_positions, num_positions) 
      TYPE(hfx_libheric_accelerator_context) :: ctx
      REAL(dp)                               :: raw_positions
      INTEGER                                :: num_positions
   END SUBROUTINE hfx_libheric_update_positions

   SUBROUTINE hfx_libheric_update_cell(ctx, raw_cell_mat) 
      TYPE(hfx_libheric_accelerator_context)    :: ctx
      REAL(dp)        :: raw_cell_mat
   END SUBROUTINE hfx_libheric_update_cell

   SUBROUTINE hfx_libheric_optimize_screening_function(x_data, y_data, n_data, y_cutoff, model, n_model)
      REAL(dp) :: x_data
      REAL(dp) :: y_data
      INTEGER  :: n_data
      REAL(dp) :: y_cutoff
      REAL(dp) :: model
      INTEGER  :: n_model
   END SUBROUTINE
   
#endif
END MODULE hfx_libheric_wrapper
