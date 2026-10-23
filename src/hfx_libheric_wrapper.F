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
#endif
   USE particle_types, ONLY: particle_type
   USE kinds, ONLY: dp, &
                    int_8
    
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
      TYPE(particle_type), DIMENSION(:), POINTER              :: particle_set => NULL()
      REAL(dp), DIMENSION(:, :), POINTER                      :: cell_mat => NULL()
   END TYPE hfx_libheric_accelerator_context

   TYPE hfx_libheric_acc_device
      CHARACTER(kind=C_CHAR, len=LIBHERIC_FSTR_LEN) :: name
      CHARACTER(kind=C_CHAR, len=LIBHERIC_FSTR_LEN) :: bdf
      CHARACTER(kind=C_CHAR, len=LIBHERIC_FSTR_LEN) :: cpu_name
      INTEGER                                       :: num_cores
      TYPE(acc_kind)                                :: kind
   END TYPE hfx_libheric_acc_device

#endif
CONTAINS
   SUBROUTINE hfx_libheric_create_accelerator_contexts(max_set, nkind, l_max, basis_parameter_l, &
                                                       cell_mat, acc_dev, particle_set, & 
                                                       mem_available, compute_raw_cartesian, ctx_handles) 
      INTEGER, INTENT(IN)                                 :: max_set, nkind, l_max
      TYPE(hfx_libheric_basis_type), &
         DIMENSION(nkind, 0:l_max), TARGET, INTENT(IN)    :: basis_parameter_l
      REAL(dp), DIMENSION(3, 3), TARGET, INTENT(IN)       :: cell_mat
      TYPE(hfx_libheric_acc_device), INTENT(IN)           :: acc_dev
      TYPE(particle_type), DIMENSION(:), &
         TARGET, INTENT(IN)                               :: particle_set
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
      INTEGER :: i
      ALLOCATE(ctx_handles(acc_dev%num_cores))
      DO i = 1, SIZE(ctx_handles)
         ctx_handles(i)%basis_parameter_l => basis_parameter_l
         ctx_handles(i)%particle_set => particle_set
         ctx_handles(i)%cell_mat => cell_mat
      END DO
#endif
   END SUBROUTINE hfx_libheric_create_accelerator_contexts
#if !defined(__LIBHERIC)

   SUBROUTINE hfx_set_quartet_class(ctx, la, lb, lc, ld)
      TYPE(hfx_libheric_accelerator_context) :: ctx
      INTEGER                                :: la, lb, lc, ld
   END SUBROUTINE
   
#endif
END MODULE hfx_libheric_wrapper
