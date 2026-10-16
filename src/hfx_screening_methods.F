!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Several screening methods used in HFX calcualtions
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
MODULE hfx_screening_methods
   USE ao_util,                         ONLY: exp_radius_very_extended
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE hfx_libint_interface,            ONLY: evaluate_eri_screen
   USE hfx_types,                       ONLY: hfx_type, &
                                              hfx_basis_type,&
                                              hfx_p_kind,&
                                              hfx_potential_type,&
                                              hfx_screen_coeff_type,&
                                              log_zero,&
                                              powell_min_log, &
                                              hfx_libheric_basis_type,&
                                              hfx_libheric_screening_model_type,&
                                              hfx_libheric_screening_models, &
                                              hfx_libheric_create_accelerator_contexts
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE libint_wrapper,                  ONLY: cp_libint_t
   USE libheric_interface,              ONLY: libheric_optimize_screening_function => optimize_screening_function, &
                                              query_accelerator_devices, &
                                              free_accelerator_contexts, &
                                              accelerator_context_is_functional, &
                                              accelerator_context_set_potential_operator, &
                                              accelerator_context_set_quartet_class, &
                                              accelerator_context_set_atom_ids, &
                                              accelerator_context_set_set_ids, &
                                              accelerator_context_set_shell_ids, &
                                              accelerator_context_set_B_offset, &
                                              accelerator_context_set_D_offset, &
                                              accelerator_context_set_exponent_ids, &
                                              accelerator_context_set_Q_offset, &
                                              accelerator_context_submit_peri, &
                                              accelerator_context_finish_quartet_class, &
                                              accelerator_context_get_next_staged_atom_quartet, &
                                              accelerator_context_update_cell, &
                                              flush_logger, &
                                              atom_quartet_type, &
                                              acc_kind, acc_device
   USE machine,                         ONLY: default_output_unit
   USE orbital_pointers,                ONLY: ncoset, nco, nso
   USE particle_types,                  ONLY: particle_type
   USE powell,                          ONLY: opt_state_type,&
                                              powell_optimize
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE, INTRINSIC :: iso_c_binding,     ONLY: C_PTR, C_BOOL
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: update_pmax_mat, &
             calc_screening_functions, &
             calc_pair_dist_radii, &
             hfx_libheric_init_screening_models, &
             hfx_libheric_free_screening_models, &
             hfx_libheric_eval_screening_model

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_screening_methods'

!***

CONTAINS

! **************************************************************************************************
!> \brief calculates max values of two-electron integrals in a quartet/shell
!>      w.r.t. different zetas using the library lib_int
!> \param lib ...
!> \param ra position
!> \param rb position
!> \param zeta zeta
!> \param zetb zeta
!> \param la_min angular momentum
!> \param la_max angular momentum
!> \param lb_min angular momentum
!> \param lb_max angular momentum
!> \param npgfa number of primitive cartesian gaussian in actual shell
!> \param npgfb number of primitive cartesian gaussian in actual shell
!> \param max_val schwarz screening value
!> \param potential_parameter contains info for libint
!> \param tmp_R_1 pgf_based screening coefficients
!> \param rab2 Squared Distance of centers ab
!> \param p_work ...
!> \par History
!>     03.2007 created [Manuel Guidon]
!>     02.2009 refactored [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE screen4(lib, ra, rb, zeta, zetb, &
                      la_min, la_max, lb_min, lb_max, &
                      npgfa, npgfb, &
                      max_val, potential_parameter, tmp_R_1, &
                      rab2, p_work)

      TYPE(cp_libint_t)                                  :: lib
      REAL(dp), INTENT(IN)                               :: ra(3), rb(3)
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zeta, zetb
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max, npgfa, &
                                                            npgfb
      REAL(dp), INTENT(INOUT)                            :: max_val
      TYPE(hfx_potential_type)                           :: potential_parameter
      TYPE(hfx_screen_coeff_type), DIMENSION(:, :), &
         POINTER                                         :: tmp_R_1
      REAL(dp)                                           :: rab2
      REAL(dp), DIMENSION(:), POINTER                    :: p_work

      INTEGER                                            :: ipgf, jpgf, la, lb
      REAL(dp)                                           :: max_val_temp, R1

      max_val_temp = max_val
      DO ipgf = 1, npgfa
         DO jpgf = 1, npgfb
            R1 = MAX(0.0_dp, tmp_R_1(jpgf, ipgf)%x(1)*rab2 + tmp_R_1(jpgf, ipgf)%x(2))
            DO la = la_min, la_max
               DO lb = lb_min, lb_max
                  !Build primitives
                  CALL evaluate_eri_screen(lib, ra, rb, ra, rb, &
                                           zeta(ipgf), zetb(jpgf), zeta(ipgf), zetb(jpgf), &
                                           la, lb, la, lb, &
                                           max_val_temp, potential_parameter, R1, R1, p_work)

                  max_val = MAX(max_val, max_val_temp)
               END DO !lb
            END DO !la
         END DO !jpgf
      END DO !ipgf
   END SUBROUTINE screen4

! **************************************************************************************************
!> \brief updates the maximum of the density matrix in compressed form for screening purposes
!> \param pmax_set buffer to store matrix
!> \param map_atom_to_kind_atom ...
!> \param set_offset ...
!> \param atomic_block_offset ...
!> \param pmax_atom ...
!> \param full_density_alpha ...
!> \param full_density_beta ...
!> \param natom ...
!> \param kind_of ...
!> \param basis_parameter ...
!> \param nkind ...
!> \param is_assoc_atomic_block ...
!> \par History
!>     09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      - updates for each pair of shells the maximum absolute value of p
! **************************************************************************************************
   SUBROUTINE update_pmax_mat(pmax_set, map_atom_to_kind_atom, set_offset, atomic_block_offset, &
                              pmax_atom, full_density_alpha, full_density_beta, natom, &
                              kind_of, basis_parameter, &
                              nkind, is_assoc_atomic_block)

      TYPE(hfx_p_kind), DIMENSION(:), POINTER            :: pmax_set
      INTEGER, DIMENSION(:), POINTER                     :: map_atom_to_kind_atom
      INTEGER, DIMENSION(:, :, :, :), POINTER            :: set_offset
      INTEGER, DIMENSION(:, :), POINTER                  :: atomic_block_offset
      REAL(dp), DIMENSION(:, :), POINTER                 :: pmax_atom, full_density_alpha, &
                                                            full_density_beta
      INTEGER, INTENT(IN)                                :: natom
      INTEGER                                            :: kind_of(*)
      TYPE(hfx_basis_type), DIMENSION(:), POINTER        :: basis_parameter
      INTEGER                                            :: nkind
      INTEGER, DIMENSION(:, :)                           :: is_assoc_atomic_block

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'update_pmax_mat'

      INTEGER :: act_atomic_block_offset, act_set_offset, handle, i, img, jatom, jkind, jset, &
         katom, kind_kind_idx, kkind, kset, mb, mc, nimg, nsetb, nsetc
      INTEGER, DIMENSION(:), POINTER                     :: nsgfb, nsgfc
      REAL(dp)                                           :: pmax_tmp

      CALL timeset(routineN, handle)

      DO i = 1, SIZE(pmax_set)
         pmax_set(i)%p_kind = 0.0_dp
      END DO

      pmax_atom = log_zero

      nimg = SIZE(full_density_alpha, 2)

      DO jatom = 1, natom
         jkind = kind_of(jatom)
         nsetb = basis_parameter(jkind)%nset
         nsgfb => basis_parameter(jkind)%nsgf

         DO katom = 1, natom
            IF (.NOT. is_assoc_atomic_block(jatom, katom) >= 1) CYCLE
            kkind = kind_of(katom)
            IF (kkind < jkind) CYCLE
            nsetc = basis_parameter(kkind)%nset
            nsgfc => basis_parameter(kkind)%nsgf
            act_atomic_block_offset = atomic_block_offset(katom, jatom)
            DO jset = 1, nsetb
               DO kset = 1, nsetc
                  IF (katom >= jatom) THEN
                     pmax_tmp = 0.0_dp
                     act_set_offset = set_offset(kset, jset, kkind, jkind)
                     DO img = 1, nimg
                        i = act_set_offset + act_atomic_block_offset - 1
                        DO mc = 1, nsgfc(kset)
                           DO mb = 1, nsgfb(jset)
                              pmax_tmp = MAX(pmax_tmp, ABS(full_density_alpha(i, img)))
                              IF (ASSOCIATED(full_density_beta)) THEN
                                 pmax_tmp = MAX(pmax_tmp, ABS(full_density_beta(i, img)))
                              END IF
                              i = i + 1
                           END DO
                        END DO
                     END DO
                     IF (pmax_tmp == 0.0_dp) THEN
                        pmax_tmp = log_zero
                     ELSE
                        pmax_tmp = LOG10(pmax_tmp)
                     END IF
                     kind_kind_idx = INT(get_1D_idx(jkind, kkind, INT(nkind, int_8)))
                     pmax_set(kind_kind_idx)%p_kind(jset, &
                                                    kset, &
                                                    map_atom_to_kind_atom(jatom), &
                                                    map_atom_to_kind_atom(katom)) = pmax_tmp
                  ELSE
                     pmax_tmp = 0.0_dp
                     act_set_offset = set_offset(jset, kset, jkind, kkind)
                     DO img = 1, nimg
                        DO mc = 1, nsgfc(kset)
                           i = act_set_offset + act_atomic_block_offset - 1 + mc - 1
                           DO mb = 1, nsgfb(jset)
                              pmax_tmp = MAX(pmax_tmp, ABS(full_density_alpha(i, img)))
                              IF (ASSOCIATED(full_density_beta)) THEN
                                 pmax_tmp = MAX(pmax_tmp, ABS(full_density_beta(i, img)))
                              END IF
                              i = i + nsgfc(kset)
                           END DO
                        END DO
                     END DO
                     IF (pmax_tmp == 0.0_dp) THEN
                        pmax_tmp = log_zero
                     ELSE
                        pmax_tmp = LOG10(pmax_tmp)
                     END IF
                     kind_kind_idx = INT(get_1D_idx(jkind, kkind, INT(nkind, int_8)))
                     pmax_set(kind_kind_idx)%p_kind(jset, &
                                                    kset, &
                                                    map_atom_to_kind_atom(jatom), &
                                                    map_atom_to_kind_atom(katom)) = pmax_tmp
                  END IF
               END DO
            END DO
         END DO
      END DO

      DO jatom = 1, natom
         jkind = kind_of(jatom)
         nsetb = basis_parameter(jkind)%nset
         DO katom = 1, natom
            IF (.NOT. is_assoc_atomic_block(jatom, katom) >= 1) CYCLE
            kkind = kind_of(katom)
            IF (kkind < jkind) CYCLE
            nsetc = basis_parameter(kkind)%nset
            pmax_tmp = log_zero
            DO jset = 1, nsetb
               DO kset = 1, nsetc
                  kind_kind_idx = INT(get_1D_idx(jkind, kkind, INT(nkind, int_8)))
                  pmax_tmp = MAX(pmax_set(kind_kind_idx)%p_kind(jset, &
                                                                kset, &
                                                                map_atom_to_kind_atom(jatom), &
                                                                map_atom_to_kind_atom(katom)), pmax_tmp)
               END DO
            END DO
            pmax_atom(jatom, katom) = pmax_tmp
            pmax_atom(katom, jatom) = pmax_tmp
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE update_pmax_mat

! **************************************************************************************************
!> \brief calculates screening functions for schwarz screening
!> \param qs_env qs_env
!> \param basis_parameter ...
!> \param lib structure to libint
!> \param potential_parameter contains infos on potential
!> \param coeffs_set set based coefficients
!> \param coeffs_kind kind based coefficients
!> \param coeffs_pgf pgf based coefficients
!> \param radii_pgf coefficients for long-range screening
!> \param max_set Maximum Number of basis set sets in the system
!> \param max_pgf Maximum Number of basis set pgfs in the system
!> \param n_threads ...
!> \param i_thread Thread ID of current task
!> \param p_work ...
!> \par History
!>     02.2009 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      This routine calculates (ab|ab) for different distances Rab = |a-b|
!>      and uses the powell optimiztion routines in order to fit the results
!>      in the following form:
!>
!>                 (ab|ab) = (ab|ab)(Rab) = c2*Rab^2 + c0
!>
!>      The missing linear term assures that the functions is monotonically
!>      decaying such that c2 can be used as upper bound when applying the
!>      Minimum Image Convention in the periodic case. Furthermore
!>      it seems to be a good choice to fit the logarithm of the (ab|ab)
!>      The fitting takes place at several levels: kind, set and pgf within
!>      the corresponding ranges of the prodiuct charge distributions
!>      Doing so, we only need arrays of size nkinds^2*2 instead of big
!>      screening matrices
! **************************************************************************************************

   SUBROUTINE calc_screening_functions(qs_env, basis_parameter, lib, potential_parameter, &
                                       coeffs_set, coeffs_kind, coeffs_pgf, radii_pgf, &
                                       max_set, max_pgf, n_threads, i_thread, &
                                       p_work)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_basis_type), DIMENSION(:), POINTER        :: basis_parameter
      TYPE(cp_libint_t)                                  :: lib
      TYPE(hfx_potential_type)                           :: potential_parameter
      TYPE(hfx_screen_coeff_type), &
         DIMENSION(:, :, :, :), POINTER                  :: coeffs_set
      TYPE(hfx_screen_coeff_type), DIMENSION(:, :), &
         POINTER                                         :: coeffs_kind
      TYPE(hfx_screen_coeff_type), &
         DIMENSION(:, :, :, :, :, :), POINTER            :: coeffs_pgf, radii_pgf
      INTEGER, INTENT(IN)                                :: max_set, max_pgf, n_threads, i_thread
      REAL(dp), DIMENSION(:), POINTER                    :: p_work

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_screening_functions'

      INTEGER                                            :: handle, i, ikind, ipgf, iset, jkind, &
                                                            jpgf, jset, la, lb, ncoa, ncob, nkind, &
                                                            nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      REAL(dp) :: kind_radius_a, kind_radius_b, max_contraction_a, max_contraction_b, max_val, &
         max_val_temp, R1, ra(3), radius, rb(3), x(2)
      REAL(dp), DIMENSION(:), POINTER                    :: set_radius_a, set_radius_b
      REAL(dp), DIMENSION(:, :), POINTER                 :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      REAL(dp), SAVE                                     :: DATA(2, 0:100)
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis
      TYPE(hfx_screen_coeff_type), DIMENSION(:, :), &
         POINTER                                         :: tmp_R_1
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

!$OMP MASTER
      CALL timeset(routineN, handle)
!$OMP END MASTER

      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set)

      nkind = SIZE(qs_kind_set, 1)

!$OMP MASTER
      ALLOCATE (coeffs_pgf(max_pgf, max_pgf, max_set, max_set, nkind, nkind))

      DO ikind = 1, nkind
         DO jkind = 1, nkind
            DO iset = 1, max_set
               DO jset = 1, max_set
                  DO ipgf = 1, max_pgf
                     DO jpgf = 1, max_pgf
                        coeffs_pgf(jpgf, ipgf, jset, iset, jkind, ikind)%x(:) = 0.0_dp
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

!$OMP END MASTER
!$OMP BARRIER
      ra = 0.0_dp
      rb = 0.0_dp
      DO ikind = 1, nkind
         NULLIFY (qs_kind, orb_basis)
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=orb_basis)
         NULLIFY (la_max, la_min, npgfa, zeta)

         la_max => basis_parameter(ikind)%lmax
         la_min => basis_parameter(ikind)%lmin
         npgfa => basis_parameter(ikind)%npgf
         nseta = basis_parameter(ikind)%nset
         zeta => basis_parameter(ikind)%zet
         set_radius_a => basis_parameter(ikind)%set_radius
         first_sgfa => basis_parameter(ikind)%first_sgf
         sphi_a => basis_parameter(ikind)%sphi
         nsgfa => basis_parameter(ikind)%nsgf
         rpgfa => basis_parameter(ikind)%pgf_radius

         DO jkind = 1, nkind
            NULLIFY (qs_kind, orb_basis)
            qs_kind => qs_kind_set(jkind)
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=orb_basis)
            NULLIFY (lb_max, lb_min, npgfb, zetb)

            lb_max => basis_parameter(jkind)%lmax
            lb_min => basis_parameter(jkind)%lmin
            npgfb => basis_parameter(jkind)%npgf
            nsetb = basis_parameter(jkind)%nset
            zetb => basis_parameter(jkind)%zet
            set_radius_b => basis_parameter(jkind)%set_radius
            first_sgfb => basis_parameter(jkind)%first_sgf
            sphi_b => basis_parameter(jkind)%sphi
            nsgfb => basis_parameter(jkind)%nsgf
            rpgfb => basis_parameter(jkind)%pgf_radius

            DO iset = 1, nseta
               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               max_contraction_a = MAXVAL([(SUM(ABS(sphi_a(1:ncoa, i))), i=sgfa, sgfa + nsgfa(iset) - 1)])
               DO jset = 1, nsetb
                  ncob = npgfb(jset)*ncoset(lb_max(jset))
                  sgfb = first_sgfb(1, jset)
                  max_contraction_b = MAXVAL([(SUM(ABS(sphi_b(1:ncob, i))), i=sgfb, sgfb + nsgfb(jset) - 1)])
                  radius = set_radius_a(iset) + set_radius_b(jset)
                  DO ipgf = 1, npgfa(iset)
                     DO jpgf = 1, npgfb(jset)
                        radius = rpgfa(ipgf, iset) + rpgfb(jpgf, jset)
                        DO i = i_thread, 100, n_threads
                           rb(1) = 0.0_dp + REAL(i, dp)*0.01_dp*radius
                           max_val = 0.0_dp
                           R1 = MAX(0.0_dp, radii_pgf(jpgf, ipgf, jset, iset, jkind, ikind)%x(1)*rb(1)**2 + &
                                    radii_pgf(jpgf, ipgf, jset, iset, jkind, ikind)%x(2))
                           DO la = la_min(iset), la_max(iset)
                              DO lb = lb_min(jset), lb_max(jset)
                                 !Build primitives
                                 max_val_temp = 0.0_dp
                                 CALL evaluate_eri_screen(lib, ra, rb, ra, rb, &
                                                          zeta(ipgf, iset), zetb(jpgf, jset), zeta(ipgf, iset), zetb(jpgf, jset), &
                                                          la, lb, la, lb, &
                                                          max_val_temp, potential_parameter, R1, R1, p_work)
                                 max_val = MAX(max_val, max_val_temp)
                              END DO !lb
                           END DO !la
                           max_val = SQRT(max_val)
                           max_val = max_val*max_contraction_a*max_contraction_b
                           DATA(1, i) = rb(1)
                           IF (max_val == 0.0_dp) THEN
                              DATA(2, i) = powell_min_log
                           ELSE
                              DATA(2, i) = LOG10((max_val))
                           END IF
                        END DO
!$OMP BARRIER
!$OMP MASTER
                        CALL optimize_it(DATA, x, powell_min_log)
                        coeffs_pgf(jpgf, ipgf, jset, iset, jkind, ikind)%x = x
!$OMP END MASTER
!$OMP BARRIER

                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

!$OMP MASTER
      ALLOCATE (coeffs_set(max_set, max_set, nkind, nkind))

      DO ikind = 1, nkind
         DO jkind = 1, nkind
            DO iset = 1, max_set
               DO jset = 1, max_set
                  coeffs_set(jset, iset, jkind, ikind)%x(:) = 0.0_dp
               END DO
            END DO
         END DO
      END DO
!$OMP END MASTER
!$OMP BARRIER
      ra = 0.0_dp
      rb = 0.0_dp
      DO ikind = 1, nkind
         NULLIFY (qs_kind, orb_basis)
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=orb_basis)
         NULLIFY (la_max, la_min, npgfa, zeta)
!      CALL get_gto_basis_set(gto_basis_set=orb_basis,&
!                             lmax=la_max,&
!                             lmin=la_min,&
!                             npgf=npgfa,&
!                             nset=nseta,&
!                             zet=zeta,&
!                             set_radius=set_radius_a,&
!                             first_sgf=first_sgfa,&
!                             sphi=sphi_a,&
!                             nsgf_set=nsgfa)
         la_max => basis_parameter(ikind)%lmax
         la_min => basis_parameter(ikind)%lmin
         npgfa => basis_parameter(ikind)%npgf
         nseta = basis_parameter(ikind)%nset
         zeta => basis_parameter(ikind)%zet
         set_radius_a => basis_parameter(ikind)%set_radius
         first_sgfa => basis_parameter(ikind)%first_sgf
         sphi_a => basis_parameter(ikind)%sphi
         nsgfa => basis_parameter(ikind)%nsgf

         DO jkind = 1, nkind
            NULLIFY (qs_kind, orb_basis)
            qs_kind => qs_kind_set(jkind)
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=orb_basis)
            NULLIFY (lb_max, lb_min, npgfb, zetb)

            lb_max => basis_parameter(jkind)%lmax
            lb_min => basis_parameter(jkind)%lmin
            npgfb => basis_parameter(jkind)%npgf
            nsetb = basis_parameter(jkind)%nset
            zetb => basis_parameter(jkind)%zet
            set_radius_b => basis_parameter(jkind)%set_radius
            first_sgfb => basis_parameter(jkind)%first_sgf
            sphi_b => basis_parameter(jkind)%sphi
            nsgfb => basis_parameter(jkind)%nsgf

            DO iset = 1, nseta
               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               max_contraction_a = MAXVAL([(SUM(ABS(sphi_a(1:ncoa, i))), i=sgfa, sgfa + nsgfa(iset) - 1)])
               DO jset = 1, nsetb
                  ncob = npgfb(jset)*ncoset(lb_max(jset))
                  sgfb = first_sgfb(1, jset)
                  max_contraction_b = MAXVAL([(SUM(ABS(sphi_b(1:ncob, i))), i=sgfb, sgfb + nsgfb(jset) - 1)])
                  radius = set_radius_a(iset) + set_radius_b(jset)
                  tmp_R_1 => radii_pgf(:, :, jset, iset, jkind, ikind)
                  DO i = i_thread, 100, n_threads
                     rb(1) = 0.0_dp + REAL(i, dp)*0.01_dp*radius
                     max_val = 0.0_dp
                     CALL screen4(lib, ra, rb, &
                                  zeta(:, iset), zetb(:, jset), &
                                  la_min(iset), la_max(iset), lb_min(jset), lb_max(jset), &
                                  npgfa(iset), npgfb(jset), &
                                  max_val, potential_parameter, tmp_R_1, rb(1)**2, p_work)
                     max_val = SQRT(max_val)
                     max_val = max_val*max_contraction_a*max_contraction_b
                     DATA(1, i) = rb(1)
                     IF (max_val == 0.0_dp) THEN
                        DATA(2, i) = powell_min_log
                     ELSE
                        DATA(2, i) = LOG10((max_val))
                     END IF
                  END DO
!$OMP BARRIER
!$OMP MASTER
                  CALL optimize_it(DATA, x, powell_min_log)
                  coeffs_set(jset, iset, jkind, ikind)%x = x
!$OMP END MASTER
!$OMP BARRIER
               END DO
            END DO

         END DO
      END DO

      ! ** now kinds
!$OMP MASTER
      ALLOCATE (coeffs_kind(nkind, nkind))

      DO ikind = 1, nkind
         DO jkind = 1, nkind
            coeffs_kind(jkind, ikind)%x(:) = 0.0_dp
         END DO
      END DO
!$OMP END MASTER
      ra = 0.0_dp
      rb = 0.0_dp
      DO ikind = 1, nkind
         NULLIFY (qs_kind, orb_basis)
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=orb_basis)
         NULLIFY (la_max, la_min, npgfa, zeta)

         la_max => basis_parameter(ikind)%lmax
         la_min => basis_parameter(ikind)%lmin
         npgfa => basis_parameter(ikind)%npgf
         nseta = basis_parameter(ikind)%nset
         zeta => basis_parameter(ikind)%zet
         kind_radius_a = basis_parameter(ikind)%kind_radius
         first_sgfa => basis_parameter(ikind)%first_sgf
         sphi_a => basis_parameter(ikind)%sphi
         nsgfa => basis_parameter(ikind)%nsgf

         DO jkind = 1, nkind
            NULLIFY (qs_kind, orb_basis)
            qs_kind => qs_kind_set(jkind)
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=orb_basis)
            NULLIFY (lb_max, lb_min, npgfb, zetb)

            lb_max => basis_parameter(jkind)%lmax
            lb_min => basis_parameter(jkind)%lmin
            npgfb => basis_parameter(jkind)%npgf
            nsetb = basis_parameter(jkind)%nset
            zetb => basis_parameter(jkind)%zet
            kind_radius_b = basis_parameter(jkind)%kind_radius
            first_sgfb => basis_parameter(jkind)%first_sgf
            sphi_b => basis_parameter(jkind)%sphi
            nsgfb => basis_parameter(jkind)%nsgf

            radius = kind_radius_a + kind_radius_b
            DO iset = 1, nseta
               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               max_contraction_a = MAXVAL([(SUM(ABS(sphi_a(1:ncoa, i))), i=sgfa, sgfa + nsgfa(iset) - 1)])
               DO jset = 1, nsetb
                  ncob = npgfb(jset)*ncoset(lb_max(jset))
                  sgfb = first_sgfb(1, jset)
                  max_contraction_b = MAXVAL([(SUM(ABS(sphi_b(1:ncob, i))), i=sgfb, sgfb + nsgfb(jset) - 1)])
                  DO i = i_thread, 100, n_threads
                     rb(1) = 0.0_dp + REAL(i, dp)*0.01_dp*radius
                     max_val = 0.0_dp
                     tmp_R_1 => radii_pgf(:, :, jset, iset, jkind, ikind)
                     CALL screen4(lib, ra, rb, &
                                  zeta(:, iset), zetb(:, jset), &
                                  la_min(iset), la_max(iset), lb_min(jset), lb_max(jset), &
                                  npgfa(iset), npgfb(jset), &
                                  max_val, potential_parameter, tmp_R_1, rb(1)**2, p_work)
                     DATA(1, i) = rb(1)
                     max_val = SQRT(max_val)
                     max_val = max_val*max_contraction_a*max_contraction_b
                     IF (max_val == 0.0_dp) THEN
                        DATA(2, i) = MAX(powell_min_log, DATA(2, i))
                     ELSE
                        DATA(2, i) = MAX(LOG10(max_val), DATA(2, i))
                     END IF
                  END DO
               END DO
            END DO
!$OMP BARRIER
!$OMP MASTER
            CALL optimize_it(DATA, x, powell_min_log)
            coeffs_kind(jkind, ikind)%x = x
!$OMP END MASTER
!$OMP BARRIER
         END DO
      END DO

!$OMP MASTER
      CALL timestop(handle)
!$OMP END MASTER

   END SUBROUTINE calc_screening_functions

! **************************************************************************************************
!> \brief calculates radius functions for longrange screening
!> \param qs_env qs_env
!> \param basis_parameter ...
!> \param radii_pgf pgf based coefficients
!> \param max_set Maximum Number of basis set sets in the system
!> \param max_pgf Maximum Number of basis set pgfs in the system
!> \param eps_schwarz ...
!> \param n_threads ...
!> \param i_thread Thread ID of current task
!> \par History
!>     02.2009 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      This routine calculates the pair-distribution radius of a product
!>      Gaussian and uses the powell optimiztion routines in order to fit
!>      the results in the following form:
!>
!>                 (ab| = (ab(Rab) = c2*Rab^2 + c0
!>
! **************************************************************************************************

   SUBROUTINE calc_pair_dist_radii(qs_env, basis_parameter, &
                                   radii_pgf, max_set, max_pgf, eps_schwarz, &
                                   n_threads, i_thread)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_basis_type), DIMENSION(:), POINTER        :: basis_parameter
      TYPE(hfx_screen_coeff_type), &
         DIMENSION(:, :, :, :, :, :), POINTER            :: radii_pgf
      INTEGER, INTENT(IN)                                :: max_set, max_pgf
      REAL(dp)                                           :: eps_schwarz
      INTEGER, INTENT(IN)                                :: n_threads, i_thread

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_pair_dist_radii'

      INTEGER                                            :: handle, i, ikind, ipgf, iset, jkind, &
                                                            jpgf, jset, la, lb, ncoa, ncob, nkind, &
                                                            nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      REAL(dp) :: cutoff, ff, max_contraction_a, max_contraction_b, prefactor, R1, R_max, ra(3), &
         rab(3), rab2, radius, rap(3), rb(3), rp(3), x(2), zetp
      REAL(dp), DIMENSION(:, :), POINTER                 :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      REAL(dp), SAVE                                     :: DATA(2, 0:100)
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

!$OMP MASTER
      CALL timeset(routineN, handle)
!$OMP END MASTER
      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set)

      nkind = SIZE(qs_kind_set, 1)
      ra = 0.0_dp
      rb = 0.0_dp
!$OMP MASTER
      ALLOCATE (radii_pgf(max_pgf, max_pgf, max_set, max_set, nkind, nkind))
      DO ikind = 1, nkind
         DO jkind = 1, nkind
            DO iset = 1, max_set
               DO jset = 1, max_set
                  DO ipgf = 1, max_pgf
                     DO jpgf = 1, max_pgf
                        radii_pgf(jpgf, ipgf, jset, iset, jkind, ikind)%x(:) = 0.0_dp
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      DATA = 0.0_dp
!$OMP END MASTER
!$OMP BARRIER

      DO ikind = 1, nkind
         NULLIFY (qs_kind, orb_basis)
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=orb_basis)
         NULLIFY (la_max, la_min, npgfa, zeta)

         la_max => basis_parameter(ikind)%lmax
         la_min => basis_parameter(ikind)%lmin
         npgfa => basis_parameter(ikind)%npgf
         nseta = basis_parameter(ikind)%nset
         zeta => basis_parameter(ikind)%zet
         first_sgfa => basis_parameter(ikind)%first_sgf
         sphi_a => basis_parameter(ikind)%sphi
         nsgfa => basis_parameter(ikind)%nsgf
         rpgfa => basis_parameter(ikind)%pgf_radius

         DO jkind = 1, nkind
            NULLIFY (qs_kind, orb_basis)
            qs_kind => qs_kind_set(jkind)
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=orb_basis)
            NULLIFY (lb_max, lb_min, npgfb, zetb)

            lb_max => basis_parameter(jkind)%lmax
            lb_min => basis_parameter(jkind)%lmin
            npgfb => basis_parameter(jkind)%npgf
            nsetb = basis_parameter(jkind)%nset
            zetb => basis_parameter(jkind)%zet
            first_sgfb => basis_parameter(jkind)%first_sgf
            sphi_b => basis_parameter(jkind)%sphi
            nsgfb => basis_parameter(jkind)%nsgf
            rpgfb => basis_parameter(jkind)%pgf_radius

            DO iset = 1, nseta
               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               max_contraction_a = MAXVAL([(SUM(ABS(sphi_a(1:ncoa, i))), i=sgfa, sgfa + nsgfa(iset) - 1)])
               DO jset = 1, nsetb
                  ncob = npgfb(jset)*ncoset(lb_max(jset))
                  sgfb = first_sgfb(1, jset)
                  max_contraction_b = MAXVAL([(SUM(ABS(sphi_b(1:ncob, i))), i=sgfb, sgfb + nsgfb(jset) - 1)])
                  DO ipgf = 1, npgfa(iset)
                     DO jpgf = 1, npgfb(jset)
                        radius = rpgfa(ipgf, iset) + rpgfb(jpgf, jset)
                        DO i = i_thread, 100, n_threads
                           rb(1) = 0.0_dp + 0.01_dp*radius*i
                           R_max = 0.0_dp
                           DO la = la_min(iset), la_max(iset)
                              DO lb = lb_min(jset), lb_max(jset)
                                 zetp = zeta(ipgf, iset) + zetb(jpgf, jset)
                                 ff = zetb(jpgf, jset)/zetp
                                 rab = 0.0_dp
                                 rab(1) = rb(1)
                                 rab2 = rb(1)**2
                                 prefactor = EXP(-zeta(ipgf, iset)*ff*rab2)
                                 rap(:) = ff*rab(:)
                                 rp(:) = ra(:) + rap(:)
                                 rb(:) = ra(:) + rab(:)
                                 cutoff = 1.0_dp
                                 R1 = exp_radius_very_extended( &
                                      la, la, lb, lb, ra=ra, rb=rb, rp=rp, &
                                      zetp=zetp, eps=eps_schwarz, prefactor=prefactor, cutoff=cutoff, epsabs=1.0E-12_dp)
                                 R_max = MAX(R_max, R1)
                              END DO
                           END DO
                           DATA(1, i) = rb(1)
                           DATA(2, i) = R_max
                        END DO
                        ! the radius can not be negative, we take that into account in the code as well by using a MAX
                        ! the functional form we use for fitting does not seem particularly accurate
!$OMP BARRIER
!$OMP MASTER
                        CALL optimize_it(DATA, x, 0.0_dp)
                        radii_pgf(jpgf, ipgf, jset, iset, jkind, ikind)%x = x
!$OMP END MASTER
!$OMP BARRIER
                     END DO !jpgf
                  END DO !ipgf
               END DO
            END DO
         END DO
      END DO
!$OMP MASTER
      CALL timestop(handle)
!$OMP END MASTER
   END SUBROUTINE calc_pair_dist_radii

!
!
! little driver routine for the powell minimizer
! data is the data to fit, x is of the form (x(1)*DATA(1)**2+x(2))
! only values of DATA(2) larger than fmin are taken into account
! it constructs an approximate upper bound of the fitted function
!
!
! **************************************************************************************************
!> \brief ...
!> \param DATA ...
!> \param x ...
!> \param fmin ...
! **************************************************************************************************
   SUBROUTINE optimize_it(DATA, x, fmin)

      REAL(KIND=dp), INTENT(IN)                          :: DATA(2, 0:100)
      REAL(KIND=dp), INTENT(OUT)                         :: x(2)
      REAL(KIND=dp), INTENT(IN)                          :: fmin

      INTEGER                                            :: i, k
      REAL(KIND=dp)                                      :: f, large_weight, small_weight, weight
      TYPE(opt_state_type)                               :: opt_state

! initial values

      x(1) = 0.0_dp
      x(2) = 0.0_dp

      ! we go in two steps, first we do the symmetric weight to get a good, unique initial guess
      ! we restart afterwards for the assym.
      ! the assym function appears to have several local minima, depending on the data to fit
      ! the loop over k can make the switch gradual, but there is not much need, seemingly
      DO k = 0, 4, 4

         small_weight = 1.0_dp
         large_weight = small_weight*(10.0_dp**k)

         ! init opt run
         opt_state%state = 0
         opt_state%nvar = 2
         opt_state%iprint = 3
         opt_state%unit = default_output_unit
         opt_state%maxfun = 100000
         opt_state%rhobeg = 0.1_dp
         opt_state%rhoend = 0.000001_dp

         DO

            ! compute function
            IF (opt_state%state == 2) THEN
               opt_state%f = 0.0_dp
               DO i = 0, 100
                  f = x(1)*DATA(1, i)**2 + x(2)
                  IF (f > DATA(2, i)) THEN
                     weight = small_weight
                  ELSE
                     weight = large_weight
                  END IF
                  IF (DATA(2, i) > fmin) opt_state%f = opt_state%f + weight*(f - DATA(2, i))**2
               END DO
            END IF

            IF (opt_state%state == -1) EXIT
            CALL powell_optimize(opt_state%nvar, x, opt_state)
         END DO

         ! dealloc mem
         opt_state%state = 8
         CALL powell_optimize(opt_state%nvar, x, opt_state)

      END DO

   END SUBROUTINE optimize_it

! **************************************************************************************************
!> \brief Given a 2d index pair, this function returns a 1d index pair for
!>        a symmetric upper triangle NxN matrix
!>        The compiler should inline this function, therefore it appears in
!>        several modules
!> \param i 2d index
!> \param j 2d index
!> \param N matrix size
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
!> \brief Evaluate screening function 
!> \param model - min(max(y_data), ax^8 + bx^4 + cx^2 + d); max(y_data) is the maximum value of the
!>                sampled data used to create the model
!> \param x2 - distance squared
!> \param y_min - default value if x2 is out of valid range i.e. [0, model%cutoff2]
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
PURE FUNCTION hfx_libheric_eval_screening_model(model, x2) RESULT(f)
   TYPE(hfx_libheric_screening_model_type), INTENT(IN) :: model
   REAL(dp), INTENT(IN)                                :: x2
   REAL(dp)                                            :: f
   
   REAL(dp)                                       :: x4, x8
   x4 = x2 * x2
   x8 = x4 * x4
   f = MERGE(MIN(model%max, model%x(1) + x2 * model%x(2) + x4 * model%x(3) + x8 * model%x(4)), model%y_min, x2 .LT. model%cutoff2)
END FUNCTION hfx_libheric_eval_screening_model

! **************************************************************************************************
!> \brief Compare error statistics (mean absolute error, mean squared error, max error) of new model
!>        to old model, and check if new model is a true upper bound to the fitted data
!> \param model_str - model identifier i.e. (la, lb, ikind, jkind, [iset, jset, [ipgf, jpgf]])
!> \param model - see hfx_libheric_eval_screen_model docs
!> \param ref_model - old model
!> \param n_data - number of data samples used to create model
!> \param x_data, y_data - data samples used to create model
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
SUBROUTINE hfx_libheric_display_screening_model_error(model_str, model, ref_model, n_data, x_data, y_data)
      CHARACTER(LEN=*), INTENT(IN)                        :: model_str
      TYPE(hfx_libheric_screening_model_type), INTENT(IN) :: model
      TYPE(hfx_screen_coeff_type), INTENT(IN)             :: ref_model
      INTEGER, INTENT(IN)                                 :: n_data
      REAL(dp), DIMENSION(:), INTENT(IN)                  :: x_data, y_data
      
      REAL(dp) :: mae_err, mse_err, inf_err, f, &
                  mae_ref, mse_ref, inf_ref, f_ref, &
                  ref_bad_fit_severity, new_bad_fit_severity, &
                  y_range
      INTEGER  :: i, ref_bad_fit_idx, new_bad_fit_idx
      CHARACTER(:), ALLOCATABLE :: error_str 


      mae_err = 0.0_dp
      mae_ref = 0.0_dp
      mse_err = 0.0_dp
      mse_ref = 0.0_dp
      inf_err = 0.0_dp
      inf_ref = 0.0_dp
      
      ref_bad_fit_idx = -1
      ref_bad_fit_severity = 0.0_dp
      new_bad_fit_idx = -1
      new_bad_fit_severity = 0.0_dp

      y_range = MAXVAL(y_data(1:n_data)) - MINVAL(y_data(1:n_data))
      DO i = 1, n_data
         f = hfx_libheric_eval_screening_model(model, x_data(i) ** 2)
         mae_err = mae_err + ABS(y_data(i) - f)
         mse_err = mse_err + (y_data(i) - f) ** 2
         inf_err = MAX(inf_err, ABS(y_data(i) - f))
         IF((y_data(i) - f) / y_range > new_bad_fit_severity) THEN
            new_bad_fit_severity = (y_data(i) - f) / y_range
            new_bad_fit_idx = i
         END IF

         f_ref = ref_model%x(1) * x_data(i) ** 2 + ref_model%x(2)
         mae_ref = mae_ref + ABS(y_data(i) - f_ref)
         mse_ref = mse_ref + (y_data(i) - f_ref) ** 2
         inf_ref = MAX(inf_ref, ABS(y_data(i) - f_ref))
         IF((y_data(i) - f) / y_range > ref_bad_fit_severity) THEN
            ref_bad_fit_severity = (y_data(i) - f) / y_range
            ref_bad_fit_idx = i
         END IF 
      END DO

      mae_err = mae_err / (y_range * n_data) 
      mae_ref = mae_ref / (y_range * n_data) 
      mse_err = mse_err / (y_range * n_data) 
      mse_ref = mse_ref / (y_range * n_data) 
      inf_err = mse_err / y_range 
      inf_ref = mse_ref / y_range

      IF(new_bad_fit_severity .GT. 1e-6) THEN
         ALLOCATE(CHARACTER(2 * n_data * 30 + 4096) :: error_str)
         WRITE(error_str, "(A,*(ES28.20,2X))") "x = vec([", x_data
         WRITE(error_str, "(A,/)") "])"
         WRITE(error_str, "(A,*(ES28.20,2X))") "y = vec([", y_data
         WRITE(error_str, "(A,/)") "])"
         WRITE(error_str, "(A,/,T3,5(A,ES16.8),A,/)") &
            TRIM(model_str), "model = min(" ,model%max, ", ", &
            model%x(4), " * x^8 + ", model%x(3), " * x^4 + ", model%x(2), " * x^2 + ", model%x(1), ")"
         WRITE(error_str, "(T6,A,I3,2(A,ES14.6),/)") "new model is bad fit for y(", new_bad_fit_idx, &
                                           ") with severity = ", new_bad_fit_severity
         WRITE(error_str, "(T6,3(A,ES14.6),/)") "mae_err = ", mae_err, ", mse_err = ", mse_err, ", inf_err = ", inf_err
         CALL cp_abort(__LOCATION__, TRIM(error_str))
      ELSE IF(mae_err .GT. mae_ref .OR. mse_err .GT. mse_ref .OR. inf_err .GT. inf_ref) THEN
         ALLOCATE(CHARACTER(2 * n_data * 30 + 4096) :: error_str)
         WRITE(error_str, "(A,*(ES28.20,2X))") "x = vec([", x_data
         WRITE(error_str, "(A,/)") "])"
         WRITE(error_str, "(A,*(ES28.20,2X))") "y = vec([", y_data
         WRITE(error_str, "(A,/)") "])"
         WRITE(error_str, "(A,/,T3,5(A,ES16.8),A,/)") &
            TRIM(model_str), "model = min(" ,model%max, ", ", &
            model%x(4), " * x^8 + ", model%x(3), " * x^4 + ", model%x(2), " * x^2 + ", model%x(1), ")"
         WRITE(error_str, "(T3,2(A,ES16.8,/))") "ref_model = ", ref_model%x(1), " * x^2 + ", ref_model%x(2)
         WRITE(error_str, "(T3,A,/)") "Reference model is better fit in at least one metric"
         IF(ref_bad_fit_severity > 1e-6) &
            WRITE(*, "(T6,A,I3,2(A,ES14.6),/)") "reference model is bad fit for y(", ref_bad_fit_idx, &
                                              ") with severity = ", ref_bad_fit_severity
         WRITE(error_str, "(T6,3(A,ES14.6),/)") "mae_err = ", mae_err, ", mse_err = ", mse_err, ", inf_err = ", inf_err
         WRITE(error_str, "(T6,3(A,ES14.6),/)") "mae_ref = ", mae_ref, ", mse_ref = ", mse_ref, ", inf_ref = ", inf_ref
         CALL cp_warn(__LOCATION__, TRIM(error_str))
      END IF
      
   END SUBROUTINE hfx_libheric_display_screening_model_error
   
! **************************************************************************************************
!> \brief Initialize screening parameters per angular momentum pair (la, lb); the model is documented
!>        in hfx_libheric_eval_screen_model
!> \param x_data - contains old models to compare the models to
!> \param max_set - maximum number of shellsets in any basis set
!> \param nkind - number unique basis sets
!> \param l_max - maximum angular momentum
!> \param basis_parameter_l - all hfx_libheric basis sets
!> \param screen - all screening parameters, i.e. for (R_pgf, pgf, sets, kinds)
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_init_screening_models(x_data, max_set, nkind, l_max, basis_parameter_l, screen)
      TYPE(hfx_type), INTENT(IN)                  :: x_data
      INTEGER, INTENT(IN)                         :: max_set, nkind, l_max
      TYPE(hfx_libheric_basis_type), &
          DIMENSION(1:nkind,0:l_max), INTENT(IN)  :: basis_parameter_l
      TYPE(hfx_libheric_screening_models), &
          INTENT(OUT)                             :: screen

      
      INTEGER, PARAMETER                          :: GRID_SIZE = 100

      REAL(dp), DIMENSION(0:GRID_SIZE, 2)         :: opt_data_kind, opt_data_set, opt_data_pgf!, &
      !                                               ref_data_kind, ref_data_set, ref_data_pgf
      !REAL(dp), DIMENSION(:, :), ALLOCATABLE      :: set_integrals, kind_integrals
      REAL(dp)                                    :: max_val, log_max_val, & 
                                                     ra(3), rb(3), rp(3), chi, ff, rab2, prefactor, &
                                                     max_contraction_a, max_contraction_b, &
                                                     pgf_radius, set_radius, kind_radius, cell_mat(3, 3)
      INTEGER                                     :: la, lb, ikind, jkind, iset, jset, ncoa, ncob, sgfa, sgfb, &
                                                     ipgf, jpgf, i, j, n_data, handle, nints, quartet_ctr
      INTEGER(int_8)                              :: mem_available
      CHARACTER(len=256)                          :: model_params
      TYPE(particle_type), DIMENSION(:), &
          ALLOCATABLE                             :: particle_set
      TYPE(acc_kind), DIMENSION(:), ALLOCATABLE   :: acc_kinds
      TYPE(acc_device), DIMENSION(:), ALLOCATABLE :: acc_devices
      TYPE(C_PTR), DIMENSION(:), ALLOCATABLE      :: ctx_handles
      TYPE(C_PTR)                                 :: local_ctx
      TYPE(atom_quartet_type)                     :: atom_quartet
      !REAL(dp), DIMENSION(:), POINTER             :: p_work
      !LOGICAL                                     :: assert_flag
      !REAL(dp)                                    :: R1
      
      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_libheric_init_screening_models'

      CALL timeset(routineN, handle)

      mem_available = 1024 * 1024 * 1024 ! 1GiB
      cell_mat(:, :) = 0.0_dp

      ALLOCATE(particle_set(nkind))
      !ALLOCATE(set_integrals(nso(l_max) ** 4, 0:GRID_SIZE), kind_integrals(nso(l_max) ** 4, 0:GRID_SIZE))
      DO ikind = 1, nkind
         particle_set(ikind)%r(:) = 0.0_dp
         particle_set(ikind)%atom_index = ikind
         ALLOCATE(particle_set(ikind)%atomic_kind)
         particle_set(ikind)%atomic_kind%kind_number = ikind
      END DO
      
      ALLOCATE(acc_kinds(0))
      ! Search for CPU accelerator is implicit
      CALL query_accelerator_devices(acc_kinds, acc_devices)
      CALL hfx_libheric_create_accelerator_contexts(max_set, nkind, l_max, basis_parameter_l, &
                                                    cell_mat, acc_devices(1), particle_set, &
                                                    mem_available, LOGICAL(.TRUE., KIND=C_BOOL), ctx_handles)
      local_ctx = ctx_handles(1)
      CALL accelerator_context_set_potential_operator(local_ctx, x_data%potential_parameter%potential_type, &
                                                      x_data%potential_parameter%omega)

      ALLOCATE(screen%coeffs(nkind, nkind, 0:l_max, 0:l_max))

      ra(:) = 0.0_dp
      rb(:) = 0.0_dp
      rp(:) = 0.0_dp

      DO la = 0, l_max
      DO lb = 0, l_max
         nints = (nco(la) * nco(lb)) ** 2
         DO ikind = 1, nkind
         DO jkind = 1, nkind
            opt_data_kind = -HUGE(max_val)
            !ref_data_kind = -HUGE(max_val)
            ASSOCIATE(bas_lai => basis_parameter_l(ikind, la), bas_lbj => basis_parameter_l(jkind, lb), &
                      coeff => screen%coeffs(jkind, ikind, lb, la))
            ALLOCATE(coeff%set(bas_lbj%nset, bas_lai%nset))
            ALLOCATE(coeff%pgf  (MAXVAL(bas_lbj%npgf), MAXVAL(bas_lai%npgf), bas_lbj%nset, bas_lai%nset))
            ALLOCATE(coeff%R_pgf(MAXVAL(bas_lbj%npgf), MAXVAL(bas_lai%npgf), bas_lbj%nset, bas_lai%nset))
            kind_radius = bas_lai%kind_radius + bas_lbj%kind_radius
            DO iset = 1, bas_lai%nset
            DO jset = 1, bas_lbj%nset
               ASSOCIATE(isetg => bas_lai%global_set(iset), jsetg => bas_lbj%global_set(jset)) 
               set_radius = bas_lai%set_radius(iset) + bas_lbj%set_radius(jset)
               opt_data_set = -HUGE(max_val)
               !ref_data_set = -HUGE(max_val)
               ncoa = bas_lai%npgf(iset) * nco(la)
               ncob = bas_lbj%npgf(jset) * nco(lb)
               sgfa = bas_lai%first_sgf(1, iset)
               sgfb = bas_lbj%first_sgf(1, jset)
               max_contraction_a = MAXVAL([(SUM(ABS(bas_lai%sphi(1:ncoa, i))), i=sgfa, sgfa + bas_lai%nsgf(iset) - 1)])
               max_contraction_b = MAXVAL([(SUM(ABS(bas_lbj%sphi(1:ncob, i))), i=sgfb, sgfb + bas_lbj%nsgf(jset) - 1)])
               !set_integrals(:, :) = 0.0_dp
               !kind_integrals(:, :) = 0.0_dp
               DO ipgf = 1, bas_lai%npgf(iset)
               DO jpgf = 1, bas_lbj%npgf(jset)
                  pgf_radius = bas_lai%pgf_radius(ipgf, iset) + bas_lbj%pgf_radius(jpgf, jset)
                  opt_data_pgf = -HUGE(max_val)
                  DO i = 0, GRID_SIZE 
                     !! ra is set to the origin, and rb only moves along the x-axis
                     rb(1) = (pgf_radius * i) / REAL(GRID_SIZE, KIND=dp)
                     chi = bas_lai%zet(ipgf, iset) + bas_lbj%zet(jpgf, jset)
                     ff = bas_lbj%zet(jpgf, jset) / chi
                     rab2 = rb(1) ** 2
                     prefactor = EXP(-bas_lai%zet(ipgf, iset) * ff * rab2)
                     rp(1) = ff*rb(1)
                     opt_data_pgf(i, 1) = rb(1)
                     opt_data_pgf(i, 2) = exp_radius_very_extended(la, la, lb, lb, ra=ra, rb=rb, rp=rp, &
                                                                   zetp=chi, eps=x_data%screening_parameter%eps_schwarz, &
                                                                   prefactor=prefactor, cutoff=1.0_dp, epsabs=1.0E-12_dp)

                  END DO
                  n_data = GRID_SIZE + 1
                  CALL libheric_optimize_screening_function(opt_data_pgf(0, 1), opt_data_pgf(0, 2), n_data, 0.0_dp, &
                                                            coeff%R_pgf(jpgf, ipgf, jset, iset)%x(1), 4)
                  coeff%R_pgf(jpgf, ipgf, jset, iset)%cutoff2 = ((pgf_radius * n_data) / REAL(GRID_SIZE, KIND=dp)) ** 2
                  coeff%R_pgf(jpgf, ipgf, jset, iset)%max     = MAXVAL(opt_data_pgf(:, 2)) 
                  coeff%R_pgf(jpgf, ipgf, jset, iset)%y_min   = 0.0_dp
                  WRITE(model_params, "(A,4(I1,A),6(I2,A))") "(la=", la, ", lb=", lb, ", ikind=", ikind, ", jkind=", jkind, &
                                                             ", iset=", iset, "(", isetg, &
                                                             "), jset=", jset, "(", jsetg, &
                                                             "), ipgf=", ipgf, ", jpgf=", jpgf,")"

                  !! Reference screening functions are only available if libint is available
#if(__LIBINT)
                  CALL hfx_libheric_display_screening_model_error(model_params, coeff%R_pgf(jpgf, ipgf, jset, iset), &
                                                                  x_data%pair_dist_radii_pgf(jpgf, ipgf, jsetg, &
                                                                                             isetg, jkind, ikind), &
                                                                  n_data, opt_data_pgf(:, 1), opt_data_pgf(:, 2)) 
#endif
                  cell_mat(1, 1) = pgf_radius / REAL(GRID_SIZE, KIND=dp)
                  cell_mat(1, 2) = set_radius / REAL(GRID_SIZE, KIND=dp)
                  cell_mat(1, 3) = kind_radius / REAL(GRID_SIZE, KIND=dp)

                  CALL accelerator_context_update_cell(local_ctx, cell_mat(1, 1))
                  quartet_ctr = 0
                  opt_data_pgf = -HUGE(max_val)
                  CALL accelerator_context_set_quartet_class(local_ctx, la, lb, la, lb)
                  IF(.NOT. accelerator_context_is_functional(local_ctx)) &
                      CPABORT("ACC not functional")
                  DO i = 0, GRID_SIZE
                     CALL accelerator_context_set_atom_ids(local_ctx, ikind - 1, jkind - 1, ikind - 1, jkind - 1)
                     CALL accelerator_context_set_set_ids(local_ctx, iset - 1, jset - 1, iset - 1, jset - 1, 0.0_dp)
                     CALL accelerator_context_set_shell_ids(local_ctx, 0, 0, 0, 0)
                     CALL accelerator_context_set_B_offset(local_ctx, i, 0, 0)
                     CALL accelerator_context_set_D_offset(local_ctx, i, 0, 0)
                     CALL accelerator_context_set_exponent_ids(local_ctx, ipgf - 1, jpgf - 1, ipgf - 1, jpgf - 1)
                     CALL accelerator_context_set_Q_offset(local_ctx, 0, 0, 0)
                     CALL accelerator_context_submit_peri(local_ctx)
                     
                     CALL accelerator_context_set_atom_ids(local_ctx, ikind - 1, jkind - 1, ikind - 1, jkind - 1)
                     CALL accelerator_context_set_set_ids(local_ctx, iset - 1, jset - 1, iset - 1, jset - 1, 0.0_dp)
                     CALL accelerator_context_set_shell_ids(local_ctx, 0, 0, 0, 0)
                     CALL accelerator_context_set_B_offset(local_ctx, 0, i, 0)
                     CALL accelerator_context_set_D_offset(local_ctx, 0, i, 0)
                     CALL accelerator_context_set_exponent_ids(local_ctx, ipgf - 1, jpgf - 1, ipgf - 1, jpgf - 1)
                     CALL accelerator_context_set_Q_offset(local_ctx, 0, 0, 0)
                     CALL accelerator_context_submit_peri(local_ctx)
                     
                     CALL accelerator_context_set_atom_ids(local_ctx, ikind - 1, jkind - 1, ikind - 1, jkind - 1)
                     CALL accelerator_context_set_set_ids(local_ctx, iset - 1, jset - 1, iset - 1, jset - 1, 0.0_dp)
                     CALL accelerator_context_set_shell_ids(local_ctx, 0, 0, 0, 0)
                     CALL accelerator_context_set_B_offset(local_ctx, 0, 0, i)
                     CALL accelerator_context_set_D_offset(local_ctx, 0, 0, i)
                     CALL accelerator_context_set_exponent_ids(local_ctx, ipgf - 1, jpgf - 1, ipgf - 1, jpgf - 1)
                     CALL accelerator_context_set_Q_offset(local_ctx, 0, 0, 0)
                     CALL accelerator_context_submit_peri(local_ctx)

                     DO WHILE(accelerator_context_get_next_staged_atom_quartet(local_ctx, atom_quartet))
                        j = quartet_ctr / 3
                        max_val = SQRT(MAX(MAXVAL(atom_quartet%integrals(1:nints, 1)), &
                                          -MINVAL(atom_quartet%integrals(1:nints, 1)))) * max_contraction_a * max_contraction_b
                        log_max_val = MERGE(powell_min_log, LOG10(max_val), max_val .EQ. 0.0_dp)
                        log_max_val = MAX(log_max_val, powell_min_log)
                        IF(MODULO(quartet_ctr, 3) .EQ. 0) THEN
                           opt_data_pgf(j, 1) = pgf_radius * j / REAL(GRID_SIZE, KIND=dp)
                           opt_data_pgf(j, 2) = log_max_val
                        ELSEIF(MODULO(quartet_ctr, 3) .EQ. 1) THEN
                           opt_data_set(j, 1) = set_radius * j / REAL(GRID_SIZE, KIND=dp)
                           opt_data_set(j, 2) = MAX(opt_data_set(j, 2), log_max_val)
                           !set_integrals(1:nints, j) = set_integrals(1:nints, j) + atom_quartet%integrals(1:nints, 1)
                        ELSE
                           opt_data_kind(j, 1) = kind_radius * j / REAL(GRID_SIZE, KIND=dp)
                           opt_data_kind(j, 2) = MAX(opt_data_kind(j, 2), log_max_val)
                           !kind_integrals(1:nints, j) = kind_integrals(1:nints, j) + atom_quartet%integrals(1:nints, 1)
                        END IF
                        quartet_ctr = quartet_ctr + 1
                     END DO
                  END DO
                  CALL accelerator_context_finish_quartet_class(local_ctx)
                  DO WHILE(accelerator_context_get_next_staged_atom_quartet(local_ctx, atom_quartet))
                     j = quartet_ctr / 3
                     max_val = SQRT(MAX(MAXVAL(atom_quartet%integrals(1:nints, 1)), &
                                       -MINVAL(atom_quartet%integrals(1:nints, 1)))) * max_contraction_a * max_contraction_b
                     log_max_val = MERGE(powell_min_log, LOG10(max_val), max_val .EQ. 0.0_dp)
                     log_max_val = MAX(log_max_val, powell_min_log)
                     IF(MODULO(quartet_ctr, 3) .EQ. 0) THEN
                        opt_data_pgf(j, 1) = pgf_radius * j / REAL(GRID_SIZE, KIND=dp)
                        opt_data_pgf(j, 2) = log_max_val
                     ELSEIF(MODULO(quartet_ctr, 3) .EQ. 1) THEN
                        opt_data_set(j, 1) = set_radius * j / REAL(GRID_SIZE, KIND=dp)
                        opt_data_set(j, 2) = MAX(opt_data_set(j, 2), log_max_val)
                        !set_integrals(1:nints, j) = set_integrals(1:nints, j) + atom_quartet%integrals(1:nints, 1)
                     ELSE
                        opt_data_kind(j, 1) = kind_radius * j / REAL(GRID_SIZE, KIND=dp)
                        opt_data_kind(j, 2) = MAX(opt_data_kind(j, 2), log_max_val)
                        !kind_integrals(1:nints, j) = kind_integrals(1:nints, j) + atom_quartet%integrals(1:nints, 1)
                     END IF
                     quartet_ctr = quartet_ctr + 1
                  END DO
                  !ref_data_pgf = -HUGE(max_val)
                  !ASSOCIATE(rad_x => coeff%R_pgf(jpgf, ipgf, jset, iset)%x, &
                  !          zeta => bas_lai%zet(ipgf, iset), zetb => bas_lbj%zet(jpgf, jset))
                  !DO i = 0, GRID_SIZE 
                  !   !! ra is set to the origin, and rb only moves along the x-axis
                  !   !! Collect data for pgf_coeff optimization
                  !   rb(1) = (pgf_radius * i) / REAL(GRID_SIZE, KIND=dp)
                  !   R1 = MAX(0.0_dp, rad_x(1) * rb(1) ** 2 + rad_x(2))
                  !   max_val = 0.0_dp 
                  !   CALL evaluate_eri_screen(x_data%lib, ra, rb, ra, rb, zeta, zetb, zeta, zetb, &
                  !                            la, lb, la, lb, max_val, x_data%potential_parameter, R1, R1, p_work)
                  !   max_val = SQRT(max_val) * max_contraction_a * max_contraction_b
                  !   log_max_val = MERGE(powell_min_log, LOG10(max_val), max_val .EQ. 0.0_dp)
                  !   log_max_val = MAX(log_max_val, powell_min_log)
                  !   ref_data_pgf(i, 1) = rb(1)
                  !   ref_data_pgf(i, 2) = log_max_val

                  !   !! Collect data for set_coeff optimization
                  !   rb(1) = REAL(i, dp) / REAL(GRID_SIZE, dp) * set_radius
                  !   R1 = MAX(0.0_dp, rad_x(1) * rb(1) ** 2 + rad_x(2))
                  !   max_val = 0.0_dp 
                  !   CALL evaluate_eri_screen(x_data%lib, ra, rb, ra, rb, zeta, zetb, zeta, zetb, &
                  !                            la, lb, la, lb, max_val, x_data%potential_parameter, R1, R1, p_work)
                  !   max_val = SQRT(max_val) * max_contraction_a * max_contraction_b
                  !   log_max_val = MERGE(powell_min_log, LOG10(max_val), max_val .EQ. 0.0_dp)
                  !   log_max_val = MAX(log_max_val, powell_min_log)
                  !   ref_data_set(i, 1) = rb(1)
                  !   ref_data_set(i, 2) = MAX(log_max_val, ref_data_set(i, 2))
                  !   
                  !   !! Collect data for kind_coeff optimization
                  !   rb(1) = (kind_radius * i) / REAL(GRID_SIZE, KIND=dp)
                  !   R1 = MAX(0.0_dp, rad_x(1) * rb(1) ** 2 + rad_x(2))
                  !   max_val = 0.0_dp 
                  !   CALL evaluate_eri_screen(x_data%lib, ra, rb, ra, rb, zeta, zetb, zeta, zetb, &
                  !                            la, lb, la, lb, max_val, x_data%potential_parameter, R1, R1, p_work)
                  !   max_val = SQRT(max_val) * max_contraction_a * max_contraction_b
                  !   log_max_val = MERGE(powell_min_log, LOG10(max_val), max_val .EQ. 0.0_dp)
                  !   log_max_val = MAX(log_max_val, powell_min_log)
                  !   ref_data_kind(i, 1) = rb(1)
                  !   ref_data_kind(i, 2) = MAX(log_max_val, ref_data_kind(i, 2))
                  !END DO
                  !END ASSOCIATE
                 
                  !assert_flag = ALL([(ABS(ref_data_pgf(i, 2) - opt_data_pgf(i, 2)) .LT. 1e-9, i=0,GRID_SIZE)])
                  !CPASSERT(assert_flag)
                   
                  n_data = GRID_SIZE + 1
                  CALL libheric_optimize_screening_function(opt_data_pgf(0, 1), opt_data_pgf(0, 2), n_data, powell_min_log, &
                                                            coeff%pgf(jpgf, ipgf, jset, iset)%x(1), 4)
                  coeff%pgf(jpgf, ipgf, jset, iset)%cutoff2 = ((pgf_radius * n_data) / REAL(GRID_SIZE, KIND=dp)) ** 2
                  coeff%pgf(jpgf, ipgf, jset, iset)%max     = MAXVAL(opt_data_pgf(:, 2)) 
                  coeff%pgf(jpgf, ipgf, jset, iset)%y_min   = powell_min_log
                  !! Reference screening functions are only available if libint is available
#if(__LIBINT)
                  CALL hfx_libheric_display_screening_model_error(model_params, coeff%pgf(jpgf, ipgf, jset, iset), &
                                                                  x_data%screen_funct_coeffs_pgf(jpgf, ipgf, jsetg, &
                                                                                                 isetg, jkind, ikind), &
                                                                  n_data, opt_data_pgf(:, 1), opt_data_pgf(:, 2)) 
#endif
               END DO
               END DO
               !assert_flag = ALL([(ABS(ref_data_set(i, 2) - opt_data_set(i, 2)) .LT. 1e-9, i=0,GRID_SIZE)])
               !CPASSERT(assert_flag)

               n_data = GRID_SIZE + 1
               CALL libheric_optimize_screening_function(opt_data_set(0, 1), opt_data_set(0, 2), n_data, powell_min_log, &
                                                             coeff%set(jset, iset)%x(1), 4)
               coeff%set(jset, iset)%cutoff2 = ((set_radius * n_data) / REAL(GRID_SIZE, KIND=dp)) ** 2
               coeff%set(jset, iset)%max     = MAXVAL(opt_data_set(:, 2)) 
               coeff%set(jset, iset)%y_min   = powell_min_log
               WRITE(model_params, "(A,4(I1,A),4(I2,A))") "(la=", la, ", lb=", lb, ", ikind=", ikind, ", jkind=", jkind, &
                                                          ", iset=", iset, "(", isetg, &
                                                          "), jset=", jset, "(", jsetg, &
                                                          "))"
               !! Reference screening functions are only available if libint is available
#if(__LIBINT)
               CALL hfx_libheric_display_screening_model_error(model_params, coeff%set(jset, iset), &
                                                               x_data%screen_funct_coeffs_set(jsetg, isetg, jkind, ikind), &
                                                               n_data, opt_data_set(:, 1), opt_data_set(:, 2)) 
#endif
               END ASSOCIATE
            END DO
            END DO
            !assert_flag = ALL([(ABS(ref_data_kind(i, 2) - opt_data_kind(i, 2)) .LT. 1e-9, i=0,GRID_SIZE)])
            !CPASSERT(assert_flag)
            n_data = GRID_SIZE + 1
            CALL libheric_optimize_screening_function(opt_data_kind(0, 1), opt_data_kind(0, 2), n_data, powell_min_log, &
                                             coeff%kind%x(1), 4)
            coeff%kind%cutoff2 = ((kind_radius * n_data) / REAL(GRID_SIZE, KIND=dp)) ** 2
            coeff%kind%max     = MAXVAL(opt_data_kind(:, 2)) 
            coeff%kind%y_min   = powell_min_log
            WRITE(model_params, "(A,4(I1,A))") "(la=", la, ", lb=", lb, ", ikind=", ikind, ", jkind=", jkind, ")"
            !! Reference screening functions are only available if libint is available
#if(__LIBINT)
            CALL hfx_libheric_display_screening_model_error(model_params, coeff%kind, &
                                                            x_data%screen_funct_coeffs_kind(jkind, ikind), &
                                                            n_data, opt_data_kind(:, 1), opt_data_kind(:, 2)) 
#endif

            END ASSOCIATE
         END DO
         END DO
     END DO
     END DO
     CALL free_accelerator_contexts(ctx_handles(1), SIZE(ctx_handles))
     CALL timestop(handle)
   END SUBROUTINE hfx_libheric_init_screening_models
   
! **************************************************************************************************
!> \brief Deallocate memory of all screening parameters  
!> \param screen - all screening parameters 
!> \par History
!>      09.2026 created [Johannes Menzel]
!> \author Johannes Menzel
! **************************************************************************************************
   SUBROUTINE hfx_libheric_free_screening_models(screen)
      TYPE(hfx_libheric_screening_models), INTENT(INOUT) :: screen
      INTEGER :: la, lb, ikind, jkind
      IF(ALLOCATED(screen%coeffs)) THEN 
         DO la = 0, UBOUND(screen%coeffs, 4)
         DO lb = 0, UBOUND(screen%coeffs, 3)
         DO ikind = 1, UBOUND(screen%coeffs, 2)
         DO jkind = 1, UBOUND(screen%coeffs, 1)
            IF(ALLOCATED(screen%coeffs(jkind, ikind, lb, la)%set))   DEALLOCATE(screen%coeffs(jkind, ikind, lb, la)%set)
            IF(ALLOCATED(screen%coeffs(jkind, ikind, lb, la)%pgf))   DEALLOCATE(screen%coeffs(jkind, ikind, lb, la)%pgf)
            IF(ALLOCATED(screen%coeffs(jkind, ikind, lb, la)%R_pgf)) DEALLOCATE(screen%coeffs(jkind, ikind, lb, la)%R_pgf)
         END DO
         END DO
         END DO
         END DO
         DEALLOCATE(screen%coeffs)
      END IF
   END SUBROUTINE hfx_libheric_free_screening_models

END MODULE hfx_screening_methods
