# 1 "elphon_selfen_ph.F"
# 1 "./symbol.inc" 1 

!-------- to be costumized by user (usually done in the makefile)-------
!#define vector              compile for vector machine
!#define essl                use ESSL instead of LAPACK
!#define single_BLAS         use single prec. BLAS

!#define wNGXhalf            gamma only wavefunctions (X-red)
!#define wNGZhalf            gamma only wavefunctions (Z-red)

!#define NGXhalf             charge stored in REAL array (X-red)
!#define NGZhalf             charge stored in REAL array (Z-red)
!#define NOZTRMM             replace ZTRMM by ZGEMM
!#define 1                 compile for parallel machine with 1
!------------- end of user part --------------------------------
# 17

# 62

# 91

!
!   charge density: full grid mode
!










# 113

!
!   charge density complex
!






# 133

!
!   wavefunctions: full grid mode
!

# 182

!
!   wavefunctions complex
!








































!
!   common definitions
!








!
!   mpi parallel macros
!







# 268

# 274

# 279

# 286

# 293



# 302







# 317











# 339









!
! OpenMP macros
!
# 354

# 357

# 366









!
! profiling macros
!
# 381




!
! shmem macros
!
# 390

# 393

# 396

!
! quadruple precision
!
# 414











!
! for the 1 interface
!
# 430

!
! CUDA includes
!
# 438

!
! SIMD related definitions
!
# 1 "./simd.inc" 1 
!!#if   defined(__MIC__) || defined(__AVX512F__)
!!#define SIMD512
!!#undef  SIMD256
!!#elif defined(__AVX__) || defined(__AVX2__)
!!#define SIMD256
!!#undef  SIMD512
!!#endif

# 17





# 25





# 35















# 54





# 443 "./symbol.inc" 2 
!
! memalign macros
!
# 456




!
! Macros for PGI/NV HPC compilers version specific code
!
# 466

# 473



# 485










# 497


# 501


!
! OpenACC macros
!
# 527



































# 568












































!
! combined OpenMP and OpenACC macros
!
# 617



!
! routines replaced in LAPACK >=3.6
!
# 625


!
! Macros for HDF5 error check
!


# 634


!
! Macros for GNU version specific code
!
# 641


!
! For machine learning
!




!
! Extra safe initializations (+ overflow protections)
!
# 655




!
! Macros for memory estimation
!




!
! Offloading related macros
!
# 671


# 695

! Line is included when  !OFFLOADING

! Line is a comment when !OFFLOADING



!replace blas and lapack wrapper calls with
!normal calls if no offloading is used:





































# 747









# 759

! Line is a comment when !_OPENACC or   ACC_OFFLOAD



# 769

! Line is a comment when !ACC_OFFLOAD



# 779

! Line is included when  !_OPENACC

! Line is a comment when !_OPENACC



# 789

! Line is a comment when !_OPENMP or   OMP_OFFLOAD



# 799

! Line is a comment when !OMP_OFFLOAD



# 809

! Line is included when  !_OPENMP

! Line is a comment when !_OPENMP


# 2 "elphon_selfen_ph.F" 2 

module elphon_selfen_ph
    use prec
    use constant
    implicit none

    private
    public electron_phonon_selfen_ph_driver
    public elphon_prepare_selfenergy_accumulators_ph
    public elph_accumulators_init_ph
    public elph_accumulators_selfen_ph
    public elph_mels_cache_prefill_ph
    public elph_mels_wf_redistribute_ph
    public elph_accumulators_finalize_ph
    public elph_accumulators_print_ph
    public elph_accumulators_free_ph

    type elph_accumulator_selfen_ph
        integer :: nbands !< number of bands at which to compute the self-energy
        integer :: nbands_sum !< number of bands to sum over
        integer :: natoms !< number of atoms in the system
        integer :: nkpoints !< number of kpoints
        integer :: nspin !< number of spin channels
        integer :: nw !< number of energies at which to evaluate the self-energy
        real(q) :: enwin !< range of energies
        real(q) :: delta !< small imaginary part to add to the self-energy
        real(q), allocatable :: carrier_per_cell(:) !< number of electrons per cell (for the corresponding temperatures)
        real(q) :: carrier_per_cell0 !< number of electrons per cell (without doping)
        real(q), allocatable :: efermi(:) !< fermi energies (for the corresponding temperatures)
        real(q) :: emin !< energy of the lowest band which will be used in the sum-over-states
        real(q) :: emax !< energy of the highest band which will be used in the sum-over-states
        integer :: band_start !< index of the band at which to start computing the self-energy
        integer :: band_stop  !< index of the band at which to stop computing the self-energy
        logical :: tetrahedron !< use tetrahedron integration method
        real(q) :: cbm !< conduction band minimum
        real(q) :: vbm !< valence band maximum
        logical :: ismetal !< logical to know wether the material is a metal or not
!> Energy cutoff (eV) for choosing which states to compute.
!> There are two cases:
!>  - metal: the self-energy for all the states within
!>           [efermi-select_energy_window(1):efermi+select_energy_window(2)] is computed
!>  - gaped system: all the states within [vbm-select_energy_window(1):cbm+select_energy_window(1)]
        real(q) :: select_energy_window(2)
        integer :: energy_grid_mode !< definition of the energy grid: 0-centered around KS eigenvalue, 1-in the range of total energies
!> energy grids at which to evaluate the self-energy in eV (nw,nbks)
!> where nbks is a combined index for bands, kpoints and spin (see bks_idx)
        real(q), allocatable :: energies(:, :)
!> store the fan self-energy. The dimensions are (ntemps,nw,nbks)
!> where nbks is a combined index for bands, kpoints and spin (see bks_idx)
        complex(q), allocatable :: t(:, :, :, :, :, :)
        integer :: ntemps !< number of temperatures
        real(q), allocatable :: tempsk(:) !< array of temperatures in kelvin
        logical, allocatable :: compute_mask(:, :, :) !< integer controlling whether to compute self-energy at ikpt and band
        integer :: nbks !< Number of combined bands, kpoints and spin indexes (see bks_idx)
        integer, allocatable :: bks_idx(:, :, :) !< Combine band, kpoint and spin index into (1._q,0._q)
!> displacements of 24 tetrahedra that contribute to a k-point
        integer :: kcut24(3, 4, 24)
!> index of the corner of each of the 24 tetrahedra for which the k-point has 0 displacement
        integer :: kcut24_origin(24)
        logical :: kcut24_initialized
!> store the phonon self-energy. The dimensions are (ntemps, nw, n-irkpoints, n-phbands)
!> where nbks is a combined index for bands, kpoints and spin (see bks_idx)
        complex(q), allocatable :: selfen_ph(:, :, :, :)
!> phonon energies at which we compute the phonon self energy (nw,nb_q_ir , 3 * nbatoms)
        real(q), allocatable :: phonon_energies(:, :, :)
    end type elph_accumulator_selfen_ph

    type elph_accumulators_selfen_ph
        integer :: naccumulators_ph
!> global index of the local bands
        integer, allocatable :: band_idx(:)
!> phonon self-energy accumulators
        type(elph_accumulator_selfen_ph), allocatable :: selfen_ph(:)
    end type elph_accumulators_selfen_ph

contains

!> @brief Private subroutine to initialize the self-energy accumulators
    subroutine elphon_prepare_selfenergy_accumulators_ph(selfen, latt_cur, kpoints_inter, elph_set, wdes_dense_ibz, &
                                                         w_dense_ibz, t_info, info, chempot, io)
        use base, only: in_struct, info_struct
        use tutor, only: vtutor
        use string, only: str
        use mkpoints_struct_def, only: kpoints_struct
        use mkpoints, only: kpoint_integers, kpoints_along_path, kpoints_along_path_write, kint_kindex
        use lattice, only: latt
        use poscar, only: type_info
        use wave_struct_def, only: wavedes, wavespin
        use vhdf5
        use elphon_base
        use electron_transport, only: chemical_potential
        implicit none

        type(elph_accumulator_selfen_ph), intent(inout) :: selfen(:)
        type(latt), intent(in) :: latt_cur
        type(kpoints_struct), intent(in) :: kpoints_inter
        type(elph_settings), intent(inout) :: elph_set
        type(wavedes), intent(in) :: wdes_dense_ibz
        type(wavespin), target, intent(inout) :: w_dense_ibz
        type(type_info), intent(in) :: t_info
        type(info_struct), intent(in) :: info
        type(chemical_potential), intent(in) :: chempot
        type(in_struct), intent(in) :: io

! local variables
        integer :: idx, id, i, iden, ikpt, ikfbz
        integer :: kint(3)
        logical :: kpoints_line_exists
        real(q), allocatable :: kpoints_dists(:)
        real(q), allocatable :: kpoints_vkpt(:, :)

        integer(HID_T) :: groupid


        

! Check if ikpts is allocated, if not use all the points
        if (.not. allocated(elph_set%selfen_ikpt)) then
            allocate (elph_set%selfen_ikpt(wdes_dense_ibz%nkpts))
            elph_set%selfen_ikpt = [(i, i=1, wdes_dense_ibz%nkpts, 1)]
        end if

! Check if ikpts values are within range
        if (any(elph_set%selfen_ikpt <= 0 .or. elph_set%selfen_ikpt > wdes_dense_ibz%nkpts)) then
            call vtutor%error('Invalid IKPT indexes in ELPH_SELFEN_IKPT: '//str(elph_set%selfen_ikpt)// &
                              '. Should be between 1 and '//str(wdes_dense_ibz%nkpts))
        END IF

! Now if ELPH_SELFEN_KPTS is set then we use these k-point coordinates to set the SELFEN_IKPT array
        if (allocated(elph_set%selfen_kpts)) then
            deallocate (elph_set%selfen_ikpt)
            allocate (elph_set%selfen_ikpt(size(elph_set%selfen_kpts, 2)))
            do i = 1, size(elph_set%selfen_kpts, 2)
! Since the coordinates are input by the user, use a lower threshold for comparing points
                call kpoint_integers(kpoints_inter, latt_cur, elph_set%selfen_kpts(:, i), kint, eps=1e-3_q)
                if (all(kint == -1)) call vtutor%error('Could not find '//str(elph_set%selfen_kpts(:, i))//' in the regular grid')
                ikfbz = kint_kindex(kpoints_inter, kint)
                elph_set%selfen_ikpt(i) = kpoints_inter%fbz2ibz(1, ikfbz)
            end do
! some output
            if (io%iu6 > 0) then
                write (io%iu6, *) ' SELFEN_KPTS               IKPT  MATCHED TO'
                do i = 1, size(elph_set%selfen_kpts, 2)
                    ikpt = elph_set%selfen_ikpt(i)
                    write (io%iu6, '(3f8.4,a,i6,3f8.4)') elph_set%selfen_kpts(:, i), '->', ikpt, kpoints_inter%vkpt(:, ikpt)
                end do
            end if
        end if

        inquire (file="KPOINTS_LINE", exist=kpoints_line_exists)
        if (kpoints_line_exists) then
            deallocate (elph_set%selfen_ikpt)
            call kpoints_along_path(kpoints_inter, latt_cur, kpoints_dists, kpoints_vkpt, elph_set%selfen_ikpt)
            call kpoints_along_path_write(elph_set%selfen_ikpt, kpoints_dists, kpoints_vkpt, io%iu6)

            call vh5_error(vh5_group_open_or_create(ih5outfileid, trim(grp_results), groupid),"elphon_selfen_ph.F",157)
            call kpoints_along_path_hdf5write(groupid, elph_set%selfen_ikpt, kpoints_dists, kpoints_vkpt, group='electron_phonon')

        end if

        idx = 0
! loop over broadening
        do id = 1, size(elph_set%selfen_delta)
! loop over number of bands to sum
            do i = 1, size(elph_set%nbands_sum)
! loop over chemical potentials
                do iden = 1, size(chempot%selfen_muij, 2)
                    idx = idx + 1
                    if (idx > elph_set%selfen_naccumulators) then
                        CALL vtutor%bug('Index is larger than number of accumulators', "elphon_selfen_ph.F", 171)
                    end if
                    call elph_accumulator_selfen_init_ph(selfen(idx), wdes_dense_ibz%nb_tot, &
                                                         elph_set%nbands_sum(i), wdes_dense_ibz%nkpts, &
                                                         wdes_dense_ibz%ispin, elph_set%selfen_nw, &
                                                         elph_set%selfen_wrange, &
                                                         elph_set%selfen_delta(id), elph_set%selfen_ntemps, &
                                                         elph_set%selfen_temps, &
                                                         elph_set%selfen_band_start, elph_set%selfen_band_stop, &
                                                         t_info%nions, w_dense_ibz%celtot, wdes_dense_ibz%nb_totk, &
                                                         chempot%selfen_muij(:, iden), chempot%selfen_nij(:, iden), info%nelect)
                    call elph_accumulator_selfen_init_tetrahedron_ph(selfen(idx), kpoints_inter)
                    if (size(elph_set%selfen_ikpt) > 0) then
                        call elph_accumulator_selfen_selector_ikpt_ph(selfen(idx), elph_set%selfen_ikpt)
                    end if
                    if (any(elph_set%selfen_energy_window > 0)) then
                        call elph_accumulator_selfen_selector_energy_window_ph( &
                            selfen(idx), w_dense_ibz, elph_set%selfen_energy_window)
                    end if
                    if (elph_set%selfen_gaps) then
                        call elph_accumulator_selfen_selector_gaps_ph(selfen(idx), w_dense_ibz)
                    end if
                    call elph_accumulator_selfen_sparsify_ph(selfen(idx), w_dense_ibz%celtot)
                end do ! loop over chemical potentials
            end do ! loop over number of bands
        end do ! loop over broadening

! write some information to the output
        do i = 1, elph_set%selfen_naccumulators
            call elph_accumulator_selfen_print_ph(selfen(i), i, io%iu6)
        end do
        

    end subroutine elphon_prepare_selfenergy_accumulators_ph

!> redistribute the k-points of the wfs according to the work that needs to be performed
    subroutine elph_mels_wf_redistribute_ph(self, elph_accumulators, kpoints)
        use elphon_mels, only: wf_redistribute_init, wf_redistribution, elph_mels_type, &
                               elph_mels_get_ikpfbz, wf_redistribute_write_workload_cpu, &
                               wf_redistribute_write_workload_kpoint, wf_redistribute_workload, &
                               wf_redistribute_wfs_window, elph_mels_windows_init, wf_redistribute_free
        use mkpoints_struct_def, only: kpoints_struct
        implicit none

        type(elph_mels_type) :: self
        type(elph_accumulators_selfen_ph) :: elph_accumulators
        type(kpoints_struct) :: kpoints
!local variables
        type(wf_redistribution) :: wf_redis

        if (.not. self%elph_set%wf_redistribute) return

        
        call wf_redistribute_init(wf_redis, kpoints%nkpts, self%comm_kbatch)
        call elph_mels_get_ikpfbz(self, kpoints, self%ikp_fbz_list)
        call wf_redistribute_compute_workload_ph( &
            wf_redis, self%elph_set, elph_accumulators, self%ikp_fbz_list, self%fbz2ibz, self%w_dense_ibz, self%wdes_dense_ibz, &
            kpoints, self%phonon_freqs_ibz, self%latt_cur, self%wave_map, self%comm_kbatch)
        call wf_redistribute_write_workload_cpu('current', wf_redis%workload_orig, wf_redis%nkpts_batch_orig, self%io%iu0)
        call wf_redistribute_write_workload_kpoint(wf_redis, self%io%iu0)
        call wf_redistribute_workload(wf_redis, self%wdes_dense_ibz, kpoints, self%comm_kbatch)
        call wf_redistribute_write_workload_cpu('new', wf_redis%workload, wf_redis%nkpts_batch, self%io%iu0)

        call wf_redistribute_wfs_window(wf_redis, self%nrspinors, self%nproj, self%wavwin_cw, self%wavwin_cproj, &
                                 self%latt_cur, self%latt_ini, self%wdes_dense_ibz, self%wave_map, self%wdes_kbatch, &
                                 self%w_kbatch, self%dproj_kbatch)

        call elph_mels_windows_init(self)

        call wf_redistribute_free(wf_redis)
        
    end subroutine elph_mels_wf_redistribute_ph

!> compute workload given the current WF distribution
    subroutine wf_redistribute_compute_workload_ph(self, elph_set, elph_accumulators, ikp_fbz_list, fbz2ibz, &
                                                   w_dense_ibz, wdes_dense_ibz, kpoints, phonon_freqs_ibz, &
                                                   latt_cur, wave_map, comm_kbatch)
        use mkpoints_struct_def, only: kpoints_struct
        use triplets, only: triplet_type, kpoint_init_kindex, kpoint_init_kibz, kpoint_init_kint
        use elphon_base, only: elph_settings
        use elphon_mels, only: wf_redistribution
        use wave_rotate, only: fbz2ibz_stars
        use wave_struct_def, only: wavespin, wavedes
        use poscar_struct_def, only: latt
        use wave_interpolate, only: wave_mapper
        use mpimy, only: communic
        implicit none

        type(wf_redistribution) :: self
        type(elph_settings) :: elph_set
        type(elph_accumulators_selfen_ph) :: elph_accumulators
        integer, intent(in) :: ikp_fbz_list(:)
        type(fbz2ibz_stars), intent(in) :: fbz2ibz !< contains mapping from IBZ to FBZ using stars
        type(wavespin), intent(in) :: w_dense_ibz !< wavefunction storage in the dense ibz
        type(wavedes), intent(in) :: wdes_dense_ibz !< wavefunction descriptor for the dense ibz
        type(kpoints_struct), intent(in) :: kpoints
        real(q), intent(in) :: phonon_freqs_ibz(:, :)
        type(latt), intent(in) :: latt_cur
        type(wave_mapper), intent(in) :: wave_map !< map local WFs to points in the ibz
        type(communic), intent(in) :: comm_kbatch
! local variables
        type(triplet_type) :: triplet
        integer :: ik, ik_ibz, isp, workload
        integer :: ikp_fbz
        logical :: delta_is_zero

        
!do_iqfbz = 0
        do isp = 1, wdes_dense_ibz%ispin
! #ifdef OLD_SELECTION
!             ! HM: here we only check conbinations of k' and k with both k and k' in the IBZ.
!             ! I believe this is not coorect. The right thing is to check for k' in the FBZ and k in the IBZ as is 1._q below
!             ! This version is naturally faster (depending on the symmetry of the system).
!             ! I did not observe large changes in the final results but I cannot explain why this would not make a difference.
!             ! loop over the local ibz k' points
!             do ik = wave_map%nkpts_orig + 1, wave_map%nkpts_batch
!                 ikp_ibz = wave_map%kpoints_index(ik)
!                 call kpoint_init_kibz(triplet%kp, kpoints, latt_cur, ikp_ibz)
! #else
! loop over all local fbz k' points
            do ik = 1, size(ikp_fbz_list)
                ikp_fbz = ikp_fbz_list(ik)
                call kpoint_init_kindex(triplet%kp, kpoints, ikp_fbz)
! #endif
! loop over all the k points
                do ik_ibz = 1, wdes_dense_ibz%nkpts
                    if (.not. elph_accumulators_shallido_ph(elph_accumulators, ik_ibz, isp)) cycle

                    call kpoint_init_kibz(triplet%k, kpoints, latt_cur, ik_ibz)
                    call kpoint_init_kint(triplet%q, kpoints, triplet%k%kint - triplet%kp%kint)

                    delta_is_zero = elph_accumulators_delta_is_zero_ph( &
                                    elph_accumulators, elph_set, w_dense_ibz, phonon_freqs_ibz, &
                                    kpoints, triplet%k%ibz, isp, triplet%kp%kint, triplet%q%kint)
                    if (delta_is_zero) cycle

!do_iqfbz(triplet%q%kindex) = &
!do_iqfbz(triplet%q%kindex) + 1

! #ifdef OLD_SELECTION
!                     workload = fbz2ibz%star_size(ikp_ibz)
! #else
                    workload = 1
! #endif
                    self%workload_orig(comm_kbatch%node_me) = &
                        self%workload_orig(comm_kbatch%node_me) + workload
                    self%workload_ikp(triplet%kp%ibz) = &
                        self%workload_ikp(triplet%kp%ibz) + workload
                end do !k
            end do !k'
        end do !isp

!write(*,*) 'will do', count(do_iqfbz/=0), '/', size(do_iqfbz), &
!           'max:', maxval(do_iqfbz), &
!           'avg:', 1.0_q*sum(do_iqfbz)/count(do_iqfbz/=0), &
!           'potential savings:', 100.0_q*count(do_iqfbz/=0)/sum(do_iqfbz)
        self%nkpts_batch_orig(comm_kbatch%node_me) = wave_map%nkpts_batch

        CALL M_sum_i(comm_kbatch, self%workload_orig, size(self%workload_orig))
        CALL M_sum_i(comm_kbatch, self%nkpts_batch_orig, size(self%nkpts_batch_orig))
        CALL M_sum_i(comm_kbatch, self%workload_ikp, size(self%workload_ikp))
        
    end subroutine wf_redistribute_compute_workload_ph

!> Fill in the WF cache with the states that should be present on all the nodes
    subroutine elph_mels_cache_prefill_ph(self, elph_accumulators)
        use wave_windower, only: wave_window_get, wave_cache_set, wave_window_flush_local_all
        use elphon_base, only: find_value
        use elphon_mels, only: elph_mels_type, elph_mels_windows_free, elph_mels_windows_init
        implicit none

        type(elph_mels_type) :: self
        type(elph_accumulators_selfen_ph), intent(in) :: elph_accumulators
! local variables
        integer :: nbk_global, isp, ikp_fbz, ik_ibz, idx
        logical :: has_wfk_locally
        complex(q), allocatable :: cw1(:)
        complex(q), allocatable :: cproj1(:)

        if (.not. self%elph_set%wf_cache_prefill) return

        
        allocate (cw1(self%wdes_dense_ibz%nrplwv))
        allocate (cproj1(self%nproj*4*self%nrspinors))
! Gather all the required WFs
        do isp = 1, self%wdes_dense_ibz%ispin
            ikp_fbz = 1
            do ik_ibz = 1, self%wdes_dense_ibz%nkpts
                if (.not. elph_accumulators_shallido_ph(elph_accumulators, ik_ibz, isp)) cycle
                idx = find_value(self%wave_map%kpoints_index, ik_ibz)
                has_wfk_locally = idx /= 0
! if I have these orbitals locally then skip
! otherwise get orbitals for this k-point
                do nbk_global = self%band_start_k, self%band_stop_k
! check if I have the bands locally otherwise get them from the node that has them
                    if (mod(nbk_global - 1, self%wdes_dense_ibz%nb_par) + 1 == self%wdes_dense_ibz%nb_low .and. &
                        has_wfk_locally) cycle
                    call wave_window_get(self%wavwin_cw, nbk_global, ik_ibz, isp, cw1)
                    call wave_cache_set(self%cache_cw, nbk_global, ik_ibz, isp, cw1)
                    call wave_window_get(self%wavwin_cproj, nbk_global, ik_ibz, isp, cproj1)
                    call wave_cache_set(self%cache_cproj, nbk_global, ik_ibz, isp, cproj1)
                end do
            end do
        end do
! TODO: this is needed otherwise the communication does not happen
        call elph_mels_windows_free(self)
        call elph_mels_windows_init(self)
!call wave_window_flush_local_all(self%wavwin_cw)
!call wave_window_flush_local_all(self%wavwin_cproj)
        

    end subroutine elph_mels_cache_prefill_ph

!> @brief Intiialize the accumulator object for the self-energy
    subroutine elph_accumulator_selfen_init_ph(self, nbands, nbands_sum, nkpoints, nspin, nw, enwin, delta, &
                                               ntemps, temps, band_start, band_stop, &
                                               natoms, eigenvalues, nb_totk, efermi, carrier_per_cell, carrier_per_cell0)
        use constant
        use ini, only: register_allocate
        use elphon_base, only: celtot_max
        type(elph_accumulator_selfen_ph) :: self
        integer, intent(in) :: nw
        integer, intent(in) :: nbands
        integer, intent(in) :: nbands_sum
        integer, intent(in) :: nkpoints
        integer, intent(in) :: nspin
        integer, intent(in) :: natoms
        real(q), intent(in) :: enwin
        real(q), intent(in) :: delta
        real(q), intent(in) :: efermi(ntemps)
        real(q), intent(in) :: carrier_per_cell(ntemps) !< number of carriers per cell (with doping)
        real(q), intent(in) :: carrier_per_cell0        !< number of carriers per cell (without doping)
        integer, intent(in) :: ntemps
        integer, intent(in) :: band_start
        integer, intent(in) :: band_stop
!> temperatures at which the self-energy is to be evaluated (in K)
        real(q), intent(in) :: temps(ntemps)
        complex(q), intent(in) :: eigenvalues(:, :, :)
        integer, intent(in) :: nb_totk(:, :)
! Local variables
        integer :: spindeg, nelect

        self%nbands = nbands
        self%nbands_sum = nbands; if (nbands_sum > 0 .and. nbands_sum < nbands) self%nbands_sum = nbands_sum
        self%tetrahedron = delta <= 0
        self%kcut24_initialized = .false.
        self%nkpoints = nkpoints
        self%nspin = nspin
        if (nw > 0) then
            self%nw = nw + mod(nw + 1, 2)
            self%energy_grid_mode = 0
        else
            self%energy_grid_mode = 1
            self%nw = max(2, abs(nw))
        end if
        self%emax = celtot_max(eigenvalues, nb_totk) + 5.0_q
        self%emin = minval(real(eigenvalues, q)) - 5.0_q
        self%delta = abs(delta)
        self%enwin = enwin
        self%band_start = 1; if (band_start > 0 .and. band_start < nbands) self%band_start = band_start
        self%band_stop = nbands; if (band_stop > 0 .and. band_stop < nbands) self%band_stop = band_stop
        self%ntemps = ntemps
        allocate (self%tempsk(self%ntemps))
        self%tempsk = temps
        self%natoms = natoms
! logicals wether to perform fan and dw calculation
        allocate (self%t(self%ntemps, 3, self%natoms, 3, self%natoms, 3*self%natoms))
        self%t = cmplx(0.0_q, 0.0_q, q)
        allocate (self%efermi(self%ntemps))
        self%efermi = efermi
        allocate (self%carrier_per_cell(self%ntemps))
        self%carrier_per_cell = carrier_per_cell
        self%carrier_per_cell0 = carrier_per_cell0

! compute CBM and VBM and detect wether we are dealing with a metal or not at 0K
        spindeg = 2; if (self%nspin == 2) spindeg = 1
        nelect = nint(self%carrier_per_cell0)
        self%cbm = minval(real(eigenvalues(nelect/spindeg + 1, :, :)))
        self%vbm = maxval(real(eigenvalues(nelect/spindeg, :, :)))
        self%ismetal = .false.; if (self%vbm > self%cbm) self%ismetal = .true.

! allocate compute_mask array with the information of whether a particular
! k-point should be computed
        allocate (self%compute_mask(self%band_stop, self%nkpoints, self%nspin))
        call register_allocate(0.125_q*STORAGE_SIZE(self%compute_mask)*SIZE(self%compute_mask,KIND=qi8), 'elph_accumulator_selfen')
        self%compute_mask = .false.
    end subroutine elph_accumulator_selfen_init_ph

!> Initialize some arrays needed for using the tetrahedron method
    subroutine elph_accumulator_selfen_init_tetrahedron_ph(self, kpoints)
        use mkpoints, only: kpoints_struct, get_tetdisp_24tet
        type(elph_accumulator_selfen_ph), intent(inout) :: self
        type(kpoints_struct), intent(in) :: kpoints
! local variables
        integer :: ic, itet

! get tetrahedron displacements
        call get_tetdisp_24tet(kpoints%b, self%kcut24)
! get [0,0,0] point of each kcut
        do itet = 1, 24
            do ic = 1, 4
                if (all(self%kcut24(:, ic, itet) == [0, 0, 0])) exit
            end do
            self%kcut24_origin(itet) = ic
        end do
! transform kcut
        do itet = 1, 24
            do ic = 1, 4
                self%kcut24(:, ic, itet) = matmul(kpoints%ilsnf, self%kcut24(:, ic, itet))
            end do
        end do
        self%kcut24_initialized = .true.
    end subroutine elph_accumulator_selfen_init_tetrahedron_ph

!> Choose the k-points to compute explicitely
    subroutine elph_accumulator_selfen_selector_ikpt_ph(self, ikpts)
        type(elph_accumulator_selfen_ph) :: self
        integer :: ikpts(:)
! local variables
        integer :: ikpt, ikibz, isp
! loop over ikpts and set compute_mask
        do isp = 1, self%nspin
            do ikpt = 1, size(ikpts)
                ikibz = ikpts(ikpt)
! Here I set both spin components
! TODO: setting each spin channel separately is not implemented for now
                self%compute_mask(:, ikibz, isp) = .true.
            end do
        end do
    end subroutine elph_accumulator_selfen_selector_ikpt_ph

!> Select k-points based on an energy window
!> on output compute_mask is true for the states for which we want to compute the self-energy
    subroutine elph_accumulator_selfen_selector_energy_window_ph(self, w, select_energy_window)
        use wave_struct_def, only: wavespin
        use wave, only: ismetal
        type(elph_accumulator_selfen_ph) :: self
        type(wavespin), intent(in) :: w
        real(q), intent(in) :: select_energy_window(2)
! local variables
        integer :: ispin, ikibz, ibibz
        real(q) :: enkibz
        self%select_energy_window = select_energy_window

! find vbm and cbm
        self%compute_mask = .false.
        do ispin = 1, self%nspin
            do ikibz = 1, self%nkpoints
                do ibibz = 1, self%band_stop
                    enkibz = real(w%celtot(ibibz, ikibz, ispin), q)
                    self%compute_mask(ibibz, ikibz, ispin) = &
                        self%compute_mask(ibibz, ikibz, ispin) .or. &
                        (enkibz < self%cbm + self%select_energy_window(2) .and. &
                         enkibz > self%vbm - self%select_energy_window(1))
                end do ! bands
            end do ! kpoints
        end do ! spin
    end subroutine elph_accumulator_selfen_selector_energy_window_ph

!> Find the direct and indirect gaps and select those k-points to be computed
    subroutine elph_accumulator_selfen_selector_gaps_ph(self, w)
        use wave_struct_def, only: wavespin
        use wave, only: ismetal
        use elphon_accumulators, only: el_tol
        type(elph_accumulator_selfen_ph) :: self
        type(wavespin) :: w
! local variables
        integer :: ib, ispin, ikibz, nelect, spindeg
        integer :: ikpt_cbm, ikpt_vbm
        real(q) :: cbmk, vbmk
        real(q) :: e, direct_gap, gap

! find direct gaps
        self%compute_mask = .false.
        nelect = nint(self%carrier_per_cell0)
        spindeg = 2; if (self%nspin == 2) spindeg = 1
        do ispin = 1, self%nspin
            direct_gap = 1e8
            do ikibz = 1, self%nkpoints
                vbmk = real(w%celtot(nelect/spindeg, ikibz, ispin))
                cbmk = real(w%celtot(nelect/spindeg + 1, ikibz, ispin))
                gap = cbmk - vbmk
! store this point if the gap is the lowest
                if (direct_gap > gap) then
                    direct_gap = gap
                end if
            end do
! now we know what the value of the direct gap is we find if there is more than (1._q,0._q) point with this value
            do ikibz = 1, self%nkpoints
                vbmk = real(w%celtot(nelect/spindeg, ikibz, ispin))
                cbmk = real(w%celtot(nelect/spindeg + 1, ikibz, ispin))
                gap = cbmk - vbmk
                if (abs(direct_gap - gap) < el_tol) then
! select degenerate valence bands
                    do ib = 1, nelect/spindeg
                        e = real(w%celtot(ib, ikibz, ispin))
                        if (abs(e - vbmk) < el_tol) self%compute_mask(ib, ikibz, ispin) = .true.
                    end do
! select degenerate conduction bands
                    do ib = nelect/spindeg + 1, size(self%compute_mask, 1)
                        e = real(w%celtot(ib, ikibz, ispin))
                        if (abs(e - cbmk) < el_tol) self%compute_mask(ib, ikibz, ispin) = .true.
                    end do
                end if
            end do
        end do

! find indirect gaps
        do ispin = 1, self%nspin
! find k-points where vbm and cbm are
            do ikibz = 1, self%nkpoints
                vbmk = real(w%celtot(nelect/spindeg, ikibz, ispin))
                cbmk = real(w%celtot(nelect/spindeg + 1, ikibz, ispin))
                if (abs(vbmk - self%vbm) < el_tol) ikpt_vbm = ikibz
                if (abs(cbmk - self%cbm) < el_tol) ikpt_cbm = ikibz
            end do
! select degenerate valence bands
            do ib = 1, nelect/spindeg
                e = real(w%celtot(ib, ikpt_vbm, ispin))
                if (abs(e - self%vbm) < el_tol) self%compute_mask(ib, ikpt_vbm, ispin) = .true.
            end do
! select degenerate condunction bands
            do ib = nelect/spindeg + 1, self%nbands
                e = real(w%celtot(ib, ikpt_cbm, ispin))
                if (abs(e - self%cbm) < el_tol) self%compute_mask(ib, ikpt_cbm, ispin) = .true.
            end do
        end do

! set band_start and band_stop
        self%band_start = 1
        self%band_stop = self%nbands
        ispin = 1
! set band_start
        do ib = 1, size(self%compute_mask, 1)
            if (any(self%compute_mask(ib, :, ispin))) exit
        end do
        self%band_start = ib
! set band_stop
        do ib = size(self%compute_mask, 1), 1, -1
            if (any(self%compute_mask(ib, :, ispin))) exit
        end do
        self%band_stop = ib
    end subroutine elph_accumulator_selfen_selector_gaps_ph

!> Use information of compute_mask to only store the states for which
!> the FAN and DW self-energy are computed
    subroutine elph_accumulator_selfen_sparsify_ph(self, eigenvalues)
        use ini, only: register_allocate
        type(elph_accumulator_selfen_ph) :: self
        complex(q), intent(in) :: eigenvalues(:, :, :)
! local variables
        integer :: isp, ikibz, iband, ibibz, ispin, ibks
        integer :: nbks, iw
        real(q) :: enkibz

! set index of each element
        nbks = 0
        allocate (self%bks_idx(self%band_start:self%band_stop, self%nkpoints, self%nspin))
        self%bks_idx = 0
        do isp = 1, self%nspin
            do ikibz = 1, self%nkpoints
                do iband = self%band_start, self%band_stop
                    if (.not. self%compute_mask(iband, ikibz, isp)) cycle
                    nbks = nbks + 1
                    self%bks_idx(iband, ikibz, isp) = nbks
                end do
            end do
        end do

! create energy grids
        self%nbks = nbks
        allocate (self%energies(self%nw, self%nbks))
        do ispin = 1, self%nspin
            do ikibz = 1, self%nkpoints
                do ibibz = self%band_start, self%band_stop
                    enkibz = real(eigenvalues(ibibz, ikibz, ispin), q)
                    ibks = self%bks_idx(ibibz, ikibz, ispin)
                    if (ibks < 1) cycle
                    do iw = 1, self%nw
                        if (self%energy_grid_mode == 0) then
                            self%energies(iw, ibks) = &
                                enkibz + ((real(iw, q) - 0.5_q)/self%nw - 0.5_q)*self%enwin
                        else
                            self%energies(iw, ibks) = &
                                self%emin + (self%emax - self%emin)/(self%nw - 1)*(iw - 1)
                        end if
                    end do ! frequencies
                end do ! bands
            end do ! kpoints
        end do ! spin

! allocate arrays
        allocate (self%selfen_ph(self%ntemps, self%nw, self%nkpoints, self%natoms*3))
        call register_allocate(0.125_q*STORAGE_SIZE(self%selfen_ph)*SIZE(self%selfen_ph,KIND=qi8) + &
                               0.125_q*STORAGE_SIZE(self%bks_idx)*SIZE(self%bks_idx,KIND=qi8), 'elph_accumulator_selfen_ph')
        self%selfen_ph(:, :, :, :) = cmplx(0.0_q, 0.0_q, q)
    end subroutine elph_accumulator_selfen_sparsify_ph

!> @brief Write information about the computation of the self-enrgy
    subroutine elph_accumulator_selfen_print_ph(self, naccumulator, iu)
        use string, only: str
        type(elph_accumulator_selfen_ph), intent(in) :: self
        integer, intent(in) :: naccumulator !< index of the accumulator
        integer, intent(in) :: iu
! local variables
        integer :: ikpt, nkpoints_compute, ib, itemp
        logical, allocatable :: compute_band(:)
        if (iu < 0) return
        write (iu, '(A,I3)') 'Self-energy accumulator N =', naccumulator
        write (iu, '(A)') '----------------------------------'
        write (iu, '(A,I8,A,I8,A)') ' Band range:               [', self%band_start, ':', self%band_stop, ']'
        write (iu, '(A,I6)') ' Number of bands to sum over:', self%nbands_sum
        write (iu, '(A,F6.3)') ' Valence band maximum:   ', self%vbm
        write (iu, '(A,F6.3)') ' Conduction band minimum:', self%cbm
        write (iu, '(A,F6.3)') ' Complex imaginary shift (delta):', self%delta
        write (iu, '(A)') ' Chemical potentials mu(T):'
        write (iu, '(2A20)') 'Temperature (K)', 'mu (eV)'
        do itemp = 1, self%ntemps
            write (iu, '(2F20.8)') self%tempsk(itemp), self%efermi(itemp)
        end do
        if (.not. allocated(self%compute_mask)) return
        nkpoints_compute = 0
        do ikpt = 1, self%nkpoints
            if (any(self%compute_mask(:, ikpt, 1))) then
                nkpoints_compute = nkpoints_compute + 1
            end if
        end do
        write (iu, '(A,I6,A,I6,A)') ' Selected k-points at which to compute the self-energy: [', &
            nkpoints_compute, ' / ', self%nkpoints, ']'
        write (iu, '(A)') ' Selected bands and k-points for computation:'

        write (iu, '(A)', advance='no') '      '
        allocate (compute_band(lbound(self%compute_mask, 1):ubound(self%compute_mask, 1)))
        compute_band = .false.
        do ib = lbound(self%compute_mask, 1), ubound(self%compute_mask, 1)
            if (.not. any(self%compute_mask(ib, :, 1))) cycle
            compute_band(ib) = .true.
            write (iu, '(I4)', advance='no') ib
        end do
        write (iu, *)

        do ikpt = 1, self%nkpoints
            if (.not. any(self%compute_mask(:, ikpt, 1))) cycle
            write (iu, '(I6)', advance='no') ikpt
            do ib = lbound(self%compute_mask, 1), ubound(self%compute_mask, 1)
                if (.not. compute_band(ib)) cycle
                write (iu, '(L4)', advance='no') self%compute_mask(ib, ikpt, 1)
            end do
            write (iu, *)
        end do
        write (iu, *)
    end subroutine elph_accumulator_selfen_print_ph

!> @brief Sum results across the nodes
    subroutine elph_accumulator_selfen_sum_ph(self, comm_batch)
        use mpimy
        type(elph_accumulator_selfen_ph) :: self
        type(communic) :: comm_batch
        
        CALL m_sum_z(comm_batch, self%selfen_ph, size(self%selfen_ph))
        
    end subroutine elph_accumulator_selfen_sum_ph


!> @brief Write the information of the phonon self-energy accumulator to a file
    subroutine elph_accumulator_selfen_hdf5write_ph(self, phonon_freqs_ibz, kpoints, fileid, idx, group, subgroup)
        use vhdf5
        use string, only: str
        use mkpoints_struct_def, only: kpoints_struct
        implicit none

        type(elph_accumulator_selfen_ph) :: self
        real(q), intent(in) :: phonon_freqs_ibz(:, :)
        type(kpoints_struct), intent(in) :: kpoints
        integer, intent(in) :: idx
        integer(HID_T), intent(in) :: fileid
        character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
        character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_DOS)
! local variables
        integer(HID_T) :: groupid, subgroupid
        character(len=MAX_LEN_GROUP) :: my_group
        character(len=MAX_LEN_GROUP) :: my_subgroup

        my_group = GRP_RESULTS; if (present(group)) my_group = trim(group)
        my_subgroup = GRP_ELPHON; if (present(subgroup)) my_subgroup = trim(subgroup)

! open group
        call vh5_error(vh5_group_open_or_create(fileid, trim(my_group)//'/'//trim(my_subgroup)//'/phonons', groupid),"elphon_selfen_ph.F",758)
        call vh5_error(vh5_group_open_or_create(groupid, 'self_energy_'//str(idx), subgroupid),"elphon_selfen_ph.F",759)

! write important metadata
        call vh5_error(vh5_write(subgroupid, 'nbands', self%nbands),"elphon_selfen_ph.F",762)
        call vh5_error(vh5_write(subgroupid, 'nbands_sum', self%nbands_sum),"elphon_selfen_ph.F",763)
        call vh5_error(vh5_write(subgroupid, 'delta', self%delta),"elphon_selfen_ph.F",764)
        call vh5_error(vh5_write(subgroupid, 'band_start', self%band_start),"elphon_selfen_ph.F",765)
        call vh5_error(vh5_write(subgroupid, 'band_stop', self%band_stop),"elphon_selfen_ph.F",766)
        call vh5_error(vh5_write(subgroupid, 'nw', self%nw),"elphon_selfen_ph.F",767)
        call vh5_error(vh5_write(subgroupid, 'enwin', self%enwin),"elphon_selfen_ph.F",768)
        call vh5_error(vh5_write(subgroupid, 'efermi', self%efermi),"elphon_selfen_ph.F",769)
        call vh5_error(vh5_write(subgroupid, 'tetrahedron', self%tetrahedron),"elphon_selfen_ph.F",770)
        call vh5_error(vh5_write(subgroupid, 'select_energy_window', self%select_energy_window),"elphon_selfen_ph.F",771)
        call vh5_error(vh5_write(subgroupid, 'bks_idx', self%bks_idx),"elphon_selfen_ph.F",772)

! write main data
        call vh5_error(vh5_write(subgroupid, 'selfen_ph', self%selfen_ph),"elphon_selfen_ph.F",775)
        call vh5_error(vh5_write(subgroupid, 'phonon_freqs_ibz', phonon_freqs_ibz),"elphon_selfen_ph.F",776)
        call vh5_error(vh5_write(subgroupid, "kpoint_generating_vectors", kpoints%B),"elphon_selfen_ph.F",777)
        call vh5_error(vh5_write(subgroupid, "kpoint_coords", kpoints%vkpt),"elphon_selfen_ph.F",778)
        call vh5_error(vh5_write(subgroupid, "kpoints_symmetry_weight", kpoints%wtkpt),"elphon_selfen_ph.F",779)
        call vh5_error(vh5_write(subgroupid, 'temps', self%tempsk),"elphon_selfen_ph.F",780)
        call vh5_error(vh5_write(subgroupid, 'energies', self%energies),"elphon_selfen_ph.F",781)

! close group
        call vh5_error(vh5_group_close_writing(subgroupid),"elphon_selfen_ph.F",784)
        call vh5_error(vh5_group_close_writing(groupid),"elphon_selfen_ph.F",785)
    end subroutine elph_accumulator_selfen_hdf5write_ph


!> @brief Free the accumulator structure
    subroutine elph_accumulator_selfen_free_ph(self)
        use ini, only: deregister_allocate
        type(elph_accumulator_selfen_ph) :: self
        deallocate (self%tempsk)
        deallocate (self%energies)
        deallocate (self%efermi)
        call deregister_allocate(0.125_q*STORAGE_SIZE(self%selfen_ph)*SIZE(self%selfen_ph,KIND=qi8) + &
                                 0.125_q*STORAGE_SIZE(self%bks_idx)*SIZE(self%bks_idx,KIND=qi8) + &
                                 0.125_q*STORAGE_SIZE(self%compute_mask)*SIZE(self%compute_mask,KIND=qi8), 'elph_accumulator_selfen_ph')
        deallocate (self%selfen_ph)
        deallocate (self%compute_mask)
        deallocate (self%bks_idx)
    end subroutine elph_accumulator_selfen_free_ph

! from here on its about routines that deal with a set of accumulators
    subroutine elph_accumulators_init_ph(self, elph_set, wdes)
        use elphon_base, only: elph_settings
        use wave_struct_def, only: wavedes
        type(elph_accumulators_selfen_ph) :: self
        type(elph_settings), intent(in) :: elph_set
        type(wavedes), intent(in) :: wdes
! local variables
        integer :: nb_global, nb
        self%naccumulators_ph = elph_set%selfen_naccumulators
        allocate (self%selfen_ph(self%naccumulators_ph))
! keep track of the index of the bands
        allocate (self%band_idx(wdes%nbands))
        do nb_global = 1, wdes%nb_tot
            if (mod(nb_global - 1, wdes%nb_par) + 1 == wdes%nb_low) then
                nb = 1 + (nb_global - 1)/wdes%nb_par
                self%band_idx(nb) = nb_global
            end if
        end do ! loop over bands
    end subroutine elph_accumulators_init_ph

    subroutine elph_accumulators_free_ph(self)
        type(elph_accumulators_selfen_ph) :: self
! local variables
        integer :: i
        do i = 1, self%naccumulators_ph
            call elph_accumulator_selfen_free_ph(self%selfen_ph(i))
        end do
    end subroutine elph_accumulators_free_ph

!> Finalize calculation of the self-energy
    subroutine elph_accumulators_finalize_ph(self, comm_global)
        use mpimy, only: communic
        type(elph_accumulators_selfen_ph) :: self
        type(communic), intent(in) :: comm_global
! local variables
        integer :: i
        do i = 1, self%naccumulators_ph
            call elph_accumulator_selfen_sum_ph(self%selfen_ph(i), comm_global)
        end do
    end subroutine elph_accumulators_finalize_ph

!> Write results to hdf5 file
    subroutine elph_accumulators_print_ph(self, w_dense_ibz, phonon_freqs_ibz, kpoints, io)
        use base, only: in_struct
        use wave_struct_def, only: wavespin
        use mkpoints_struct_def, only: kpoints_struct

        USE vhdf5, ONLY: IH5OUTFILEID

        implicit none

        type(elph_accumulators_selfen_ph) :: self
        type(wavespin), intent(in) :: w_dense_ibz
        type(kpoints_struct), intent(in) :: kpoints
        real(q), intent(in) :: phonon_freqs_ibz(:, :)
        type(in_struct), intent(in) :: io
! local variables
        integer :: i
        do i = 1, self%naccumulators_ph

            call elph_accumulator_selfen_hdf5write_ph(self%selfen_ph(i), phonon_freqs_ibz, kpoints, ih5outfileid, i)

            call elph_accumulator_selfen_print_result_ph(self%selfen_ph(i), i, w_dense_ibz%celtot, io%iu6)
        end do
    end subroutine elph_accumulators_print_ph

!> Write the results from the calculation of the self-energy to the outcar
    subroutine elph_accumulator_selfen_print_result_ph(self, i, celtot, iu)
        use string, only: str
        type(elph_accumulator_selfen_ph) :: self
        integer, intent(in) :: i
        complex(q), intent(in) :: celtot(:, :, :) !< eigenvalues in the IBZ
        integer, intent(in) :: iu
! local variables
        integer :: isp, ikpt, ib
        integer :: ibks, itemp
        real(q) :: reph, imph

        if (iu < 0) return

        WRITE (iu, '(A,I3)') 'Phonon self-energy accumulator N=', i
        if (self%nw > 1) then
            write (iu, *) 'nw > 1 so the self-energy is not written'
            return
        end if

        do itemp = 1, self%ntemps
            write (iu, '(A,F8.0,A)') ' T= ', self%tempsk(itemp), ' K'
            write (iu, '(3A6,3A12)') 'ispin', 'iqpt', 'iband', 'KS eV', 're(ph) eV', 'im(ph) eV'
            do isp = 1, self%nspin
                do ikpt = 1, self%nkpoints
                    if (.not. any(self%compute_mask(:, ikpt, isp))) cycle
                    do ib = 1, self%natoms*3
                        ibks = self%bks_idx(ib, ikpt, isp)
                        if (ibks < 1) cycle
                        reph = real(self%selfen_ph(itemp, 1, ikpt, ib), q)
                        imph = aimag(self%selfen_ph(itemp, 1, ikpt, ib))
                        write (iu, '(3I6,5F12.6)') isp, ikpt, ib, &
                            real(celtot(ib, ikpt, isp), q), &
                            reph, imph
                    end do
                end do !kpoints
            end do !spin
            write (iu, *)
        end do !itemp
        write (iu, *)
    end subroutine elph_accumulator_selfen_print_result_ph

    subroutine electron_phonon_selfen_ph_driver(elph_set, elph_mels, elph_accumulators_ph, kpoints)
        use elphon_base, only: elph_settings
        use elphon_mels, only: elph_mels_type, elph_mels_compute_g
        use triplets, only: triplet_type, kpoint_init_kindex, kpoint_init_kint
        use mkpoints, only: vkpt_kindex, kpoints_struct
        implicit none

        type(elph_settings), intent(in) :: elph_set
        type(elph_mels_type), intent(in) :: elph_mels
        type(elph_accumulators_selfen_ph), intent(inout) :: elph_accumulators_ph
        type(kpoints_struct), intent(in) :: kpoints
! local variables
        type(triplet_type) :: triplet
        integer :: isp
        integer :: ikp_fbz, iq_fbz
        integer :: ik_ibz
        integer :: ikp
        integer :: ibz2kindex(kpoints%nkpts)
        complex(q), allocatable :: gcart_kkp(:, :, :, :)
        complex(q), allocatable :: gmode_kkp(:, :, :)

        if (elph_set%dry_run) return

! create mapping from the ibz to kindex
        do ik_ibz = 1, kpoints%nkpts
            ibz2kindex(ik_ibz) = vkpt_kindex(kpoints, elph_mels%latt_cur, kpoints%vkpt(:, ik_ibz))
        end do

        allocate (gcart_kkp(elph_mels%nbands_k, elph_mels%nbands, 3, elph_mels%natoms_pc))
        allocate (gmode_kkp(elph_mels%nbands_k, elph_mels%nbands, 3*elph_mels%natoms_pc))

        

        do isp = 1, elph_mels%ispin
! loop over the q points in the fbz.
! Those of (q, k', ir-k) are selected.
            do iq_fbz=1,kpoints%nkfull
                if (elph_mels%io%iu6>=0) then
                    if (iq_fbz==1 .or. iq_fbz==kpoints%nkfull .or. mod(iq_fbz,100)==0) then
                        write(*,*) 'ibz qpoint [',iq_fbz,'/',kpoints%nkfull,']'
                    endif
                endif
                call kpoint_init_kindex(triplet%q, kpoints, iq_fbz)

! loop over k' in each distribution bin
                do ikp = 1, size(elph_mels%ikp_fbz_list)
                    ikp_fbz = elph_mels%ikp_fbz_list(ikp)
                    call kpoint_init_kindex(triplet%kp, kpoints, ikp_fbz)
                    call kpoint_init_kint(triplet%k, kpoints, triplet%q%kint + triplet%kp%kint)

! Pass through only k in irreducible k-points.
                    if (triplet%k%kindex/=ibz2kindex(triplet%k%ibz)) cycle

! contributes = elph_accumulators_contributes_ph(elph_accumulators_ph, elph_set, &
!                                                elph_mels%w_dense_ibz, elph_mels%phonon_freqs_ibz, &
!                                                elph_mels%latt_cur, kpoints, triplet, isp)
! if (.not. contributes) cycle

! compute the electron-phonon matrix elements
                    call elph_mels_compute_g(elph_mels,kpoints,triplet%k%kindex,ikp_fbz,iq_fbz,isp,gcart_kkp,gmode_kkp)

! use the electron-phonon matrix elements to compute the electron self-energy
                    call elph_accumulators_accumulate_ph( &
                        elph_accumulators_ph, elph_mels%band_start_k, &
                        triplet, isp, kpoints, &
                        elph_mels%w_dense_ibz, elph_mels%wdes_dense_ibz, &
                        elph_mels%phonon_freqs_ibz, gmode_kkp)
                end do !loop over k'
            end do !loop over q
        end do ! loop over spin
        

    end subroutine electron_phonon_selfen_ph_driver

!> checck if the k, q and kp scattering event contributes to the any of the self-energy accumulators
    function elph_accumulators_contributes_ph(self, elph_set, w_dense_ibz, phonon_freqs_ibz, kpoints, &
                                              triplet, isp) result(contributes)
        use elphon_base, only: elph_settings
        use wave_struct_def, only: wavespin
        use mkpoints_struct_def, only: kpoints_struct
        use mkpoints, only: kindex_vkpt, ikibz_kindex
        use triplets, only: triplet_type
        type(elph_accumulators_selfen_ph) :: self
        type(elph_settings) :: elph_set
        type(wavespin) :: w_dense_ibz
        real(q) :: phonon_freqs_ibz(:, :)
        type(kpoints_struct) :: kpoints
        type(triplet_type) :: triplet
        integer :: isp
        logical :: contributes

        
        contributes = .true.
        if (.not. elph_accumulators_shallido_ph(self, triplet%k%ibz, isp)) then
            contributes = .false.
            
            return
        end if
        if (elph_accumulators_delta_is_zero_ph(self, elph_set, w_dense_ibz, phonon_freqs_ibz, &
                                               kpoints, triplet%k%ibz, isp, triplet%kp%kint, &
                                               triplet%q%kint)) then
            contributes = .false.
            
            return
        end if
        
    end function elph_accumulators_contributes_ph

! Check if any of the delta functions for the accumulators is not (0._q,0._q)
    function elph_accumulators_delta_is_zero_ph(self, elph_set, w_dense_ibz, phonon_freqs_ibz, &
                                                kpoints_inter, ik_ibz, isp, kpint, qint) &
        result(delta_is_zero)
        use elphon_base, only: elph_settings
        use wave_struct_def, only: wavespin
        use mkpoints_struct_def, only: kpoints_struct
        type(elph_accumulators_selfen_ph), intent(in) :: self
        type(elph_settings), intent(in) :: elph_set
        type(wavespin), intent(in) :: w_dense_ibz
        real(q), intent(in) :: phonon_freqs_ibz(:, :)
        type(kpoints_struct), intent(in) :: kpoints_inter
        integer, intent(in) :: ik_ibz
        integer, intent(in) :: isp
        integer, intent(in) :: kpint(3)
        integer, intent(in) :: qint(3)
        logical :: delta_is_zero
! local variables
        integer :: i
        delta_is_zero = .false.
! next check is specific for transport computations
        if (elph_set%selfen_imag_skip) then
        do i = 1, self%naccumulators_ph
        if (elph_accumulator_selfen_tet_contributes_low_ph(self%selfen_ph(i), kpoints_inter, &
                                                           ik_ibz, isp, kpint, qint, w_dense_ibz%celen(:, :, isp), &
                                                           phonon_freqs_ibz)) then
            return
        end if
        end do
        delta_is_zero = .true.
        end if
    end function elph_accumulators_delta_is_zero_ph

!> Check wether this q-point contributes to the imaginary part of the self-energy
!> of this k-point
    function elph_accumulator_selfen_tet_contributes_low_ph(self, kpoints, ikibz, ispin, kpint, qint, &
                                                            eleig_ibz, pheig_ibz) result(contributes)
        use constant
        use mkpoints, only: kpoints_struct, kpoint_integers, get_indexes, kint_mod, kint_kindex
        use elphon_accumulators, only: tetrahedron_kcontributes, get_energies, get_energies_complex, thztoev
        type(elph_accumulator_selfen_ph), intent(in) :: self
        type(kpoints_struct), intent(in) :: kpoints !< kpoints structure
        integer, intent(in) :: ikibz !< index of the k in the ibz
        integer, intent(in) :: ispin !< spin index
        integer, intent(in) :: kpint(3) !< coordinates of k'
        integer, intent(in) :: qint(3) !< coordinates of q
        complex(q), intent(in) :: eleig_ibz(:, :) !< electron eigenvalues in the IBZ
        real(q), intent(in) :: pheig_ibz(:, :) !< phonon eigenvalues in the IBZ
        logical contributes
! local variables
        integer :: imode, ibibz, ibfbz, ibks
        integer :: idx_kp(4, 24)
        integer :: idx_q(4, 24)
        real(q) :: eleig_tetra(4, 24)
        real(q) :: pheig_tetra(4, 24, self%natoms*3)
        real(q) :: w(self%nw)

        

! find the index of the k' point
        call get_indexes(kpoints, kpint, self%kcut24, idx_kp)

! find the index of the q  point
        call get_indexes(kpoints, qint, self%kcut24, idx_q)

! loop over phonon modes
        do imode = 1, self%natoms*3
! get the phonon frequencies on the 24 tetrahedron for q
            call get_energies(idx_q, imode, pheig_ibz, pheig_tetra(:, :, imode))
        end do
        pheig_tetra = pheig_tetra*thztoev

! loop over bands at k'
        do ibfbz = 1, size(eleig_ibz, 1) !self%nbands_sum
! get the electron energies on the 24 corners for k'
            call get_energies_complex(idx_kp, ibfbz, eleig_ibz, eleig_tetra)

! loop over bands at k
            do ibibz = self%band_start, self%band_stop
                ibks = self%bks_idx(ibibz, ikibz, ispin)
                if (ibks < 1) cycle
                w = self%energies(:, ibks)
                do imode = 1, self%natoms*3
! check wether the electron_energies+phonon_energies contribute
                    if (tetrahedron_kcontributes(eleig_tetra + pheig_tetra(:, :, imode), w)) then
                        contributes = .true.
                        
                        return
                    end if
                end do !phonon modes
            end do !bands at k
        end do !bands at k'
        contributes = .false.
        
    end function elph_accumulator_selfen_tet_contributes_low_ph

!> return fractional coordinates of kindex
!> The idea behind this function is slightly different from kindex_kint function.
    function get_kpoint_coordinates_from_ikfbz(kpoints, kindex, igrpop) result(vk)
        use mkpoints, only: kpoints_fbz_to_ibz
        use wave_rotate, only: rotate_vkpt_irot
        use mkpoints_struct_def, only: kpoints_struct
        implicit none

        type(kpoints_struct), target, intent(in) :: kpoints
        integer, intent(in) :: kindex !< index of the k-point in the fbz
        integer, intent(in) :: igrpop(3, 3, 48) !< reciprocal space symmetry operations

        integer :: ikibz !< index of the k-point in the ibz
        integer :: irot !< index of the rotation that expands the k-point from the ibz to the fbz
        real(q) :: vk(3)

        call kpoints_fbz_to_ibz(kpoints, kindex, ikibz, irot)
        vk = rotate_vkpt_irot(kpoints%vkpt(:, ikibz), igrpop, irot)
    end function get_kpoint_coordinates_from_ikfbz

    function elph_accumulators_shallido_ph(self, ik_ibz, isp) result(shallido)
        type(elph_accumulators_selfen_ph) :: self
        integer, intent(in) :: ik_ibz
        integer, intent(in) :: isp
        logical :: shallido
! local variables
        integer :: i
        shallido = .false.
! first check is wether this k-point is marked for computation
        do i = 1, self%naccumulators_ph
            if (elph_accumulator_selfen_shallido_ph(self%selfen_ph(i), ik_ibz, isp)) then
                shallido = .true.
                return
            end if
        end do
    end function elph_accumulators_shallido_ph

!> @brief Determine whether we should compute this k and k' combination
    pure logical function elph_accumulator_selfen_shallido_ph(self, ikibz, ispin) result(shallido)
        type(elph_accumulator_selfen_ph), intent(in) :: self
        integer, intent(in) :: ikibz
        integer, intent(in) :: ispin
        shallido = any(self%compute_mask(:, ikibz, ispin))
    end function elph_accumulator_selfen_shallido_ph

    subroutine elph_accumulators_accumulate_ph(self, band_start_k, triplet, isp, &
                                               kpoints, w_dense_ibz, wdes_dense_ibz, &
                                               phonon_freqs_ibz, gmode_kkp)
        use mkpoints_struct_def, only: kpoints_struct
        use wave_struct_def, only: wavespin, wavedes
        use wave_rotate, only: wave_rotator
        use triplets, only: triplet_type
        implicit none

        type(elph_accumulators_selfen_ph), intent(inout) :: self
        integer, intent(in) :: band_start_k
        type(triplet_type), intent(in) :: triplet
        integer, intent(in) :: isp
        type(kpoints_struct), intent(in) :: kpoints
        type(wavespin), intent(in) :: w_dense_ibz
        type(wavedes), intent(in) :: wdes_dense_ibz
        real(q), intent(in) :: phonon_freqs_ibz(:, :)
        complex(q), intent(in) :: gmode_kkp(:, :, :)
! local variables
        integer :: i

        
        do i = 1, self%naccumulators_ph
            call elph_accumulator_selfen_tet_only_ph(self%selfen_ph(i), kpoints, isp, triplet, &
                                                     w_dense_ibz%celen(:, :, isp), w_dense_ibz%celtot(:, :, isp), &
                                                     self%band_idx, phonon_freqs_ibz(:, triplet%q%ibz), &
                                                     band_start_k, gmode_kkp, wdes_dense_ibz%rspin, wdes_dense_ibz%comm_inb)
        end do
        
    end subroutine elph_accumulators_accumulate_ph

!> Accumulate the contribution to phonon self-energy.
!> This is a high lever wrapper for the lorentzian (not implemented) and tetrahedron method.
    subroutine elph_accumulator_selfen_tet_only_ph(self, kpoints, ispin, triplet, &
                                                   eleig_kp_celen, eleig_k_celtot, enfbz_idx, eph_thz, band_start_k, g, &
                                                   rspin, comm)
        use mkpoints, only: kpoints_struct
        use mpimy, only: communic
        use triplets, only: triplet_type
        implicit none

        type(elph_accumulator_selfen_ph), intent(inout) :: self
        type(kpoints_struct), intent(in) :: kpoints !< kpoints structure
        integer, intent(in) :: ispin !< spin index
        type(triplet_type), intent(in) :: triplet
        integer, intent(in) :: enfbz_idx(:) !< index of the bands
        real(q), intent(in) :: eph_thz(:) !< phonon frequencies in THz with 2pi of an ir-qpoint
        complex(q), intent(in) :: eleig_kp_celen(:, :) !< electron eigenvalues in the IBZ
        complex(q), intent(in) :: eleig_k_celtot(:, :) !< electron eigenvalues for all bands in the IBZ
        integer, intent(in) :: band_start_k !< first band index of k in the array with g matrix elements
        complex(q), intent(in) :: g(:, :, :) ! (band at k, band at k', phonon mode)
        real(q), intent(in) :: rspin

        type(communic) :: comm

        if (self%tetrahedron) then
            call elph_accumulator_selfen_tet_ph(self, kpoints, ispin, triplet, &
                                                eleig_kp_celen, eleig_k_celtot, enfbz_idx, eph_thz, band_start_k, &
                                                g, rspin, comm)
        end if
    end subroutine elph_accumulator_selfen_tet_only_ph

!> @brief Accumulate the contribution to the phonon self-energy using the tetrahedron method.
!> Contributions from q-points in FBZ are accumurated to those of corresponding ir-q-points.
!> This works because only ir-k-points come in this subroutine. The idea is:
!> (ir-k, kp, q) is mapped to (R.k, R.kp, R.q=ir-q)
!> implicitly under the invariance of g for rotation.
    subroutine elph_accumulator_selfen_tet_ph(self, kpoints, ispin, triplet, &
                                              eleig_kp_celen, eleig_k_celtot, enfbz_idx, eph_thz, &
                                              band_start_k, g, rspin, comm)
        use constant, only: bolkev
        use mkpoints, only: kpoints_struct, get_indexes
        use tet_macdonald, only: weight_24tetra_delta
        use elphon_base, only: fermi_dirac
        use elphon_accumulators, only: phon_tol, thztoev, get_energies_complex
        use mpimy, only: communic
        use triplets, only: triplet_type

        implicit none

        type(elph_accumulator_selfen_ph), intent(inout) :: self
        type(kpoints_struct), intent(in) :: kpoints !< kpoints structure
        integer, intent(in) :: ispin !< spin index
        type(triplet_type), intent(in) :: triplet
        integer, intent(in) :: enfbz_idx(:) !< index of the bands
        real(q), intent(in) :: eph_thz(:)
        complex(q), intent(in) :: eleig_kp_celen(:, :) !< electron eigenvalues in the IBZ
        complex(q), intent(in) :: eleig_k_celtot(:, :) !< electron eigenvalues for all bands in the IBZ
        integer, intent(in) :: band_start_k !< first band index of k in the array with g matrix elements
        complex(q), intent(in) :: g(:, :, :) ! (band at k, band at k', phonon mode)
        real(q), intent(in) :: rspin
        type(communic) :: comm !< communicator to paralelize the calculation

! local variables
        integer :: imode, iw, ibfbz, ibibz, ibks, itemp
        integer :: iqibz !< index of the q in the ibz
        integer :: ikibz !< index of the k in the ibz
        integer :: ikpibz !< index of the k' in the ibz
        integer :: idx_k(4, 24)
        integer :: idx_kp(4, 24)
        real(q) :: weight = 1.0
        real(q) :: eleig_tetra_k(4, 24)
        real(q) :: eleig_tetra_kp(4, 24)
        real(q) :: deltamw(self%natoms*3)
        real(q) :: fk(self%ntemps)
        real(q) :: fkp(self%ntemps)
        real(q) :: temps(self%ntemps)
        real(q) :: eph(self%natoms*3)
        complex(q) :: selfen_ph(self%ntemps, self%nw)

        
        eph = eph_thz*thztoev
        temps = self%tempsk*bolkev ! convert to eV

        ikibz = triplet%k%ibz
        ikpibz = triplet%kp%ibz
        iqibz = triplet%q%ibz

        call get_indexes(kpoints, triplet%k%kint, self%kcut24, idx_k)
        call get_indexes(kpoints, triplet%kp%kint, self%kcut24, idx_kp)

! loop over bands at k'
        do ibfbz = 1, size(eleig_kp_celen, 1)
! when ncore/=1 we can still use this comunicator to paralelize the self-energy computation
            if (mod(ibfbz, comm%ncpu) /= comm%node_me - 1) cycle
            if (enfbz_idx(ibfbz) > self%nbands_sum) cycle

            call get_energies_complex(idx_kp, ibfbz, eleig_kp_celen, eleig_tetra_kp)
            do itemp = 1, self%ntemps
                fkp(itemp) = fermi_dirac(real(eleig_kp_celen(ibfbz, ikpibz), q) - self%efermi(itemp), temps(itemp))
            end do

! loop over bands at k
            do ibibz = self%band_start, self%band_stop
                call get_energies_complex(idx_k, ibibz, eleig_k_celtot, eleig_tetra_k)
                do itemp = 1, self%ntemps
                    fk(itemp) = fermi_dirac(real(eleig_k_celtot(ibibz, ikibz), q) - self%efermi(itemp), temps(itemp))
                end do
                ibks = self%bks_idx(ibibz, ikibz, ispin)
                if (ibks < 1) cycle
                selfen_ph(:, :) = cmplx(0.0_q, 0.0_q, q)
                deltamw = weight_24tetra_delta(kpoints%volwgt, self%kcut24_origin, &
                                               self%natoms*3, eph, eleig_tetra_k - eleig_tetra_kp)

! loop over phonon modes
                do imode = 1, self%natoms*3
                    if (abs(eph(imode)) < phon_tol) cycle
! loop over energies (this loop will be removed since nw=1 for phonon selfenergy.)
                    do iw = 1, self%nw
! First dim of selfen_ph counts temperatures.
                        selfen_ph(:, iw) = selfen_ph(:, iw) + rspin*weight &
                                           *g(ibibz - band_start_k + 1, ibfbz, imode) &
                                           *conjg(g(ibibz - band_start_k + 1, ibfbz, imode)) &
                                           *(fkp - fk) *deltamw(imode) *cmplx(0.0_q, pi, q)
                    end do
                    self%selfen_ph(:, :, iqibz, imode) = self%selfen_ph(:, :, iqibz, imode) + selfen_ph(:, :)
                end do ! modes
            end do ! k bands
        end do ! k' bands
        
    end subroutine elph_accumulator_selfen_tet_ph

end module elphon_selfen_ph
