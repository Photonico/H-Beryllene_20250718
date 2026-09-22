# 1 "elphon_driver.F"
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


# 2 "elphon_driver.F" 2 

! Unified driver to compute electron and phonon self-energies due to electron-phonon coupling
! elphon_mels_sandwich routine should have a caching mechanism for the potential
! two consecutive calls with the same q-point skip the interpolation
module elphon_driver
 use prec
 use mkpoints_struct_def, only : kpoints_struct
 use elphon_base, only: elph_settings
 use elphon_mels, only: elph_mels_type, elph_mels_init, elph_mels_compute_wfs, &
                        elph_mels_pot_init, elph_mels_pot_read_hdf5, elph_mels_pot_setup, &
                        elph_mels_ifc_init
 use elphon_accumulators, only: elph_accumulators_selfen, elph_accumulators_g
 implicit none

 contains
 subroutine electron_phonon_driver(ELPH_SET, HAMILTONIAN, KPOINTS, GRID, LATT_CUR, LATT_INI, &
            T_INFO, NONLR_S, NONL_S, W, LMDIM, P, SV, CQIJ, CDIJ, SYMM, INFO, IO, &
            GRIDC, GRIDUS, C_TO_US, IRDMAX)
!! USE moffload_struct_def
    USE base, ONLY: info_struct, in_struct, symmetry
    USE hamil_struct_def, ONLY: ham_handle
    USE tutor, ONLY: vtutor
    USE string, ONLY: str
    USE mpimy, ONLY : communic
    USE mgrid_struct_def, ONLY: grid_3d, transit
    USE wave, ONLY: wdes_from_wdes, wdes_set_npro
    USE lattice, ONLY: latt
    USE nonl_struct_def, ONLY: nonl_struct
    USE nonlr_struct_def, ONLY: nonlr_struct
    USE poscar, ONLY: type_info
    USE pseudo, ONLY: potcar
    USE wave_struct_def, ONLY: wavespin, wavedes
    USE ini, ONLY: dump_allocate
    USE elphon_kgrid, ONLY: elph_kgrid
    USE elphon_base, ONLY: elph_driver_el, elph_driver_ph, elph_driver_mels, elph_header_write, elph_footer_write, &
                           elph_write_settings, elph_scattering_approx_mrta_tau, elph_scattering_approx_mrta_lambda, &
                           elph_scattering_approx_erta_tau, elph_scattering_approx_erta_lambda, elph_scattering_approx_crta
    USE elphon_mels, ONLY: elph_mels_ifc_init, elph_mels_pot_init, elph_mels_get_w_dense, elph_mels_ifc_init, &
                           elph_mels_ifc_read_hdf5, elph_mels_ifc_setup, elph_mels_ifc_compute_ibz_phonons, elph_mels_ifc_free, elph_mels_pot_free, &
                           elph_mels_wfs_free, elph_mels_free, elph_mels_read_dimensions, elph_mels_set_nbandsk, &
                           elph_mels_cache_init, elph_mels_cache_free, elph_mels_windows_init, elph_mels_wave_caches_init, elph_mels_windows_free, &
                           elph_mels_wave_caches_free, elph_mels_cache_prefill, elph_mels_cache_prefill_bcast, elph_mels_barrier, elph_mels_get_ikpfbz, &
                           elph_mels_wf_redistribute, elph_mels_ifc_compute_dispersion, elph_mels_pot_dispersion, elph_mels_wf_write
    USE elphon_common, ONLY: elphon_comm_init, elphon_get_kpoints_inter, elphon_prepare_selfenergy_accumulators, elphon_set_nbands, &
                             elphon_get_velocity, elphon_write_velocity, transport_elphon
    USE elphon_accumulators, ONLY: elph_accumulators_init, elph_accumulators_finalize, elph_accumulators_print, elph_accumulators_free
    use elphon_selfen_ph, only: elphon_prepare_selfenergy_accumulators_ph, elph_accumulators_init_ph, elph_accumulators_selfen_ph, &
                                electron_phonon_selfen_ph_driver, elph_mels_wf_redistribute_ph, elph_mels_cache_prefill_ph, &
                                elph_accumulators_finalize_ph, elph_accumulators_print_ph, elph_accumulators_free_ph
    use electron_transport, only: chemical_potential
    use bandgap_struct, only: bandgap_info
    TYPE (elph_settings)  :: ELPH_SET !< information for the electron-phonon calculation
    TYPE (ham_handle)     :: HAMILTONIAN
    TYPE (kpoints_struct) :: KPOINTS
    TYPE (grid_3d),TARGET :: GRID
    TYPE (latt)           :: LATT_CUR
    TYPE (latt)           :: LATT_INI
    TYPE (type_info)      :: T_INFO
    TYPE (nonlr_struct)   :: NONLR_S
    TYPE (nonl_struct)    :: NONL_S
    TYPE (potcar)         :: P(:)
    TYPE (grid_3d)        :: GRIDC
    TYPE (grid_3d)        :: GRIDUS
    TYPE (transit)        :: C_TO_US
    TYPE (info_struct)    :: INFO
    TYPE (in_struct)      :: IO
    TYPE (symmetry)       :: SYMM
    INTEGER               :: LMDIM, IRDMAX
    TYPE (wavespin), TARGET :: W
    COMPLEX(q) :: CQIJ (:,:,:,:)
    COMPLEX(q) :: CDIJ (:,:,:,:)
    COMPLEX(q)   :: SV(:,:)
! local variables
    type(communic) :: comm_kbatch !< Communicator for each batch of kpoints
    type(communic) :: comm_global !< Communicator for each batch of kpoints
    type(elph_mels_type) :: elph_mels
    type(kpoints_struct) :: kpoints_inter
    type(elph_accumulators_selfen) :: elph_accumulators !< accumulators for electron self-energy
    type(wavedes)  :: wdes_dense_ibz
    type(elph_kgrid) :: kgrid ! kspclib related variables
    logical :: qpoints_exists
! velocity matrix elements
    real(q), allocatable :: velocity(:,:,:,:)
    type(chemical_potential) :: chempot !< chemical potentials computed at different temperatures and doping levels
    type(bandgap_info),allocatable :: bandgap(:)

! local variables (phonon selfenergy)
    type(elph_accumulators_selfen_ph) :: elph_accumulators_ph

# 96


    call elph_header_write(io)
    call elph_write_settings(elph_set,io%iu6)

    
    call wdes_from_wdes(wdes_dense_ibz,w%wdes)
    comm_global = w%wdes%comm
    if (elph_set%ncore/=w%wdes%comm_inb%ncpu) then
       call vtutor%alert('The NCORE='//str(elph_set%ncore)//' setting in the INCAR was not satisfied, try to continue using NCORE='//str(w%wdes%comm_inb%ncpu))
       elph_set%ncore = w%wdes%comm_inb%ncpu
    endif
    call elphon_comm_init(wdes_dense_ibz, elph_set%kpar, elph_set%ncore, comm_kbatch, comm_global, io)
! Check if the NCORE setting is compatible with the (1._q,0._q) in the SCF calculation
    if (wdes_dense_ibz%comm_inb%ncpu/=w%wdes%comm_inb%ncpu) then
       call vtutor%error('The NCORE='//str(w%wdes%comm_inb%ncpu)//' setting in the INCAR cannot be used in the electron-phonon calculation. '//&
                         'You might need to change ELPH_KPAR <= Number of MPI ranks / NCORE.')
    endif
! TODO: change this
    elph_mels%comm_kbatch = comm_kbatch
    elph_mels%comm_global = comm_global

! set number of bands based on the input by user or using number of plane-waves
    call elphon_set_nbands(wdes_dense_ibz,elph_set,info%nelect,t_info,maxval(w%wdes%nplwkp_tot))
    call wdes_set_npro(wdes_dense_ibz,t_info,p,info%loverl)

! read kpoint_dense file
    call elphon_get_kpoints_inter(kpoints_inter, kgrid, kpoints, elph_set, latt_cur, t_info, wdes_dense_ibz, comm_kbatch, symm%isym, io)

! initialize elphon object
    call elph_mels_init(elph_mels,elph_set,wdes_dense_ibz,lmdim,grid,latt_ini,latt_cur,info,p,symm,t_info,io)

! read dimensions from phelel_params.hdf5
    if (any(elph_set%selfen_approx(:)%scattering_approx/=elph_scattering_approx_crta)) then
        call elph_mels_read_dimensions(elph_mels,filename='phelel_params.hdf5')
    endif

! prepare WF descriptor in the dense IBZ
    call elph_mels_get_w_dense(elph_mels,kpoints_inter,w%wdes%ngdim)

! compute KS states on dense grid
    call elph_mels_compute_wfs(elph_mels,kpoints,kpoints_inter,kgrid,w,hamiltonian,nonlr_s,nonl_s,&
                               cqij,cdij,sv,chempot,bandgap)

    if (elph_set%driver ==  elph_driver_ph) then
! initialize phonon accumulators
        call elph_accumulators_init_ph(elph_accumulators_ph,elph_set,elph_mels%wdes_dense_ibz)

        call elphon_prepare_selfenergy_accumulators_ph(&
                        elph_accumulators_ph%selfen_ph, latt_cur, kpoints_inter, elph_set, &
                        elph_mels%wdes_dense_ibz, elph_mels%w_dense_ibz, t_info, info, chempot, io)
    else
! initialize accumulators
        call elph_accumulators_init(elph_accumulators,elph_set,elph_mels%wdes_dense_ibz)

        call elphon_prepare_selfenergy_accumulators(elph_accumulators%selfen_el, elph_accumulators%id_name, elph_accumulators%id_size, latt_cur, kpoints_inter, elph_set, elph_mels%wdes_dense_ibz, &
                                                    elph_mels%w_dense_ibz, t_info, info, chempot, bandgap, io)
    endif

    if (elph_set%transport.or.any(elph_set%selfen_approx(:)%scattering_approx==elph_scattering_approx_mrta_lambda)&
                          .or.any(elph_set%selfen_approx(:)%scattering_approx==elph_scattering_approx_mrta_tau)&
                          .or.any(elph_set%selfen_approx(:)%scattering_approx==elph_scattering_approx_erta_lambda)&
                          .or.any(elph_set%selfen_approx(:)%scattering_approx==elph_scattering_approx_erta_tau)) then
        call elphon_get_velocity(velocity, p, lmdim, elph_mels%nbands, elph_mels%wf_ncdij, latt_cur, latt_ini, symm, nonlr_s, nonl_s, &
                                 t_info, elph_set, elph_mels%wdes_dense_ibz, elph_mels%w_dense_ibz, &
                                 elph_mels%wave_map, elph_mels%wdes_kbatch, elph_mels%w_kbatch, sv, cqij, cdij, &
                                 grid, gridc, gridus, c_to_us, irdmax, info, io, comm_kbatch)
        call elphon_write_velocity(elph_set, velocity, elph_mels%w_dense_ibz%celtot, io)
    endif
    deallocate(elph_mels%w_kbatch%CPROJ)

! if only CRTA is requested we can skip the computation of the electron self-energy altogether
    if (elph_set%transport.and.all(elph_set%selfen_approx(:)%scattering_approx==elph_scattering_approx_crta)) then
        call transport_elphon(elph_set,elph_accumulators%selfen_el,velocity,latt_cur,kgrid,kpoints_inter,&
                              elph_mels%w_dense_ibz,elph_mels%wdes_dense_ibz,bandgap,symm,info,io,elph_mels%comm_global)
        
        return
    endif

! initialize the force constants from phelel_params.hdf5
    call elph_mels_ifc_init(elph_mels)
    call elph_mels_ifc_read_hdf5(elph_mels,filename='phelel_params.hdf5')
    call elph_mels_ifc_setup(elph_mels)
    call elph_mels_ifc_compute_dispersion(elph_mels,qpoints_exists,filename='QPOINTS')
    if (qpoints_exists) then
        
        return
    endif
    call elph_mels_ifc_compute_ibz_phonons(elph_mels,kpoints_inter)

! initialize the electron-phonon potential from phelel_params.hdf5
    call elph_mels_pot_init(elph_mels,w%wdes%COMM_intra_node)
    call elph_mels_pot_read_hdf5(elph_mels,filename='phelel_params.hdf5')
    call elph_mels_pot_setup(elph_mels)
    call elph_mels_pot_dispersion(elph_mels,qpoints_exists,filename='QPOINTS_POT')
    if (qpoints_exists) then
        
        return
    endif

! allocate caches
    if (elph_set%driver == elph_driver_ph) then
        call elph_mels_set_nbandsk(elph_mels,minval(elph_accumulators_ph%selfen_ph(:)%band_start),maxval(elph_accumulators_ph%selfen_ph(:)%band_stop))
    else
        call elph_mels_set_nbandsk(elph_mels,minval(elph_accumulators%selfen_el(:)%band_start),maxval(elph_accumulators%selfen_el(:)%band_stop))
    end if
    call elph_mels_wf_write(elph_mels)
    call elph_mels_cache_init(elph_mels)

    call elph_mels_windows_init(elph_mels)
    call elph_mels_wave_caches_init(elph_mels)

    if (elph_set%driver == elph_driver_ph) then
        call elph_mels_wf_redistribute_ph(elph_mels,elph_accumulators_ph,kpoints_inter)
        call elph_mels_cache_prefill_ph(elph_mels,elph_accumulators_ph)
    else
        call elph_mels_wf_redistribute(elph_mels,elph_accumulators,kpoints_inter)
        call elph_mels_cache_prefill(elph_mels,elph_accumulators)
    end if

! generate integer list of all the k' points that need to be handled locally
    call elph_mels_get_ikpfbz(elph_mels,kpoints_inter,elph_mels%ikp_fbz_list)

! report memory usage
    call dump_allocate(io%iu6)
    call wforce(io%iu6)

! main loop where electron-phonon matrix elements are computed and accumulated
    select case (elph_set%driver)
        case(elph_driver_el)
            call electron_phonon_selfen_el_driver(elph_set,elph_mels,elph_accumulators,kpoints_inter,velocity)
        case(elph_driver_ph)
            call electron_phonon_selfen_ph_driver(elph_set,elph_mels,elph_accumulators_ph,kpoints_inter)
        case(elph_driver_mels)
            call electron_phonon_mels_driver(elph_set,elph_mels,kpoints_inter)
    end select

    call elph_mels_barrier(elph_mels)
    call elph_mels_windows_free(elph_mels)
    call elph_mels_wave_caches_free(elph_mels)

! accumulators
    if (elph_set%driver == elph_driver_ph) then
        call elph_accumulators_finalize_ph(elph_accumulators_ph,comm_global)
        call elph_accumulators_print_ph(elph_accumulators_ph,elph_mels%w_dense_ibz,elph_mels%phonon_freqs_ibz,kpoints_inter,io)
    else
        call elph_accumulators_finalize(elph_accumulators,elph_mels%w_dense_ibz,comm_global)
        call elph_accumulators_print(elph_accumulators,elph_mels%w_dense_ibz,io)
    end if

! transport calculation
    if (elph_set%transport) then
!call elphon_get_velocity(velocity, p, lmdim, elph_mels%nbands, elph_mels%wf_ncdij, latt_cur, latt_ini, symm, nonlr_s, nonl_s, &
!                         t_info, elph_set, elph_mels%wdes_dense_ibz, elph_mels%w_dense_ibz, &
!                         elph_mels%wave_map, elph_mels%wdes_kbatch, elph_mels%w_kbatch, sv, cqij, cdij, &
!                         grid, gridc, gridus, c_to_us, irdmax, info, io, comm_kbatch, comm_global)
        call transport_elphon(elph_set,elph_accumulators%selfen_el,velocity,latt_cur,kgrid,kpoints_inter,&
                              elph_mels%w_dense_ibz,elph_mels%wdes_dense_ibz,bandgap,symm,info,io,elph_mels%comm_global)
    endif

! free objects
    if (elph_set%driver == elph_driver_ph) then
        call elph_accumulators_free_ph(elph_accumulators_ph)
    else
        call elph_accumulators_free(elph_accumulators)
    end if
    call elph_mels_cache_free(elph_mels)
    call elph_mels_ifc_free(elph_mels)
    call elph_mels_pot_free(elph_mels)
    call elph_mels_wfs_free(elph_mels,kpoints,w,nonl_s)
    call elph_mels_free(elph_mels)
    

    call elph_footer_write(io)

# 273


 end subroutine electron_phonon_driver

!> Main loop for the computation of the electron-phonon self-energy
 subroutine electron_phonon_selfen_el_driver(elph_set,elph_mels,elph_accumulators,kpoints,velocity)
    USE elphon_base, ONLY: find_value, elph_progress
    USE elphon_mels, ONLY: elph_mels_compute_g, elph_mels_get_ikpfbz
    USE elphon_accumulators, ONLY: elph_accumulators_accumulate, elph_accumulators_contributes, elph_accumulators_prepare_dw,&
                                   elph_accumulators_shallido, elph_accumulators_delta_is_zero
    USE triplets, ONLY: triplet_type, triplet_init_kp_and_q, kpoint_init_kindex, kpoint_init_kibz, kpoint_init_vkpt,&
                        kpoint_init_kint
    USE mkpoints, ONLY: kindex_from_kp_and_q, ikfbz_kindex, vkpt_kindex, kindex_kint, kint_mod, kint_kindex
    type(elph_settings) :: elph_set
    type(elph_mels_type) :: elph_mels
    type(elph_accumulators_selfen) :: elph_accumulators
    type(kpoints_struct) :: kpoints
    real(q), allocatable, intent(in) :: velocity(:,:,:,:)
! local variables
    type(triplet_type) :: triplet
    logical :: contributes
    integer :: isp
    integer :: iq_fbz,ikp_fbz
    integer :: ik_ibz
    integer :: ik
    integer :: ibz2kindex(kpoints%nkpts)
    complex(q),allocatable :: gcart_kkp(:,:,:,:)
    complex(q),allocatable :: gmode_kkp(:,:,:)

    if (elph_set%dry_run) return

! create mapping from the ibz to kindex
    do ik_ibz=1,kpoints%nkpts
        ibz2kindex(ik_ibz) = vkpt_kindex(kpoints,elph_mels%latt_cur,kpoints%vkpt(:,ik_ibz))
    enddo

    allocate(gcart_kkp(elph_mels%nbands_k,elph_mels%nbands_local,3,elph_mels%natoms_pc))
    allocate(gmode_kkp(elph_mels%nbands_k,elph_mels%nbands_local,3*elph_mels%natoms_pc))

! prepare computation of the debye waller term
    call elph_accumulators_prepare_dw(elph_accumulators,elph_set,elph_mels%fbz2ibz,elph_mels%elph_ifc,elph_mels%wdes_dense_ibz,elph_mels%comm_global)
    
! loop for electron self-energy
    do isp=1,elph_mels%ispin


! loop over the q points in the fbz
        do iq_fbz=1,kpoints%nkfull
            call elph_progress('fbz qpoint',iq_fbz,kpoints%nkfull,elph_mels%io%iu0)
            call kpoint_init_kindex(triplet%q, kpoints, iq_fbz)

! loop over k'
            do ik=1,size(elph_mels%ikp_fbz_list)
                ikp_fbz = elph_mels%ikp_fbz_list(ik)
                call kpoint_init_kindex(triplet%kp, kpoints, ikp_fbz)

                call kpoint_init_kint(triplet%k, kpoints, triplet%q%kint+triplet%kp%kint)

! skip if not ibz
                if (triplet%k%kindex/=ibz2kindex(triplet%k%ibz)) cycle

! check if delta(ek-ek'+wq) or delta(ek-ek'-wq) function is (0._q,0._q)
                contributes = elph_accumulators_contributes(elph_accumulators,elph_set,elph_mels%w_dense_ibz,elph_mels%phonon_freqs_ibz,&
                                                            elph_accumulators%band_idx,kpoints,triplet,isp)
                if (.not.contributes) cycle
! compute the electron-phonon matrix elements
                call elph_mels_compute_g(elph_mels,kpoints,triplet%k%kindex,ikp_fbz,iq_fbz,isp,gcart_kkp,gmode_kkp)
! use the electron-phonon matrix elements to compute the electron self-energy
                call elph_accumulators_accumulate(elph_accumulators,elph_set,elph_mels%band_start_k,elph_mels%band_start_kp,triplet,isp,kpoints,elph_mels%waverot,&
                                                  elph_mels%w_dense_ibz,elph_mels%wdes_dense_ibz,elph_mels%phonon_freqs_ibz,elph_mels%latt_cur,gmode_kkp,gcart_kkp, &
                                                  velocity)
            enddo !loop over k'
        enddo !loop over q
# 371

    enddo ! loop over spin
    
 end subroutine electron_phonon_selfen_el_driver

!> Main loop for the computation of electron-phonon matrix elements
 subroutine electron_phonon_mels_driver(elph_set,elph_mels,kpoints)

    USE vhdf5, ONLY: vh5_write_lattice_ion
    USE elphon_base, ONLY: vh5_write_elph_kpoints_regular, vh5_write_elph_kpoints_list, vh5_write_elph_eigenvalues_list

    USE elphon_base, ONLY: find_value, &
                           get_elph_kpoints_eigenvals_list, get_elph_nkpoints_local, &
                           elph_progress, wdes_get_global_band_index
    USE elphon_mels, ONLY: elph_mels_compute_g, elph_mels_compute_g_ibz, elph_mels_get_ikpfbz
    USE elphon_accumulators, ONLY: elph_accumulator_g, elph_accumulators_g_init, elph_accumulators_g_accumulate, &
                                   elph_accumulators_g_free
    USE triplets, ONLY: triplet_type, triplet_init_kp_and_q, kpoint_init_kindex, kpoint_init_kibz, kpoint_init_vkpt,&
                        kpoint_init_kint
    USE mkpoints_struct_def, ONLY: LineMode, ExplicitList
    USE mkpoints, ONLY: kindex_from_kp_and_q, ikfbz_kindex, vkpt_kindex, kindex_kint, kint_mod, kint_kindex
    USE wave_rotate, ONLY: kpoint_star, wave_rotator_get_star
    USE wave_struct_def, ONLY: wavedes1
    USE wave, ONLY: setwdes
    USE vhdf5
    type(elph_settings) :: elph_set
    type(elph_mels_type) :: elph_mels
    type(kpoints_struct) :: kpoints
! local variables
    type(triplet_type) :: triplet
    type(elph_accumulators_g) :: elph_accumulators
    type(wavedes1) :: wdes1
    type(kpoint_star) :: star

    integer(HID_T) :: fileid

    logical :: lhdf5_mpiio
    integer :: isp, ikp_progress
    integer :: nkpts_k_local, nkpts_kp_local
    integer :: ik_ibz
    integer :: ikp_ibz, ikp_write
    integer :: istar, ik, ikp
    integer :: nkpts_k, nkpts_kp
    integer,allocatable :: kpoint_index(:) ! indexes of the k-points to compute
    integer,allocatable :: band_idx(:)
    real(q),allocatable :: wtkpt_k(:)
    real(q),allocatable :: vkpt_k(:,:)
    real(q),allocatable :: vkpt_kp(:,:)
    real(q),allocatable :: eig_k(:,:,:)
    real(q),allocatable :: eig_kp(:,:,:)
    complex(q),allocatable :: gcart_kkp(:,:,:,:)
    complex(q),allocatable :: gmode_kkp(:,:,:)

    allocate(gcart_kkp(elph_mels%nbands_k,elph_mels%nbands_local,3,elph_mels%natoms_pc))
    allocate(gmode_kkp(elph_mels%nbands_k,elph_mels%nbands_local,3*elph_mels%natoms_pc))

! select k-points
    if (allocated(elph_set%selfen_ikpt)) then
        allocate(kpoint_index(size(elph_set%selfen_ikpt)))
        kpoint_index = elph_set%selfen_ikpt
    endif


    call get_elph_kpoints_eigenvals_list(kpoints,elph_mels%fbz2ibz,kpoint_index,elph_mels%w_dense_ibz%celtot,&
                                         nkpts_k,nkpts_kp,wtkpt_k,vkpt_k,vkpt_kp,eig_k,eig_kp,lwriteh5)

    call get_elph_nkpoints_local(elph_mels%wave_map,elph_mels%waverot,elph_mels%fbz2ibz,kpoints,nkpts_k,nkpts_k_local,nkpts_kp_local)
    call wdes_get_global_band_index(elph_mels%wdes_dense_ibz,band_idx)

    lhdf5_mpiio = elph_mels%elph_set%usehdf5mpiio
    call elph_accumulators_g_init(elph_accumulators,elph_mels%band_start_k,elph_mels%band_stop_k,elph_mels%nbands_k,&
                                                    elph_mels%band_start_kp,elph_mels%band_stop_kp,elph_mels%nbands_kp,&
                                                    nkpts_k,nkpts_k_local,nkpts_kp,nkpts_kp_local,elph_mels%ispin,elph_mels%natoms_pc,band_idx,lhdf5_mpiio,elph_mels%comm_global)

! need wdes1 to write matrix elements block
    call setwdes(elph_mels%wdes_dense_ibz,wdes1,1)

    
    if (kpoints%mode==LineMode.or.kpoints%mode==ExplicitList) then

! fix k and k' is k-points

! loop for computing electron-phonon matrix elements along a path
        do isp=1,elph_mels%ispin
! loop over k' points in the local ibz
            do ikp=1,elph_mels%wave_map%nkpts_batch
                call elph_progress("k'-point",ikp,elph_mels%wave_map%nkpts_batch,elph_mels%io%iu0)
                ikp_ibz = elph_mels%wave_map%kpoints_index(ikp+elph_mels%wave_map%nkpts_orig)
                do ik=1,size(kpoint_index)
                    ik_ibz = kpoint_index(ik)
! compute the electron-phonon matrix elements
                    call elph_mels_compute_g_ibz(elph_mels,kpoints,ik_ibz,ikp_ibz,isp,gcart_kkp,gmode_kkp)
! accumulate matrix elements
                    call elph_accumulators_g_accumulate(elph_accumulators,wdes1,ik,ikp_ibz,isp,gmode_kkp,elph_mels%w_q)
                enddo !loop over k
            enddo !loop over k'
        enddo ! loop over spin

! TODO fix q k is kpoints and k'

    else
! loop for computing electron-phonon matrix elements on a regular k-point grid
        do isp=1,elph_mels%ispin
! loop over k' points in the local ibz
!do ik=1,size(elph_mels%ikp_fbz_list)
!    ikp_fbz = elph_mels%ikp_fbz_list(ik)
!    call kpoint_init_kindex(triplet%kp, kpoints, elph_mels%latt_cur, ikp_fbz)

! loop over k' points in the local ibz
            ikp_progress = 0
            do ikp=1,elph_mels%wave_map%nkpts_batch
                ikp_ibz = elph_mels%wave_map%kpoints_index(ikp+elph_mels%wave_map%nkpts_orig)
! loop over points in the star of k'
                call wave_rotator_get_star(elph_mels%waverot,kpoints%vkpt(:,ikp_ibz),star)
                do istar=1,star%nkpts
                    call kpoint_init_vkpt(triplet%kp, kpoints, elph_mels%latt_cur, star%vkpt(:,istar))
                    ikp_write = sum(elph_mels%fbz2ibz%star_size(1:ikp_ibz-1)) + istar
                    ikp_progress = ikp_progress + 1
                    call elph_progress("k' point",ikp_progress,nkpts_kp_local,elph_mels%io%iu0)

! loop over k in the ibz
                    do ik=1,size(kpoint_index)
                        ik_ibz = kpoint_index(ik)
                        call kpoint_init_kibz(triplet%k, kpoints, elph_mels%latt_cur, ik_ibz)
! compute q=k-kp
                        call kpoint_init_kint(triplet%q, kpoints, triplet%k%kint-triplet%kp%kint)

! compute the electron-phonon matrix elements
                        call elph_mels_compute_g(elph_mels,kpoints,triplet%k%kindex,triplet%kp%kindex,triplet%q%kindex,isp,gcart_kkp,gmode_kkp)

! accumulate matrix elements
                        call elph_accumulators_g_accumulate(elph_accumulators,wdes1,ik_ibz,ikp_write,isp,gmode_kkp,elph_mels%w_q)
                    enddo !loop over k
                enddo
            enddo !loop over k'
        enddo ! loop over spin
    endif
    


    call elph_accumulators_g_free(elph_accumulators)

! write data to HDF5 file
    call vh5_error(vh5_file_open_readwrite('vaspelph.h5',fileid),"elphon_driver.F",514)
    call vh5_write_lattice_ion(fileid, elph_mels%info%sznam1, elph_mels%t_info, elph_mels%latt_cur%scale, elph_mels%latt_cur%a,&
                               .false., elph_mels%elph_pot%primitive_positions)
    call vh5_write_elph_kpoints_regular(fileid,elph_mels%waverot%nrotk,elph_mels%waverot%igrpop,&
                                elph_mels%fbz2ibz%irot_fbz2ibz,elph_mels%fbz2ibz%indx_fbz2ibz)
    call vh5_write_elph_kpoints_list(fileid,wtkpt_k,vkpt_k,vkpt_kp)
    if (allocated(vkpt_k))  deallocate(vkpt_k)
    if (allocated(vkpt_kp)) deallocate(vkpt_kp)

! write eigenvalues
    call vh5_write_elph_eigenvalues_list(fileid,eig_k,eig_kp)
    if (allocated(eig_k))  deallocate(eig_k)
    if (allocated(eig_kp)) deallocate(eig_kp)
    if (lwriteh5) then
        call vh5_error(vh5_file_close(fileid),"elphon_driver.F",528)
    endif


 end subroutine electron_phonon_mels_driver

end module elphon_driver
