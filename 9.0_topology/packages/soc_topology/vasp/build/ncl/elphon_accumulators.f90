# 1 "elphon_accumulators.F"
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


# 2 "elphon_accumulators.F" 2 
!> This module implements structures and routines to compute
!> quantities using electron-phonon matrix elements
module elphon_accumulators
    use prec
    use constant
    use mpimy, only: communic
    use bandgap_struct, only: bandgap_info

    use vhdf5

# 14

    USE elphon_base, ONLY : elph_selfen_approx_t
    implicit none

!> @brief Datatype holding the information of the self-energy
    type elph_accumulator_selfen
        integer :: nbands !< number of bands at which to compute the self-energy
        integer :: nbands_sum !< number of bands to sum over
        integer :: natoms !< number of atoms in the system
        integer :: nkpoints !< number of kpoints
        integer :: nspin !< number of spin channels
!> spin degeneracy: 2 for spin unpolarized, 1 for collinear and non collinear spin calculations
        integer :: spindeg !< spin degeneracy
        integer :: nw !< number of energies at which to evaluate the self-energy
        real(q) :: enwin !< range of energies
        real(q) :: delta !< small imaginary part to add to the self-energy
        real(q),allocatable :: carrier_per_cell(:) !< number of electrons per cell (for the corresponding temperatures)
        real(q) :: nelect !< number of electrons per cell (without doping)
        real(q),allocatable :: efermi(:) !< fermi energies (for the corresponding temperatures)
        real(q) :: emin !< energy of the lowest band which will be used in the sum-over-states
        real(q) :: emax !< energy of the highest band which will be used in the sum-over-states
        integer :: band_start !< index of the band at which to start computing the self-energy
        integer :: band_stop  !< index of the band at which to stop computing the self-energy
        type(elph_selfen_approx_t) :: approx !< type of approximation used to compute the fan self-energy
        integer,allocatable :: id_idx(:)  !< indexes for which this accumulator was generated
        integer,allocatable :: id_size(:)  !< size of the indexes for which this accumulator was generated
        character(len=15),allocatable :: id_name(:) !< name of the fields for which this accumulator was generated
!> used to determine the energy beyond which the delta functions are considered to be (0._q,0._q) in broadening methods
        real(q) :: broad_tol
        logical :: transport !< compute transport properties
        real(q) :: transport_emin !< minimum energy of window for transport computation
        real(q) :: transport_emax !< maximum energy of window for transport computation
        integer :: transport_nedos !< number of energy steps in the energy window
        logical :: tetrahedron !< use tetrahedron integration method
        logical :: fan !< store the Fan contribution to the self-energy
        logical :: dfan !< store the Fan contribution to the self-energy
        logical :: dw !< store the Debye-Waller contribution to the self-energy
!> Energy cutoff (eV) for choosing which states to compute.
!> There are two cases:
!>  - metal: the self-energy for all the states within
!>           [efermi-select_energy_window(1):efermi+select_energy_window(2)] is computed
!>  - gaped system: all the states within [vbm-select_energy_window(1):cbm+select_energy_window(1)]
        real(q) :: select_energy_window(2)
        integer :: energy_grid_mode !< definition of the energy grid: 0-centered around KS eigenvalue, 1-in the range of total energies
!> energy grids at which to evaluate the self-energy in eV (nw,nbks)
!> where nbks is a combined index for bands, kpoints and spin (see bks_idx)
        real(q),allocatable :: energies(:,:)
!> store the fan self-energy. The dimensions are (ntemps,nw,nbks)
!> where nbks is a combined index for bands, kpoints and spin (see bks_idx)
        complex(q),allocatable :: selfen_fan(:,:,:)
!> store the derivative of the fan self-energy w.r.t w. The dimensions are (ntemps,nw,nbks)
!> where nbks is a combined index for bands, kpoints and spin (see bks_idx)
        complex(q),allocatable :: selfen_dfan(:,:,:)
!> store the Debye-Waller self-energy. The dimensions are (ntemps,nbks)
!> where nbks is a combined index for bands, kpoints and spin (see bks_idx)
        real(q),allocatable :: selfen_dw(:,:)
!> auxiliary array for the computation of the debye-waller
!> term using the rigid-ion approximation
        complex(q),allocatable :: t(:,:,:,:,:,:)
        integer :: ntemps !< number of temperatures
        real(q),allocatable :: tempsk(:) !< array of temperatures in kelvin
        logical,allocatable :: compute_mask(:,:,:) !< integer controlling whether to compute self-energy at ikpt and band
        integer :: nbks !< Number of combined bands, kpoints and spin indexes (see bks_idx)
        integer,allocatable :: bks_idx(:,:,:) !< Combine band, kpoint and spin index into (1._q,0._q)
!> displacements of 24 tetrahedra that contribute to a k-point
        integer :: kcut24(3,4,24)
!> index of the corner of each of the 24 tetrahedra for which the k-point has 0 displacement
        integer :: kcut24_origin(24)
        logical :: kcut24_initialized
!$ integer(kind=omp_lock_kind),allocatable :: omp_fan_lock(:)
! to compute the bandgaps
        logical :: report_gaps !< report bandgap and renormalization after the electron-phonon calculation
        type(bandgap_info),allocatable :: bandgap(:) !< bandgap information
        real(q),allocatable :: direct_gap(:)
        real(q),allocatable :: direct_gap_renorm(:,:)
        real(q),allocatable :: fundamental_gap(:)
        real(q),allocatable :: fundamental_gap_renorm(:,:)
    end type elph_accumulator_selfen

!> container for g-matrix elements which will be stored in memory in a distrbuted fashion
    type elph_accumulator_g
        integer :: nkpts_k
        integer :: nkpts_kp
        integer :: nkpts_k_local
        integer :: nkpts_kp_local
        integer :: nbands_k
        integer :: band_start_k
        integer :: band_stop_k
        integer :: nbands_kp
        integer :: nbands_kp_local
        integer :: band_start_kp
        integer :: band_stop_kp
        integer :: nspin
        integer :: natoms
        integer :: nkpts
        integer :: ikkp_counter
        integer :: ikkp_max
!> electron-phonon matrix elements
        complex(q),allocatable :: g(:,:,:,:,:)
!> phonon frequencies
        real(q),allocatable :: w_q(:,:)
        integer,allocatable :: idx_k(:)
        integer,allocatable :: idx_kp(:)
! paralelism
        type(communic) :: comm
        integer,allocatable :: band_idx(:)
    end type elph_accumulator_g

!> @brief Container for information related to writting the electron-phonon data to file
    type elph_accumulator_g_hdf5mpiio

        integer(hid_t) :: fileid !< The fileid where the data is going to be written
        integer(hid_t) :: mels_groupid !< The group id of the matrix elements group
        integer(hid_t) :: phonons_groupid !< The group id of the phonons group
        integer(hid_t) :: plist_id

        logical :: lwritesave !< Save what is the master CPU
        integer :: natoms
        integer :: nkpts_k !< Number of points in the IBZ
        integer :: nkpts_kp !< Number of points in the FBZ
        integer :: nbands_k !< Number of bands at k
        integer :: nbands_kp !< Number of bands at k'
        integer :: nrotk !< Number of rotation operations
    end type elph_accumulator_g_hdf5mpiio

!> contains a list of accumulators for electron self-energy
    type elph_accumulators_selfen
!> number of accumulators_el
        integer :: naccumulators_el
!> name of the id of each accumulator
        character(len=15) :: id_name(4)
!> number of possible elements for each id
        integer :: id_size(4)
!> global index of the local bands
        integer,allocatable :: band_idx(:) !< global index of the bands that are present locally
!> electron self-energy accumulators
        type(elph_accumulator_selfen),allocatable :: selfen_el(:)
    end type elph_accumulators_selfen

!> contains a list of accumulators for g matrix elements
    type elph_accumulators_g
!> logical to determine wether to use hdf5 with mpiio support
        logical :: lhdf5_mpiio
!> electron-phonon matrix element accumulator
        type(elph_accumulator_g) :: g_store
!> electron-phonon matrix element accumulator to hdf5mpiio
        type(elph_accumulator_g_hdf5mpiio) :: g_hdf5mpiio
    end type elph_accumulators_g

    real(q), parameter :: PHON_TOL=1e-4 !< phonon tolerance to skip the calculation (eV)
    real(q), parameter :: EL_TOL=4e-5 !< electron tolerance to detect degenerate eigenvalues (eV)
    real(q), parameter :: PLANCK_SI = 6.626070040e-34_q
    real(q), parameter :: THZTOEV = 1e12 * PLANCK_SI / EVTOJ / 2 / PI
    contains

!> @brief Intiialize the accumulator object for the self-energy
    subroutine elph_accumulator_selfen_init(self,nbands,nbands_sum,nkpoints,nspin,rspin,nw,enwin,delta,&
                                            ntemps,temps,band_start,band_stop,approx,broad_tol,transport,&
                                            natoms,eigenvalues,nb_totk,efermi,carrier_per_cell,nelect,&
                                            fan,dfan,dw,id_idx,id_size,id_name)
        use constant
        use ini, only: register_allocate
        use elphon_base, only: celtot_max
        type(elph_accumulator_selfen) :: self
        integer,intent(in) :: nw
        integer,intent(in) :: nbands
        integer,intent(in) :: nbands_sum
        integer,intent(in) :: nkpoints
        integer,intent(in) :: nspin
        real(q),intent(in) :: rspin
        integer,intent(in) :: natoms
        real(q),intent(in) :: enwin
        real(q),intent(in) :: delta
        real(q),intent(in) :: efermi(ntemps)
        real(q),intent(in) :: carrier_per_cell(ntemps) !< number of carriers per cell (with doping)
        real(q),intent(in) :: nelect                   !< number of carriers per cell (without doping)
        integer,intent(in) :: ntemps
        integer,intent(in) :: band_start
        integer,intent(in) :: band_stop
        type(elph_selfen_approx_t),intent(in) :: approx !< use static approximation for the electron self-energy
        real(q),intent(in) :: broad_tol
        logical,intent(in) :: transport
        logical,intent(in) :: fan !< store the Fan contribution to the self-energy
        logical,intent(in) :: dfan !< store the derivative of the Fan contribution to the self-energy
        logical,intent(in) :: dw !< store the Debye-Waller contribution to the self-energy
        integer,intent(in) :: id_idx(:) !< integer indentifying this electron self-energy accumulator
        integer,intent(in) :: id_size(:) !< size of the integer indentifying this electron self-energy accumulator
        character(len=*),intent(in) :: id_name(:) !< list of names indentifying this electron self-energy accumulator
!> temperatures at which the self-energy is to be evaluated (in K)
        real(q),intent(in) :: temps(ntemps)
        complex(q),intent(in) :: eigenvalues(:,:,:)
        integer,intent(in) :: nb_totk(:,:)

        self%nbands = nbands
        self%nbands_sum = nbands; if (nbands_sum>0.and.nbands_sum<nbands) self%nbands_sum = nbands_sum
        self%approx = approx
        self%id_idx = id_idx
        self%id_size = id_size
        self%id_name = id_name
        self%broad_tol = broad_tol
        self%transport = transport
        self%tetrahedron = delta<=0
        self%kcut24_initialized = .false.
        self%nkpoints = nkpoints
        self%nspin = nspin
        self%spindeg = nint(rspin)
        if (nw>0) then
            self%nw = nw+mod(nw+1,2)
            self%energy_grid_mode=0
        else
            self%energy_grid_mode=1
            self%nw = max(2,abs(nw))
        endif
        self%emax=celtot_max(eigenvalues,nb_totk)+5.0_q
        self%emin = minval(real(eigenvalues,q))-5.0_q
        self%delta = abs(delta)
        self%enwin = enwin
        self%band_start = 1; if (band_start>0.and.band_start<nbands) self%band_start = band_start
        self%band_stop = nbands; if (band_stop>0.and.band_stop<nbands) self%band_stop = band_stop
        self%ntemps = ntemps
        allocate(self%tempsk(self%ntemps))
        self%tempsk = temps
        self%natoms = natoms
! logicals wether to perform fan and dw calculation
        self%fan = fan
        self%dfan = dfan
        self%dw = dw
        allocate(self%t(self%ntemps,3,self%natoms,3,self%natoms,3*self%natoms))
        self%t=cmplx(0.0_q,0.0_q,q)
        allocate(self%efermi(self%ntemps))
        self%efermi=efermi
        allocate(self%carrier_per_cell(self%ntemps))
        self%carrier_per_cell=carrier_per_cell
        self%nelect=nelect
        self%report_gaps=.false.

! allocate compute_mask array with the information of whether a particular
! k-point should be computed
        allocate(self%compute_mask(self%band_stop,self%nkpoints,self%nspin))
        call register_allocate(0.125_q*STORAGE_SIZE(self%compute_mask)*SIZE(self%compute_mask,KIND=qi8),'elph_accumulator_selfen')
        self%compute_mask = .false.
    end subroutine elph_accumulator_selfen_init

    subroutine elph_accumulator_selfen_set_bandgap(self,bandgap)
        use bandgap_struct
        type(elph_accumulator_selfen) :: self
        type(bandgap_info), intent(in) :: bandgap(:)

        allocate(self%bandgap(size(bandgap)))
        self%bandgap(:) = bandgap(:)
    end subroutine elph_accumulator_selfen_set_bandgap

!> Select bands and k-point for which to compute the self-energy
!> based on the interval from vbm-elph_selfen_energy_window(1) up to to vbm
!>            and from cbm to cbm+elph_selfen_energy_window(2)
!> On output compute_mask is true for the states for which we want to compute the self-energy
    subroutine elph_accumulator_selfen_selector_energy_window(self,w,select_energy_window)
        use wave_struct_def, only: wavespin
        type(elph_accumulator_selfen) :: self
        type(wavespin),intent(in) :: w
        real(q),intent(in) :: select_energy_window(2)
! local variables
        integer :: ispin,ikibz,ibibz,iw,ikpt
        integer :: nelect
        real(q) :: enk, vbm, cbm
        self%select_energy_window = select_energy_window

        self%compute_mask = .false.
! spin unpolarized and non-collinear spin case
        vbm = self%bandgap(1)%fundamental_valence%eigenvalue
        cbm = self%bandgap(1)%fundamental_conduction%eigenvalue
        do ikibz=1,self%nkpoints
            do ibibz=1,self%band_stop
                enk = real(w%celtot(ibibz,ikibz,1),q)
                self%compute_mask(ibibz,ikibz,1) = &
                self%compute_mask(ibibz,ikibz,1) .or. &
                ( enk < cbm + self%select_energy_window(2) .and. &
                  enk > vbm - self%select_energy_window(1) )
            enddo ! bands
        enddo ! kpoints

! for collinear spin case we might select a few more states
! such as to compute the spin-dependent gaps as well
        if (self%nspin==2) then
          do ispin=1,self%nspin
            vbm = self%bandgap(1+ispin)%fundamental_valence%eigenvalue
            cbm = self%bandgap(1+ispin)%fundamental_conduction%eigenvalue
            do ikibz=1,self%nkpoints
                do ibibz=1,self%band_stop
                    enk = real(w%celtot(ibibz,ikibz,ispin),q)
                    self%compute_mask(ibibz,ikibz,ispin) = &
                    self%compute_mask(ibibz,ikibz,ispin) .or. &
                    ( enk < cbm + self%select_energy_window(2) .and. &
                      enk > vbm - self%select_energy_window(1) )
                enddo ! bands
            enddo ! kpoints
          enddo !ispin
        endif
    end subroutine elph_accumulator_selfen_selector_energy_window

!> Choose the k-points to compute explicitely
    subroutine elph_accumulator_selfen_selector_ikpt(self,ikpts)
        type(elph_accumulator_selfen) :: self
        integer :: ikpts(:)
! local variables
        integer :: ikpt, ikibz, isp
! loop over ikpts and set compute_mask
        do isp=1,self%nspin
            do ikpt=1,size(ikpts)
                ikibz = ikpts(ikpt)
! Here I set both spin components
! TODO: setting each spin channel separately is not implemented for now
                self%compute_mask(:,ikibz,isp) = .true.
            enddo
        enddo
    end subroutine elph_accumulator_selfen_selector_ikpt

!> Find the direct and indirect gaps and select those k-points to be computed
    subroutine elph_accumulator_selfen_selector_gaps(self,w)
        use wave_struct_def, only: wavespin
        type(elph_accumulator_selfen) :: self
        type(wavespin) :: w
! local variables
        integer :: ib, ispin, ikibz, nelect
        integer :: ikpt_cbm, ikpt_vbm, ikpt_direct
        integer :: ispin_vbm_dummy, ispin_cbm_dummy
        real(q) :: e, vbm, cbm

        self%report_gaps = .true.

! find direct gaps
        self%compute_mask = .false.
        do ispin=1,self%nspin
            call update_mask(self%bandgap(1)%direct_valence)
            call update_mask(self%bandgap(1)%direct_conduction)
        enddo

! Also select to compute the states for the per-spin gaps
        if (self%nspin==2) then
          do ispin=1,self%nspin
            call update_mask(self%bandgap(1+ispin)%direct_valence)
            call update_mask(self%bandgap(1+ispin)%direct_conduction)
          enddo !ispin
        endif

! find fundamental gaps
        do ispin=1,self%nspin
! find k-points where vbm and cbm are
            call update_mask(self%bandgap(1)%fundamental_valence)
            call update_mask(self%bandgap(1)%fundamental_conduction)
        enddo

        if (self%nspin==2) then
          do ispin=1,self%nspin
            call update_mask(self%bandgap(1+ispin)%fundamental_valence)
            call update_mask(self%bandgap(1+ispin)%fundamental_conduction)
          enddo !ispin
        endif

! set band_start and band_stop
        self%band_start=1
        self%band_stop=self%nbands
        ispin=1
! set band_start
        do ib=1,size(self%compute_mask,1)
            if (any(self%compute_mask(ib,:,ispin))) exit
        enddo
        self%band_start=ib
! set band_stop
        do ib=size(self%compute_mask,1),1,-1
            if (any(self%compute_mask(ib,:,ispin))) exit
        enddo
        self%band_stop=ib
        contains
        subroutine update_mask(op)
          use bandgap_struct, only: orbital_property
          type(orbital_property),intent(in) :: op
!local variables
          integer :: ib
          do ib=1,self%nbands
! select degenerate bands
              e = real(w%celtot(ib,op%kpoint,op%spin))
              if (abs(e-op%eigenvalue)>EL_TOL) cycle
              self%compute_mask(ib,op%kpoint,op%spin) = .true.
          enddo
        end subroutine
    end subroutine elph_accumulator_selfen_selector_gaps

!> Use transport criteria to select which states need to be computed
!> This goes trough all the tetrahedrons and checks if they contribute to any of the energies
!> Edits the compute_maks, band_start and band_stop of the electron-phonon accumulator structure
    subroutine elph_accumulator_selfen_selector_transport(self,w,kpoints,driver,emin_plot,emax_plot,do_plot,dfermi_tol,transport_nedos, emin_compute, emax_compute)
        use constant, only : bolkev
        use wave_struct_def, only: wavespin
        use mkpoints_struct_def, only : kpoints_struct
        use gauss_quad, only : gauss_legendre
        use elphon_base, only : transport_emin_emax
        use tutor, only: vtutor
        use string, only: str
        type(elph_accumulator_selfen) :: self
        type(wavespin) :: w
        type(kpoints_struct),intent(in) :: kpoints
        integer,intent(in) :: driver !< the type of driver used for transport calculation (1 linear grid and simpson rule, 2 gauss-legendre integration)
        real(q),intent(in) :: emin_plot !< min energy for plotting transport function from the INCAR
        real(q),intent(in) :: emax_plot !< max energy for plotting transport function from the INCAR
        logical,intent(in) :: do_plot !< wether to plot the transport function
        real(q),intent(in) :: dfermi_tol !< tolerance for the fermi occupations
        integer,intent(in) :: transport_nedos !< number of points for transport computation
        real(q),intent(in) :: emin_compute !< min energy for computing transport function from the INCAR
        real(q),intent(in) :: emax_compute !< max energy for computing transport function from the INCAR
! local variables
        integer :: itemp, itet
        integer :: isp, ik, ib
        integer :: kindexes(4)
        real(q) :: smearing, ewmin, ewmax, etmin, etmax
        real(q) :: eig(4)
        real(q), allocatable :: energies(:)
        real(q), allocatable :: weights(:)
        real(q), allocatable :: x(:)

! default energy range
        self%transport_emin = minval(self%efermi)
        self%transport_emax = maxval(self%efermi)
! if we selected that we want to plot the transport function
! then set transport_emin and transport_emax to values of transport_emin_plot and transport_emin_plot in the INCAR file
        if (do_plot) then
            self%transport_emin=emin_plot
            self%transport_emax=emax_plot
! if ranges were not set in the incar file then use min and max eigenval
            if (self%transport_emin == huge(self%transport_emin)) then
               self%transport_emin = self%emin
            endif
! if ranges were not set in the incar file then use min and max eigenval
            if (self%transport_emax == huge(self%transport_emax)) then
               self%transport_emax = self%emax
            endif
        endif

        self%compute_mask=.false.
        self%transport_nedos = transport_nedos

! readjust transport_emin and transport_emax according to necessary energy windows for transport calculation
        select case(driver)
        case (1)
           if (0.0_q<dfermi_tol .and. dfermi_tol<=1.0_q) then
! find max and min energies that need to be computed based on the smearings and temperatures
               do itemp=1,self%ntemps
                   smearing = bolkev*self%tempsk(itemp)
                   call transport_emin_emax(self%efermi(itemp),smearing,dfermi_tol,emin_compute, emax_compute, ewmin,ewmax)
                   self%transport_emin = min(ewmin,self%transport_emin)
                   self%transport_emax = max(ewmax,self%transport_emax)
                enddo
           else
               call vtutor%error('Invalid value for transport_dfermi_tol. It must be in the interval ]0:1] ')
           endif

        case(2)
           allocate(energies(transport_nedos))
           allocate(weights(transport_nedos))
           allocate(x(transport_nedos))
           call gauss_legendre(-1.0_q,1.0_q,x,weights,transport_nedos,eps=1e-12_q)
! determine max energy window for all the temperatures
           do itemp=1,self%ntemps
               smearing = bolkev*self%tempsk(itemp)
               energies = smearing*log((1+x)/(1-x)) + self%efermi(itemp)
               self%transport_emin = min(self%transport_emin,minval(energies))
               self%transport_emax = max(self%transport_emax,maxval(energies))
           enddo
           deallocate(energies)
           deallocate(weights)
           deallocate(x)
        case default
           call vtutor%error('elph_accumulator_selfen_selector_transport: Transport driver '//str(driver)//' not implemented')
        end select

! check which k-points need to be computed given these energy ranges
        self%band_start = w%wdes%nb_tot
        self%band_stop  = 1
        if (kpoints%ntet==0) then
            call vtutor%alert('The tetrahedron method was not initialized so the transport window selection will not work.')
            self%compute_mask=.true.
        endif
        do itet=1,kpoints%ntet
           kindexes = kpoints%idtet(1:4,itet)
           do isp=1,self%nspin
              do ib=1,self%nbands
! get energies on the 4 corners of the tetrahedron
                 eig(:) = real(w%celtot(ib,kindexes(:),isp),q)
! max and min energy of the tetrahedra
                 etmin = minval(eig)
                 etmax = maxval(eig)
! determine wether this tetrahedron can contibute to the
! energy grid above (compare only min and max value)
                 if (self%transport_emax<etmin.and.self%transport_emax<etmax) cycle ! tetrahedron above max
                 if (self%transport_emin>etmin.and.self%transport_emin>etmax) cycle ! tetrahedron below min
! select k-points
                 do ik=1,4
                    self%compute_mask(ib,kindexes(ik),isp)=.true.
                 enddo
! select bands
                 self%band_start = min(self%band_start,ib)
                 self%band_stop  = max(self%band_stop, ib)
              enddo
           enddo
        enddo
        if (count(self%compute_mask)==0) then
           call vtutor%error('No states were selected for transport computation. &
                             &This happens because no KS states were found in the emin: '//str(self%transport_emin)//&
                             &' emax: '//str(self%transport_emax)//' energy range around the fermi level. &
                             &Change the carrier density, temperature or number of integration points.')
        endif
     end subroutine elph_accumulator_selfen_selector_transport

!> Use information of compute_mask to only store the states for which
!> the FAN and DW self-energy are computed
    subroutine elph_accumulator_selfen_sparsify(self,eigenvalues)
        use ini, only: register_allocate
        type(elph_accumulator_selfen) :: self
        complex(q),intent(in) :: eigenvalues(:,:,:)
! local variables
        integer :: isp, ikibz, iband, ibibz, ispin, ibks
        integer :: nbks, iw, n
        real(q) :: enkibz

! set index of each element
        nbks=0
        allocate(self%bks_idx(self%band_start:self%band_stop,self%nkpoints,self%nspin))
        self%bks_idx=0
        do isp=1,self%nspin
            do ikibz=1,self%nkpoints
                do iband=self%band_start,self%band_stop
                    if (.not.self%compute_mask(iband,ikibz,isp)) cycle
                    nbks = nbks+1
                    self%bks_idx(iband,ikibz,isp) = nbks
                enddo
            enddo
        enddo

! create energy grids
        self%nbks = nbks
        allocate(self%energies(self%nw,self%nbks))
        do ispin=1,self%nspin
            do ikibz=1,self%nkpoints
                do ibibz=self%band_start,self%band_stop
                    enkibz = real(eigenvalues(ibibz,ikibz,ispin),q)
                    ibks = self%bks_idx(ibibz,ikibz,ispin)
                    if (ibks<1) cycle
                    do iw=1,self%nw
                        if (self%energy_grid_mode==0) then
                            self%energies(iw,ibks) = &
                            enkibz + ((real(iw,q)-0.5_q)/self%nw-0.5_q)*self%enwin
                        else
                            self%energies(iw,ibks) = &
                            self%emin + (self%emax-self%emin)/(self%nw-1)*(iw-1)
                        endif
                    enddo ! frequencies
                enddo ! bands
            enddo ! kpoints
        enddo ! spin

! allocate arrays
        allocate(self%selfen_fan(self%ntemps,self%nw,nbks))
        allocate(self%selfen_dw(self%ntemps,nbks))
!$ allocate(self%omp_fan_lock(nbks))
!$ do n=1,nbks
!$    call omp_init_lock(self%omp_fan_lock(n))
!$ enddo
        call register_allocate(0.125_q*STORAGE_SIZE(self%selfen_dw)*SIZE(self%selfen_dw,KIND=qi8)+&
                               0.125_q*STORAGE_SIZE(self%selfen_fan)*SIZE(self%selfen_fan,KIND=qi8)+&
                               0.125_q*STORAGE_SIZE(self%bks_idx)*SIZE(self%bks_idx,KIND=qi8),'elph_accumulator_selfen')
        self%selfen_fan=cmplx(0.0_q,0.0_q,q)
        self%selfen_dw=0.0_q
        if (self%dfan) then
            allocate(self%selfen_dfan(self%ntemps,self%nw,nbks))
            call register_allocate(0.125_q*STORAGE_SIZE(self%selfen_dfan)*SIZE(self%selfen_dfan,KIND=qi8),'elph_accumulator_selfen')
            self%selfen_dfan=0.0_q
        endif
    end subroutine elph_accumulator_selfen_sparsify

!> @brief Write information about the computation of the self-enrgy
    subroutine elph_accumulator_selfen_print(self,naccumulator,iu)
        use string, only: str
        use elphon_base, only: elph_selfen_approx_str
        type(elph_accumulator_selfen),intent(in) :: self
        integer,intent(in) :: naccumulator !< index of the accumulator
        integer,intent(in) :: iu
! local variables
        integer :: ikpt, nkpoints_compute, ib, itemp
        logical,allocatable :: compute_band(:)
        if (iu<0) return
        write(iu,'(A,I3)') 'Electron self-energy accumulator N =',naccumulator
        write(iu,'(A)')    '----------------------------------'
        write(iu,'(A,I8,A,I8,A)') ' Band range:               [',self%band_start,':',self%band_stop,']'
        write(iu,'(A,I6)')        ' Number of bands to sum over:',self%nbands_sum
        write(iu,'(A,A)')         ' Scattering approximation: ', elph_selfen_approx_str(self%approx,verbose=.true.)
        write(iu,'(A,L)')         ' Static self-energy: ', self%approx%static
!write(iu,'(A,F6.3)')      ' Valence band maximum:   ', self%vbm
!write(iu,'(A,F6.3)')      ' Conduction band minimum:', self%cbm
        write(iu,'(A,F6.3)')      ' Complex imaginary shift (delta):',self%delta
        write(iu,'(A)')        ' Chemical potentials mu(T):'
        write(iu,'(2A20)')     'Temperature (K)', 'mu (eV)'
        do itemp=1,self%ntemps
            write(iu, '(2F20.8)') self%tempsk(itemp), self%efermi(itemp)
        enddo
        if (self%transport) then
            write(iu,'(A,F8.3,A,F8.3,A,F8.3,A)') ' Transport energy range:   [',self%transport_emin,':',&
                                                                      self%transport_emax,']  wich corresponds to ',&
                                                                      (self%transport_emax-self%transport_emin), ' eV'
        endif
        if (.not.allocated(self%compute_mask)) return
        nkpoints_compute = 0
        do ikpt=1,self%nkpoints
            if (any(self%compute_mask(:,ikpt,1))) then
                nkpoints_compute = nkpoints_compute + 1
            endif
        enddo
        write(iu,'(A,I6,A,I6,A)') ' Number of selected k-points at which to compute the self-energy: [',&
                   nkpoints_compute, ' / ', self%nkpoints,']'
        write(iu,'(A)') ' Selected bands and k-points at which to compute the self-energy:'

        write(iu,'(A)',advance='no') '      '
        allocate(compute_band(lbound(self%compute_mask,1):ubound(self%compute_mask,1)))
        compute_band = .false.
        do ib=lbound(self%compute_mask,1),ubound(self%compute_mask,1)
            if (.not.any(self%compute_mask(ib,:,1))) cycle
            compute_band(ib) = .true.
            write(iu,'(I4)',advance='no') ib
        enddo
        write(iu,*)

        do ikpt=1,self%nkpoints
            if (.not.any(self%compute_mask(:,ikpt,1))) cycle
            write(iu,'(I6)',advance='no') ikpt
            do ib=lbound(self%compute_mask,1),ubound(self%compute_mask,1)
                if (.not.compute_band(ib)) cycle
                write(iu,'(L4)',advance='no') self%compute_mask(ib,ikpt,1)
            enddo
            write(iu,*)
        enddo
        write(iu,*)
    end subroutine elph_accumulator_selfen_print

!> Write the results from the calculation of the self-energy to the outcar
    subroutine elph_accumulator_selfen_print_result(self,i,celtot,iu)
        use string, only: str
        type(elph_accumulator_selfen) :: self
        integer,intent(in) :: i
        complex(q),intent(in) :: celtot(:,:,:) !< eigenvalues in the IBZ
        integer,intent(in) :: iu
! local variables
        integer :: isp, ikpt, ib
        integer :: ibks, itemp
        real(q) :: refan, imfan, redw, dfan, z

        if (iu<0) return

        WRITE(iu,'(A,I3)') 'Electron self-energy accumulator N=',i
        if (self%nw>1) then
            write(iu,*) 'nw > 1 so the self-energy is not written'
            return
        endif

        do itemp=1,self%ntemps
            write(iu,'(A,F8.0,A)') ' T= ',self%tempsk(itemp),' K'
            if (self%dfan) then
                write(iu,'(3A6,6A12)') 'ispin', 'ikpt', 'iband', 'KS eV', 're(Fan) eV', 'im(Fan) eV', 'DW eV', 'Fan+DW eV', 'Z'
            else
                write(iu,'(3A6,6A12)') 'ispin', 'ikpt', 'iband', 'KS eV', 're(Fan) eV', 'im(Fan) eV', 'DW eV', 'Fan+DW eV'
            endif
            do isp=1,self%nspin
                do ikpt=1,self%nkpoints
                    if (.not.any(self%compute_mask(:,ikpt,isp))) cycle
                    do ib=lbound(self%bks_idx,1),ubound(self%bks_idx,1)
                        ibks = self%bks_idx(ib,ikpt,isp)
                        if (ibks<1) cycle
                        refan = real(self%selfen_fan(itemp,1,ibks),q)
                        imfan = aimag(self%selfen_fan(itemp,1,ibks))
                        redw = self%selfen_dw(itemp,ibks)
                        if (self%dfan) then
                            dfan = self%selfen_dfan(itemp,1,ibks)
                            z = 1.0_q/(1.0_q-dfan)
                            write(iu,'(3I6,6F12.6)') isp,ikpt,ib,&
                                                     real(celtot(ib,ikpt,isp),q),&
                                                     refan, imfan, redw, refan+redw, z
                        else
                            write(iu,'(3I6,5F12.6)') isp,ikpt,ib,&
                                                     real(celtot(ib,ikpt,isp),q),&
                                                     refan, imfan, redw, refan+redw
                        endif
                    enddo
                enddo !kpoints
            enddo !spin
            write(iu,*)
        enddo !itemp
        write(iu,*)
    end subroutine elph_accumulator_selfen_print_result

!> get compressed ibks index from index of the band, kpoint and spin
    function elph_accumulator_selfen_get_ibks(self,ib,ik,ispin) result(ibks)
        use tutor, only: vtutor
        use string, only: str
        type(elph_accumulator_selfen),intent(in) :: self
        integer,intent(in) :: ib !< index of the band
        integer,intent(in) :: ik !< index of the k-point
        integer,intent(in) :: ispin !< index of the spin
        integer :: ibks
        if (ib<self%band_start) call vtutor%bug('ib indx out of bounds',"elphon_accumulators.F",720)
        if (ib>self%band_stop) call vtutor%bug('ib index out of bounds',"elphon_accumulators.F",721)
        ibks = self%bks_idx(ib,ik,ispin)
    end function elph_accumulator_selfen_get_ibks

!> get Fan and Debye Waller self-energy from the index of the band, kpoint and spin
    subroutine elph_accumulator_selfen_get_fan_dw_bks(self,ib,ik,ispin,iw,itemp,e0,fan,dw)
        use tutor, only: vtutor
        use string, only: str
        type(elph_accumulator_selfen),intent(in) :: self
        integer,intent(in) :: ib !< index of the band
        integer,intent(in) :: ik !< index of the k-point
        integer,intent(in) :: ispin !< index of the spin
        integer,intent(in) :: iw !< index of the frequency point
        integer,intent(in) :: itemp !< temperature index
        real(q),intent(out) :: e0
        complex(q),intent(out) :: fan
        real(q),intent(out) :: dw
!local variables
        integer :: ibks
        ibks = elph_accumulator_selfen_get_ibks(self,ib,ik,ispin)
        if (ibks<1) call vtutor%bug('the KS state ib: '//str(ib)//&
                                                ' ik: '//str(ik)//&
                                             ' ispin: '//str(ispin)//' was not selected for computation',"elphon_accumulators.F",743)
        e0 = self%energies(iw,ibks)
        fan = self%selfen_fan(itemp,iw,ibks)
        dw = self%selfen_dw(itemp,ibks)
    end subroutine elph_accumulator_selfen_get_fan_dw_bks

!> If computation of the band gaps was selected, compute them
    subroutine elph_accumulator_selfen_gaps_compute(self)
        use tutor, only: vtutor
        use string, only: str
        type(elph_accumulator_selfen) :: self
! local variables
        integer :: itemp, ispin
        integer :: ik, ibks
        integer :: ikpt_vbm, ikpt_cbm
        integer :: iband_vbm, iband_cbm
        integer :: ispin_vbm, ispin_cbm
        real(q) :: vbm_e0, cbm_e0
        complex(q) :: vbm_fan, cbm_fan
        real(q) :: vbm_dw, cbm_dw
        real(q) :: vbm_qp, cbm_qp
        if (.not.self%report_gaps) return
        allocate(self%direct_gap(size(self%bandgap)))
        allocate(self%fundamental_gap(size(self%bandgap)))
        allocate(self%direct_gap_renorm(self%ntemps,size(self%bandgap)))
        allocate(self%fundamental_gap_renorm(self%ntemps,size(self%bandgap)))
! get direct gaps
        do itemp=1,self%ntemps
! the direct gap is by definition at the same k-point so we
! can get the index from the conduction or the valuence,
! the result must be the same
          ikpt_vbm  = self%bandgap(1)%direct_valence%kpoint
          iband_vbm = self%bandgap(1)%direct_valence%band
          ispin_vbm = self%bandgap(1)%direct_valence%spin
          ikpt_cbm  = self%bandgap(1)%direct_conduction%kpoint
          iband_cbm = self%bandgap(1)%direct_conduction%band
          ispin_cbm = self%bandgap(1)%direct_conduction%spin
! get valence band maximum correction
          call elph_accumulator_selfen_get_fan_dw_bks(self,iband_vbm,ikpt_vbm,ispin_vbm,1,itemp,vbm_e0,vbm_fan,vbm_dw)
! get conduction band minimum correction
          call elph_accumulator_selfen_get_fan_dw_bks(self,iband_cbm,ikpt_cbm,ispin_cbm,1,itemp,cbm_e0,cbm_fan,cbm_dw)
! get direct gap
          self%direct_gap = self%bandgap(1)%direct_conduction%eigenvalue-&
                            self%bandgap(1)%direct_valence%eigenvalue
          self%direct_gap = cbm_e0-vbm_e0
! get renormalization of direct gap
          vbm_qp = vbm_e0+real(vbm_fan,q)+vbm_dw
          cbm_qp = cbm_e0+real(cbm_fan,q)+cbm_dw
          self%direct_gap_renorm(itemp,1) = cbm_qp-vbm_qp
        enddo !itemp

        if (self%nspin==2) then
          do ispin=1,self%nspin
            do itemp=1,self%ntemps
              ikpt_vbm  = self%bandgap(1+ispin)%direct_valence%kpoint
              iband_vbm = self%bandgap(1+ispin)%direct_valence%band
              ispin_vbm = self%bandgap(1+ispin)%direct_valence%spin
              ikpt_cbm  = self%bandgap(1+ispin)%direct_conduction%kpoint
              iband_cbm = self%bandgap(1+ispin)%direct_conduction%band
              ispin_cbm = self%bandgap(1+ispin)%direct_conduction%spin
! get valence band maximum correction
              call elph_accumulator_selfen_get_fan_dw_bks(self,iband_vbm,ikpt_vbm,ispin_vbm,1,itemp,vbm_e0,vbm_fan,vbm_dw)
! get conduction band minimum correction
              call elph_accumulator_selfen_get_fan_dw_bks(self,iband_cbm,ikpt_cbm,ispin_cbm,1,itemp,cbm_e0,cbm_fan,cbm_dw)
! get direct gap
              self%direct_gap = self%bandgap(1+ispin)%direct_conduction%eigenvalue-&
                                self%bandgap(1+ispin)%direct_valence%eigenvalue
              self%direct_gap = cbm_e0-vbm_e0
! get renormalization of direct gap
              vbm_qp = vbm_e0+real(vbm_fan,q)+vbm_dw
              cbm_qp = cbm_e0+real(cbm_fan,q)+cbm_dw
              self%direct_gap_renorm(itemp,1+ispin) = cbm_qp-vbm_qp
            enddo !itemp
          enddo !spin
        endif

! get indirect gaps
        do itemp=1,self%ntemps
          ikpt_vbm  = self%bandgap(1)%fundamental_valence%kpoint
          iband_vbm = self%bandgap(1)%fundamental_valence%band
          ispin_vbm = self%bandgap(1)%fundamental_valence%spin
          ikpt_cbm  = self%bandgap(1)%fundamental_conduction%kpoint
          iband_cbm = self%bandgap(1)%fundamental_conduction%band
          ispin_cbm = self%bandgap(1)%fundamental_conduction%spin
! get valence band maximum correction
          call elph_accumulator_selfen_get_fan_dw_bks(self,iband_vbm,ikpt_vbm,ispin_vbm,1,itemp,vbm_e0,vbm_fan,vbm_dw)
! get conduction band minimum correction
          call elph_accumulator_selfen_get_fan_dw_bks(self,iband_cbm,ikpt_cbm,ispin_cbm,1,itemp,cbm_e0,cbm_fan,cbm_dw)
! get fundamental gap
!self%fundamental_gap = self%bandgap(1)%fundamental_conduction%eigenvalue-&
!                       self%bandgap(1)%fundamental_valence%eigenvalue
          self%fundamental_gap(1) = cbm_e0-vbm_e0
! get renormalization of direct gap
          vbm_qp = vbm_e0+real(vbm_fan,q)+vbm_dw
          cbm_qp = cbm_e0+real(cbm_fan,q)+cbm_dw
          self%fundamental_gap_renorm(itemp,1) = cbm_qp-vbm_qp
        enddo !itemp

        if (self%nspin==2) then
          do ispin=1,self%nspin
            do itemp=1,self%ntemps
              ikpt_vbm  = self%bandgap(1)%fundamental_valence%kpoint
              iband_vbm = self%bandgap(1)%fundamental_valence%band
              ispin_vbm = self%bandgap(1)%fundamental_valence%spin
              ikpt_cbm  = self%bandgap(1)%fundamental_conduction%kpoint
              iband_cbm = self%bandgap(1)%fundamental_conduction%band
              ispin_cbm = self%bandgap(1)%fundamental_conduction%spin
! get valence band maximum correction
              call elph_accumulator_selfen_get_fan_dw_bks(self,iband_vbm,ikpt_vbm,ispin_vbm,1,itemp,vbm_e0,vbm_fan,vbm_dw)
! get conduction band minimum correction
              call elph_accumulator_selfen_get_fan_dw_bks(self,iband_cbm,ikpt_cbm,ispin_cbm,1,itemp,cbm_e0,cbm_fan,cbm_dw)
! get direct gap
!self%fundamental_gap = self%bandgap(1+ispin)%fundamental_conduction%eigenvalue-&
!                       self%bandgap(1+ispin)%fundamental_valence%eigenvalue
              self%fundamental_gap(1+ispin) = cbm_e0-vbm_e0
! get renormalization of direct gap
              vbm_qp = vbm_e0+real(vbm_fan,q)+vbm_dw
              cbm_qp = cbm_e0+real(cbm_fan,q)+cbm_dw
              self%fundamental_gap_renorm(itemp,1+ispin) = cbm_qp-vbm_qp
            enddo !itemp
          enddo !spin
        endif

    end subroutine elph_accumulator_selfen_gaps_compute

!> If computation of the band gaps was selected, compute them
!> and write them to the iu unit
    subroutine elph_accumulator_selfen_gaps_write(self,i,iu)
        use tutor, only: vtutor
        use string, only: str
        type(elph_accumulator_selfen),intent(in) :: self
        integer,intent(in) :: i !< index of the accumulator
        integer,intent(in) :: iu
! local variables
        integer :: itemp, ispin
        integer :: ik, ibks
        integer :: ikpt_vbm, ikpt_cbm
        integer :: iband_vbm, iband_cbm
        integer :: ispin_vbm, ispin_cbm
        real(q) :: vbm_e0, cbm_e0
        complex(q) :: vbm_fan, cbm_fan
        real(q) :: vbm_dw, cbm_dw
        real(q) :: vbm_qp, cbm_qp
        if (.not.self%report_gaps) return
        if (iu<0) return
! get direct gaps
        write(iu,'(A)')
        write(iu,'(A,I3)') 'Electron self-energy accumulator N=',i
        write(iu,'(A)') 'Direct gap'

        if (self%nspin==2) write(iu,'(A)') 'spin independent'
        write(iu,'(4A20)') 'Temperature (K)', 'KS gap (eV)', 'QP gap (eV)', 'KS-QP gap (meV)'
        do itemp=1,self%ntemps
          write(iu,'(4F20.6)') self%tempsk(itemp), self%direct_gap, &
                               self%direct_gap_renorm(itemp,1), (self%direct_gap_renorm(itemp,1)-self%direct_gap)*1000
        enddo !itemp

        if (self%nspin==2) then
          do ispin=1,self%nspin
            if (self%nspin==2) write(iu,'(A,I3)') 'spin component ',ispin
            write(iu,'(4A20)') 'Temperature (K)', 'KS gap (eV)', 'QP gap (eV)', 'KS-QP gap (meV)'
            do itemp=1,self%ntemps
              write(iu,'(4F20.6)') self%tempsk(itemp), self%direct_gap, &
                                   self%direct_gap_renorm(itemp,1+ispin), (self%direct_gap_renorm(itemp,1+ispin)-self%direct_gap)*1000
            enddo !itemp
          enddo !spin
        endif

! get indirect gaps
        write(iu,'(A)')
        write(iu,'(A)') 'Fundamental gap'
        if (self%nspin==2) write(iu,'(A)') 'spin independent'
        write(iu,'(4A20)') 'Temperature (K)', 'KS gap (eV)', 'QP gap (eV)', 'KS-QP gap (meV)'
        do itemp=1,self%ntemps
          write(iu,'(4F20.6)') self%tempsk(itemp), self%fundamental_gap(1), &
                               self%fundamental_gap_renorm(itemp,1), (self%fundamental_gap_renorm(itemp,1)-self%fundamental_gap(1))*1000
        enddo !itemp

        if (self%nspin==2) then
          do ispin=1,self%nspin
            if (self%nspin==2) write(iu,'(A,I3)') 'spin component ',ispin
            write(iu,'(4A20)') 'Temperature (K)', 'KS gap (eV)', 'QP gap (eV)', 'KS-QP gap (meV)'
            do itemp=1,self%ntemps
              write(iu,'(4F20.6)') self%tempsk(itemp), self%fundamental_gap(1+ispin), &
                                   self%fundamental_gap_renorm(itemp,1+ispin), (self%fundamental_gap_renorm(itemp,1+ispin)-self%fundamental_gap(1+ispin))*1000
            enddo !itemp
          enddo !spin
        endif
    end subroutine elph_accumulator_selfen_gaps_write

!> @brief Determine whether we should compute this k and k' combination
    pure logical function elph_accumulator_selfen_shallido(self,ikibz,ispin) result(shallido)
        type(elph_accumulator_selfen),intent(in) :: self
        integer,intent(in) :: ikibz
        integer,intent(in) :: ispin
        shallido = any(self%compute_mask(:,ikibz,ispin))
    end function elph_accumulator_selfen_shallido

!> Accumulate the contribution of of a q-point to the Fan self-energy
!> This is a high lever wrapper for the lorentzian and fan variants
    subroutine elph_accumulator_selfen_fan(self,kpoints,latt_cur,ikibz,ispin,vkptp,vqpt,&
                           eleig_ibz,pheig_ibz,enfbz,band_idx,eph_thz,band_start_k,band_start_kp,&
                           g,qweight,comm,velocity,igrpop)
        use mkpoints, only: kpoints_struct
        use lattice, only: latt
        type(elph_accumulator_selfen),intent(inout) :: self
        type(kpoints_struct),intent(in) :: kpoints !< kpoints structure
        type(latt),intent(in) :: latt_cur !< lattice structure
        integer,intent(in) :: ikibz !< index of the k in the ibz
        integer,intent(in) :: ispin !< spin index
        real(q),intent(in) :: vkptp(3) !< coordinates of k'
        real(q),intent(in) :: vqpt(3) !< coordinates of q
        real(q),intent(in) :: enfbz(:) !< eigenvalues at k' (which is usually in the FBZ)
        integer,intent(in) :: band_idx(:) !< global index of the bands that are present locally
        real(q),intent(in) :: eph_thz(:) !< phonon eigenenergies
        complex(q),intent(in) :: eleig_ibz(:,:) !< electron eigenvalues in the IBZ
        real(q),intent(in) :: pheig_ibz(:,:) !< phonon eigenvalues in the IBZ
        integer,intent(in) :: band_start_k !< first band index of k in the array with g matrix elements
        integer,intent(in) :: band_start_kp !< first band index of k in the array with g matrix elements
        complex(q),intent(in) :: g(:,:,:) ! (band at k, band at k', phonon mode)
        real(q),intent(in) :: qweight
        type(communic) :: comm
        real(q), allocatable, intent(in) :: velocity(:,:,:,:)
        integer, intent(in) :: igrpop(3,3,48)
        if (self%tetrahedron) then
            call elph_accumulator_selfen_fan_tet(self,kpoints,latt_cur,ikibz,ispin,vkptp,vqpt,&
                               eleig_ibz,pheig_ibz,enfbz,band_idx,eph_thz,band_start_k,band_start_kp,g,comm,velocity,igrpop)
        else
            call elph_accumulator_selfen_fan_lorentzian(self,ikibz,ispin,enfbz,band_idx,eph_thz,band_start_k,band_start_kp,g,qweight,comm)
        endif
    end subroutine elph_accumulator_selfen_fan

!> @brief Accumulate the contribution of of a q-point to the Fan self-energy
    subroutine elph_accumulator_selfen_fan_lorentzian(self,ikibz,ispin,enfbz,band_idx,eph_thz,band_start_k,band_start_kp,g,qweight,comm)
        use constant
        use tutor
        use elphon_base, only: fermi_dirac, bose_einstein
        use tutor, only: vtutor
        type(elph_accumulator_selfen),intent(inout) :: self
        integer,intent(in) :: ikibz
        integer,intent(in) :: ispin
        real(q),intent(in) :: enfbz(:)
        integer,intent(in) :: band_idx(:)
        real(q),intent(in) :: eph_thz(self%natoms*3)
        integer,intent(in) :: band_start_k !< first band index of k in the array with g matrix elements
        integer,intent(in) :: band_start_kp !< first band index of k' in the array with g matrix elements
        complex(q),intent(in) :: g(:,:,:) !< (band at k, band at k', phonon mode)
        real(q),intent(in) :: qweight
        type(communic) :: comm !< communicator to paralelize the calculation
! Local variables
        integer :: ibfbz, ibibz, ibks, iw, imode, itemp
        real(q) :: n(self%ntemps)
        real(q) :: n_modes(self%ntemps,self%natoms*3)
        real(q) :: f(self%ntemps)
        real(q) :: temps(self%ntemps)
        real(q) :: eph(self%natoms*3)
        real(q) :: eph_denom(self%natoms*3)
        complex(q) :: fan(self%ntemps,self%nw)
        real(q) :: z(self%ntemps,self%nw)
        complex(q) :: idelta, g2
        real(q) :: enkfbz, w, denommw, denompw

        
        eph = eph_thz*THZTOEV
        if (self%approx%static) then
            eph_denom = 0
        else
            eph_denom = eph
        endif
        temps = self%tempsk*BOLKEV ! convert to eV

        idelta = cmplx(0.0_q,self%delta,q)

! loop over phonon modes
        do imode=1,self%natoms*3
! precompute phonon occupations
            do itemp=1,self%ntemps
                n_modes(itemp,imode) = bose_einstein(eph(imode),temps(itemp))
            enddo
        enddo

! loop over bands at k'
!$omp parallel do default(none) private(ibfbz,ibibz,ibks,itemp,f,n,w,g2,enkfbz,fan,z,denommw,denompw) &
!$omp& shared(comm,g,ikibz,ispin,temps,band_idx,enfbz,self,band_start_k,eph,eph_denom,n_modes,idelta,qweight,band_start_kp)
        do ibfbz=band_start_kp,size(enfbz) !self%nbands_sum
! when ncore/=1 we can still use this comunicator to paralelize the self-energy computation
            if (mod(ibfbz,comm%ncpu)/=comm%node_me-1) cycle
            if (band_idx(ibfbz)>self%nbands_sum) cycle
            enkfbz = enfbz(ibfbz)

! precompute electron occupations
            do itemp=1,self%ntemps
                f(itemp) = fermi_dirac(enkfbz-self%efermi(itemp),temps(itemp))
            enddo

! loop over bands at k
            do ibibz=self%band_start,self%band_stop
                ibks = self%bks_idx(ibibz,ikibz,ispin)
                if (ibks<1) cycle
                fan(:,:) = cmplx(0.0_q,0.0_q,q)
                z(:,:) = 0.0_q

! loop over phonon modes
                do imode=1,self%natoms*3
                    if (abs(eph(imode))<PHON_TOL) cycle
                    n(:) = n_modes(:,imode)

! loop over energies
                    do iw=1,self%nw
                        w = self%energies(iw,ibks)
                        g2 = g(ibibz-band_start_k+1,ibfbz-band_start_kp+1,imode)*&
                       conjg(g(ibibz-band_start_k+1,ibfbz-band_start_kp+1,imode))
! compute for all temperatures
                        fan(:,iw) = fan(:,iw) + &
                                        g2*((n(:)+1.0_q-f(:))/(w-enkfbz-eph_denom(imode)+idelta)+&
                                            (n(:)+f(:)      )/(w-enkfbz+eph_denom(imode)+idelta))*qweight
! calculate contribution to dfan
                        if (.not.self%dfan) cycle
                        denommw = (w-enkfbz-eph_denom(imode))
                        denompw = (w-enkfbz+eph_denom(imode))
                        z(:,iw) = z(:,iw) + &
                                        g2*((n(:)+1.0_q-f(:))*(self%delta**2-denommw**2)/(self%delta**2+denommw**2)**2+&
                                            (n(:)+f(:)      )*(self%delta**2-denompw**2)/(self%delta**2+denompw**2)**2)*qweight
                    enddo ! energies
                enddo ! modes

! add to array
!$ call omp_set_lock(self%omp_fan_lock(ibks))
                self%selfen_fan(:,:,ibks) = self%selfen_fan(:,:,ibks) + fan(:,:)
                if (self%dfan) self%selfen_dfan(:,:,ibks) = self%selfen_dfan(:,:,ibks) + z(:,:)
!$ call omp_unset_lock(self%omp_fan_lock(ibks))
            enddo ! k bands
        enddo ! k' bands
!$omp end parallel do
!if (any(self%selfen_fan/=self%selfen_fan)) call vtutor%error('NAN present in selfen_fan_lorentzian')
        
    end subroutine elph_accumulator_selfen_fan_lorentzian

!> @brief Accumulate the contribution of of a q-point to the Fan self-energy
    subroutine elph_accumulator_selfen_fans(self,ikibz,ispin,eph_thz,g,gs,qweight)
        use constant
        use tutor
        use elphon_base, only: fermi_dirac, bose_einstein
        use tutor, only: vtutor
        use ieee_arithmetic, only: ieee_is_nan
        type(elph_accumulator_selfen),intent(inout) :: self
        integer,intent(in) :: ikibz
        integer,intent(in) :: ispin
        real(q),intent(in) :: eph_thz(self%natoms*3)
        complex(q),intent(in) :: g(:,:,:) ! (band at k, band at k', phonon mode)
        complex(q),intent(in) :: gs(:,:,:) ! (band at k, band at k', phonon mode)
        real(q),intent(in) :: qweight
! Local variables
        integer :: ibibz, iw, imode, itemp
        integer :: ibks
        real(q) :: n(self%ntemps)
        real(q) :: temps(self%ntemps)
        real(q) :: eph(self%natoms*3)
        complex(q) :: idelta
        real(q) :: w

        
        temps = self%tempsk*BOLKEV ! convert to eV
        eph = eph_thz*THZTOEV

        idelta = cmplx(0.0_q,self%delta,q)
! loop over phonon modes
        do imode=1,self%natoms*3
            if (abs(eph(imode))<PHON_TOL) cycle
! precompute phonon occupations
            do itemp=1,self%ntemps
                n(itemp) = bose_einstein(eph(imode),temps(itemp))
            enddo

! loop over bands at k
            do ibibz=self%band_start,self%band_stop
                ibks = self%bks_idx(ibibz,ikibz,ispin)
                if (ibks<1) cycle

! loop over energies
                do iw=1,self%nw
                    w = self%energies(iw,ibks)
! compute for all temperatures
                    self%selfen_fan(:,iw,ibks) = &
                    self%selfen_fan(:,iw,ibks) + &
                                        gs(ibibz,ibibz,imode)*conjg(g(ibibz,ibibz,imode))*&
                                       (2*n+1.0_q)*qweight
                enddo ! energies
            enddo ! k bands
        enddo ! modes
! if (any(self%selfen_fan/=self%selfen_fan)) call vtutor%error('NAN present in selfen_fan')
        
    end subroutine elph_accumulator_selfen_fans

!> Initialize some arrays needed for using the tetrahedron method
    subroutine elph_accumulator_selfen_init_tetrahedron(self,kpoints)
        use mkpoints, only: kpoints_struct, get_tetdisp_24tet
        type(elph_accumulator_selfen),intent(inout) :: self
        type(kpoints_struct),intent(in) :: kpoints
! local variables
        integer :: ic, itet

! get tetrahedron displacements
        call get_tetdisp_24tet(kpoints%b,self%kcut24)
! get [0,0,0] point of each kcut
        do itet=1,24
            do ic=1,4
                if( all(self%kcut24(:,ic,itet)==[0,0,0])) exit
            enddo
            self%kcut24_origin(itet)=ic
        enddo
! transform kcut
        do itet=1,24
          do ic=1,4
            self%kcut24(:,ic,itet) = matmul(kpoints%ilsnf,self%kcut24(:,ic,itet))
          enddo
        enddo
        self%kcut24_initialized = .true.
    end subroutine elph_accumulator_selfen_init_tetrahedron

!> Get the electron self-energy in a bks array
    subroutine elph_accumulator_selfen_get(self,itemp,iw,self_energy)
        type(elph_accumulator_selfen) :: self
        integer :: itemp
        integer :: iw
        complex(q),allocatable,intent(inout) :: self_energy(:,:,:)
! local variables
        integer :: isp, ikibz, iband, ibks

        self_energy = 0.0_q
        do isp=1,self%nspin
          do ikibz=1,self%nkpoints
            do iband=self%band_start,self%band_stop
              ibks = self%bks_idx(iband,ikibz,isp)
              if (ibks<1) cycle
              self_energy(iband,ikibz,isp)=self%selfen_fan(itemp,iw,ibks)
            enddo
          enddo
        enddo
    end subroutine elph_accumulator_selfen_get

!> @brief Accumulate the contribution of of a q-point to the Fan self-energy
!> using the tetrahedron method
    subroutine elph_accumulator_selfen_fan_tet(self,kpoints,latt_cur,ikibz,ispin,vkptp,vqpt,&
                           eleig_ibz,pheig_ibz,enfbz,band_idx,eph_thz,band_start_k,band_start_kp,g,comm,velocity,igrpop)
        use constant
        use mkpoints, only: kpoints_struct, kpoint_integers, kint_mod, get_tetdisp_24tet, get_indexes
        use tet_macdonald, only: weight_24tetra_delta
        use lattice, only: latt
        use elphon_base, only: fermi_dirac, bose_einstein
        type(elph_accumulator_selfen),intent(inout) :: self
        type(kpoints_struct),intent(in) :: kpoints !< kpoints structure
        type(latt),intent(in) :: latt_cur !< lattice structure
        integer,intent(in) :: ikibz !< index of the k in the ibz
        integer,intent(in) :: ispin !< spin index
        real(q),intent(in) :: vkptp(3) !< coordinates of k'
        real(q),intent(in) :: vqpt(3) !< coordinates of q
        real(q),intent(in) :: enfbz(:)
        integer,intent(in) :: band_idx(:) !< index of the bands
        real(q),intent(in) :: eph_thz(:)
        complex(q),intent(in) :: eleig_ibz(:,:) !< electron eigenvalues in the IBZ
        real(q),intent(in) :: pheig_ibz(:,:) !< phonon eigenvalues in the IBZ
        integer,intent(in) :: band_start_k !< first band index of k in the array with g matrix elements
        integer,intent(in) :: band_start_kp !< first band index of k' in the array with g matrix elements
        complex(q),intent(in) :: g(:,:,:) ! (band at k, band at k', phonon mode)
        type(communic) :: comm !< communicator to paralelize the calculation
        real(q), allocatable, intent(in) :: velocity(:,:,:,:)
        integer, intent(in) :: igrpop(3,3,48)
! local variables
        integer :: kpint(3) !< coordinates of k'
        integer :: qint(3) !< coordinates of q

        call kpoint_integers(kpoints,latt_cur,vkptp,kpint)
        call kpoint_integers(kpoints,latt_cur,vqpt,qint)

        call elph_accumulator_selfen_fan_tet_low(self,kpoints,latt_cur,ikibz,ispin,kpint,qint,&
                           eleig_ibz,pheig_ibz,enfbz,band_idx,eph_thz,band_start_k,band_start_kp,&
                           g,comm,velocity,igrpop)
    end subroutine

!> @brief Accumulate the contribution of of a q-point to the Fan self-energy
!> using the tetrahedron method
    subroutine elph_accumulator_selfen_fan_tet_low(self,kpoints,latt_cur,ikibz,ispin,kpint,qint,&
                           eleig_ibz,pheig_ibz,enfbz,band_idx,eph_thz,band_start_k,band_start_kp,g,comm,velocity,igrpop)
        use constant
        use lattice, only: latt
        use mkpoints, only: kpoints_struct, kpoint_integers, kint_mod, get_tetdisp_24tet, get_indexes, kint_kindex, kpoints_fbz_to_ibz
        use wave_rotate, only: rotate_vkpt_irot
        use tet_macdonald, only: weight_24tetra_delta
        use tutor, only: vtutor
        use elphon_base, only: fermi_dirac, bose_einstein, &
                               elph_scattering_approx_erta_lambda, elph_scattering_approx_erta_tau, &
                               elph_scattering_approx_mrta_lambda, elph_scattering_approx_mrta_tau
        type(elph_accumulator_selfen),intent(inout) :: self
        type(kpoints_struct),intent(in) :: kpoints !< kpoints structure
        type(latt),intent(in) :: latt_cur !< lattice structure
        integer,intent(in) :: ikibz !< index of the k in the ibz
        integer,intent(in) :: ispin !< spin index
        integer,intent(in) :: kpint(3) !< integer coordinates of k'
        integer,intent(in) :: qint(3) !< integer coordinates of q
        real(q),intent(in) :: enfbz(:)
        integer,intent(in) :: band_idx(:) !< index of the bands
        real(q),intent(in) :: eph_thz(:)
        complex(q),intent(in) :: eleig_ibz(:,:) !< electron eigenvalues in the IBZ
        real(q),intent(in) :: pheig_ibz(:,:) !< phonon eigenvalues in the IBZ
        integer,intent(in) :: band_start_k !< first band index of k in the array with g matrix elements
        integer,intent(in) :: band_start_kp !< first band index of k in the array with g matrix elements
        complex(q),intent(in) :: g(:,:,:) ! (band at k, band at k', phonon mode)
        type(communic) :: comm !< communicator to paralelize the calculation
        real(q), allocatable, intent(in) :: velocity(:,:,:,:) ! (nb,ik,isp,idir)
        integer, intent(in) :: igrpop(3,3,48)
! local variables
        integer :: imode, iw, ibfbz, ibibz, ibks, itemp
        integer :: idx_kp(4,24)
        integer :: idx_q(4,24)
        real(q) :: eleig_tetra(4,24)
        real(q) :: pheig_tetra(4,24,self%natoms*3)
        real(q) :: enkfbz
        real(q) :: w(self%nw)
        real(q) :: deltamw(self%nw)
        real(q) :: deltapw(self%nw)
        real(q) :: n(self%ntemps)
        real(q) :: n_modes(self%ntemps,self%natoms*3)
        real(q) :: f(self%ntemps)
        real(q) :: temps(self%ntemps)
        real(q) :: eph(self%natoms*3)
        integer ikpfbz,ikpibz,irotp
        complex(q) :: fan(self%ntemps,self%nw)
        real(q) :: erta(self%ntemps,self%nw)
        real(q), allocatable :: rot_v(:,:)
        real(q) :: vmkp(3),vnk(3),alpha
        real(q) :: vnk2, dot_vnk_vmkp, norm_vmkp, norm_vnk

        
        eph = eph_thz*THZTOEV
        temps = self%tempsk*BOLKEV ! convert to eV

        if (.not.self%kcut24_initialized) then
            call vtutor%bug('Must call elph_accumulator_selfen_init_tetrahedron before elph_accumulator_selfen_fan_tet',"elphon_accumulators.F",1281)
        endif

! velocity at kp in fbz
        if ((self%approx%scattering_approx == elph_scattering_approx_mrta_tau) .or. &
            (self%approx%scattering_approx == elph_scattering_approx_mrta_lambda) .or. &
            (self%approx%scattering_approx == elph_scattering_approx_erta_tau) .or. &
            (self%approx%scattering_approx == elph_scattering_approx_erta_lambda)) then
            ikpfbz=kint_kindex(kpoints,kpint)
            call kpoints_fbz_to_ibz(kpoints,ikpfbz,ikpibz,irotp)
            allocate(rot_v(1:size(enfbz),1:3))
            do ibfbz=1,size(enfbz)
               vmkp=matmul(velocity(band_idx(ibfbz),ikpibz, ispin,:),latt_cur%A(:,:)) ! velocity in the basis of reciprocal lattice vecors
               vmkp=rotate_vkpt_irot(vmkp,igrpop,irotp) ! rotate velocity to fbz point
               rot_v(ibfbz,1:3)=matmul(latt_cur%B,vmkp) ! rotated velocity in cartesian basis
            enddo
        endif

! find the index of the k' point
        call get_indexes(kpoints,kpint,self%kcut24,idx_kp)
! find the index of the q  point
        call get_indexes(kpoints,qint,-self%kcut24,idx_q)

! loop over phonon modes
        do imode=1,self%natoms*3
! get the phonon frequencies on the 24 tetrahedron for q
            call get_energies(idx_q,imode,pheig_ibz,pheig_tetra(:,:,imode))

! precompute phonon occupations
            do itemp=1,self%ntemps
                n_modes(itemp,imode) = bose_einstein(eph(imode),temps(itemp))
            enddo
        enddo
        pheig_tetra=THZTOEV*pheig_tetra

! loop over bands at k'
!$omp parallel do default(none) private(ibfbz,ibibz,ibks,itemp,f,n,deltapw,deltamw,w,enkfbz,fan,eleig_tetra,alpha,vnk,vmkp,norm_vnk,norm_vmkp,dot_vnk_vmkp,vnk2,erta) &
!$omp& shared(comm,g,ikibz,ispin,temps,band_idx,enfbz,self,eleig_ibz,band_start_k,kpoints,pheig_tetra,eph,n_modes,idx_kp,velocity,rot_v,band_start_kp)
        do ibfbz=band_start_kp,size(enfbz) !self%nbands_sum
! when ncore/=1 we can still use this comunicator to paralelize the self-energy computation
            if (mod(ibfbz,comm%ncpu)/=comm%node_me-1) cycle
            if (band_idx(ibfbz)>self%nbands_sum) cycle
            enkfbz = enfbz(ibfbz)

! get the electron energies on the 24 corners for k'
            call get_energies_complex(idx_kp,ibfbz,eleig_ibz,eleig_tetra)

! precompute electron occupations
            do itemp=1,self%ntemps
                f(itemp) = fermi_dirac(enkfbz-self%efermi(itemp),temps(itemp))
            enddo

! loop over bands at k
            do ibibz=self%band_start,self%band_stop
                ibks=self%bks_idx(ibibz,ikibz,ispin)
                if (ibks<1) cycle
                w = self%energies(:,ibks)
                fan(:,:) = cmplx(0.0_q,0.0_q,q)

! compute the weights for all the energies in the adiabatic approx.
                if (self%approx%static) then
                   deltamw = weight_24tetra_delta(kpoints%volwgt,self%kcut24_origin,self%nw,w,eleig_tetra)
                   deltapw = deltamw
                endif

! loop over phonon modes
                do imode=1,self%natoms*3
                   if (abs(eph(imode))<PHON_TOL) cycle
                   n = n_modes(:,imode)

! compute the weights for all the energies in non-adiabatic case
                   if (.not.self%approx%static) then
!\delta[w-(e_{mk+q}-\omega_{q\nu})] = \delta[w-e_{mk+q}+\omega_{q\nu})
                       deltapw = weight_24tetra_delta(kpoints%volwgt,self%kcut24_origin,self%nw,w,eleig_tetra-pheig_tetra(:,:,imode))
!\delta[w-(e_{mk+q}+\omega_{q\nu})] = \delta[w-e_{mk+q}-\omega_{q\nu})
                       deltamw = weight_24tetra_delta(kpoints%volwgt,self%kcut24_origin,self%nw,w,eleig_tetra+pheig_tetra(:,:,imode))
                   endif

! loop over energies
                   do iw=1,self%nw
! compute for all temperatures
                       fan(:,iw) = fan(:,iw) +&
                                g(ibibz-band_start_k+1,ibfbz-band_start_kp+1,imode)*&
                          conjg(g(ibibz-band_start_k+1,ibfbz-band_start_kp+1,imode))*&
                               ((n(:)+1.0_q-f(:))*cmplx(0.0_q,-pi,q)*deltamw(iw)+&
                                (n(:)+f(:)      )*cmplx(0.0_q,-pi,q)*deltapw(iw))
                   enddo ! energies
                enddo ! modes
                if ((self%approx%scattering_approx == elph_scattering_approx_mrta_tau) .or. &
                    (self%approx%scattering_approx == elph_scattering_approx_mrta_lambda) .or. &
                    (self%approx%scattering_approx == elph_scattering_approx_erta_tau) .or. &
                    (self%approx%scattering_approx == elph_scattering_approx_erta_lambda)) then
                    vnk(:) = velocity(ibibz,ikibz,ispin,:)
                    vmkp(:) = rot_v(ibfbz,:)
                    dot_vnk_vmkp=dot_product(vnk,vmkp)
                    vnk2=dot_product(vnk,vnk)
                    do iw=1,self%nw
                       erta(:,iw) = abs((enkfbz-self%efermi(:))/(w(iw)-self%efermi(:)))
                    enddo
                endif
                if ((self%approx%scattering_approx == elph_scattering_approx_mrta_lambda) .or. &
                    (self%approx%scattering_approx == elph_scattering_approx_erta_lambda)) then
                   norm_vnk=sqrt(vnk2)
                   norm_vmkp=sqrt(dot_product(vmkp,vmkp))
                   alpha = 1.0_q
                   if ((norm_vmkp*norm_vnk)>1e-6) alpha = (1.0_q-dot_vnk_vmkp/(norm_vmkp*norm_vnk))
!alpha = (1.0_q-dot_vnk_vmkp/(norm_vmkp*norm_vnk))
                   if     (self%approx%scattering_approx == elph_scattering_approx_mrta_lambda) then
                      fan(:,:)=fan(:,:)*alpha
                   elseif (self%approx%scattering_approx == elph_scattering_approx_erta_lambda) then
                      fan(:,:)=fan(:,:)*alpha*erta(:,:)
                   endif
                endif
                if ((self%approx%scattering_approx == elph_scattering_approx_mrta_tau) .or. &
                    (self%approx%scattering_approx == elph_scattering_approx_erta_tau)) then
                   alpha = 1.0_q
                   if (vnk2>1e-6) alpha = (1.0_q-dot_vnk_vmkp/vnk2)
!alpha = (1.0_q-dot_vnk_vmkp/vnk2)
                   if     (self%approx%scattering_approx == elph_scattering_approx_mrta_tau) then
                      fan(:,:)=fan(:,:)*alpha
                   elseif (self%approx%scattering_approx == elph_scattering_approx_erta_tau) then
                      fan(:,:)=fan(:,:)*alpha*erta(:,:)
                   endif
                endif
! add to array
!$ call omp_set_lock(self%omp_fan_lock(ibks))
                self%selfen_fan(:,:,ibks) = self%selfen_fan(:,:,ibks) + fan(:,:)
!$ call omp_unset_lock(self%omp_fan_lock(ibks))
            enddo ! k bands
        enddo ! k' bands
!$omp end parallel do
!if (any(self%selfen_fan/=self%selfen_fan)) call vtutor%error('NAN present in selfen_fan_tet')
        
    end subroutine elph_accumulator_selfen_fan_tet_low

!> Check wether this q-point contributes to the imaginary part of the self-energy
    function elph_accumulator_selfen_fan_contributes(self,kpoints,ikibz,ispin,kpint,qint,&
                           band_idx,eleig_ibz,pheig_ibz) result(contributes)
        use mkpoints, only: kpoints_struct
        type(elph_accumulator_selfen),intent(in) :: self
        type(kpoints_struct),intent(in) :: kpoints !< kpoints structure
        integer,intent(in) :: ikibz !< index of the k in the ibz
        integer,intent(in) :: ispin !< spin index
        integer,intent(in) :: kpint(3) !< integer coordinates of k'
        integer,intent(in) :: qint(3) !< integer coordinates of q
        integer,intent(in) :: band_idx(:) !< global index of the bands that are present locally
        complex(q),intent(in) :: eleig_ibz(:,:) !< electron eigenvalues in the IBZ
        real(q),intent(in) :: pheig_ibz(:,:) !< phonon eigenvalues in the IBZ
        logical contributes
        if (self%tetrahedron) then
            contributes = elph_accumulator_selfen_fan_tet_contributes_low(self,kpoints,ikibz,ispin,kpint,qint,&
                           band_idx,eleig_ibz,pheig_ibz)
        else
            contributes = elph_accumulator_selfen_fan_broadening_contributes_low(self,kpoints,ikibz,ispin,kpint,qint,&
                           band_idx,eleig_ibz,pheig_ibz)
        endif
    end function

!> Check wether this q-point contributes to the imaginary part of the self-energy
!> of this k-point
    function elph_accumulator_selfen_fan_tet_contributes(self,kpoints,latt_cur,ikibz,ispin,vkptp,vqpt,&
                           band_idx,eleig_ibz,pheig_ibz) result(contributes)
        use constant
        use lattice, only: latt
        use mkpoints, only: kpoints_struct, kpoint_integers, get_indexes, kint_mod, kint_kindex
        type(elph_accumulator_selfen),intent(in) :: self
        type(kpoints_struct),intent(in) :: kpoints !< kpoints structure
        type(latt),intent(in) :: latt_cur !< lattice structure
        integer,intent(in) :: ikibz !< index of the k in the ibz
        integer,intent(in) :: ispin !< spin index
        real(q),intent(in) :: vkptp(3) !< coordinates of k'
        real(q),intent(in) :: vqpt(3) !< coordinates of q
        integer,intent(in) :: band_idx(:) !< global index of the bands present locally
        complex(q),intent(in) :: eleig_ibz(:,:) !< electron eigenvalues in the IBZ
        real(q),intent(in) :: pheig_ibz(:,:) !< phonon eigenvalues in the IBZ
        logical contributes
! local variables
        integer :: kpint(3) !< coordinates of k'
        integer :: qint(3) !< coordinates of q

        call kpoint_integers(kpoints,latt_cur,vkptp,kpint)
        call kpoint_integers(kpoints,latt_cur,vqpt,qint)

        contributes = elph_accumulator_selfen_fan_tet_contributes_low(self,kpoints,ikibz,ispin,kpint,qint,&
                           band_idx,eleig_ibz,pheig_ibz)
    end function elph_accumulator_selfen_fan_tet_contributes

!> Check wether this q-point contributes to the imaginary part of the self-energy
!> of this k-point when using the tetrahedron method
    function elph_accumulator_selfen_fan_tet_contributes_low(self,kpoints,ikibz,ispin,kpint,qint,&
                           band_idx,eleig_ibz,pheig_ibz) result(contributes)
        use constant
        use mkpoints, only: kpoints_struct, kpoint_integers, get_indexes, kint_mod, kint_kindex
        type(elph_accumulator_selfen),intent(in) :: self
        type(kpoints_struct),intent(in) :: kpoints !< kpoints structure
        integer,intent(in) :: ikibz !< index of the k in the ibz
        integer,intent(in) :: ispin !< spin index
        integer,intent(in) :: kpint(3) !< integer coordinates of k'
        integer,intent(in) :: qint(3) !< integer coordinates of q
        integer,intent(in) :: band_idx(:) !< global index of the band that is present locally
        complex(q),intent(in) :: eleig_ibz(:,:) !< electron eigenvalues in the IBZ
        real(q),intent(in) :: pheig_ibz(:,:) !< phonon eigenvalues in the IBZ
        logical contributes
! local variables
        integer :: imode, ibibz, ibfbz, ibks
        integer :: idx_kp(4,24)
        integer :: idx_q(4,24)
        real(q) :: eleig_tetra(4,24)
        real(q) :: pheig_tetra(4,24,self%natoms*3)
        real(q) :: w(self%nw)

! find the index of the k' point
        call get_indexes(kpoints,kpint,self%kcut24,idx_kp)
! find the index of the q  point
        call get_indexes(kpoints,qint,-self%kcut24,idx_q)

! loop over phonon modes
        do imode=1,self%natoms*3
! get the phonon frequencies on the 24 tetrahedron for q
            call get_energies(idx_q,imode,pheig_ibz,pheig_tetra(:,:,imode))
        enddo
        pheig_tetra=pheig_tetra*THZTOEV

! loop over bands at k'
        do ibfbz=1,size(eleig_ibz,1) !self%nbands_sum
            if (band_idx(ibfbz)>self%nbands_sum) cycle
! get the electron energies on the 24 corners for k'
            call get_energies_complex(idx_kp,ibfbz,eleig_ibz,eleig_tetra)

! loop over bands at k
            do ibibz=self%band_start,self%band_stop
                ibks = self%bks_idx(ibibz,ikibz,ispin)
                if (ibks<1) cycle
                w = self%energies(:,ibks)
!if (.not.self%compute_mask(ibibz,ikibz,ispin)) cycle
!if (all(w(:)<self%transport_emin.or.&
!             self%transport_emax<w(:))) then
!    cycle
!endif

                do imode=1,self%natoms*3
! check wether the electron_energies+phonon_energies contribute
                    if (tetrahedron_kcontributes(eleig_tetra+pheig_tetra(:,:,imode),w)) then
                        contributes = .true.
                        return
                    endif
! check wether the electron_energies-phonon_energies contribute
                    if (tetrahedron_kcontributes(eleig_tetra-pheig_tetra(:,:,imode),w)) then
                        contributes = .true.
                        return
                    endif
                enddo !phonon modes
            enddo !bands at k
        enddo !bands at k'
        contributes = .false.
    end function elph_accumulator_selfen_fan_tet_contributes_low

!> Check wether this q-point contributes to the imaginary part of the self-energy
!> of this k-point when using lorentzian smearing
    function elph_accumulator_selfen_fan_broadening_contributes_low(self,kpoints,ikibz,ispin,kpint,qint,&
                           band_idx,eleig_ibz,pheig_ibz) result(contributes)
        use constant
        use mkpoints, only: kpoints_struct, kpoint_integers, get_indexes, kint_mod, kint_kindex
        type(elph_accumulator_selfen),intent(in) :: self
        type(kpoints_struct),intent(in) :: kpoints !< kpoints structure
        integer,intent(in) :: ikibz !< index of the k in the ibz
        integer,intent(in) :: ispin !< spin index
        integer,intent(in) :: kpint(3) !< integer coordinates of k'
        integer,intent(in) :: qint(3) !< integer coordinates of q
        integer,intent(in) :: band_idx(:) !< global index of the bands that are present locally
        complex(q),intent(in) :: eleig_ibz(:,:) !< electron eigenvalues in the IBZ
        real(q),intent(in) :: pheig_ibz(:,:) !< phonon eigenvalues in the IBZ
        logical contributes
! local variables
        integer :: imode, ibibz, ibfbz, ibks, ibz_q, ibz_kp
        integer :: qindex, kpkindex
        real(q) :: enfbz
        real(q) :: energy_range
        real(q) :: w(self%nw)
        real(q) :: eph_denom(self%natoms*3)

! early exit if we don't exclude anything
        if (self%broad_tol==0.0_q) then
            contributes = .true.
            return
        endif

! early exit if we exclude everything
        if (self%broad_tol==1.0_q) then
            contributes = .false.
            return
        endif

! find the index of the k' point
        kpkindex = kint_kindex(kpoints,kpint)
        ibz_kp = kpoints%fbz2ibz(1,kpkindex)

! find the index of the q  point
        qindex = kint_kindex(kpoints,qint)
        ibz_q = kpoints%fbz2ibz(1,qindex)

! get frequencies of phonon modes at q
        if (self%approx%static) then
            eph_denom = 0.0_q
        else
            eph_denom = pheig_ibz(:,ibz_q)*THZTOEV
        endif

! we choose wether this delta function contributes based on
! the percentage of the integral of the function that we want to include.
! For lorentzian broadening we have the form:
! f(x) = delta/(delta^2+x^2)
! the integral between -y and y is given by
! int(delta/(delta^2+x^2),(x,-y,y)) = -2.atan(y/delta)
! the total integral y->infty is
! -pi
! We want to find the value of y such that the we obtain a percentage given by
! broad_tol of the total integral
! -2.atan(y/delta) = -pi*(1-broad_tol)
! solving for y yields
! y = delta*tan(pi*(1-broad_tol)/2)
! with broad_tol [0:1[
        energy_range = self%delta*tan(PI*(1.0_q-self%broad_tol)/2.0_q)

! loop over bands at k'
        do ibfbz=1,size(eleig_ibz,1) !self%nbands_sum
            if (band_idx(ibfbz)>self%nbands_sum) cycle
            enfbz = real(eleig_ibz(ibfbz,ibz_kp),q)

! loop over bands at k
            do ibibz=self%band_start,self%band_stop
                ibks = self%bks_idx(ibibz,ikibz,ispin)
                if (ibks<1) cycle
                w = self%energies(:,ibks)

                do imode=1,self%natoms*3
! check wether the electron_energies+phonon_energies contribute
                    if (any(energy_range>abs(w-enfbz-eph_denom(imode)))) then
                        contributes = .true.
                        return
                    endif
! check wether the electron_energies-phonon_energies contribute
                    if (any(energy_range>abs(w-enfbz+eph_denom(imode)))) then
                        contributes = .true.
                        return
                    endif
                enddo !phonon modes
            enddo !bands at k
        enddo !bands at k'
        contributes = .false.
    end function elph_accumulator_selfen_fan_broadening_contributes_low

!> Determine whether a certain k-point contributes to a certain energy
!> using the tetrahedron method. This amounts to looking wether the energy is
!> inside any of the 24 tetrahedra around this k-point
    pure function tetrahedron_kcontributes(energies_tetra,energies)
!> energies at the corners of the 24 tetrahedra around e_k
        real(q),intent(in) :: energies_tetra(4,24)
        real(q),intent(in) :: energies(:)
        logical :: tetrahedron_kcontributes
! local variables
        real(q) :: emin, emax
! get max and min energy of every 4 corners of the 24 tetrahedra
        emax=maxval(energies_tetra(:,:))
        emin=minval(energies_tetra(:,:))
! if the energy sits in between the min and max energies then it is likely
! that the contribution is non-(0._q,0._q)
        tetrahedron_kcontributes = any(emin < energies(:) .and. energies(:) < emax)
    end function tetrahedron_kcontributes

!> @brief Get energies corresponding to the corners of the tetrahedron
    subroutine get_energies(idxs,iband,energies_ibz,energies)
        integer,intent(in) :: idxs(4,24)
        integer,intent(in) :: iband
        real(q),intent(in) :: energies_ibz(:,:)
        real(q),intent(inout) :: energies(4,24) !< energies at the corners of the tetrahedron
! local variables
        integer :: itet,ic,ikibz
        do itet=1,24
            do ic=1,4
                ikibz = idxs(ic,itet)
                energies(ic,itet) = energies_ibz(iband,ikibz)
            enddo
        enddo
    end subroutine get_energies

!> @brief Get energies corresponding to the corners of the tetrahedron
    subroutine get_energies_complex(idxs,iband,energies_ibz,energies)
        integer,intent(in) :: idxs(4,24)
        integer,intent(in) :: iband
        complex(q),intent(in) :: energies_ibz(:,:)
        real(q),intent(inout) :: energies(4,24) !< energies at the corners of the tetrahedron
! local variables
        integer :: itet,ic,ikibz
        do itet=1,24
            do ic=1,4
                ikibz = idxs(ic,itet)
                energies(ic,itet) = real(energies_ibz(iband,ikibz))
            enddo
        enddo
    end subroutine get_energies_complex

!> @brief Accumulate the contribution of a q-point to the Debye Waller energy
!> I follow the notation and conventions in Giustino's RMP (Eq. 195-198)
!> Here I compute the t matrix which is then used to compute the DW self-energy
!>
!> \f(
!> t^\nu_{\kappa\alpha,\kappa'\alpha'}=
!> \int \frac{d\mathbf{q}}{\Omega_\text{BZ}}
!> \frac{1}
!> {2\omega_{\mathbf{q}\nu}}
!> \left[
!> \frac{
!>     e_{\kappa\alpha\nu}(\mathbf{q})
!>     e^*_{\kappa\alpha\nu}(\mathbf{q})
!> }{M_\kappa}+
!> \frac{
!>     e_{\kappa'\alpha\nu}(\mathbf{q})
!>     e^*_{\kappa'\alpha\nu}(\mathbf{q})
!> }{M_\kappa'}
!> \right]
!> (2n_{\mathbf{q}\nu} + 1)
!> \f)
    subroutine elph_accumulator_selfen_dw_t(self,eph_thz,eph_eiv,qweight)
        use constant
        use elphon_base, only: bose_einstein
        type(elph_accumulator_selfen) :: self
        real(q),intent(in) :: eph_thz(self%natoms*3)
        real(q),intent(in) :: qweight
!> phonon eigenvectors
        complex(q),intent(in) :: eph_eiv(3*self%natoms,3*self%natoms)
! Local variables
        integer :: iatom1, iatom2, imode1, imode2
        integer :: idir1, idir2
        integer :: imode, itemp
        complex(q) :: idelta
        real(q) :: n(self%ntemps)
        real(q) :: temps(self%ntemps)
        real(q) :: eph(self%natoms*3)

        
        temps = self%tempsk*BOLKEV ! convert to eV
        eph = eph_thz*THZTOEV
        idelta = cmplx(0.0_q,self%delta,q)

! Compute debye waller factor t
! Loop over phonon modes
!self%t = cmplx(0.0_q,0.0_q,q)
        do imode=1,self%natoms*3
            if (abs(eph(imode))<PHON_TOL) cycle
! precompute phonon occupations for different temperature
            do itemp=1,self%ntemps
                n(itemp) = bose_einstein(eph(imode),temps(itemp))
            enddo

            do iatom1=1,self%natoms
            do idir1=1,3
                imode1=(iatom1-1)*3+idir1
                do iatom2=1,self%natoms
                do idir2=1,3
                    imode2=(iatom2-1)*3+idir2
! note that eph_eiv is scalled by 1/sqrt(2*m*w)
                    self%t(:,idir1,iatom1,idir2,iatom2,imode) = &
                    self%t(:,idir1,iatom1,idir2,iatom2,imode)+&
                         (eph_eiv((iatom1-1)*3+idir1,imode)*&
                    conjg(eph_eiv((iatom1-1)*3+idir2,imode))+&
                          eph_eiv((iatom2-1)*3+idir1,imode)*&
                    conjg(eph_eiv((iatom2-1)*3+idir2,imode)))*&
                          (2.0_q*n+1.0_q)*qweight
                enddo
                enddo
            enddo
            enddo
        enddo

        
    end subroutine elph_accumulator_selfen_dw_t

    subroutine elph_accumulator_selfen_dw_sum_t(self,comm_batch)
       use mpimy
       type(elph_accumulator_selfen) :: self
       type(communic) :: comm_batch
       CALL m_sum_z( comm_batch, self%t, size(self%t))
    end subroutine elph_accumulator_selfen_dw_sum_t

!> @brief Compute Debye Waller matrix elements
!> I follow the notation and conventions in Giustino's RMP (Eq. 195-198)
!> Here these two formulas are implemented:
!>
!> \f(
!> \Sigma^\text{DW}_{nn\mathbf{k}}=
!> -\sum'_{m}
!> \frac{g^{2,\text{DW}}_{mn}(\mathbf{k})}
!> {{e_{n\mathbf{k}}-e_{m\mathbf{k}}}}
!> \f)
!>
!> With the second-order electron-phonon matrix element
!> expressed in terms of first order matrix elements:
!>
!> \f(
!> g^{2,\text{DW}}_{mn}(\mathbf{k})=
!> \sum_{\kappa\alpha,\kappa'\alpha',\nu}
!> h^*_{mn,\kappa\alpha}(\mathbf{k})
!> h_{mn,\kappa'\alpha'}(\mathbf{k})
!> t^\nu_{\kappa\alpha,\kappa'\alpha'}
!> \f)
    subroutine elph_accumulator_selfen_dw_g(self,ikibz,ispin,ek,ekp,ekp_idx,band_start_k,band_start_kp,g_atom_gamma,comm)
        type(elph_accumulator_selfen),target :: self
        integer,intent(in) :: ispin
        real(q),intent(in) :: ek(:) !< eigenvalues for k
        real(q),intent(in) :: ekp(:) !< eigenvalues for k'
        integer,intent(in) :: ekp_idx(:) !< number of the bands at k'
        integer,intent(in) :: band_start_k !< first band index of k in the array with g matrix elements
        integer,intent(in) :: band_start_kp !< first band index of k' in the array with g matrix elements
        complex(q),intent(in) :: g_atom_gamma(:,:,:,:)
        type(communic) :: comm !< communicator to paralelize the calculation
! Local variables
        integer :: iatom1, iatom2
        integer :: ikibz
        integer :: imode
        integer :: ibmk, ibnk
        integer :: idir1, idir2
        integer :: ibks
        real(q) :: enk, emk
        complex(q) :: ediff
        complex(q) :: idelta
        complex(q) :: dw2(self%ntemps,self%band_start:self%band_stop,size(ekp))
        idelta = cmplx(0.0_q,self%delta,q)
! TODO: make sure the users are aware of this
        if (self%delta==0.0_q) idelta = cmplx(0.0_q,0.01,q)

        

        
! Compute debye waller matrix elements
        dw2=cmplx(0.0_q,0.0_q,q)
        do imode=1,self%natoms*3
        do iatom2=1,self%natoms
        do idir2=1,3
            do iatom1=1,self%natoms
            do idir1=1,3
                do ibmk=band_start_kp,size(ekp) !self%nbands_sum
! when ncore/=1 we can still use this comunicator to paralelize the self-energy computation
                if (mod(ibmk,comm%ncpu)/=comm%node_me-1) cycle
                if (ekp_idx(ibmk)>self%nbands_sum) cycle
                do ibnk=self%band_start,self%band_stop
                    dw2(:,ibnk,ibmk) = dw2(:,ibnk,ibmk)+&
                                      self%t(:,idir1,iatom1,idir2,iatom2,imode)/2.0*&
                                      g_atom_gamma(ibnk-band_start_k+1,ibmk-band_start_kp+1,idir1,iatom1)*&
                                conjg(g_atom_gamma(ibnk-band_start_k+1,ibmk-band_start_kp+1,idir2,iatom2))
                enddo ! ibnk
                enddo ! ibmk
            enddo ! idir2
            enddo ! iatom2
        enddo ! idir1
        enddo ! iatom1
        enddo ! imode
        

        
! loop over nk bands
        do ibnk=self%band_start,self%band_stop
            ibks=self%bks_idx(ibnk,ikibz,ispin)
            if (ibks<1) cycle
            enk=ek(ibnk)
! loop over mk bands
            do ibmk=1,size(ekp) !self%nbands_sum
! when ncore/=1 we can still use this comunicator to paralelize the self-energy computation
                if (mod(ibmk,comm%ncpu)/=comm%node_me-1) cycle
                if (ekp_idx(ibmk)>self%nbands_sum) cycle
                emk=ekp(ibmk)
                ediff = enk-emk+idelta
!if(abs(ediff)<EL_TOL) cycle
                self%selfen_dw(:,ibks) = &
                self%selfen_dw(:,ibks)&
                -real(dw2(:,ibnk,ibmk)/ediff,q)
            enddo ! k' bands
        enddo ! k bands
        
        
    end subroutine elph_accumulator_selfen_dw_g


!> Compute additional terms contributing to the electron-phonon self-energy in the static approximation
    subroutine elph_accumulator_selfen_dws_g(self,ikibz,ispin,band_start_k,g_atom_gamma,gs_atom_gamma)
        type(elph_accumulator_selfen),target :: self
        integer,intent(in) :: ispin
        integer,intent(in) :: band_start_k !< first band index of k in the array with g matrix elements
        complex(q),intent(in) :: g_atom_gamma(:,:,:,:)
        complex(q),intent(in) :: gs_atom_gamma(:,:,:,:)
! Local variables
        integer :: iatom1, iatom2
        integer :: ikibz
        integer :: imode
        integer :: ibmk, ibnk
        integer :: idir1, idir2
        integer :: ibks
        complex(q) :: idelta
        complex(q) :: dw2(self%ntemps,self%band_start:self%band_stop,self%nbands)
        idelta = cmplx(0.0_q,self%delta,q)
! TODO: make sure the users are aware of this
        if (self%delta==0.0_q) idelta = cmplx(0.0_q,0.01,q)

        

! Compute debye waller matrix elements
        dw2=cmplx(0.0_q,0.0_q,q)
        do iatom1=1,self%natoms
        do idir1=1,3
            do iatom2=1,self%natoms
            do idir2=1,3
                do ibmk=self%band_start,self%band_stop
                do ibnk=self%band_start,self%band_stop
                do imode=1,self%natoms*3
                    dw2(:,ibnk,ibmk) = dw2(:,ibnk,ibmk)+&
                                      self%t(:,idir1,iatom1,idir2,iatom2,imode)/2.0*&
                                      gs_atom_gamma(ibnk-band_start_k+1,ibmk,idir1,iatom1)*&
                                      g_atom_gamma(ibnk-band_start_k+1,ibmk,idir2,iatom2)
                enddo ! imode
                enddo ! ibnk
                enddo ! ibmk
            enddo ! idir2
            enddo ! iatom2
        enddo ! idir1
        enddo ! iatom1

! loop over nk bands
        do ibnk=self%band_start,self%band_stop
            ibks=self%bks_idx(ibnk,ikibz,ispin)
            if (ibks<1) cycle
            self%selfen_dw(:,ibks) = &
            self%selfen_dw(:,ibks)&
            +real(dw2(:,ibnk,ibnk))
        enddo ! k bands
        
    end subroutine elph_accumulator_selfen_dws_g

!> @brief Sum results across the nodes
    subroutine elph_accumulator_selfen_sum(self,comm_batch)
       use mpimy
       type(elph_accumulator_selfen) :: self
       type(communic) :: comm_batch
       
       CALL m_sum_z( comm_batch, self%selfen_fan, size(self%selfen_fan))
       CALL m_sum_d( comm_batch, self%selfen_dw, size(self%selfen_dw))
       if (self%dfan) then
            CALL m_sum_d( comm_batch, self%selfen_dfan, size(self%selfen_dfan))
       endif
       
    end subroutine elph_accumulator_selfen_sum

!> Compute the real part of the self-energy using the Kramers-Konig relation
    subroutine elph_accumulator_selfen_kk(self)
       use constant, only: pi
       type(elph_accumulator_selfen) :: self
! local variables
       integer :: i,j,ibks,itemp
       real(q) :: x(self%nw)
       real(q) :: fan_re(self%nw), fan_im(self%nw)

! only use Kramers-Konig in the tetrahedron case
       if (.not.self%tetrahedron) return
       if (self%nw == 1) return
       
       do ibks=1,self%nbks
           x = self%energies(:,ibks)
           do itemp=1,self%ntemps
               fan_re=0.0d0
               fan_im=aimag(self%selfen_fan(itemp,:,ibks))
               do i=2,self%nw-1
                  do j=1,self%nw-1
                     if (j==i-1) cycle
                     if (j==i) cycle
                     fan_re(i)=fan_re(i)+(fan_im(j)/(x(i)-x(j))+fan_im(j+1)/(x(i)-x(j+1)))*(x(j+1)-x(j))/2
                  enddo
               enddo
               self%selfen_fan(itemp,:,ibks) = -fan_re/pi+(0.0_q,1.0_q)*fan_im
           enddo !itemp
       enddo
       
    end subroutine elph_accumulator_selfen_kk

!> find degenerate eigenvalues and average the FAN and DW energies
    subroutine elph_accumulator_selfen_degavg(self,w)
       use wave_struct_def, only: wavespin
       type(elph_accumulator_selfen) :: self
       type(wavespin),intent(in) :: w
! local variables
       integer :: ib,ik,ispin
       integer :: ibstart,ibstop,ibks
       integer :: ndeg
       complex(q) :: fan_avg(self%ntemps,self%nw)
       complex(q) :: dfan_avg(self%ntemps,self%nw)
       complex(q) :: dw_avg(self%ntemps)
       
       if (self%dfan) dfan_avg = 0
! find degenerate states
       do ispin=1,self%nspin
           do ik=1,self%nkpoints
               ib = self%band_start
               do while (ib<self%band_stop)
                   ibstart = ib
                   ibstop  = ib
                   do while( abs(w%celtot(ibstop+1,ik,ispin)-w%celtot(ibstart,ik,ispin))<EL_TOL )
                      ibstop=ibstop+1
                      if (ibstop>=self%band_stop) exit
                   enddo

! debug
!write(*,'(a,i3,a)',advance='no') 'ndeg:',ndeg, ' nbands:'
!do ib=ibstart,ibstop
!   write(*,'(i3)',advance='no') ib
!enddo
!write(*,*)

! sum all the elements
                   fan_avg=0
                   dw_avg=0
                   ndeg = 0
                   do ib=ibstart,ibstop
                       ibks=self%bks_idx(ib,ik,ispin)
                       if (ibks<1) cycle
                       ndeg=ndeg+1
                       fan_avg(:,:) = fan_avg(:,:) + self%selfen_fan(:,:,ibks)
                       dw_avg(:) = dw_avg(:) + self%selfen_dw(:,ibks)
                       if (self%dfan) dfan_avg(:,:) = dfan_avg(:,:) + self%selfen_dfan(:,:,ibks)
                   enddo
! and set them
                   do ib=ibstart,ibstop
                       ibks=self%bks_idx(ib,ik,ispin)
                       if (ibks<1) cycle
                       self%selfen_fan(:,:,ibks) = fan_avg(:,:)/ndeg
                       self%selfen_dw(:,ibks) = dw_avg(:)/ndeg
                       if (self%dfan) self%selfen_dfan(:,:,ibks) = dfan_avg(:,:)/ndeg
                   enddo
! move to next eigenvalue
                   ib=ibstop+1
               enddo
           enddo
       enddo
       
    end subroutine elph_accumulator_selfen_degavg


!> @brief Write the information of the electron self-energy accumulator to a file
    subroutine elph_accumulator_selfen_hdf5write(self,fileid,idx,group,subgroup)
        use vhdf5
        use string, only: str
        use elphon_base, only: elph_selfen_approx_hdf5write
        type(elph_accumulator_selfen) :: self
        integer,intent(in) :: idx
        integer(HID_T),intent(in) :: fileid
        character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
        character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_DOS)
! local variables
        integer(HID_T) :: groupid, subgroupid
        character(len=MAX_LEN_GROUP) :: my_group
        character(len=MAX_LEN_GROUP) :: my_subgroup

        my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
        my_subgroup = GRP_ELPHON;   if (present(subgroup)) my_subgroup = trim(subgroup)

! open group
        call vh5_error(vh5_group_open_or_create(fileid, trim(my_group)//'/'//trim(my_subgroup)//'/electrons', groupid),"elphon_accumulators.F",2044)
        call vh5_error(vh5_group_open_or_create(groupid, 'self_energy_'//str(idx), subgroupid),"elphon_accumulators.F",2045)

! write important metadata
        call vh5_error(vh5_write(subgroupid,'nbands',self%nbands),"elphon_accumulators.F",2048)
        call vh5_error(vh5_write(subgroupid,'nbands_sum',self%nbands_sum),"elphon_accumulators.F",2049)
        call vh5_error(vh5_write(subgroupid,'delta',self%delta),"elphon_accumulators.F",2050)
        call vh5_error(vh5_write(subgroupid,'band_start',self%band_start),"elphon_accumulators.F",2051)
        call vh5_error(vh5_write(subgroupid,'band_stop',self%band_stop),"elphon_accumulators.F",2052)
        call vh5_error(vh5_write(subgroupid,'nw',self%nw),"elphon_accumulators.F",2053)
        call vh5_error(vh5_write(subgroupid,'enwin',self%enwin),"elphon_accumulators.F",2054)
        call vh5_error(vh5_write(subgroupid,'efermi',self%efermi),"elphon_accumulators.F",2055)
        call vh5_error(vh5_write(subgroupid,'carrier_per_cell',self%carrier_per_cell),"elphon_accumulators.F",2056)
        call vh5_error(vh5_write(subgroupid,'tetrahedron',self%tetrahedron),"elphon_accumulators.F",2057)
        call vh5_error(vh5_write(subgroupid,'select_energy_window',self%select_energy_window),"elphon_accumulators.F",2058)
!call vh5_error(vh5_write(subgroupid,'compute_ikpt',any(self%compute_mask,1)),"elphon_accumulators.F",2059)
!call vh5_error(vh5_write(subgroupid,'compute_mask',self%compute_mask),"elphon_accumulators.F",2060)
        call vh5_error(vh5_write(subgroupid,'bks_idx',self%bks_idx),"elphon_accumulators.F",2061)

        call vh5_error(vh5_write(subgroupid,'id_idx',self%id_idx),"elphon_accumulators.F",2063)
        call vh5_error(vh5_write(subgroupid,'id_name',self%id_name),"elphon_accumulators.F",2064)

        call elph_selfen_approx_hdf5write(self%approx,subgroupid)

! write bandgap renormalization to hdf5 file
        if (self%report_gaps) then
            call vh5_error(vh5_write(subgroupid,'fundamental_gap',self%fundamental_gap),"elphon_accumulators.F",2070)
            call vh5_error(vh5_write(subgroupid,'direct_gap',self%direct_gap),"elphon_accumulators.F",2071)
            call vh5_error(vh5_write(subgroupid,'fundamental_gap_renorm',self%fundamental_gap_renorm),"elphon_accumulators.F",2072)
            call vh5_error(vh5_write(subgroupid,'direct_gap_renorm',self%direct_gap_renorm),"elphon_accumulators.F",2073)
        endif

! write main data
        call vh5_error(vh5_write(subgroupid,'selfen_fan',self%selfen_fan),"elphon_accumulators.F",2077)
        if (self%dfan) call vh5_error(vh5_write(subgroupid,'selfen_dfan',self%selfen_dfan),"elphon_accumulators.F",2078)
        call vh5_error(vh5_write(subgroupid,'selfen_dw',self%selfen_dw),"elphon_accumulators.F",2079)
        call vh5_error(vh5_write(subgroupid,'temps',self%tempsk),"elphon_accumulators.F",2080)
        call vh5_error(vh5_write(subgroupid,'energies',self%energies),"elphon_accumulators.F",2081)

! close group
        call vh5_error(vh5_group_close_writing(subgroupid),"elphon_accumulators.F",2084)
        call vh5_error(vh5_group_close_writing(groupid),"elphon_accumulators.F",2085)
    end subroutine elph_accumulator_selfen_hdf5write

!> @brief Write the information of the electron self-energy accumulator to a file
    subroutine elph_accumulators_selfen_hdf5write(self,fileid,group)
        use vhdf5
        use string, only: str
        use elphon_base, only: elph_selfen_approx_hdf5write
        type(elph_accumulators_selfen) :: self
        integer(HID_T),intent(in) :: fileid
        character(len=*), intent(in), optional :: group
!local variables
        integer(HID_T) :: groupid
        character(len=MAX_LEN_GROUP) :: my_group

        my_group    = GRP_RESULTS//'/'//GRP_ELPHON//'/electrons/self_energy_meta'
        if (present(group))    my_group    = trim(group)

! open group
        call vh5_error(vh5_group_open_or_create(fileid, my_group, groupid),"elphon_accumulators.F",2104)

        call vh5_error(vh5_write(groupid,'ncalculators',product(self%id_size)),"elphon_accumulators.F",2106)
        call vh5_error(vh5_write(groupid,'id_size',self%id_size),"elphon_accumulators.F",2107)
        call vh5_error(vh5_write(groupid,'id_name',self%id_name),"elphon_accumulators.F",2108)

! close group
        call vh5_error(vh5_group_close_writing(groupid),"elphon_accumulators.F",2111)
    end subroutine elph_accumulators_selfen_hdf5write


!> @brief Free the accumulator structure
    subroutine elph_accumulator_selfen_free(self)
        use ini, only: deregister_allocate
        type(elph_accumulator_selfen) :: self
        deallocate(self%tempsk)
        deallocate(self%energies)
        deallocate(self%efermi)
        deallocate(self%carrier_per_cell)

        call deregister_allocate(0.125_q*STORAGE_SIZE(self%selfen_dw)*SIZE(self%selfen_dw,KIND=qi8)+&
                                 0.125_q*STORAGE_SIZE(self%selfen_fan)*SIZE(self%selfen_fan,KIND=qi8)+&
                                 0.125_q*STORAGE_SIZE(self%bks_idx)*SIZE(self%bks_idx,KIND=qi8)+&
                                 0.125_q*STORAGE_SIZE(self%compute_mask)*SIZE(self%compute_mask,KIND=qi8),'elph_accumulator_selfen')
        deallocate(self%selfen_fan)
!$ deallocate(self%omp_fan_lock)
        deallocate(self%selfen_dw)
        deallocate(self%compute_mask)
        deallocate(self%bks_idx)
    end subroutine elph_accumulator_selfen_free

! from here on its about routines that deal with a set of accumulators
    subroutine elph_accumulators_init(self,elph_set,wdes)
        use elphon_base, only: elph_settings, wdes_get_global_band_index
        use wave_struct_def, only: wavedes
        type(elph_accumulators_selfen) :: self
        type(elph_settings),intent(in) :: elph_set
        type(wavedes),intent(in) :: wdes
! local variables
        self%naccumulators_el = elph_set%selfen_naccumulators
        allocate(self%selfen_el(self%naccumulators_el))
! keep track of the index of the bands
        call wdes_get_global_band_index(wdes,self%band_idx)
    end subroutine elph_accumulators_init

    subroutine elph_accumulators_free(self)
        type(elph_accumulators_selfen) :: self
! local variables
        integer :: i
        do i=1,self%naccumulators_el
            call elph_accumulator_selfen_free(self%selfen_el(i))
        enddo
    end subroutine elph_accumulators_free

!> checck if the k, q and kp scattering event contributes to the any of the self-energy accumulators
    function elph_accumulators_contributes(self,elph_set,w_dense_ibz,phonon_freqs_ibz,band_idx,kpoints,&
                                           triplet,isp) result(contributes)
        use elphon_base, only: elph_settings
        use lattice, only: latt
        use wave_struct_def, only: wavespin
        use mkpoints_struct_def, only: kpoints_struct
        use mkpoints, only: kindex_vkpt, ikibz_kindex
        use triplets, only: triplet_type
        type(elph_accumulators_selfen) :: self
        type(elph_settings) :: elph_set
        type(wavespin) :: w_dense_ibz
        real(q) :: phonon_freqs_ibz(:,:)
        integer,intent(in) :: band_idx(:)
        type(kpoints_struct) :: kpoints
        type(triplet_type) :: triplet
        integer :: isp
        logical :: contributes

        contributes = .true.
        if (.not.elph_accumulators_shallido(self,triplet%k%ibz,isp)) then
            contributes = .false.
            return
        endif
        if (elph_accumulators_delta_is_zero(self,elph_set,w_dense_ibz,phonon_freqs_ibz,&
                                            band_idx,kpoints,triplet%k%ibz,isp,triplet%kp%kint,triplet%q%kint)) then
            contributes = .false.
            return
        endif
    end function elph_accumulators_contributes

    function elph_accumulators_shallido(self,ik_ibz,isp) result(shallido)
        type(elph_accumulators_selfen) :: self
        integer,intent(in) :: ik_ibz
        integer,intent(in) :: isp
        logical :: shallido
! local variables
        integer :: i
        shallido=.false.
! first check is wether this k-point is marked for computation
        do i=1,self%naccumulators_el
           if (elph_accumulator_selfen_shallido(self%selfen_el(i),ik_ibz,isp)) then
              shallido=.true.
              return
           endif
        enddo
    end function elph_accumulators_shallido

! Check if any of the delta functions for the accumulators is not (0._q,0._q)
    function elph_accumulators_delta_is_zero(self,elph_set,w_dense_ibz,phonon_freqs_ibz,band_idx,kpoints_inter,&
                                             ik_ibz,isp,kpint,qint) result(delta_is_zero)
        use elphon_base, only: elph_settings
        use wave_struct_def, only: wavespin
        use mkpoints_struct_def, only: kpoints_struct
        type(elph_accumulators_selfen),intent(in) :: self
        type(elph_settings),intent(in) :: elph_set
        type(wavespin),intent(in) :: w_dense_ibz
        integer,intent(in) :: band_idx(:)
        real(q),intent(in) :: phonon_freqs_ibz(:,:)
        type(kpoints_struct),intent(in) :: kpoints_inter
        integer,intent(in) :: ik_ibz
        integer,intent(in) :: isp
        integer,intent(in) :: kpint(3)
        integer,intent(in) :: qint(3)
        logical :: delta_is_zero
! local variables
        integer :: i
        delta_is_zero = .false.
! next check is specific for transport computations
        if (elph_set%selfen_imag_skip) then
           do i=1,self%naccumulators_el
              if (elph_accumulator_selfen_fan_contributes(self%selfen_el(i),kpoints_inter,&
                         ik_ibz,isp,kpint,qint,band_idx,w_dense_ibz%celen(:,:,isp),phonon_freqs_ibz)) then
                 return
              endif
           enddo
           delta_is_zero = .true.
        endif
    end function elph_accumulators_delta_is_zero

!> Receive the matrix elements and loop over all the accumulators storing what is required
    subroutine elph_accumulators_accumulate(self,elph_set,band_start_k,band_start_kp,triplet,isp,&
                                            kpoints,waverot,w_dense_ibz,wdes_dense_ibz,phonon_freqs_ibz,&
                                            latt_cur,gmode_kkp,gcart_kkp,velocity)
        use lattice, only: latt
        use elphon_base, only: elph_settings
        use mkpoints_struct_def, only: kpoints_struct
        use mkpoints, only: kindex_vkpt, kpoints_fbz_to_ibz
        use wave_struct_def, only: wavespin, wavedes
        use wave_rotate, only: wave_rotator, rotate_vkpt_irot
        use triplets, only: triplet_type
        type(elph_accumulators_selfen) :: self
        type(elph_settings) :: elph_set
        integer,intent(in) :: band_start_k
        integer,intent(in) :: band_start_kp
        type(triplet_type),intent(in) :: triplet
        integer,intent(in) :: isp
        type(kpoints_struct) :: kpoints
        type(wave_rotator) :: waverot
        type(wavespin) :: w_dense_ibz
        type(wavedes) :: wdes_dense_ibz
        real(q) :: phonon_freqs_ibz(:,:)
        type(latt) :: latt_cur
        complex(q),intent(in) :: gmode_kkp(:,:,:)
        complex(q),intent(in) :: gcart_kkp(:,:,:,:)
        real(q), allocatable, intent(in) :: velocity(:,:,:,:)
! local variables
        integer :: i, ikibz, ikpibz, iqibz, irotp
        real(q) :: vk(3), vkp(3), vqpt(3)

        
        ikibz = triplet%k%ibz
        ikpibz= triplet%kp%ibz
        iqibz = triplet%q%ibz

        ikibz = triplet%k%ibz
        vk = kpoints%vkpt(:,ikibz)
        call kpoints_fbz_to_ibz(kpoints,triplet%kp%kindex,ikpibz,irotp)
        vkp = rotate_vkpt_irot(kpoints%vkpt(:,ikpibz),waverot%igrpop,irotp)
        vqpt = vk-vkp

        if (elph_set%selfen_fan) then
           do i=1,self%naccumulators_el
              call elph_accumulator_selfen_fan(self%selfen_el(i),kpoints,latt_cur,&
                           ikibz,isp,vkp,vqpt,&
                           w_dense_ibz%celen(:,:,isp),phonon_freqs_ibz,&
                           real(w_dense_ibz%celen(:,ikpibz,isp),q),self%band_idx,phonon_freqs_ibz(:,iqibz),&
                           band_start_k,band_start_kp,gmode_kkp,1.0_q/kpoints%nkfull,wdes_dense_ibz%comm_inb,velocity,waverot%igrpop)
           enddo
        endif

! check if q=0
        if (ikibz==ikpibz.and.irotp==1) then
            if (elph_set%selfen_dw) then
                do i=1,self%naccumulators_el
                   call elph_accumulator_selfen_dw_g(self%selfen_el(i),ikibz,isp,&
                                                     real(w_dense_ibz%celtot(:,ikpibz,isp),q),&
                                                     real(w_dense_ibz%celen(:,ikpibz,isp),q),self%band_idx,band_start_k,band_start_kp,gcart_kkp,wdes_dense_ibz%comm_inb)
                enddo
            endif
        endif
        
    end subroutine elph_accumulators_accumulate

!> Prepare the computation of the Debye-Waller term
    subroutine elph_accumulators_prepare_dw(self,elph_set,fbz2ibz,elph_ifc,wdes_dense_ibz,comm_global)
        use elphon_base, only: elph_settings
        use wave_struct_def, only: wavedes
        use wave_rotate, only: kpoint_star, fbz2ibz_stars, kstar_from_fbz2ibz
        use elphon_potential, only: phonon_ifc, phonon_ifc_get_phonons_zero_disp
        type(elph_accumulators_selfen) :: self
        type(elph_settings) :: elph_set
        type(fbz2ibz_stars) :: fbz2ibz
        type(phonon_ifc) :: elph_ifc
        type(wavedes) :: wdes_dense_ibz
        type(communic) :: comm_global
! local variables
        type(kpoint_star) :: kstar
        integer :: ikp_ibz, ikp_star, i
        real(q) :: vkp(3), vk(3), vqpt(3)
        real(q),allocatable :: w_q(:)
        complex(q),allocatable :: e_q(:,:)
! Compute phonon modes in the FBZ and accumulate for computation of the Debye-Waller term
        if (elph_set%selfen_dw.or.elph_set%selfen_dws) then
           allocate(w_q(1:3*elph_ifc%natoms_pc))
           allocate(e_q(1:3*elph_ifc%natoms_pc,1:3*elph_ifc%natoms_pc))

           do ikp_ibz=1,wdes_dense_ibz%nkpts

              if (mod(ikp_ibz-1,comm_global%ncpu).ne.comm_global%node_me-1) cycle

              call kstar_from_fbz2ibz(kstar, wdes_dense_ibz%vkpt(:,ikp_ibz), &
                                                      fbz2ibz, ikp_ibz)
              do ikp_star=1,kstar%nkpts
                 vk = [0.0_q,0.0_q,0.0_q]
                 vkp = kstar%vkpt(:,ikp_star)
                 vqpt=vk-vkp
                 call phonon_ifc_get_phonons_zero_disp(elph_ifc,vqpt,w_q,e_q)
                 do i=1,elph_set%selfen_naccumulators
                    call elph_accumulator_selfen_dw_t(self%selfen_el(i),&
                                                      w_q,e_q,1.0_q/fbz2ibz%nkfbz)
                 enddo
              end do !loop over k' star
           end do !loop over k'
! sum t matrix
           do i=1,elph_set%selfen_naccumulators
              call elph_accumulator_selfen_dw_sum_t(self%selfen_el(i),comm_global)
           enddo
        endif
    end subroutine elph_accumulators_prepare_dw

!> Finalize calculation of the self-energy
    subroutine elph_accumulators_finalize(self,w_dense_ibz,comm_global)
        use wave_struct_def, only: wavespin
        type(elph_accumulators_selfen) :: self
        type(wavespin),intent(in) :: w_dense_ibz
        type(communic),intent(in) :: comm_global
! local variables
        integer :: i
        do i=1,self%naccumulators_el
            call elph_accumulator_selfen_sum(self%selfen_el(i),comm_global)
            call elph_accumulator_selfen_kk(self%selfen_el(i))
            call elph_accumulator_selfen_degavg(self%selfen_el(i),w_dense_ibz)
        enddo
    end subroutine elph_accumulators_finalize

!> Write results to hdf5 file
    subroutine elph_accumulators_print(self,w_dense_ibz,io)
        use base, only: in_struct
        use wave_struct_def, only: wavespin

        USE vhdf5, ONLY: IH5OUTFILEID

        type(elph_accumulators_selfen) :: self
        type(wavespin),intent(in) :: w_dense_ibz
        type(in_struct),intent(in) :: io
! local variables
        integer :: i

!write metadata about electron-phonon self-energy accumulators

        call elph_accumulators_selfen_hdf5write(self,ih5outfileid)


        do i=1,self%naccumulators_el
            call elph_accumulator_selfen_print_result(self%selfen_el(i),i,w_dense_ibz%celtot,io%iu6)
            call elph_accumulator_selfen_gaps_compute(self%selfen_el(i))
            call elph_accumulator_selfen_gaps_write(self%selfen_el(i),i,io%iu6)
            call elph_accumulator_selfen_gaps_write(self%selfen_el(i),i,io%iu0)

            call elph_accumulator_selfen_hdf5write(self%selfen_el(i),ih5outfileid,i)

        enddo
    end subroutine elph_accumulators_print

!> intialize structure to store the electron-phonon matrix elements locally
    subroutine elph_accumulator_g_init(self,nkpts_k,nkpts_kp,nkpts_k_local,nkpts_kp_local,natoms,nspin,&
                                       band_start_k,band_stop_k,nbands_k,band_start_kp,band_stop_kp,nbands_kp,band_idx,comm)
        type(elph_accumulator_g) :: self
        integer, intent(in) :: nkpts_k
        integer, intent(in) :: nkpts_kp
        integer, intent(in) :: nkpts_k_local
        integer, intent(in) :: nkpts_kp_local
        integer, intent(in) :: natoms
        integer, intent(in) :: nspin
        integer, intent(in) :: band_start_k
        integer, intent(in) :: band_stop_k
        integer, intent(in) :: nbands_k
        integer, intent(in) :: band_start_kp
        integer, intent(in) :: band_stop_kp
        integer, intent(in) :: nbands_kp
        integer, intent(in) :: band_idx(:)
        type(communic),intent(in) :: comm

        self%nkpts_k = nkpts_k
        self%nkpts_kp = nkpts_kp
        self%nkpts_k_local = nkpts_k_local
        self%natoms = natoms
        self%nspin = nspin
        self%band_start_k = band_start_k
        self%band_stop_k = band_stop_k
        self%nbands_k = nbands_k
        self%band_start_kp = band_start_kp
        self%band_stop_kp = band_stop_kp
        self%nbands_kp = nbands_kp
        self%band_idx = band_idx
!consistency check
        self%nbands_kp_local = size(self%band_idx)
! communicator
        self%comm = comm
! reset k-point counter
        self%ikkp_counter = 0
        self%ikkp_max = nkpts_k_local*nkpts_kp_local
! big allocations
        allocate(self%g(self%nbands_k,self%nbands_kp_local,self%natoms*3,self%ikkp_max,self%nspin))
        allocate(self%w_q(self%natoms*3,self%ikkp_max))
        allocate(self%idx_k(self%ikkp_max))
        allocate(self%idx_kp(self%ikkp_max))
    end subroutine elph_accumulator_g_init

!> accumulate g matrix elements for a certain k and k'
    subroutine elph_accumulator_g_accumulate(self,ik,ikp,ispin,g,w_q)
        use tutor, only: vtutor
        use string, only: str
        type(elph_accumulator_g) :: self
        integer, intent(in) :: ik !< index of the k point
        integer, intent(in) :: ikp !< index of the k' point
        integer, intent(in) :: ispin
        complex(q),intent(in) :: g(:,:,:)
        real(q),intent(in) :: w_q(:)
! consistency check
        if (size(g,1) /= self%nbands_k) then
            call vtutor%bug("Inconsistent dimension 1 of the electron-phonon matrix elements expected: "//&
                            str(self%nbands_k)//" got: "//str(size(g,1)),"elphon_accumulators.F",2451)
        endif
        if (size(g,2) /= self%nbands_kp_local) then
            call vtutor%bug("Inconsistent dimension 2 of the electron-phonon matrix elements expected: "//&
                            str(self%nbands_kp_local)//" got: "//str(size(g,2)),"elphon_accumulators.F",2455)
        endif
        if (size(g,3) /= self%natoms*3) then
            call vtutor%bug("Inconsistent dimension 3 of the electron-phonon matrix elements expected: "//&
                            str(self%natoms*3)// " got: "//str(size(g,3)),"elphon_accumulators.F",2459)
        endif
        self%ikkp_counter = self%ikkp_counter + 1
        self%g(:,:,:,self%ikkp_counter,ispin) = g(:,:,:)
        self%w_q(:,self%ikkp_counter) = w_q
        self%idx_k(self%ikkp_counter) = ik
        self%idx_kp(self%ikkp_counter) = ikp
    end subroutine elph_accumulator_g_accumulate


!> write matrix elements to an hdf5 file and free memory
    subroutine elph_accumulator_g_vh5write_free(self,filename)
        use vhdf5
        use tutor, only: vtutor
        type(elph_accumulator_g) :: self
        character(len=*),optional :: filename
! local variables
        character(:),allocatable :: my_filename
        integer(HID_T) :: fileid, mels_groupid, phonons_groupid
! local variables
        integer :: irank, irank_remote
        integer(size_t) :: dims(6)

! check if the g array is complete
        if (self%ikkp_counter/=self%ikkp_max) then
           call vtutor%bug('elph_accumulator_g_vh5write does not contain all the k-points yet',"elphon_accumulators.F",2484)
        endif

        irank=self%comm%node_me
        if (irank==1) then
! create an hdf5 file to write the matrix elements to
           my_filename = "vaspelph.h5"; if (present(filename)) my_filename = filename
           call vh5_error(vh5_file_create_or_overwrite(my_filename,fileid),"elphon_accumulators.F",2491)

! open matrix elements group
           call vh5_error(vh5_group_open_or_create(fileid, 'matrix_elements', mels_groupid),"elphon_accumulators.F",2494)
! open the phonon group
!call vh5_error(vh5_group_open_or_create(fileid, GRP_PHONONS, phonons_groupid),"elphon_accumulators.F",2496)
! for now we keep compatibility
           phonons_groupid = mels_groupid

! prepare the electron-phonon array for writting in parallel
! (this is a complex array but I can create it as real)
           dims = [self%nbands_k,self%nbands_kp,3*self%natoms,self%nkpts_k,self%nkpts_kp,self%nspin]
           call vh5_error(vh5_create_double_complex_array_nd(mels_groupid, 'elph', 6, dims),"elphon_accumulators.F",2503)

! prepare the phonon frequencies array
           dims(:3) = [3*self%natoms,self%nkpts_k,self%nkpts_kp]
           call vh5_error(vh5_create_double_array_nd(phonons_groupid, 'phonon_eigenvalues', 3, dims),"elphon_accumulators.F",2507)

! write some dimensions to file
           call vh5_error(vh5_write(mels_groupid, 'band_start_k', self%band_start_k),"elphon_accumulators.F",2510)
           call vh5_error(vh5_write(mels_groupid, 'nbands_k', self%nbands_k),"elphon_accumulators.F",2511)
           call vh5_error(vh5_write(mels_groupid, 'band_start_kp', self%band_start_kp),"elphon_accumulators.F",2512)
           call vh5_error(vh5_write(mels_groupid, 'nbands_kp', self%nbands_kp),"elphon_accumulators.F",2513)
           call vh5_error(vh5_write(mels_groupid, 'natoms', self%natoms),"elphon_accumulators.F",2514)
           call vh5_error(vh5_write(mels_groupid, 'nspin', self%nspin),"elphon_accumulators.F",2515)
           call vh5_error(vh5_write(mels_groupid, 'nkpts_k', self%nkpts_k),"elphon_accumulators.F",2516)
           call vh5_error(vh5_write(mels_groupid, 'nkpts_kp', self%nkpts_kp),"elphon_accumulators.F",2517)

! receive data from rank 2 up to NCPU
           do irank_remote=1,self%comm%ncpu
              if (irank_remote>1) then
! receive number of k-points
                 CALL M_recv_i(self%comm, irank_remote, self%ikkp_counter, 1)
! allocate and receive k-point lists
                 allocate(self%idx_k(self%ikkp_counter))
                 allocate(self%idx_kp(self%ikkp_counter))
                 CALL M_recv_i(self%comm, irank_remote, self%idx_k, self%ikkp_counter)
                 CALL M_recv_i(self%comm, irank_remote, self%idx_kp, self%ikkp_counter)
! allocate and receive band indexes
                 allocate(self%band_idx(self%nbands_kp_local))
                 CALL M_recv_i(self%comm, irank_remote, self%band_idx, self%nbands_kp_local)
! allocate and receive g matrix elements
                 allocate(self%g(self%nbands_k,self%nbands_kp_local,self%natoms*3,self%ikkp_counter,self%nspin))
                 CALL M_recv_z(self%comm, irank_remote, self%g, size(self%g))
! allocate and receive phonon frequencies
                 allocate(self%w_q(self%natoms*3,self%ikkp_counter))
                 CALL M_recv_d(self%comm, irank_remote, self%w_q, size(self%w_q))
              endif
! write local information and deallocate
              call elph_accumulator_g_vh5write_block(self,phonons_groupid,mels_groupid,self%ikkp_counter,self%idx_k,self%idx_kp,self%band_idx,self%g,self%w_q)
              deallocate(self%idx_k)
              deallocate(self%idx_kp)
              deallocate(self%band_idx)
              deallocate(self%g)
              deallocate(self%w_q)
           enddo
! close file and quit
           if (lwriteh5) then
!call vh5_error(vh5_group_close(phonons_groupid),"elphon_accumulators.F",2549)
              call vh5_error(vh5_group_close(mels_groupid),"elphon_accumulators.F",2550)
              call vh5_error(vh5_file_close(fileid),"elphon_accumulators.F",2551)
           endif
        else
! send number of k-points
           CALL M_send_i(self%comm, 1, self%ikkp_counter, 1)
! send k-point lists
           CALL M_send_i(self%comm, 1, self%idx_k, self%ikkp_counter)
           CALL M_send_i(self%comm, 1, self%idx_kp, self%ikkp_counter)
! send band indexes
           CALL M_send_i(self%comm, 1, self%band_idx, self%nbands_kp_local)
! send g matrix elements
           CALL M_send_z(self%comm, 1, self%g, size(self%g))
! send phonon frequencies
           CALL M_send_d(self%comm, 1, self%w_q, size(self%w_q))
           deallocate(self%idx_k)
           deallocate(self%idx_kp)
           deallocate(self%band_idx)
           deallocate(self%g)
           deallocate(self%w_q)
        endif
    end subroutine elph_accumulator_g_vh5write_free

!> Write a block of g matrix elements
    subroutine elph_accumulator_g_vh5write_block(self,phonons_groupid,mels_groupid,ikkp_counter,idx_k,idx_kp,band_idx,g,w_q)
        use vhdf5
        use tutor, only: vtutor
        type(elph_accumulator_g) :: self
        integer,intent(in) :: ikkp_counter
        integer,intent(in) :: idx_k(:)
        integer,intent(in) :: idx_kp(:)
        integer,intent(in) :: band_idx(:)
        complex(q),intent(in) :: g(:,:,:,:,:)
        real(q),intent(in) :: w_q(:,:)
        integer(HID_T),intent(in) :: phonons_groupid
        integer(HID_T),intent(in) :: mels_groupid
! local variables
        integer :: isp, imode, ikp, ik, ikkp, ib_global
        integer :: ib_local
        integer(SIZE_T) :: start(6)
        integer(SIZE_T) :: count(6)

        if (ikkp_counter /= size(idx_k)) then
            call vtutor%bug('Incorrect dimensions of idx_k',"elphon_accumulators.F",2593)
        endif

        do ikkp=1,ikkp_counter
            ik = idx_k(ikkp)
            ikp = idx_kp(ikkp)
! write phonon frquencies
            start(:3) = [1,ik,ikp]
            count(:3) = [self%natoms*3,1,1]
            call vh5_error(vh5_write_double_subarray_nd(phonons_groupid, "phonon_eigenvalues", 3, start, count, w_q(:,ikkp)),"elphon_accumulators.F",2602)
! write g matrix elements
            do isp=1,self%nspin
               do imode=1,self%natoms*3
                  do ib_local=1,self%nbands_kp_local
                      ib_global = band_idx(ib_local)
                      start = [1,ib_global,imode,ik,ikp,isp]
                      count = [self%nbands_k,1,1,1,1,1]
                      call vh5_error(vh5_write_double_complex_subarray_nd(mels_groupid, "elph", 6, start, count, g(:,ib_local,imode,ikkp,isp)),"elphon_accumulators.F",2610)
                  enddo
               enddo
            enddo
        enddo
    end subroutine elph_accumulator_g_vh5write_block


!> @brief create a the file containing the electron-phonon matrix elements
!> and intialize with the basic structure.
!> This file is open in parallel so that each CPU in the batch can write its own data
!>
!> The basic structure of the file will be:
!> - vkpt_ibz : contains the reduced coordinates of the k-points for which the
!>              electron-phonon matrix elements were computed
!> - indx_fbz2ibz : contains the mapping between the FBZ and the IBZ
!> - irot_ibz2fbz : contains an index with the operation that rotates the IBZ point to the FBZ,
!>                  negative index means that that inversion symmetry is applied as well
!> - elph : contains the the electron-phonon matrix elements (n1,n2,ph_mode,nk1,nk2)
!> - eigenvalues : contains the ibz eigenvalues
    subroutine elph_accumulator_g_hdf5mpiio_init(self,nb_tot,band_start_k,nbands_k,band_start_kp,nbands_kp,natoms,nkpts_k,nkpts_kp,nspin,comm)
        use mpimy
        use vhdf5
        type(elph_accumulator_g_hdf5mpiio) :: self
        integer :: nb_tot !< Total number of bands computed in the electron-phonon driver
        integer :: band_start_k !< Number of bands for k
        integer :: nbands_k !< Number of bands for k
        integer :: band_start_kp !< Number of bands for k
        integer :: nbands_kp !< Number of bands for k'
        integer :: nkpts_k !< Number of points for k
        integer :: nkpts_kp !< Number of points for k'
        integer :: nspin !< Number of spin channels
        integer :: natoms !< Number of atoms
        type(communic) :: comm

! local variables
        integer(size_t) :: dims(6)

! create a file with MPIIO access
        call vh5_error(vh5_file_create_or_overwrite('vaspelph.h5',self%fileid,comm=comm%mpi_comm),"elphon_accumulators.F",2649)
! From now on all the CPUS write everything.
! This must be undone once we are finished and before calling any other HDF5 routines
! TODO: this should be 1._q in another way
        self%lwritesave = lwriteh5
        self%natoms = natoms
        self%nkpts_k = nkpts_k
        self%nkpts_kp = nkpts_kp
        self%nbands_k = nbands_k
        self%nbands_kp = nbands_kp
        lwriteh5 = .true.
        self%plist_id = vh5_data_access_mpiio()

! open group
        call vh5_error(vh5_group_open_or_create(self%fileid, 'matrix_elements', self%mels_groupid),"elphon_accumulators.F",2663)
!call vh5_error(vh5_group_open_or_create(self%fileid, GRP_PHONONS, self%phonons_groupid),"elphon_accumulators.F",2664)
! for now we keep compatibility
        self%phonons_groupid = self%mels_groupid

! write some dimensions to file
        call vh5_error(vh5_write(self%mels_groupid, 'nbtot', nb_tot),"elphon_accumulators.F",2669)
        call vh5_error(vh5_write(self%mels_groupid, 'band_start_k', band_start_k),"elphon_accumulators.F",2670)
        call vh5_error(vh5_write(self%mels_groupid, 'nbands_k', nbands_k),"elphon_accumulators.F",2671)
        call vh5_error(vh5_write(self%mels_groupid, 'band_start_kp', band_start_kp),"elphon_accumulators.F",2672)
        call vh5_error(vh5_write(self%mels_groupid, 'nbands_kp', nbands_kp),"elphon_accumulators.F",2673)
        call vh5_error(vh5_write(self%mels_groupid, 'natoms', self%natoms),"elphon_accumulators.F",2674)
        call vh5_error(vh5_write(self%mels_groupid, 'nspin', nspin),"elphon_accumulators.F",2675)
        call vh5_error(vh5_write(self%mels_groupid, 'nkpts_k', self%nkpts_k),"elphon_accumulators.F",2676)
        call vh5_error(vh5_write(self%mels_groupid, 'nkpts_kp', nkpts_kp),"elphon_accumulators.F",2677)

! prepare the electron-phonon array for writting in parallel
        dims = [nbands_k,nbands_kp,3*self%natoms,nkpts_k,nkpts_kp,nspin]
        call vh5_error(vh5_create_double_complex_array_nd(self%mels_groupid, 'elph', 6, dims),"elphon_accumulators.F",2681)
        dims(:3) = [3*self%natoms,self%nkpts_k,nkpts_kp]
        call vh5_error(vh5_create_double_array_nd(self%phonons_groupid, 'phonon_eigenvalues', 3, dims),"elphon_accumulators.F",2683)

    end subroutine elph_accumulator_g_hdf5mpiio_init

!> @brief write a block of electron-phonon matrix elements
    subroutine elph_accumulator_g_hdf5mpiio_accumulate(self,ik,wdesk2,ikp,isp,elph,w_q)
        use vhdf5
        use tutor, only: vtutor
        use string, only: str
        use wave_struct_def, only: wavedes1
        use wave, only: nb_local
        type(elph_accumulator_g_hdf5mpiio),intent(in) :: self
        integer,intent(in) :: ik !< index of the k point
        type(wavedes1),intent(in) :: wdesk2 !< wavefunction descriptor of k' (can be distributed)
        integer,intent(in) :: ikp !< index of the k' point
        integer,intent(in) :: isp !< index of the spin channel to be written
        complex(q),intent(in) :: elph(:,:,:) !< Matrix block containing the electron-phonon matrix elements
        real(q),intent(in) :: w_q(:)

! Local variables
        integer :: ib, ib_global, imode
        integer(size_t) :: start(6), count(6)
        complex(q),allocatable :: elph_tmp(:,:,:)

! Some safety checks!
        if (ik > self%nkpts_k) then
            call vtutor%bug("Invalid k-point index: "//str(ik)//&
                            " should be <= "//str(self%nkpts_k),"elphon_accumulators.F",2710)
        endif
        if (ikp > self%nkpts_kp) then
            call vtutor%bug("Invalid k'-point index: "//str(ikp)//&
                            " should be <= "//str(self%nkpts_kp),"elphon_accumulators.F",2714)
        endif

        if (size(elph,1) /= self%nbands_k) then
            call vtutor%bug('Invalid first dimension of ELPH array found: '//str(size(elph,1))//&
                            ' should be: '//str(self%nbands_k),"elphon_accumulators.F",2719)
        endif
        if (size(elph,2) /= wdesk2%nbands) then
            call vtutor%bug('Invalid second dimension of ELPH array found: '//str(size(elph,2))//&
                            ' should be: '//str(wdesk2%nbands),"elphon_accumulators.F",2723)
        endif
        if (self%nbands_kp /= wdesk2%nb_tot) then
            call vtutor%bug('Invalid value for nbands_kp found: '//str(self%nbands_kp)//&
                            ' should be: '//str(wdesk2%nb_tot),"elphon_accumulators.F",2727)
        endif

        
! Only (1._q,0._q) of the nodes in the intra-band group needs to write because the results were gathered
        if (wdesk2%comm_inb%node_me==1) then
! Only (1._q,0._q) of the nodes sharing the band needs to write because the results were gathered
        if (wdesk2%comm_inter%ncpu==1) then
! All the bands are present locally so we can write in (1._q,0._q) big block
            start = [1,1,1,ik,ikp,isp]
            count = [self%nbands_k,self%nbands_kp,3*self%natoms,1,1,1]
            call vh5_error(vh5_write_double_complex_subarray_nd(self%mels_groupid, "elph", 6, start, count, elph, xfer_prp = self%plist_id ),"elphon_accumulators.F",2738)
        else
# 2755

! gather all information onto (1._q,0._q) node and write
            allocate(elph_tmp(self%nbands_k,self%nbands_kp,3*self%natoms))
            elph_tmp = 0
            do imode=1,3*self%natoms
                do ib_global=1,wdesk2%nb_tot
                    ib = nb_local(ib_global,wdesk2)
                    if (ib==0) cycle
                    elph_tmp(:,ib_global,imode) = elph(:,ib,imode)
                enddo
            enddo
            CALL M_sum_z( wdesk2%comm_inter, elph_tmp, size(elph_tmp) )
            start = [1,1,1,ik,ikp,isp]
            count = [self%nbands_k,self%nbands_kp,3*self%natoms,1,1,1]
            call vh5_error(vh5_write_double_complex_subarray_nd(self%mels_groupid, "elph", 6, start, count, elph_tmp, xfer_prp = self%plist_id ),"elphon_accumulators.F",2769)
            deallocate(elph_tmp)

        endif
        endif

        start(:3) = [1,ik,ikp]
        count(:3) = [3*self%natoms,1,1]
        call vh5_error(vh5_write_double_subarray_nd(self%phonons_groupid, "phonon_eigenvalues", 3, start, count, w_q, xfer_prp = self%plist_id ),"elphon_accumulators.F",2777)
        

    end subroutine elph_accumulator_g_hdf5mpiio_accumulate

!> @brief Close the file containing the electron-phonon matrix elements
    subroutine elph_accumulator_g_hdf5mpiio_free(self)
        use vhdf5
        type(elph_accumulator_g_hdf5mpiio) :: self


! close the matrix elements group
        call vh5_error(vh5_group_close_writing(self%mels_groupid),"elphon_accumulators.F",2789)
!call vh5_error(vh5_group_close_writing(self%phonons_groupid),"elphon_accumulators.F",2790)

! close the file
        call vh5_error(vh5_file_close(self%fileid),"elphon_accumulators.F",2793)

! restore status where only master writes
        lwriteh5 = self%lwritesave

! close plist
        call vh5_error(vh5_plist_close(self%plist_id),"elphon_accumulators.F",2799)

    end subroutine elph_accumulator_g_hdf5mpiio_free

!> Initialize accumulators for matrix elements
    subroutine elph_accumulators_g_init(self,band_start_k,band_stop_k,nbands_k,&
                                        band_start_kp,band_stop_kp,nbands_kp,nkpts_k,nkpts_k_local,nkpts_kp,nkpts_kp_local,ispin,natoms_pc,band_idx,lhdf5_mpiio,comm_global)
        type(elph_accumulators_g) :: self
        integer,intent(in) :: band_start_k
        integer,intent(in) :: band_stop_k
        integer,intent(in) :: nbands_k
        integer,intent(in) :: band_start_kp
        integer,intent(in) :: band_stop_kp
        integer,intent(in) :: nbands_kp
        integer,intent(in) :: nkpts_k
        integer,intent(in) :: nkpts_k_local
        integer,intent(in) :: nkpts_kp
        integer,intent(in) :: nkpts_kp_local
        integer,intent(in) :: ispin
        integer,intent(in) :: natoms_pc
        integer,intent(in) :: band_idx(:)
        logical,intent(in) :: lhdf5_mpiio
        type(communic) :: comm_global
        self%lhdf5_mpiio = lhdf5_mpiio
        if (lhdf5_mpiio) then
! open an HDF5 file in parallel
            call elph_accumulator_g_hdf5mpiio_init(self%g_hdf5mpiio,nbands_k,band_start_k,nbands_k,band_start_kp,nbands_kp,&
                                      natoms_pc,nkpts_k,&
                                      nkpts_kp,ispin,comm_global)
        else
! store electron-phonon matrix elements in memory and then write to file
            call elph_accumulator_g_init(self%g_store,nkpts_k,nkpts_kp,nkpts_k_local,nkpts_kp_local,natoms_pc,ispin,&
                                band_start_k,band_stop_k,nbands_k,band_start_kp,band_stop_kp,nbands_kp,band_idx,comm_global)
        endif
    end subroutine elph_accumulators_g_init

!> Accumulate matrix elements in the g matrix elements accumlators
    subroutine elph_accumulators_g_accumulate(self,wdes1,ik,ikp_ibz,isp,gmode_kkp,w_q)
        use wave_struct_def, only: wavedes1
        type(elph_accumulators_g) :: self
        type(wavedes1) :: wdes1
        integer,intent(in) :: ik
        integer,intent(in) :: ikp_ibz
        integer,intent(in) :: isp
        complex(q),intent(in) :: gmode_kkp(:,:,:)
        real(q),intent(in) :: w_q(:)
        if (self%lhdf5_mpiio) then
! write the matrix elements to hdf5 file
            call elph_accumulator_g_hdf5mpiio_accumulate(self%g_hdf5mpiio, ik, &
                          wdes1, ikp_ibz, isp, gmode_kkp, w_q)
        else
            call elph_accumulator_g_accumulate(self%g_store,ik,ikp_ibz,isp,gmode_kkp,w_q)
        endif
    end subroutine elph_accumulators_g_accumulate

!> Free electron-phonon g matirx elements accumulators
    subroutine elph_accumulators_g_free(self)
        type(elph_accumulators_g) :: self

        if (self%lhdf5_mpiio) then
            call elph_accumulator_g_hdf5mpiio_free(self%g_hdf5mpiio)
        else
            call elph_accumulator_g_vh5write_free(self%g_store)
        endif

    end subroutine elph_accumulators_g_free

end module
