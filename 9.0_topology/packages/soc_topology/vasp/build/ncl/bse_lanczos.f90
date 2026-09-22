# 1 "bse_lanczos.F"
!pm:  implementation of the Lanczos algorithm for computing the absorption spectrum
! Implementation is based on the reformulation of the problem of computing the dielectric function, \epsilon_{ij}(\omega) as a continued fraction
! For the diagonal elements, {ii},  we have
!                                                    |u_0^{ii}|^2
! \epsilon_{ii}(\omega) = 1 - ----------------------------------------------------------        (1)
!                                                                    b_1
!                             \omega  - a_1 + i\eta + ----------------------------------
!                                                                                b_2
!                                                     \omega - a_2 + i\eta + -----------
!                                                                                ...
!
! where:
! i)   |u_0^{ii}|^2 is the amplitude of the initial guess vector, before being normalised to 1. We take the dipole moments as the initial guess,
! from OPTMAT, so {ii} is the direction of the dipole.
! ii)  a_i and b_i are the Lanczos coefficients, and are real numbers.
! iii) \omega is the frequency, and \eta is the width, controlled by CSHIFT inside the INCAR.
!
! For the off-diagonal elements, this has to be modified. The assumption at the start is that \epsilon has an operational form like
!                                 1
! \epsilon(\omega) = 1 - <P|------------|P> (2),
!                            \omega - H
!
! where H is the BSE hamiltonian, and <P|P> = |u_0|^2 is real. But if the left and right vectors are not the same, i.e. they are dipole moments
! along different directions, we cannot follow this and there is no guarantee that <P|P> is real. We solve this by rotating the matrix 45 degrees and using
!                                         1
! \bar\epsilon_{m}(\omega) = - <R_m|------------|R_m> (3),
!                                    H - \omega
! where |R_m> is the rotated dipole vector. We will have
! |R_1> = 1/\sqrt{2}(|P_x> + |P_y>)   (4a)
! |R_2> = 1/\sqrt{2}(|P_y> + |P_z>)   (4b)
! |R_3> = 1/\sqrt{2}(|P_z> + |P_x>)   (4c)
!
! This way the vectors on the left and right side of (3) are the same, and <R_m|R_m> is a real number. At the end the three \bar\epsilon auxiliar functions ! are rotated back so that we can recover the {xy}, {xz}, and {yz} components
!
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


# 36 "bse_lanczos.F" 2 
# 1 "./bse.inc" 1 
! currently the default is to use double precision (flag double_prec_bse);
! the flag single_prec_bse allows to select single precision storage of the BSE
! matrix and thus saves a factor 2 in storage (but usually also in CPU time!!);
! important: use of single precision is absolutely safe and highly recommended!!
! if you want to use single precision define precompiler flag single_prec_bse !

# 9


! the flag USE_ZHEEVX determines whether ZHEEVX is used; (1._q,0._q) has to know that
! ZHEEVX is about twice to three times faster than ZHEEV(D), but also requires
! twice as much memory -- that's the downside of the coin and maybe a problem
!
! alternatively, (1._q,0._q) can also use ZHEEVD by defining precompiler flag USE_ZHEEVD
! which uses "divide-and-conquer" algorithms for diagonalization and could be
! a bit faster than ZHEEV for larger matrices (if eigenvectors are computed);
! currently USE_ZHEEVX is the default if no different precompiler flag is given
!
! for those which like to use/test ZHEEV another (undocument) flag BSE_ZHEEV
! may be set to cancel the definition of USE_ZHEEVX (and/or USE_ZHEEVD ...);
! this only affects the BSE routines (and no use of USE_ZHEEVX/D elsewhere ...);

# 43


# 64


# 85





























# 116




# 122




# 144



! Debugging
!#define bsedebug

# 159



# 164





# 37 "bse_lanczos.F" 2 
module bse_lanczos

  use prec     
  use constant 
  use wave_high     
  use bse_te
  use bse_struct
  use vaspxml
  use iso_c_binding

  implicit none

!> (1) bse, (0) ip -> according to settings of ladder and lhartree (and ltriplet)
      integer,private  :: exciton
!> (1) going beyond ("antires==2"), (0) sticking with tamm-dancoff approximation
      integer,private  :: beyondtd
!> number of spin components (2 for collinear spin-polarization and 1 else); just needed for a norm factor (containing a spin-degeneracy factor ...)
      integer,private  :: ispin
!> number of k-points -> what is nkpts in vasp ...
      integer,private  :: maxk
!  some further (local) variables

!> matrix dimension used for the calculation (-> rank of bse matrix from bse::calculate_bse)
      integer,private  :: matdim
!> number of rows at each node (bse matrix distributed for 1 jobs ... !)
      integer,private, allocatable  :: matdim_node(:)
!> corresponding "global index" of first local matrix element at each node
      integer,private, allocatable  :: matstart(:)

!> the imaginary unit and complex (0._q,0._q) are always needed
      complex(mq),private, parameter  :: imun     = (0.0_mq,1.0_mq)
      complex(mq),private, parameter  :: czero    = (0.0_mq,0.0_mq)

      integer,private :: node_lead=0    !< number of the leading process.
      integer,private :: node_me        !< number of my process (set by vasp).
      integer,private :: nranks         !< total number of processes (set by vasp).
      type(communic),private :: my_comm !< communicator to be used (set by vasp).



contains

!****************** subroutine prepare_lanczos_paralell    ******************************
!
!> sets up the parallelization structure to be used with the lanczos algorithm
!
!*************************************************************************************

  subroutine prepare_lanczos_parallel(whf, io, ncv_)

!! use moffload_struct_def

    use base
    use ini
    use string, only: str
    use tutor, only: vtutor
    type (wavespin)      whf
    type (in_struct)     io
    integer :: ncv_,ncv         !< total number of pair states
    logical :: lbeyondtd        !< beyond tamm-dancoff or not
    integer :: nkpts,nbtd
    integer :: hstrip
    integer :: n

    integer :: i
    integer  :: istat

    integer(qi8) :: size8,size8a,size8b
    integer :: nprow, npcol, myrow, mycol, np, nq

! argument NCV_ will be set negative if BSEMATRIX shall not be allocated (for calculation of IP spectra we do not need it, just the distribution and LDAGW)
    ncv=abs(ncv_)
! the strip size is ceil(ncv/nranks), some nodes have less or no data
    node_me = whf%wdes%comm%node_me-1
    nranks = whf%wdes%comm%ncpu

! column of ranks
    call BLACS_GRIDINFO(bse_desc(2), nprow, npcol, myrow, mycol)
    np = numroc(ncv,bse_desc(5),myrow,0,nprow)
    nq = numroc(ncv,bse_desc(6),mycol,0,npcol)

! determine the local number of rows
    hstrip = max(1,numroc(ncv,bse_desc(5),0,0,nprow))
    ncv_local = np

    i=ncv_local

! simple check (should never fail -- unless 1 communication has a problem ...)
    CALL m_sum_i(whf%wdes%comm, i, 1 )
    if (i /= ncv) then
       call vtutor%bug("internal error in prepare_lanczos_parallel: dimension wrong " // str(I) // " " // str(NCV) // " " // str(NCV_LOCAL), &
           "bse_lanczos.F", 128)
    endif

! allocate arrays matdim_node and matstart here and get some mpi parameters
    allocate(matdim_node(0:nranks-1))
    allocate(matstart(0:nranks))

! now store the local dimension ncv_local in matdim_node(node_me) ...
    matdim_node=0
    do i=0,nranks-1
       if (i==node_me) matdim_node(i)=numroc(ncv,bse_desc(5),i,0,nprow)
    enddo
! ... and synchronize all nodes (here by summing instead of broadcasting ;-))
    CALL m_sum_i(whf%wdes%comm, matdim_node, whf%wdes%nb_par )
! now we calculate the corresponding global index (minus (1._q,0._q)) of the first element stored on each node and store it in MATSTART(MY_NODE) [it's basically
! an "offset" with respect to the first element in the   global   array ...]
    matstart(0)=0

    nq=0
    do i=1,nranks
       nq=nq+numroc(ncv,bse_desc(5),i-1,0,nprow)
       matstart(i)=nq
    enddo

! calculate the base (first index) on each node and check against MATSTART
    first_row_index=1
    do i=1,nranks
       if (i==node_me+1) exit
! add number of data on node i
       first_row_index=min(ncv,node_me*hstrip+1)
    enddo
    if (first_row_index/=min(ncv,node_me*hstrip+1)) then
       call vtutor%bug("internal error in prepare_lanczos_parallel: first_row_index not consistent " // "with matstart", "bse_lanczos.F", 160)
    endif

! copy the 1 communicator for later use
    my_comm=whf%wdes%comm

! now allocate the rest ...
    nbtd=1

    size8a=ncv_local
    size8b=ncv * nbtd
    size8=size8a*size8b

    allocate(optmat(ncv, 3),ldagw(ncv))
! if input argument ncv_ was negative do not fully allocate bsematrix ...
    if (ncv_ > 0) then
       if(ncv_local > 0) then
          allocate(bsematrix(ncv_local, ncv*nbtd), stat = istat)
       else
          allocate(bsematrix(1,1), stat = istat)
       endif
# 190

       if ( istat/=0 ) then
          call vtutor%error( "prepare_lanczos_parallel (bsematrix) is not able to allocate "//&
            str(8._q*SIZE8*2_qi8/1024_qi8) //" kB of data on MPI rank 0." )
       endif
       if (io%iu0>=0) write(io%iu0,'(a,f8.3,a,i7)') ' bse lanczos double prec attempting allocation of',8.e-9_q*size8*2_qi8,' gbyte  rank=',ncv
       if (io%iu6>=0) write(io%iu0,'(a,f8.3,a,i7)') ' bse lanczos double prec attempting allocation of',8.e-9_q*size8*2_qi8,' gbyte  rank=',ncv
       call register_allocate(8._q*2_qi8* size8 , "bse")

    else  ! ncv_ < 0
! if not needed do just a dummy allocation of bsematrix
       allocate(bsematrix(1,1))
    endif

!! 

! for a simple check whether the BSE matrix is set properly (and in particular
! completely) later preset it with some non-(0._q,0._q) value (defined in the head)

!$acc enter data create(bsematrix) 
!$acc kernels present(bsematrix) 
    bsematrix=0
!$acc end kernels

!! 

  end subroutine prepare_lanczos_parallel
 

!****************** subroutine calculate_bse_lanczos    ******************************
!
!> main driver to setup the haydock-lanczos solver
!> it should check wheter we are using the td approximation or not and send the
!> necessary parameters
!
!*************************************************************************************

  subroutine calculate_bse_lanczos(shift,omega,ispin_in,maxk_in,ncv,kptweight_in, encut_sp_in, lexciton, lbeyondtd, nedos_in, conv, io, nbseeig)

    use base
    use constant
    use reader_tags
    implicit none 

    integer :: ispin_in,maxk_in,ncv
    real(q) :: shift
    real(q) :: omega
    real(q) :: kptweight_in(maxk_in)
    real(q) :: encut_sp_in
    real(q) :: conv
    real(q) :: bse_prec
    logical :: lexciton
    logical :: lbeyondtd, lopen
    integer :: nedos_in, madtdin, nbseeig, ierr
    character(40):: prec_str
    type (in_struct) :: io

    integer :: it

!> precision control section
!> set default threshold for epsilon
    bse_prec = 1.d-3  !> this is the equivalent of BSEPREC=(M)EDIUM/(m)edium
!> check precision set in INCAR
    call open_incar_if_found(io%iu5, lopen)
    call process_incar(lopen, io%iu0, io%iu5, 'BSEPREC', prec_str, 40, ierr, lwritexml=.false., lcontinue=.true.)

!> check if INCAR existed and if BSEPREC was read
    if ((ierr/=0).and.(ierr/=3)) then
      call vtutor%warning("Error while reading BSEPREC tag from INCAR. Default value ('Normal') will be used!")
      prec_str = 'Normal'
    endif
    if ((prec_str(1:1)=='L') .or. (prec_str(1:1)=='l') .or. (prec_str(1:1)=='F') .or. (prec_str(1:1)=='f')) then
      bse_prec = 1.d-2
    elseif ((prec_str(1:1)=='H') .or. (prec_str(1:1)=='h')) then
      bse_prec = 1.d-4
    elseif ((prec_str(1:1)=='A') .or. (prec_str(1:1)=='a')) then
      bse_prec = 1.d-5
    elseif ((prec_str(1:1)=='M') .or. (prec_str(1:1)=='m') .or. (prec_str(1:1)=='N') .or. (prec_str(1:1)=='n')) then
      bse_prec = 1.d-3
    else 
      call vtutor%warning("Undefined value for tag BSEPREC in INCAR. Default value ('Normal') will be used!")
    endif

    call close_incar_if_found(io%iu5)

!> go beyond tamm-dancoff approximation (antires==2) or not (antires/=2)?
    if (lbeyondtd) then
      if (io%iu0>=0) write(io%iu0,*) 'ANTIRES was set to 2 at some point during execution.'
      lbeyondtd = .false.
      beyondtd = 0
    else
      if (io%iu0>=0) write(io%iu0,*) 'Calling Haydock-Lanczos solver to compute optical spectra.'
    endif

    matdim=ncv
    maxk=maxk_in
!> calculate excitonic spectra or independent-particle spectra only?
    if (lexciton) then
      exciton=1
    else
      exciton=0
    endif

    ispin = ispin_in

    
    if (io%iu0>=0) write(io%iu0,*) 'starting the solver.' 
    if (exciton>=0 .and. beyondtd==0) call lanczos_solver(matdim,omega,encut_sp_in,nedos_in,shift,conv,io,nbseeig,ispin,bse_prec)

    

    if (io%iu0>=0) write(io%iu0,*) 'end of calculate_bse_lanczos'


  end subroutine calculate_bse_lanczos

!****************** subroutine lanczos_bse_solver    ******************************
!
!> computes the haydock-lanczos coefficients and dielectric function
!
!******************************************************************************

  subroutine lanczos_solver(fulldim,omega,omegamax,nedos_in,shift,conv_thr,io,nbseeig,ispin,prec_eps)

    use base
    use string, only: str
    use chi_glb, only: lgwlf

    use vhdf5

    implicit none
    type (in_struct)        :: io
    integer                 :: fulldim

!> like in the time-evolution code, every step computes a new vector from two previous iterations. these are
!> stored in vector_0, vector_1, and v_tmp_1
# 331

    complex(mq), allocatable :: vector_0(:,:)   !< initial and (n-1)th iteration of the Lanczos vectors
    complex(mq), allocatable :: vector_1(:,:)   !< second and n-th iteration of the Lanczos vectors
    complex(mq), allocatable :: v_tmp_1(:,:), v_tmp_2(:,:)   !< auxiliar vector to store vector_1 before copying it to vector_0
    complex(mq) :: ortho_test(1:6,1:6) !< used for all dot procducts between vectors

    real(q), allocatable :: a(:,:,:), b(:,:,:)
    complex(mq), pointer :: spectra1(:,:), spectra2(:,:)
    
    integer  :: id1, id2, it, iw, maxiter, itercount, ierr, nbseeig, ispin, nedos_in, np, nq, cols, mycol, myrow, npcol, nprow, nomega
    real(mq) :: err(1:6), pref, norm_v0(1:6), zz, z1
    real(q)  :: omega, omegamax, shift, conv_thr

    character(len=:), allocatable :: prefix

    complex(q), allocatable :: eps(:,:,:)
    character(len=10), dimension(6) :: eps_dir
# 350

    real(q), parameter  :: prec_dot_prod=1.d-8

    real(q) :: prec_eps
!> debugging option in case bsematrix has a finite imaginary component in the diagonal elements
!   complex(mq), pointer :: bsematrix_p(:) => NULL()

    eps_dir = [character(len=10) :: 'x', 'y', 'z', 'xy', 'xz', 'yz']
    prec_eps=1.d-3

!> because we use VASP native units, the prefactor is the same as in the diagonalisation solver, i.e. 4\pi e^2/volume
    pref = edeps/(omega*real(maxk,mq))
    if (ispin==1) pref=2._mq*pref

    allocate(vector_0(1:fulldim,1:6),vector_1(1:fulldim,1:6),a(1:fulldim,1:6,1:6),b(0:fulldim,1:6,1:6))
    allocate(v_tmp_1(1:matdim_node(node_me),1:6),v_tmp_2(1:fulldim,1:6))
!> the a and b coefficients have to be stored in order to reproduce the full spectra using the continued fraction
!> but the maximum rank of their array is the rank of the bse matrix, as the number of iterations cannot exceed the
!> matrix's rank

    vector_0 = (0._q,0._q); vector_1 = (0._q,0._q)
    v_tmp_1 = (0._q,0._q); ortho_test = (0._q,0._q); v_tmp_2 = (0._q,0._q)

    a=0._mq; b=0._mq

! spectra1 and spectra2 do not need to be store on each GPU so they are kept at CPU level
! here I allocate only them for the sparser frequency grid that will be used for convergence
!> frequency grid for convergence
    nomega = int(sqrt(real(nedos_in)))
    allocate(spectra1(1:nomega,1:6),spectra2(1:nomega,1:6))
!> the spectra is always a complex array
    spectra1 = czero; spectra2 = czero
    
    norm_v0 = 0._mq

    vector_0(1:fulldim,1:3) = optmat(1:fulldim,1:3)

!> debug options in case bsematrix has a non-(0._q,0._q) imaginary component in the diagonal elements
!   call c_f_pointer(c_loc(bsematrix), bsematrix_p, [size(bsematrix)])
!   call hermicity(fulldim, bsematrix_p, bse_desc)

!> off-diagonal components xy, yz, zx, in that order
     vector_0(1:fulldim,4) = (optmat(1:fulldim,1) + optmat(1:fulldim,2))/sqrt(2._mq)  !> \epsilon_{xy}
     vector_0(1:fulldim,5) = (optmat(1:fulldim,2) + optmat(1:fulldim,3))/sqrt(2._mq)  !> \epsilon_{yz}
     vector_0(1:fulldim,6) = (optmat(1:fulldim,1) + optmat(1:fulldim,3))/sqrt(2._mq)  !> \epsilon_{zx}

!> compute the norm along each direction and normalise the vector
    do id1 = 1, 6
# 401

      norm_v0(id1) = sqrt(real(dot_product(vector_0(1:fulldim,id1), vector_0(1:fulldim,id1)),mq))

      vector_1(1:fulldim,id1) = vector_0(1:fulldim,id1)/norm_v0(id1)
    enddo

    call blacs_gridinfo(vec_desc(2), nprow, npcol, myrow, mycol)
    np = numroc(fulldim,bse_desc(5), myrow, 0, nprow)
    nq = numroc(fulldim,bse_desc(6), mycol, 0, npcol)

    cols = 6
    vec_desc=bse_desc
    dvec_desc=bse_desc
    vec_desc(4)=cols
    vec_desc(6)=cols
    vec_desc(3)=fulldim
    vec_desc(5)=fulldim
    if(my_comm%node_me == 1) then
       vec_desc(9)=fulldim
    else
       vec_desc(9)=1
    endif
    dvec_desc(4)=cols
    dvec_desc(6)=cols

!> get time from cpu
    if (node_me==0) then
      call cpu_time(zz)
      if (io%iu0>=0) write(io%iu0,'(a,f8.3,a)') ' Lanczos setup complete after: ', zz, ' sec'
      if (io%iu6>=0) write(io%iu6,'(a,f8.3,a)') ' Lanczos setup complete after: ', zz, ' sec'
      if (io%iu0>=0) write(io%iu0,'(a)') ' Starting Lanczos loop '
      if (io%iu6>=0) write(io%iu6,'(a)') ' Starting lanczos loop '
    endif

!> starting the loop
    it = 1
    itercount = 1
    do while (it <= fulldim)

!> temporary vectors, set to czero for precaution
      v_tmp_1 = (0._q,0._q); v_tmp_2 = (0._q,0._q)

! |v_tmp_1> = bsematrix|vector_1>
      if (matdim_node(node_me)>0) &
      call ZGEMM('n', 'n', matdim_node(node_me), 6, fulldim, (1._q,0._q), bsematrix, matdim_node(node_me), vector_1(1,1), fulldim, &
&                    (0._q,0._q), v_tmp_1(1,1), matdim_node(node_me))

!> for 1 v_tmp_1 (disordered) needs to be copied into vector_1 (ordered) (copy+order all to taks 1, then broadcast from task 1 to all)
      call PZGEMR2D(fulldim, 6, v_tmp_1, 1, 1, dvec_desc, v_tmp_2, 1, 1, vec_desc, dvec_desc(2))
      call M_bcast_z_from(MY_COMM, v_tmp_2, size(v_tmp_2), 1)

! a(it) = <vector_1|bsematrix|vector_1> = <vector_1|v_tmp_1>
      if (matdim_node(node_me)>0) &
      call ZGEMM('C','n', 6, 6, fulldim, (1._q,0._q), vector_1(1,1), fulldim, v_tmp_2(1,1), fulldim, (0._q,0._q), ortho_test(1,1), 6)
      a(it,1:6,1:6) = real(ortho_test(1:6,1:6),mq)

! v_tmp_1 <- bsematrix|vector_1> - a(it)|vector_1> - b(it-1)|v0> = |v_tmp_1> - a(it)|vector_1> - b(it-1)|vector_0>
      do id1 = 1,6
        v_tmp_2(:,id1) = v_tmp_2(:,id1) - a(it,id1,id1)*vector_1(:,id1)
        if (it > 1) then
          v_tmp_2(:,id1) = v_tmp_2(:,id1) - sqrt(b(it-1,id1,id1))*vector_0(:,id1)
        endif
      enddo

!> b(it)^2 = <v_tmp_2|v_tmp_2>
      if (matdim_node(node_me)>0) &
      call ZGEMM('C','n', 6, 6, fulldim, (1._q,0._q), v_tmp_2(1,1), fulldim, v_tmp_2(1,1), fulldim, (0._q,0._q), ortho_test(1,1), 6)
      b(it,1:6,1:6) = real(ortho_test(1:6,1:6),mq)

!> |vector_0> <- |vector_1>
      vector_0 = vector_1

!> |vector_1> <- |v_tmp_1>/sqrt(b(it))
!> old vector_1 is copied to vector_0, but now I need to reassign the distributed vector v_tmp_1 to vector_1
!> seems easier to (0._q,0._q) vector_1, reassing the column-stripes contained in v_tmp_1, and then summ and bcast them
      vector_1 = (0._q,0._q)
      do id1=1,6
        vector_1(:,id1) = v_tmp_2(:,id1)/sqrt(b(it,id1,id1))
      enddo 

!> orthogonality test between |vector_0> and |vector_1>, <vector_0|vector_1>
      if (matdim_node(node_me)>0) &
      call ZGEMM('C','n', 6, 6, fulldim, (1._q,0._q), vector_0(1,1), fulldim, vector_1(1,1), fulldim, (0._q,0._q), ortho_test(1,1), 6)
      do id1 = 1,6
        if (abs(ortho_test(id1,id1)) > prec_dot_prod) then
          call vtutor%error("At iteration " // str(it) // " the orthogonality condition between the two consequtive Lanczos vector, |0> and |1>, was &
               broken. It is not certain that this is recoverable, so the code will now stop! &
               Most likely this comes from having a non-hermitian BSE hamiltonian. Please try running the DFT steps with &
               LORBITALREAl=.TRUE. in the INCAR.")
        endif
      enddo

      if (mod(it,10) == 0) then
!> lanczos spectra at current iteration
        call lanczos_spectra(a(1:fulldim,1:6,1:6), b(1:fulldim,1:6,1:6), it  , shift, 0, omegamax, nomega, spectra2, norm_v0, 0)
!> lanczos spectra at previous iteration
        call lanczos_spectra(a(1:fulldim,1:6,1:6), b(1:fulldim,1:6,1:6), it-1, shift, 0, omegamax, nomega, spectra1, norm_v0, 0)
!> rms between spectra(it)=spectra2 and spectra(it-1) = spectra1
        call spectra_convergence(spectra1(1:nomega,1:6), spectra2(1:nomega,1:6), nomega, err(1:6))        
!> if error condiction is satisfied set make the cycle break
        if (sum(err(1:3))/3._q < prec_eps) then
          it = fulldim
        endif
      endif

!> write to OUTCAR and stdout only every 100 iterations
      if (mod(it,100) == 0) then
!> get time
        if (node_me == 0) call cpu_time(z1)
!> write information about convergence to OUTCAR and stdout
        if (io%iu0>=0) write(io%iu0,'(a,i6,a,e14.8,a,f8.3,a)') ' Lanczos iteration ', it ,' done: RMS ', sum(err(1:3))/3._q, '; Time ', z1, ' sec'
        if (io%iu6>=0) write(io%iu6,'(a,i6,a,e14.8,a,f8.3,a)') ' Lanczos iteration ', it ,' done: RMS ', sum(err(1:3))/3._q, '; Time ', z1, ' sec'
      endif

      it = it + 1
      itercount = itercount + 1

    enddo !> it

!> since the next iteration is not 1._q once convergence is reached
    itercount = itercount - 1

!> deallocate and realocate spectra1 and spectra2 to the full frequency grid
    deallocate(spectra1,spectra2)
    allocate(spectra1(1:nedos_in,1:6),spectra2(1:nedos_in,1:6))
    spectra1 = czero; spectra2 = czero

!> spectra1 is used to compute the resonant part
    call lanczos_spectra(a(1:fulldim,1:6,1:6), b(1:fulldim,1:6,1:6), itercount, shift, 0, omegamax, nedos_in, spectra1(1:nedos_in,1:6), norm_v0, 0)
!> now spectra2 is used to compute the anti-resonant part
    call lanczos_spectra(a(1:fulldim,1:6,1:6), b(1:fulldim,1:6,1:6), itercount, shift, 0, omegamax, nedos_in, spectra2(1:nedos_in,1:6), norm_v0, 1)
!> add the anti-resonant to the resonant part, multiply by the prefactor
    spectra1 = pref*(spectra1 + conjg(spectra2))

!> diagonal components are dealt with here
    allocate(eps(1:nedos_in,1:3,1:3))
    eps= czero
    
    do id1 = 1,3
      eps(1:nedos_in,id1,id1) = 1._q + spectra1(1:nedos_in,id1)
    enddo !> id1

!> k-index for spectra1(:,k) is
!> id1 = 1, id2 = 2 => k = 4
!> id1 = 1, id2 = 3 => k = 6
!> id1 = 2, id2 = 3 => k = 5
!> I rotate the auxilar \bar\epsilon here back to the xy, xz, and yz axis according to
!> \epsilon_{ij} = \bar\epsilon_{ij} - 0.5[\epsilon_{ii} + \epsilon_{jj} - 2], i\neq j

    eps(:,1,2) = spectra1(:,4) - 0.5_q*(eps(:,1,1) + eps(:,2,2) - 2._q)
    eps(:,1,3) = spectra1(:,6) - 0.5_q*(eps(:,1,1) + eps(:,3,3) - 2._q)
    eps(:,2,3) = spectra1(:,5) - 0.5_q*(eps(:,2,2) + eps(:,3,3) - 2._q)
    eps(:,2,1) = eps(:,1,2); eps(:,3,2) = eps(:,2,3); eps(:,3,1) = eps(:,1,3)

!> write the dielectric function in the xml file
    call xml_epsilon_w(omegamax/(nedos_in-1),real(eps),aimag(eps),nedos_in)

    deallocate(vector_0)
    deallocate(vector_1) 
    deallocate(v_tmp_1)
    deallocate(v_tmp_2)
    deallocate(spectra1,spectra2)
    
    call xml_lanczos_w(itercount,a(1:itercount,1:6,1:6),b(1:itercount,1:6,1:6))


    if (lgwlf) then
      prefix = "bse"
    else
      prefix = "tdhf"
    end if
    call vh5_write_dielectric_energies(ih5outfileid, omegamax/(nedos_in-1), nedos_in, prefix)
    call vh5_write_dielectric_dynamic(ih5outfileid, eps, prefix)


    deallocate(eps)

    if (io%iu6>0) then
      write(io%iu6,*) ""
      write(io%iu6,*) "======================================================================="
      write(io%iu6,*) ""
      write(io%iu6,*) "Parameters set for Lanczos algorithm"
      write(io%iu6,'(a,f8.3)') "   Threshold: ", prec_eps
      write(io%iu6,'(a,f8.3)') "   Maximum energy: ", omegamax
      write(io%iu6,*) "   Number of energy steps: ", nedos_in
      write(io%iu6,*) "   Number of iterations needed to reach set accuracy:", itercount
      write(io%iu6,*) "   Final accuracy along all directions:"
      do id1 = 1, 6
        write(io%iu6,'(2a,e14.7)') "      ", eps_dir(id1), err(id1)
      enddo
      write(io%iu6,*) ""
      write(io%iu6,*) "======================================================================="
    endif

    deallocate(a,b)

  end subroutine lanczos_solver


!****************** subroutine lanczos_spectra    *****************************
!
!> computes the spectra from haydock-lanczos method using a continued fraction
!
!******************************************************************************

  subroutine lanczos_spectra(a, b, iter,  shift, term, omegamax, nomega, spectra, u0, res)
    implicit none

    real(q), intent(in) :: a(:,:,:), b(:,:,:)
    real(mq), intent(in) :: u0(:) 
    integer, intent(in) :: iter, term, res, nomega
    real(q),intent(in) :: omegamax, shift
    complex(mq), intent(out) :: spectra(:,:)
!> workspace
    integer  :: iw, it, id1,id2
    real(q) :: deltaw
!> terminator parameters
    complex(mq) :: c1, c2, f, g

    c1 = czero
    c2 = czero
    f = czero
    g = czero

    deltaw = omegamax/(nomega-1)
    spectra = czero
    do id1 = 1,6
      do iw = 0, nomega-1
        if (term > 0) then
          c1 = sum(a(1:iter,id1,id1))/iter
          c1 = sum(b(1:iter,id1,id1))/iter
          f = (iw*deltaw)**2 - c1**2 + c2**2
          g = 2.0*iw*deltaw*c2**2
          spectra(iw+1,id1) = (f**2 - 2.0*iw*deltaw*g)/g
        endif
        do it = 0, iter-1, 1
          if (res == 0) spectra(iw+1,id1) =   1.0/( iw*deltaw - a(iter-it,id1,id1) + imun*shift - b(iter-it,id1,id1)*spectra(iw+1,id1))
          if (res == 1) spectra(iw+1,id1) =   1.0/(-iw*deltaw - a(iter-it,id1,id1) + imun*shift - b(iter-it,id1,id1)*spectra(iw+1,id1))
        enddo
      enddo
      spectra(:,id1) = -u0(id1)**2*spectra(:,id1)
    enddo
  end subroutine lanczos_spectra

!****************** subroutine spectra_convergence    *************************
!
!> computes the RMS error between two different iterations of the spectra from haydock-lanczos
!
!******************************************************************************
  subroutine spectra_convergence(spectra1, spectra2, nomega, error)
    complex(mq), intent(in) :: spectra1(:,:), spectra2(:,:)
    integer, intent(in) :: nomega
    real(mq) :: error(:)
!> workspace
    integer :: iw, id

    error = 0._mq 
    do id = 1,6
      do iw = 1, nomega
        error(id) = error(id) + (aimag(spectra1(iw,id)) - aimag(spectra2(iw,id)))**2
      enddo
      error(id) = sqrt(error(id)/nomega)
    enddo

  end subroutine spectra_convergence

!****************** subroutine hermiticity    *********************************
!
!> enforces hermiticity of a matrix by removing the non-(0._q,0._q) imaginary part in the diagonal elements
!
!******************************************************************************

  subroutine hermicity(n, a, desca ) 

      implicit none
      integer,parameter :: ctxt_=2
!> the blocking factor used to distribute the rows of the array.
      integer,parameter :: mb_=5
!> the blocking factor used to distribute the columns of the array.
      integer,parameter :: nb_=6
!> the leading dimension of the local array.  lld_a >= max(1,locr(m_a)).
      integer,parameter :: lld_=9
!> dimension of desca
      integer,parameter :: dlen_=9

      integer desca( dlen_ ) !< distributed matrix descriptor array
      complex(q) a(:)        !< distributed matrix
      integer n              !< dimension of the matrix
! local
      integer np,nq
      integer i1res,j1res,jcol
      integer myrow, mycol, nprow, npcol
      integer i1,i2,j1,j2

      intrinsic  min
      integer    numroc,irow_jcol
      external   numroc,blacs_gridinfo


      call blacs_gridinfo(desca(ctxt_),nprow,npcol,myrow,mycol)
      np = numroc(n,desca(mb_),myrow,0,nprow)
      nq = numroc(n,desca(nb_),mycol,0,npcol)

      jcol = 0
! loop over local columns
      do j1=1,nq,desca(nb_)
        j1res=min(desca(nb_),nq-j1+1)
        do j2=1,j1res
          irow_jcol = jcol

! loop over local rows
          do i1=1,np,desca(mb_)
            i1res=min(desca(mb_),np-i1+1)
            do i2=1,i1res
              irow_jcol = irow_jcol+1
              if(  (desca(mb_)*myrow+nprow*(i1-1)+i2) == &
                   (desca(nb_)*mycol+npcol*(j1-1)+j2) ) then
                a(irow_jcol) = real(a(irow_jcol),mq)
                goto 100 ! next column, since we have found the right row
              endif
            enddo
          enddo
100       jcol = jcol + desca(lld_)
        enddo
      enddo

      return

  end subroutine hermicity

end module bse_lanczos

