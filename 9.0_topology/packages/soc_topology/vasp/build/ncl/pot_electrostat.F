#include "symbol.inc"

module pot_electrostat

GPACC use moffload_struct_def
GPACC use moffload
      use prec
      use mpimy
      use mgrid
      use lattice
      use constant
      use pot_struct_def
      use coulomb_cutoff
      use pseudo_struct_def, only: NPSPTS, potcar
      use poscar_struct_def, only: type_info
      use charge, only: rhoadd

      implicit none

contains

    subroutine pothar(gridc, latt_cur, coulomb_pot, charge_density, &
                      potential, denc)
        type(grid_3d), intent(in)            :: gridc
        type(latt), intent(in)               :: latt_cur
        type(coulomb_potential), intent(in)  :: coulomb_pot
        complex(q)                           :: charge_density(gridc%rc%np)
        complex(q), intent(out)              :: potential(gridc%rc%np)
        real(q), intent(out)                 :: denc

        PROFILING_START('pothar')
        !
        if (coulomb_pot%kernel_truncate%active) then
            !
            call cutoff_hartree_potential(gridc, latt_cur, &
                                          coulomb_pot, &
                                          charge_density, potential)
        else
            !
            call potential_from_charge_density(gridc, latt_cur, &
                                               coulomb_pot, &
                                               charge_density, &
                                               potential)
        end if
        !
        call electrostatic_energy(gridc, charge_density, potential, &
                                  denc)
        !
        PROFILING_STOP('pothar')
    end subroutine

    subroutine potion(gridc, p, latt_cur, t_info, coulomb_pot, &
                      ion_potential, structure_factor, pscenc, &
                      over_ions)
      type (grid_3d), intent(in)           :: gridc
      type (type_info), intent(in)         :: t_info
      type (potcar), intent(in)            :: p (t_info%ntyp)
      type (latt), intent(in)              :: latt_cur
      type (coulomb_potential), intent(in) :: coulomb_pot
      complex(q), intent(in)               :: structure_factor(gridc%mplwv,t_info%ntyp)
      complex(q), intent(out)              :: ion_potential(gridc%rc%np)
      real(q), intent(out)                 :: pscenc       ! non-electrostatic contributions to energy
      integer, optional                    :: over_ions(:) ! if less than the total number of ions is required

      call pscenc_energy(p,t_info,latt_cur,coulomb_pot,pscenc)
      !
      associate(is_truncate => coulomb_pot%kernel_truncate%active)
      !
      !$acc kernels present(ion_potential) __IF_ASYNC__
      ion_potential = 0._q
      !$acc end kernels
      !
      if (is_truncate) then
         !
         call cutoff_ion_potential(p, t_info, gridc, latt_cur, &
                                   coulomb_pot, structure_factor, &
                                   ion_potential, &
                                   use_grad_kernels=.false., &
                                   over_ions=over_ions)
      else
         !
         call potential_from_ion_charge(p, t_info, gridc, latt_cur, &
                                        coulomb_pot, structure_factor, &
                                        ion_potential, &
                                        use_grad_kernels=.false., &
                                        over_ions=over_ions)
      end if
      end associate
    end subroutine

    subroutine grad_potion(gridc, p, latt_cur, t_info, coulomb_pot, &
                           grad_ion_potential, structure_factor)
      type (grid_3d), intent(in)           :: gridc
      type (type_info), intent(in)         :: t_info
      type (potcar), intent(in)            :: p (t_info%ntyp)
      type (latt), intent(in)              :: latt_cur
      type (coulomb_potential), intent(in) :: coulomb_pot
      complex(q), intent(in)               :: structure_factor(gridc%mplwv,t_info%ntyp)
      complex(q), intent(out)              :: grad_ion_potential(gridc%rc%np)

      associate(is_truncate => coulomb_pot%kernel_truncate%active)
      !
      !$acc kernels present(grad_ion_potential) __IF_ASYNC__
      grad_ion_potential = 0._q
      !$acc end kernels
      !
      if (is_truncate) then
         !
         call cutoff_ion_potential(p, t_info, gridc, latt_cur, &
                                   coulomb_pot, structure_factor, &
                                   grad_ion_potential, &
                                   use_grad_kernels=.true.)
      else
         !
         call potential_from_ion_charge(p, t_info, gridc, latt_cur, &
                                   coulomb_pot, structure_factor, &
                                   grad_ion_potential, &
                                   use_grad_kernels=.true.)
      end if
      end associate
      end subroutine

      subroutine potential_from_charge_density(gridc, latt_cur, &
                                               coulomb_pot, &
                                               charge_density, &
                                               potential)
          type (grid_3d), intent(in)           :: gridc
          type (latt), intent(in)              :: latt_cur
          type(coulomb_potential), intent(in)  :: coulomb_pot
          complex(q), intent(in)               :: charge_density(gridc%rc%np)
          complex(q), intent(inout)            :: potential(gridc%rc%np)
          integer                              :: idx(3)
          real(q)                              :: g_vector(3)
          real(q)                              :: kernel
          integer                              :: n1, n2, n3, ni, nc

!$acc routine (index_to_gvector) seq
!$acc routine (kernel_for_gvector) seq
!$acc routine (potential_from_kernel_for_gvector) seq
          !
!$acc parallel loop collapse(2) private(n2,n3,ni,idx,g_vector,kernel) &
!$acc present(gridc,latt_cur,coulomb_pot,charge_density,potential) __IF_ASYNC__
          do nc=1,gridc%rc%ncol
__CPU     n2= gridc%rc%i2(nc)
__CPU     n3= gridc%rc%i3(nc)
          do n1=1,gridc%rc%nrow
            !
__GPU       n2= gridc%rc%i2(nc)
__GPU       n3= gridc%rc%i3(nc)
            !
            ni = (nc - 1) * gridc%rc%nrow + n1
            !
            idx = [gridc%lpctx(n1), gridc%lpcty(n2), gridc%lpctz(n3)]
            !
            call index_to_gvector(latt_cur, idx, g_vector)
            !
            call kernel_for_gvector(coulomb_pot, &
                                    g_vector, kernel, all(idx==0))
            !
            call potential_from_kernel_for_gvector(kernel, &
                                                   charge_density(ni), &
                                                   latt_cur%omega, &
                                                   potential(ni))
          enddo
          enddo
          !
          call setunb(potential,gridc)
      end subroutine

      subroutine potential_from_kernel_for_gvector(kernel, &
                                                        charge_density,&
                                                        volume, &
                                                        potential)
!$acc routine seq
         real(q), intent(in)     :: kernel
         complex(q), intent(in)  :: charge_density
         real(q), intent(in)     :: volume
         complex(q), intent(out) :: potential

         potential = kernel * charge_density / volume
      end subroutine

      subroutine electrostatic_energy(gridc, charge_density, potential,&
                                      energy)
         type(grid_3d), intent(in) :: gridc
         complex(q), intent(in)    :: charge_density(gridc%rc%np)
         complex(q), intent(in)    :: potential(gridc%rc%np)
         real(q), intent(out)      :: energy
         integer :: nc, n1, n2, n3, ni
         real(q) :: factm, dum

         energy = 0._q
         !
         !$acc parallel loop collapse(2) reduction(+:energy) private(n2,n3,ni,dum,factm) &
         !$acc present(gridc,potential,charge_density) __IF_ASYNC__
         do nc = 1, gridc%rc%ncol
__CPU    n2 = gridc%rc%i2(nc)
__CPU    n3 = gridc%rc%i3(nc)
         do n1 = 1, gridc%rc%nrow
__GPU      n2 = gridc%rc%i2(nc)
__GPU      n3 = gridc%rc%i3(nc)
           !
           ni = (nc - 1) * gridc%rc%nrow + n1
           SETFACT1
           SETFACT
           !
           dum = MULFACT potential(ni) * conjg(charge_density(ni))
           energy = energy + dum
         enddo
         enddo
         !
         energy = -energy / 2._q
         CALLMPI( M_sum_d(gridc%comm, energy, 1))
      end subroutine

      subroutine iteration_information(iterating_over_ions, &
                                       iterate_over_array, &
                                       total_number_types, &
                                       over_ions)
         logical, intent(out)              :: iterating_over_ions
         integer, allocatable, intent(out) :: iterate_over_array(:)
         integer, intent(in)               :: total_number_types
         integer, optional, intent(in)     :: over_ions(:)
         integer                           :: nt

         if (present(over_ions)) then
            !
            iterating_over_ions = .true.
         else
            !
            iterating_over_ions = .false.
         end if
         !
         if (iterating_over_ions) then
             !
             allocate(iterate_over_array, source=over_ions)
         else
             !
             allocate(iterate_over_array(total_number_types))
             !
             do nt = 1, total_number_types
                !
                iterate_over_array(nt) = nt
             end do
         end if
      end subroutine

      subroutine potential_from_ion_charge(p, t_info, gridc, latt_cur, &
                                           coulomb_pot, &
                                           structure_factor, &
                                           potential, &
                                           use_grad_kernels, &
                                           over_ions)
         type(grid_3d), intent(in)           :: gridc
         type(type_info), intent(in)         :: t_info
         type(potcar), intent(in)            :: p (t_info%ntyp)
         type (latt), intent(in)             :: latt_cur
         type(coulomb_potential), intent(in) :: coulomb_pot
         complex(q), intent(in)              :: structure_factor(gridc%mplwv,t_info%ntyp)
         complex(q), intent(out)             :: potential(gridc%rc%np)
         logical, intent(in)                 :: use_grad_kernels
         integer, optional                   :: over_ions(:)
         integer                             :: nc, ni, n1, n2, n3
         integer                             :: j, idx_type
         real(q)                             :: g_vector(3)
         real(q)                             :: g_norm
         real(q)                             :: psgma2, zz
         real(q)                             :: kernel, ps_kernel
         real(q)                             :: summed_kernel
         integer, allocatable                :: iterate_over(:)
         integer                             :: idx(3)
         logical                             :: iterating_over_ions
         complex(q)                          :: charge
         complex(q)                          :: potential_
         logical                             :: is_g0
!$acc routine (index_to_gvector) seq
!$acc routine (kernel_for_gvector) seq
!$acc routine (kernel_gradient_for_gvector) seq
!$acc routine (pseudo_kernel_for_gvector) seq
!$acc routine (pseudo_kernel_gradient_for_gvector) seq
!$acc routine (phase_factor_for_index) seq
!$acc routine (potential_from_kernel_for_gvector) seq

         !$acc enter data copyin(t_info) __IF_ASYNC__
         !$acc enter data copyin(t_info%vca,t_info%posion) __IF_ASYNC__
         !
         !$acc kernels present(potential) __IF_ASYNC__
         potential = 0._q
         !$acc end kernels
         !
         call iteration_information(iterating_over_ions, iterate_over, &
                                    t_info%ntyp, over_ions)
         !
         do j = 1, size(iterate_over)
           !
           if (iterating_over_ions) then
               !
               idx_type = t_info%ityp(iterate_over(j))
           else
               !
               idx_type = iterate_over(j)
           end if
           !
           psgma2 = p(idx_type)%psgmax - p(idx_type)%psgmax / npspts
           !
           zz = -p(idx_type)%zvalf
           !
           !
!$acc parallel loop collapse(2) private(n2,n3,ni,idx,g_vector,g_norm,is_g0,charge,kernel,ps_kernel,summed_kernel,potential_) &
!$acc& present(gridc,p,t_info,latt_cur,coulomb_pot,structure_factor,potential) __IF_ASYNC__
           do nc = 1, gridc%rc%ncol
__CPU      n2 = gridc%rc%i2(nc)
__CPU      n3 = gridc%rc%i3(nc)
           do n1 = 1, gridc%rc%nrow
__GPU        n2 = gridc%rc%i2(nc)
__GPU        n3 = gridc%rc%i3(nc)
             !
             ni = (nc - 1) * gridc%rc%nrow + n1
             !
             idx = [gridc%lpctx(n1), gridc%lpcty(n2), gridc%lpctz(n3)]
             !
             is_g0 = all(idx==0)
             !
             call index_to_gvector(latt_cur, idx, g_vector)
             !
             g_norm = sqrt(  g_vector(1)**2._q &
                           + g_vector(2)**2._q &
                           + g_vector(3)**2._q )
             !
             if ( g_norm < psgma2 )  then
                 !
                 if (use_grad_kernels) then
                    !
                     call kernel_gradient_for_gvector(&
                                     coulomb_pot, &
                                     g_vector, kernel, is_g0)
                     !
                     call pseudo_kernel_gradient_for_gvector(&
                                     p, t_info, coulomb_pot, &
                                     g_vector, idx_type, ps_kernel, &
                                     is_g0)
                else
                     !
                     call kernel_for_gvector(&
                                     coulomb_pot, &
                                     g_vector, kernel, is_g0)
                     !
                     call pseudo_kernel_for_gvector(&
                                     p, t_info, coulomb_pot, &
                                     g_vector, idx_type, ps_kernel, &
                                     is_g0)
                end if
                !
                if (iterating_over_ions) then
                   !
                   call phase_factor_for_index(idx, &
                                   t_info%posion(:,iterate_over(j)), &
                                   t_info%vca(idx_type), &
                                   charge)
                else
                  !
                   charge = structure_factor(ni, idx_type)
                end if
                !
                summed_kernel = ps_kernel + zz * kernel
                !
                call potential_from_kernel_for_gvector(summed_kernel, &
                                                       charge, &
                                                       latt_cur%omega,&
                                                       potential_)
                potential(ni) = potential(ni) + potential_
             else
                potential(ni) = 0._q
             endif
           enddo
           enddo
           enddo
           !$acc exit data delete(t_info) __IF_ASYNC__
           !$acc exit data delete(t_info%vca,t_info%posion) __IF_ASYNC__
           call setunb(potential, gridc)
      end subroutine

!> calculate the contribution to the total energy from the non-coulomb
!> part of the g=0 component of the pseudopotential and the force on the
!> unit cell due to the change in this energy as the size of the cell
!> changes
      subroutine pscenc_energy(p, t_info, latt_cur, coulomb_pot, pscenc)
          type(type_info), intent(in)         :: t_info
          type(potcar), intent(in)            :: p (t_info%ntyp)
          type(latt), intent(in)              :: latt_cur
          type(coulomb_potential), intent(in) :: coulomb_pot
          real(q), intent(out)                :: pscenc
          real(q)                             :: zvsum
          integer :: nt

          if (coulomb_pot%kernel_truncate%active) then
              !
              pscenc = 0._q
          else
              !
              zvsum=0
              do nt=1,t_info%ntyp
                 !
                 zvsum=zvsum+p(nt)%zvalf*t_info%nityp(nt)*t_info%vca(nt)
              enddo
              !
              pscenc = 0._q
              do nt=1,t_info%ntyp
                 !
                 pscenc=pscenc + p(nt)%pscore * t_info%vca(nt) &
                                              * t_info%nityp(nt) &
                                              * (zvsum / latt_cur%omega)
              enddo
          end if
      end subroutine

     subroutine cutoff_hartree_potential(gridc, latt_cur, coulomb_pot, &
                                         charge_density, potential)
        type(grid_3d), intent(in)            :: gridc
        type(latt), intent(in)               :: latt_cur
        type(coulomb_potential), intent(in)  :: coulomb_pot
        complex(q)                           :: charge_density(gridc%rc%np)
        complex(q), intent(out)              :: potential(gridc%rc%np)
        type(latt)                           :: latt_cur_padded
        type(coulomb_potential)              :: periodic_coulomb_pot
        complex(q), allocatable              :: padded_charge_density(:)
        complex(q), allocatable              :: coarse_charge_density(:)
        complex(q), allocatable              :: padded_potential(:)
        complex(q), allocatable              :: coarse_open_potential(:)
        complex(q), allocatable              :: coarse_periodic_potential(:)
        complex(q), allocatable              :: periodic_potential(:)
        real(q)                              :: rescale

        call rescaling_factor(coulomb_pot%kernel_truncate, rescale)
        !
        associate(kernel_truncate    => coulomb_pot%kernel_truncate, &
                  coarsen_before_pad => coulomb_pot%kernel_truncate%coarsen_before_pad, &
                  fine_grid          => coulomb_pot%kernel_truncate%fine_grid, &
                  coarse_grid        => coulomb_pot%kernel_truncate%coarse_grid, &
                  translation_table  => coulomb_pot%kernel_truncate%translation_table, &
                  padding_factor     => coulomb_pot%kernel_truncate%padding_factor)
        !
GPACC   call copyin_typed_var(periodic_coulomb_pot)
        !
        allocate(padded_charge_density(fine_grid%mplwv))
        allocate(padded_potential(fine_grid%mplwv))
        !
        padded_charge_density = 0._q
        padded_potential = 0._q
        !
        !$acc enter data copyin(padded_charge_density) __IF_ASYNC__
        !$acc enter data copyin(padded_potential) __IF_ASYNC__
        !
        call generate_padded_lattice(kernel_truncate, &
                                     latt_cur, latt_cur_padded)
        !
        !$acc enter data copyin(latt_cur_padded) __IF_ASYNC__
        !
        if (coarsen_before_pad) then
            !
            !$acc kernels present(potential) __IF_ASYNC__
            potential = 0._q
            !$acc end kernels
            !
            allocate(coarse_charge_density(coarse_grid%mplwv))
            coarse_charge_density = 0._q
            !$acc enter data copyin(coarse_charge_density) __IF_ASYNC__
            !
            allocate(coarse_open_potential(coarse_grid%mplwv))
            coarse_open_potential = 0._q
            !$acc enter data copyin(coarse_open_potential) __IF_ASYNC__
            !
            allocate(coarse_periodic_potential(coarse_grid%mplwv))
            coarse_periodic_potential = 0._q
            !$acc enter data copyin(coarse_periodic_potential) __IF_ASYNC__
            !
            allocate(periodic_potential(gridc%mplwv))
            periodic_potential = 0._q
            !
            !$acc enter data copyin(periodic_potential) __IF_ASYNC__
            call cp_grid(coulomb_pot%kernel_truncate%fine_grid, &
                         coulomb_pot%kernel_truncate%coarse_grid, &
                         coulomb_pot%kernel_truncate%translation_table, &
                         charge_density, coarse_charge_density)
            !
            call setunb(coarse_charge_density, &
                        coulomb_pot%kernel_truncate%coarse_grid)
            !
            call pad3d(kernel_truncate, &
                       coarse_charge_density, &
                       padded_charge_density, 1, rescale)
            !
            call setunb(padded_charge_density, &
                        coulomb_pot%kernel_truncate%fine_grid)
            !
            call potential_from_charge_density(fine_grid,&
                                latt_cur_padded,coulomb_pot,&
                                padded_charge_density, padded_potential)
            !
            call pad3d(kernel_truncate, coarse_open_potential, &
                       padded_potential, -1)
            !
            call setunb(coarse_open_potential, coarse_grid)
            !
            call potential_from_charge_density(coarse_grid, latt_cur, &
                                periodic_coulomb_pot, &
                                coarse_charge_density, &
                                coarse_periodic_potential)
            !
            call rc_add(coarse_open_potential, 1._q, &
                        coarse_periodic_potential, -1._q, &
                        coarse_open_potential, coarse_grid)
            !
            call cpb_grid(fine_grid, coarse_grid, translation_table, &
                          coarse_open_potential, potential)
            !
            call setunb(potential, fine_grid)
            !
            call potential_from_charge_density(gridc, latt_cur, &
                                periodic_coulomb_pot, charge_density, &
                                periodic_potential)
            !
            call rc_add(potential, 1._q, &
                        periodic_potential, 1._q, &
                        potential, gridc)
            !$acc exit data delete(coarse_charge_density) __IF_ASYNC__
            !$acc exit data delete(coarse_open_potential) __IF_ASYNC__
            !$acc exit data delete(coarse_periodic_potential) __IF_ASYNC__
            !$acc exit data delete(periodic_potential) __IF_ASYNC__
        else
            !
            call pad3d(kernel_truncate, charge_density, &
                       padded_charge_density, 1, rescale)
            !
            call setunb(padded_charge_density, fine_grid)
            !
            call potential_from_charge_density(fine_grid, &
                                latt_cur_padded,coulomb_pot, &
                                padded_charge_density,padded_potential)
            !
            call pad3d(kernel_truncate, potential, padded_potential, -1)
        end if
        !$acc exit data delete(padded_charge_density) __IF_ASYNC__
        !$acc exit data delete(padded_potential) __IF_ASYNC__
        !$acc exit data delete(latt_cur_padded) __IF_ASYNC__
        !
        call setunb(potential, gridc)
GPACC   call delete_typed_var(periodic_coulomb_pot)
        end associate
     end subroutine

     subroutine cutoff_ion_potential(p, t_info, gridc, latt_cur, &
                                     coulomb_pot, structure_factor, &
                                     potential, use_grad_kernels, &
                                     over_ions)
        type (grid_3d), intent(in)           :: gridc
        type (type_info), intent(in)         :: t_info
        type (potcar), intent(in)            :: p (t_info%ntyp)
        type (latt), intent(in)              :: latt_cur
        type (coulomb_potential), intent(in) :: coulomb_pot
        complex(q), intent(in)               :: structure_factor(gridc%mplwv,t_info%ntyp)
        complex(q), intent(inout)            :: potential(gridc%rc%np)
        logical, intent(in)                  :: use_grad_kernels
        integer, optional                    :: over_ions(:)
        type(type_info)                      :: t_info_pad
        type(latt)                           :: latt_cur_pad
        complex(q), allocatable              :: structure_factor_pad(:,:)
        complex(q), allocatable              :: structure_factor_coarse(:,:)
        complex(q), allocatable              :: potential_pad(:)
        complex(q), allocatable              :: potential_coarse(:)
        complex(q), allocatable              :: potential_coarse_periodic(:)
        complex(q), allocatable              :: potential_periodic(:)
        type(coulomb_potential)              :: coulomb_pot_periodic
        real(q), allocatable                 :: padded_posion(:,:)
        integer :: nt

        associate(kernel_truncate    => coulomb_pot%kernel_truncate, &
                  coarsen_before_pad => coulomb_pot%kernel_truncate%coarsen_before_pad, &
                  fine_grid          => coulomb_pot%kernel_truncate%fine_grid, &
                  coarse_grid        => coulomb_pot%kernel_truncate%coarse_grid, &
                  translation_table  => coulomb_pot%kernel_truncate%translation_table, &
                  padding_factor     => coulomb_pot%kernel_truncate%padding_factor)
        !
GPACC   call copyin_typed_var(coulomb_pot_periodic)
        !
        allocate(structure_factor_pad(fine_grid%mplwv, t_info%ntyp))
        allocate(structure_factor_coarse(coarse_grid%mplwv, t_info%ntyp))
        allocate(potential_pad(fine_grid%mplwv))
        allocate(potential_coarse(coarse_grid%mplwv))
        allocate(potential_coarse_periodic(coarse_grid%mplwv))
        allocate(potential_periodic(gridc%mplwv))
        !
        call generate_padded_lattice(kernel_truncate, latt_cur, &
                                     latt_cur_pad)
        !$acc enter data copyin(latt_cur_pad) __IF_ASYNC__
        !
        call generate_padded_t_info(kernel_truncate, t_info, &
                                    t_info_pad, padded_posion)
        !
        structure_factor_pad = 0._q
        call stufak(fine_grid, t_info_pad, structure_factor_pad)
        !$acc enter data copyin(structure_factor_pad) __IF_ASYNC__
        !
        potential_pad = 0._q
        !$acc enter data copyin(potential_pad) __IF_ASYNC__
        call potential_from_ion_charge(p, t_info_pad,&
                                       fine_grid,latt_cur_pad, &
                                       coulomb_pot, &
                                       structure_factor_pad, &
                                       potential_pad, &
                                       use_grad_kernels=.false., &
                                       over_ions=over_ions)
        !
        if (coarsen_before_pad) then
            !
            potential_coarse = 0._q
            !$acc enter data copyin(potential_coarse) __IF_ASYNC__
            !
            call pad3d(kernel_truncate, potential_coarse, &
                       potential_pad, -1)
            !
            call setunb(potential_coarse, coarse_grid)
            !
            structure_factor_coarse = 0._q
            !$acc enter data copyin(structure_factor_coarse) __IF_ASYNC__
            !
            do nt = 1, t_info%ntyp
                !
                call cp_grid(fine_grid, coarse_grid, &
                             translation_table, &
                             structure_factor(:,nt), &
                             structure_factor_coarse(:,nt))
                !
                call setunb(structure_factor_coarse(:,nt), &
                            coarse_grid)
            end do
            !
            potential_coarse_periodic = 0._q
            !$acc enter data copyin(potential_coarse_periodic) __IF_ASYNC__
            !
            call potential_from_ion_charge(p, t_info, coarse_grid, &
                                   latt_cur, &
                                   coulomb_pot_periodic, &
                                   structure_factor_coarse,&
                                   potential_coarse_periodic, &
                                   use_grad_kernels=use_grad_kernels, &
                                   over_ions=over_ions)
            !
            call rc_add(potential_coarse, 1.0_q, &
                        potential_coarse_periodic, -1.0_q, &
                        potential_coarse, coarse_grid)
            !$acc kernels present(potential) __IF_ASYNC__
            potential = 0._q
            !$acc end kernels
            !
            call cpb_grid(fine_grid, coarse_grid, translation_table, &
                          potential_coarse, potential)
            !
            call setunb(potential, fine_grid)
            !
            potential_periodic = 0._q
            !$acc enter data copyin(potential_periodic) __IF_ASYNC__
            !
            call potential_from_ion_charge(p, t_info, gridc, latt_cur,&
                                   coulomb_pot_periodic, &
                                   structure_factor,&
                                   potential_periodic, &
                                   use_grad_kernels=use_grad_kernels, &
                                   over_ions=over_ions)
            !
            call rc_add(potential, 1.0_q, &
                        potential_periodic, 1.0_q, &
                        potential, gridc)
            !$acc exit data delete(potential_coarse) __IF_ASYNC__
            !$acc exit data delete(structure_factor_coarse) __IF_ASYNC__
            !$acc exit data delete(potential_coarse_periodic) __IF_ASYNC__
            !$acc exit data delete(potential_periodic) __IF_ASYNC__
        else
            !
            call pad3d(kernel_truncate, potential, potential_pad, -1)
        end if
        !
        !$acc exit data delete(structure_factor_pad) __IF_ASYNC__
        !$acc exit data delete(potential_pad) __IF_ASYNC__
        !$acc exit data delete(latt_cur_pad) __IF_ASYNC__
GPACC   call delete_typed_var(coulomb_pot_periodic)
        !
        call setunb(potential,gridc)
        end associate
   end subroutine

   subroutine potion_particle_mesh(gridc, p, latt_cur, t_info, &
                                   coulomb_pot, potential, &
                                   pscenc, energy, force)
       type(grid_3d), intent(in)           :: gridc
       type(type_info), intent(in)         :: t_info
       type(potcar)                        :: p(t_info%ntyp)
       type(latt), intent(in)              :: latt_cur
       type(coulomb_potential), intent(in) :: coulomb_pot
       real(q), intent(out)                :: pscenc
       real(q), intent(out)                :: energy
       real(q), optional                   :: force(:,:)
       real(q), allocatable                :: qtot(:)
       complex(q), allocatable             :: charge(:)
       complex(q), intent(out)             :: potential(:)
       integer                             :: n_type
       real(q) :: zvsum,g,gx,gy,gz,tewen0
       integer :: nt,nis,ni,n,n1,n2,n3,nc,factm

       allocate(charge(gridc%MPLWV))
       allocate(qtot(t_info%nions))
       !
       call pscenc_energy(p, t_info, latt_cur, coulomb_pot, pscenc)
       !
       charge = 0._q; potential = 0._q; qtot = 0._q
       !
       do n_type = 1, t_info%ntyp
          p(n_type)%usespl => p(n_type)%rhospl
          p(n_type)%usez = -p(n_type)%zvalf
          p(n_type)%usecut = p(n_type)%rcutrho
       end do
       !
       call rhoadd(t_info, latt_cur, p, gridc, charge, qtot)
       !
       call pothar(gridc, latt_cur, coulomb_pot, charge, potential, &
                   energy)
       energy = -1_q * energy - pscenc
       do nt=1,t_info%ntyp
         energy = energy - p(nt)%eself*t_info%vca(nt)*t_info%nityp(nt)
       enddo
       !
       call FFT3D(potential, gridc, 1)
       !
       ! calculate force on each ion
       !
       ! (dE/dR_i,a) = e^2 \int d^3 r (drho_i(|r-R_i|)/dR_i,a) V(r)
       ! i=1..NION
       ! a=1..3 (x,y,z)
       !
       ! drho_i(|r-R_i|)/dR_i,a = -rho'(|r-R_i|) (r_a-R_i,a)/|r-R_i|
       !
       ! Evaluate f_i(r)=d rho/dr*(r-R_i)/|r-R_i| on the grid for each ion
       ! and integrate F_i = \int dr f_i(r) V(r)
       !
       nis = 1
       !
       if (present(force)) then
           type: do nt = 1, t_info%ntyp
              if (.not.associated(p(nt)%usespl)) then
                 nis = nis+t_info%nityp(nt); cycle type
              endif
              ions: do ni=nis,t_info%nityp(nt)+nis-1
                 call rhoder(t_info,latt_cur,p(nt),ni,gridc,qtot,potential,force(1,ni))
              enddo ions
              nis = nis+t_info%nityp(nt)
           enddo type
           force=force*latt_cur%omega
       end if
       ! take the potential back to reciprocal space
       call rl_add(potential,1._q/gridc%nplwv,potential,0.0_q,potential,gridc)
       call FFT3D(potential,gridc,-1)
       call setunb_compat(potential,gridc)
   end subroutine

   subroutine rescaling_factor(kernel_truncate, rescale)
       type(kernel_truncation) :: kernel_truncate
       real(q) :: rescale

       associate(dimensionality => kernel_truncate%dimensionality, &
                 padding_factor => kernel_truncate%padding_factor)
       if (dimensionality == 0) then
           !
           rescale = real(padding_factor,q)**3.0_q
       else if (dimensionality == 2) then
           !
           rescale = real(padding_factor,q)
       end if
       end associate
   end subroutine

end module
