#include "symbol.inc"

module coulomb_cutoff
      use prec, only: q
      use mgrid, only: grid_3d, inilgrd, gen_rc_sub_grid, gen_rc_grid
      use lattice, only: latt
      use pot_struct_def
      use pseudo_struct_def, only: potcar
      use poscar_struct_def, only: type_info
      use mgrid_struct_def, only: transit
      use base, only: info_struct
      use constant, only: EDEPS, PI, TPI, FELECT, CITPI
      use lattice, only: lattic
      use tutor, only: vtutor
      use mpimy, only: communic
      use pseudo_struct_def, only: NPSPTS, potcar
GPACC use moffload

    implicit none

 contains

     subroutine index_to_gvector(latt_cur, idx, g_vector)
!$acc routine seq
      type(latt), intent(in) :: latt_cur
      integer, intent(in)    :: idx(3)
      real(q), intent(out)   :: g_vector(3)

      g_vector(1) = idx(1) * latt_cur%B(1,1) + idx(2) * latt_cur%B(1,2) + idx(3) * latt_cur%B(1,3)
      g_vector(2) = idx(1) * latt_cur%B(2,1) + idx(2) * latt_cur%B(2,2) + idx(3) * latt_cur%B(2,3)
      g_vector(3) = idx(1) * latt_cur%B(3,1) + idx(2) * latt_cur%B(3,2) + idx(3) * latt_cur%B(3,3)
      !
      g_vector = g_vector * TPI
     end subroutine

     subroutine periodic_kernel(g_vector, kernel, calculate_g0)
!$acc routine seq
       real(q), intent(in)     :: g_vector(3)
       real(q), intent(out)    :: kernel
       logical, intent(in)     :: calculate_g0
       real(q)                 :: g_square

       if (calculate_g0) then
          kernel = 0.0_q
       else
          g_square = g_vector(1)**2._q &
                   + g_vector(2)**2._q &
                   + g_vector(3)**2._q
          kernel = 2._q * tpi * felect / g_square
       end if
     end subroutine

     subroutine periodic_kernel_gradient(g_vector, kernel_gradient, &
                                              calculate_g0)
!$acc routine seq
       real(q), intent(in)     :: g_vector(3)
       real(q), intent(out)    :: kernel_gradient
       logical, intent(in)     :: calculate_g0
       real(q)                 :: g_cube
       !
       if (calculate_g0) then
          kernel_gradient = 0.0_q
       else
          g_cube = ( g_vector(1)**2._q &
                   + g_vector(2)**2._q &
                   + g_vector(3)**2._q )**(3._q / 2._q)
          kernel_gradient = -4._q * tpi * felect / g_cube
       end if
     end subroutine

     subroutine spherical_truncated_kernel(kernel_truncate, &
                                                g_vector, kernel, &
                                                calculate_g0)
!$acc routine seq
       type(kernel_truncation), intent(in) :: kernel_truncate
       real(q), intent(in)     :: g_vector(3)
       real(q), intent(out)    :: kernel
       logical, intent(in)     :: calculate_g0
       real(q)                 :: g_norm
       !
       !
       associate(rc => kernel_truncate%truncation_length)
       if (calculate_g0) then
           kernel = tpi * felect * rc**2._q
       else
          g_norm = sqrt(  g_vector(1)**2._q &
                        + g_vector(2)**2._q &
                        + g_vector(3)**2._q )
          kernel = 2._q * tpi * felect / g_norm**2._q &
                 * ( 1._q - cos(g_norm * rc) )
       end if
       end associate
     end subroutine

     subroutine get_surface_index(normal_idx, surface_idx)
!$acc routine seq
         integer, intent(in) :: normal_idx
         integer, intent(out) :: surface_idx(2)
         integer :: i, idx_
         !
         idx_ = 1
#ifdef debug_coulomb_cutoff
         if (normal_idx < 1 .or. normal_idx > 3) then
             call vtutor%error("Index of the surface normal component is &
                       & less than 1 (x-axis) and greater than 3 &
                       & (z-axis).")
         end if
#endif
         do i = 1, 3
            if (i == normal_idx) cycle
            surface_idx(idx_) = i
            idx_ = idx_ + 1
         end do
     end subroutine get_surface_index

     function surface_and_normal_both_zero(g_para, g_perp) &
                                                result(gohere)
!$acc routine seq
         real(q), intent(in) :: g_para, g_perp
         logical :: gohere
         !
         gohere = .false.
         if (g_para < 1e-10_q .and. abs(g_perp) < 1e-10_q) then
            gohere = .true.
         end if
     end function

     function surface_zero_but_normal_nonzero(g_para, g_perp) &
                                                   result(gohere)
!$acc routine seq
         real(q), intent(in) :: g_para, g_perp
         logical :: gohere
         !
         gohere = .false.
         if (g_para < 1e-10_q .and. abs(g_perp) > 1e-10_q) then
            gohere = .true.
         end if
     end function

     function surface_nonzero(g_para, g_perp) result(gohere)
!$acc routine seq
         real(q), intent(in) :: g_para, g_perp
         logical :: gohere
         !
         gohere = .false.
         if (g_para > 1e-10_q) then
             gohere = .true.
         end if
     end function

     subroutine slab_truncated_kernel(kernel_truncate, g_vector, &
                                           kernel, calculate_g0)
!$acc routine seq
       type(kernel_truncation), intent(in) :: kernel_truncate
       real(q), intent(in)     :: g_vector(3)
       real(q), intent(out)    :: kernel
       logical, intent(in)     :: calculate_g0
       integer                 :: surface_index(2)
       real(q)                 :: g_perp, g_para, g_norm

!$acc routine(get_surface_index) seq
!$acc routine(surface_and_normal_both_zero) seq
!$acc routine(surface_zero_but_normal_nonzero) seq
!$acc routine(surface_nonzero) seq
       associate(&
           normal_index => kernel_truncate%surface_normal_direction, &
           rc => kernel_truncate%truncation_length, &
           use_simplified => kernel_truncate%use_simplified_2d_kernel)
           !
           call get_surface_index(normal_index, surface_index)
           !
           g_perp = g_vector(normal_index)
           !
           g_para = sqrt(  g_vector(surface_index(1))**2._q &
                         + g_vector(surface_index(2))**2._q )
           !
           g_norm = sqrt(  g_vector(1)**2._q &
                         + g_vector(2)**2._q &
                         + g_vector(3)**2._q )
           !
           if (surface_and_normal_both_zero(g_para, g_perp) &
               .or. calculate_g0) then
               !
               kernel = -1._q * tpi * felect * rc**2
           else if (surface_zero_but_normal_nonzero(g_para, g_perp) &
                    .and. (.not. use_simplified) ) then
               !
               kernel = 2.0_q * tpi * felect / g_perp**2 &
                      * (  1 - cos(g_perp * rc) - g_perp * rc &
                      * sin(g_perp * rc) )
           else if (surface_nonzero(g_para, g_perp) .and. &
                    (.not. use_simplified) ) then
               !
               kernel = 2.0_q * tpi * felect / g_norm**2 &
                      * ( 1 + exp(-g_para * rc) &
                        * ( g_perp / g_para * sin(g_perp * Rc) &
                          - cos(g_perp * rc)) )
           else if (use_simplified) then
               !
               kernel = 2._q * tpi * felect / g_norm**2 &
                      * ( 1 - exp(-g_para * rc) * cos(g_perp * rc) )
#ifdef debug_coulomb_cutoff
           else
               call vtutor%error("misplaced g-vector in the &
                  &slab_truncated_kernel subroutine. This subroutine &
                  &uses the absolute value of the g-vector to &
                  &determine the conditionals.")
#endif
           end if
       end associate
     end subroutine

     function calculate_g0_component(is_g0) result(calculate)
!$acc routine seq
       logical, intent(in), optional       :: is_g0
       logical :: calculate

       if (present(is_g0)) then
           calculate = is_g0
       else
           calculate = .false.
       end if
     end function

     subroutine kernel_for_gvector(coulomb_pot, &
                                   g_vector, kernel, is_g0)
!$acc routine seq
       type(coulomb_potential), intent(in) :: coulomb_pot
       real(q), intent(in)                 :: g_vector(3)
       real(q), intent(out)                :: kernel
       logical, intent(in), optional       :: is_g0
       logical                             :: calculate_g0
!$acc routine(calculate_g0_component) seq
!$acc routine(spherical_truncated_kernel) seq
!$acc routine(slab_truncated_kernel) seq
!$acc routine(periodic_kernel) seq
       associate(kernel_truncate => coulomb_pot%kernel_truncate)
       !
       calculate_g0 = calculate_g0_component(is_g0)
       !
       if (kernel_truncate%active) then
           !
           if (kernel_truncate%dimensionality == 0) then
           !
           call spherical_truncated_kernel(kernel_truncate, g_vector, &
                                           kernel, calculate_g0)
           else if (kernel_truncate%dimensionality == 2) then
           !
           call slab_truncated_kernel(kernel_truncate, g_vector, &
                                      kernel, calculate_g0)
           end if
       else
           !
           call periodic_kernel(g_vector, kernel, calculate_g0)
       end if
       !
       end associate
     end subroutine

     subroutine pseudo_kernel_for_gvector(p, t_info, coulomb_pot, &
                                     g_vector, idx_type, &
                                     kernel, is_g0)
!$acc routine seq
       type(type_info), intent(in)         :: t_info
       type(potcar), intent(in)            :: p(t_info%ntyp)
       type(coulomb_potential), intent(in) :: coulomb_pot
       real(q), intent(in)                 :: g_vector(3)
       integer, intent(in)                 :: idx_type
       real(q), intent(out)                :: kernel
       logical, intent(in), optional       :: is_g0
       logical                             :: calculate_g0
       integer                             :: i
       real(q)                             :: argsc
       real(q)                             :: g_norm, rem

!$acc routine(calculate_g0_component) seq
       associate(kernel_truncate => coulomb_pot%kernel_truncate)
       !
       calculate_g0 = calculate_g0_component(is_g0)
       !
       argsc = npspts / p(idx_type)%psgmax
       !
       g_norm = sqrt(  g_vector(1)**2._q &
                     + g_vector(2)**2._q &
                     + g_vector(3)**2._q )
       !
       if (calculate_g0) then
          !
          if (kernel_truncate%active) then
              !
              kernel = p(idx_type)%pscore
          else
              !
              kernel = 0._q
          end if
       else
!
! convert the magnitude of the reciprocal lattice vector to a position
! in the pseudopotential arrays and interpolate the pseudopotential and
! its derivative
!
          i  = int(g_norm * argsc)+1
          rem = g_norm - p(idx_type)%psp(i,1)
          !
          kernel = ( p(idx_type)%psp(i,2) &
                 +   rem * (  p(idx_type)%psp(i,3) &
                            + rem * ( p(idx_type)%psp(i,4) &
                                     + rem * p(idx_type)%psp(i,5) )))
       end if
       !
       end associate
     end subroutine

     subroutine pseudo_kernel_gradient_for_gvector(p, t_info, &
                                        coulomb_pot, &
                                        g_vector, idx_type, &
                                        kernel_gradient, is_g0)
!$acc routine seq
       type(type_info), intent(in)         :: t_info
       type(potcar), intent(in)            :: p(t_info%ntyp)
       type(coulomb_potential), intent(in) :: coulomb_pot
       real(q), intent(in)                 :: g_vector(3)
       real(q), intent(out)                :: kernel_gradient
       integer, intent(in)                 :: idx_type
       logical, intent(in), optional       :: is_g0
       logical                             :: calculate_g0
       real(q)                             :: g_norm, rem
       integer                             :: i
       real(q)                             :: argsc

!$acc routine(calculate_g0_component) seq
       !
       calculate_g0 = calculate_g0_component(is_g0)
       !
       argsc = npspts / p(idx_type)%psgmax
       !
       g_norm = sqrt(  g_vector(1)**2._q &
                     + g_vector(2)**2._q &
                     + g_vector(3)**2._q )
       !
       if (calculate_g0) then
          !
          kernel_gradient = 0._q
       else
          !
          i  = int(g_norm * argsc) + 1
          rem = g_norm - p(idx_type)%psp(i,1)
          !
          kernel_gradient = p(idx_type)%psp(i,3) &
                          + rem * (  2._q * p(idx_type)%psp(i,4) &
                                   + 3._q * rem* p(idx_type)%psp(i,5))
       end if
     end subroutine

     subroutine kernel_gradient_for_gvector(coulomb_pot, &
                                            g_vector, kernel_gradient, &
                                            is_g0)
!$acc routine seq
       type(coulomb_potential), intent(in) :: coulomb_pot
       real(q), intent(in)                 :: g_vector(3)
       real(q), intent(out)                :: kernel_gradient
       logical, intent(in), optional       :: is_g0
       logical                             :: calculate_g0

!$acc routine(calculate_g0_component) seq
!$acc routine(periodic_kernel_gradient) seq
       !
       calculate_g0 = calculate_g0_component(is_g0)
       !
       call periodic_kernel_gradient(g_vector, kernel_gradient, &
                                     calculate_g0)
     end subroutine

     subroutine phase_factor_for_index(idx, ion_positions, &
                                            virtual_crystal_approx, &
                                            phase_factor)
!$acc routine seq
        integer, intent(in)     :: idx(3)
        real(q), intent(in)     :: ion_positions(3)
        real(q), intent(in)     :: virtual_crystal_approx
        complex(q), intent(out) :: phase_factor

        phase_factor = exp(-citpi * (  idx(1) * ion_positions(1) &
                                     + idx(2) * ion_positions(2) &
                                     + idx(3) * ion_positions(3) )) &
                     * virtual_crystal_approx
     end subroutine

     subroutine change_grid_for_kernel_truncation(grid, kernel_truncate)
         type(grid_3d), intent(inout)        :: grid
         type(kernel_truncation), intent(in) :: kernel_truncate
         integer, allocatable :: grid_dim(:)
         integer, allocatable :: remainder(:)
         integer :: idx
         logical :: correct

         allocate(grid_dim, source=grid%ngptar)
         allocate(remainder, mold=grid%ngptar)
         !
         associate(dimensionality => kernel_truncate%dimensionality, &
                   padding_factor => kernel_truncate%padding_factor, &
                   surface_normal_direction => &
                              kernel_truncate%surface_normal_direction)
         !
         correct = .false.
         remainder = 0
         do while (.not. correct)
            !
            call FFTCHK(grid_dim)
            if (dimensionality == 0) then
                !
                remainder(1) = modulo(grid_dim(1), padding_factor)
                remainder(2) = modulo(grid_dim(2), padding_factor)
                remainder(3) = modulo(grid_dim(3), padding_factor)
            else if (dimensionality == 2) then
                !
                if (surface_normal_direction == 1) then
                    remainder(1) = modulo(grid_dim(1), padding_factor)
                else if (surface_normal_direction == 2) then
                    remainder(2) = modulo(grid_dim(2), padding_factor)
                else if (surface_normal_direction == 3) then
                    remainder(3) = modulo(grid_dim(3), padding_factor)
                end if
            end if
            !
            if (all(remainder == 0)) then
               correct = .true.
            end if
            !
            if (.not. correct) then
                if (dimensionality == 0) then
                    !
                    grid_dim(1) = grid_dim(1) + 1
                    grid_dim(2) = grid_dim(2) + 1
                    grid_dim(3) = grid_dim(3) + 1
                else if (dimensionality == 2) then
                    !
                    if (surface_normal_direction == 1) then
                        grid_dim(1) = grid_dim(1) + 1
                    else if (surface_normal_direction == 2) then
                        grid_dim(2) = grid_dim(2) + 1
                    else if (surface_normal_direction == 3) then
                        grid_dim(3) = grid_dim(3) + 1
                    end if
                end if
            end if
         end do
         end associate
         !
         grid%ngptar(1) = grid_dim(1)
         grid%ngptar(2) = grid_dim(2)
         grid%ngptar(3) = grid_dim(3)
     end subroutine

     subroutine create_grids_for_kernel_truncation(gridc, &
                                                   kernel_truncate, &
                                                   comm_kin, &
                                                   coarse_grid, &
                                                   fine_grid, &
                                                   translation_table)
         type(grid_3d), intent(in), target      :: gridc
         type(kernel_truncation), intent(inout) :: kernel_truncate
         type(grid_3d), target :: coarse_grid, fine_grid
         type(transit), target :: translation_table
         type(communic), target :: comm_kin
         integer :: padding_factor
         integer :: fine_ngx, fine_ngy, fine_ngz
         integer :: coarse_ngx, coarse_ngy, coarse_ngz

         associate(coarsen_before_pad => &
                              kernel_truncate%coarsen_before_pad, &
                   dimensionality => kernel_truncate%dimensionality, &
                   padding_factor => kernel_truncate%padding_factor, &
                   surface_normal_direction => &
                              kernel_truncate%surface_normal_direction,&
                   ngx => gridc%NGX, ngy => gridc%NGY, ngz => gridc%NGZ)
         !
         if (coarsen_before_pad) then
             !
             kernel_truncate%fine_grid => gridc
             fine_grid = gridc
             !
             if (dimensionality == 0) then
                 !
                 coarse_ngx = ngx / padding_factor
                 coarse_ngy = ngy / padding_factor
                 coarse_ngz = ngz / padding_factor
             elseif (dimensionality == 2) then
                 !
                 if (surface_normal_direction == 1) then
                    coarse_ngx = ngx / padding_factor
                    coarse_ngy = ngy
                    coarse_ngz = ngz
                 elseif (surface_normal_direction == 2) then
                    coarse_ngx = ngx
                    coarse_ngy = ngy / padding_factor
                    coarse_ngz = ngz
                 elseif (surface_normal_direction == 3) then
                    coarse_ngx = ngx
                    coarse_ngy = ngy
                    coarse_ngz = ngz / padding_factor
                end if
            end if
            !
            call inilgrd(coarse_ngx, coarse_ngy, coarse_ngz, &
                         coarse_grid)
            coarse_grid%comm => COMM_KIN
            call gen_rc_sub_grid(coarse_grid, gridc, &
                                 translation_table, &
                                 .true., .true.)
            CALLMPI( mapset(coarse_grid))
            call fftgridplan(coarse_grid)
            !
            kernel_truncate%coarse_grid => coarse_grid
            kernel_truncate%translation_table => translation_table
         !
         else
             !
             kernel_truncate%coarse_grid => gridc
             coarse_grid = gridc
             !
             if (dimensionality == 0) then
                 !
                 fine_ngx = ngx * padding_factor
                 fine_ngy = ngy * padding_factor
                 fine_ngz = ngz * padding_factor
             elseif (dimensionality == 2) then
                 !
                 if (surface_normal_direction == 1) then
                     fine_ngx = ngx * padding_factor
                     fine_ngy = ngy
                     fine_ngz = ngz
                 elseif (surface_normal_direction == 2) then
                     fine_ngx = ngx
                     fine_ngy = ngy * padding_factor
                     fine_ngz = ngz
                 elseif (surface_normal_direction == 3) then
                     fine_ngx = ngx
                     fine_ngy = ngy
                     fine_ngz = ngz * padding_factor
                 end if
             end if
             !
             call inilgrd(fine_ngx, fine_ngy, fine_ngz, fine_grid)
             fine_grid%comm => comm_kin
             fine_grid%comm_kin => comm_kin
             call gen_rc_grid(fine_grid)
             CALLMPI( mapset(fine_grid))
             call fftgridplan(fine_grid)
             !
             kernel_truncate%fine_grid => fine_grid
         end if
     end associate
     end subroutine

     subroutine set_truncation_length(latt_cur, kernel_truncate)
         type(latt), intent(in) :: latt_cur
         type(kernel_truncation), intent(inout) :: kernel_truncate

         associate(truncation_factor => &
                         kernel_truncate%truncation_factor, &
                   surface_normal_direction => &
                         kernel_truncate%surface_normal_direction, &
                   dimensionality => kernel_truncate%dimensionality)

         if (dimensionality == 0) then
             !
             kernel_truncate%truncation_length = truncation_factor &
                                               * maxval(latt_cur%ANORM)
         else if (dimensionality == 2) then
             !
             kernel_truncate%truncation_length = truncation_factor &
                       * latt_cur%ANORM(surface_normal_direction)
         end if
         end associate
     end subroutine

     subroutine setup_padding_and_truncation_factor(latt_cur, kernel_truncate)
         type(latt), intent(in) :: latt_cur
         type(kernel_truncation), intent(inout) :: kernel_truncate
         logical :: set_padding, set_truncation

         associate(is_truncated => kernel_truncate%active, &
                   dimensionality => kernel_truncate%dimensionality, &
                   coarsen_before_pad => &
                             kernel_truncate%coarsen_before_pad, &
                   use_simplified => &
                              kernel_truncate%use_simplified_2d_kernel)
         !
         if (kernel_truncate%padding_factor == huge(1)) then
             !
             set_padding = .true.
         else
             !
             set_padding = .false.
         end if
         !
         if (kernel_truncate%truncation_factor == huge(1._q)) then
             !
             set_truncation = .true.
         else
             !
             set_truncation = .false.
         end if

         !
         if (is_truncated .and. dimensionality == 0) then
             !
             if (coarsen_before_pad) then
                !
                if (set_padding) kernel_truncate%padding_factor = 3
                if (set_truncation) kernel_truncate%truncation_factor = sqrt(3.0_q)
             else
                !
                if (set_padding) kernel_truncate%padding_factor = 3
                if (set_truncation) kernel_truncate%truncation_factor = sqrt(3.0_q)
             end if
         elseif (is_truncated .and. dimensionality == 2) then
            !
            if (use_simplified) then
                !
                if ((.not. set_padding) .or. (.not. set_truncation)) then
                    !
                    call vtutor%error("If using a simplified 2D kernel &
                        &INCAR cannot set the padding and truncation.")
                end if
                !
                kernel_truncate%padding_factor = 2
                kernel_truncate%truncation_factor = 1._q
            else
                !
                if (set_padding) kernel_truncate%padding_factor = 2
                if (set_truncation) kernel_truncate%truncation_factor = 1._q
            end if
         else
            call vtutor%error("Currently only a dimensionality of 0/2 &
                           &are suppored.")
         end if
         !
         call set_truncation_length(latt_cur, kernel_truncate)
         !
         end associate
     end subroutine

     subroutine validate_kernel_dimensionality(kernel_truncate)
         type(kernel_truncation), intent(in) :: kernel_truncate
         associate(dimensionality => kernel_truncate%dimensionality)
         if ((dimensionality < 0) .or. (dimensionality > 3)) then
             !
             call vtutor%error("Please choose a dimensionality value&
                        & of between 0 (atomic or molecular system), &
                        & 1 (nanorods or wires), 2 (surfaces), &
                        & 3 (bulk)")
         endif
         end associate
     end subroutine

     subroutine validate_surface_normal_direction(kernel_truncate)
         type(kernel_truncation), intent(in) :: kernel_truncate

         associate(dimensionality => kernel_truncate%dimensionality, &
                   surface_normal_direction => &
                            kernel_truncate%surface_normal_direction)
         !
         if (dimensionality /= 2 .and. &
             surface_normal_direction /= -1) then
             !
             call vtutor%error("Set the surface normal direction only &
                 &when you have a 2D or quasi-2D system.")
         end if
         if (dimensionality == 2 .and. &
             ( surface_normal_direction > 3 .or. &
               surface_normal_direction < 1 )) then
             !
             call vtutor%error("Since you are using a two dimensional &
                 &electrostatic correction, make sure that the &
                 &ISURFACE is set to between 1-3, where 1-> x, &
                 &2-> y and 3-> z")
         end if
         end associate
     end subroutine

     subroutine generate_padded_lattice(kernel_truncate, latt_cur, &
                                        latt_cur_padded)
        type(kernel_truncation), intent(in) :: kernel_truncate
        type(latt), intent(in)              :: latt_cur
        type(latt), intent(inout)           :: latt_cur_padded

        associate(dimensionality => kernel_truncate%dimensionality, &
                  padding_factor => kernel_truncate%padding_factor, &
                  surface_normal_direction => &
                  kernel_truncate%surface_normal_direction)
        latt_cur_padded%scale = latt_cur%scale
        if (kernel_truncate%active) then
            !
            if (dimensionality == 0) then
                !
                latt_cur_padded%A = latt_cur%A * padding_factor
            else if (dimensionality == 2) then
                !
                latt_cur_padded%A = latt_cur%A
                latt_cur_padded%A(surface_normal_direction,:) = &
                        latt_cur%A(surface_normal_direction,:) &
                        * padding_factor
            end if
        else
            !
            latt_cur_padded%A = latt_cur%A
        end if
        !
        call lattic(latt_cur_padded)
        end associate
     end subroutine


!> This routine generates the padded version of t_info for the
!> coulomb cutoff routines. posion are the positions of the ions
!> in direct coordinates. Keep this variable in scope if you
!> want t_info_padded to work as expected. This requirement
!> stems from the fact that posion in t_info is only a pointer
!> to the positions
     subroutine generate_padded_t_info(kernel_truncate, t_info, &
                                       t_info_padded, posion)
        type(kernel_truncation), intent(in)        :: kernel_truncate
        type(type_info), intent(in)                :: t_info
        type(type_info), intent(out)               :: t_info_padded
        real(q), allocatable, target, intent(out)  :: posion(:,:)

        associate(dimensionality => kernel_truncate%dimensionality, &
                  padding_factor => kernel_truncate%padding_factor, &
                  surface_normal_direction => &
                  kernel_truncate%surface_normal_direction)
        !
        t_info_padded = t_info
        !
        if (kernel_truncate%active) then
            !
            allocate(posion(3,t_info%nions))
            !
            posion = t_info%posion
            !
            if (dimensionality == 0) then
                !
                posion = posion / padding_factor
            else if (dimensionality == 2) then
                !
                posion(surface_normal_direction,:) = &
                    posion(surface_normal_direction,:) / padding_factor
            end if
            !
            t_info_padded%posion => posion
        end if
        end associate
     end subroutine

     subroutine pad3d(kernel_truncate, quantity, padded_quantity, &
                      direction, real_space_rescale)
        type(kernel_truncation), intent(in) :: kernel_truncate
        complex(q), intent(inout)  :: quantity(:)
        complex(q), intent(inout)  :: padded_quantity(:)
        integer, intent(in)        :: direction
        real(q), optional          :: real_space_rescale
        complex(q), allocatable    :: quantity_coarse(:)
        complex(q), allocatable    :: quantity_fine(:)
        RGRID, allocatable         :: real_quantity_coarse(:)
        RGRID, allocatable         :: real_quantity_fine(:)
        integer :: nc, n1, n2, n3, ni, ncp, pad_ni

        !$acc enter data copyin(direction) __IF_ASYNC__
        !
        associate(coarse_grid => kernel_truncate%coarse_grid, &
                  fine_grid => kernel_truncate%fine_grid, &
                  translation_table => &
                  kernel_truncate%translation_table, &
                  coarsen_before_pad => &
                  kernel_truncate%coarsen_before_pad)
        !
        if (direction == 1) then
            allocate(quantity_coarse(coarse_grid%MPLWV))
            quantity_coarse = 0._q
            !$acc enter data copyin(quantity_coarse) __IF_ASYNC__
        else if (direction == -1) then
            allocate(quantity_fine(fine_grid%MPLWV))
            quantity_fine = 0._q
            !$acc enter data copyin(quantity_fine) __IF_ASYNC__
        end if
        !
        allocate(real_quantity_coarse(coarse_grid%RL%NP))
        allocate(real_quantity_fine(fine_grid%RL%NP))
        !
        real_quantity_coarse = 0._q
        real_quantity_fine = 0._q
        !$acc enter data copyin(real_quantity_coarse) __IF_ASYNC__
        !$acc enter data copyin(real_quantity_fine) __IF_ASYNC__
        if (direction == 1) then
            call rc_add(quantity, 1.0_q, &
                        quantity, 0.0_q, &
                        quantity_coarse, coarse_grid)
        else if (direction == -1) then
            call rc_add(padded_quantity, 1._q, &
                        padded_quantity, 0._q, &
                        quantity_fine, fine_grid)
        end if
        !
        if (direction == 1) then
            call FFT3D(quantity_coarse, coarse_grid, 1)
            call rl_add(quantity_coarse, 1._q, &
                        quantity_coarse, 0._q, &
                        real_quantity_coarse, coarse_grid)
        else if (direction == -1) then
            !
            call FFT3D(quantity_fine, fine_grid, 1)
            call rl_add(quantity_fine, 1._q, &
                        quantity_fine, 0._q, &
                        real_quantity_fine, fine_grid)
        end if
        !$acc parallel loop collapse(2) private(n2,n3,ncp,ni,pad_ni) &
        !$acc present(kernel_truncate,real_quantity_fine,real_quantity_coarse) __IF_ASYNC__
        do nc = 1, coarse_grid%RL%NCOL
__CPU   n2 = coarse_grid%RL%I2(nc)
__CPU   n3 = coarse_grid%RL%I3(nc)
__CPU   ncp = fine_grid%RL%INDEX(n2-1, n3-1)
        do n1 = 1, coarse_grid%RL%NROW
__GPU      n2 = coarse_grid%RL%I2(nc)
__GPU      n3 = coarse_grid%RL%I3(nc)
__GPU      ncp = fine_grid%RL%INDEX(n2-1, n3-1)
           !
           ni = (nc-1) * coarse_grid%RL%NROW + n1
           pad_ni = (ncp-1) * fine_grid%RL%NROW + n1
           !
           if (direction == 1) then
               !
               real_quantity_fine(pad_ni) = real_quantity_coarse(ni)
           else if (direction == -1) then
               !
               real_quantity_coarse(ni) = real_quantity_fine(pad_ni)
           end if
           end do
        end do
        !
        if (direction == 1) then
            !
            if (present(real_space_rescale)) then
                !$acc enter data copyin(real_space_rescale)
                call rl_add(real_quantity_fine, real_space_rescale, &
                            real_quantity_fine, 0._q, &
                            real_quantity_fine, fine_grid)
                !$acc exit data delete(real_space_rescale)
            end if
            call fft_rc_scale(real_quantity_fine, padded_quantity, &
                              fine_grid)
        else if (direction == -1) then
            !
            if (present(real_space_rescale)) then
                !$acc enter data copyin(real_space_rescale)
                call rl_add(real_quantity_coarse, real_space_rescale, &
                            real_quantity_coarse, 0._q, &
                            real_quantity_coarse, coarse_grid)
                !$acc exit data delete(real_space_rescale)
            end if
            call fft_rc_scale(real_quantity_coarse, quantity, coarse_grid)
        end if
        !
        if (direction == 1) then
            !$acc exit data delete(quantity_coarse) __IF_ASYNC__
        else if (direction == -1) then
            !$acc exit data delete(quantity_fine) __IF_ASYNC__
        end if
        !$acc exit data delete(real_quantity_coarse) __IF_ASYNC__
        !$acc exit data delete(real_quantity_fine) __IF_ASYNC__
        !$acc exit data delete(direction) __IF_ASYNC__
        end associate
     end subroutine

     subroutine write_coulomb_cutoff(unit_idx, kernel_truncate)
         integer, intent(in)                 :: unit_idx
         type(kernel_truncation), intent(in) :: kernel_truncate

         if (unit_idx > 0) then
           write(unit_idx, 1406) kernel_truncate%active, &
                                 kernel_truncate%dimensionality, &
                                 kernel_truncate%surface_normal_direction, &
                                 kernel_truncate%coarsen_before_pad
         end if
         !
1406 FORMAT(&
        'Coulomb truncation method'/ &
        '  LTRUNCATE          = ', L6, '  Coulomb kernel truncation to be used' / &
        '  IDIMENSIONALITY    = ', I6, '  Dimensionality of cutoff (0 molecules / 2 surfaces)' / &
        '  ISURFACE           = ', I6, '  Direction of surface normal'/ &
        '  LCOARSEN           = ', L6, '  Coarsen the charge grid before padding, preferred'/ &
    )
     end subroutine

#ifdef debug_coulomb_cutoff
! ------------- DEBUG ROUTINES -------------

     subroutine debug_writer_potential(pgridc, t_info, latt_cur, &
                                       filename, potential, description)
         use fileio, only: outpot
         use poscar, only: outpos
         type(grid_3d), intent(in) :: pgridc
         character(len=*), intent(in) :: filename
         complex(q), intent(in) :: potential(:)
         type (type_info), intent(in) :: t_info
         type(latt), intent(in) :: latt_cur
         character(len=*) :: description
         character(len=40) :: charac
         complex(q), allocatable :: real_potential(:)
! work
         integer :: iupad
         integer :: ierr
         integer :: node_me, ionode
         node_me = pgridc%comm%node_me
         ionode  = pgridc%comm%ionode
         allocate(real_potential(pgridc%MPLWV))
         call rc_add(potential, 1._q, &
                     potential, 0._q, &
                     real_potential, pgridc)
         call FFT3D(real_potential, pgridc, 1)
         charac = " "
         !
         io_begin
         open(newunit=iupad, file=filename, iostat=ierr, status="replace")
         call outpos(iupad, .false., charac, t_info, &
                     latt_cur%scale, &
                     latt_cur%A, &
                     .false., t_info%posion)
         io_end
         call outpot(pgridc, iupad, .true., real_potential)
         io_begin
         close(iupad)
         io_end
         call sum_rl(description // " (real space sum):", &
                     potential, pgridc)
     end subroutine debug_writer_potential

     subroutine debug_write_charge(pgridc, t_info, latt_cur, &
                                   filename, charge_density, &
                                   description)
         use fileio, only: outchg
         use poscar, only: outpos
         type(grid_3d), intent(in) :: pgridc
         character(len=*), intent(in) :: filename
         complex(q), intent(in) :: charge_density(:)
         type (type_info), intent(in) :: t_info
         type(latt), intent(in) :: latt_cur
         character(len=*) :: description
         character(len=40) :: charac
! work
         integer :: iupad
         integer :: ierr
         integer :: node_me, ionode
         node_me = pgridc%comm%node_me
         ionode  = pgridc%comm%ionode
         charac = " "
         io_begin
         open(newunit=iupad, file=filename, iostat=ierr, status="replace")
         call outpos(iupad, .false., charac, t_info, &
                     latt_cur%scale, &
                     latt_cur%A, &
                     .false., t_info%posion)
         io_end
         call outchg(pgridc, iupad, .true., charge_density)
         io_begin
         close(iupad)
         io_end
         if (node_me==1) then
            print *, description,  " (g=0 component):", charge_density(1)
         end if
     end subroutine debug_write_charge

! positive for open, negative for close
     subroutine debug_toggle_file(filename, iupad, open_or_close)
         character(len=*), intent(in) :: filename
         integer, intent(inout) :: iupad
         integer, intent(in) :: open_or_close
         integer :: ierr

         if (open_or_close > 0) then
             open(newunit=iupad, file=filename, iostat=ierr, status="replace")
         else
             close(iupad)
         end if
     end subroutine debug_toggle_file

#endif

end module coulomb_cutoff
