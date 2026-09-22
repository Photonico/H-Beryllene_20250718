#include "symbol.inc"

module coulomb_cutoff_gradients
GPACC use moffload_struct_def
GPACC use moffload
      use prec, only: q
      use mgrid, only: grid_3d, inilgrd, gen_rc_sub_grid, gen_rc_grid
      use lattice, only: latt, lattic, dirkar
      use pot_electrostat, only: potion
      use pot_struct_def
      use pseudo_struct_def, only: potcar
      use poscar_struct_def, only: type_info
      use mgrid_struct_def, only: transit
      use base, only: info_struct
      use constant, only: EDEPS, PI, TPI, FELECT, CITPI

      implicit none

contains
    subroutine truncate_forloc(gridc, p, t_info, latt_cur, &
                               coulomb_pot, charge_density, &
                               structure_factor, forces)
       type(grid_3d), intent(in) :: gridc
       type(latt), intent(in) :: latt_cur
       type(type_info), intent(in) :: t_info
       type(potcar), intent(in) :: p(t_info%ntyp)
       type(coulomb_potential), intent(in) :: coulomb_pot
       complex(q) :: charge_density(:)
       complex(q) :: structure_factor(:,:)
       real(q), intent(out) :: forces(:,:)
       integer :: i, nis, nt, n1, n2, n3, nc, ni, nion, ni_add
       complex(q) :: for
       real(q) :: for1, for2, for3
       real(q) :: dummy, factm
       complex(q) :: charge
       complex(q), allocatable :: ion_pot(:)

       PROFILING_START('truncate_forloc')
       !
       allocate(ion_pot(gridc%MPLWV))
       !$acc enter data copyin(ion_pot) __IF_ASYNC__
       !
       nis = 1
       !
       do nt = 1, t_info%ntyp
       !
       ni_add = t_info%NITYP(nt)
       !
       do nion = nis, ni_add + nis - 1
           !acc kernels present(ion_pot) __IF_ASYNC__
           ion_pot = 0.0_q
           !acc end kernels
           !
           call potion(gridc, p, latt_cur, t_info, coulomb_pot, &
                       ion_pot, structure_factor, dummy, [nion])
           !
           for1 = 0.0_q; for2 = 0.0_q; for3 = 0.0_q
           !
           !$acc parallel loop reduction(+:for1,for2,for3) &
           !$acc private(n2,n3,ni,factm,for) &
           !$acc present(gridc,charge_density,ion_pot)
           do nc = 1, gridc%RC%NCOL
__CPU      n2 = gridc%RC%I2(nc)
__CPU      n3 = gridc%RC%I3(nc)
           do n1 = 1, gridc%RC%NROW
__GPU          n2 = gridc%RC%I2(nc)
__GPU          n3 = gridc%RC%I3(nc)
               !
               ni = (nc-1) * gridc%RC%NROW + n1
               !
               SETFACT1
               SETFACT
               !
               for = cmplx(0.0_q, -1.0_q) * ion_pot(ni) &
                   * MULFACT conjg(charge_density(ni))
               !
               for1 = for1 - gridc%LPCTX_(n1) * real(for,q)
               for2 = for2 - gridc%LPCTY_(n2) * real(for,q)
               for3 = for3 - gridc%LPCTZ_(n3) * real(for,q)
           end do
           end do
           !
           forces(1,nion) = for1 * TPI
           forces(2,nion) = for2 * TPI
           forces(3,nion) = for3 * TPI
           enddo
           !
           nis = nis + ni_add
           !
       enddo
       !
       CALLMPI( M_sum_d(gridc%comm, forces(1,1), t_info%nions*3))
       !
       call dirkar(t_info%nions, forces, latt_cur%b)
       !$acc exit data delete(ion_pot) __IF_ASYNC__
       !
       PROFILING_STOP('truncate_forloc')
    end subroutine
end module
