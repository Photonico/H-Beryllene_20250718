#include "symbol.inc"
!************************************************************************
! RCS:  $Id: constant.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
!> Physical constants and unit conversion used in VASP
!>
!> VASP uses mostly atomic units (a.u.) except:
!> + length - in Angstroem instead of Bohr
!> + energy - in eV instead of Hartree (2*Rydberg)
!>
!> Symbols:
!> + \f( e   \f) - electron charge
!> + \f( m_e \f) - electron mass
!> + \f( c   \f) - speed of light
!> + \f( \epsilon_0 \f) - permitivity of free space
!
!***********************************************************************
MODULE constant
      USE prec

      ! Some important Parameters, to convert to a.u.
      REAL(q), PARAMETER :: AUTOA=0.529177249_q !< 1 a.u. of length (Bohr) in Angstroem
      REAL(q), PARAMETER :: RYTOEV=13.605826_q !< 1 Ry in Ev
      REAL(q), PARAMETER :: CLIGHT=137.037  !< speed of light in a.u.
      REAL(q), PARAMETER :: EVTOJ=1.60217733E-19_q !< 1 eV in Joule
      REAL(q), PARAMETER :: AMTOKG=1.6605402E-27_q !< 1 atomic mass unit ("proton mass") in kg
      REAL(q), PARAMETER :: BOLKEV=8.6173857E-5_q !< Boltzmanns constant in eV/K
      REAL(q), PARAMETER :: BOLK=BOLKEV*EVTOJ !< Boltzmanns constant in Joule/K

      REAL(q), PARAMETER :: EVTOKCAL=23.06

      REAL(q),   PARAMETER  :: PI =3.141592653589793238_q !< Pi
      REAL(q),   PARAMETER  :: TPI=2*PI !< 2 Pi
      COMPLEX(q),PARAMETER  :: CITPI = (0._q,1._q)*TPI !< 2 i Pi

!> electron charge divided by 4*pi times the permittivity of free space
!>         in atomic units this is just e^2
!>
!> \f( \frac{e^2}{4 \pi \epsilon_0} \f)
      REAL(q), PARAMETER :: FELECT = 2*RYTOEV*AUTOA
!> electron charge divided by the permittivity of free space
!>         in atomic units this is just 4 pi e^2
!>
!> \f( \frac{e^2}{\epsilon_0} \f)
      REAL(q), PARAMETER :: EDEPS = 4*PI*FELECT
!> reduced planck constant squared divided by 2 times the electron mass
!>
!> \f( \frac{\hbar^2}{2 m_e} \f)
      REAL(q), PARAMETER :: HSQDTM = RYTOEV*AUTOA*AUTOA

!> vector field A times momentum times e/ (2 m_e c) is an energy
!> magnetic moments are supplied in Bohr magnetons
!>~~~
!> e / (2 m_e c) A(r) p(r) = energy
!> e / (2 m_e c) m_s x ( r - r_s) / (r-r_s)^3 hbar nabla =
!> e^2 hbar^2 / (2 m_e^2 c^2) 1/ lenght^3 = energy
!>~~~
!> conversion factor from magnetic moment to energy
!> checked independently in SI by Gilles de Wijs
      REAL(q), PARAMETER :: MAGMOMTOENERGY=1/CLIGHT**2*AUTOA**3*RYTOEV

!> dimensionless number connecting input and output magnetic moments
!> AUTOA e^2 (2 m_e c^2)
      REAL(q), PARAMETER :: MOMTOMOM=AUTOA/CLIGHT/CLIGHT/2

      REAL(q), PARAMETER :: AUTOA2=AUTOA *AUTOA
      REAL(q), PARAMETER :: AUTOA3=AUTOA2*AUTOA
      REAL(q), PARAMETER :: AUTOA4=AUTOA2*AUTOA2
      REAL(q), PARAMETER :: AUTOA5=AUTOA3*AUTOA2

      CONTAINS

      SUBROUTINE NUP_CONST()
      END SUBROUTINE

END MODULE
