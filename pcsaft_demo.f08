!***********************************************************************************************************************************
!
!                                       -- PC-SAFT FORTRAN 2008 implementation demo program --
!
!***********************************************************************************************************************************
! Daniel Cordier, Institut UTINAM, Besançon, France (until september, 2015).
!             GSMA, Reims, France (from september, 2015).
!
!             - daniel.cordier@univ-reims.fr
!             - Research Gate ...: https://www.researchgate.net/profile/Daniel_Cordier
!             - OrCID ...........: http://orcid.org/0000-0003-4515-6271
!
program pcsaft_demo
  use utils_DC
  
  use mod_pcsaft
  
  implicit none
  
  character(len=*), parameter :: subname = 'Main Program: "pcsaft_demo"'
  
  character(len=1) :: state
  integer :: nc ! Number of species in a given mixture.
  
  integer :: k
  character(len=comp_name_ncmax), allocatable, dimension(:) :: compound !
  real(dp), allocatable, dimension(:) :: x, x0  ! mole fractions of the various species taken into account.
  real(dp)                  :: P      ! the pressure (in Pa).
  real(dp)                  :: T      ! the temperature (in K).
  real(dp)                  :: rho, rho_MKSA    ! the density of the mixture (kg.m^-3).
  real(dp), allocatable, dimension(:) :: muskT_k
  real(dp), allocatable, dimension(:) :: ln_phi,  ln_phi0! the activity coefficients in neperian log (no unit).
  
  real(dp), allocatable, dimension(:) :: x_weigth ! mass fractions of the various species taken into account.
  real(dp), allocatable, dimension(:) :: M_mol    ! molar mass
  real(dp)                  :: M_bar, Tm, Delta_Hm, f
  real(dp)                  :: alpha_V
  
  integer  :: i, Npt, ok
  real(dp) :: Pmin, Pmax, Delta_P, P_rhoT, T_Prho, sum, Vmol, Nmolecules, a, Navogadro
  
  logical :: underflow_support, gradual
  character(len=100) :: temperature

  ! ===============================================================================================================
  ! Pour la tabulation de la densité des mélanges N2/CO2 (cf. collaboration avec Sébastien Lebonnois) :
  real(dp) :: xN2, xN2_min, xN2_max, Delta_xN2 
  ! ===============================================================================================================
  
  !----------------------------------------------------------------------------------------------------------------
  write(6,*) ''
  call wf(' --------------------------------------------------------------------', 'bright red')
  call wf(' -                     PC_SAFT demo program                         -', 'bright red')
  call wf(' --------------------------------------------------------------------', 'bright red')
  call wf(' - D. Cordier, 2013-2021', 'bright red')
  call wf('   https://www.researchgate.net/profile/Daniel_Cordier', 'bright blue')
  call wf('   http://orcid.org/0000-0003-4515-6271', 'bright blue')
  write(6,*) '' 
  
  ! -------------------------------------------------------------------------------------------------------------------------------
  ! Properties of the compilation:  
  call display_compil_numpres (rho_MKSA)
  
  !----------------------------------------------------------------------------------------------------------------
  ! 1- Vapor mixture of N2 and CH4:
  
  write(6,*) ''
  call wf(' --------------------------------------------------------------------', 'bright red')
  call wf(' > Example: a vapor containing N2 and CH4 ---------------------------', 'bright red')
  nc = 2
  allocate(compound(1:nc), x(1:nc), x0(1:nc), muskT_k(1:nc), ln_phi(1:nc),  ln_phi0(1:nc), x_weigth(1:nc), M_mol(1:nc), stat=ok)
  
  state = 'V' ! Physical state.
  compound(1)= 'N2'  ! species 1.
  compound(2)= 'CH4' ! species 2.
  ! Chemical composition in mole fraction:
  x(1)= 0.20_dp      ! N2
  x(2)= 1._dp - x(1) ! CH4
  
  T= 90._dp ! Temperature (K).
  P= 1.5_dp * onebar_in_Pa ! Pressure (Pa).
  
  call pcsaft_PT(subname,1,state,compound,x, P,T,rho,muskT_k,ln_phi) ! PC-SAFT computation.
  
  write(6,'(A18,A1)') '  - State       = ', state
  write(6,'(A18,I3)') '  - N. species  = ', nc
  do i= 1, nc
     write(6,'(A12,I3,A3,A3)') '  - Species ', i ,' = ', compound(i)
  end do
  write(6,'(A18,F7.3,A4)') '  - P           = ', P/onebar_in_Pa, ' bar'
  write(6,'(A18,F7.3,A2)') '  - T           = ', T, ' K'
  write(6,'(A18,F7.3,A8)') '  - rho         = ', rho, ' kg m^-3'
  write(6,*) ''
  
  deallocate(compound, x, x0, muskT_k, ln_phi,  ln_phi0, x_weigth, M_mol)

  !----------------------------------------------------------------------------------------------------------------
  ! 2- Liquid mixture of N2 and CH4:
  write(6,*) ''
  call wf(' --------------------------------------------------------------------', 'bright red')
  call wf(' > Example: a liquid containing N2 and CH4 --------------------------', 'bright red')
  nc = 2
  allocate(compound(1:nc), x(1:nc), x0(1:nc), muskT_k(1:nc), ln_phi(1:nc),  ln_phi0(1:nc), x_weigth(1:nc), M_mol(1:nc), stat=ok)
  
  state = 'L' ! Physical state.
  compound(1)= 'N2'  ! species 1.
  compound(2)= 'CH4' ! species 2.
  ! Chemical composition in mole fraction:
  x(1)= 0.20_dp      ! N2
  x(2)= 1._dp - x(1) ! CH4
  
  T= 90._dp ! Temperature (K).
  P= 1.5_dp * onebar_in_Pa ! Pressure (Pa).
  
  call pcsaft_PT(subname,1,state,compound,x, P,T,rho,muskT_k,ln_phi) ! PC-SAFT computation.
  
  write(6,'(A18,A1)') '  - State       = ', state
  write(6,'(A18,I3)') '  - N. species  = ', nc
  do i= 1, nc
     write(6,'(A12,I3,A3,A3)') '  - Species ', i ,' = ', compound(i)
  end do
  write(6,'(A18,F7.3,A4)') '  - P           = ', P/onebar_in_Pa, ' bar'
  write(6,'(A18,F7.3,A2)') '  - T           = ', T, ' K'
  write(6,'(A18,F7.3,A8)') '  - rho         = ', rho, ' kg m^-3'
  write(6,*) ''
  
  deallocate(compound, x, x0, muskT_k, ln_phi,  ln_phi0, x_weigth, M_mol)
  
  !----------------------------------------------------------------------------------------------------------------
  ! 3- Liquid mixture of N2, CH4 and C2H6:
  write(6,*) ''
  call wf(' --------------------------------------------------------------------', 'bright red')
  call wf(' > Example: a liquid containing N2, CH4 and C2H6 --------------------', 'bright red')
  nc = 3
  allocate(compound(1:nc), x(1:nc), x0(1:nc), muskT_k(1:nc), ln_phi(1:nc),  ln_phi0(1:nc), x_weigth(1:nc), M_mol(1:nc), stat=ok)
  
  state = 'L' ! Physical state.
  compound(1)= 'N2'   ! species 1.
  compound(2)= 'CH4'  ! species 2.
  compound(3)= 'C2H6' ! species 3.
  ! Chemical composition in mole fraction:
  x(1)= 0.20_dp             ! N2
  x(2)= 0.20_dp             ! CH4
  x(3)= 1._dp - x(1) - x(2) ! C2H6
  
  T= 90._dp ! Temperature (K).
  P= 1.5_dp * onebar_in_Pa ! Pressure (Pa).
  
  call pcsaft_PT(subname,1,state,compound,x, P,T,rho,muskT_k,ln_phi) ! PC-SAFT computation.
  
  write(6,'(A18,A1)') '  - State       = ', state
  write(6,'(A18,I3)') '  - N. species  = ', nc
  do i= 1, nc
     write(6,'(A12,I3,A3,A3)') '  - Species ', i ,' = ', compound(i)
  end do
  write(6,'(A18,F7.3,A4)') '  - P           = ', P/onebar_in_Pa, ' bar'
  write(6,'(A18,F7.3,A2)') '  - T           = ', T, ' K'
  write(6,'(A18,F7.3,A8)') '  - rho         = ', rho, ' kg m^-3'
  write(6,*) '' 
  
  !deallocate(compound, x, x0, muskT_k, ln_phi,  ln_phi0, x_weigth, M_mol)

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Appel test de PC-SAFT avec comme variables indépendantes P et rho :
  call pcsaft_Prho(subname,1,state,compound,x,P,rho,T_Prho,muskT_k,ln_phi)
  
  write(6,*) ' > T= ', T, ' -- T_Prho= ', T_Prho
  
  write(6,*) ''
  call wf(' > Computation done!', 'bright red')
  write(6,*) ''
  
end program pcsaft_demo
