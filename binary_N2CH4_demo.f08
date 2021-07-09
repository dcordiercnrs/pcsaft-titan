!***********************************************************************************************************************************
!
!                           -- Demo program of a binary system simulation using PC-SAFT --
!
!***********************************************************************************************************************************
! Daniel Cordier, Institut UTINAM, Besançon, France (until september, 2015).
!                 GSMA, Reims, France (from september, 2015).
!
!             - daniel.cordier@univ-reims.fr
!             - Research Gate ...: https://www.researchgate.net/profile/Daniel_Cordier
!             - OrCID ...........: http://orcid.org/0000-0003-4515-6271
!
program binary_N2CH4_demo
  use utils_dc
  use mod_pcsaft
  
  implicit none
  
  character(len=*), parameter :: subname = 'Main Program: "binary_N2CH4_demo"'
  
  integer, parameter :: nc = 2 ! Number of involved species (here we have a _binary_ mixture).
  
  integer :: k
  character(len=comp_name_ncmax), dimension(1:nc) :: compound !
  real(dp), dimension(1:nc) :: x, y      ! mole fractions of the various species taken into account.
  real(dp)                  :: P         ! the pressure (in Pa).
  real(dp)                  :: T         ! the temperature (in K).
  real(dp)                  :: rho, rho_MKSA    ! the density of the mixture (kg.m^-3).
  real(dp), dimension(1:nc) :: muskT_k
  real(dp), dimension(1:nc) :: ln_phi, F, delta_2D
  real(dp), dimension(1:nc,1:nc) :: JFn
  
  real(dp), parameter :: pd = 1.E-4_dp       ! Parameter for the computation of derivative by finite differences.
  real(dp), parameter :: precis_NR= 1.E-6_dp ! Precision for the Newton-Raphson algorithm.
  
  integer :: n_Nr
  integer, parameter :: n_NR_max = 100 ! Allowed maximum number of iterations in the Newton-Raphson algo.
  
  real(dp) :: x_N2_1, x_N2_2, Delta_x_N2, delta_2D_max
  real(dp) :: Phi_N2_L, Phi_C2H6_L, Phi_N2_V, Phi_C2H6_V, F1_1, F2_1, F1_2, F2_2
   
  integer  :: iN2, N_N2, i, j, uo, nNR
  character(len=8)  :: date
  character(len=10) :: time
  
  character(len=*), parameter :: file_out_L = 'binary_diag_N2CH4_L.dat' 
  character(len=*), parameter :: file_out_V = 'binary_diag_N2CH4_V.dat' 
  
  !----------------------------------------------------------------------------------------------------------------
  write(6,*) ''
  call wf(' --------------------------------------------------------------------', 'bright red')
  call wf(' -             Simulation of the binary diagram N2-CH4              -', 'bright red')
  call wf(' --------------------------------------------------------------------', 'bright red')
  write(6,*) '' 

  !----------------------------------------------------------------------------------------------------------------
  ! We fix the temperature of the system:
  T= 110._dp ! in K
  !T= 115._dp ! in K
  !T= 120._dp ! in K
  
  compound(1)= 'N2'   ! species 1.
  compound(2)= 'CH4'  ! species 2.

  x_N2_1= 0.0001_dp
  x_N2_2= 0.999_dp
  
  N_N2= 3000
  Delta_x_N2= abs(x_N2_2-x_N2_1)/(N_N2-1)
  
  x(1)= x_N2_1
  x(2)= 1._dp - x(1)
  
  ! Estimation of mole fractions in the vapor phase:
  y(1)= 0.99999_dp
  y(2)= 1._dp - y(1)
  
  P= 0.001_dp * onebar_in_Pa ! First guess for the pressure.
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Initialization of PC-SAFT :
  call wf(' --------------------------------------------------------------------', &
                'bright red')
  call wf(' >>> Initialization of PC-SAFT', 'bright red')
  call pcsaft_PT(subname,1,'L',compound,x, P,T,rho,muskT_k,ln_phi)
  !write(6,*) ' > rho= ', rho
  call wf(' --------------------------------------------------------------------', &
                'bright red')
                
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Computation : 
  open(unit=1,status='unknown',form='formatted',file=trim(file_out_L))
  open(unit=2,status='unknown',form='formatted',file=trim(file_out_V))
  
  call wf(' > We write the output in the file: ', 'bright red', trim(file_out_L), 'bright blue')
  call wf(' > We write the output in the file: ', 'bright red', trim(file_out_V), 'bright blue')
  
  write(6,*) ''
  
  call date_and_time(date,time)
  do uo = 1, 2
        write(uo,'(A)')   '# PC-SAFT-Titan -- Daniel Cordier'
        write(uo,'(A)')   '# https://www.researchgate.net/profile/Daniel_Cordier'
        write(uo,'(A)')   '# http://orcid.org/0000-0003-4515-6271'
        write(uo,'(A)')   '#'
        write(uo,'(A20,A2,A1,A2,A1,A4,A4,A2,a2,A2,a4,A2,a1)') '# Computation date: ',date(7:8),'/',date(5:6),'/', &
              date(1:4), ' -- ', time(1:2),'h ',time(3:4),'min ',time(5:6),'s'
        write(uo,'(A15,F8.2,A2)')   '# Temperature: ', T, ' K'
        write(uo,'(A13,A6)')   '# Species 1: ', compound(1)
        write(uo,'(A13,A6)')   '# Species 2: ', compound(2)
        write(uo,'(A)')   '#'
  end do
  
  write(1,'(2A15)') '# xN2 (liq)    ', '  P (bar)      '
  write(2,'(2A15)') '# yN2 (vap)    ', '  P (bar)      '
     
  do iN2= 1, N_N2-1
     ! 2D Newton-Raphson algo to solve the equations of the liquid-vapor equilibrium:
     delta_2D_max= 1._dp
     !write(6,*) ' > iN2= ', iN2
     nNR = 1
     do while (delta_2D_max > precis_NR)
        call pcsaft_PT(subname,2,'L',compound,x, P,T,rho,muskT_k,ln_phi)
        Phi_N2_L  = exp(ln_phi(1))
        Phi_C2H6_L= exp(ln_phi(2))
        call pcsaft_PT(subname,3,'V',compound,y, P,T,rho,muskT_k,ln_phi)
        Phi_N2_V  = exp(ln_phi(1))
        Phi_C2H6_V= exp(ln_phi(2))
     
        F(1)= Phi_N2_L   * x(1) - Phi_N2_V   * y(1)
        F(2)= Phi_C2H6_L * x(2) - Phi_C2H6_V * (1._dp-y(1))
     
        ! -------- Derivatives / P -----------------------------------------
        P= P * (1._dp-pd)     
        call pcsaft_PT(subname,4,'L',compound,x, P,T,rho,muskT_k,ln_phi)
        Phi_N2_L  = exp(ln_phi(1))
        Phi_C2H6_L= exp(ln_phi(2))
        call pcsaft_PT(subname,5,'V',compound,y, P,T,rho,muskT_k,ln_phi)
        Phi_N2_V  = exp(ln_phi(1))
        Phi_C2H6_V= exp(ln_phi(2))
     
        F1_1= Phi_N2_L   * x(1) - Phi_N2_V   * y(1)
        F2_1= Phi_C2H6_L * x(2) - Phi_C2H6_V * (1._dp-y(1))
     
        P= P / (1._dp-pd) * (1._dp+pd)
        call pcsaft_PT(subname,6,'L',compound,x, P,T,rho,muskT_k,ln_phi)
        Phi_N2_L  = exp(ln_phi(1))
        Phi_C2H6_L= exp(ln_phi(2))
        call pcsaft_PT(subname,7,'V',compound,y, P,T,rho,muskT_k,ln_phi)
        Phi_N2_V  = exp(ln_phi(1))
        Phi_C2H6_V= exp(ln_phi(2))
     
        F1_2= Phi_N2_L   * x(1) - Phi_N2_V   * y(1)
        F2_2= Phi_C2H6_L * x(2) - Phi_C2H6_V * (1._dp-y(1))
     
        P= P / (1._dp+pd) ! Back to initial value.
     
        JFn(1,1)= (F1_2-F1_1)/2._dp/P/pd
        JFn(2,1)= (F2_2-F2_1)/2._dp/P/pd
        ! -------- End compu. derivatives / P -----------------------------------------
     
        ! -------- Begin comp. derivative / y(1) -----------------------------------------
        y(1)= y(1) * (1._dp-pd)
        y(2)= 1._dp - y(1)
        call pcsaft_PT(subname,8,'L',compound,x, P,T,rho,muskT_k,ln_phi)
        Phi_N2_L  = exp(ln_phi(1))
        Phi_C2H6_L= exp(ln_phi(2))
        call pcsaft_PT(subname,9,'V',compound,y, P,T,rho,muskT_k,ln_phi)
        Phi_N2_V  = exp(ln_phi(1))
        Phi_C2H6_V= exp(ln_phi(2))
     
        F1_1= Phi_N2_L   * x(1) - Phi_N2_V   * y(1)
        F2_1= Phi_C2H6_L * x(2) - Phi_C2H6_V * (1._dp-y(1))
     
        y(1)= y(1) / (1._dp-pd) * (1._dp+pd)
        y(2)= 1._dp - y(1)
        call pcsaft_PT(subname,10,'L',compound,x, P,T,rho,muskT_k,ln_phi)
        Phi_N2_L  = exp(ln_phi(1))
        Phi_C2H6_L= exp(ln_phi(2))
        call pcsaft_PT(subname,11,'V',compound,y, P,T,rho,muskT_k,ln_phi)
        Phi_N2_V  = exp(ln_phi(1))
        Phi_C2H6_V= exp(ln_phi(2))
     
        F1_2= Phi_N2_L   * x(1) - Phi_N2_V   * y(1)
        F2_2= Phi_C2H6_L * x(2) - Phi_C2H6_V * (1._dp-y(1))
     
        y(1)= y(1) / (1._dp+pd) ! On revient à la valeur de départ.
        y(2)= 1._dp - y(1)
        
        JFn(1,2)= (F1_2-F1_1)/2._dp/y(1)/pd
        JFn(2,2)= (F2_2-F2_1)/2._dp/y(1)/pd
        ! -------- End comp. derivative / y(1) -----------------------------------------

        ! We solve de 2x2 system:
        call resol_sys2x2 (subname, 1, JFn, -F, delta_2D)
     
        delta_2D_max= abs(delta_2D(1)/P)
        if ( abs(delta_2D(2)/y(1)) > delta_2D_max ) then
           delta_2D_max= abs(delta_2D(2)/y(1))
        end if
        
        ! Mise à jour des variables :
        P   = P    + delta_2D(1)
        y(1)= y(1) + delta_2D(2)
        y(2)= 1._dp - y(1)
        
        nNR= nNR+1 ! NUmber of Newton-Raphson iterations.
        if ( nNR > n_NR_max ) then
           write(6,*) ''
           write(6,*) ' ----------------------------------------------------------------------------------'
           write(6,*) ' > The maximum number of iterations in the Newton-Raphson algo is reached, we stop.'
           write(6,*) '   Please, try to fix the problem.'
           write(6,*) ' ----------------------------------------------------------------------------------'
           write(6,*) ''
           stop
        end if
        
     end do ! "do while" Newton-Raphson algo.
        
     write(1,'(2ES15.7)') x(1), P/onebar_in_Pa ! sortie pour le fichier des données de la phase LIQUIDE.
     write(2,'(2ES15.7)') y(1), P/onebar_in_Pa ! sortie pour le fichier des données de la phase VAPEUR.
     
     x(1)= x(1) + Delta_x_N2
     x(2)= 1._dp - x(1)
     
  end do
  
  close(unit=1)
  close(unit=2)
  close(unit=3)
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Good end message :   
  write(6,*) ''
  call wf(' --------------------------------------------------------------------', &
          'bright red')
  call wf(' > Computation done!', 'bright red')
  write(6,*) ''
   
  stop
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contains

!===================================================================================================================================
subroutine resol_sys2x2 (sub, ncall, A, B, X)
! D. Cordier - July 2021, 7th.
! Solve a 2x2 linear system.
  use utils_dc    
  implicit none
  character(len=*), parameter :: subname = 'resol_sys2x2'
  
  character(len=*),       intent(in) :: sub
  integer,                intent(in) :: ncall
  real(dp), dimension(1:2,1:2), intent(in) :: A
  real(dp), dimension(1:2),     intent(in) :: B
  real(dp), dimension(1:2),    intent(out) :: X
  
  real(dp) :: det
  
  det = A(1,1)*A(2,2) - A(2,1)*A(1,2)

  if ( det == 0._dp ) then
     write(6,*) ' > In the routine "', subname(1:len_trim(subname)), '" det == 0.!'
     stop
  end if
  
  X(1) = (B(1)*A(2,2) - B(2)*A(1,2)) / det
  X(2) = (A(1,1)*B(2) - A(2,1)*B(1)) / det

  return
      
end subroutine resol_sys2x2

end program binary_N2CH4_demo
