!-----------------------------------------------------------------------------------------------------------------------------------
! 
!                                     -- CSTES, TYPE definitions for the PC-SAFT EOS --
!
!-----------------------------------------------------------------------------------------------------------------------------------
! D. Cordier, Institut UTINAM, Besançon, France (until september, 2015).
!             GSMA, Reims, France (from september, 2015).
!
!             - daniel.cordier@univ-reims.fr
!             - Research Gate ...: https://www.researchgate.net/profile/Daniel_Cordier
!             - OrCID ...........: http://orcid.org/0000-0003-4515-6271
!
! -  3 novembre 2013 --: écriture première version.
! - 21  janvier 2014 --: introduction du fichier contenant les donnes pour les paramtres d'interation "k_ij".
! - 18 novembre 2019 --: passage des constantes de "1.3806488d-23"  "1.3806488E-23_dp".
! -  7 mai 2020 -------: introduction d'une nouvelle subroutine 'extract_species_prop' qui permet d'extraire directement
!            dans la base de données 'PC-SAFT' certaines propriétés des espèces étudiées. Cette souroutine permet d'utiliser
!            ces données depuis un programme appelant, avant tout appel à PC-SAFT stricto-sensus.
! - 12 mai 2020 -------: - mise en place de 'ares' comme sortie optionelles de 'Zandmu_pcsaft_PT'.
!                        - mise en place de la sortie de 'h^res/RT' de 'pcsaft_PT'.
! - 30 juin 2020 ------: - mise en place du paramètre logique 'normized_xi_sum' qui contrôle la condition de normalisation des fractions molaires 
!                          à l'entrée de certaines subroutines. Cela permet, par exemple, de calculer des dérivées numériques, par différences
!                          finie, sans se préoccuper de cette condition de normalisation.
! - 11 janvier 2021 ---: - introduction de la subroutine 'alphaV', qui provide the volumetric thermal expansion coefficient (K^-1).
!
module mod_pcsaft
  use utils_DC
  
  ! =========================================================================================================================================
  ! Paramètre qui permet de contrôler la condition de normalisation des fractions molaires à l'entrée de certaines subroutines,
  ! Cela permet, par exemple, de calculer des dérivées numériques, par différences finie, sans se préoccuper de cette condition de normalisation.
  logical, parameter :: normized_xi_sum = .true.  ! In this case, we force the sum of mole fraction to equal the unity !
  !logical, parameter :: normized_xi_sum = .false. !   In this case, we _DO NOT_ force the sum of mole fraction to equal the unity !
  ! =========================================================================================================================================
  
  character(len=*), private, parameter :: mod_name = 'mod_pcsaft' ! Name of this module (PRIVATE).
  
  character(len=*), parameter :: datafilename       = './DATABASE_PCSAFT/COMPOUNDS_DATA_PC-SAFT_withASSO.data' ! PC-SAFT parameters database file name.

                    ! in this file there are:  Molar_mass, Num._seg., Seg._diam., Seg._engergy and References.
  character(len=*), parameter :: datafilename_inter = './DATABASE_PCSAFT/COMPOUNDS_DATA_PC-SAFT-PARAM_INTERAC.data' ! PC-SAFT interaction 
                    ! parameters database file name, in this file one can find the k_ij's and/or the parameters "a" and "b" 
                    ! required to compute the "k_ij" at a given temperature.

  ! Physical constants used this implementation of PC-SAFT:
  integer,  parameter :: comp_name_ncmax = 15   ! Max. number of characters for compound names.
  integer,  parameter :: comp_ref_ncmax  = 100  ! Max. number of characters for compound reference (for PC-SAFT parameters).
  
  real(dp), parameter ::              kB = 1.3806488E-23_dp  ! Boltzmann's constant in J.K^-1
  real(dp), parameter ::  angstrom_meter = 1.E-10_dp         ! one Angstrom in meter
  real(dp), parameter ::            Navo = 6.0221413E+23_dp  ! Avogadro's number
  real(dp), parameter ::    onebar_in_Pa = 1.0000E+5_dp      ! One bar in pascal
  real(dp), parameter :: onekgm3_in_gcm3 = 1.0000E-3_dp      ! One kg.m^-3 in g.cm^-3
  real(dp), parameter ::           R_gas = 8.3144621_dp      ! The gas contant in J.K^-1.mol^-1
  
  ! Definitions of derived types associated to chemical species properties:
  type compdata
     character(len=comp_name_ncmax) :: name      ! The name of the considered compound.
     real(dp)                       :: molmass   ! The molar mass (in kg.mol^-1) for the considered compound.
     real(dp)                       :: numseg    ! The number of segment for the considered compound.
     real(dp)                       :: segdiam   ! The segment diameter.
     real(dp)                       :: segenerg  ! The segment energy
     integer                        :: assoflag  ! Flag: 0=non associative species, 1= associative species
     real(dp)                       :: assenerg  ! The association energy
     real(dp)                       :: assvol    ! The association volume
     character(len=comp_ref_ncmax)  :: ref       ! The reference in which these parameters have been taken.
  end type compdata

  type compdataINTER
     character(len=comp_name_ncmax) :: name1      ! The name of the FIRST  considered compound.
     character(len=comp_name_ncmax) :: name2      ! The name of the SECOND considered compound.
     real(dp)                       :: a          ! The "a" coefficient in Eq. (8) of Tan et al. (2013).
     real(dp)                       :: b          ! The "b" coefficient in Eq. (8) of Tan et al. (2013).
     character(len=comp_ref_ncmax)  :: ref        ! The reference in which these parameters have been taken.
  end type compdataINTER

  type compprop
     character(len=comp_name_ncmax) :: name      ! The name of the considered compound.
     real(dp)                       :: x         ! Mole fraction
     real(dp)                       :: molmass   ! The molar mass (in kg.mol^-1) for the considered compound.
     real(dp)                       :: numseg    ! The number of segment for the considered compound.
     real(dp)                       :: segdiam   ! The segment diameter.
     real(dp)                       :: segenerg  ! The segment energy
     integer                        :: assoflag  ! Flag: 0=non associative species, 1= associative species
     real(dp)                       :: assenerg  ! The association energy
     real(dp)                       :: assvol    ! The association volume
  end type compprop
   
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Subroutines available in this module:
  public  :: read_database_pcsaft, & ! Read the 'PC-SAFT database', i.e. files containing the parameters for numerous species.
             extract_species_prop, & ! Extract species properties from 'PC-SAFT database', useful if we need these propertie outside PC-SAFT.
             P_pcsaft_rho, & ! Provides 'P' when temperature T and density rho are given.
             pcsaft_rhoT, pcsaft_PT, pcsaft_Prho, &
             alphaV, & ! Provide the volumetric thermal expansion coefficient (K^-1)
             cp_solid, cpcv_ssound, cp0_idealgas_joback
  
  !  --- PC-SAFT internal subroutines:
  private :: P_pcsaft_eta, Zandmu_pcsaft_PT, Zandmu_pcsaft_Prho, univ_cst
  private :: association_terms
  private :: intext_linear ! Subroutines utilitaires, également fournie par le module '../../prog_Fortran/INTERPOL/interpol_mod.f08',
                           ! avec lequel elle ne rentre pas en conflit du fait de sa déclaration en 'private' ici.

!===================================================================================================================================
!===================================================================================================================================
contains

!===================================================================================================================================
subroutine extract_species_prop (path, ncall, compound, species_prop)
! Goal of this subroutine: extract the properties of the species specified in 'compound' and leave them in 'species_prop'.
! D. Cordier, 7 mai 2020.
  use utils_dc
  
  implicit none
  character(len=*), parameter :: subname = 'extract_species_prop'
  
  ! Inputs/outputs of the program:
  character(len=*), intent(in) :: path
  integer, intent(in)          :: ncall
  character(len=len_trim(path)+len_trim(subname)+1) :: subtemp
   
  character(len=comp_name_ncmax), dimension(:), intent(in) :: compound     ! List of the names of the involved species.
  type(compprop),                dimension(:), intent(out) :: species_prop ! proprieties associated to involved species.
  
  integer(si) :: i, j
  logical :: first = .true.
  logical :: present
  character(len=6) :: temp_string
  character(len=comp_name_ncmax) :: temp_string2
  type(compdata),      dimension(:), allocatable, save :: compdatabase
  type(compdataINTER), dimension(:), allocatable, save :: compdatabase_interac
  character(len=comp_name_ncmax) :: name1, name2, name3, name4
  
  ! -------------------------------------------------------------------------------------------------------------------------------- 
  ! We read the database of PC-SAFT parameters:
  if (first) then
     call read_database_pcsaft(subname,1,compdatabase,compdatabase_interac)
     first= .false.
     write(6,*) ''
     write(temp_string,'(I4)') size(compdatabase)
     call wf(' > Size of the database: ', 'bright red', &
                        temp_string(1:len_trim(temp_string)), 'bright black', &
                        ' species:', 'bright red')
     do i= 1, size(compdatabase)
        temp_string2= compdatabase(i)%name
        call wf('   - ', 'bright red', &
                        temp_string2(1:len_trim(temp_string2)), 'bright black' )
     end do
     write(6,*) ''
     call wf(' > Species taken into account in this computation: ', 'bright red')
     do i= 1, size(compound)
        temp_string2= compound(i)
        call wf('   - ', 'bright red', &
                        temp_string2(1:len_trim(temp_string2)), 'bright blue' )
     end do
     write(6,*) ''
  end if
  
  do i= 1, size(compound)
     name1  = trim(ADJUSTL(compound(i))) ! Enlève les blancs en début et fin de chaîne.
     present= .false.
     do j= 1, size(compdatabase)
        name2= trim(ADJUSTL(compdatabase(j)%name))
        !write(6,*) name2
        if ( name1(1:len_trim(name1)) == name2(1:len_trim(name2)) ) then
           present= .true.
           species_prop(i)%name     =  name1 ! Name
           species_prop(i)%x        =  -99.  ! Mole fraction, UNDEFINED HERE -> -99. magic number!
           species_prop(i)%molmass  =  compdatabase(j)%molmass  ! The molar mass (in kg.mol^-1) for the considered compound.
           species_prop(i)%numseg   =  compdatabase(j)%numseg   ! The number of segment for the considered compound.
           species_prop(i)%segdiam  =  compdatabase(j)%segdiam  ! The segment diameter.
           species_prop(i)%segenerg =  compdatabase(j)%segenerg ! The segment energy  
           ! Next a series of parameters UNDEFINED here -> magic number.
           species_prop(i)%assoflag =  -99._dp ! Flag: 0=non associative species, 1= associative species
           species_prop(i)%assenerg =  -99._dp ! The association energy
           species_prop(i)%assvol   =  -99._dp ! The association volume
           exit
        end if
     end do
     if ( .NOT. present ) then
        write(6,*) ''
        call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
        write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
        write(6,*) '   the species "', name1(1:len_trim(name1)), '" is not present in the database!'
        write(6,*) ''
        stop
     end if
  end do
  
  return
  
end subroutine extract_species_prop

!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 
!                  -- Computation of "Z" (compressibility factor) and "mu^res" in the frame of the PC-SAFT theory --
!
!-----------------------------------------------------------------------------------------------------------------------------------
! D. Cordier, Institut UTINAM, Besançon, France.
! 1er novembre 2013 -- daniel.cordier@obs-besancon.fr
!
! References: * Sunil Kumar Maity, "Modeling and simulation of solid-liquid equilibrium by
!               perturbed-chain statistical associating fluid theory", Master of Technology in Chemical Engineering,
!               2003,
!               Indian Institute of Technology, Kharagpur, India.
!             * Abdelkrim Belkadi, Thse de l'Universit de Toulouse, 2008.
!             * Gross & Sadowski, Ind. Eng. Chem. Res. 2001, 40, 1244-1260.
!
! - June 2sd 2014: change the name 'Zandmu_pcsaft' -> 'Zandmu_pcsaft_PT'.
! - 27 novembre 2014 : correction du nom de la routine dans 'subname'.
!
! - 15 avril 2020 : - passage des constantes numriques de l'criture '1.d0' vers l'criture '1._dp'
!                   - mise dans le module 'NEW_MOD_PCSAFT'.
!
! - 12 mai 2020   : - introduction d'une sortie optionelle : 'ares'.
!
subroutine Zandmu_pcsaft_PT(sub, ncall, state, species_prop, k_ij, P, T, Z, muskT_k, eta, rho, ares)
  use utils_DC
  
  implicit none
  character(len=*), parameter :: subname = 'Zandmu_pcsaft_PT'
  
  ! Inputs/outputs of the program:
  character(len=*), intent(in) :: sub
  integer, intent(in)          :: ncall
  character(len=len_trim(sub)+len_trim(subname)+1) :: subtemp

  character(len=1),             intent(in) :: state     ! the physical state of the phase:  - state='L' = liquid phase
                                                           !                                   - state='V' = vapor phase.
  type(compprop), dimension(:), intent(in) :: species_prop ! The properties of studied species.
  real(dp),     dimension(:,:), intent(in) :: k_ij    ! interaction parameters at the given temperature T
  real(dp), intent(in)                     :: P, T    ! Pressure (in Pa) and temperature (in K).
  real(dp), intent(out)                    :: Z       ! Compressibility factor (no unit)
  real(dp), dimension(:),      intent(out) :: muskT_k ! The residual chemical potential divided by kT (no unit)
  real(dp), intent(out)                    :: eta     ! the reduced density (no unit)
  real(dp), intent(out)                    :: rho     ! the density (kg.m^-3) of the system at P and T
  real(dp), optional, intent(out)          :: ares    ! Residual (compare to ideal gas value) Helmoltz energy (J.mol^-1).
  
  integer             :: ntour_NR
  integer, parameter  :: ntour_NR_max = 100      ! Max. number of tours in the Newton-Raphson's algorithm.
  real(dp), parameter :: delta_NR_max = 1.E-8_dp ! Convergence criterion for the Newton-Raphson's algorithm.
  real(dp), parameter :: pd = 1.E-5_dp
  real(dp) :: Psaft, Psaft1, Psaft2, Psaft_0, eta1, eta2, rho1, rho2, Z1, Z2, dP_pcsaftsdeta, Delta_P
  real(dp) :: ln_eta, ln_eta1, ln_eta2, Delta_P1, Delta_P2
  logical  :: increase

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Call path of this subroutine:
  subtemp= sub(1:len_trim(sub))//'/'//subname(1:len_trim(subname))
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! We determine the value of "eta" that matches the pressure of the system P (i.e. to get Psaft=P):
  !   
  if ( state == 'L' ) then
     eta= 0.5_dp
  end if
  if ( state == 'V' ) then
     eta= 1.E-10_dp
  end if
 
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Recherche par Newton-Raphson :
  ! Initialisation de la routine 'P_pcsaft_eta' :

  call P_pcsaft_eta(subtemp,1,species_prop,k_ij,T,eta,Psaft_0,rho,Z,muskT_k)
   
  ! The Newton-Raphson's algorithm itself:
  Delta_P= Psaft_0/P-1._dp
  ntour_NR= 1
  do while ( abs(Delta_P) > delta_NR_max)
     call P_pcsaft_eta(subtemp,2,species_prop,k_ij,T,eta,Psaft_0,rho,Z,muskT_k)  
     Delta_P= Psaft_0/P-1._dp
     
     eta1= eta * (1._dp-pd)
     call P_pcsaft_eta(subtemp,3,species_prop,k_ij,T,eta1,Psaft1,rho1,Z1,muskT_k)     
     eta2= eta * (1._dp+pd)
     call P_pcsaft_eta(subtemp,4,species_prop,k_ij,T,eta2,Psaft2,rho2,Z2,muskT_k)
     
     dP_pcsaftsdeta= (Psaft2-Psaft1)/(eta2-eta1)/P
  
     eta= eta - Delta_P/dP_pcsaftsdeta

     !if ( eta < 0.d0 ) then
     !   print*, ' eta < 0.d0 on impose eta= 1.d-10 ! (Dans Zandmu_pcsaft)'
     !   eta= 1.d-20 
     !end if
     
     ntour_NR= ntour_NR + 1
     if ( ntour_NR > ntour_NR_max ) then
        write(6,*) ''
        call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
        write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
        write(6,*) '   the number of iteration in the Newton-Raphson''s algorithm is larger than the allowed maximum value!'
        write(6,*) '   ntour_NR     = ', ntour_NR
        write(6,*) '   ntour_NR_max = ', ntour_NR_max
        write(6,*) ''
        stop
     end if
  end do
  !
  ! Fin recherche par Newton-Raphson
  !---------------------------------------------------------------------------------------------------------------------------------
   
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Recherche par dichotomie :
  !!
  !!call P_pcsaft_eta(subname,1,species_prop,T,eta1,Psaft1,rho,Z,muskT_k)
  !!Delta_P1= Psaft1/P-1.d0
  !!
  !!call P_pcsaft_eta(subtemp,2,species_prop,T,eta2,Psaft2,rho,Z,muskT_k)
  !!Delta_P2= Psaft2/P-1.d0
  !!
  !!if ( Delta_P1 * Delta_P2 >= 0.d0 ) then
  !!   write(6,*) ''
  !!   write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
  !!   write(6,*) '   Delta_P1 * Delta_P2 => 0.d0'
  !!   write(6,*) '   the bisection method can be applied here, please check the values of "eta1" and "eta2"'
  !!   write(6,*) ''
  !!   stop
  !!end if
  !!
  !!if ( (Delta_P1 < 0.d0) .AND. (Delta_P2 > 0.d0) ) then
  !!   increase= .true.
  !!else
  !!   increase= .false.
  !!end if
  !!
  !!ntour_NR= 0
  !!eta= (eta1+eta2)/2.d0
  !!do while (abs(eta1-eta2)/eta > 1.d-6)
  !!   eta= (eta1+eta2)/2.d0
  !!   call P_pcsaft_eta(subtemp,3,species_prop,T,eta,Psaft_0,rho,Z,muskT_k)
  !!   Delta_P= Psaft_0/P-1.d0
  !!   if (increase) then
  !!      if (Delta_P > 0.d0) then
  !!         eta2= eta
  !!      else
  !!         eta1= eta
  !!      end if
  !!   else
  !!      if (Delta_P > 0.d0) then
  !!         eta1= eta
  !!      else
  !!         eta2= eta
  !!      end if
  !!   end if
  !!   !print*, ' eta= ', eta
  !!   !read(5,*)
  !!   ntour_NR= ntour_NR + 1
  !!   if ( ntour_NR > ntour_NR_max ) then
  !!     write(6,*) ''
  !!      write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', sub(1:len_trim(sub)), ''''
  !!      write(6,*) '   the number of iteration in the bisection method''s algorithm is larger than the allowed maximum value!'
  !!      write(6,*) '   ntour_NR     = ', ntour_NR
  !!      write(6,*) '   ntour_NR_max = ', ntour_NR_max
  !!      write(6,*) ''
  !!      stop
  !!   end if
  !!end do
  
  ! Fin recherche par dichotomie.
  !---------------------------------------------------------------------------------------------------------------------------------
   
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Final values of "rho", "Z", "muskT_k":
  if (present(ares)) then
    call P_pcsaft_eta(subtemp, 5, species_prop, k_ij, T, eta, Psaft_0, rho, Z, muskT_k, ares)   
  else
    call P_pcsaft_eta(subtemp, 6, species_prop, k_ij, T, eta, Psaft_0, rho, Z, muskT_k)
  end if
  
  !---------------------------------------------------------------------------------------------------------------------------------
    
  return

end subroutine Zandmu_pcsaft_PT

!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 
!                  -- Computation of "Z" (compressibility factor) and "mu^res" in the frame of the PC-SAFT theory --
!
!-----------------------------------------------------------------------------------------------------------------------------------
! D. Cordier, Institut UTINAM, Besançon, France.
! June 2sd 2014 -- daniel.cordier@obs-besancon.fr
!
! References: * Sunil Kumar Maity, "Modeling and simulation of solid-liquid equilibrium by
!               perturbed-chain statistical associating fluid theory", Master of Technology in Chemical Engineering,
!               2003,
!               Indian Institute of Technology, Kharagpur, India.
!             * Abdelkrim Belkadi, Thse de l'Universit de Toulouse, 2008.
!             * Gross & Sadowski, Ind. Eng. Chem. Res. 2001, 40, 1244-1260.
!
! - 27 novembre 2014 : correction du nom de la routine dans 'subname'.
!
! - 15 avril 2020 : - passage des constantes numriques de l'criture '1.d0' vers l'criture '1._dp'
!                   - mise dans le module 'NEW_MOD_PCSAFT'.
!
subroutine Zandmu_pcsaft_Prho(sub,ncall,state,species_prop,k_ij,P,rho,Z,muskT_k,T,ares)
  use utils_DC
  
  implicit none
  character(len=*), parameter :: subname = 'Zandmu_pcsaft_Prho'
  
  ! Inputs/outputs of the program:
  character(len=*), intent(in) :: sub
  integer, intent(in)          :: ncall
  character(len=len_trim(sub)+len_trim(subname)+1) :: subtemp
  
  character(len=1),             intent(in) :: state     ! the physical state of the phase:  - state='L' = liquid phase
                                                           !                                   - state='V' = vapor phase.
  type(compprop), dimension(:), intent(in) :: species_prop ! The properties of studied species.
  real(dp),     dimension(:,:), intent(in) :: k_ij    ! interaction parameters at the given temperature T
  real(dp), intent(in)                     :: P, rho    ! Pressure (in Pa) and density (in kg.m^-3).
  real(dp), intent(out)                    :: Z       ! Compressibility factor (no unit)
  real(dp), dimension(:),      intent(out) :: muskT_k ! The residual chemical potential divided by kT (no unit)
  real(dp), intent(out)                    :: T       ! the temperature (in K) 
  real(dp),          intent(out), optional :: ares    ! Residual (compare to ideal gas value) Helmoltz energy (in J.mol^-1)
  
  integer             :: ntour_NR
  integer, parameter  :: ntour_NR_max = 100000   ! Max. number of tours in the Newton-Raphson's algorithm.
  real(dp), parameter :: delta_NR_max = 1.E-6_dp ! Convergence criterion for the Newton-Raphson's algorithm.
  real(dp), parameter :: pd = 1.E-5_dp
  real(dp) :: Psaft, Psaft1, Psaft2, Psaft_0, temp, temp1, temp2, Z1, Z2, dP_pcsaftsdtemp, Delta_P

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Call path of this subroutine:
  subtemp= sub(1:len_trim(sub))//'/'//subname(1:len_trim(subname))
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! We determine the value of "eta" that matches the pressure of the system P (i.e. to get Psaft=P):
  !   
  !if ( state == 'L' ) then
  !   eta= 0.5
  !end if
  !if ( state == 'V' ) then
  !   eta= 1.d-10
  !end if
 
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Recherche par Newton-Raphson :
  ! Initialisation de la routine 'P_pcsaft_rho' :

  temp= 273._dp ! On initialise la valeur de la température.
  call P_pcsaft_rho(subtemp,1,species_prop,k_ij,temp,rho,Psaft_0,Z,muskT_k)
   
  ! The Newton-Raphson's algorithm itself:
  Delta_P= Psaft_0/P-1._dp
  ntour_NR= 1

  do while ( abs(Delta_P) > delta_NR_max)
     call P_pcsaft_rho(subtemp,2,species_prop,k_ij,temp,rho,Psaft_0,Z,muskT_k)
     Delta_P= Psaft_0/P-1._dp
     
     temp1= temp * (1._dp-pd)
     call P_pcsaft_rho(subtemp,3,species_prop,k_ij,temp1,rho,Psaft1,Z,muskT_k)     
     
     temp2= temp * (1._dp+pd)
     call P_pcsaft_rho(subtemp,4,species_prop,k_ij,temp2,rho,Psaft2,Z,muskT_k)
     
          
     dP_pcsaftsdtemp= (Psaft2-Psaft1)/(temp2-temp1)/temp
  
     temp= temp - Delta_P/dP_pcsaftsdtemp     

     !print*, ' temp= ', temp
     !read(5,*)
     
     ntour_NR= ntour_NR + 1
     if ( ntour_NR > ntour_NR_max ) then
        write(6,*) ''
        call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
        write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
        write(6,*) '   the number of iteration in the Newton-Raphson''s algorithm is larger than the allowed maximum value!'
        write(6,*) '   ntour_NR     = ', ntour_NR
        write(6,*) '   ntour_NR_max = ', ntour_NR_max
        write(6,*) ''
        stop
     end if
  end do
  !
  ! Fin recherche par Newton-Raphson
  !---------------------------------------------------------------------------------------------------------------------------------
   
   
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Final values of "rho", "Z", "muskT_k":
  if (present(ares)) then
 	 call P_pcsaft_rho(subtemp,5,species_prop,k_ij,temp,rho,Psaft2,Z,muskT_k,ares)
  else
     call P_pcsaft_rho(subtemp,6,species_prop,k_ij,temp,rho,Psaft2,Z,muskT_k)
  end if
  T= temp
  
  !---------------------------------------------------------------------------------------------------------------------------------
    
  return

end subroutine Zandmu_pcsaft_Prho

!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 
!                                   -- PC-SAFT parameters: reading of data files --
!
!-----------------------------------------------------------------------------------------------------------------------------------
! D. Cordier, Institut UTINAM, Besançon, France.
! daniel.cordier@obs-besancon.fr
!
! -  3 novembre 2013  : criture de la premire version.
! - 21  janvier 2014  : mise en place de la lecture du deuxme fichier (celui contenant les donnes des paramtres d'interaction
!                       "k_ij".
! - 18 septembre 2014 : scurit : dtection des doublons dans le fichier de 'k_ij'.
! - 27 septembre 2016 : rcriture des "do while" pour les lectures de fichiers de longueurs indtermines.
!
! - 15 avril 2020 : - mise dans le module 'NEW_MOD_PCSAFT'.
!
subroutine read_database_pcsaft(sub, ncall, compdatabase, compdatabase_interac)
  use utils_DC
  
  implicit none
  character(len=*), parameter :: subname = 'read_database_pcsaft'
  
  ! Inputs/outputs of the program:
  character(len=*), intent(in) :: sub
  integer, intent(in)          :: ncall

  type(compdata),      dimension(:), allocatable, intent(out) :: compdatabase
  type(compdataINTER), dimension(:), allocatable, intent(out) :: compdatabase_interac
  
  integer :: nl, ndieze, ntot, ok, i, j, iostatus
  character(len=1)   :: line1
  character(len=300) :: line
  logical :: here
  !

!-----------------------------------------------------------------------------------------------------------------------------------
! Readin of the FIRST data file:
  
  ! We check that the PC-SAFT parameters database is really present in the directory of work:
  inquire(file=datafilename(1:len_trim(datafilename)),exist=here)
  if ( .NOT. here ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', sub(1:len_trim(sub))
     write(6,*) '   the file ''', datafilename(1:len_trim(datafilename)), ''' is not present'
     write(6,*) ''
     stop
  end if
  
  ! We count the number of species taken into account by the database:
  nl=0; ndieze=0
  open(unit=1,status='unknown',form='formatted',file=datafilename(1:len_trim(datafilename)))
  do
     read(1,'(A)',iostat=iostatus) line1
     if (iostatus < 0) then
        exit
     end if
     if ( line1 /= '#' ) then
        nl=nl+1
     else
        ndieze=ndieze+1
     end if
  end do
  close(unit=1)

  ! We allocate the memory for the database:
  allocate(compdatabase(1:nl),stat=ok);       call alloc_error(subname,'compdatabase',ok)
  
  ! We read the database file:
  open(unit=1,status='unknown',form='formatted',file=datafilename(1:len_trim(datafilename)))
  ntot=nl+ndieze
  nl=0
  do i= 1, ntot
     read(1,'(A)') line
     if ( line(1:1) /= '#' ) then
        nl=nl+1
!#            |   (kg.mol^-1)  |    "m"    | "sig" (Ang) | e/kB (K)     | 
!N2                28.01340E-03      1.2414        3.2992        89.2230   NIST Used by Tan et al. (2013)
        read(line,'(A15,ES17.5,F12.4,F14.4,F15.4,A100)') compdatabase(nl)%name, compdatabase(nl)%molmass, &
                                                        compdatabase(nl)%numseg, compdatabase(nl)%segdiam, &
                                                        compdatabase(nl)%segenerg, compdatabase(nl)%ref
     end if
  end do
  close(unit=1)

!-----------------------------------------------------------------------------------------------------------------------------------
! Readin of the SECOND data file:
                       
  ! We check that the PC-SAFT parameters database is really present in the directory of work:
  inquire(file=datafilename_inter(1:len_trim(datafilename_inter)),exist=here)
  if ( .NOT. here ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', sub(1:len_trim(sub))
     write(6,*) '   the file ''', datafilename_inter(1:len_trim(datafilename_inter)), ''' is not present'
     write(6,*) ''
     stop
  end if
  
  ! We count the number of species taken into account by the database:
  nl=0; ndieze=0
  open(unit=1,status='unknown',form='formatted',file=datafilename_inter(1:len_trim(datafilename_inter)))
  do
     read(1,'(A)',iostat=iostatus) line1
     if(iostatus < 0) then
       exit
     end if
     if ( line1 /= '#' ) then
        nl=nl+1
     else
        ndieze=ndieze+1
     end if
  end do
  close(unit=1)

  ! We allocate the memory for the database:
  allocate(compdatabase_interac(1:nl),stat=ok);     call alloc_error(subname,'compdatabase_interac',ok)

  ! We read the database file:
  open(unit=1,status='unknown',form='formatted',file=datafilename_inter(1:len_trim(datafilename_inter)))
  ntot=nl+ndieze
  nl=0
  do i= 1, ntot
     read(1,'(A)') line
     if ( line(1:1) /= '#' ) then
        nl=nl+1
        read(line,'(2A15,2F12.4,A100)') compdatabase_interac(nl)%name1, compdatabase_interac(nl)%name2, &
                                        compdatabase_interac(nl)%a, compdatabase_interac(nl)%b, &
                                        compdatabase_interac(nl)%ref
        !write(6,'(2A15,2F12.4,A100)') compdatabase_interac(nl)%name1, compdatabase_interac(nl)%name2, &
        !                              compdatabase_interac(nl)%a, compdatabase_interac(nl)%b, &
        !                              compdatabase_interac(nl)%ref
     end if
  end do  
  close(unit=1)
  
  ! We search for duplicates:
  do i= 1, nl
     !write(6,*) '''', compdatabase_interac(i)%name1, '''', '|', '''', compdatabase_interac(i)%name2, ''''
     do j= 1, nl
        ! Cas o deux lignes se rptent simplement :
        if ( (compdatabase_interac(i)%name1 .eq. compdatabase_interac(j)%name1) .AND. (i /= j) ) then
           if ( compdatabase_interac(i)%name2 .eq. compdatabase_interac(j)%name2 ) then
              write(6,*) ''
              call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
              write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', sub(1:len_trim(sub))
              write(6,*) '   there is a duplicate:'
              write(6,*) '   ', i, compdatabase_interac(i)%name1, compdatabase_interac(i)%name2
              write(6,*) '   ', j, compdatabase_interac(j)%name1, compdatabase_interac(j)%name2
              write(6,*) ''
              stop
           end if
        end if
        ! 
        if ( (compdatabase_interac(i)%name1 .eq. compdatabase_interac(j)%name2) .AND. (i /= j) ) then
           if ( compdatabase_interac(i)%name2 .eq. compdatabase_interac(j)%name1 ) then
              write(6,*) ''
              call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
              write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', sub(1:len_trim(sub))
              write(6,*) '   there is a duplicate:'
              write(6,*) '   ', i, compdatabase_interac(i)%name1, compdatabase_interac(i)%name2
              write(6,*) '   ', j, compdatabase_interac(j)%name1, compdatabase_interac(j)%name2
              write(6,*) ''
              stop
           end if
        end if
     end do
  end do
  
  return
  
end subroutine read_database_pcsaft

!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 
!                           -- Computation of the pressure P in the frame of the PC-SAFT theory --
!
!                                       >>>    AS A FUNCTION OF "rho" and "T"   <<<
!
!-----------------------------------------------------------------------------------------------------------------------------------
! D. Cordier, Institut UTINAM, Besançon, France.
! June 2sd 2014 -- daniel.cordier@obs-besancon.fr
!
! References: * Sunil Kumar Maity, "Modeling and simulation of solid-liquid equilibrium by
!               perturbed-chain statistical associating fluid theory", Master of Technology in Chemical Engineering,
!               2003,
!               Indian Institute of Technology, Kharagpur, India.
!               WARNING: a lot of mistakes/typos in this references! 
!             * Abdelkrim Belkadi, Thse de l'Universit de Toulouse, 2008.
!             * Gross & Sadowski, Ind. Eng. Chem. Res. 2001, 40, 1244-1260.
!
! - 15 avril 2020 : - passage des constantes numriques de l'criture '1.d0' vers l'criture '1._dp'
!                   - mise dans le module 'NEW_MOD_PCSAFT'.
!
subroutine P_pcsaft_rho(sub, ncall, species_prop, k_ij, T, rho_MKSA, Psaft, Z, muskT_k, ares)
  use utils_DC
  
  implicit none
  character(len=*), parameter :: subname = 'P_pcsaft_rho'
  
  ! Inputs/outputs of the program:
  character(len=*), intent(in) :: sub
  integer, intent(in)          :: ncall
  character(len=len_trim(sub)+len_trim(subname)+1) :: subtemp

  type(compprop), dimension(:), intent(in) :: species_prop ! The properties of studied species.
  real(dp),     dimension(:,:), intent(in) :: k_ij         ! interaction parameters at the given temperature T
  real(dp), intent(in)                     :: T            ! Temperature (in K).
  real(dp), intent(in)                     :: rho_MKSA     ! the density (kg.m^-3) 
  real(dp), intent(out)                    :: Psaft        ! The pressure (in Pa) computed in the frame of the PC-SAFT
  real(dp), intent(out)                    :: Z            ! The compressibility factor (no unit)
  real(dp), dimension(:), intent(out)      :: muskT_k      ! The residual chemical potential divided by kT (no unit)
  real(dp),               intent(out), optional :: ares    ! Residual (compare to ideal gas value) Helmoltz energy (J.mol^-1).
  
  real(dp) :: rho ! density in "number of molecules / angstrom^3 (NMA)".
  real(dp) :: eta ! reduced density (no unit).
  
  logical :: opt_outputs
  
  integer  :: i, j, k, n, ok
  real(dp), save :: m_bar
  real(dp), dimension(:), allocatable :: d
  real(dp), dimension(0:3) :: zeta
  real(dp), dimension(:,:), allocatable :: zeta_xk
  
  real(dp), save :: sum_ximidi0, sum_ximidi1, sum_ximidi2, sum_ximidi3

  real(dp), dimension(:),   allocatable :: gii_hs
  real(dp), dimension(:),   allocatable :: rho_dgii_hs
  real(dp), dimension(:,:), allocatable :: dghs_ii_dxk

  real(dp) :: C_1
  real(dp), dimension(:,:), allocatable :: sigma, epsilon
  real(dp), dimension(0:2,0:6) :: a, b
  real(dp), dimension(0:6) :: a_i, b_i
  real(dp) :: I_1, I_2, detaI1sdeta, detaI2sdeta
  real(dp) :: C_2
  real(dp) :: m2_epsilon_sigma3_bar, m2_epsilon2_sigma3_bar
  real(dp) :: Z_hs, Z_hc, Z_disp, Zassoc
  real(dp) :: molar_mass_bar
  real(dp) :: a_res, a_hc, a_hs, a_disp ! reduced Helmholtz energy (no unit)
  real(dp) :: sum_xi_mim1_ln_gii, sum_xi_mim1_ghsiim1_dghsiisxk, sum_xj_dares_dxj
  real(dp), dimension(:),   allocatable :: da_hs_s_dxk, dahc_dxk, dares_dxk
  real(dp), dimension(:),   allocatable :: m2_epsilon_sigma3_bar_xk, m2_epsilon2_sigma3_bar_xk
  real(dp), dimension(:),   allocatable :: C_1_xk
  real(dp), dimension(:,:), allocatable :: ai_xk, bi_xk
  real(dp), dimension(:),   allocatable :: I_1_xk, I_2_xk, dadisp_xk
  real(dp), dimension(:),   allocatable :: muskT_k_assoc

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Call path of this subroutine:
  subtemp= sub(1:len_trim(sub))//'/'//subname(1:len_trim(subname))
  
  !---------------------------------------------------------------------------------------------------------------------------------
  if ( rho_MKSA <= 0._dp ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the considered density "rho" is <= 0.!'
     write(6,*) '   We have: "rho" = ', rho_MKSA
     write(6,*) ''
     stop
  end if
    
  !---------------------------------------------------------------------------------------------------------------------------------
  ! We compute the mean number of segment "m_bar" and various quantities used in the following:
     m_bar= 0._dp ! cf. Gross & Sadowski (2001)
     do i= 1, size(species_prop)
        m_bar= m_bar + species_prop(i)%x * species_prop(i)%numseg
     end do
     
     ! Computation of the d_i's:
     allocate(d(1:size(species_prop)),stat=ok);                       call alloc_error(subname,'d',ok)
     
     do i= 1, size(species_prop)
        ! Eq. (4.2.9) de Maity (2003) et aussi Eq. [II-136] p. 50 de Belkadi (2008)
        d(i)= species_prop(i)%segdiam * (1._dp - 0.12_dp * exp(-3._dp * species_prop(i)%segenerg / T)) ! En Angstrom
        !write(6,*) ' d(',i,')= ', d(i)
     end do

     sum_ximidi0= 0._dp
     do i= 1, size(species_prop)
        ! cf. Eq. (4.2.8) p. 26 Maity (2003) and Eq. [II-131] p. 49 Belkadi (2008)
        sum_ximidi0= sum_ximidi0 + species_prop(i)%x * species_prop(i)%numseg
     end do
     
     sum_ximidi1= 0._dp
     do i= 1, size(species_prop)
        ! cf. Eq. (4.2.8) p. 26 Maity (2003) and Eq. [II-131] p. 49 Belkadi (2008)
        sum_ximidi1= sum_ximidi1 + species_prop(i)%x * species_prop(i)%numseg * d(i)
     end do
     
     sum_ximidi2= 0._dp
     do i= 1, size(species_prop)
        ! cf. Eq. (4.2.8) p. 26 Maity (2003) and Eq. [II-131] p. 49 Belkadi (2008)
        sum_ximidi2= sum_ximidi2 + species_prop(i)%x * species_prop(i)%numseg * d(i)**2
     end do
     
     sum_ximidi3= 0._dp
     do i= 1, size(species_prop)
        ! cf. Eq. (4.2.20) Maity (2003) and Eq. [II-124] p. 48 Belkadi (2008)
        sum_ximidi3= sum_ximidi3 + species_prop(i)%x * species_prop(i)%numseg * d(i)**3
     end do

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Computation of the density "rho" from "MKSA" to "number of molecules / angstrom^3 (NMA)" :
  molar_mass_bar= 0._dp
  do i= 1, size(species_prop)
     molar_mass_bar= molar_mass_bar + species_prop(i)%x * species_prop(i)%molmass ! Mean molar mass in kg.mol^-1
  end do

  rho= rho_MKSA / molar_mass_bar * Navo * angstrom_meter**3._dp ! density in "number of molecules / angstrom^3 (NMA)"
  
  ! Reduced density (no unit):
  eta= rho/6._dp*pi * sum_ximidi3 ! cf. Eq. (4.2.20) Maity (2003) and Eq. [II-124] p. 48 Belkadi (2008)
                                 ! unit: number of molecules / angstrom^3 ! V
                                 
  !---------------------------------------------------------------------------------------------------------------------------------
  
  allocate(gii_hs(1:size(species_prop)),stat=ok);                           call alloc_error(subname,'gii_hs',ok)
  allocate(rho_dgii_hs(1:size(species_prop)),stat=ok);                      call alloc_error(subname,'rho_dgii_hs',ok)
  allocate(sigma(1:size(species_prop),1:size(species_prop)),stat=ok);       call alloc_error(subname,'sigma',ok)
  allocate(epsilon(1:size(species_prop),1:size(species_prop)),stat=ok);     call alloc_error(subname,'epsilon',ok)
  allocate(zeta_xk(0:3,1:size(species_prop)),stat=ok);                      call alloc_error(subname,'zeta_xk',ok)
  allocate(da_hs_s_dxk(1:size(species_prop)),stat=ok);                      call alloc_error(subname,'da_hs_s_dxk',ok)
  allocate(dahc_dxk(1:size(species_prop)),stat=ok);                         call alloc_error(subname,'dahc_dxk',ok)
  allocate(dares_dxk(1:size(species_prop)),stat=ok);                        call alloc_error(subname,'dares_dxk',ok)
  allocate(dghs_ii_dxk(1:size(species_prop),1:size(species_prop)),stat=ok); call alloc_error(subname,'dghs_ii_dxk',ok)
  allocate(m2_epsilon_sigma3_bar_xk(1:size(species_prop)),stat=ok);         call alloc_error(subname,'m2_epsilon_sigma3_bar_xk',ok)
  allocate(m2_epsilon2_sigma3_bar_xk(1:size(species_prop)),stat=ok);        call alloc_error(subname,'m2_epsilon2_sigma3_bar_xk',ok)
  allocate(C_1_xk(1:size(species_prop)),stat=ok);                           call alloc_error(subname,'C_1_xk',ok)
  allocate(ai_xk(0:6,1:size(species_prop)),stat=ok);                        call alloc_error(subname,'ai_xk',ok)
  allocate(bi_xk(0:6,1:size(species_prop)),stat=ok);                        call alloc_error(subname,'bi_xk',ok)
  allocate(I_1_xk(1:size(species_prop)),stat=ok);                           call alloc_error(subname,'I_1_xk',ok)
  allocate(I_2_xk(1:size(species_prop)),stat=ok);                           call alloc_error(subname,'I_2_xk',ok)
  allocate(dadisp_xk(1:size(species_prop)),stat=ok);                        call alloc_error(subname,'dadisp_xk',ok)
    
  !---------------------------------------------------------------------------------------------------------------------------------
  !!! Beginning of the computation of the compressibility factor "Z":
  !!!rho= 6.d0/pi * eta / sum_ximidi3  ! cf. Eq. (4.2.20) Maity (2003) and Eq. [II-124] p. 48 Belkadi (2008)
  !!!                                  ! unit: number of molecules / angstrom^3 ! V
    
  ! Voir Eq. (9) de Gross & Sadowski (2001)
  zeta(0)= pi/6._dp * rho * sum_ximidi0 ! Angstrom^-3   ! V                             
  zeta(1)= pi/6._dp * rho * sum_ximidi1 ! Angstrom^-2   ! V
  zeta(2)= pi/6._dp * rho * sum_ximidi2 ! Angstrom^-1   ! V
  zeta(3)= pi/6._dp * rho * sum_ximidi3 ! no unit       ! V

  ! Eq. (8) de Gross & Sadowski (2001) et ((4.2.26 de MAity (2003) :
  Z_hs= zeta(3)/(1._dp-zeta(3)) + &
        3._dp*zeta(1)*zeta(2)/zeta(0)/(1._dp-zeta(3))**2 + &
        (3._dp*zeta(2)**3-zeta(3)*zeta(2)**3)/zeta(0)/(1._dp-zeta(3))**3 ! V
           
  ! We compute the term "g_ii^hs", cf. for instance (4.2.7) p. 26 of Maity (2003):
  do i= 1, size(species_prop)
     ! ATTENTION: cette formule est (4.2.7) p. 26 de Maity (2003) en accord avec [II-129] p. 49 de Belkadi (2008) mais
     !            dans le papier original de Gross & Sadowski (2001) il n'y a pas de carr au dnominateur du deuxime terme.
     !            en fait il y a de carr dans la formule (A.7) de G&S.
     gii_hs(i)= 1._dp/(1._dp-zeta(3)) + d(i)*d(i)/(d(i)+d(i)) * 3._dp*zeta(2)/(1._dp-zeta(3))**2 + (d(i)*d(i)/(d(i)+d(i)))**2 * &
             2._dp*zeta(2)**2/(1._dp-zeta(3))**3 ! V
  end do  
  
  do i= 1, size(species_prop)
        ! Voir (4.2.27) p. 29 de Maity (2003) et (A.27) de G&S:
        j=i
        rho_dgii_hs(i)= zeta(3)/(1._dp-zeta(3))**2 + d(i)*d(j)/(d(i)+d(j)) * &
                         (3._dp*zeta(2)/(1._dp-zeta(3))**2 + 6._dp*zeta(2)*zeta(3)/(1._dp-zeta(3))**3) + &
                         (d(i)*d(j)/(d(i)+d(j)))**2 * &
                         (4._dp*zeta(2)**2/(1._dp-zeta(3))**3 + 6._dp*zeta(2)**2*zeta(3)/(1._dp-zeta(3))**4) ! V
  end do
  
  ! On calcule Z_hc:
  ! Eq. (4.2.25) p. 28 de Maity (2003) et aussi Eq. (5) de Gross & Sadowski (2001):
  Z_hc= m_bar* Z_hs
  do i= 1, size(species_prop)
     Z_hc= Z_hc - species_prop(i)%x * (species_prop(i)%numseg-1._dp) /gii_hs(i) * rho_dgii_hs(i) ! V
  end do
  
  ! We compute Z^disp:
  ! Eq. (4.2.11) p. 26 de Maity (2003), ATTENTION dans G&S il manque l'exposant <<-1>> !
  C_1=  1._dp/(1._dp + m_bar * (8._dp*eta-2._dp*eta**2)/(1._dp-eta)**4 + (1._dp-m_bar) * &
        (20._dp*eta-27._dp*eta**2+12._dp*eta**3-2._dp*eta**4)/(1._dp-eta)**2/(2._dp-eta)**2) ! V
        
  do i= 1, size(species_prop)
     do j= 1, size(species_prop)
        sigma(i,j)  = 1._dp/2._dp * (species_prop(i)%segdiam + species_prop(j)%segdiam)
        ! ATTENTION 1 : on met une valeur fixe pour k_ij --->>> il faudrait faire une table initiale o on a tous les k_ij
        !               ce que je mets l est destin aux premiers tests.
        ! ATTENTION 2 : - la formule utlise ici est celle de Maity (2003) (4.2.15) p. 27
        !               - Belkadi (2008) en propose une autre (cf. p. 38)
        !               - Gross & Sadowski (2001) Eq. (24) est encore autre chose !

        ! ATTENTION : les calculs des 'k_ij'  la temprature de travail 'T' se fait dans la routine 'pcsaft'
        epsilon(i,j)= sqrt(species_prop(i)%segenerg * species_prop(j)%segenerg) * (1._dp - k_ij(i,j))
        !print*, ' k_ij(i,j)= ', k_ij(i,j)
     end do
  end do
  
  call univ_cst(subname,1,a,b)
  
  do i= 0, 6
     ! Eq. (4.2.18) p. 27 de Maity (2003), aussi Eq. (18) de Gross & Sadowski (2001) et galement [II-145] et [II-146] p. 51 de Belkadi
     a_i(i)= a(0,i) + (m_bar-1._dp)/m_bar * a(1,i) + (m_bar-1._dp)/m_bar * (m_bar-2._dp)/m_bar * a(2,i) ! V
     b_i(i)= b(0,i) + (m_bar-1._dp)/m_bar * b(1,i) + (m_bar-1._dp)/m_bar * (m_bar-2._dp)/m_bar * b(2,i) ! V
  end do
  
  ! the I1 and I2:
  I_1= 0._dp
  I_2= 0._dp
  do i= 0, 6
     ! Eq. (4.2.16) (qui est fausse cf. (4.1.20) p. 21 et (4.2.15) de Maity (2003); et aussi Eq. (16) & (17) de Gross & Sadowski (2001)
     ! [II-143]et [II-144] p. 51 de Belkadi
     I_1= I_1 + a_i(i) * eta**i
     I_2= I_2 + b_i(i) * eta**i
  end do  
  
  detaI1sdeta= 0._dp
  detaI2sdeta= 0._dp
  do i= 0, 6
     ! Eq. (4.2.29) & (4.2.30) de Maity (2003)
     detaI1sdeta= detaI1sdeta + a_i(i) * (i+1._dp) * eta**i ! V
     detaI2sdeta= detaI2sdeta + b_i(i) * (i+1._dp) * eta**i ! V
  end do
  
  ! Eq. (4.2.31) p. 29 de Maity (2003)
  C_2= -C_1**2 * (m_bar * (-4._dp*eta**2+20._dp*eta+8._dp)/(1._dp-eta)**5 + (1._dp-m_bar) * &
        (2._dp*eta**3+12._dp*eta**2-48._dp*eta+40._dp) &
        / (1._dp-eta)**3 / (2._dp-eta)**3 ) ! V
  
  ! Calcul de <m^2 Epsilon sigma^3> :
  m2_epsilon_sigma3_bar= 0._dp
  do i= 1, size(species_prop)
     do j= 1, size(species_prop)
        ! Eq. (4.2.12) p. 26 de Maity (2003)
        m2_epsilon_sigma3_bar= m2_epsilon_sigma3_bar + species_prop(i)%x * species_prop(j)%x * &
                                                       species_prop(i)%numseg * species_prop(j)%numseg * &
                                                       epsilon(i,j)/T * sigma(i,j)**3 ! V
     end do
  end do
  
  ! Calcul de <m^2 Epsilon^2 sigma^3> :
  m2_epsilon2_sigma3_bar= 0._dp
  do i= 1, size(species_prop)
     do j= 1, size(species_prop)
        ! Eq. (4.2.13) p. 27 de Maity (2003)
        m2_epsilon2_sigma3_bar= m2_epsilon2_sigma3_bar + species_prop(i)%x * species_prop(j)%x * &
                                                         species_prop(i)%numseg * species_prop(j)%numseg * &
                                                         epsilon(i,j)**2/T**2 * sigma(i,j)**3 ! V
     end do
  end do
  
  ! On finalise le calcul de Z^disp:
  ! Eq. (4.2.28) p. 29 de Maity (2003) et (A.28) de G&S :
  
  !print*, ' Dans Zdisp :'
  !print*, ' > rho                   = ', rho 
  !print*, ' > detaI1sdeta           = ', detaI1sdeta
  !print*, ' > m2_epsilon_sigma3_bar = ', m2_epsilon_sigma3_bar
  !print*, ' > C_1                   = ', C_1
  !print*, ' > detaI2sdeta           = ', detaI2sdeta
  !print*, ' > C_2                   = ', C_2
  !print*, ' > eta                   = ', eta
  !print*, ' > I_2                   = ', I_2
  !print*, ' > m2_epsilon2_sigma3_bar= ', m2_epsilon2_sigma3_bar   
  
  Z_disp= -2._dp*pi*rho*detaI1sdeta*m2_epsilon_sigma3_bar - pi*rho*m_bar*(C_1*detaI2sdeta+C_2*eta*I_2)*m2_epsilon2_sigma3_bar ! V
  
  ! Calcul final de Z :
  ! Eq. (4.2.24) p. 28 de Maity (2003) et (A.24) de G&S :
  Z= 1._dp + Z_hc + Z_disp
  ! End of the computation of the compressibility factor "Z":
  !---------------------------------------------------------------------------------------------------------------------------------
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Beginning of the computation of the pressure "Psaft":
  Psaft= (1._dp + Z_hc + Z_disp) * kB * T * rho / angstrom_meter**3.0_dp ! with rho in "number of molecules / angstrom^3", Psaft in pascal.
  
  ! End of the computation of the pressure "Psaft".
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Beginning of the computation of the chemical potential "mu_res/kT" (Eq. 4.2.33 of Maity, 2003):
  !
  a_hs= 1._dp/zeta(0) * ( 3._dp*zeta(1)*zeta(2)/(1._dp-zeta(3))  + &              ! Eq. (A.6) de Gross & Sadowski (2001)
                         zeta(2)**3/zeta(3)/(1._dp-zeta(3))**2 + &              ! et aussi Eq. (4.2.6) de Maity (2003)
                         (zeta(2)**3/zeta(3)**2-zeta(0)) * log(1._dp-zeta(3)) & !
                         ) ! V                                                    !
  
  sum_xi_mim1_ln_gii= 0._dp
  do i= 1, size(species_prop)
     sum_xi_mim1_ln_gii= sum_xi_mim1_ln_gii + species_prop(i)%x *               & ! 
                                              (species_prop(i)%numseg - 1._dp) * & ! 
                                              log(gii_hs(i)) ! V 
  end do
  
  a_hc= m_bar * a_hs - sum_xi_mim1_ln_gii ! cf. Eq. (A.4) p. 1256 de Gross & Sadowski (2001) et aussi Eq. (4.2.4) de Maity (2003) ! V
  
  ! Eq. (A.10) p. 1256 de Gross & Sadowski (2001) et aussi (4.2.10) de Maity (2003):
  a_disp= -2._dp * pi * rho * I_1 * m2_epsilon_sigma3_bar - pi * rho * m_bar * C_1 * I_2 * m2_epsilon2_sigma3_bar ! V
  !write(6,*) rho, I_1, m2_epsilon_sigma3_bar, m_bar, C_1, I_2, m2_epsilon2_sigma3_bar
  
  a_res= a_hc + a_disp ! Eq. (A.3) p. 1256 de G.&S. (2001) et aussi (4.2.3) p. 25 de Maity (2003) ! V
  
  if (present(ares)) then
     ares= R_gas * T * a_res
  end if
  
  do k= 1, size(species_prop)
     zeta_xk(0,k)= pi/6._dp * rho * species_prop(k)%numseg 
     zeta_xk(1,k)= pi/6._dp * rho * species_prop(k)%numseg * d(k)
     zeta_xk(2,k)= pi/6._dp * rho * species_prop(k)%numseg * d(k)**2
     zeta_xk(3,k)= pi/6._dp * rho * species_prop(k)%numseg * d(k)**3
  end do
  
  ! Calcul des (dg^hs(i,i)/dx_k) : cf. Eq. (A.37) de G.&S. (2001) et aussi (4.2.37) p. 30 de Maity (2003).
  do i= 1, size(species_prop)
     do k= 1, size(species_prop)
        dghs_ii_dxk(i,k)= zeta_xk(3,k) / (1._dp - zeta(3))**2 + &
                          d(i)*d(i)/(d(i)+d(i)) *              &
                          ( 3._dp*zeta_xk(2,k)/(1._dp-zeta(3))**2 + 6._dp*zeta(2)*zeta_xk(3,k)/(1._dp-zeta(3))**3) + &
                          ( d(i)*d(i)/(d(i)+d(i)) )**2 *       &
                          ( 4._dp*zeta(2)*zeta_xk(2,k)/(1._dp-zeta(3))**3 + 6._dp*zeta(2)**2*zeta_xk(3,k)/(1._dp-zeta(3))**4 ) ! V
     end do
  end do
  
  ! Calcul de (dahs/dxk), cf. Eq. (A.36) de G.&S. (2001) et aussi (4.2.36) p. 30 de Maity (2001):
  do k= 1, size(species_prop)
     da_hs_s_dxk(k)= - zeta_xk(0,k)/zeta(0) * a_hs + 1._dp/zeta(0) *     &
                         (                                               &
                           3._dp * (zeta_xk(1,k)*zeta(2)+zeta(1)*zeta_xk(2,k)) / (1._dp-zeta(3)) +  &
                           3._dp*zeta(1)*zeta(2)*zeta_xk(3,k) / (1._dp-zeta(3))**2               +  &
                           3._dp*zeta(2)**2*zeta_xk(2,k) / zeta(3) / (1._dp-zeta(3))**2          +  &
                           zeta(2)**3*zeta_xk(3,k)*(3._dp*zeta(3)-1._dp) / zeta(3)**2 / (1._dp-zeta(3))**3 + &
                           ( (3._dp*zeta(2)**2*zeta_xk(2,k)*zeta(3)-2._dp*zeta(2)**3*zeta_xk(3,k))/zeta(3)**3 - zeta_xk(0,k) ) * &
                           log(1._dp - zeta(3)) +                                                  &
                           ( zeta(0) - zeta(2)**3/zeta(3)**2 ) * zeta_xk(3,k) / (1._dp-zeta(3))    &
                         )  ! V
  end do
  
  ! We compute (dahc/dxk), cf. Eq. (A.35) p. 1258 de Gross & Sadowski (2001) and also (4.2.35) p. 30 de Maity (2001):
  do k= 1, size(species_prop)
     sum_xi_mim1_ghsiim1_dghsiisxk= 0._dp
     do i= 1, size(species_prop)
        sum_xi_mim1_ghsiim1_dghsiisxk= sum_xi_mim1_ghsiim1_dghsiisxk + &
                                           species_prop(i)%x * (species_prop(i)%numseg-1._dp) / gii_hs(i) * dghs_ii_dxk(i,k)
     end do
     dahc_dxk(k)= species_prop(k)%numseg * a_hs + &
                  m_bar * da_hs_s_dxk(k) - sum_xi_mim1_ghsiim1_dghsiisxk & ! ATTENTION typo de (4.2.35) de Maity (2003)
                  - (species_prop(k)%numseg-1._dp) * log(gii_hs(k))         ! La formule n'est pas cohrente avec (4.2.4) p. 25 de Maity (2003
                                                                           ! et diffre de (A.35) de G&S (2001)
                  ! ATTENTION : le terme en "- (species_prop(k)%numseg-1.d0) * log(gii_hs(k)) " n'est pas dans Maity (2003) (4.2.35)
                  !             et n'est pas non plus dans G&S (2001) (A.35) alors qu'il est prsent dans le code de Gross et
                  !             qu'on le trouve en drivant (4.2.4) p. 25 de Maity (2003), ((A.4) de G&S 2001).
  end do
  
  !-------------------------------------------
  ! Beginning of computation of (da^disp/dxk):
  
  ! Use of Eq. (A.39) of G&S (2001)
  do k= 1, size(species_prop)
     m2_epsilon_sigma3_bar_xk(k)= 0._dp
     do j= 1, size(species_prop)
        m2_epsilon_sigma3_bar_xk(k)= m2_epsilon_sigma3_bar_xk(k) + species_prop(j)%x * species_prop(j)%numseg * epsilon(k,j)/T * &
                                                                   sigma(k,j)**3 ! V
     end do
     m2_epsilon_sigma3_bar_xk(k)= 2._dp * species_prop(k)%numseg * m2_epsilon_sigma3_bar_xk(k) ! V
  end do

  ! Use of Eq. (A.40) of G&S (2001)
  do k= 1, size(species_prop)
     m2_epsilon2_sigma3_bar_xk(k)= 0._dp
     do j= 1, size(species_prop)
        m2_epsilon2_sigma3_bar_xk(k)= m2_epsilon2_sigma3_bar_xk(k) + species_prop(j)%x * species_prop(j)%numseg * &
                                                                     epsilon(k,j)**2/T**2 * sigma(k,j)**3 ! V
     end do
     m2_epsilon2_sigma3_bar_xk(k)= 2._dp * species_prop(k)%numseg * m2_epsilon2_sigma3_bar_xk(k) ! V
  end do

  ! Use of Eq. (A.41) of G.&S. (2001):
  do k= 1, size(species_prop)
     C_1_xk(k)= C_2 * zeta_xk(3,k) - C_1**2 * (species_prop(k)%numseg * (8._dp*eta-2._dp*eta**2)/(1._dp-eta)**4 - &
                          species_prop(k)%numseg * (20._dp*eta - 27._dp*eta**2 + 12._dp * eta**3 - 2._dp*eta**4 ) &
                                               / (1._dp-eta)**2 / (2._dp-eta)**2 ) ! V
  end do
  
  ! Use of Eq. (A.44) and (A.45) from G&S (2001):
  do i= 0, 6
     do k= 1, size(species_prop)
        ai_xk(i,k)= species_prop(k)%numseg/m_bar**2 * a(1,i) + species_prop(k)%numseg/m_bar**2 * (3._dp - 4._dp/m_bar) * a(2,i) ! V
        bi_xk(i,k)= species_prop(k)%numseg/m_bar**2 * b(1,i) + species_prop(k)%numseg/m_bar**2 * (3._dp - 4._dp/m_bar) * b(2,i) ! V
     end do
  end do
  
  ! Use of Eq. (A.42) and Eq. (A.43) of G&S (2001): 
  do k= 1, size(species_prop)
     I_1_xk(k)= 0._dp
     I_2_xk(k)= 0._dp
     do i= 0, 6
        I_1_xk(k)= I_1_xk(k) + a_i(i) * i * zeta_xk(3,k) * eta**(i-1) + ai_xk(i,k) * eta**i ! V
        I_2_xk(k)= I_2_xk(k) + b_i(i) * i * zeta_xk(3,k) * eta**(i-1) + bi_xk(i,k) * eta**i ! V
     end do
  end do
  
  ! Use Eq. (A.38) of G.&S. (2001):
  do k= 1, size(species_prop)
     dadisp_xk(k)= -2._dp*pi*rho * (I_1_xk(k) * m2_epsilon_sigma3_bar + I_1 * m2_epsilon_sigma3_bar_xk(k)) - &
                    pi * rho * ( (species_prop(k)%numseg*C_1*I_2 + m_bar*C_1_xk(k)*I_2 + m_bar*C_1*I_2_xk(k)) * &
                               m2_epsilon2_sigma3_bar + &
                               m_bar*C_1*I_2*m2_epsilon2_sigma3_bar_xk(k) ) ! V sauf que a^disp dpend de rho (cf. 4.2.10 p. 26)
                               ! et rho dpend de xk !! -> pb ?
    ! TEST TEST TEST en ajoutant le terme de drive de 'rho'
    !dadisp_xk(k)= dadisp_xk(k) - 6.d0/pi*eta*species_prop(k)%numseg * d(k)**3/sum_ximidi3**2 * a_disp/rho
        
  end do
  ! End of computation of (da^disp/dxk):
  !-------------------------------------------

  ! Finally we get the derivatives (da^res/dxk):
  do k= 1, size(species_prop)
     dares_dxk(k)= dahc_dxk(k) + dadisp_xk(k) ! V voir par ex. (4.2.3) de Maity (2001).
  end do
 
  ! Computation of the chemical potential following the Eq. (A.33) of G&S (2001) and also Maity (2001) Eq. (4.2.33) p. 29
  do k= 1, size(species_prop)
     sum_xj_dares_dxj= 0._dp
     do j= 1, size(species_prop)
        sum_xj_dares_dxj= sum_xj_dares_dxj + species_prop(j)%x * dares_dxk(j)
     end do 
     muskT_k(k)= a_res + (Z-1._dp) + dares_dxk(k) - sum_xj_dares_dxj ! V   
     !! write(6,*) '>>>', a_res, (Z-1.d0), dares_dxk(k), sum_xj_dares_dxj
     !!write(6,*) ' > muskT_k(',k,')=', muskT_k(k)
  end do    
  

  !---------------------------------------------------------------------------------------------------------------------------------
  ! EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER 
  ! Computation of association terms:
  !allocate(muskT_k_assoc(1:size(species_prop)),stat=ok);                           call alloc_error(subname,'muskT_k_assoc',ok)
  !
  !call association_terms(subname,1,species_prop,rho,T,Zassoc,muskT_k_assoc)
  !
  !deallocate(muskT_k_assoc)
  ! EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER 
  !---------------------------------------------------------------------------------------------------------------------------------
  
  deallocate(d)
  deallocate(gii_hs)
  deallocate(rho_dgii_hs)
  deallocate(sigma)
  deallocate(epsilon)
  deallocate(zeta_xk)
  deallocate(da_hs_s_dxk)
  deallocate(dghs_ii_dxk)
  deallocate(dahc_dxk)
  deallocate(dares_dxk)
  deallocate(m2_epsilon_sigma3_bar_xk)
  deallocate(m2_epsilon2_sigma3_bar_xk)
  deallocate(C_1_xk)
  deallocate(ai_xk)
  deallocate(bi_xk)
  deallocate(I_1_xk)
  deallocate(I_2_xk)
  deallocate(dadisp_xk)

  return
  
end subroutine P_pcsaft_rho

!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 
!                           -- Computation of the pressure P in the frame of the PC-SAFT theory --
!
!                              >>>    AS A FUNCTION OF "T" and "eta" the reduced density  <<<
!
!-----------------------------------------------------------------------------------------------------------------------------------
! D. Cordier, Institut UTINAM, Besançon, France.
! 1er novembre 2013 -- daniel.cordier@obs-besancon.fr
! 10 mai 2020 -- daniel.cordier@univ-reims.fr
!
! References: * Sunil Kumar Maity, "Modeling and simulation of solid-liquid equilibrium by
!               perturbed-chain statistical associating fluid theory", Master of Technology in Chemical Engineering,
!               2003,
!               Indian Institute of Technology, Kharagpur, India.
!               WARNING: a lot of mistakes/typos in this references! 
!             * Abdelkrim Belkadi, Thèse de l'Universit de Toulouse, 2008.
!             * Gross & Sadowski, Ind. Eng. Chem. Res. 2001, 40, 1244-1260.
!
! - 15 avril 2020 : - passage des constantes numriques de l'écriture '1.d0' vers l'écriture '1._dp'
!                   - mise dans le module 'NEW_MOD_PCSAFT'.
!
subroutine P_pcsaft_eta(sub,ncall,species_prop,k_ij,T,eta,Psaft,rho,Z,muskT_k,ares)
  use utils_DC
  
  implicit none
  character(len=*), parameter :: subname = 'P_pcsaft_eta'
  
  ! Inputs/outputs of the program:
  character(len=*), intent(in) :: sub
  integer, intent(in)          :: ncall
  character(len=len_trim(sub)+len_trim(subname)+1) :: subtemp

  type(compprop), dimension(:), intent(in) :: species_prop ! The properties of studied species.
  real(dp),     dimension(:,:), intent(in) :: k_ij    ! interaction parameters at the given temperature T
  real(dp), intent(in)                     :: T     ! Temperature (in K).
  real(dp), intent(in)                     :: eta   ! The reduced density (no unit)
  real(dp), intent(out)                    :: Psaft ! The pressure (in Pa) computed in the frame of the PC-SAFT
  real(dp), intent(out)                    :: rho   ! the density (kg.m^-3) of the system at P and T
  real(dp), intent(out)                    :: Z     ! The compressibility factor (no unit)
  real(dp), dimension(:), intent(out)      :: muskT_k ! The residual chemical potential divided by kT (no unit)
  real(dp),          intent(out), optional :: ares ! Residual (compare to ideal gas value) Helmoltz energy (J.mol^-1).
  
  logical :: opt_outputs
  
  integer  :: i, j, k, n, ok
  real(dp), save :: m_bar
  real(dp), dimension(:), allocatable :: d
  real(dp), dimension(0:3) :: zeta
  real(dp), dimension(:,:), allocatable :: zeta_xk
  
  real(dp), save :: sum_ximidi0, sum_ximidi1, sum_ximidi2, sum_ximidi3

  real(dp), dimension(:),   allocatable :: gii_hs
  real(dp), dimension(:),   allocatable :: rho_dgii_hs
  real(dp), dimension(:,:), allocatable :: dghs_ii_dxk

  real(dp) :: C_1
  real(dp), dimension(:,:), allocatable :: sigma, epsilon
  real(dp), dimension(0:2,0:6) :: a, b
  real(dp), dimension(0:6) :: a_i, b_i
  real(dp) :: I_1, I_2, detaI1sdeta, detaI2sdeta
  real(dp) :: C_2
  real(dp) :: m2_epsilon_sigma3_bar, m2_epsilon2_sigma3_bar
  real(dp) :: Z_hs, Z_hc, Z_disp, Zassoc
  real(dp) :: molar_mass_bar
  real(dp) :: a_res, a_hc, a_hs, a_disp ! reduced Helmholtz energy (no unit)
  real(dp) :: sum_xi_mim1_ln_gii, sum_xi_mim1_ghsiim1_dghsiisxk, sum_xj_dares_dxj
  real(dp), dimension(:),   allocatable :: da_hs_s_dxk, dahc_dxk, dares_dxk
  real(dp), dimension(:),   allocatable :: m2_epsilon_sigma3_bar_xk, m2_epsilon2_sigma3_bar_xk
  real(dp), dimension(:),   allocatable :: C_1_xk
  real(dp), dimension(:,:), allocatable :: ai_xk, bi_xk
  real(dp), dimension(:),   allocatable :: I_1_xk, I_2_xk, dadisp_xk
  real(dp), dimension(:),   allocatable :: muskT_k_assoc

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Call path of this subroutine:
  subtemp= sub(1:len_trim(sub))//'/'//subname(1:len_trim(subname))
  
  !---------------------------------------------------------------------------------------------------------------------------------
  if ( eta <= 0._dp ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the considered reduced density "eta" is <= 0.!'
     write(6,*) '   We have: "eta" = ', eta
     write(6,*) ''
     stop
  end if
  
  !write(6,*) ' > In "P_pcsaft": ', species_prop(1)%name, species_prop(1)%assoflag, species_prop(1)%assenerg, &
  !                                 species_prop(1)%assvol 
  !stop
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! We compute the mean number of segment "m_bar" and various quantities used in the following:
     m_bar= 0._dp ! cf. Gross & Sadowski (2001)
     do i= 1, size(species_prop)
        m_bar= m_bar + species_prop(i)%x * species_prop(i)%numseg
     end do
     
     ! Computation of the d_i's:
     allocate(d(1:size(species_prop)),stat=ok);                       call alloc_error(subname,'d',ok)
     
     do i= 1, size(species_prop)
        ! Eq. (4.2.9) de Maity (2003) et aussi Eq. [II-136] p. 50 de Belkadi (2008)
        d(i)= species_prop(i)%segdiam * (1._dp - 0.12_dp * exp(-3._dp * species_prop(i)%segenerg / T)) ! En Angstrom
        !write(6,*) ' d(',i,')= ', d(i)
     end do

     sum_ximidi0= 0._dp
     do i= 1, size(species_prop)
        ! cf. Eq. (4.2.8) p. 26 Maity (2003) and Eq. [II-131] p. 49 Belkadi (2008)
        sum_ximidi0= sum_ximidi0 + species_prop(i)%x * species_prop(i)%numseg
     end do
     
     sum_ximidi1= 0._dp
     do i= 1, size(species_prop)
        ! cf. Eq. (4.2.8) p. 26 Maity (2003) and Eq. [II-131] p. 49 Belkadi (2008)
        sum_ximidi1= sum_ximidi1 + species_prop(i)%x * species_prop(i)%numseg * d(i)
     end do
     
     sum_ximidi2= 0._dp
     do i= 1, size(species_prop)
        ! cf. Eq. (4.2.8) p. 26 Maity (2003) and Eq. [II-131] p. 49 Belkadi (2008)
        sum_ximidi2= sum_ximidi2 + species_prop(i)%x * species_prop(i)%numseg * d(i)**2
     end do
     
     sum_ximidi3= 0._dp
     do i= 1, size(species_prop)
        ! cf. Eq. (4.2.20) Maity (2003) and Eq. [II-124] p. 48 Belkadi (2008)
        sum_ximidi3= sum_ximidi3 + species_prop(i)%x * species_prop(i)%numseg * d(i)**3
     end do
  !---------------------------------------------------------------------------------------------------------------------------------
  
  allocate(gii_hs(1:size(species_prop)),stat=ok);                           call alloc_error(subname,'gii_hs',ok)
  allocate(rho_dgii_hs(1:size(species_prop)),stat=ok);                      call alloc_error(subname,'rho_dgii_hs',ok)
  allocate(sigma(1:size(species_prop),1:size(species_prop)),stat=ok);       call alloc_error(subname,'sigma',ok)
  allocate(epsilon(1:size(species_prop),1:size(species_prop)),stat=ok);     call alloc_error(subname,'epsilon',ok)
  allocate(zeta_xk(0:3,1:size(species_prop)),stat=ok);                      call alloc_error(subname,'zeta_xk',ok)
  allocate(da_hs_s_dxk(1:size(species_prop)),stat=ok);                      call alloc_error(subname,'da_hs_s_dxk',ok)
  allocate(dahc_dxk(1:size(species_prop)),stat=ok);                         call alloc_error(subname,'dahc_dxk',ok)
  allocate(dares_dxk(1:size(species_prop)),stat=ok);                        call alloc_error(subname,'dares_dxk',ok)
  allocate(dghs_ii_dxk(1:size(species_prop),1:size(species_prop)),stat=ok); call alloc_error(subname,'dghs_ii_dxk',ok)
  allocate(m2_epsilon_sigma3_bar_xk(1:size(species_prop)),stat=ok);         call alloc_error(subname,'m2_epsilon_sigma3_bar_xk',ok)
  allocate(m2_epsilon2_sigma3_bar_xk(1:size(species_prop)),stat=ok);        call alloc_error(subname,'m2_epsilon2_sigma3_bar_xk',ok)
  allocate(C_1_xk(1:size(species_prop)),stat=ok);                           call alloc_error(subname,'C_1_xk',ok)
  allocate(ai_xk(0:6,1:size(species_prop)),stat=ok);                        call alloc_error(subname,'ai_xk',ok)
  allocate(bi_xk(0:6,1:size(species_prop)),stat=ok);                        call alloc_error(subname,'bi_xk',ok)
  allocate(I_1_xk(1:size(species_prop)),stat=ok);                           call alloc_error(subname,'I_1_xk',ok)
  allocate(I_2_xk(1:size(species_prop)),stat=ok);                           call alloc_error(subname,'I_2_xk',ok)
  allocate(dadisp_xk(1:size(species_prop)),stat=ok);                        call alloc_error(subname,'dadisp_xk',ok)
    
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Beginning of the computation of the compressibility factor "Z":
  rho= 6._dp/pi * eta / sum_ximidi3  ! cf. Eq. (4.2.20) Maity (2003) and Eq. [II-124] p. 48 Belkadi (2008)
                                    ! unit: number of molecules / angstrom^3 ! V
  
  ! Voir Eq. (9) de Gross & Sadowski (2001)
  zeta(0)= pi/6._dp * rho * sum_ximidi0 ! Angstrom^-3   ! V                             
  zeta(1)= pi/6._dp * rho * sum_ximidi1 ! Angstrom^-2   ! V
  zeta(2)= pi/6._dp * rho * sum_ximidi2 ! Angstrom^-1   ! V
  zeta(3)= pi/6._dp * rho * sum_ximidi3 ! no unit       ! V

  ! Eq. (8) de Gross & Sadowski (2001) et ((4.2.26 de MAity (2003) :
  Z_hs= zeta(3)/(1._dp-zeta(3)) + &
        3._dp*zeta(1)*zeta(2)/zeta(0)/(1._dp-zeta(3))**2 + &
        (3._dp*zeta(2)**3-zeta(3)*zeta(2)**3)/zeta(0)/(1._dp-zeta(3))**3 ! V
           
  ! We compute the term "g_ii^hs", cf. for instance (4.2.7) p. 26 of Maity (2003):
  do i= 1, size(species_prop)
     ! ATTENTION: cette formule est (4.2.7) p. 26 de Maity (2003) en accord avec [II-129] p. 49 de Belkadi (2008) mais
     !            dans le papier original de Gross & Sadowski (2001) il n'y a pas de carr au dnominateur du deuxime terme.
     !            en fait il y a de carr dans la formule (A.7) de G&S.
     gii_hs(i)= 1._dp/(1._dp-zeta(3)) + d(i)*d(i)/(d(i)+d(i)) * 3._dp*zeta(2)/(1._dp-zeta(3))**2 + (d(i)*d(i)/(d(i)+d(i)))**2 * &
             2._dp*zeta(2)**2/(1._dp-zeta(3))**3 ! V
  end do  
  
  do i= 1, size(species_prop)
        ! Voir (4.2.27) p. 29 de Maity (2003) et (A.27) de G&S:
        j=i
        rho_dgii_hs(i)= zeta(3)/(1._dp-zeta(3))**2 + d(i)*d(j)/(d(i)+d(j)) * &
                         (3._dp*zeta(2)/(1._dp-zeta(3))**2 + 6._dp*zeta(2)*zeta(3)/(1._dp-zeta(3))**3) + &
                         (d(i)*d(j)/(d(i)+d(j)))**2 * &
                         (4._dp*zeta(2)**2/(1._dp-zeta(3))**3 + 6._dp*zeta(2)**2*zeta(3)/(1._dp-zeta(3))**4) ! V
  end do
  
  ! On calcule Z_hc:
  ! Eq. (4.2.25) p. 28 de Maity (2003) et aussi Eq. (5) de Gross & Sadowski (2001):
  Z_hc= m_bar* Z_hs
  do i= 1, size(species_prop)
     Z_hc= Z_hc - species_prop(i)%x * (species_prop(i)%numseg-1._dp) /gii_hs(i) * rho_dgii_hs(i) ! V
  end do
  
  ! We compute Z^disp:
  ! Eq. (4.2.11) p. 26 de Maity (2003), ATTENTION dans G&S il manque l'exposant <<-1>> !
  C_1=  1._dp/(1._dp + m_bar * (8_dp*eta-2._dp*eta**2)/(1._dp-eta)**4 + (1._dp-m_bar) * &
        (20._dp*eta-27._dp*eta**2+12._dp*eta**3-2._dp*eta**4)/(1._dp-eta)**2/(2._dp-eta)**2) ! V
        
  do i= 1, size(species_prop)
     do j= 1, size(species_prop)
        sigma(i,j)  = 1._dp/2._dp * (species_prop(i)%segdiam + species_prop(j)%segdiam)
        ! ATTENTION 1 : on met une valeur fixe pour k_ij --->>> il faudrait faire une table initiale o on a tous les k_ij
        !               ce que je mets l est destin aux premiers tests.
        ! ATTENTION 2 : - la formule utlise ici est celle de Maity (2003) (4.2.15) p. 27
        !               - Belkadi (2008) en propose une autre (cf. p. 38)
        !               - Gross & Sadowski (2001) Eq. (24) est encore autre chose !

        ! ATTENTION : les calculs des 'k_ij'  la temprature de travail 'T' se fait dans la routine 'pcsaft'
        epsilon(i,j)= sqrt(species_prop(i)%segenerg * species_prop(j)%segenerg) * (1._dp - k_ij(i,j))
        !print*, ' k_ij(i,j)= ', k_ij(i,j)
     end do
  end do
  
  call univ_cst(subname,1,a,b)
  
  do i= 0, 6
     ! Eq. (4.2.18) p. 27 de Maity (2003), aussi Eq. (18) de Gross & Sadowski (2001) et galement [II-145] et [II-146] p. 51 de Belkadi
     a_i(i)= a(0,i) + (m_bar-1._dp)/m_bar * a(1,i) + (m_bar-1._dp)/m_bar * (m_bar-2._dp)/m_bar * a(2,i) ! V
     b_i(i)= b(0,i) + (m_bar-1._dp)/m_bar * b(1,i) + (m_bar-1._dp)/m_bar * (m_bar-2._dp)/m_bar * b(2,i) ! V
  end do
  
  ! the I1 and I2:
  I_1= 0._dp
  I_2= 0._dp
  do i= 0, 6
     ! Eq. (4.2.16) (qui est fausse cf. (4.1.20) p. 21 et (4.2.15) de Maity (2003); et aussi Eq. (16) & (17) de Gross & Sadowski (2001)
     ! [II-143]et [II-144] p. 51 de Belkadi
     I_1= I_1 + a_i(i) * eta**i
     I_2= I_2 + b_i(i) * eta**i
  end do  
  
  detaI1sdeta= 0._dp
  detaI2sdeta= 0._dp
  do i= 0, 6
     ! Eq. (4.2.29) & (4.2.30) de Maity (2003)
     detaI1sdeta= detaI1sdeta + a_i(i) * (i+1._dp) * eta**i ! V
     detaI2sdeta= detaI2sdeta + b_i(i) * (i+1._dp) * eta**i ! V
  end do
  
  ! Eq. (4.2.31) p. 29 de Maity (2003)
  C_2= -C_1**2 * (m_bar * (-4._dp*eta**2+20._dp*eta+8._dp)/(1._dp-eta)**5 + (1._dp-m_bar) * &
        (2._dp*eta**3+12._dp*eta**2-48._dp*eta+40._dp) &
        / (1._dp-eta)**3 / (2._dp-eta)**3 ) ! V
  
  ! Calcul de <m^2 Epsilon sigma^3> :
  m2_epsilon_sigma3_bar= 0._dp
  do i= 1, size(species_prop)
     do j= 1, size(species_prop)
        ! Eq. (4.2.12) p. 26 de Maity (2003)
        m2_epsilon_sigma3_bar= m2_epsilon_sigma3_bar + species_prop(i)%x * species_prop(j)%x * &
                                                       species_prop(i)%numseg * species_prop(j)%numseg * &
                                                       epsilon(i,j)/T * sigma(i,j)**3 ! V
     end do
  end do
  
  ! Calcul de <m^2 Epsilon^2 sigma^3> :
  m2_epsilon2_sigma3_bar= 0._dp
  do i= 1, size(species_prop)
     do j= 1, size(species_prop)
        ! Eq. (4.2.13) p. 27 de Maity (2003)
        m2_epsilon2_sigma3_bar= m2_epsilon2_sigma3_bar + species_prop(i)%x * species_prop(j)%x * &
                                                         species_prop(i)%numseg * species_prop(j)%numseg * &
                                                         epsilon(i,j)**2/T**2 * sigma(i,j)**3 ! V
     end do
  end do
  
  ! On finalise le calcul de Z^disp:
  ! Eq. (4.2.28) p. 29 de Maity (2003) et (A.28) de G&S :
  
  !print*, ' Dans Zdisp :'
  !print*, ' > rho                   = ', rho 
  !print*, ' > detaI1sdeta           = ', detaI1sdeta
  !print*, ' > m2_epsilon_sigma3_bar = ', m2_epsilon_sigma3_bar
  !print*, ' > C_1                   = ', C_1
  !print*, ' > detaI2sdeta           = ', detaI2sdeta
  !print*, ' > C_2                   = ', C_2
  !print*, ' > eta                   = ', eta
  !print*, ' > I_2                   = ', I_2
  !print*, ' > m2_epsilon2_sigma3_bar= ', m2_epsilon2_sigma3_bar   
  
  Z_disp= -2._dp*pi*rho*detaI1sdeta*m2_epsilon_sigma3_bar - pi*rho*m_bar*(C_1*detaI2sdeta+C_2*eta*I_2)*m2_epsilon2_sigma3_bar ! V
  
  ! Calcul final de Z :
  ! Eq. (4.2.24) p. 28 de Maity (2003) et (A.24) de G&S :
  Z= 1._dp + Z_hc + Z_disp
  ! End of the computation of the compressibility factor "Z":
  !---------------------------------------------------------------------------------------------------------------------------------
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Beginning of the computation of the pressure "Psaft":
  Psaft= (1._dp + Z_hc + Z_disp) * kB * T * rho / angstrom_meter**3.0_dp ! with rho in "number of molecules / angstrom^3", Psaft in pascal.
  
  ! End of the computation of the pressure "Psaft".
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Beginning of the computation of the chemical potential "mu_res/kT" (Eq. 4.2.33 of Maity, 2003):
  !
  a_hs= 1._dp/zeta(0) * ( 3._dp*zeta(1)*zeta(2)/(1._dp-zeta(3))  + &              ! Eq. (A.6) de Gross & Sadowski (2001)
                         zeta(2)**3/zeta(3)/(1._dp-zeta(3))**2   + &              ! et aussi Eq. (4.2.6) de Maity (2003)
                         (zeta(2)**3/zeta(3)**2-zeta(0)) * log(1._dp-zeta(3)) & !
                         ) ! V                                                    !
  
  sum_xi_mim1_ln_gii= 0._dp
  do i= 1, size(species_prop)
     sum_xi_mim1_ln_gii= sum_xi_mim1_ln_gii + species_prop(i)%x *                & ! 
                                              (species_prop(i)%numseg - 1._dp) * & ! 
                                              log(gii_hs(i)) ! V 
  end do
  
  a_hc= m_bar * a_hs - sum_xi_mim1_ln_gii ! cf. Eq. (A.4) p. 1256 de Gross & Sadowski (2001) et aussi Eq. (4.2.4) de Maity (2003) ! V
  
  ! Eq. (A.10) p. 1256 de Gross & Sadowski (2001) et aussi (4.2.10) de Maity (2003):
  a_disp= -2._dp * pi * rho * I_1 * m2_epsilon_sigma3_bar - pi * rho * m_bar * C_1 * I_2 * m2_epsilon2_sigma3_bar ! V
  !write(6,*) rho, I_1, m2_epsilon_sigma3_bar, m_bar, C_1, I_2, m2_epsilon2_sigma3_bar
  
  a_res= a_hc + a_disp ! Eq. (A.3) p. 1256 de G.&S. (2001) et aussi (4.2.3) p. 25 de Maity (2003) ! V
  
  if (present(ares)) then
     ares= R_gas * T * a_res
  end if
  
  do k= 1, size(species_prop)
     zeta_xk(0,k)= pi/6._dp * rho * species_prop(k)%numseg 
     zeta_xk(1,k)= pi/6._dp * rho * species_prop(k)%numseg * d(k)
     zeta_xk(2,k)= pi/6._dp * rho * species_prop(k)%numseg * d(k)**2
     zeta_xk(3,k)= pi/6._dp * rho * species_prop(k)%numseg * d(k)**3
  end do
  
  ! Calcul des (dg^hs(i,i)/dx_k) : cf. Eq. (A.37) de G.&S. (2001) et aussi (4.2.37) p. 30 de Maity (2003).
  do i= 1, size(species_prop)
     do k= 1, size(species_prop)
        dghs_ii_dxk(i,k)= zeta_xk(3,k) / (1._dp - zeta(3))**2 + &
                          d(i)*d(i)/(d(i)+d(i)) *              &
                          ( 3._dp*zeta_xk(2,k)/(1._dp-zeta(3))**2 + 6._dp*zeta(2)*zeta_xk(3,k)/(1._dp-zeta(3))**3) + &
                          ( d(i)*d(i)/(d(i)+d(i)) )**2 *       &
                          ( 4._dp*zeta(2)*zeta_xk(2,k)/(1._dp-zeta(3))**3 + 6._dp*zeta(2)**2*zeta_xk(3,k)/(1._dp-zeta(3))**4 ) ! V
     end do
  end do
  
  ! Calcul de (dahs/dxk), cf. Eq. (A.36) de G.&S. (2001) et aussi (4.2.36) p. 30 de Maity (2001):
  do k= 1, size(species_prop)
     da_hs_s_dxk(k)= - zeta_xk(0,k)/zeta(0) * a_hs + 1._dp/zeta(0) *     &
                         (                                              &
                           3._dp * (zeta_xk(1,k)*zeta(2)+zeta(1)*zeta_xk(2,k)) / (1._dp-zeta(3)) +  &
                           3._dp*zeta(1)*zeta(2)*zeta_xk(3,k) / (1._dp-zeta(3))**2               +  &
                           3._dp*zeta(2)**2*zeta_xk(2,k) / zeta(3) / (1._dp-zeta(3))**2          +  &
                           zeta(2)**3*zeta_xk(3,k)*(3._dp*zeta(3)-1._dp) / zeta(3)**2 / (1._dp-zeta(3))**3 + &
                           ( (3._dp*zeta(2)**2*zeta_xk(2,k)*zeta(3)-2._dp*zeta(2)**3*zeta_xk(3,k))/zeta(3)**3 - zeta_xk(0,k) ) * &
                           log(1._dp - zeta(3)) +                                                  &
                           ( zeta(0) - zeta(2)**3/zeta(3)**2 ) * zeta_xk(3,k) / (1._dp-zeta(3))    &
                         )  ! V
  end do
  
  ! We compute (dahc/dxk), cf. Eq. (A.35) p. 1258 de Gross & Sadowski (2001) and also (4.2.35) p. 30 de Maity (2001):
  do k= 1, size(species_prop)
     sum_xi_mim1_ghsiim1_dghsiisxk= 0._dp
     do i= 1, size(species_prop)
        sum_xi_mim1_ghsiim1_dghsiisxk= sum_xi_mim1_ghsiim1_dghsiisxk + &
                                           species_prop(i)%x * (species_prop(i)%numseg-1._dp) / gii_hs(i) * dghs_ii_dxk(i,k)
     end do
     dahc_dxk(k)= species_prop(k)%numseg * a_hs + &
                  m_bar * da_hs_s_dxk(k) - sum_xi_mim1_ghsiim1_dghsiisxk & ! ATTENTION typo de (4.2.35) de Maity (2003)
                  - (species_prop(k)%numseg-1._dp) * log(gii_hs(k))         ! La formule n'est pas cohrente avec (4.2.4) p. 25 de Maity (2003
                                                                           ! et diffre de (A.35) de G&S (2001)
                  ! ATTENTION : le terme en "- (species_prop(k)%numseg-1.d0) * log(gii_hs(k)) " n'est pas dans Maity (2003) (4.2.35)
                  !             et n'est pas non plus dans G&S (2001) (A.35) alors qu'il est prsent dans le code de Gross et
                  !             qu'on le trouve en drivant (4.2.4) p. 25 de Maity (2003), ((A.4) de G&S 2001).
  end do
  
  !-------------------------------------------
  ! Beginning of computation of (da^disp/dxk):
  
  ! Use of Eq. (A.39) of G&S (2001)
  do k= 1, size(species_prop)
     m2_epsilon_sigma3_bar_xk(k)= 0._dp
     do j= 1, size(species_prop)
        m2_epsilon_sigma3_bar_xk(k)= m2_epsilon_sigma3_bar_xk(k) + species_prop(j)%x * species_prop(j)%numseg * epsilon(k,j)/T * &
                                                                   sigma(k,j)**3 ! V
     end do
     m2_epsilon_sigma3_bar_xk(k)= 2._dp * species_prop(k)%numseg * m2_epsilon_sigma3_bar_xk(k) ! V
  end do

  ! Use of Eq. (A.40) of G&S (2001)
  do k= 1, size(species_prop)
     m2_epsilon2_sigma3_bar_xk(k)= 0._dp
     do j= 1, size(species_prop)
        m2_epsilon2_sigma3_bar_xk(k)= m2_epsilon2_sigma3_bar_xk(k) + species_prop(j)%x * species_prop(j)%numseg * &
                                                                     epsilon(k,j)**2/T**2 * sigma(k,j)**3 ! V
     end do
     m2_epsilon2_sigma3_bar_xk(k)= 2._dp * species_prop(k)%numseg * m2_epsilon2_sigma3_bar_xk(k) ! V
  end do

  ! Use of Eq. (A.41) of G.&S. (2001):
  do k= 1, size(species_prop)
     C_1_xk(k)= C_2 * zeta_xk(3,k) - C_1**2 * (species_prop(k)%numseg * (8._dp*eta-2._dp*eta**2)/(1._dp-eta)**4 - &
                          species_prop(k)%numseg * (20._dp*eta - 27._dp*eta**2 + 12._dp * eta**3 - 2._dp*eta**4 ) &
                                               / (1._dp-eta)**2 / (2._dp-eta)**2 ) ! V
  end do
  
  ! Use of Eq. (A.44) and (A.45) from G&S (2001):
  do i= 0, 6
     do k= 1, size(species_prop)
        ai_xk(i,k)= species_prop(k)%numseg/m_bar**2 * a(1,i) + species_prop(k)%numseg/m_bar**2 * (3._dp - 4._dp/m_bar) * a(2,i) ! V
        bi_xk(i,k)= species_prop(k)%numseg/m_bar**2 * b(1,i) + species_prop(k)%numseg/m_bar**2 * (3._dp - 4._dp/m_bar) * b(2,i) ! V
     end do
  end do
  
  ! Use of Eq. (A.42) and Eq. (A.43) of G&S (2001): 
  do k= 1, size(species_prop)
     I_1_xk(k)= 0._dp
     I_2_xk(k)= 0._dp
     do i= 0, 6
        I_1_xk(k)= I_1_xk(k) + a_i(i) * i * zeta_xk(3,k) * eta**(i-1) + ai_xk(i,k) * eta**i ! V
        I_2_xk(k)= I_2_xk(k) + b_i(i) * i * zeta_xk(3,k) * eta**(i-1) + bi_xk(i,k) * eta**i ! V
     end do
  end do
  
  ! Use Eq. (A.38) of G.&S. (2001):
  do k= 1, size(species_prop)
     dadisp_xk(k)= -2._dp*pi*rho * (I_1_xk(k) * m2_epsilon_sigma3_bar + I_1 * m2_epsilon_sigma3_bar_xk(k)) - &
                    pi * rho * ( (species_prop(k)%numseg*C_1*I_2 + m_bar*C_1_xk(k)*I_2 + m_bar*C_1*I_2_xk(k)) * &
                               m2_epsilon2_sigma3_bar + &
                               m_bar*C_1*I_2*m2_epsilon2_sigma3_bar_xk(k) ) ! V sauf que a^disp dpend de rho (cf. 4.2.10 p. 26)
                               ! et rho dpend de xk !! -> pb ?
    ! TEST TEST TEST en ajoutant le terme de drive de 'rho'
    !dadisp_xk(k)= dadisp_xk(k) - 6.d0/pi*eta*species_prop(k)%numseg * d(k)**3/sum_ximidi3**2 * a_disp/rho
        
  end do
  ! End of computation of (da^disp/dxk):
  !-------------------------------------------

  ! Finally we get the derivatives (da^res/dxk):
  do k= 1, size(species_prop)
     dares_dxk(k)= dahc_dxk(k) + dadisp_xk(k) ! V voir par ex. (4.2.3) de Maity (2001).
  end do
 
  ! Computation of the chemical potential following the Eq. (A.33) of G&S (2001) and also Maity (2001) Eq. (4.2.33) p. 29
  do k= 1, size(species_prop)
     sum_xj_dares_dxj= 0._dp
     do j= 1, size(species_prop)
        sum_xj_dares_dxj= sum_xj_dares_dxj + species_prop(j)%x * dares_dxk(j)
     end do 
     muskT_k(k)= a_res + (Z-1._dp) + dares_dxk(k) - sum_xj_dares_dxj ! V   
     !! write(6,*) '>>>', a_res, (Z-1.d0), dares_dxk(k), sum_xj_dares_dxj
     !!write(6,*) ' > muskT_k(',k,')=', muskT_k(k)
  end do    
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Computation of the density "rho" in MKSA:

  molar_mass_bar= 0._dp
  do i= 1, size(species_prop)
     molar_mass_bar= molar_mass_bar + species_prop(i)%x * species_prop(i)%molmass ! Mean molar mass in kg.mol^-1
  end do

  rho= rho * molar_mass_bar / Navo / angstrom_meter**3.0_dp ! rho in kg.m^-3

  !---------------------------------------------------------------------------------------------------------------------------------
  ! EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER 
  ! Computation of association terms:
  !allocate(muskT_k_assoc(1:size(species_prop)),stat=ok);                           call alloc_error(subname,'muskT_k_assoc',ok)
  !
  !call association_terms(subname,1,species_prop,rho,T,Zassoc,muskT_k_assoc)
  !
  !deallocate(muskT_k_assoc)
  ! EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER -- EN CHANTIER 
  !---------------------------------------------------------------------------------------------------------------------------------
  
  deallocate(d)
  deallocate(gii_hs)
  deallocate(rho_dgii_hs)
  deallocate(sigma)
  deallocate(epsilon)
  deallocate(zeta_xk)
  deallocate(da_hs_s_dxk)
  deallocate(dghs_ii_dxk)
  deallocate(dahc_dxk)
  deallocate(dares_dxk)
  deallocate(m2_epsilon_sigma3_bar_xk)
  deallocate(m2_epsilon2_sigma3_bar_xk)
  deallocate(C_1_xk)
  deallocate(ai_xk)
  deallocate(bi_xk)
  deallocate(I_1_xk)
  deallocate(I_2_xk)
  deallocate(dadisp_xk)

  return
  
end subroutine P_pcsaft_eta

!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 
!                       -- Computation of activity coefficient and density following the PC-SAFT theory --
!                                *** Independant variables: rho (density) and T (temperature) ***
!
!-----------------------------------------------------------------------------------------------------------------------------------
! D. Cordier, Institut UTINAM, Besançon, France.
! June 2sd 2014 -- daniel.cordier@obs-besancon.fr
! D. Cordier, CNRS, GSMA, Reims, France.
! 25/02/2019.
!
! - The independant variables are rho (density) and T (temperature).
!
!   25/02/2019 : on modifie le contrle de la normalisation des fractions molaires en tenant compte de la prcision  laquelle le 
!                calcul est fait, on utilise :
!                if (sum > 1._dp + epsilon(sum)) then
!
! - 15 avril 2020 : - passage des constantes numriques de l'criture '1.d0' vers l'criture '1._dp'
!                   - mise dans le module 'NEW_MOD_PCSAFT'.
!
subroutine pcsaft_rhoT(sub,ncall,state,compound,x,rho_MKSA,T,P,muskT_k,ln_phi_k,a_res)
  use utils_DC
  
  implicit none
  character(len=*), parameter :: subname = 'pcsaft_rhoT'
  
  ! Inputs/outputs of the program:
  character(len=*), intent(in) :: sub
  integer, intent(in)          :: ncall
  character(len=len_trim(sub)+len_trim(subname)+1) :: subtemp

  character(len=1),                intent(in) :: state    ! the physical state of the phase:  - state='L' = liquid phase
                                                          !                                   - state='V' = vapor phase.
  character(len=comp_name_ncmax), dimension(:), intent(in) :: compound ! List of the name of compound taken into account.
  real(dp),          dimension(:), intent(in) :: x        ! mole fractions of the various species taken into account.
  real(dp),                        intent(in) :: rho_MKSA ! the density (in kg.m^-3).
  real(dp),                        intent(in) :: T        ! the temperature (in K).
  
  real(dp),                        intent(out) :: P        ! the pressure (in Pa).
  real(dp), dimension(:),          intent(out) :: muskT_k  ! The residual chemical potential divided by kT (no unit)
  real(dp), dimension(:),          intent(out) :: ln_phi_k ! the activity coefficients in neperian log (no unit).
  real(dp),              intent(out), optional :: a_res    ! Residual (compare to ideal gas value) Helmoltz energy (J.mol^-1).
  
  type(compdata),      dimension(:), allocatable, save :: compdatabase
  type(compdataINTER), dimension(:), allocatable, save :: compdatabase_interac
  
  integer :: i, j, k, ok
  logical :: here, trouve
  logical :: first = .true.
  
  character(len=6) :: temp_string
  character(len=comp_name_ncmax) :: temp_string2
  character(len=comp_name_ncmax) :: name1, name2, name3, name4
  real(dp) :: sum, Z, eta
  real(dp), dimension(:,:), allocatable :: k_ij
  
  type(compprop), dimension(:), allocatable :: species_prop

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Call path of this subroutine:
  subtemp= sub(1:len_trim(sub))//'/'//subname(1:len_trim(subname))
  
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! We check various obvious things:
  if ((state /= 'L') .AND. (state /= 'V')) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the considered physical state is neither Liquid nor Vapor!'
     write(6,*) '   We have: "state" = ', state
     write(6,*) ''
     stop
  end if
  
  if ( size(compound) /= size(x) ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the array "compound" (names of considered species) and "x" (their mole fractions) have different sizes!'
     write(6,*) '   Size of "compound" ---: ', size(compound)
     write(6,*) '   Size of "x" ----------: ', size(x)
     write(6,*) ''
     stop
  end if

  if ( size(ln_phi_k) /= size(x) ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the array "ln_phi_k" (names of considered species) and "x" (their mole fractions) have different sizes!'
     write(6,*) '   Size of "ln_phi_k" ------: ', size(ln_phi_k)
     write(6,*) '   Size of "x" -------------: ', size(x)
     write(6,*) ''
     stop
  end if

  if ((rho_MKSA < 0._dp) .OR. (T < 0._dp)) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the density and/or the temperature is/are negative!'
     write(6,*) '   rho = ', rho_MKSA
     write(6,*) '   T   = ', T
     write(6,*) ''
     stop
  end if

  sum= 0._dp
  do i=1, size(x)
     sum= sum + x(i)
  end do
  if ( (sum > 1._dp + epsilon(sum)) .AND. normized_xi_sum ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the sum of mole fractions is larger than 1!'
     write(6,*) '   Sum x(i) = ', sum
     write(6,*) ''
  !   stop
  end if    
  
  ! -------------------------------------------------------------------------------------------------------------------------------- 
  !
  allocate(species_prop(1:size(compound)),stat=ok);            call alloc_error(subname,'species_prop',ok)
  allocate(k_ij(1:size(compound),1:size(compound)),stat=ok);   call alloc_error(subname,'k_ij',ok)
  
  ! -------------------------------------------------------------------------------------------------------------------------------- 
  ! We read the database of PC-SAFT parameters:
  if (first) then
     call read_database_pcsaft(subname,1,compdatabase,compdatabase_interac)
     first= .false.
     write(6,*) ''
     write(temp_string,'(I4)') size(compdatabase)
     call wf(' > Size of the database: ', 'bright red', &
             temp_string(1:len_trim(temp_string)), 'bright black', &
             ' species:', 'bright red')
     do i= 1, size(compdatabase)
        temp_string2= compdatabase(i)%name
        call wf('   - ', 'bright red', &
                        temp_string2(1:len_trim(temp_string2)), 'bright black' )
     end do
     write(6,*) ''
     call wf(' > Species taken into account in this computation: ', 'bright red')
     do i= 1, size(compound)
        temp_string2= compound(i)
        call wf('   - ', 'bright red', &
                        temp_string2(1:len_trim(temp_string2)), 'bright blue' )
     end do
     write(6,*) ''
  end if
  
  !do j= 1, size(compdatabase)
  !   write(6,'(A13,ES17.5,F12.4,F14.4,F15.4,A100)') compdatabase(j)%name,     compdatabase(j)%molmass, &
  !                                                  compdatabase(j)%numseg,   compdatabase(j)%segdiam, &
  !                                                  compdatabase(j)%segenerg, compdatabase(j)%ref
  !end do
  
  !do i= 1, size(x)
  !   write(6,*) '|',trim(ADJUSTL(compound(i))),'|'
  !end do
  
  ! Construction du tableau des paramtres d'interaction : -----------------------------
  do i= 1, size(compound)
     k_ij(i,i)= 0._dp ! On met des zros sur la diagonale
  end do
  
  do i= 1, size(compound)
     name1  = trim(ADJUSTL(compound(i)))
     do j= i+1, size(compound)
        name2  = trim(ADJUSTL(compound(j)))
        trouve = .false.
        do k= 1, size(compdatabase_interac)
           name3  = trim(ADJUSTL(compdatabase_interac(k)%name1))
           name4  = trim(ADJUSTL(compdatabase_interac(k)%name2))
           if (  ( (name1(1:len_trim(name1)) == name3(1:len_trim(name3))) .AND. &
                   (name2(1:len_trim(name2)) == name4(1:len_trim(name4)))   ) .OR. &
                 ( (name2(1:len_trim(name2)) == name3(1:len_trim(name3))) .AND. &
                   (name1(1:len_trim(name1)) == name4(1:len_trim(name4)))   )           ) then
              trouve = .true.
              !print*, compdatabase_interac(k)%a
              !print*, compdatabase_interac(k)%b
              !stop
              k_ij(i,j)= compdatabase_interac(k)%a + compdatabase_interac(k)%b * 1.d-4 * T 
              k_ij(j,i)= k_ij(i,j)
           end if
        end do
        if ( .NOT. trouve) then
           write(6,*) ''
           write(6,*) ' -----------------------------------------------------------------------------------------------------------'
           call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
           write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', sub(1:len_trim(sub)), ''''
           write(6,*) '   we did not find interaction parameter "k_ij" for the couple:'
           write(6,*) '   - ', name1
           write(6,*) '   - ', name2
           write(6,*) '   Please check the database file called: "', datafilename_inter(1:len_trim(datafilename_inter)), '" '
           write(6,*) ' -----------------------------------------------------------------------------------------------------------'
           write(6,*) ''
           stop
        end if
     end do
  end do
  ! Fin de la construction du tableau des paramtres d'interaction : -------------------
  
  do i= 1, size(x)
     name1  = trim(ADJUSTL(compound(i))) ! Enlve les blancs en dbut et fin de chane.
     here= .false.
     do j= 1, size(compdatabase)
        name2= trim(ADJUSTL(compdatabase(j)%name))
        !write(6,*) name2
        if ( name1(1:len_trim(name1)) == name2(1:len_trim(name2)) ) then
           here= .true.
           species_prop(i)%name     =  name1 ! Name
           species_prop(i)%x        =  x(i)  ! Mole fraction
           species_prop(i)%molmass  =  compdatabase(j)%molmass  ! The molar mass (in kg.mol^-1) for the considered compound.
           species_prop(i)%numseg   =  compdatabase(j)%numseg   ! The number of segment for the considered compound.
           species_prop(i)%segdiam  =  compdatabase(j)%segdiam  ! The segment diameter.
           species_prop(i)%segenerg =  compdatabase(j)%segenerg ! The segment energy  
           
           species_prop(i)%assoflag =  compdatabase(j)%assoflag ! Flag: 0=non associative species, 1= associative species
           species_prop(i)%assenerg =  compdatabase(j)%assenerg ! The association energy
           species_prop(i)%assvol   =  compdatabase(j)%assvol   ! The association volume
           exit
        end if
     end do
     if ( .NOT. here ) then
        write(6,*) ''
        call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
        write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', sub(1:len_trim(sub)), ''''
        write(6,*) '   the species "', name1(1:len_trim(name1)), '" is not present in the database!'
        write(6,*) ''
        stop
     end if
  end do
  
  !do i= 1, size(x)
  !   write(6,'(A13,ES15.5,ES17.5,F12.4,F14.4,F15.4)') species_prop(i)%name, species_prop(i)%x, species_prop(i)%molmass, &
  !                                                    species_prop(i)%numseg, species_prop(i)%segdiam,                  &
  !                                                    species_prop(i)%segenerg
  !end do

  ! -------------------------------------------------------------------------------------------------------------------------------- 
  ! In this subroutine, no need to go through a 'Zandmu_pcsaft' like subroutine, we call directly a 'P_pcsaft_eta' like
  ! subroutine:
  !write(6,*) ' > Rho_MKSA= ', rho_MKSA
  
  if(PRESENT(a_res))then
     call P_pcsaft_rho(subtemp,1,species_prop,k_ij,T,rho_MKSA,P,Z,muskT_k,a_res)
  else
     call P_pcsaft_rho(subtemp,2,species_prop,k_ij,T,rho_MKSA,P,Z,muskT_k)
  end if

  ! We then derive the fugacity coefficients "phi_k":
  !write(6,'(4(A,ES15.5))') 'P=', P, ' T= ', T, ' eta= ', eta, ' rho= ', rho
  do k= 1, size(species_prop)
     ln_phi_k(k)= muskT_k(k) - log(Z)
     !write(6,*) ' ln_phi_k(k)= ', ln_phi_k(k), muskT_k(k), log(Z)
  end do
  
  return
  
end subroutine pcsaft_rhoT

!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 
!                      -- Computation of activity coefficient and density following the PC-SAFT theory --
!                                *** Independant variables: P (pressure) and T (temperature) ***
!
!-----------------------------------------------------------------------------------------------------------------------------------
! D. Cordier, CNRS, Institut UTINAM, Besançon, France.
! 1er novembre 2013 -- daniel.cordier@obs-besancon.fr
! D. Cordier, CNRS, GSMA, Reims, France.
! 25/02/2019.
!
! - The independant variables are P (pressure) and T (temperature).
!
!
! - 2 juin 2014: we change the name: 'pcsaft' -> 'pcsaft_PT'
!   25/02/2019 : on modifie le contrle de la normalisation des fractions molaires en tenant compte de la prcision  laquelle le 
!                calcul est fait, on utilise :
!                if (sum > 1._dp + epsilon(sum)) then
!
! - 15 avril 2020 : - mise dans le module 'NEW_MOD_PCSAFT'.
! -    7 mai 2020 : - introduction de l'argument optional 're_initialize'.
!                   - changement de la variable 'present' en 'present_' pour éviter le conflit avec le mot clé 'preent'
!
subroutine pcsaft_PT(sub, ncall, state, compound, x, P, T, rho, muskT_k, ln_phi_k, re_initialize, hressRT)
  use utils_DC
  
  implicit none
  character(len=*), parameter :: subname = 'pcsaft_PT'
  
  ! Inputs/outputs of the program:
  character(len=*), intent(in) :: sub
  integer, intent(in)          :: ncall
  character(len=len_trim(sub)+len_trim(subname)+1) :: subtemp
  
  character(len=1),                intent(in) :: state    ! the physical state of the phase:  - state='L' = liquid phase
                                                          !                                   - state='V' = vapor phase.
  character(len=comp_name_ncmax),  dimension(:), intent(in) :: compound ! List of the name of compound taken into account.
  real(dp),           dimension(:), intent(in) :: x        ! mole fractions of the various species taken into account.
  real(dp),                         intent(in) :: P        ! the pressure (in Pa).
  real(dp),                         intent(in) :: T        ! the temperature (in K).
  logical,                optional, intent(in) :: re_initialize ! To re-initialize the routine.
  
  real(dp),                         intent(out) :: rho      ! the density of the mixture (kg.m^-3).
  real(dp), dimension(:),           intent(out) :: muskT_k  ! The residual chemical potential divided by kT (no unit)
  real(dp), dimension(:),           intent(out) :: ln_phi_k ! the fugacity coefficients in neperian log (no unit).
  real(dp),               optional, intent(out) :: hressRT  ! the residual molar enthalpy /RT
  
  type(compdata),      dimension(:), allocatable, save :: compdatabase
  type(compdataINTER), dimension(:), allocatable, save :: compdatabase_interac
  
  integer  :: i, j, k, ok
  logical  :: present_, trouve
  logical :: first = .true.
  
  character(len=6) :: temp_string
  character(len=comp_name_ncmax) :: temp_string2
  character(len=comp_name_ncmax) :: name1, name2, name3, name4
  real(dp) :: sum, Z, eta
  real(dp), dimension(:,:), allocatable :: k_ij
  
  type(compprop), dimension(:), allocatable :: species_prop
  real(dp), parameter :: pd = 1.E-4_dp ! For the numerical estimation of derivatives.
  real(dp) :: T1, T2, ares1, ares2, daressdT
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Call path of this subroutine:
  subtemp= sub(1:len_trim(sub))//'/'//subname(1:len_trim(subname))
  
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! We check various obvious things:
  if ((state /= 'L') .AND. (state /= 'V')) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the considered physical state is neither Liquid nor Vapor!'
     write(6,*) '   We have: "state" = ', state
     write(6,*) ''
     stop
  end if
  
  if ( size(compound) /= size(x) ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the array "compound" (names of considered species) and "x" (their mole fractions) have different sizes!'
     write(6,*) '   Size of "compound" ---: ', size(compound)
     write(6,*) '   Size of "x" ----------: ', size(x)
     write(6,*) ''
     stop
  end if

  if ( size(ln_phi_k) /= size(x) ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the array "ln_phi_k" (names of considered species) and "x" (their mole fractions) have different sizes!'
     write(6,*) '   Size of "ln_phi_k" ------: ', size(ln_phi_k)
     write(6,*) '   Size of "x" -------------: ', size(x)
     write(6,*) ''
     stop
  end if

  if ((P < 0._dp) .OR. (T < 0._dp)) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the pressure and/or the temperature is/are negative!'
     write(6,*) '   P = ', P
     write(6,*) '   T = ', T
     write(6,*) ''
     stop
  end if

  sum= 0._dp
  do i= 1, size(x)
     sum= sum + x(i)
  end do
  if ( (sum > 1._dp + epsilon(sum)) .AND. normized_xi_sum ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the sum of mole fractions is larger than 1!'
     write(6,*) '   Sum x(i) = ', sum
     write(6,*) '   ncall    = ', ncall
     write(6,*) ''
!     stop
  end if  

  ! -------------------------------------------------------------------------------------------------------------------------------- 
  !
  allocate(species_prop(1:size(compound)),stat=ok);            call alloc_error(subname,'species_prop',ok)
  allocate(k_ij(1:size(compound),1:size(compound)),stat=ok);   call alloc_error(subname,'k_ij',ok)
  
  ! -------------------------------------------------------------------------------------------------------------------------------- 
  ! We read the database of PC-SAFT parameters:
  if (present(re_initialize)) then
     if (re_initialize) then
        first = re_initialize
     end if
  end if
  if (first) then
     call read_database_pcsaft(subname,1,compdatabase,compdatabase_interac)
     first= .false.
     write(6,*) ''
     write(temp_string,'(I4)') size(compdatabase)
     call wf(' > Size of the database: ', 'bright red', &
                        temp_string(1:len_trim(temp_string)), 'bright black', &
                        ' species:', 'bright red')
     do i= 1, size(compdatabase)
        temp_string2= compdatabase(i)%name
        call wf('   - ', 'bright red', &
                        temp_string2(1:len_trim(temp_string2)), 'bright black' )
     end do
     write(6,*) ''
     call wf(' > Species taken into account in this computation: ', 'bright red')
     do i= 1, size(compound)
        temp_string2= compound(i)
        call wf('   - ', 'bright red', &
                        temp_string2(1:len_trim(temp_string2)), 'bright blue' )
     end do
     write(6,*) ''
  end if
  
  !do j= 1, size(compdatabase)
  !   write(6,'(A13,ES17.5,F12.4,F14.4,F15.4,A100)') compdatabase(j)%name,     compdatabase(j)%molmass, &
  !                                                  compdatabase(j)%numseg,   compdatabase(j)%segdiam, &
  !                                                  compdatabase(j)%segenerg, compdatabase(j)%ref
  !end do
  
  !do i= 1, size(x)
  !   write(6,*) '|',trim(ADJUSTL(compound(i))),'|'
  !end do
  
  ! Construction du tableau des paramètres d'interaction : -----------------------------
  do i= 1, size(compound)
     k_ij(i,i)= 0._dp ! On met des zéros sur la diagonale
  end do
  
  do i= 1, size(compound)
     name1  = trim(ADJUSTL(compound(i)))
     do j= i+1, size(compound)
        name2  = trim(ADJUSTL(compound(j)))
        trouve = .false.
        do k= 1, size(compdatabase_interac)
           name3  = trim(ADJUSTL(compdatabase_interac(k)%name1))
           name4  = trim(ADJUSTL(compdatabase_interac(k)%name2))
           if (  ( (name1(1:len_trim(name1)) == name3(1:len_trim(name3))) .AND. &
                   (name2(1:len_trim(name2)) == name4(1:len_trim(name4)))   ) .OR. &
                 ( (name2(1:len_trim(name2)) == name3(1:len_trim(name3))) .AND. &
                   (name1(1:len_trim(name1)) == name4(1:len_trim(name4)))   )           ) then
              trouve = .true.
              !print*, compdatabase_interac(k)%a
              !print*, compdatabase_interac(k)%b
              !stop
              k_ij(i,j)= compdatabase_interac(k)%a + compdatabase_interac(k)%b * 1.d-4 * T 
              k_ij(j,i)= k_ij(i,j)
           end if
        end do
        if ( .NOT. trouve) then
           write(6,*) ''
           write(6,*) ' -----------------------------------------------------------------------------------------------------------'
           call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
           write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
           write(6,*) '   we did not find interaction parameter "k_ij" for the couple:'
           write(6,*) '   - ', name1
           write(6,*) '   - ', name2
           write(6,*) '   Please check the database file called: "', datafilename_inter(1:len_trim(datafilename_inter)), '" '
           write(6,*) ' -----------------------------------------------------------------------------------------------------------'
           write(6,*) ''
           stop
        end if
     end do
  end do
  ! Fin de la construction du tableau des paramètres d'interaction : -------------------
  
  do i= 1, size(x)
     name1   = trim(ADJUSTL(compound(i))) ! Enlève les blancs en début et fin de chaîne.
     present_= .false.
     do j= 1, size(compdatabase)
        name2= trim(ADJUSTL(compdatabase(j)%name))
        !write(6,*) name2
        if ( name1(1:len_trim(name1)) == name2(1:len_trim(name2)) ) then
           present_= .true.
           species_prop(i)%name     =  name1 ! Name
           species_prop(i)%x        =  x(i)  ! Mole fraction
           species_prop(i)%molmass  =  compdatabase(j)%molmass  ! The molar mass (in kg.mol^-1) for the considered compound.
           species_prop(i)%numseg   =  compdatabase(j)%numseg   ! The number of segment for the considered compound.
           species_prop(i)%segdiam  =  compdatabase(j)%segdiam  ! The segment diameter.
           species_prop(i)%segenerg =  compdatabase(j)%segenerg ! The segment energy  
           
           species_prop(i)%assoflag =  compdatabase(j)%assoflag ! Flag: 0=non associative species, 1= associative species
           species_prop(i)%assenerg =  compdatabase(j)%assenerg ! The association energy
           species_prop(i)%assvol   =  compdatabase(j)%assvol   ! The association volume
           exit
        end if
     end do
     if ( .NOT. present_ ) then
        write(6,*) ''
        call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
        write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
        write(6,*) '   the species "', name1(1:len_trim(name1)), '" is not present in the database!'
        write(6,*) ''
        stop
     end if
  end do
  
  !do i= 1, size(x)
  !   write(6,'(A13,ES15.5,ES17.5,F12.4,F14.4,F15.4)') species_prop(i)%name, species_prop(i)%x, species_prop(i)%molmass, &
  !                                                    species_prop(i)%numseg, species_prop(i)%segdiam,                  &
  !                                                    species_prop(i)%segenerg
  !end do

  ! -------------------------------------------------------------------------------------------------------------------------------- 
  ! We compute "Z", the "mu_k/kT"'s and the density "rho" in the frame of the PC-SAFT theory:
  if (present(hressRT)) then
     T1 = T * (1._dp-pd)
     call Zandmu_pcsaft_PT(subtemp,1,state,species_prop,k_ij,P,T1,Z,muskT_k,eta,rho,ares1)
     T2 = T * (1._dp+pd)
     call Zandmu_pcsaft_PT(subtemp,2,state,species_prop,k_ij,P,T2,Z,muskT_k,eta,rho,ares2)

     daressdT = (ares2 - ares1) / (T2 - T1) / R_gas / T !!! ATTENTION : dans '' ares = R*T*ares !!!!
     
     call Zandmu_pcsaft_PT(subtemp,2,state,species_prop,k_ij,P,T,Z,muskT_k,eta,rho)
     hressRT = -T * daressdT + Z - 1._dp ! cf. Eq. 27 p. 871 in Gord et al. (2013)
                                         !     also Eq. 4.2.46 p. 32 in Kumar Maity (MSc Report).
  else
     call Zandmu_pcsaft_PT(subtemp,2,state,species_prop,k_ij,P,T,Z,muskT_k,eta,rho)
  end if
  
  ! We then derive the fugacity coefficients "phi_k":
  !write(6,'(4(A,ES15.5))') 'P=', P, ' T= ', T, ' eta= ', eta, ' rho= ', rho
  do k= 1, size(species_prop)
     ln_phi_k(k)= muskT_k(k) - log(Z)
     !write(6,*) ' ln_phi_k(k)= ', ln_phi_k(k), muskT_k(k), log(Z)
  end do
  
  return
  
end subroutine pcsaft_PT

!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 
!                      -- Computation of activity coefficient and density following the PC-SAFT theory --
!                                *** Independant variables: P (pressure) and rho (density) ***
!
!-----------------------------------------------------------------------------------------------------------------------------------
! D. Cordier, Institut UTINAM, Besançon, France.
! June 2sd -- daniel.cordier@obs-besancon.fr
! D. Cordier, CNRS, GSMA, Reims, France.
! 25/02/2019.
!
! - The independant variables are rho (density) and T (temperature).
!
! - 25/02/2019 : on modifie le contrle de la normalisation des fractions molaires en tenant compte de la prcision  laquelle le 
!                calcul est fait, on utilise :
!                if (sum > 1._dp + epsilon(sum)) then
!
! - 15 avril 2020 : - passage des constantes numriques de l'criture '1.d0' vers l'criture '1._dp'
!                   - mise dans le module 'NEW_MOD_PCSAFT'.
!
subroutine pcsaft_Prho(sub,ncall,state,compound,x,P,rho,T,muskT_k,ln_phi_k,ares)
  use utils_DC
  
  implicit none
  character(len=*), parameter :: subname = 'pcsaft_Prho'
  
  ! Inputs/outputs of the program:
  character(len=*), intent(in) :: sub
  integer, intent(in)          :: ncall
  character(len=len_trim(sub)+len_trim(subname)+1) :: subtemp
  
  character(len=1),                intent(in) :: state    ! the physical state of the phase:  - state='L' = liquid phase
                                                          !                                   - state='V' = vapor phase.
  character(len=comp_name_ncmax), dimension(:), intent(in) :: compound ! List of the name of compound taken into account.
  real(dp),          dimension(:), intent(in) :: x        ! mole fractions of the various species taken into account.
  real(dp),                        intent(in) :: P        ! the pressure (in Pa).
  real(dp),                        intent(in) :: rho      ! the density (in kg.m^-3).
  
  real(dp),                        intent(out) :: T       ! the temperature (in K).
  real(dp), dimension(:),          intent(out) :: muskT_k ! The residual chemical potential divided by kT (no unit)
  real(dp), dimension(:),          intent(out) :: ln_phi_k ! the activity coefficients in neperian log (no unit).
  real(dp),               intent(out), optional :: ares    ! Residual (compare to ideal gas value) Helmoltz energy (J.mol^-1).

  type(compdata),      dimension(:), allocatable, save :: compdatabase
  type(compdataINTER), dimension(:), allocatable, save :: compdatabase_interac
  
  integer  :: i, j, k, ok
  logical  :: here, trouve
  logical :: first = .true.
  
  character(len=6) :: temp_string
  character(len=comp_name_ncmax) :: temp_string2
  character(len=comp_name_ncmax) :: name1, name2, name3, name4
  real(dp) :: sum, Z, eta
  real(dp), dimension(:,:), allocatable :: k_ij
  
  type(compprop), dimension(:), allocatable :: species_prop

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Call path of this subroutine:
  subtemp= sub(1:len_trim(sub))//'/'//subname(1:len_trim(subname))
  
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! We check various obvious things:
  if ((state /= 'L') .AND. (state /= 'V')) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the considered physical state is neither Liquid nor Vapor!'
     write(6,*) '   We have: "state" = ', state
     write(6,*) ''
     stop
  end if
  
  if ( size(compound) /= size(x) ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the array "compound" (names of considered species) and "x" (their mole fractions) have different sizes!'
     write(6,*) '   Size of "compound" ---: ', size(compound)
     write(6,*) '   Size of "x" ----------: ', size(x)
     write(6,*) ''
     stop
  end if

  if ( size(ln_phi_k) /= size(x) ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the array "ln_phi_k" (names of considered species) and "x" (their mole fractions) have different sizes!'
     write(6,*) '   Size of "ln_phi_k" ------: ', size(ln_phi_k)
     write(6,*) '   Size of "x" -------------: ', size(x)
     write(6,*) ''
     stop
  end if

  if ((P < 0._dp) .OR. (rho < 0._dp)) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the pressure and/or the density is/are negative!'
     write(6,*) '   P   = ', P
     write(6,*) '   rho = ', rho
     write(6,*) ''
     stop
  end if

  sum= 0._dp
  do i=1, size(x)
     sum= sum + x(i)
  end do
  if ( (sum > 1._dp + epsilon(sum)) .AND. normized_xi_sum ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the sum of mole fractions is larger than 1!'
     write(6,*) '   Sum x(i) = ', sum
     write(6,*) ''
     stop
  end if  

  ! -------------------------------------------------------------------------------------------------------------------------------- 
  !
  allocate(species_prop(1:size(compound)),stat=ok);            call alloc_error(subname,'species_prop',ok)
  allocate(k_ij(1:size(compound),1:size(compound)),stat=ok);   call alloc_error(subname,'k_ij',ok)
  
  ! -------------------------------------------------------------------------------------------------------------------------------- 
  ! We read the database of PC-SAFT parameters:
  if (first) then
     call read_database_pcsaft(subname,1,compdatabase,compdatabase_interac)
     first= .false.
     write(6,*) ''
     write(temp_string,'(I4)') size(compdatabase)
     call wf(' > Size of the database: ', 'bright red', &
                        temp_string(1:len_trim(temp_string)), 'bright black', &
                        ' species:', 'bright red')
     do i= 1, size(compdatabase)
        temp_string2= compdatabase(i)%name
        call wf('   - ', 'bright red', &
                        temp_string2(1:len_trim(temp_string2)), 'bright black' )
     end do
     write(6,*) ''
     call wf(' > Species taken into account in this computation: ', 'bright red')
     do i= 1, size(compound)
        temp_string2= compound(i)
        call wf('   - ', 'bright red', &
                        temp_string2(1:len_trim(temp_string2)), 'bright blue' )
     end do
     write(6,*) ''
  end if
  
  !do j= 1, size(compdatabase)
  !   write(6,'(A13,ES17.5,F12.4,F14.4,F15.4,A100)') compdatabase(j)%name,     compdatabase(j)%molmass, &
  !                                                  compdatabase(j)%numseg,   compdatabase(j)%segdiam, &
  !                                                  compdatabase(j)%segenerg, compdatabase(j)%ref
  !end do
  
  !do i= 1, size(x)
  !   write(6,*) '|',trim(ADJUSTL(compound(i))),'|'
  !end do
  
  ! Construction du tableau des paramtres d'interaction : -----------------------------
  do i= 1, size(compound)
     k_ij(i,i)= 0._dp ! On met des zros sur la diagonale
  end do
  
  do i= 1, size(compound)
     name1  = trim(ADJUSTL(compound(i)))
     do j= i+1, size(compound)
        name2  = trim(ADJUSTL(compound(j)))
        trouve = .false.
        do k= 1, size(compdatabase_interac)
           name3  = trim(ADJUSTL(compdatabase_interac(k)%name1))
           name4  = trim(ADJUSTL(compdatabase_interac(k)%name2))
           if (  ( (name1(1:len_trim(name1)) == name3(1:len_trim(name3))) .AND. &
                   (name2(1:len_trim(name2)) == name4(1:len_trim(name4)))   ) .OR. &
                 ( (name2(1:len_trim(name2)) == name3(1:len_trim(name3))) .AND. &
                   (name1(1:len_trim(name1)) == name4(1:len_trim(name4)))   )           ) then
              trouve = .true.
              !print*, compdatabase_interac(k)%a
              !print*, compdatabase_interac(k)%b
              !stop
              k_ij(i,j)= compdatabase_interac(k)%a + compdatabase_interac(k)%b * 1.E-4_dp * T 
              k_ij(j,i)= k_ij(i,j)
           end if
        end do
        if ( .NOT. trouve) then
           write(6,*) ''
           write(6,*) ' -----------------------------------------------------------------------------------------------------------'
           call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
           write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
           write(6,*) '   we did not find interaction parameter "k_ij" for the couple:'
           write(6,*) '   - ', name1
           write(6,*) '   - ', name2
           write(6,*) '   Please check the database file called: "', datafilename_inter(1:len_trim(datafilename_inter)), '" '
           write(6,*) ' -----------------------------------------------------------------------------------------------------------'
           write(6,*) ''
           stop
        end if
     end do
  end do
  ! Fin de la construction du tableau des paramtres d'interaction : -------------------
  
  do i= 1, size(x)
     name1  = trim(ADJUSTL(compound(i))) ! Enlve les blancs en dbut et fin de chane.
     here= .false.
     do j= 1, size(compdatabase)
        name2= trim(ADJUSTL(compdatabase(j)%name))
        !write(6,*) name2
        if ( name1(1:len_trim(name1)) == name2(1:len_trim(name2)) ) then
           here= .true.
           species_prop(i)%name     =  name1 ! Name
           species_prop(i)%x        =  x(i)  ! Mole fraction
           species_prop(i)%molmass  =  compdatabase(j)%molmass  ! The molar mass (in kg.mol^-1) for the considered compound.
           species_prop(i)%numseg   =  compdatabase(j)%numseg   ! The number of segment for the considered compound.
           species_prop(i)%segdiam  =  compdatabase(j)%segdiam  ! The segment diameter.
           species_prop(i)%segenerg =  compdatabase(j)%segenerg ! The segment energy  
           
           species_prop(i)%assoflag =  compdatabase(j)%assoflag ! Flag: 0=non associative species, 1= associative species
           species_prop(i)%assenerg =  compdatabase(j)%assenerg ! The association energy
           species_prop(i)%assvol   =  compdatabase(j)%assvol   ! The association volume
           exit
        end if
     end do
     if ( .NOT. here ) then
        write(6,*) ''
        call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
        write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
        write(6,*) '   the species "', name1(1:len_trim(name1)), '" is not present in the database!'
        write(6,*) ''
        stop
     end if
  end do
  
  !do i= 1, size(x)
  !   write(6,'(A13,ES15.5,ES17.5,F12.4,F14.4,F15.4)') species_prop(i)%name, species_prop(i)%x, species_prop(i)%molmass, &
  !                                                    species_prop(i)%numseg, species_prop(i)%segdiam,                  &
  !                                                    species_prop(i)%segenerg
  !end do

  ! -------------------------------------------------------------------------------------------------------------------------------- 
  ! We compute "Z", the "mu_k/kT"'s and the density "rho" in the frame of the PC-SAFT theory:

!!  call Zandmu_pcsaft(sub(1:len_trim(sub))//'/'//subname,1,state,species_prop,k_ij,P,T,Z,muskT_k,eta,rho)
  
  if (present(ares) ) then
     call Zandmu_pcsaft_Prho(subtemp,1,state,species_prop,k_ij, &
                             P,rho,Z,muskT_k,T,ares)
  else
     call Zandmu_pcsaft_Prho(subtemp,2,state,species_prop,k_ij,P,rho,Z,muskT_k,T)
  end if
  
  ! We then derive the fugacity coefficients "phi_k":
  !write(6,'(4(A,ES15.5))') 'P=', P, ' T= ', T, ' eta= ', eta, ' rho= ', rho
  do k= 1, size(species_prop)
     ln_phi_k(k)= muskT_k(k) - log(Z)
     !write(6,*) ' ln_phi_k(k)= ', ln_phi_k(k), muskT_k(k), log(Z)
  end do
  
  return
  
end subroutine pcsaft_Prho

!===================================================================================================================================
! D. Cordier, CNRS, France.
! http://orcid.org/0000-0003-4515-6271 
! 11 janvier 2021.
! Estimation of the "volumetric coefficient of thermal expansion" (cf. https://en.wikipedia.org/wiki/Thermal_expansion)
!
subroutine alphaV(sub, ncall, state, compound, x, P, T, alpha_V)
  use utils_DC
  
  implicit none
  character(len=*), parameter :: subname = 'alphaV'
  
  ! Inputs/outputs of the program:
  character(len=*), intent(in) :: sub
  integer, intent(in)          :: ncall
  character(len=len_trim(sub)+len_trim(subname)+1) :: subtemp
  
  character(len=1),                intent(in) :: state    ! the physical state of the phase:  - state='L' = liquid phase
                                                          !                                   - state='V' = vapor phase.
  character(len=comp_name_ncmax),  dimension(:), intent(in) :: compound ! List of the name of compound taken into account.
  real(dp),           dimension(:), intent(in) :: x        ! mole fractions of the various species taken into account.
  real(dp),                         intent(in) :: P        ! the pressure (in Pa).
  real(dp),                         intent(in) :: T        ! the temperature (in K).
  real(dp),                        intent(out) :: alpha_V      ! volumetric coefficient of thermal expansion (K^-1).

  real(dp), parameter :: pd = 1.E-4_dp ! For the numerical estimation of derivatives.
  real(dp) :: T1, T2, ares1, ares2, daressdT
  integer  :: i, ok
  real(dp) :: sum, rho, rho_1, rho_2, T_1, T_2
  real(dp), allocatable, dimension(:) :: muskT_k  ! The residual chemical potential divided by kT (no unit)
  real(dp), allocatable, dimension(:) :: ln_phi_k ! the fugacity coefficients in neperian log (no unit).
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Call path of this subroutine:
  subtemp= sub(1:len_trim(sub))//'/'//subname(1:len_trim(subname))
  
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! We check various obvious things:
  if ((state /= 'L') .AND. (state /= 'V')) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the considered physical state is neither Liquid nor Vapor!'
     write(6,*) '   We have: "state" = ', state
     write(6,*) ''
     stop
  end if
  
  if ( size(compound) /= size(x) ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the array "compound" (names of considered species) and "x" (their mole fractions) have different sizes!'
     write(6,*) '   Size of "compound" ---: ', size(compound)
     write(6,*) '   Size of "x" ----------: ', size(x)
     write(6,*) ''
     stop
  end if

  if ((P < 0._dp) .OR. (T < 0._dp)) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the pressure and/or the temperature is/are negative!'
     write(6,*) '   P = ', P
     write(6,*) '   T = ', T
     write(6,*) ''
     stop
  end if

  sum= 0._dp
  do i= 1, size(x)
     sum= sum + x(i)
  end do
  if ( (sum > 1._dp + epsilon(sum)) .AND. normized_xi_sum ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
     write(6,*) '   the sum of mole fractions is larger than 1!'
     write(6,*) '   Sum x(i) = ', sum
     write(6,*) '   ncall    = ', ncall
     write(6,*) ''
!     stop
  end if  
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! 
  allocate(muskT_k(1:size(compound)),stat=ok);   call alloc_error(subtemp,'muskT_k',ok)
  allocate(ln_phi_k(1:size(compound)),stat=ok);  call alloc_error(subtemp,'ln_phi_k',ok)
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! We compute the value of the density 'rho' (kg m^-3) for the given compisiton 'x' (mole fractions), pressure 'P' (Pa) and
  ! temperature 'T' (K):
  call pcsaft_PT(subtemp, 1, state, compound, x, P, T, rho, muskT_k, ln_phi_k)
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! We compute the finite difference require for the estmation of 'alpha_V':
  T_1 = T * (1._dp - pd)
  call pcsaft_PT(subtemp, 2, state, compound, x, P, T_1, rho_1, muskT_k, ln_phi_k)
  
  T_2 = T * (1._dp + pd)
  call pcsaft_PT(subtemp, 3, state, compound, x, P, T_2, rho_2, muskT_k, ln_phi_k)
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Finally, 'alpha_V' (K^-1): 
  alpha_V = -1._dp/rho * (rho_2 - rho_1) / (T_2 - T_1)
  
  deallocate(muskT_k,ln_phi_k)
  
  return

end subroutine alphaV

!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 
!                      -- Experimental isobaric heat capacity Cp for some organic solids --
!
!-----------------------------------------------------------------------------------------------------------------------------------
! D. Cordier, Institut UTINAM, Besançon, France.
! June 10th 2014 -- daniel.cordier@obs-besancon.fr
!
! - 15 avril 2020 : - passage des constantes numriques de l'criture '1.d0' vers l'criture '1._dp'
!                   - mise dans le module 'NEW_MOD_PCSAFT'.
!
subroutine cp_solid(sub,ncall,compound,T,Cpsolid)
  use utils_DC
  
  implicit none
  character(len=*), parameter :: subname = 'cp_solid'
  
  ! Inputs/outputs of the program:
  character(len=*), intent(in) :: sub
  integer, intent(in)          :: ncall
  character(len=len_trim(sub)+len_trim(subname)+1) :: subtemp
  
  character(len=comp_name_ncmax), intent(in) :: compound   ! List of the name of compound taken into account.
  real(dp),                        intent(in) :: T         ! the temperature (in K).
  
  real(dp),                        intent(out) :: Cpsolid ! heat capacity at constant pressure (J.K^-1.mol^-1).

  real(dp), parameter :: one_cal_in_joule = 4.184_dp
  
  integer :: ok, Npt, flag, i
  real(dp), dimension(:), allocatable :: Temp, Cp_exp
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Call path of this subroutine:
  subtemp= sub(1:len_trim(sub))//'/'//subname(1:len_trim(subname))

  !---------------------------------------------------------------------------------------------------------------------------------
  ! 
  flag= 0
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! 'C4H10' -- butane:
  if (compound(1:len_trim(compound)) .eq. 'C4H10') then
     flag = 1
     Npt= 49 ! Number of considered groups
     allocate(Temp(1:Npt),stat=ok);   call alloc_error(subtemp,'Temp',ok)
     allocate(Cp_exp(1:Npt),stat=ok); call alloc_error(subtemp,'Cp_exp',ok)
     
     Temp(1) =  12.73_dp; Cp_exp(1) =  0.484_dp;
     Temp(2) =  15.05_dp; Cp_exp(2) =  0.793_dp;
     Temp(3) =  18.44_dp; Cp_exp(3) =  1.510_dp;
     Temp(4) =  22.36_dp; Cp_exp(4) =  2.455_dp;     
     Temp(5) =  26.71_dp; Cp_exp(5) =  3.419_dp;
     Temp(6) =  31.53_dp; Cp_exp(6) =  4.724_dp;
     Temp(7) =  36.89_dp; Cp_exp(7) =  5.862_dp;
     Temp(8) =  42.03_dp; Cp_exp(8) =  6.980_dp;
     Temp(9) =  48.87_dp; Cp_exp(9) =  8.426_dp;
     Temp(10)=  54.59_dp; Cp_exp(10)=  9.412_dp;
     Temp(11)=  59.62_dp; Cp_exp(11)= 10.276_dp;
     Temp(12)=  60.72_dp; Cp_exp(12)= 10.412_dp;
     Temp(13)=  65.35_dp; Cp_exp(13)= 11.242_dp;
     Temp(14)=  66.27_dp; Cp_exp(14)= 11.313_dp;
     Temp(15)=  72.30_dp; Cp_exp(15)= 12.109_dp;
     Temp(16)=  74.07_dp; Cp_exp(16)= 12.412_dp;
     Temp(17)=  79.32_dp; Cp_exp(17)= 13.161_dp;
     Temp(18)=  85.76_dp; Cp_exp(18)= 14.052_dp;
     Temp(19)=  88.39_dp; Cp_exp(19)= 14.367_dp;
     Temp(20)=  92.28_dp; Cp_exp(20)= 14.927_dp;
     Temp(21)=  96.29_dp; Cp_exp(21)= 15.41_dp;
     Temp(22)=  98.17_dp; Cp_exp(22)= 15.72_dp;
     Temp(23)= 103.55_dp; Cp_exp(23)= 16.39_dp;
     Temp(24)= 113.77_dp; Cp_exp(24)= 19.92_dp;
     Temp(25)= 119.80_dp; Cp_exp(25)= 20.13_dp;
     Temp(26)= 126.78_dp; Cp_exp(26)= 20.47_dp;
     Temp(27)= 130.17_dp; Cp_exp(27)= 20.74_dp;
     Temp(28)= 133.39_dp; Cp_exp(28)= 21.19_dp;
     Temp(29)= 139.88_dp; Cp_exp(29)= 27.07_dp;     
     Temp(30)= 142.22_dp; Cp_exp(30)= 27.12_dp;
     Temp(31)= 149.50_dp; Cp_exp(31)= 27.42_dp;
     Temp(32)= 156.03_dp; Cp_exp(32)= 27.73_dp;
     Temp(33)= 162.29_dp; Cp_exp(33)= 27.71_dp;
     Temp(34)= 168.49_dp; Cp_exp(34)= 27.92_dp;
     Temp(35)= 174.59_dp; Cp_exp(35)= 27.92_dp;
     Temp(36)= 180.83_dp; Cp_exp(36)= 27.99_dp;
     Temp(37)= 186.63_dp; Cp_exp(37)= 28.22_dp;
     Temp(38)= 191.76_dp; Cp_exp(38)= 28.32_dp;
     Temp(39)= 196.83_dp; Cp_exp(39)= 28.51_dp;
     Temp(40)= 203.24_dp; Cp_exp(40)= 28.71_dp;
     Temp(41)= 209.00_dp; Cp_exp(41)= 28.78_dp;
     Temp(42)= 215.94_dp; Cp_exp(42)= 29.11_dp;
     Temp(43)= 222.71_dp; Cp_exp(43)= 29.29_dp;
     Temp(44)= 230.81_dp; Cp_exp(44)= 29.62_dp;
     Temp(45)= 238.43_dp; Cp_exp(45)= 29.84_dp;
     Temp(46)= 247.01_dp; Cp_exp(46)= 30.36_dp;
     Temp(47)= 251.43_dp; Cp_exp(47)= 30.89_dp;
     Temp(48)= 262.04_dp; Cp_exp(48)= 31.30_dp;
     Temp(49)= 268.14_dp; Cp_exp(49)= 31.46_dp;
        
     if ( ( T < Temp(1)) .OR. (T > Temp(size(Temp))) ) then
        write(6,*) ''
        call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
        write(6,*) ' > in subroutine "', subname(1:len_trim(subname)), '" called by: "', subtemp(1:len_trim(subtemp)), '"'
        write(6,*) '   T < Tmin or T > Tmax'
        write(6,*) '   T = ', T
        write(6,*) ''
        stop
     end if
     
     ! Conversion to J.K^-1.mol^-1 :   
     do i= 1, Npt
        Cp_exp(i)= Cp_exp(i) * one_cal_in_joule
     end do
     
     call intext_linear(subtemp,1,Temp,Cp_exp,T,Cpsolid)

     deallocate(Temp,Cp_exp)
  end if 

  !---------------------------------------------------------------------------------------------------------------------------------
  ! 'C6H6' -- benzene (Oliver et al., 1948):
  ! Tm = 279.1 K
  if (compound(1:len_trim(compound)) .eq. 'C6H6') then
     flag = 1
     Npt= 39 ! Number of considered groups
     allocate(Temp(1:Npt),stat=ok);   call alloc_error(subtemp,'Temp',ok)
     allocate(Cp_exp(1:Npt),stat=ok); call alloc_error(subtemp,'Cp_exp',ok)
     
     Temp(1) = 13._dp; Cp_exp(1) = 0.685_dp;
     Temp(2) = 14._dp; Cp_exp(2) = 0.830_dp;
     Temp(3) = 15._dp; Cp_exp(3) = 0.995_dp;
     Temp(4) = 20._dp; Cp_exp(4) = 2.000_dp;     
     Temp(5) = 25._dp; Cp_exp(5) = 3.145_dp;
     Temp(6) = 30._dp; Cp_exp(6) = 4.300_dp;
     Temp(7) = 35._dp; Cp_exp(7) = 5.385_dp;
     Temp(8) = 40._dp; Cp_exp(8) = 6.340_dp;
     Temp(9) = 45._dp; Cp_exp(9) = 7.165_dp;
     Temp(10)= 50._dp; Cp_exp(10)= 7.885_dp;
     Temp(11)= 55._dp; Cp_exp(11)= 8.505_dp;
     Temp(12)= 60._dp; Cp_exp(12)= 9.065_dp;
     Temp(13)= 65._dp; Cp_exp(13)= 9.540_dp;
     Temp(14)= 70._dp; Cp_exp(14)= 9.975_dp;
     Temp(15)= 75._dp; Cp_exp(15)=10.375_dp;
     Temp(16)= 80._dp; Cp_exp(16)=10.750_dp;
     Temp(17)= 85._dp; Cp_exp(17)=11.105_dp;
     Temp(18)= 90._dp; Cp_exp(18)=11.430_dp;
     Temp(19)= 95._dp; Cp_exp(19)=11.745_dp;
     Temp(20)=100._dp; Cp_exp(20)=12.050_dp;
     Temp(21)=110._dp; Cp_exp(21)=12.670_dp;
     Temp(22)=120._dp; Cp_exp(22)=13.310_dp;
     Temp(23)=130._dp; Cp_exp(23)=14.000_dp;
     Temp(24)=140._dp; Cp_exp(24)=14.700_dp;
     Temp(25)=150._dp; Cp_exp(25)=15.450_dp;
     Temp(26)=160._dp; Cp_exp(26)=16.230_dp;
     Temp(27)=170._dp; Cp_exp(27)=17.090_dp;
     Temp(28)=180._dp; Cp_exp(28)=18.020_dp;
     Temp(29)=190._dp; Cp_exp(29)=18.980_dp;     
     Temp(30)=200._dp; Cp_exp(30)=20.010_dp;
     Temp(31)=210._dp; Cp_exp(31)=21.140_dp;
     Temp(32)=220._dp; Cp_exp(32)=22.320_dp;
     Temp(33)=230._dp; Cp_exp(33)=23.550_dp;
     Temp(34)=240._dp; Cp_exp(34)=24.880_dp;
     Temp(35)=250._dp; Cp_exp(35)=26.300_dp;
     Temp(36)=260._dp; Cp_exp(36)=27.760_dp;
     Temp(37)=270._dp; Cp_exp(37)=29.310_dp;
     Temp(38)=278.69_dp; Cp_exp(38)=30.760_dp;
     Temp(39)=280._dp; Cp_exp(39)=31.59_dp;
     
     if ( ( T < Temp(1)) .OR. (T > Temp(size(Temp))) ) then
        write(6,*) ''
        call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
        write(6,*) ' > in subroutine "', subname(1:len_trim(subname)), '" called by: "', subtemp(1:len_trim(subtemp)), '"'
        write(6,*) '   T < Tmin or T > Tmax'
        write(6,*) '   T = ', T
        write(6,*) ''
        stop
     end if
     
     ! Conversion to J.K^-1.mol^-1 :   
     do i= 1, Npt
        Cp_exp(i)= Cp_exp(i) * one_cal_in_joule
     end do
     
     call intext_linear(subtemp,1,Temp,Cp_exp,T,Cpsolid)

     deallocate(Temp,Cp_exp)
  end if  
  

  !---------------------------------------------------------------------------------------------------------------------------------
  ! 'CH3CN' -- acetonitrile (Putnam et al. 1965):
  ! Tm = 229.3  K
  if (compound(1:len_trim(compound)) .eq. 'CH3CN') then
     flag = 1
     Npt= 44 ! Number of considered groups
     allocate(Temp(1:Npt),stat=ok);   call alloc_error(subtemp,'Temp',ok)
     allocate(Cp_exp(1:Npt),stat=ok); call alloc_error(subtemp,'Cp_exp',ok)
     
     Temp(1) =  20._dp; Cp_exp(1) =  0.882_dp;
     Temp(2) =  25._dp; Cp_exp(2) =  1.595_dp;
     Temp(3) =  30._dp; Cp_exp(3) =  2.407_dp;
     Temp(4) =  35._dp; Cp_exp(4) =  3.253_dp;     
     Temp(5) =  40._dp; Cp_exp(5) =  4.075_dp; 
     Temp(6) =  45._dp; Cp_exp(6) =  4.846_dp;
     Temp(7) =  50._dp; Cp_exp(7) =  5.570_dp;
     Temp(8) =  55._dp; Cp_exp(8) =  6.247_dp;
     Temp(9) =  60._dp; Cp_exp(9) =  6.868_dp;
     Temp(10)=  65._dp; Cp_exp(10)=  7.438_dp;
     Temp(11)=  70._dp; Cp_exp(11)=  7.959_dp;
     Temp(12)=  75._dp; Cp_exp(12)=  8.456_dp;
     Temp(13)=  80._dp; Cp_exp(13)=  8.933_dp;
     Temp(14)=  85._dp; Cp_exp(14)=  9.386_dp;
     Temp(15)=  90._dp; Cp_exp(15)=  9.817_dp;
     Temp(16)=  95._dp; Cp_exp(16)= 10.22_dp;
     Temp(17)= 100._dp; Cp_exp(17)= 10.62_dp;
     Temp(18)= 105._dp; Cp_exp(18)= 11.00_dp;
     Temp(19)= 110._dp; Cp_exp(19)= 11.38_dp;
     Temp(20)= 115._dp; Cp_exp(20)= 11.76_dp;
     Temp(21)= 120._dp; Cp_exp(21)= 12.14_dp;
     Temp(22)= 125._dp; Cp_exp(22)= 12.51_dp;
     Temp(23)= 130._dp; Cp_exp(23)= 12.88_dp;
     Temp(24)= 135._dp; Cp_exp(24)= 13.24_dp;
     Temp(25)= 140._dp; Cp_exp(25)= 13.60_dp;
     Temp(26)= 145._dp; Cp_exp(26)= 13.98_dp;
     Temp(27)= 150._dp; Cp_exp(27)= 14.34_dp;
     Temp(28)= 155._dp; Cp_exp(28)= 14.71_dp;
     Temp(29)= 160._dp; Cp_exp(29)= 15.08_dp;     
     Temp(30)= 165._dp; Cp_exp(30)= 15.45_dp;
     Temp(31)= 170._dp; Cp_exp(31)= 15.82_dp;
     Temp(32)= 175._dp; Cp_exp(32)= 16.19_dp;
     Temp(33)= 180._dp; Cp_exp(33)= 16.55_dp; 
     Temp(34)= 185._dp; Cp_exp(34)= 16.92_dp;
     Temp(35)= 190._dp; Cp_exp(35)= 17.30_dp;
     Temp(36)= 195._dp; Cp_exp(36)= 17.66_dp;
     Temp(37)= 200._dp; Cp_exp(37)= 18.04_dp;
     Temp(38)= 205._dp; Cp_exp(38)= 18.41_dp;
     Temp(39)= 210._dp; Cp_exp(39)= 18.78_dp;
     Temp(40)= 215._dp; Cp_exp(40)= 19.15_dp;
     Temp(41)= 220._dp; Cp_exp(41)= 17.73_dp;
     Temp(42)= 225._dp; Cp_exp(42)= 18.10_dp;
     Temp(43)= 230._dp; Cp_exp(43)= 20.91_dp;
     Temp(44)= 235._dp; Cp_exp(44)= 20.96_dp;

     
     if ( ( T < Temp(1)) .OR. (T > Temp(size(Temp))) ) then
        write(6,*) ''
        call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
        write(6,*) ' > in subroutine "', subname(1:len_trim(subname)), '" called by: "', subtemp(1:len_trim(subtemp)), '"'
        write(6,*) '   T < Tmin or T > Tmax'
        write(6,*) '   T = ', T
        write(6,*) ''
        stop
     end if
     
     ! Conversion to J.K^-1.mol^-1 :   
     do i= 1, Npt
        Cp_exp(i)= Cp_exp(i) * one_cal_in_joule
     end do
     
     call intext_linear(subtemp,1,Temp,Cp_exp,T,Cpsolid)

     deallocate(Temp,Cp_exp)
  end if  
  

  !---------------------------------------------------------------------------------------------------------------------------------
  ! 'HCN' -- hydrogen cyanide (Giauque & Ruehrwein):
  !  Tm = 260. K
  if (compound(1:len_trim(compound)) .eq. 'HCN') then
     flag = 1
     Npt= 49 ! Number of considered groups
     allocate(Temp(1:Npt),stat=ok);   call alloc_error(subtemp,'Temp',ok)
     allocate(Cp_exp(1:Npt),stat=ok); call alloc_error(subtemp,'Cp_exp',ok)
     
     Temp(1) =  79.79_dp; Cp_exp(1) = 6.795_dp;
     Temp(2) =  84.93_dp; Cp_exp(2) = 7.158_dp;
     Temp(3) =  89.95_dp; Cp_exp(3) = 7.456_dp;
     Temp(4) =  94.97_dp; Cp_exp(4) = 7.718_dp;     
     Temp(5) = 102.06_dp; Cp_exp(5) = 8.126_dp;
     Temp(6) = 106.94_dp; Cp_exp(6) = 8.349_dp;
     Temp(7) = 111.60_dp; Cp_exp(7) = 8.599_dp;
     Temp(8) = 116.85_dp; Cp_exp(8) = 8.833_dp;
     Temp(9) = 122.48_dp; Cp_exp(9) = 9.034_dp;
     Temp(10)= 124.67_dp; Cp_exp(10)= 9.121_dp;
     Temp(11)= 127.73_dp; Cp_exp(11)= 9.244_dp;
     Temp(12)= 133.19_dp; Cp_exp(12)= 9.453_dp;
     Temp(13)= 133.34_dp; Cp_exp(13)= 9.468_dp;
     Temp(14)= 138.42_dp; Cp_exp(14)= 9.631_dp;
     Temp(15)= 141.10_dp; Cp_exp(15)= 9.780_dp;
     Temp(16)= 143.62_dp; Cp_exp(16)= 9.882_dp;
     Temp(17)= 148.75_dp; Cp_exp(17)=10.040_dp;
     Temp(18)= 149.34_dp; Cp_exp(18)=10.120_dp;
     Temp(19)= 154.04_dp; Cp_exp(19)=10.270_dp;
     Temp(20)= 157.89_dp; Cp_exp(20)=10.420_dp;
     Temp(21)= 159.29_dp; Cp_exp(21)=10.440_dp;
     Temp(22)= 163.70_dp; Cp_exp(22)=10.630_dp;
     Temp(23)= 164.43_dp; Cp_exp(23)=10.620_dp;
     Temp(24)= 166.30_dp; Cp_exp(24)=10.700_dp;
     Temp(25)= 166.40_dp; Cp_exp(25)=10.720_dp;
     Temp(26)= 168.23_dp; Cp_exp(26)=10.810_dp;
     Temp(27)= 169.55_dp; Cp_exp(27)=11.930_dp;
     Temp(28)= 170.20_dp; Cp_exp(28)=12.190_dp;
     Temp(29)= 171.27_dp; Cp_exp(29)=11.680_dp;     
     Temp(30)= 172.31_dp; Cp_exp(30)=11.230_dp;
     Temp(31)= 172.73_dp; Cp_exp(31)=10.930_dp;
     Temp(32)= 174.42_dp; Cp_exp(32)=10.950_dp;
     Temp(33)= 175.81_dp; Cp_exp(33)=10.960_dp;
     Temp(34)= 177.04_dp; Cp_exp(34)=11.010_dp;
     Temp(35)= 177.58_dp; Cp_exp(35)=11.040_dp;
     Temp(36)= 181.51_dp; Cp_exp(36)=11.190_dp;
     Temp(37)= 183.61_dp; Cp_exp(37)=11.250_dp;
     Temp(38)= 190.30_dp; Cp_exp(38)=11.520_dp;
     Temp(39)= 196.32_dp; Cp_exp(39)=11.760_dp;
     
     Temp(40)= 202.10_dp; Cp_exp(40)=12.020_dp;
     Temp(41)= 207.76_dp; Cp_exp(41)=12.220_dp;
     Temp(42)= 213.43_dp; Cp_exp(42)=12.490_dp;
     Temp(43)= 224.74_dp; Cp_exp(43)=13.080_dp;
     Temp(44)= 230.94_dp; Cp_exp(44)=13.420_dp;
     Temp(45)= 236.66_dp; Cp_exp(45)=13.770_dp;
     Temp(46)= 242.51_dp; Cp_exp(46)=14.090_dp;
     Temp(47)= 248.43_dp; Cp_exp(47)=14.470_dp;
     Temp(48)= 254.14_dp; Cp_exp(48)=14.790_dp;
     Temp(49)= 266.57_dp; Cp_exp(49)=16.920_dp;
     
        
     if ( ( T < Temp(1)) .OR. (T > Temp(size(Temp))) ) then
        write(6,*) ''
        call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
        write(6,*) ' > in subroutine "', subname(1:len_trim(subname)), '" called by: "', subtemp(1:len_trim(subtemp)), '"'
        write(6,*) '   T < Tmin or T > Tmax'
        write(6,*) '   T = ', T
        write(6,*) ''
        stop
     end if
     
     ! Conversion to J.K^-1.mol^-1 :   
     do i= 1, Npt
        Cp_exp(i)= Cp_exp(i) * one_cal_in_joule
     end do
     
     call intext_linear(subtemp,1,Temp,Cp_exp,T,Cpsolid)

     deallocate(Temp,Cp_exp)
  end if  

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Security:
  if ( flag == 0 ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine "', subname(1:len_trim(subname)), '" called by: "', subtemp(1:len_trim(subtemp)), '"'
     write(6,*) '  the compound: ', compound(1:len_trim(compound)), ' is not present!'
     write(6,*) ''
     stop
  end if
  return
  
end subroutine cp_solid

!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 
!               -- Computation of heat capacity Cp, Cv and the speed of sound in the frame of the PC-SAFT theory  and
!                                                Joback's group-contribution method --
!
!-----------------------------------------------------------------------------------------------------------------------------------
! D. Cordier, Institut UTINAM, Besançon, France.
! June 4th 2014 -- daniel.cordier@obs-besancon.fr
!
! - 15 avril 2020 : - passage des constantes numriques de l'criture '1.d0' vers l'criture '1._dp'
!                   - mise dans le module 'NEW_MOD_PCSAFT'.
!
! Background: - for the equations giving specific heat and the speed of sound: see Diamantonis & Economou (2011)
!             - for some species, the ideal gas contribution to Cp (isochoric heat capacity) is computed thanks to the Joback's 
!               group-contribution method (see routine 'cp0_idealgas_joback' for details).
!
subroutine cpcv_ssound(sub,ncall,state,compound,x,P,T,rho,Cp,Cv,ssound)
  use utils_DC
  
  implicit none
  character(len=*), parameter :: subname = 'cpcv_ssound'
  
  ! Inputs/outputs of the program:
  character(len=*), intent(in) :: sub
  integer, intent(in)          :: ncall
  character(len=len_trim(sub)+len_trim(subname)+1) :: subtemp
  
  character(len=1),                intent(in) :: state    ! the physical state of the phase:  - state='L' = liquid phase
                                                          !                                   - state='V' = vapor phase.
  character(len=comp_name_ncmax), dimension(:), intent(in) :: compound ! List of the name of compound taken into account.
  real(dp),          dimension(:), intent(in) :: x        ! mole fractions of the various species taken into account.
  real(dp),                        intent(in) :: P        ! the pressure (in Pa).
  real(dp),                        intent(in) :: T        ! the temperature (in K).
  
  real(dp),                        intent(out) :: rho      ! the density of the mixture (kg.m^-3).
  real(dp),                        intent(out) :: Cp       ! heat capacity at constant pressure (J.K^-1.mol^-1).
  real(dp),                        intent(out) :: Cv       ! heat capacity at constant volume (J.K^-1.mol^-1).
  real(dp),                        intent(out) :: ssound   ! speed of sound (m.s^-1).
  
  logical :: first = .true.
  logical :: here
  integer :: i, j, ok
  real(dp), parameter :: pd= 1.E-3_dp
  real(dp) :: ares, ares1, ares2, aP_TT, aP_rho, aP_rhorho, aP_rhoT1, aP_rhoT2, aP_rhoT, Cp0_idealgas, Cv0_idealgas
  real(dp) :: Ptemp, Ttemp, P1, P2, k_T, alpha, Cpmoy
  real(dp), dimension(1:size(compound)) :: muskT_k, ln_phi_k

  type(compdata),      dimension(:), allocatable, save :: compdatabase
  type(compdataINTER), dimension(:), allocatable, save :: compdatabase_interac
  character(len=comp_name_ncmax) :: name1, name2, name3, name4
  real(dp), dimension(:), allocatable, save :: molmass
  real(dp) :: molmass_average
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Call path of this subroutine:
  subtemp= sub(1:len_trim(sub))//'/'//subname(1:len_trim(subname))

  ! -------------------------------------------------------------------------------------------------------------------------------- 
  ! We read the database of PC-SAFT parameters in order to get the molecular weights:
  if (first) then
     call read_database_pcsaft(subname,1,compdatabase,compdatabase_interac)
     allocate(molmass(1:size(x)),stat=ok); call alloc_error(subtemp,'molmass',ok)
     first= .false.
  end if

  do i= 1, size(x)
     name1  = trim(ADJUSTL(compound(i))) ! Enlve les blancs en dbut et fin de chane.
     here= .false.
     do j= 1, size(compdatabase)
        name2= trim(ADJUSTL(compdatabase(j)%name))
        !write(6,*) name2
        if ( name1(1:len_trim(name1)) == name2(1:len_trim(name2)) ) then
           here= .true.
           molmass(i)= compdatabase(j)%molmass  ! The molar mass (in kg.mol^-1) for the considered compound.
           exit
        end if
     end do
     if ( .NOT. here ) then
        write(6,*) ''
        call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
        write(6,*) ' > In subroutine ''', subname(1:len_trim(subname)), ''' called by ''', subtemp(1:len_trim(subtemp)), ''''
        write(6,*) '   the species "', name1(1:len_trim(name1)), '" is not present in the database!'
        write(6,*) ''
        stop
     end if
  end do

  ! -------------------------------------------------------------------------------------------------------------------------------- 
  ! Average molar mass:    
  molmass_average= 0._dp
  do i= 1, size(x)
     molmass_average= molmass_average + molmass(i)
  end do
  molmass_average= molmass_average / size(x)
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Computation of 'rho' at the given pressure 'P' and temperature 'T' :
  call pcsaft_PT(subtemp,1,state,compound,x,P,T,rho,muskT_k,ln_phi_k)
  
  !print*, ' Check 0: T= ', T, ' -- P= ', P, ' rho= ', rho

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Computation of the SECOND derivative of 'ares' T :  
  
  ! Premire drive par rapport  rho :
  call pcsaft_rhoT(subtemp,1,state,compound,x,rho,T,Ptemp,muskT_k,ln_phi_k,ares)
  
  call pcsaft_rhoT(subtemp,1,state,compound,x,rho,T*(1._dp-pd),Ptemp,muskT_k,ln_phi_k,ares1)
  call pcsaft_rhoT(subtemp,2,state,compound,x,rho,T*(1._dp+pd),Ptemp,muskT_k,ln_phi_k,ares2)
  
  aP_TT= (ares2-2._dp*ares+ares1)/T**2/pd**2 
  
  !print*, ' Check 1: P= ', P, ' -- Ptemp= ', Ptemp, ' rho= ', rho
  !print*, ' >>> aP_TT= ', aP_TT, ' -T * d2ares/dT2= ', -T*aP_TT

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Computation of the derivative '(dP/drho)_T' and 'k_T ':  
  
  call pcsaft_rhoT(subtemp,3,state,compound,x,rho*(1._dp-pd),T,P1,muskT_k,ln_phi_k)
  call pcsaft_rhoT(subtemp,4,state,compound,x,rho*(1._dp+pd),T,P2,muskT_k,ln_phi_k)
  
  k_T= (rho * (P2-P1)/2._dp/rho/pd )**(-1._dp)
  
  !print*, ' >>> k_T= ', k_T
  
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Computation of the derivative '(dP/dT)_v' and 'alpha':
  
  call pcsaft_rhoT(subtemp,5,state,compound,x,rho,T*(1._dp-pd),P1,muskT_k,ln_phi_k)
  call pcsaft_rhoT(subtemp,6,state,compound,x,rho,T*(1._dp+pd),P2,muskT_k,ln_phi_k)
  
  alpha= k_T * (P2-P1)/2._dp/T/pd
  
  !print*, ' >>> alpha= ', alpha, T*alpha**2/k_T/rho * molmass(1)
  
  Cpmoy= 0._dp
  do i= 1, size(x)
     call cp0_idealgas_joback(subtemp,1,compound(i),T,Cp0_idealgas)
     Cpmoy= Cpmoy + x(i) * Cp0_idealgas
  end do
  Cp0_idealgas= Cpmoy
  
  !print*, ' >>> Cp0_idealgas= ', Cp0_idealgas
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Mayer's relation:
  Cv0_idealgas= Cp0_idealgas - R_gas

  !print*, ' Cv0_idealgas= ', Cv0_idealgas
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! PC-SAFT isochoric heat capacity:  
  Cv= Cv0_idealgas - T*aP_TT ! cf. Eq. (19) Diamantonis & Economou (2011).

  !print*, 'c1 = Cv-Cp0_idealgas= ', Cv-Cp0_idealgas
  !stop
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! PC-SAFT isobaric heat capacity:    
  Cp = Cv + T*alpha**2/k_T/rho * molmass_average

  !print*, 'c2 = Cp-Cp0_idealgas= ', Cp-Cp0_idealgas
  !stop
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! PC-SAFT speed of sound:
  
  ssound= sqrt(Cp/Cv /k_T /rho) ! m.s^-1
  
  !print*, ' T= ', T
  !print*, ' P= ', P
  !print*, ' k_T= ', k_T 
  !print*, ' rho= ', rho
  !stop
  
  !print*, ''
  !print*, ' >>>>> Cp ( T= ', T, ' )= ', Cp0_idealgas -T*aP_TT + T*alpha**2/k_T/rho * molmass_average
  !print*, ''
  
  !print*, ''
  !print*, ' >>>>> ssound ( T= ', T, ' )= ', ssound, ' m.s^-1'
  !print*, ''
 
  return
  
end subroutine cpcv_ssound

!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 
!                          -- Computation of ideal gas heat capacity Cp with the Joback method --
!
!-----------------------------------------------------------------------------------------------------------------------------------
! D. Cordier, Institut UTINAM, Besançon, France.
! June 4th 2014 -- daniel.cordier@obs-besancon.fr
!
! - 15 avril 2020 : - passage des constantes numriques de l'criture '1.d0' vers l'criture '1._dp'
!                   - mise dans le module 'NEW_MOD_PCSAFT'.
!
subroutine cp0_idealgas_joback(sub,ncall,compound,T,Cp0_idealgas)
  use utils_DC
  
  implicit none
  character(len=*), parameter :: subname = 'cp0_idealgas_joback'
  
  ! Inputs/outputs of the program:
  character(len=*), intent(in) :: sub
  integer, intent(in)          :: ncall
  character(len=len_trim(sub)+len_trim(subname)+1) :: subtemp
  
  character(len=comp_name_ncmax),  intent(in)  :: compound     ! the name of the compound 
  real(dp),                        intent(in)  :: T            ! the temperature (in K).
  real(dp),                        intent(out) :: Cp0_idealgas ! heat capacity at constant pressure (J.K^-1.mol^-1).

  integer :: i, flag
  integer :: Ng, ok
  real(dp), dimension(:),  allocatable :: a, b, c, d
  integer,  dimension(:),  allocatable :: Nk
  real(dp),  dimension(:), allocatable :: n, Theta
  real(dp) :: sum_A, sum_B, sum_C, sum_D
  real(dp) :: sum
  real(dp) :: Cp1, Cp2, T1, T2, aa, bb
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Call path of this subroutine:
  subtemp= sub(1:len_trim(sub))//'/'//subname(1:len_trim(subname))
  
  !---------------------------------------------------------------------------------------------------------------------------------
  flag= 0

  !---------------------------------------------------------------------------------------------------------------------------------
  ! 'N2' -- nitrogen:
  if (compound(1:len_trim(compound)) .eq. 'N2') then
     flag = 1     
     Cp0_idealgas=  7._dp/2._dp * R_gas ! Cela fait environ 29 J.K-1.mol-1
     ! La valeur "gaz parfait" ci-dessus ne conduit pas  un rsultat satisfaisant  77K pour la vitesse du son publie dans
     ! le rapport NASA de Zuckerwar & Mazel (1985), en prenant 859 m.s^-1  77 K on trouve :
     Cp0_idealgas=  16.8_dp ! J.K^-1.mol^-1
  end if 

  !---------------------------------------------------------------------------------------------------------------------------------
  ! 'CH4' -- methane (cf. Setzmann & Wagner, 1991, p. 1086, Eq. 3.8):
!  if (compound(1:len_trim(compound)) .eq. 'CH4') then
!     flag = 1    
!     allocate(n(0:5),stat=ok);     call alloc_error(subtemp,'n',ok)
!     allocate(Theta(1:5),stat=ok); call alloc_error(subtemp,'Theta',ok)
!     
!     n(0)= 4.0016d0
!     n(1)= 0.008449d0
!     n(2)= 4.6942d0
!     n(3)= 3.4865d0
!     n(4)= 1.6572d0
!     n(5)= 1.4115d0
!     
!     Theta(1)=   648.d0
!     Theta(2)=  1957.d0
!     Theta(3)=  3895.d0
!     Theta(4)=  5705.d0
!     Theta(5)= 15080.d0
!     
!     sum= n(0)
!     do i= 1, 5
!        sum= sum + n(i) * Theta(i)**2 /T**2 * exp(Theta(i)/T) /(exp(Theta(i)/T) - 1.d0)**2
!     end do
!     
!     Cp0_idealgas=  R_gas * sum
!     
!     !print*, ' >>> Cp0_idealgas(CH4)= ', Cp0_idealgas
!     deallocate(n,Theta)
!  end if 

  ! Dtermiantion des capacites  90.729 K et 100.0 K du gaz parfait correspondant afin de retrouver les vitesses du son
  ! publies dans Table 40 p. 1126 de Setzmann & Wagner (1991) 
  if (compound(1:len_trim(compound)) .eq. 'CH4') then
     flag = 1  
     
     T1= 90.729_dp ! K
     T2= 100.0_dp  ! K
     
     Cp1= 43.5_dp ! J.mol^-1.K^-1   
     
     
     Cp2= 44.0_dp ! J.mol^-1.K^-1
          
     aa= (Cp2-Cp1)/(T2-T1)
     bb= Cp2 - aa*T2
     
     Cp0_idealgas= aa * T + bb
     
     !print*, ' >>> Cp0_idealgas(CH4)= ', Cp0_idealgas

  end if 

  !---------------------------------------------------------------------------------------------------------------------------------
  ! 'CH4' -- methane ->> deuxime version !!
!!  if (compound(1:len_trim(compound)) .eq. 'CH4') then
!!     flag = 1
!!     Ng= 1 ! Number of considered groups
!!     allocate(Nk(1:Ng),stat=ok); call alloc_error(subtemp,'d',ok)
!!     allocate(a(1:Ng),stat=ok);  call alloc_error(subtemp,'a',ok)
!!     allocate(b(1:Ng),stat=ok);  call alloc_error(subtemp,'b',ok)
!!     allocate(c(1:Ng),stat=ok);  call alloc_error(subtemp,'c',ok)
!!     allocate(d(1:Ng),stat=ok);  call alloc_error(subtemp,'d',ok)
!!     
!!     Nk(1)= 1; a(1)= 19.500d0; b(1)= -0.00808d0; c(1)= 1.53d-4;   d(1)= -9.70d-8; ! -CH3
!!          
!!     sum_A=0.d0; sum_B=0.d0; sum_C=0.d0; sum_D=0.d0;
!!     do i= 1, Ng
!!        sum_A= sum_A + Nk(i) * a(i)
!!        sum_B= sum_B + Nk(i) * b(i)
!!        sum_C= sum_C + Nk(i) * c(i)
!!        sum_D= sum_D + Nk(i) * d(i)
!!     end do
!!     
!!     Cp0_idealgas=  sum_A - 37.93d0         + &
!!                   (sum_B + 0.210d0) * T    + &
!!                   (sum_C - 3.91d-4) * T**2.d0 + &
!!                   (sum_D + 2.06d-7) * T**3.d0     ! Cp0 in J.mol^-1.K^-1
!!     deallocate(Nk,a,b,c,d)
!!     
!!     !*, ' > Cp0_idealgas=  ', Cp0_idealgas !, sum_A, sum_B, sum_C, sum_D, T
!!     !print*, sum_A - 37.93d0, (sum_B + 0.210d0) * T , (sum_B - 3.91d-4) , (sum_C + 2.06d-7) * T**3.d0
!!  end if 
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! 'C2H6' -- ethane:
!  if (compound(1:len_trim(compound)) .eq. 'C2H6') then
!     flag = 1
!     Ng= 1 ! Number of considered groups
!     allocate(Nk(1:Ng),stat=ok); call alloc_error(subtemp,'d',ok)
!     allocate(a(1:Ng),stat=ok);  call alloc_error(subtemp,'a',ok)
!     allocate(b(1:Ng),stat=ok);  call alloc_error(subtemp,'b',ok)
!     allocate(c(1:Ng),stat=ok);  call alloc_error(subtemp,'c',ok)
!     allocate(d(1:Ng),stat=ok);  call alloc_error(subtemp,'d',ok)
!     
!     Nk(1)= 2; a(1)= 19.500d0; b(1)= -0.00808d0; c(1)= 1.53d-4;   d(1)= -9.70d-8; ! -CH3
!          
!     sum_A=0.d0; sum_B=0.d0; sum_C=0.d0; sum_D=0.d0;
!     do i= 1, Ng
!        sum_A= sum_A + Nk(i) * a(i)
!        sum_B= sum_B + Nk(i) * b(i)
!        sum_C= sum_C + Nk(i) * c(i)
!        sum_D= sum_D + Nk(i) * d(i)
!     end do
!     
!     Cp0_idealgas=  sum_A - 37.93d0         + &
!                   (sum_B + 0.210d0) * T    + &
!                   (sum_C - 3.91d-4) * T**2.d0 + &
!                   (sum_D + 2.06d-7) * T**3.d0     ! Cp0 in J.mol^-1.K^-1
!     deallocate(Nk,a,b,c,d)
!     !*, ' > Cp0_idealgas=  ', Cp0_idealgas !, sum_A, sum_B, sum_C, sum_D, T
!     !print*, sum_A - 37.93d0, (sum_B + 0.210d0) * T , (sum_B - 3.91d-4) , (sum_C + 2.06d-7) * T**3.d0
!  end if 

  if (compound(1:len_trim(compound)) .eq. 'C2H6') then
     flag = 1
     ! 

     T1= 90.36_dp ! K
     T2= 100.0_dp ! K
     
     Cp1= 2.6_dp ! J.mol^-1.K^-1     
     Cp2= 4.60_dp ! J.mol^-1.K^-1
          
     aa= (Cp2-Cp1)/(T2-T1)
     bb= Cp2 - aa*T2
     
     Cp0_idealgas= aa * T + bb
     
     !print*, ' > Cp0_idealgas(C2H6)= ', Cp0_idealgas

  end if 
   
  !---------------------------------------------------------------------------------------------------------------------------------
  ! 'C3H8' -- propane:
  if (compound(1:len_trim(compound)) .eq. 'C3H8') then
     flag = 1
     ! Ajustement  partir des donnes prises sur la page 114 de Younglove & Ely (1987) :
     T1= 89.49_dp ! K
     T2= 200.0_dp ! K
     
     ! Les deux valeurs de Cp0_idealgas ajustes pour 85.48 K et 200 K (16 dc. 2015) :
     Cp1= -8.46700_dp    ! J.mol^-1.K^-1
     Cp2= 20.542187_dp ! J.mol^-1.K^-1
          
     aa= (Cp2-Cp1)/(T2-T1)
     bb= Cp2 - aa*T2
     
     Cp0_idealgas= aa * T + bb
     
     !print*, ' > Cp0_idealgas(C3H8)= ', Cp0_idealgas

  end if 
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Test computation for 'ethyphenol' (see example p. 3.6 of Poling et al.):
  if (compound(1:len_trim(compound)) .eq. 'ethylphenol') then
     flag = 1
     Ng= 5 ! Number of considered groups
     allocate(Nk(1:Ng),stat=ok); call alloc_error(subtemp,'d',ok)
     allocate(a(1:Ng),stat=ok);  call alloc_error(subtemp,'a',ok)
     allocate(b(1:Ng),stat=ok);  call alloc_error(subtemp,'b',ok)
     allocate(c(1:Ng),stat=ok);  call alloc_error(subtemp,'c',ok)
     allocate(d(1:Ng),stat=ok);  call alloc_error(subtemp,'d',ok)
     
     Nk(1)= 1; a(1)= 19.500_dp;  b(1)= -0.00808_dp; c(1)= 1.53E-4_dp;   d(1)= -9.70E-8_dp;  ! -CH3
     Nk(2)= 1; a(2)= -0.909_dp;  b(2)= 0.09500_dp;  c(2)= -0.54E-4_dp;  d(2)=  1.19E-8_dp;  ! -CH2-
     Nk(3)= 4; a(3)= -2.140_dp;  b(3)= 5.74E-2_dp;  c(3)= -1.64E-6_dp;  d(3)= -1.59E-8_dp ; ! =CH(ds)
     Nk(4)= 2; a(4)= -8.2500_dp; b(4)= 1.010E-1_dp; c(4)= -1.42E-4_dp;  d(4)= 6.78E-8_dp;   ! =C(ds)          
     Nk(5)= 1; a(5)= -2.810_dp;  b(5)= 0.11100_dp;  c(5)= -1.16E-4_dp;  d(5)=  4.94E-8_dp;  ! -ACOH
     
     sum_A=0._dp; sum_B=0._dp; sum_C=0._dp; sum_D=0._dp;
     do i= 1, Ng
        sum_A= sum_A + Nk(i) * a(i)
        sum_B= sum_B + Nk(i) * b(i)
        sum_C= sum_C + Nk(i) * c(i)
        sum_D= sum_D + Nk(i) * d(i)
     end do
     
     Cp0_idealgas=  sum_A - 37.93_dp         + &
                   (sum_B + 0.210_dp) * T    + &
                   (sum_C - 3.91E-4_dp) * T**2._dp + &
                   (sum_D + 2.06E-7_dp) * T**3._dp     ! Cp0 in J.mol^-1.K^-1
     deallocate(Nk,a,b,c,d)
     
     !print*, ' > Cp0_idealgas=  ', Cp0_idealgas !, sum_A, sum_B, sum_C, sum_D, T
     !print*, sum_A - 37.93d0, (sum_B + 0.210d0) * T , (sum_B - 3.91d-4) , (sum_C + 2.06d-7) * T**3.d0
  end if 
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! 'octane':
  if (compound(1:len_trim(compound)) .eq. 'octane') then
     flag = 1
     Ng= 2 ! Number of considered groups
     allocate(Nk(1:Ng),stat=ok); call alloc_error(subtemp,'d',ok)
     allocate(a(1:Ng),stat=ok);  call alloc_error(subtemp,'a',ok)
     allocate(b(1:Ng),stat=ok);  call alloc_error(subtemp,'b',ok)
     allocate(c(1:Ng),stat=ok);  call alloc_error(subtemp,'c',ok)
     allocate(d(1:Ng),stat=ok);  call alloc_error(subtemp,'d',ok)
     
     Nk(1)= 2; a(1)= 19.500_dp; b(1)= -0.00808_dp; c(1)= 1.53E-4_dp;   d(1)= -9.70E-8_dp; ! -CH3
     Nk(2)= 6; a(2)= -0.909_dp; b(2)= 0.09500_dp;  c(2)= -0.54E-4_dp;  d(2)=  1.19E-8_dp; ! -CH2-
          
     sum_A=0._dp; sum_B=0._dp; sum_C=0._dp; sum_D=0._dp;
     do i= 1, Ng
        sum_A= sum_A + Nk(i) * a(i)
        sum_B= sum_B + Nk(i) * b(i)
        sum_C= sum_C + Nk(i) * c(i)
        sum_D= sum_D + Nk(i) * d(i)
     end do
     
     Cp0_idealgas=  sum_A - 37.93_dp         + &
                   (sum_B + 0.210_dp) * T    + &
                   (sum_C - 3.91E-4_dp) * T**2._dp + &
                   (sum_D + 2.06E-7_dp) * T**3._dp     ! Cp0 in J.mol^-1.K^-1
     deallocate(Nk,a,b,c,d)
     
     !print*, ' > Cp0_idealgas=  ', Cp0_idealgas !, sum_A, sum_B, sum_C, sum_D, T
     !print*, sum_A - 37.93d0, (sum_B + 0.210d0) * T , (sum_B - 3.91d-4) , (sum_C + 2.06d-7) * T**3.d0
  end if 

  !---------------------------------------------------------------------------------------------------------------------------------
  ! 'C4H10' -- butane:
  if (compound(1:len_trim(compound)) .eq. 'C4H10') then
     flag = 1
     Ng= 2 ! Number of considered groups
     allocate(Nk(1:Ng),stat=ok); call alloc_error(subtemp,'d',ok)
     allocate(a(1:Ng),stat=ok);  call alloc_error(subtemp,'a',ok)
     allocate(b(1:Ng),stat=ok);  call alloc_error(subtemp,'b',ok)
     allocate(c(1:Ng),stat=ok);  call alloc_error(subtemp,'c',ok)
     allocate(d(1:Ng),stat=ok);  call alloc_error(subtemp,'d',ok)
     
     Nk(1)= 2; a(1)= 19.500_dp; b(1)= -0.00808_dp; c(1)= 1.53E-4_dp;   d(1)= -9.70E-8_dp; ! -CH3
     Nk(2)= 2; a(2)= -0.909_dp; b(2)= 0.09500_dp;  c(2)= -0.54E-4_dp;  d(2)=  1.19E-8_dp; ! -CH2-
          
     sum_A=0._dp; sum_B=0._dp; sum_C=0._dp; sum_D=0._dp;
     do i= 1, Ng
        sum_A= sum_A + Nk(i) * a(i)
        sum_B= sum_B + Nk(i) * b(i)
        sum_C= sum_C + Nk(i) * c(i)
        sum_D= sum_D + Nk(i) * d(i)
     end do
     
     Cp0_idealgas=  sum_A - 37.93_dp         + &
                   (sum_B + 0.210_dp) * T    + &
                   (sum_C - 3.91E-4_dp) * T**2._dp + &
                   (sum_D + 2.06E-7_dp) * T**3._dp     ! Cp0 in J.mol^-1.K^-1
     deallocate(Nk,a,b,c,d)
     
     !*, ' > Cp0_idealgas=  ', Cp0_idealgas !, sum_A, sum_B, sum_C, sum_D, T
     !print*, sum_A - 37.93d0, (sum_B + 0.210d0) * T , (sum_B - 3.91d-4) , (sum_C + 2.06d-7) * T**3.d0
  end if 

  !---------------------------------------------------------------------------------------------------------------------------------
  ! 'C6H6' -- benzene:
  if ( (compound(1:len_trim(compound)) .eq. 'benzene') .OR. (compound(1:len_trim(compound)) .eq. 'C6H6') ) then
     flag = 1
     Ng= 1 ! Number of considered groups
     allocate(Nk(1:Ng),stat=ok); call alloc_error(subtemp,'d',ok)
     allocate(a(1:Ng),stat=ok);  call alloc_error(subtemp,'a',ok)
     allocate(b(1:Ng),stat=ok);  call alloc_error(subtemp,'b',ok)
     allocate(c(1:Ng),stat=ok);  call alloc_error(subtemp,'c',ok)
     allocate(d(1:Ng),stat=ok);  call alloc_error(subtemp,'d',ok)
     
     ! Coefficients dans la Table C-1 de Poling et al. :
     ! 2.140 5.74E-02 1.64E-06 1.59E-08
     Nk(1)= 6; a(1)= 2.140_dp; b(1)= 5.74E-2_dp; c(1)= 1.64E-6_dp;   d(1)= 1.59E-8_dp; ! -CH(ds) (in aromatic ring)

     sum_A=0._dp; sum_B=0._dp; sum_C=0._dp; sum_D=0._dp;
     do i= 1, Ng
        sum_A= sum_A + Nk(i) * a(i)
        sum_B= sum_B + Nk(i) * b(i)
        sum_C= sum_C + Nk(i) * c(i)
        sum_D= sum_D + Nk(i) * d(i)
     end do
     
     Cp0_idealgas=  sum_A - 37.93_dp         + &
                   (sum_B + 0.210_dp) * T    + &
                   (sum_C - 3.91E-4_dp) * T**2._dp + &
                   (sum_D + 2.06E-7_dp) * T**3._dp     ! Cp0 in J.mol^-1.K^-1
     deallocate(Nk,a,b,c,d)
     
     !print*, ' > Cp0_idealgas=  ', Cp0_idealgas !, sum_A, sum_B, sum_C, sum_D, T
     !print*, sum_A - 37.93d0, (sum_B + 0.210d0) * T , (sum_B - 3.91d-4) , (sum_C + 2.06d-7) * T**3.d0
  end if 

  !---------------------------------------------------------------------------------------------------------------------------------
  ! 'CH3CN' -- acetonitrile:
  if (compound(1:len_trim(compound)) .eq. 'CH3CN') then
     flag = 1
     Ng= 2 ! Number of considered groups
     allocate(Nk(1:Ng),stat=ok); call alloc_error(subtemp,'d',ok)
     allocate(a(1:Ng),stat=ok);  call alloc_error(subtemp,'a',ok)
     allocate(b(1:Ng),stat=ok);  call alloc_error(subtemp,'b',ok)
     allocate(c(1:Ng),stat=ok);  call alloc_error(subtemp,'c',ok)
     allocate(d(1:Ng),stat=ok);  call alloc_error(subtemp,'d',ok)
     
     ! Coefficients dans la Table C-1 de Poling et al. :
     ! Pour -CN : 36.500 7.33E-02 1.84E-04 1.03E-07
     
     Nk(1)= 1; a(1)= 19.500_dp; b(1)= -0.00808_dp; c(1)= 1.53E-4_dp;   d(1)= -9.70E-8_dp; ! -CH3
     Nk(2)= 1; a(2)= 36.500_dp; b(2)= 7.33E-2_dp; c(2)=  1.84E-4_dp;   d(2)= 1.03E-7_dp; ! -CN

     sum_A=0.d0; sum_B=0.d0; sum_C=0.d0; sum_D=0.d0;
     do i= 1, Ng
        sum_A= sum_A + Nk(i) * a(i)
        sum_B= sum_B + Nk(i) * b(i)
        sum_C= sum_C + Nk(i) * c(i)
        sum_D= sum_D + Nk(i) * d(i)
     end do
     
     Cp0_idealgas=  sum_A - 37.93_dp         + &
                   (sum_B + 0.210_dp) * T    + &
                   (sum_C - 3.91E-4_dp) * T**2._dp + &
                   (sum_D + 2.06E-7_dp) * T**3._dp     ! Cp0 in J.mol^-1.K^-1
     deallocate(Nk,a,b,c,d)
     
     !print*, ' > Cp0_idealgas=  ', Cp0_idealgas !, sum_A, sum_B, sum_C, sum_D, T
     !print*, sum_A - 37.93d0, (sum_B + 0.210d0) * T , (sum_B - 3.91d-4) , (sum_C + 2.06d-7) * T**3.d0
  end if 

  !---------------------------------------------------------------------------------------------------------------------------------
  ! 'HCN' -- hydrogen cyanide:
  if (compound(1:len_trim(compound)) .eq. 'HCN') then
     flag = 1
     Ng= 1 ! Number of considered groups
     allocate(Nk(1:Ng),stat=ok); call alloc_error(subtemp,'d',ok)
     allocate(a(1:Ng),stat=ok);  call alloc_error(subtemp,'a',ok)
     allocate(b(1:Ng),stat=ok);  call alloc_error(subtemp,'b',ok)
     allocate(c(1:Ng),stat=ok);  call alloc_error(subtemp,'c',ok)
     allocate(d(1:Ng),stat=ok);  call alloc_error(subtemp,'d',ok)
     
     ! Coefficients dans la Table C-1 de Poling et al. :
     ! Pour -CN : 36.500 7.33E-02 1.84E-04 1.03E-07
     
     Nk(1)= 1; a(1)= 36.500_dp; b(1)= 7.33E-2_dp; c(1)=  1.84E-4_dp;   d(1)= 1.03E-7_dp; ! -CN

     sum_A=0._dp; sum_B=0._dp; sum_C=0._dp; sum_D=0._dp;
     do i= 1, Ng
        sum_A= sum_A + Nk(i) * a(i)
        sum_B= sum_B + Nk(i) * b(i)
        sum_C= sum_C + Nk(i) * c(i)
        sum_D= sum_D + Nk(i) * d(i)
     end do
     
     Cp0_idealgas=  sum_A - 37.93_dp         + &
                   (sum_B + 0.210_dp) * T    + &
                   (sum_C - 3.91E-4_dp) * T**2._dp + &
                   (sum_D + 2.06E-7_dp) * T**3._dp     ! Cp0 in J.mol^-1.K^-1
     deallocate(Nk,a,b,c,d)
     
     !print*, ' > Cp0_idealgas=  ', Cp0_idealgas !, sum_A, sum_B, sum_C, sum_D, T
     !print*, sum_A - 37.93d0, (sum_B + 0.210d0) * T , (sum_B - 3.91d-4) , (sum_C + 2.06d-7) * T**3.d0
  end if 

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Security:
  if ( flag == 0 ) then
     write(6,*) ''
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,*) ' > In subroutine "', subname(1:len_trim(subname)), '" called by: "', subtemp(1:len_trim(subtemp)), '"'
     write(6,*) '   the compound "', compound(1:len_trim(compound)), '" is not present!'
     write(6,*) ''
     stop
  end if
  
  return
  
end subroutine cp0_idealgas_joback

!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 
!                               -- PC-SAFT universal constant (Gross & Sadowski, 2001) --
!
!-----------------------------------------------------------------------------------------------------------------------------------
! D. Cordier, Institut UTINAM, Besançon, France.
! 1er novembre 2013 -- daniel.cordier@obs-besancon.fr
!
! - 15 avril 2020 : - passage des constantes numriques de l'criture '1.d0' vers l'criture '1._dp'
!                   - mise dans le module 'NEW_MOD_PCSAFT'.
!
subroutine univ_cst(sub,ncall,a,b)
  use utils_DC
  
  implicit none
  character(len=*), parameter :: subname = 'univ_cst'
  
  ! Inputs/outputs of the program:
  character(len=*), intent(in) :: sub
  integer, intent(in)          :: ncall

  !
  real(dp), dimension(0:2,0:6), intent(out) :: a, b
  
  ! PC-SAFT theory universal constants from Table 1 of Gross & Sadowski (2001) :
  
  !i a0i          a1i           a2i            b0i          b1i           b2i
  !0 0.9105631445 -0.3084016918 -0.0906148351 0.7240946941 -0.5755498075 0.0976883116
  !1 0.6361281449 0.1860531159 0.4527842806 2.2382791861 0.6995095521 -0.2557574982
  !2 2.6861347891 -2.5030047259 0.5962700728 -4.0025849485 3.8925673390 -9.1558561530
  !3 -26.547362491 21.419793629 -1.7241829131 -21.003576815 -17.215471648 20.642075974
  !4 97.759208784 -65.255885330 -4.1302112531 26.855641363 192.67226447 -38.804430052
  !5 -159.59154087 83.318680481 13.776631870 206.55133841 -161.82646165 93.626774077
  !6 91.297774084 -33.746922930 -8.6728470368 -355.60235612 -165.20769346 -29.666905585

  !a0i                         a1i                       a2i                      
  a(0,0)=    0.9105631445_dp; a(1,0)=  -0.3084016918_dp; a(2,0)= -0.0906148351_dp; 
  a(0,1)=    0.6361281449_dp; a(1,1)=   0.1860531159_dp; a(2,1)=  0.4527842806_dp; 
  a(0,2)=    2.6861347891_dp; a(1,2)=  -2.5030047259_dp; a(2,2)=  0.5962700728_dp; 
  a(0,3)=  -26.547362491_dp;  a(1,3)=  21.419793629_dp;  a(2,3)= -1.7241829131_dp; 
  a(0,4)=   97.759208784_dp;  a(1,4)= -65.255885330_dp;  a(2,4)= -4.1302112531_dp; 
  a(0,5)= -159.59154087_dp;   a(1,5)=  83.318680481_dp;  a(2,5)= 13.776631870_dp;  
  a(0,6)=   91.297774084_dp;  a(1,6)= -33.746922930_dp;  a(2,6)= -8.6728470368_dp; 

  ! b0i                     b1i                        b2i
  b(0,0)=   0.7240946941_dp; b(1,0)=   -0.5755498075_dp; b(2,0)=   0.0976883116_dp;
  b(0,1)=   2.2382791861_dp; b(1,1)=    0.6995095521_dp; b(2,1)=  -0.2557574982_dp;
  b(0,2)=  -4.0025849485_dp; b(1,2)=    3.8925673390_dp; b(2,2)=  -9.1558561530_dp;
  b(0,3)= -21.003576815_dp;  b(1,3)=  -17.215471648_dp;  b(2,3)=  20.642075974_dp;
  b(0,4)=  26.855641363_dp;  b(1,4)=  192.67226447_dp;   b(2,4)= -38.804430052_dp;
  b(0,5)= 206.55133841_dp;   b(1,5)= -161.82646165_dp;   b(2,5)=  93.626774077_dp;
  b(0,6)= -355.60235612_dp;  b(1,6)= -165.20769346_dp;   b(2,6)= -29.666905585_dp;

  return
  
end subroutine univ_cst

!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! 
!                                        -- Computation of the association terms --
!
!-----------------------------------------------------------------------------------------------------------------------------------
! D. Cordier, Institut UTINAM, Besançon, France.
! 11 avril 2014 -- daniel.cordier@obs-besancon.fr
!
! 16 avril 2020 : - passage de l'criture des valeurs numriques en '1.d0' vers '1._dp'
!                 - implmentation dans le module 'NEW_MOD_PCSAFT'
!
! References: * Abdelkrim Belkadi, Thse de l'Universit de Toulouse, 2008.
!             * Chapman, W. G., Gubbins, K. E., Jackson, G. and Radosz M. (1990), Ind. Eng. Chem. Res., 29, 1709-1721.
!
! Note : la thorie PC-SAFT est applicatble galement  certain fluides comme l'eau, fluides dits "associatifs", ceci moyennant
!        l'ajout de terme au facteur de compressibilit Z et au potentiel chimique.
!
subroutine association_terms(sub,ncall,species_prop,rho,T,Zassoc,muskT_k_assoc)
  use utils_DC
  
  implicit none
  character(len=*), parameter :: subname = 'association_terms'
  
  ! Inputs/outputs of the program:
  character(len=*), intent(in) :: sub
  integer, intent(in)          :: ncall

  type(compprop), dimension(:), intent(in) :: species_prop  ! The properties of studied species.
  real(dp), intent(in)                     :: rho           ! Density (in kg.m^-3)
  real(dp), intent(in)                     :: T             ! Temperature (in K).
  real(dp), intent(out)                    :: Zassoc        ! The compressibility factor (no unit)
  real(dp), dimension(:), intent(out)      :: muskT_k_assoc ! The chemical potential divided by kT (no unit)

  integer  :: i, j, ok
  real(dp) :: rho_molAA3, M_bar, epsi_k, Tr, m, fm, f, sigma
  real(dp), dimension(:), allocatable :: rho_j, d
  real(dp), dimension(:,:), allocatable :: dd, g_seg
  real(dp), dimension(0:3) :: zeta
  real(dp) :: sum_ximi, sum_ximidii, sum_ximidii2, sum_ximidii3

  !---------------------------------------------------------------------------------------------------------------------------------
  ! On liste les proprits des diffrentes espces prises en compte dans le mlange :
  write(6,*) ''
  write(6,*) ' > Hello from association term part !!!'
  write(6,*) ''
  do i= 1, size(species_prop)
     write(6,*) species_prop(i)%name, species_prop(i)%x, species_prop(i)%molmass, species_prop(i)%numseg, species_prop(i)%segdiam, &
                species_prop(i)%segenerg, species_prop(i)%assoflag, species_prop(i)%assenerg, species_prop(i)%assvol
  end do
  write(6,*) ''

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Calcul de la masse molaire moyenne :
  M_bar= 0._dp
  do i= 1, size(species_prop)
     M_bar= M_bar + species_prop(i)%x * species_prop(i)%molmass
  end do
  !write(6,*) ' > M_bar= ', M_bar
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Density en mol.AA^-3 :
  rho_molAA3= rho * M_bar * angstrom_meter**3 
  !write(6,*) ' > rho_molAA3= ', rho_molAA3, rho
  !stop
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Allocation de la mmoire :
  allocate(rho_j(1:size(species_prop)),stat=ok);                       call alloc_error(subname,'rho_j',ok)
  allocate(d(1:size(species_prop)),stat=ok);                           call alloc_error(subname,'d',ok)
  allocate(dd(1:size(species_prop),1:size(species_prop)),stat=ok);     call alloc_error(subname,'dd',ok)
  allocate(g_seg(1:size(species_prop),1:size(species_prop)),stat=ok);  call alloc_error(subname,'g_seg',ok)
  
  zeta= 0._dp ! On met tous les lments de ce vecteur  zro
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Calcul des 'rho_j' des espces (voir Eq. (23) de Chapman et al. 1990 
  do i= 1, size(species_prop)
     rho_j(i)= species_prop(i)%x * rho
  end do

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Calcul des 'd_i' et 'dd(i,j)' :
  do i= 1, size(species_prop)
     epsi_k= species_prop(i)%segenerg
     Tr    = T/epsi_k
     m     = species_prop(i)%numseg
     fm    = 0.0010477_dp + 0.025337_dp * (m-1._dp)/m                         ! Eq. (3) de Chapman et al. 1990 
     f     = (1._dp + 0.2977_dp * Tr) / (1._dp + 0.33163_dp * Tr + fm * Tr**2) ! Eq. (2) de Chapman et al. 1990 
     sigma = species_prop(i)%segdiam
     d(i)  = sigma * f                                                     ! Eq. (1a) de Chapman et al. 1990 
  end do
  
  do i= 1, size(species_prop)
     do j= 1, size(species_prop)
        dd(i,j)= (d(i)+d(j))/2._dp
     end do
  end do
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Calcul des 'zta' (cf. Eq. III-21 p. 67 de Belkadi) :
  
  sum_ximi    = 0._dp
  sum_ximidii = 0._dp
  sum_ximidii2= 0._dp
  sum_ximidii3= 0._dp
  
  do i= 1, size(species_prop)
     sum_ximi    = sum_ximi     + species_prop(i)%x * species_prop(i)%numseg
     sum_ximidii = sum_ximidii  + species_prop(i)%x * species_prop(i)%numseg * dd(i,i)
     sum_ximidii2= sum_ximidii2 + species_prop(i)%x * species_prop(i)%numseg * dd(i,i)**2
     sum_ximidii3= sum_ximidii3 + species_prop(i)%x * species_prop(i)%numseg * dd(i,i)**3
  end do
  
  zeta(0)= pi/6._dp * Navo * rho * sum_ximi     ! en particules / AA^3
  zeta(1)= pi/6._dp * Navo * rho * sum_ximidii  ! en particules / AA^2
  zeta(2)= pi/6._dp * Navo * rho * sum_ximidii2 ! en particules / AA
  zeta(3)= pi/6._dp * Navo * rho * sum_ximidii2 ! en particules (sans dim.)

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Calcul des g^seg_ij :
  
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Calcul des Delta^AiBj :
  do i= 1, size(species_prop)
     do j= 1, size(species_prop)
        g_seg(i,j)= 1._dp/(1._dp-zeta(3)) + 3._dp*dd(i,i)*dd(j,j)/(dd(i,i)+dd(j,j)) * zeta(2)/(1)**2
     end do
  end do
  
  stop

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Deallocation de la mmoire :
  deallocate(rho_j)
  deallocate(d)
  deallocate(dd)
  deallocate(g_seg)
  
  return
  
end subroutine association_terms

!===================================================================================================================================
! ----------------------------------------------------------------------------------------------------------------------------------
subroutine intext_linear(sub,ncall,x,y,x0,y0,dydx0,d2ydx20)
! Interpolation linaire.
! Daniel Cordier, 2 novembre 2011.
! daniel.cordier@obs-besancon.fr
!
! - 7 mai 2012  : mise en arguments optionnels de 'dydx0' et 'd2ydx20'
! - 5 mars 2014 : extension au cas o Nx= 2.
!
  use utils_DC
  
  implicit none
  character(len=*), parameter :: subname = 'intext_linear'
    
  character(len=*), intent(in)        :: sub   ! Name of the calling subroutine
  integer, intent(in)                 :: ncall ! Number of the call in subroutine "sub"
  real(dp), dimension(:), intent(in)  :: x, y
  real(dp), intent(in)                :: x0
  real(dp), intent(out)               :: y0
  real(dp), intent(out), optional     :: dydx0, d2ydx20
  integer :: Nx, Ny
  !real(dp), dimension(:), allocatable :: h, W, U_up, U_diag, U_down, grand_A, grand_B, S

  real(dp) :: extra, a, b
  integer :: i, ix, i0, i1, i2
  logical :: dydx0_pres, d2ydx20_pres
    
  ! Test de la prsence des drives en arguments :
  dydx0_pres   = PRESENT(dydx0)
  d2ydx20_pres = PRESENT(d2ydx20)
  
  Nx = size(x)
  Ny = size(y)

  if (Nx /= Ny) then
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,'(5A)')   ' > Problem in "',subname(1:LEN_TRIM(subname)),'" routine called from "', sub(1:LEN_TRIM(sub)), '"'
     write(6,'(A,I6)') '   call number: ', ncall
     write(6,'(A)')    ' > Nx and Ny are not equal!'
     write(6,*) ' > Nx = ', Nx
     write(6,*) ' > Ny = ', Ny
     stop
  end if  
  if ( Nx < 2 ) then
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,'(5A)')   ' > Problem in "',subname(1:LEN_TRIM(subname)),'" routine called from "', sub(1:LEN_TRIM(sub)), '"'
     write(6,'(A,I6)') '   call number: ', ncall
     write(6,'(A)')    ' > Nx and Ny are <= 2!'
     write(6,'(A,I4)')    '   Nx= ', Nx
     write(6,'(A,I4)')    '   Ny= ', Ny
     stop
  end if

  if ( x0 == x(Nx) ) then
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,'(5A)')   ' > Problem in "',subname(1:LEN_TRIM(subname)),'" routine called from "', sub(1:LEN_TRIM(sub)), '"'
     write(6,'(A,I6)') '   call number: ', ncall
     write(6,'(A)')    ' > x(1)  = x(Nx)!'
     write(6,'(A,ES15.5)') '   x(1)  = ', x(1)
     write(6,'(A,ES15.5)') '   x(Nx) = ', x(Nx) 
     stop
  end if
    
  ! Security:
  if ((x0 .lt. x(1)) .or. (x0 .gt. x(Nx))) then
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,'(5A)')   ' > Problem in "',subname(1:LEN_TRIM(subname)),'" routine called from "', sub(1:LEN_TRIM(sub)), '"'
     write(6,'(A,I6)') '   call number: ', ncall
     write(6,'(A)')    '   x0 has to be in the range x(1)--x(Nx)'
     write(6,'(A10,I5)')     '   Nx   = ', Nx
     write(6,'(A10,ES15.5)') '   x0   = ', x0
     write(6,'(A10,ES15.5)') '   x(1) = ', x(1)
     write(6,'(A10,ES15.5)') '   x(Nx)= ', x(Nx)
     if ( x0 > x(Nx) ) then
        extra= abs(x(Nx)-x0)/abs(x(Nx)-x(1)) * 100._dp
     end if
     if ( x0 < x(1) ) then
        extra= abs(x(1)-x0) /abs(x(Nx)-x(1)) * 100._dp
     end if
     !write(6,'(A,F10.3,A)') '   WARNING: extrapolation (',extra,' %)'
     call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
     write(6,'(A,F5.2,A)') '   WARNING: extrapolation (',extra,' %)'
     write(6,'(A)') ' '
     !stop
  end if
  do i= 1, Nx-1
     if ( x(i) .gt. x(i+1) ) then
        call wf(' > In module: ', 'bright red', mod_name(1:len_trim(mod_name)), 'bright yellow')
        write(6,'(4A)')   ' > Problem in "',subname(1:LEN_TRIM(subname)),'" routine called from "', sub(1:LEN_TRIM(sub)), '"'
        write(6,'(A,I6)') '   call number: ', ncall
        write(6,'(3A)')    ' > In "',subname(1:LEN_TRIM(subname)),'" the x(i)s have to be in the increasing order'
        stop
     end if
  end do
  
  ! On recherche l'emplacement 'x0' dans la grille :
  if ( x0 .lt. x(1) ) then
     ix= 1
  else
     if ( x0 .gt. x(Nx) ) then
        ix= Nx-1
     else
        if ( Nx == 2 ) then
           ix= 1
        else 
           i1= 1; i2= Nx
           do
             i0=(i1+i2)/2
             if (x0 < x(i0) ) then
                i2=i0
             else
                i1=i0
             end if
             if (i1+1 >= i2) exit
           end do
           ix=i1 !
        end if
    end if
  end if

  a= (y(ix+1)-y(ix))/(x(ix+1)-x(ix))
  b= y(ix) - a * x(ix)
  
  y0     = a * x0 + b
  if (dydx0_pres)   dydx0  = a
  if (d2ydx20_pres) d2ydx20= 0._dp
  
  return
      
end subroutine intext_linear

!===================================================================================================================================

end module mod_pcsaft
