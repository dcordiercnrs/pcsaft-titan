! ----------------------------------------------------------------------------------------------------------------------------------
! D. Cordier -- 17 avril 2020 -- Module d'utilitaires, destiné à remplacer 'nrtype.f90' et 'foul.f90', en y ajoutant quelques
!                                fonctionnalités.
!
! - 16/04/2020 (ver. 0.8.0) : mise en place en agrégeant divers éléments utilisés ça et là jusqu'alors, utilisation avec la 
!          nouvelle version modulaire de mon implémentation de PC-SAFT.
! - 28/04/2020 (ver. 0.8.1) : mise en place de 'real_to_string' du module 'foul.f90'
! -  4/05/2020 (ver. 0.8.2) : we introduce the variable type of integer, plus 'i2s'
!
module utils_dc
  
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options
  use foul, only :  wf => write_formatted   ! Module pour avoir de jolies polices de couleurs dans le terminal
  use foul, only :  i2s => integer_to_string
  use foul, only :  r2s => real_to_string
  implicit none
  
  character(len=*), private, parameter :: mod_name = 'utils_dc'
  character(len=*), private, parameter :: version  = '0.8.2'
    
  integer,  parameter :: si = 4 ! Integer type in use.
  integer,  parameter :: dp = selected_real_kind(p=14,r=300) ! ensures that the floats have a precision of at least "p" significant decimals 
                                                             ! and an exponent range of at least 10^-r to 10+r
                                                             ! cf. "Modern Fortran Explained - Incorporating Fortran 2018" 
                                                             ! by Michael Metcalf, John Reid, and Malcolm Cohen, 2018, 
                                                             ! Oxford University Press,
                                                             ! p. 15.
                                                             
  real(dp), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_dp

  public :: display_compil_numpres, alloc_error, concapath
  
! ==================================================================================================================================
CONTAINS

! ----------------------------------------------------------------------------------------------------------------------------------
subroutine display_compil_numpres (x)
  implicit none
  real(dp), intent(in) :: x
  logical              :: underflow_support, gradual
  character(len=100)   :: temp
  
  call wf(' > This file was compiled by: ','bright red', compiler_version(), 'bright yellow')
  call wf('   using the options: ', 'bright red')
  call wf('   ', 'red', compiler_options(), 'bright yellow')

  underflow_support = ieee_support_underflow_control(x)
  
  write(6,*) '  - support for underflow control --: ', underflow_support
  CALL IEEE_GET_UNDERFLOW_MODE (gradual) ! Stores underflow mode
  write(6,*) '  - gradual underflow control ------: ', gradual
  
  write(6,*) ''
  call wf('   Numerical precision: ----------------------------------', 'bright red')
  
  write(temp,*) epsilon(x)
  call wf('    - epsilon(real) .......: ', 'bright red', temp(1:len_trim(temp)), 'bright yellow')
  
  write(temp,*) huge(x)
  call wf('    - huge(real) ..........: ', 'bright red', temp(1:len_trim(temp)), 'bright yellow')
  
  write(temp,*) maxexponent(x)
  call wf('    - maxexponent(real) ...: ', 'bright red', temp(1:len_trim(temp)), 'bright yellow')
  
  write(temp,*) minexponent(x)
  call wf('    - minexponent(real) ...: ', 'bright red', temp(1:len_trim(temp)), 'bright yellow')
  
  write(temp,*) precision(x)
  call wf('    - precision(real) .....: ', 'bright red', temp(1:len_trim(temp)), 'bright yellow')
  call wf('   -------------------------------------------------------', 'bright red')
  
  return
  
end subroutine 

! ----------------------------------------------------------------------------------------------------------------------------------
subroutine alloc_error(sub,tab,nstat)
! D. Cordier, 27 octobre 2011.
! Message d'erreur en cas de non-allocation de mémoire
! pour un tableau.

  implicit none
  character(len=*), parameter :: subname = 'alloc_error'
  character(len=*), intent(in) :: sub ! Name of the calling subroutine
  character(len=*), intent(in) :: tab ! Name of the array
  integer, intent(in)          :: nstat
  character(:),allocatable :: path
  
  path= concapath(sub,subname)
 
  if ( nstat /= 0 ) then
  write(6,'(3A)')   ' > Problem in subroutine "',subname,'"'
  write(6,'(3A)')   '   call path "', path,'":'
  write(6,'(3A)')   '     - the array "',tab(1:LEN_TRIM(tab)),'" cannot be allocated'
  write(6,'(A,I4)') '     - nstat= ', nstat
  write(6,'(A)')    ''
  stop
  end if
  
end subroutine alloc_error

! ----------------------------------------------------------------------------------------------------------------------------------
function concapath(sub1,sub2) result(res)
! D. Cordier, GSMA, 5 juin 2015.
! Cette function est à l'origine 'concapath.f08' que j'avais mise au point pour ma tentative de travail sur la structure verticale
! des lacs de Titan (~/work_2015/lakes_vertical_structure/).
!
  implicit none
  character(:),allocatable :: res
  character(len=*), intent(in)  :: sub1
  character(len=*), intent(in)  :: sub2
  !character(len=:), intent(out) :: concapath
  !real :: concapath
  res= sub1(1:len_trim(sub1))//'/'//sub2(1:len_trim(sub2))
  
  return
  
end function concapath
! ----------------------------------------------------------------------------------------------------------------------------------

end module utils_DC
