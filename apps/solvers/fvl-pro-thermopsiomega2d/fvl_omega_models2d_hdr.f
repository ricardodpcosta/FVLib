! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: three-dimensional fluid models
! Modification: April, 2024

! ============================================================================
! VARIABLES
! ============================================================================

integer(kind=__fvl_integer_kind__)::omega_numpseudopatches=0
integer(kind=__fvl_integer_kind__)::omega_numpseudoboundconds=0
type(fvl_patch2d)::omega_pseudo_patch,omega_pseudopatches(100)
type(fvl_omega_pseudo_boundcond2d)::omega_pseudoboundconds(10)

abstract interface

__fvl_pure__ function fvl_omega_valuefun(coord,time) result(res)
      real(kind=__fvl_real_kind__),dimension(1:2),intent(in)::coord
      real(kind=__fvl_real_kind__),intent(in)::time
      real(kind=__fvl_real_kind__)::res
end function fvl_omega_valuefun

__fvl_pure__ subroutine fvl_omega_coefsfun(arraysize,coord,time,res)
      integer(kind=__fvl_integer_kind__),intent(in)::arraysize
      real(kind=__fvl_real_kind__),dimension(1:2),intent(in)::coord
      real(kind=__fvl_real_kind__),intent(in)::time
      real(kind=__fvl_real_kind__),dimension(1:arraysize),intent(out)::res
end subroutine fvl_omega_coefsfun

end interface

! end of file
