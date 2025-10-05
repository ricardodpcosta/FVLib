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

integer(kind=__fvl_integer_kind__)::psi_numpseudopatches=0
integer(kind=__fvl_integer_kind__)::psi_numpseudoboundconds=0
integer(kind=__fvl_integer_kind__)::psi_numsinglefields=0
type(fvl_patch2d)::psi_pseudo_patch,psi_pseudopatches(100)
type(fvl_scalarsinglefield2d)::psi_singlefields(10)
type(fvl_psi_pseudo_boundcond2d)::psi_pseudoboundconds(10)

abstract interface

__fvl_pure__ function fvl_psi_valuefun(coord,time) result(res)
      real(kind=__fvl_real_kind__),dimension(1:2),intent(in)::coord
      real(kind=__fvl_real_kind__),intent(in)::time
      real(kind=__fvl_real_kind__)::res
end function fvl_psi_valuefun

__fvl_pure__ subroutine fvl_psi_coefsfun(arraysize,coord,time,res)
      integer(kind=__fvl_integer_kind__),intent(in)::arraysize
      real(kind=__fvl_real_kind__),dimension(1:2),intent(in)::coord
      real(kind=__fvl_real_kind__),intent(in)::time
      real(kind=__fvl_real_kind__),dimension(1:arraysize),intent(out)::res
end subroutine fvl_psi_coefsfun

end interface

interface fvl_set_psi_constpseudovalue_boundcond
      procedure fvl_set_psi_constpseudovalue_boundcond1
end interface

interface fvl_set_psi_constpseudomixed_boundcond
      procedure fvl_set_psi_constpseudomixed_boundcond1
      procedure fvl_set_psi_constpseudomixed_boundcond2
end interface

interface fvl_set_psi_fixedpseudomixed_boundcond
      procedure fvl_set_psi_fixedpseudomixed_boundcond1
      procedure fvl_set_psi_fixedpseudomixed_boundcond2
      procedure fvl_set_psi_fixedpseudomixed_boundcond3
end interface

interface fvl_set_psi_constpseudocombined_boundcond
      procedure fvl_set_psi_constpseudocombined_boundcond1
      procedure fvl_set_psi_constpseudocombined_boundcond2
end interface

interface fvl_set_psi_fixedpseudocombined_boundcond
      procedure fvl_set_psi_fixedpseudocombined_boundcond1
      procedure fvl_set_psi_fixedpseudocombined_boundcond2
      procedure fvl_set_psi_fixedpseudocombined_boundcond3
end interface

! end of file
