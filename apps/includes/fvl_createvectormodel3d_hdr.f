! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: vector model
! Modification: February, 2025

! ============================================================================
! CREATE MODEL
! ============================================================================

type(fvl_vectormodel3d),target::__fvl_print_model__()_models(__fvl_num_steps__)
type(fvl_vectormodel3d),pointer::__fvl_print_model__()_model=>__fvl_print_model__()_models(1)
type(fvl_vectorfield3d_selector)::__fvl_print_model__()_fieldselector
type(fvl_vectorboundcond3d_selector)::__fvl_print_model__()_boundcondselector

abstract interface

__fvl_pure__ subroutine __fvl_print_set_model__()_valuefun(coord,time,res)
      real(kind=__fvl_real_kind__),dimension(1:3),intent(in)::coord
      real(kind=__fvl_real_kind__),intent(in)::time
      real(kind=__fvl_real_kind__),dimension(1:3),intent(out)::res
end subroutine __fvl_print_set_model__()_valuefun

__fvl_pure__ subroutine __fvl_print_set_model__()_coefsfun(arraysize,coord,time,res)
      integer(kind=__fvl_integer_kind__),intent(in)::arraysize
      real(kind=__fvl_real_kind__),dimension(1:3),intent(in)::coord
      real(kind=__fvl_real_kind__),intent(in)::time
      real(kind=__fvl_real_kind__),dimension(1:arraysize),intent(out)::res
end subroutine __fvl_print_set_model__()_coefsfun

end interface

interface __fvl_print_set_model__()_constmixed_boundcond
      procedure __fvl_print_set_model__()_constmixed_boundcond1
      procedure __fvl_print_set_model__()_constmixed_boundcond2
end interface

interface __fvl_print_set_model__()_constcombined_boundcond
      procedure __fvl_print_set_model__()_constcombined_boundcond1
      procedure __fvl_print_set_model__()_constcombined_boundcond2
end interface

interface __fvl_print_set_model__()_fixedmixed_boundcond
      procedure __fvl_print_set_model__()_fixedmixed_boundcond1
      procedure __fvl_print_set_model__()_fixedmixed_boundcond2
end interface

interface __fvl_print_set_model__()_fixedcombined_boundcond
      procedure __fvl_print_set_model__()_fixedcombined_boundcond1
      procedure __fvl_print_set_model__()_fixedcombined_boundcond2
end interface

! end of file
