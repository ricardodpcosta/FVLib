! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: module for fvl_omega_pseudo_boundcond2d object
! Modification: April, 2024

! #define __fvl_module_prefix__ fvl_omega_pseudo_boundcond2d
! #define __fvl_module_object__ fvl_omega_pseudo_boundcond2d

#include "macros.f"

module fvl_omega_pseudo_boundcond2d_mod

use fvl_lib2d

implicit none

!!dir$ attributes align:fvl_align_bytes::fvl_omega_pseudo_boundcond2d
type,extends(fvl_scalarboundcond2d)::fvl_omega_pseudo_boundcond2d
      class(fvl_scalarboundfield2d),pointer,public::staticboundfield
      class(fvl_scalarboundfield2d),pointer,public::psinnboundfield
contains
      procedure,public::updatefields=>fvl_omega_pseudo_boundcond2d_updatefields
end type

interface fvl_omega_constpseudovalue_boundcond2d
      module procedure fvl_omega_constpseudovalue_boundcond2d_init
end interface

interface fvl_omega_fixedpseudovalue_boundcond2d
      module procedure fvl_omega_fixedpseudovalue_boundcond2d_init
end interface

! interface fvl_delete
!       module procedure fvl_scalarboundcond2d_delete
! end interface
!
! interface fvl_clean
!       module procedure fvl_scalarboundcond2d_clean
! end interface
!
! interface fvl_copy
!       module procedure fvl_scalarboundcond2d_copy
! end interface
!
! interface fvl_equal
!       module procedure fvl_scalarboundcond2d_equal
! end interface

! #include "fvl_aos_hdr.f"

contains

! #include "fvl_aos_src.f"

function fvl_omega_constpseudovalue_boundcond2d_init(patch,value,psinnboundfield) result(res)
      class(fvl_patch2d),intent(in),target::patch
      real(kind=__fvl_real_kind__),intent(in)::value
      class(fvl_scalarboundfield2d),target::psinnboundfield
      type(fvl_omega_pseudo_boundcond2d)::res
      res%type=fvl_boundcond2d_value
      res%patch=>patch
      res%psinnboundfield=>psinnboundfield
      allocate(res%valuefield1,source=fvl_calc_scalarboundfield2d_init1(patch=patch))
      res%valuefield2=>null()
      allocate(res%coefsfield,source=fvl_const_darrayboundfield2d_init(patch=patch,value=[1.0d0,0.0d0]))
      allocate(res%staticboundfield,source=fvl_const_scalarboundfield2d_init(patch=patch,value=value))
end function fvl_omega_constpseudovalue_boundcond2d_init

function fvl_omega_fixedpseudovalue_boundcond2d_init(patch,valuefun,psinnboundfield) result(res)
      class(fvl_patch2d),intent(in),target::patch
      procedure(fvl_scalarboundcond2d_valuefun)::valuefun
      class(fvl_scalarboundfield2d),target::psinnboundfield
      type(fvl_omega_pseudo_boundcond2d)::res
      res%type=fvl_boundcond2d_value
      res%patch=>patch
      res%psinnboundfield=>psinnboundfield
      allocate(res%valuefield1,source=fvl_calc_scalarboundfield2d_init1(patch=patch))
      res%valuefield2=>null()
      allocate(res%coefsfield,source=fvl_const_darrayboundfield2d_init(patch=patch,value=[1.0d0,0.0d0]))
      allocate(res%staticboundfield,source=fvl_fixed_scalarboundfield2d_init(patch=patch,valuefun=valuefun))
end function fvl_omega_fixedpseudovalue_boundcond2d_init

subroutine fvl_omega_pseudo_boundcond2d_updatefields(this)
      class(fvl_omega_pseudo_boundcond2d),intent(inout)::this
      integer(kind=__fvl_integer_kind__)::i,j,k
      class(fvl_patch2d),pointer::patch
      class(fvl_mesh2d),pointer::mesh
      class(fvl_scalarboundfield2d),pointer::valuefield
      patch=>this%patch
      mesh=>patch%getmesh()
      select type(valuefield=>this%valuefield1)
            class is(fvl_calc_scalarboundfield2d)
!$omp parallel num_threads(fvl_getnumthreads()) default(shared) shared(this,patch,mesh) private(i,j,k)
!$omp do schedule(static)
                  do i=1,patch%getnumedges()
                        k=patch%getedgeindex(i)
                        do j=1,this%psinnboundfield%getnumboundpointvalues(k)
                              ! calculate boundary condition value
                              call valuefield%setboundpointvalue(j,k,this%staticboundfield%evalboundpointvalue(j,k)&
                                    -this%psinnboundfield%evalboundpointvalue(j,k))
                        end do
                  end do
!$omp end do nowait
!$omp end parallel
      end select
end subroutine fvl_omega_pseudo_boundcond2d_updatefields

end module fvl_omega_pseudo_boundcond2d_mod
! end of file
