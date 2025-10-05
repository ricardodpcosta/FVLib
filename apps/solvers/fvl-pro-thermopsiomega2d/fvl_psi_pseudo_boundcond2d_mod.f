! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: module for fvl_psi_pseudo_boundcond object
! Modification: April, 2024

! #define __fvl_module_prefix__ fvl_psi_pseudo_boundcond
! #define __fvl_module_object__ fvl_psi_pseudo_boundcond

#include "macros.f"

module fvl_psi_pseudo_boundcond2d_mod

use fvl_lib2d

implicit none

!!dir$ attributes align:fvl_align_bytes::fvl_psi_pseudo_boundcond2d
type,extends(fvl_scalarboundcond2d)::fvl_psi_pseudo_boundcond2d
      class(fvl_scalarsinglefield2d),pointer,public::psisinglefield
      class(fvl_scalarfield2d),pointer,public::omegaconvinnerfield
      class(fvl_scalarfield2d),pointer,public::omegadiffinnerfield
      class(fvl_vectormodel2d),pointer,public::forcemodel
contains
      procedure,public::updatefields=>fvl_psi_pseudo_boundcond2d_updatefields
      procedure,public::evalcirculation=>fvl_psi_pseudo_boundcond2d_evalcirculation
end type

interface fvl_psi_constpseudovalue_boundcond2d
      module procedure fvl_psi_constpseudovalue_boundcond2d_init
end interface

interface fvl_psi_constpseudomixed_boundcond2d
      module procedure fvl_psi_constpseudomixed_boundcond2d_init1
      module procedure fvl_psi_constpseudomixed_boundcond2d_init2
end interface

interface fvl_psi_fixedpseudomixed_boundcond2d
      module procedure fvl_psi_fixedpseudomixed_boundcond2d_init1
      module procedure fvl_psi_fixedpseudomixed_boundcond2d_init2
      module procedure fvl_psi_fixedpseudomixed_boundcond2d_init3
end interface

interface fvl_psi_constpseudocombined_boundcond2d
      module procedure fvl_psi_constpseudocombined_boundcond2d_init1
      module procedure fvl_psi_constpseudocombined_boundcond2d_init2
end interface

interface fvl_psi_fixedpseudocombined_boundcond2d
      module procedure fvl_psi_fixedpseudocombined_boundcond2d_init1
      module procedure fvl_psi_fixedpseudocombined_boundcond2d_init2
      module procedure fvl_psi_fixedpseudocombined_boundcond2d_init3
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

function fvl_psi_constpseudovalue_boundcond2d_init(patch,psisinglefield,omegaconvinnerfield,omegadiffinnerfield,forcemodel) result(res)
      class(fvl_patch2d),intent(in),target::patch
      class(fvl_scalarsinglefield2d),target::psisinglefield
      class(fvl_scalarfield2d),target::omegaconvinnerfield
      class(fvl_scalarfield2d),target::omegadiffinnerfield
      class(fvl_vectormodel2d),target::forcemodel
      type(fvl_psi_pseudo_boundcond2d)::res
      res%type=fvl_boundcond2d_value
      res%patch=>patch
      res%psisinglefield=>psisinglefield
      res%omegaconvinnerfield=>omegaconvinnerfield
      res%omegadiffinnerfield=>omegadiffinnerfield
      res%forcemodel=>forcemodel
      allocate(res%valuefield1,source=fvl_const_scalarboundfield2d_init(patch=patch,value=0.0d0))
      res%valuefield2=>null()
      allocate(res%coefsfield,source=fvl_const_darrayboundfield2d_init(patch=patch,value=[1.0d0,0.0d0]))
end function fvl_psi_constpseudovalue_boundcond2d_init

function fvl_psi_constpseudomixed_boundcond2d_init1(patch,value,coefs,psisinglefield,omegaconvinnerfield,omegadiffinnerfield,forcemodel) result(res)
      class(fvl_patch2d),intent(in),target::patch
      real(kind=__fvl_real_kind__),intent(in)::value
      real(kind=__fvl_real_kind__),dimension(1:2),intent(in)::coefs
      class(fvl_scalarsinglefield2d),target::psisinglefield
      class(fvl_scalarfield2d),target::omegaconvinnerfield
      class(fvl_scalarfield2d),target::omegadiffinnerfield
      class(fvl_vectormodel2d),target::forcemodel
      type(fvl_psi_pseudo_boundcond2d)::res
      res%type=fvl_boundcond2d_mixed
      res%patch=>patch
      res%psisinglefield=>psisinglefield
      res%omegaconvinnerfield=>omegaconvinnerfield
      res%omegadiffinnerfield=>omegadiffinnerfield
      res%forcemodel=>forcemodel
      allocate(res%valuefield1,source=fvl_const_scalarboundfield2d_init(patch=patch,value=0.0d0))
      allocate(res%valuefield2,source=fvl_const_scalarboundfield2d_init(patch=patch,value=value))
      allocate(res%coefsfield,source=fvl_const_darrayboundfield2d_init(patch=patch,value=coefs))
end function fvl_psi_constpseudomixed_boundcond2d_init1

function fvl_psi_constpseudomixed_boundcond2d_init2(patch,value,coefsfun,psisinglefield,omegaconvinnerfield,omegadiffinnerfield,forcemodel) result(res)
      class(fvl_patch2d),intent(in),target::patch
      real(kind=__fvl_real_kind__),intent(in)::value
      procedure(fvl_scalarboundcond2d_coefsfun)::coefsfun
      class(fvl_scalarsinglefield2d),target::psisinglefield
      class(fvl_scalarfield2d),target::omegaconvinnerfield
      class(fvl_scalarfield2d),target::omegadiffinnerfield
      class(fvl_vectormodel2d),target::forcemodel
      type(fvl_psi_pseudo_boundcond2d)::res
      res%type=fvl_boundcond2d_mixed
      res%patch=>patch
      res%psisinglefield=>psisinglefield
      res%omegaconvinnerfield=>omegaconvinnerfield
      res%omegadiffinnerfield=>omegadiffinnerfield
      res%forcemodel=>forcemodel
      allocate(res%valuefield1,source=fvl_const_scalarboundfield2d_init(patch=patch,value=0.0d0))
      allocate(res%valuefield2,source=fvl_const_scalarboundfield2d_init(patch=patch,value=value))
      allocate(res%coefsfield,source=fvl_fixed_darrayboundfield2d_init(patch=patch,arraysize=2,valuefun=coefsfun))
end function fvl_psi_constpseudomixed_boundcond2d_init2

function fvl_psi_fixedpseudomixed_boundcond2d_init1(patch,valuefun,coefs,psisinglefield,omegaconvinnerfield,omegadiffinnerfield,forcemodel) result(res)
      class(fvl_patch2d),intent(in),target::patch
      procedure(fvl_scalarboundcond2d_valuefun)::valuefun
      real(kind=__fvl_real_kind__),dimension(1:2),intent(in)::coefs
      class(fvl_scalarsinglefield2d),target::psisinglefield
      class(fvl_scalarfield2d),target::omegaconvinnerfield
      class(fvl_scalarfield2d),target::omegadiffinnerfield
      class(fvl_vectormodel2d),target::forcemodel
      type(fvl_psi_pseudo_boundcond2d)::res
      res%type=fvl_boundcond2d_mixed
      res%patch=>patch
      res%psisinglefield=>psisinglefield
      res%omegaconvinnerfield=>omegaconvinnerfield
      res%omegadiffinnerfield=>omegadiffinnerfield
      res%forcemodel=>forcemodel
      allocate(res%valuefield1,source=fvl_calc_scalarboundfield2d_init1(patch=patch))
      allocate(res%valuefield2,source=fvl_fixed_scalarboundfield2d_init(patch=patch,valuefun=valuefun))
      allocate(res%coefsfield,source=fvl_const_darrayboundfield2d_init(patch=patch,value=coefs))
end function fvl_psi_fixedpseudomixed_boundcond2d_init1

function fvl_psi_fixedpseudomixed_boundcond2d_init2(patch,value,coefsfun,psisinglefield,omegaconvinnerfield,omegadiffinnerfield,forcemodel) result(res)
      class(fvl_patch2d),intent(in),target::patch
      real(kind=__fvl_real_kind__),intent(in)::value
      procedure(fvl_scalarboundcond2d_coefsfun)::coefsfun
      class(fvl_scalarsinglefield2d),target::psisinglefield
      class(fvl_scalarfield2d),target::omegaconvinnerfield
      class(fvl_scalarfield2d),target::omegadiffinnerfield
      class(fvl_vectormodel2d),target::forcemodel
      type(fvl_psi_pseudo_boundcond2d)::res
      res%type=fvl_boundcond2d_mixed
      res%patch=>patch
      res%psisinglefield=>psisinglefield
      res%omegaconvinnerfield=>omegaconvinnerfield
      res%omegadiffinnerfield=>omegadiffinnerfield
      res%forcemodel=>forcemodel
      allocate(res%valuefield1,source=fvl_const_scalarboundfield2d_init(patch=patch,value=0.0d0))
      allocate(res%valuefield2,source=fvl_const_scalarboundfield2d_init(patch=patch,value=value))
      allocate(res%coefsfield,source=fvl_fixed_darrayboundfield2d_init(patch=patch,arraysize=2,valuefun=coefsfun))
end function fvl_psi_fixedpseudomixed_boundcond2d_init2

function fvl_psi_fixedpseudomixed_boundcond2d_init3(patch,valuefun,coefsfun,psisinglefield,omegaconvinnerfield,omegadiffinnerfield,forcemodel) result(res)
      class(fvl_patch2d),intent(in),target::patch
      procedure(fvl_scalarboundcond2d_valuefun)::valuefun
      procedure(fvl_scalarboundcond2d_coefsfun)::coefsfun
      class(fvl_scalarsinglefield2d),target::psisinglefield
      class(fvl_scalarfield2d),target::omegaconvinnerfield
      class(fvl_scalarfield2d),target::omegadiffinnerfield
      class(fvl_vectormodel2d),target::forcemodel
      type(fvl_psi_pseudo_boundcond2d)::res
      res%type=fvl_boundcond2d_mixed
      res%patch=>patch
      res%psisinglefield=>psisinglefield
      res%omegaconvinnerfield=>omegaconvinnerfield
      res%omegadiffinnerfield=>omegadiffinnerfield
      res%forcemodel=>forcemodel
      allocate(res%valuefield1,source=fvl_calc_scalarboundfield2d_init1(patch=patch))
      allocate(res%valuefield2,source=fvl_fixed_scalarboundfield2d_init(patch=patch,valuefun=valuefun))
      allocate(res%coefsfield,source=fvl_fixed_darrayboundfield2d_init(patch=patch,arraysize=2,valuefun=coefsfun))
end function fvl_psi_fixedpseudomixed_boundcond2d_init3

function fvl_psi_constpseudocombined_boundcond2d_init1(patch,value,coefs,psisinglefield,omegaconvinnerfield,omegadiffinnerfield,forcemodel) result(res)
      class(fvl_patch2d),intent(in),target::patch
      real(kind=__fvl_real_kind__),intent(in)::value
      real(kind=__fvl_real_kind__),dimension(1:2),intent(in)::coefs
      class(fvl_scalarsinglefield2d),target::psisinglefield
      class(fvl_scalarfield2d),target::omegaconvinnerfield
      class(fvl_scalarfield2d),target::omegadiffinnerfield
      class(fvl_vectormodel2d),target::forcemodel
      type(fvl_psi_pseudo_boundcond2d)::res
      res%type=fvl_boundcond2d_combined
      res%patch=>patch
      res%psisinglefield=>psisinglefield
      res%omegaconvinnerfield=>omegaconvinnerfield
      res%omegadiffinnerfield=>omegadiffinnerfield
      res%forcemodel=>forcemodel
      allocate(res%valuefield1,source=fvl_const_scalarboundfield2d_init(patch=patch,value=0.0d0))
      allocate(res%valuefield2,source=fvl_const_scalarboundfield2d_init(patch=patch,value=value))
      allocate(res%coefsfield,source=fvl_const_darrayboundfield2d_init(patch=patch,value=coefs))
end function fvl_psi_constpseudocombined_boundcond2d_init1

function fvl_psi_constpseudocombined_boundcond2d_init2(patch,value,coefsfun,psisinglefield,omegaconvinnerfield,omegadiffinnerfield,forcemodel) result(res)
      class(fvl_patch2d),intent(in),target::patch
      real(kind=__fvl_real_kind__),intent(in)::value
      procedure(fvl_scalarboundcond2d_coefsfun)::coefsfun
      class(fvl_scalarsinglefield2d),target::psisinglefield
      class(fvl_scalarfield2d),target::omegaconvinnerfield
      class(fvl_scalarfield2d),target::omegadiffinnerfield
      class(fvl_vectormodel2d),target::forcemodel
      type(fvl_psi_pseudo_boundcond2d)::res
      res%type=fvl_boundcond2d_combined
      res%patch=>patch
      res%psisinglefield=>psisinglefield
      res%omegaconvinnerfield=>omegaconvinnerfield
      res%omegadiffinnerfield=>omegadiffinnerfield
      res%forcemodel=>forcemodel
      allocate(res%valuefield1,source=fvl_const_scalarboundfield2d_init(patch=patch,value=0.0d0))
      allocate(res%valuefield2,source=fvl_const_scalarboundfield2d_init(patch=patch,value=value))
      allocate(res%coefsfield,source=fvl_fixed_darrayboundfield2d_init(patch=patch,arraysize=2,valuefun=coefsfun))
end function fvl_psi_constpseudocombined_boundcond2d_init2

function fvl_psi_fixedpseudocombined_boundcond2d_init1(patch,valuefun,coefs,psisinglefield,omegaconvinnerfield,omegadiffinnerfield,forcemodel) result(res)
      class(fvl_patch2d),intent(in),target::patch
      procedure(fvl_scalarboundcond2d_valuefun)::valuefun
      real(kind=__fvl_real_kind__),dimension(1:2),intent(in)::coefs
      class(fvl_scalarsinglefield2d),target::psisinglefield
      class(fvl_scalarfield2d),target::omegaconvinnerfield
      class(fvl_scalarfield2d),target::omegadiffinnerfield
      class(fvl_vectormodel2d),target::forcemodel
      type(fvl_psi_pseudo_boundcond2d)::res
      res%type=fvl_boundcond2d_combined
      res%patch=>patch
      res%psisinglefield=>psisinglefield
      res%omegaconvinnerfield=>omegaconvinnerfield
      res%omegadiffinnerfield=>omegadiffinnerfield
      res%forcemodel=>forcemodel
      allocate(res%valuefield1,source=fvl_const_scalarboundfield2d_init(patch=patch,value=0.0d0))
      allocate(res%valuefield2,source=fvl_fixed_scalarboundfield2d_init(patch=patch,valuefun=valuefun))
      allocate(res%coefsfield,source=fvl_const_darrayboundfield2d_init(patch=patch,value=coefs))
end function fvl_psi_fixedpseudocombined_boundcond2d_init1

function fvl_psi_fixedpseudocombined_boundcond2d_init2(patch,value,coefsfun,psisinglefield,omegaconvinnerfield,omegadiffinnerfield,forcemodel) result(res)
      class(fvl_patch2d),intent(in),target::patch
      real(kind=__fvl_real_kind__),intent(in)::value
      procedure(fvl_scalarboundcond2d_coefsfun)::coefsfun
      class(fvl_scalarsinglefield2d),target::psisinglefield
      class(fvl_scalarfield2d),target::omegaconvinnerfield
      class(fvl_scalarfield2d),target::omegadiffinnerfield
      class(fvl_vectormodel2d),target::forcemodel
      type(fvl_psi_pseudo_boundcond2d)::res
      res%type=fvl_boundcond2d_combined
      res%patch=>patch
      res%psisinglefield=>psisinglefield
      res%omegaconvinnerfield=>omegaconvinnerfield
      res%omegadiffinnerfield=>omegadiffinnerfield
      res%forcemodel=>forcemodel
      allocate(res%valuefield1,source=fvl_const_scalarboundfield2d_init(patch=patch,value=0.0d0))
      allocate(res%valuefield2,source=fvl_const_scalarboundfield2d_init(patch=patch,value=value))
      allocate(res%coefsfield,source=fvl_fixed_darrayboundfield2d_init(patch=patch,arraysize=2,valuefun=coefsfun))
end function fvl_psi_fixedpseudocombined_boundcond2d_init2

function fvl_psi_fixedpseudocombined_boundcond2d_init3(patch,valuefun,coefsfun,psisinglefield,omegaconvinnerfield,omegadiffinnerfield,forcemodel) result(res)
      class(fvl_patch2d),intent(in),target::patch
      procedure(fvl_scalarboundcond2d_valuefun)::valuefun
      procedure(fvl_scalarboundcond2d_coefsfun)::coefsfun
      class(fvl_scalarsinglefield2d),target::psisinglefield
      class(fvl_scalarfield2d),target::omegaconvinnerfield
      class(fvl_scalarfield2d),target::omegadiffinnerfield
      class(fvl_vectormodel2d),target::forcemodel
      type(fvl_psi_pseudo_boundcond2d)::res
      res%type=fvl_boundcond2d_combined
      res%patch=>patch
      res%psisinglefield=>psisinglefield
      res%omegaconvinnerfield=>omegaconvinnerfield
      res%omegadiffinnerfield=>omegadiffinnerfield
      res%forcemodel=>forcemodel
      allocate(res%valuefield1,source=fvl_const_scalarboundfield2d_init(patch=patch,value=0.0d0))
      allocate(res%valuefield2,source=fvl_fixed_scalarboundfield2d_init(patch=patch,valuefun=valuefun))
      allocate(res%coefsfield,source=fvl_fixed_darrayboundfield2d_init(patch=patch,arraysize=2,valuefun=coefsfun))
end function fvl_psi_fixedpseudocombined_boundcond2d_init3

subroutine fvl_psi_pseudo_boundcond2d_updatefields(this)
      class(fvl_psi_pseudo_boundcond2d),intent(inout)::this
      class(fvl_scalarboundfield2d),pointer::valuefield1,valuefield2
      if(this%type/=fvl_boundcond2d_mixed) then
            select type(valuefield1=>this%valuefield1)
                  class is(fvl_const_scalarboundfield2d)
                        call valuefield1%setvalue(this%psisinglefield%getvalue())
            end select
      else
            select type(valuefield2=>this%valuefield2)
                  class is(fvl_const_scalarboundfield2d)
                        select type(valuefield1=>this%valuefield1)
                              class is(fvl_const_scalarboundfield2d)
                                    call valuefield1%setvalue(this%psisinglefield%getvalue()+valuefield2%getvalue())
                        end select
                  class is(fvl_fixed_scalarboundfield2d)
                        select type(valuefield1=>this%valuefield1)
                              class is(fvl_calc_scalarboundfield2d)
                                    call valuefield1%setscalar(this%psisinglefield%getvalue())
                                    call valuefield1%addfield(this%valuefield2)
                        end select
            end select
      end if
end subroutine fvl_psi_pseudo_boundcond2d_updatefields

function fvl_psi_pseudo_boundcond2d_evalcirculation(this) result(res)
      class(fvl_psi_pseudo_boundcond2d),intent(inout)::this
      real(kind=__fvl_real_kind__)::res
      integer(kind=__fvl_integer_kind__)::i,k
      real(kind=__fvl_real_kind__),dimension(1:2)::force
      class(fvl_patch2d),pointer::patch
      class(fvl_mesh2d),pointer::mesh
      class(fvl_vectorfield2d),pointer::innerfield
#if(__fvl_debug__==1)
      if(.not. this%omegaconvinnerfield%getedgemeanfield()) then
            __fvl_log_error__("Internal error in call fvl_psi_pseudo_boundcond2d_evalcirculation: field not allocated")
      end if
      if(.not. this%omegadiffinnerfield%getedgemeanfield()) then
            __fvl_log_error__("Internal error in call fvl_psi_pseudo_boundcond2d_evalcirculation: field not allocated")
      end if
#endif
      patch=>this%patch
      mesh=>patch%getmesh()
      res=0.0d0
!$omp parallel num_threads(fvl_getnumthreads()) default(shared) shared(this,res,patch,mesh) private(i,k,force,innerfield)
!$omp do schedule(static) reduction(+:res)
      do i=1,patch%getnumedges()
            k=patch%getedgeindex(i)
            innerfield=>this%forcemodel%findinnerfield(mesh%getedgecode(k))
            if(associated(innerfield)) then
                  call innerfield%evaledgemeanvalue(k,force)
                  res=res+(this%omegaconvinnerfield%evaledgemeanvalue(k)+this%omegadiffinnerfield%evaledgemeanvalue(k)&
                        -fvl_dot2(force,mesh%getedgetangent(k)))*mesh%getedgelength(k)
            else
                  res=res+(this%omegaconvinnerfield%evaledgemeanvalue(k)+this%omegadiffinnerfield%evaledgemeanvalue(k))*mesh%getedgelength(k)
            end if
      end do
!$omp end do nowait
!$omp end parallel
end function fvl_psi_pseudo_boundcond2d_evalcirculation

end module fvl_psi_pseudo_boundcond2d_mod
! end of file
