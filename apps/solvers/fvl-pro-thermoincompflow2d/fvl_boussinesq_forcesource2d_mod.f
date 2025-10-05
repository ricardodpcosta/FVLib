! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: module for fvl_boussinesq_forcesource2d object
! Modification: February, 2025

! #define __fvl_module_prefix__ fvl_boussinesq_forcesource2d
! #define __fvl_module_object__ fvl_boussinesq_forcesource2d

#include "macros.f"

module fvl_boussinesq_forcesource2d_mod

use fvl_utils_mod
use fvl_calc_vectorfield2d_mod
use fvl_scalarmodel2d_mod
use fvl_vectormodel2d_mod

implicit none

! global variables
character(len=__fvl_character_len__),parameter,public::fvl_boussinesq_forcesource2d_label_default=fvl_char_empty
character(len=__fvl_character_len__),parameter,public::fvl_boussinesq_forcesource2d_filepath_default="."
character(len=__fvl_character_len__),parameter,public::fvl_boussinesq_forcesource2d_fileform_default="binary"
real(kind=__fvl_real_kind__),parameter,public::fvl_boussinesq_forcesource2d_reft_default=0.0d0
real(kind=__fvl_real_kind__),parameter,public::fvl_boussinesq_forcesource2d_refrho_default=0.0d0

!!dir$ attributes align:fvl_align_bytes::fvl_boussinesq_forcesource2d
type,extends(fvl_calc_vectorfield2d)::fvl_boussinesq_forcesource2d
      real(kind=__fvl_real_kind__),public::reft
      real(kind=__fvl_real_kind__),public::refrho
contains
      procedure,public::read1=>fvl_boussinesq_forcesource2d_read1
      procedure,public::read2=>fvl_boussinesq_forcesource2d_read2
      procedure,public::update=>fvl_boussinesq_forcesource2d_update
end type

interface fvl_boussinesq_forcesource2d
      module procedure fvl_boussinesq_forcesource2d_init1
      module procedure fvl_boussinesq_forcesource2d_init2
      module procedure fvl_boussinesq_forcesource2d_allocate
end interface

interface fvl_delete
      module procedure fvl_boussinesq_forcesource2d_delete
end interface

interface fvl_clean
      module procedure fvl_boussinesq_forcesource2d_clean
end interface

interface fvl_copy
      module procedure fvl_boussinesq_forcesource2d_copy
end interface

interface fvl_equal
      module procedure fvl_boussinesq_forcesource2d_equal
end interface

! #include "fvl_aos_hdr.f"

contains

! #include "fvl_aos_src.f"

function fvl_boussinesq_forcesource2d_init1(patch,type,reft,refrho,label,filepath,fileform) result(res)
      class(fvl_patch2d),intent(in),target::patch
      integer(kind=__fvl_integer_kind__),intent(in)::type
      real(kind=__fvl_real_kind__),intent(in)::reft
      real(kind=__fvl_real_kind__),intent(in)::refrho
      character(len=*),intent(in),optional::label
      character(len=*),intent(in),optional::filepath
      character(len=*),intent(in),optional::fileform
      type(fvl_boussinesq_forcesource2d)::res
#if(__fvl_debug__==1)
      if(type/=fvl_field2d_vertcentresfield .and. type/=fvl_field2d_vertboundsfield&
            .and. type/=fvl_field2d_edgecentresfield .and. type/=fvl_field2d_edgemeansfield&
            .and. type/=fvl_field2d_edgequadratsfield .and. type/=fvl_field2d_edgeboundsfield&
            .and. type/=fvl_field2d_cellcentresfield .and. type/=fvl_field2d_cellmeansfield&
            .and. type/=fvl_field2d_cellquadratsfield .and. type/=fvl_field2d_cellboundsfield&
            .and. type/=fvl_field2d_cellboundsfield) then
            __fvl_log_error__("Input error in call fvl_boussinesq_forcesource2d_init1: invalid argument value for <type>")
      end if
      if(((type==fvl_field2d_vertcentresfield .or. type==fvl_field2d_vertboundsfield)&
            .and. (patch%gettype()/=fvl_patch2d_vertspatch .and. patch%gettype()/=fvl_patch2d_boundvertspatch))&
            .or. ((type==fvl_field2d_edgecentresfield .or. type==fvl_field2d_edgemeansfield&
            .or. type==fvl_field2d_edgequadratsfield .or. type==fvl_field2d_edgeboundsfield)&
            .and. (patch%gettype()/=fvl_patch2d_edgespatch .and. patch%gettype()/=fvl_patch2d_boundedgespatch))&
            .or. ((type==fvl_field2d_cellcentresfield .or. type==fvl_field2d_cellmeansfield&
            .or. type==fvl_field2d_cellquadratsfield .or. type==fvl_field2d_cellboundsfield)&
            .and. (patch%gettype()/=fvl_patch2d_cellspatch .and. patch%gettype()/=fvl_patch2d_boundcellspatch))&
            .or. (type==fvl_field2d_boundsfield .and. patch%gettype()/=fvl_patch2d_boundspatch)) then
            __fvl_log_error__("Input error in call fvl_boussinesq_forcesource2d_init1: invalid argument value for <patch>")
      end if
      if(present(fileform)) then
            if(fvl_trim(fileform)/="ascii" .and. fvl_trim(fileform)/="binary") then
                  __fvl_log_error__("Input error in call fvl_boussinesq_forcesource2d_init1: invalid argument value for <fileform>")
            end if
      end if
#endif
      res%patch=>patch
      res%type=type
      res%reft=reft
      res%refrho=refrho
      if(present(label)) then
            res%label=label
      else
            res%label=fvl_boussinesq_forcesource2d_label_default
      end if
      if(present(filepath)) then
            res%filepath=filepath
      else
            res%filepath=fvl_boussinesq_forcesource2d_filepath_default
      end if
      if(present(fileform)) then
            res%fileform=fileform
      else
            res%fileform=fvl_boussinesq_forcesource2d_fileform_default
      end if
      call res%make()
end function fvl_boussinesq_forcesource2d_init1

function fvl_boussinesq_forcesource2d_init2(field,reft,refrho,label,filepath,fileform) result(res)
      class(fvl_vectorfield2d),intent(in)::field
      real(kind=__fvl_real_kind__),intent(in)::reft
      real(kind=__fvl_real_kind__),intent(in)::refrho
      character(len=*),intent(in),optional::label
      character(len=*),intent(in),optional::filepath
      character(len=*),intent(in),optional::fileform
      type(fvl_boussinesq_forcesource2d)::res
#if(__fvl_debug__==1)
      if(present(fileform)) then
            if(fvl_trim(fileform)/="ascii" .and. fvl_trim(fileform)/="binary") then
                  __fvl_log_error__("Input error in call fvl_boussinesq_forcesource2d_init2: invalid argument value for <fileform>")
            end if
      end if
#endif
      res%patch=>field%patch
      res%type=field%type
      res%reft=reft
      res%refrho=refrho
      if(present(label)) then
            res%label=label
      else
            res%label=fvl_boussinesq_forcesource2d_label_default
      end if
      if(present(filepath)) then
            res%filepath=filepath
      else
            res%filepath=fvl_boussinesq_forcesource2d_filepath_default
      end if
      if(present(fileform)) then
            res%fileform=fileform
      else
            res%fileform=fvl_boussinesq_forcesource2d_fileform_default
      end if
      call res%make()
      call res%setfield(field)
end function fvl_boussinesq_forcesource2d_init2

function fvl_boussinesq_forcesource2d_allocate(patch,type,label,dictfile,dictlabel) result(res)
      class(fvl_patch2d),intent(in)::patch
      integer(kind=__fvl_integer_kind__),intent(in)::type
      character(len=*),intent(in)::label
      class(fvl_dict_file),intent(in)::dictfile
      character(len=*),intent(in)::dictlabel
      class(fvl_vectorfield2d),pointer::res
      character(len=__fvl_character_len__)::filepath,fileform
      real(kind=__fvl_real_kind__)::reft,refrho
      logical(kind=__fvl_logical_kind__)::found
#if(__fvl_debug__==1)
      if(type/=fvl_field2d_vertcentresfield .and. type/=fvl_field2d_vertboundsfield&
            .and. type/=fvl_field2d_edgecentresfield .and. type/=fvl_field2d_edgemeansfield&
            .and. type/=fvl_field2d_edgequadratsfield .and. type/=fvl_field2d_edgeboundsfield&
            .and. type/=fvl_field2d_cellcentresfield .and. type/=fvl_field2d_cellmeansfield&
            .and. type/=fvl_field2d_cellquadratsfield .and. type/=fvl_field2d_cellboundsfield&
            .and. type/=fvl_field2d_cellboundsfield) then
            __fvl_log_error__("Input error in call fvl_boussinesq_forcesource2d_allocate: invalid argument value for <type>")
      end if
      if(((type==fvl_field2d_vertcentresfield .or. type==fvl_field2d_vertboundsfield)&
            .and. (patch%gettype()/=fvl_patch2d_vertspatch .and. patch%gettype()/=fvl_patch2d_boundvertspatch))&
            .or. ((type==fvl_field2d_edgecentresfield .or. type==fvl_field2d_edgemeansfield&
            .or. type==fvl_field2d_edgequadratsfield .or. type==fvl_field2d_edgeboundsfield)&
            .and. (patch%gettype()/=fvl_patch2d_edgespatch .and. patch%gettype()/=fvl_patch2d_boundedgespatch))&
            .or. ((type==fvl_field2d_cellcentresfield .or. type==fvl_field2d_cellmeansfield&
            .or. type==fvl_field2d_cellquadratsfield .or. type==fvl_field2d_cellboundsfield)&
            .and. (patch%gettype()/=fvl_patch2d_cellspatch .and. patch%gettype()/=fvl_patch2d_boundcellspatch))&
            .or. (type==fvl_field2d_boundsfield .and. patch%gettype()/=fvl_patch2d_boundspatch)) then
            __fvl_log_error__("Input error in call fvl_boussinesq_forcesource2d_allocate: invalid argument value for <patch>")
      end if
#endif
      filepath=dictfile%getvalue(dictlabel,"file_path",fvl_boussinesq_forcesource2d_filepath_default,found)
      fileform=dictfile%getvalue(dictlabel,"file_form",fvl_boussinesq_forcesource2d_fileform_default,found)
      reft=dictfile%getvalue(dictlabel,"ref_t",fvl_boussinesq_forcesource2d_reft_default)
      refrho=dictfile%getvalue(dictlabel,"ref_rho",fvl_boussinesq_forcesource2d_refrho_default)
      allocate(res,source=fvl_boussinesq_forcesource2d_init1(patch,type,reft,refrho,fvl_trim(label)//"."//fvl_trim(dictlabel),filepath,fileform))
end function fvl_boussinesq_forcesource2d_allocate

subroutine fvl_boussinesq_forcesource2d_delete(this)
      type(fvl_boussinesq_forcesource2d),intent(inout)::this
      this%label=fvl_char_empty
      this%filepath=fvl_char_empty
      this%fileform=fvl_char_empty
      this%type=0
      this%size=0
      this%numindices=0
      this%maxnumvalues=0
      this%reft=0.0d0
      this%refrho=0.0d0
      call fvl_deallocate(this%numvalues)
      call fvl_deallocate(this%values1)
      call fvl_deallocate(this%values1)
      this%patch=>null()
end subroutine fvl_boussinesq_forcesource2d_delete

subroutine fvl_boussinesq_forcesource2d_clean(this)
      type(fvl_boussinesq_forcesource2d),intent(inout)::this
      this%label=fvl_char_empty
      this%filepath=fvl_char_empty
      this%fileform=fvl_char_empty
      this%type=0
      this%size=0
      this%numindices=0
      this%maxnumvalues=0
      this%reft=0.0d0
      this%refrho=0.0d0
      call fvl_omp_clean(this%numvalues)
      call fvl_omp_clean(this%values1)
      call fvl_omp_clean(this%values2)
      this%patch=>null()
end subroutine fvl_boussinesq_forcesource2d_clean

subroutine fvl_boussinesq_forcesource2d_copy(this,that)
      type(fvl_boussinesq_forcesource2d),intent(in)::this
      type(fvl_boussinesq_forcesource2d),intent(out)::that
      that%label=this%label
      that%filepath=this%filepath
      that%fileform=this%fileform
      that%type=this%type
      that%size=this%size
      that%numindices=this%numindices
      that%maxnumvalues=this%maxnumvalues
      that%reft=this%reft
      that%refrho=this%refrho
      call fvl_omp_copy(this%numvalues,that%numvalues)
      call fvl_omp_copy(this%values1,that%values1)
      call fvl_omp_copy(this%values2,that%values2)
      that%patch=>this%patch
end subroutine fvl_boussinesq_forcesource2d_copy

function fvl_boussinesq_forcesource2d_equal(this,that) result(res)
      type(fvl_boussinesq_forcesource2d),intent(in)::this
      type(fvl_boussinesq_forcesource2d),intent(in)::that
      logical(kind=__fvl_logical_kind__)::res
      if(that%label/=this%label .or.&
            that%filepath/=this%filepath .or.&
            that%fileform/=this%fileform .or.&
            that%type/=this%type .or.&
            that%size/=this%size .or.&
            that%numindices/=this%numindices .or.&
            that%maxnumvalues/=this%maxnumvalues .or.&
            .not. fvl_equal(this%numvalues,that%numvalues) .or.&
            .not. fvl_equal(that%reft,this%reft) .or.&
            .not. fvl_equal(that%refrho,this%refrho) .or.&
            .not. fvl_equal(this%values1,that%values1) .or.&
            .not. fvl_equal(this%values2,that%values2) .or.&
            .not. associated(this%patch,that%patch)) then
            res=.false.
      else
            res=.true.
      end if
end function fvl_boussinesq_forcesource2d_equal

subroutine fvl_boussinesq_forcesource2d_read1(this)
      class(fvl_boussinesq_forcesource2d),intent(inout)::this
end subroutine fvl_boussinesq_forcesource2d_read1

subroutine fvl_boussinesq_forcesource2d_read2(this,file)
      class(fvl_boussinesq_forcesource2d),intent(inout)::this
      class(fvl_data_file),intent(inout)::file
end subroutine fvl_boussinesq_forcesource2d_read2

subroutine fvl_boussinesq_forcesource2d_update(this,tmodel,gmodel,betamodel)
      class(fvl_boussinesq_forcesource2d),intent(inout)::this
      class(fvl_scalarmodel2d),intent(in)::tmodel
      class(fvl_vectormodel2d),intent(in)::gmodel
      class(fvl_scalarmodel2d),intent(in)::betamodel
      class(fvl_scalarfield2d),pointer::tfield
      class(fvl_vectorfield2d),pointer::gfield
      class(fvl_scalarfield2d),pointer::betafield
      tfield=>tmodel%findfield(this%patch)
      gfield=>gmodel%findfield(this%patch)
      betafield=>betamodel%findfield(this%patch)
      call this%setfield(tfield)
      call this%subscalar(this%reft)
      call this%mulfield(gfield)
      call this%mulfield(betafield)
      call this%mulscalar(-this%refrho)
end subroutine fvl_boussinesq_forcesource2d_update

end module fvl_boussinesq_forcesource2d_mod
! end of file
