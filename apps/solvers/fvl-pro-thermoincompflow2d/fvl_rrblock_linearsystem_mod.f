! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: module for fvl_rblock_linearsystem object
! Modification: March, 2025

#define __fvl_module_prefix__ fvl_rblock_linearsystem
#define __fvl_module_object__ fvl_rblock_linearsystem

#include "macros.f"

module fvl_rblock_linearsystem_mod

use fvl_utils_mod
use fvl_runtime_mod
use fvl_ptr_runtime_mod
use fvl_model_mod
use fvl_ptr_model_mod
use fvl_model_mod
use fvl_ptr_model_mod
use fvl_linearsystem_mod
use fvl_linearsolver_mod

implicit none

! global variables
integer(kind=__fvl_integer_kind__),parameter,public::fvl_rblock_linearsystem_none=0
integer(kind=__fvl_integer_kind__),parameter,public::fvl_rblock_linearsystem_coupled=1
integer(kind=__fvl_integer_kind__),parameter,public::fvl_rblock_linearsystem_segregated=2
integer(kind=__fvl_integer_kind__),parameter,public::fvl_rblock_linearsystem_nopredictor=0
integer(kind=__fvl_integer_kind__),parameter,public::fvl_rblock_linearsystem_taylorexpansion=1
integer(kind=__fvl_integer_kind__),parameter,public::fvl_rblock_linearsystem_forwardeuler=2
integer(kind=__fvl_integer_kind__),parameter,public::fvl_rblock_linearsystem_method_default=1
integer(kind=__fvl_integer_kind__),parameter,public::fvl_rblock_linearsystem_predictor_default=1
integer(kind=__fvl_integer_kind__),parameter,public::fvl_rblock_linearsystem_numlevels_default=2
integer(kind=__fvl_integer_kind__),parameter,public::fvl_rblock_linearsystem_numsteps_default=1

!!dir$ attributes align:__fvl_align_bytes__::fvl_rblock_linearsystem
type,extends(fvl_linearsystem)::fvl_rblock_linearsystem
      integer(kind=__fvl_integer_kind__),public::method
      integer(kind=__fvl_integer_kind__),public::predictor
      integer(kind=__fvl_integer_kind__),public::basesize
      real(kind=__fvl_real_kind__),public::currenttimestep
      procedure(fvl_rblock_linearsystem_tfun),nopass,pointer,public::tfun1
      procedure(fvl_rblock_linearsystem_tfun),nopass,pointer,public::tfun2
      procedure(fvl_rblock_linearsystem_tfun),nopass,pointer,public::tfun3
      procedure(fvl_rblock_linearsystem_tfun),nopass,pointer,public::tfun4
      procedure(fvl_rblock_linearsystem_afun),nopass,pointer,public::afun1
      procedure(fvl_rblock_linearsystem_afun),nopass,pointer,public::afun2
      procedure(fvl_rblock_linearsystem_afun),nopass,pointer,public::afun3
      procedure(fvl_rblock_linearsystem_afun),nopass,pointer,public::afun4
      procedure(fvl_rblock_linearsystem_bfun),nopass,pointer,public::bfun1
      procedure(fvl_rblock_linearsystem_bfun),nopass,pointer,public::bfun2
      procedure(fvl_rblock_linearsystem_bfun),nopass,pointer,public::bfun3
      procedure(fvl_rblock_linearsystem_bfun),nopass,pointer,public::bfun4
      procedure(fvl_rblock_linearsystem_pfun),nopass,pointer,public::pfun1
      procedure(fvl_rblock_linearsystem_pfun),nopass,pointer,public::pfun2
      procedure(fvl_rblock_linearsystem_pfun),nopass,pointer,public::pfun3
      procedure(fvl_rblock_linearsystem_pfun),nopass,pointer,public::pfun4
      procedure(fvl_rblock_linearsystem_ufun),nopass,pointer,public::ufun1
      procedure(fvl_rblock_linearsystem_ufun),nopass,pointer,public::ufun2
      procedure(fvl_rblock_linearsystem_ufun),nopass,pointer,public::ufun3
      procedure(fvl_rblock_linearsystem_ufun),nopass,pointer,public::ufun4
      real(kind=__fvl_real_kind__),allocatable,dimension(:,:),public::rhs
      real(kind=__fvl_real_kind__),allocatable,dimension(:,:,:),public::basecoefs
      real(kind=__fvl_real_kind__),allocatable,dimension(:,:),public::basecoefsrhs
      class(fvl_linearsolver),pointer::linearsolver
contains
      procedure,public::gettfun1=>fvl_rblock_linearsystem_gettfun1
      procedure,public::getafun1=>fvl_rblock_linearsystem_getafun1
      procedure,public::getbfun1=>fvl_rblock_linearsystem_getbfun1
      procedure,public::getpfun1=>fvl_rblock_linearsystem_getpfun1
      procedure,public::getufun1=>fvl_rblock_linearsystem_getufun1
      procedure,public::gettfun2=>fvl_rblock_linearsystem_gettfun2
      procedure,public::getafun2=>fvl_rblock_linearsystem_getafun2
      procedure,public::getbfun2=>fvl_rblock_linearsystem_getbfun2
      procedure,public::getpfun2=>fvl_rblock_linearsystem_getpfun2
      procedure,public::getufun2=>fvl_rblock_linearsystem_getufun2
      procedure,public::gettfun3=>fvl_rblock_linearsystem_gettfun3
      procedure,public::getafun3=>fvl_rblock_linearsystem_getafun3
      procedure,public::getbfun3=>fvl_rblock_linearsystem_getbfun3
      procedure,public::getpfun3=>fvl_rblock_linearsystem_getpfun3
      procedure,public::getufun3=>fvl_rblock_linearsystem_getufun3
      procedure,public::gettfun4=>fvl_rblock_linearsystem_gettfun4
      procedure,public::getafun4=>fvl_rblock_linearsystem_getafun4
      procedure,public::getbfun4=>fvl_rblock_linearsystem_getbfun4
      procedure,public::getpfun4=>fvl_rblock_linearsystem_getpfun4
      procedure,public::getufun4=>fvl_rblock_linearsystem_getufun4
      procedure,public::evallhs=>fvl_rblock_linearsystem_evallhs
      procedure,public::evalrhs=>fvl_rblock_linearsystem_evalrhs
      procedure,public::precond=>fvl_rblock_linearsystem_precond
      procedure,public::update=>fvl_rblock_linearsystem_update
      procedure,public::make=>fvl_rblock_linearsystem_make
      procedure,public::solve=>fvl_rblock_linearsystem_solve
      procedure,public::solvecoupled=>fvl_rblock_linearsystem_solvecoupled
      procedure,public::solvesegregated=>fvl_rblock_linearsystem_solvesegregated
end type

interface fvl_rblock_linearsystem
      module procedure fvl_rblock_linearsystem_init1
      module procedure fvl_rblock_linearsystem_init2
      module procedure fvl_rblock_linearsystem_init3
      module procedure fvl_rblock_linearsystem_init4
      module procedure fvl_rblock_linearsystem_init5
      module procedure fvl_rblock_linearsystem_init6
end interface

interface fvl_delete
      module procedure fvl_rblock_linearsystem_delete
end interface

interface fvl_clean
      module procedure fvl_rblock_linearsystem_clean
end interface

interface fvl_copy
      module procedure fvl_rblock_linearsystem_copy
end interface

interface fvl_equal
      module procedure fvl_rblock_linearsystem_equal
end interface

abstract interface

subroutine fvl_rblock_linearsystem_tfun(step,x,res)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:size(x,1)),intent(out)::res
end subroutine fvl_rblock_linearsystem_tfun

subroutine fvl_rblock_linearsystem_afun(step,x,res)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:size(x,1)),intent(out)::res
end subroutine fvl_rblock_linearsystem_afun

subroutine fvl_rblock_linearsystem_bfun(step,res)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(inout)::res
end subroutine fvl_rblock_linearsystem_bfun

subroutine fvl_rblock_linearsystem_pfun(step,x,res)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:size(x,1)),intent(out)::res
end subroutine fvl_rblock_linearsystem_pfun

subroutine fvl_rblock_linearsystem_ufun(step)
      integer(kind=__fvl_integer_kind__),intent(in)::step
end subroutine fvl_rblock_linearsystem_ufun

end interface

#include "fvl_aos_hdr.f"

contains

#include "fvl_aos_src.f"

function fvl_rblock_linearsystem_init1(runtimes,models1,models2,models3,models4,models5,&
      tfun1,afun1,bfun1,pfun1,ufun1,tfun2,afun2,bfun2,pfun2,ufun2,tfun3,afun3,bfun3,pfun3,ufun3,tfun4,afun4,bfun4,pfun4,ufun4,&
      linearsolver,method,predictor,numsteps,numlevels) result(res)
      class(fvl_runtime),dimension(1:),contiguous,intent(in)::runtimes
      class(fvl_model),dimension(1:,1:),contiguous,intent(in),target::models1
      class(fvl_model),dimension(1:,1:),contiguous,intent(in),target::models2
      class(fvl_model),dimension(1:,1:),contiguous,intent(in),optional,target::models3
      class(fvl_model),dimension(1:,1:),contiguous,intent(in),optional,target::models4
      class(fvl_model),dimension(1:,1:),contiguous,intent(in),optional,target::models5
      procedure(fvl_rblock_linearsystem_tfun)::tfun1
      procedure(fvl_rblock_linearsystem_afun)::afun1
      procedure(fvl_rblock_linearsystem_bfun)::bfun1
      procedure(fvl_rblock_linearsystem_pfun)::pfun1
      procedure(fvl_rblock_linearsystem_ufun)::ufun1
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun2
      procedure(fvl_rblock_linearsystem_afun),optional::afun2
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun2
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun2
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun2
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun3
      procedure(fvl_rblock_linearsystem_afun),optional::afun3
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun3
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun3
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun3
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun4
      procedure(fvl_rblock_linearsystem_afun),optional::afun4
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun4
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun4
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun4
      class(fvl_linearsolver),intent(in),target::linearsolver
      integer(kind=__fvl_integer_kind__),intent(in),optional::method
      integer(kind=__fvl_integer_kind__),intent(in),optional::predictor
      integer(kind=__fvl_integer_kind__),intent(in),optional::numlevels
      integer(kind=__fvl_integer_kind__),intent(in),optional::numsteps
      type(fvl_rblock_linearsystem)::res
      integer(kind=__fvl_integer_kind__)::i,j,k,aux
      class(fvl_model),dimension(:,:),pointer::models
#if(__fvl_debug__==1)
      if(present(method)) then
            if(method/=fvl_rblock_linearsystem_none .and. method/=fvl_rblock_linearsystem_coupled .and.&
                  method/=fvl_rblock_linearsystem_segregated) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: invalid argument value for <method>")
            end if
      end if
      if(present(predictor)) then
            if(predictor/=fvl_rblock_linearsystem_nopredictor .and. predictor/=fvl_rblock_linearsystem_taylorexpansion .and.&
                  predictor/=fvl_rblock_linearsystem_forwardeuler) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: invalid argument value for <predictor>")
            end if
      end if
      if(present(numsteps)) then
            if(numsteps<1) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: invalid argument value for <numsteps>")
            end if
      end if
      if(present(numlevels)) then
            if(numlevels<2) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: invalid argument value for <numlevels>")
            end if
      end if
      if(present(numlevels)) then
            if(numlevels==3 .and. (.not. present(models3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <models3>")
            end if
            if(numlevels==3 .and. (.not. present(tfun2))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <tfun2>")
            end if
            if(numlevels==3 .and. (.not. present(afun2))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <afun2>")
            end if
            if(numlevels==3 .and. (.not. present(bfun2))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <bfun2>")
            end if
            if(numlevels==3 .and. (.not. present(pfun2))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <pfun2>")
            end if
            if(numlevels==3 .and. (.not. present(ufun2))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <ufun2>")
            end if
            if(numlevels==4 .and. (.not. present(models4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <models4>")
            end if
            if(numlevels==4 .and. (.not. present(tfun3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <tfun3>")
            end if
            if(numlevels==4 .and. (.not. present(afun3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <afun3>")
            end if
            if(numlevels==4 .and. (.not. present(bfun3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <bfun3>")
            end if
            if(numlevels==4 .and. (.not. present(pfun3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <pfun3>")
            end if
            if(numlevels==4 .and. (.not. present(ufun3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <ufun3>")
            end if
            if(numlevels==5 .and. (.not. present(models5))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <models5>")
            end if
            if(numlevels==5 .and. (.not. present(tfun4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <tfun4>")
            end if
            if(numlevels==5 .and. (.not. present(afun4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <afun4>")
            end if
            if(numlevels==5 .and. (.not. present(bfun4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <bfun4>")
            end if
            if(numlevels==5 .and. (.not. present(pfun4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <pfun4>")
            end if
            if(numlevels==5 .and. (.not. present(ufun4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init1: missing argument <ufun4>")
            end if
      end if
#endif
      ! initialize data
      res%numruntimes=size(runtimes,1)
      call fvl_allocate(res%runtimes,res%numruntimes)
      do i=1,res%numruntimes
            call res%runtimes(i)%setptr(runtimes(i))
      end do
      res%tfun1=>tfun1
      res%afun1=>afun1
      res%bfun1=>bfun1
      res%pfun1=>pfun1
      res%ufun1=>ufun1
      if(present(tfun2)) then
            res%tfun2=>tfun2
      end if
      if(present(afun2)) then
            res%afun2=>afun2
      end if
      if(present(bfun2)) then
            res%bfun2=>bfun2
      end if
      if(present(pfun2)) then
            res%pfun2=>pfun2
      end if
      if(present(ufun2)) then
            res%ufun2=>ufun2
      end if
      if(present(tfun3)) then
            res%tfun3=>tfun3
      end if
      if(present(afun3)) then
            res%afun3=>afun3
      end if
      if(present(bfun3)) then
            res%bfun3=>bfun3
      end if
      if(present(pfun3)) then
            res%pfun3=>pfun3
      end if
      if(present(ufun3)) then
            res%ufun3=>ufun3
      end if
      if(present(tfun4)) then
            res%tfun4=>tfun4
      end if
      if(present(afun4)) then
            res%afun4=>afun4
      end if
      if(present(bfun4)) then
            res%bfun4=>bfun4
      end if
      if(present(pfun4)) then
            res%pfun4=>pfun4
      end if
      if(present(ufun4)) then
            res%ufun4=>ufun4
      end if
      res%linearsolver=>linearsolver
      if(present(method)) then
            res%method=method
      else
            res%method=fvl_rblock_linearsystem_method_default
      end if
      if(present(predictor)) then
            res%predictor=predictor
      else
            res%predictor=fvl_rblock_linearsystem_predictor_default
      end if
      if(present(numlevels)) then
            res%numlevels=numlevels
      else
            res%numlevels=fvl_rblock_linearsystem_numlevels_default
      end if
      if(present(numsteps)) then
            res%numsteps=numsteps
      else
            res%numsteps=fvl_rblock_linearsystem_numsteps_default
      end if
      ! initialise data
      res%nummodels=0
      do i=1,res%numsteps
            do j=1,res%numlevels
                  select case(j)
                        case(1)
                              models=>models1
                        case(2)
                              models=>models2
                        case(3)
                              models=>models3
                        case(4)
                              models=>models4
                        case(5)
                              models=>models5
                  end select
                  res%nummodels=res%nummodels+size(models,1)
            end do
      end do
      call fvl_allocate(res%models,res%nummodels,res%numlevels,res%numsteps)
      res%pointsize=0
      res%datasize=0
      aux=0
      do i=1,res%numsteps
            do j=1,res%numlevels
                  select case(j)
                        case(1)
                              models=>models1
                        case(2)
                              models=>models2
                        case(3)
                              models=>models3
                        case(4)
                              models=>models4
                        case(5)
                              models=>models5
                  end select
                  do k=1,size(models,1)
                        call res%models(k,j,i)%setptr(models(k,i))
                        res%datasize=res%datasize+models(k,i)%getfieldssize()
                        if(i==1 .and. j==1) then
                              res%pointsize=res%pointsize+models(k,i)%getfieldssize()
                        else
                              aux=aux+models(k,i)%getfieldssize()
                        end if
                  end do
#if(__fvl_debug__==1)
                  if(i/=1 .and. j/=1) then
                        if(aux/=res%pointsize) then
                              __fvl_log_error__("Internal error in call fvl_rblock_linearsystem_init1: inconsistent model sizes")
                        end if
                  end if
                  aux=0
#endif
            end do
      end do
      call res%make()
end function fvl_rblock_linearsystem_init1

function fvl_rblock_linearsystem_init2(runtimes,models1,models2,models3,models4,models5,&
      tfun1,afun1,bfun1,pfun1,ufun1,tfun2,afun2,bfun2,pfun2,ufun2,tfun3,afun3,bfun3,pfun3,ufun3,tfun4,afun4,bfun4,pfun4,ufun4,&
      linearsolver,method,predictor,numlevels,numsteps) result(res)
      class(fvl_ptr_runtime),dimension(1:),contiguous,intent(in)::runtimes
      class(fvl_ptr_model),dimension(1:,1:),contiguous,intent(in),target::models1
      class(fvl_ptr_model),dimension(1:,1:),contiguous,intent(in),target::models2
      class(fvl_ptr_model),dimension(1:,1:),contiguous,intent(in),optional,target::models3
      class(fvl_ptr_model),dimension(1:,1:),contiguous,intent(in),optional,target::models4
      class(fvl_ptr_model),dimension(1:,1:),contiguous,intent(in),optional,target::models5
      procedure(fvl_rblock_linearsystem_tfun)::tfun1
      procedure(fvl_rblock_linearsystem_afun)::afun1
      procedure(fvl_rblock_linearsystem_bfun)::bfun1
      procedure(fvl_rblock_linearsystem_pfun)::pfun1
      procedure(fvl_rblock_linearsystem_ufun)::ufun1
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun2
      procedure(fvl_rblock_linearsystem_afun),optional::afun2
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun2
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun2
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun2
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun3
      procedure(fvl_rblock_linearsystem_afun),optional::afun3
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun3
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun3
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun3
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun4
      procedure(fvl_rblock_linearsystem_afun),optional::afun4
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun4
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun4
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun4
      class(fvl_linearsolver),intent(in),target::linearsolver
      integer(kind=__fvl_integer_kind__),intent(in),optional::method
      integer(kind=__fvl_integer_kind__),intent(in),optional::predictor
      integer(kind=__fvl_integer_kind__),intent(in),optional::numlevels
      integer(kind=__fvl_integer_kind__),intent(in),optional::numsteps
      type(fvl_rblock_linearsystem)::res
      integer(kind=__fvl_integer_kind__)::i,j,k,aux
      class(fvl_runtime),pointer::runtime
      class(fvl_model),pointer::model
      class(fvl_ptr_model),dimension(:,:),pointer::models
#if(__fvl_debug__==1)
      if(size(runtimes,1)<1) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument size for <runtimes>")
      end if
      do i=1,size(runtimes,1)
            if(.not. associated(runtimes(i)%getptr())) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument value for <runtimes>")
            end if
      end do
      if(size(models1,1)<1) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument size for <models1>")
      end if
      if(size(models1,2)<1) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument size for <models1>")
      end if
      if(size(models2,1)<1) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument size for <models2>")
      end if
      if(size(models2,2)<1) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument size for <models2>")
      end if
      if(present(models3)) then
            if(size(models3,1)<1) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument size for <models3>")
            end if
            if(size(models3,2)<1) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument size for <models3>")
            end if
      end if
      if(present(models4)) then
            if(size(models4,1)<1) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument size for <models4>")
            end if
            if(size(models4,2)<1) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument size for <models4>")
            end if
      end if
      if(present(models5)) then
            if(size(models5,1)<1) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument size for <models5>")
            end if
            if(size(models5,2)<1) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument size for <models5>")
            end if
      end if
      do i=1,size(models1,2)
            do j=1,size(models1,1)
                  if(.not. associated(models1(j,i)%getptr())) then
                        __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument value for <models1>")
                  end if
            end do
      end do
      do i=1,size(models2,2)
            do j=1,size(models2,1)
                  if(.not. associated(models2(j,i)%getptr())) then
                        __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument value for <models2>")
                  end if
            end do
      end do
      if(present(models3)) then
            do i=1,size(models3,2)
                  do j=1,size(models3,1)
                        if(.not. associated(models3(j,i)%getptr())) then
                              __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument value for <models3>")
                        end if
                  end do
            end do
      end if
      if(present(models4)) then
            do i=1,size(models4,2)
                  do j=1,size(models4,1)
                        if(.not. associated(models4(j,i)%getptr())) then
                              __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument value for <models4>")
                        end if
                  end do
            end do
      end if
      if(present(models5)) then
            do i=1,size(models5,2)
                  do j=1,size(models5,1)
                        if(.not. associated(models5(j,i)%getptr())) then
                              __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument value for <models5>")
                        end if
                  end do
            end do
      end if
      if(present(method)) then
            if(method/=fvl_rblock_linearsystem_none .and. method/=fvl_rblock_linearsystem_coupled .and.&
                  method/=fvl_rblock_linearsystem_segregated) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument value for <method>")
            end if
      end if
      if(present(predictor)) then
            if(predictor/=fvl_rblock_linearsystem_nopredictor .and. predictor/=fvl_rblock_linearsystem_taylorexpansion .and.&
                  predictor/=fvl_rblock_linearsystem_forwardeuler) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument value for <predictor>")
            end if
      end if
      if(present(numlevels)) then
            if(numlevels<2) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument value for <numlevels>")
            end if
      end if
      if(present(numsteps)) then
            if(numsteps<1) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: invalid argument value for <numsteps>")
            end if
      end if
      if(present(numlevels)) then
            if(numlevels==3 .and. (.not. present(models3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <models3>")
            end if
            if(numlevels==3 .and. (.not. present(tfun2))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <tfun2>")
            end if
            if(numlevels==3 .and. (.not. present(afun2))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <afun2>")
            end if
            if(numlevels==3 .and. (.not. present(bfun2))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <bfun2>")
            end if
            if(numlevels==3 .and. (.not. present(pfun2))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <pfun2>")
            end if
            if(numlevels==3 .and. (.not. present(ufun2))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <ufun2>")
            end if
            if(numlevels==4 .and. (.not. present(models4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <models4>")
            end if
            if(numlevels==4 .and. (.not. present(tfun3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <tfun3>")
            end if
            if(numlevels==4 .and. (.not. present(afun3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <afun3>")
            end if
            if(numlevels==4 .and. (.not. present(bfun3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <bfun3>")
            end if
            if(numlevels==4 .and. (.not. present(pfun3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <pfun3>")
            end if
            if(numlevels==4 .and. (.not. present(ufun3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <ufun3>")
            end if
            if(numlevels==5 .and. (.not. present(models5))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <models5>")
            end if
            if(numlevels==5 .and. (.not. present(tfun4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <tfun4>")
            end if
            if(numlevels==5 .and. (.not. present(afun4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <afun4>")
            end if
            if(numlevels==5 .and. (.not. present(bfun4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <bfun4>")
            end if
            if(numlevels==5 .and. (.not. present(pfun4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <pfun4>")
            end if
            if(numlevels==5 .and. (.not. present(ufun4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init2: missing argument <ufun4>")
            end if
      end if
#endif
      ! initialize data
      res%numruntimes=size(runtimes,1)
      call fvl_allocate(res%runtimes,res%numruntimes)
      do i=1,res%numruntimes
            runtime=>runtimes(i)%getptr()
            call res%runtimes(i)%setptr(runtime)
      end do
      res%tfun1=>tfun1
      res%afun1=>afun1
      res%bfun1=>bfun1
      res%pfun1=>pfun1
      res%ufun1=>ufun1
      if(present(tfun2)) then
            res%tfun2=>tfun2
      end if
      if(present(afun2)) then
            res%afun2=>afun2
      end if
      if(present(bfun2)) then
            res%bfun2=>bfun2
      end if
      if(present(pfun2)) then
            res%pfun2=>pfun2
      end if
      if(present(ufun2)) then
            res%ufun2=>ufun2
      end if
      if(present(tfun3)) then
            res%tfun3=>tfun3
      end if
      if(present(afun3)) then
            res%afun3=>afun3
      end if
      if(present(bfun3)) then
            res%bfun3=>bfun3
      end if
      if(present(pfun3)) then
            res%pfun3=>pfun3
      end if
      if(present(ufun3)) then
            res%ufun3=>ufun3
      end if
      if(present(tfun4)) then
            res%tfun4=>tfun4
      end if
      if(present(afun4)) then
            res%afun4=>afun4
      end if
      if(present(bfun4)) then
            res%bfun4=>bfun4
      end if
      if(present(pfun4)) then
            res%pfun4=>pfun4
      end if
      if(present(ufun4)) then
            res%ufun4=>ufun4
      end if
      res%linearsolver=>linearsolver
      if(present(method)) then
            res%method=method
      else
            res%method=fvl_rblock_linearsystem_method_default
      end if
      if(present(predictor)) then
            res%predictor=predictor
      else
            res%predictor=fvl_rblock_linearsystem_predictor_default
      end if
      if(present(numlevels)) then
            res%numlevels=numlevels
      else
            res%numlevels=fvl_rblock_linearsystem_numlevels_default
      end if
      if(present(numsteps)) then
            res%numsteps=numsteps
      else
            res%numsteps=fvl_rblock_linearsystem_numsteps_default
      end if
      ! initialise data
      res%nummodels=0
      do i=1,res%numsteps
            do j=1,res%numlevels
                  select case(j)
                        case(1)
                              models=>models1
                        case(2)
                              models=>models2
                        case(3)
                              models=>models3
                        case(4)
                              models=>models4
                        case(5)
                              models=>models5
                  end select
                  res%nummodels=res%nummodels+size(models,1)
            end do
      end do
      call fvl_allocate(res%models,res%nummodels,res%numlevels,res%numsteps)
      res%pointsize=0
      res%datasize=0
      aux=0
      do i=1,res%numsteps
            do j=1,res%numlevels
                  select case(j)
                        case(1)
                              models=>models1
                        case(2)
                              models=>models2
                        case(3)
                              models=>models3
                        case(4)
                              models=>models4
                        case(5)
                              models=>models5
                  end select
                  do k=1,size(models,1)
                        model=>models(k,i)%getptr()
                        call res%models(k,j,i)%setptr(model)
                        res%datasize=res%datasize+model%getfieldssize()
                        if(i==1 .and. j==1) then
                              res%pointsize=res%pointsize+model%getfieldssize()
                        else
                              aux=aux+model%getfieldssize()
                        end if
                  end do
#if(__fvl_debug__==1)
                  if(i/=1 .and. j/=1) then
                        if(aux/=res%pointsize) then
                              __fvl_log_error__("Internal error in call fvl_rblock_linearsystem_init2: inconsistent model sizes")
                        end if
                  end if
                  aux=0
#endif
            end do
      end do
      call res%make()
end function fvl_rblock_linearsystem_init2

function fvl_rblock_linearsystem_init3(runtime,models1,models2,models3,models4,models5,&
      tfun1,afun1,bfun1,pfun1,ufun1,tfun2,afun2,bfun2,pfun2,ufun2,tfun3,afun3,bfun3,pfun3,ufun3,tfun4,afun4,bfun4,pfun4,ufun4,&
      linearsolver,method,predictor,numlevels,numsteps) result(res)
      class(fvl_runtime),intent(in)::runtime
      class(fvl_model),dimension(1:),intent(in),target::models1
      class(fvl_model),dimension(1:),intent(in),target::models2
      class(fvl_model),dimension(1:),intent(in),optional,target::models3
      class(fvl_model),dimension(1:),intent(in),optional,target::models4
      class(fvl_model),dimension(1:),intent(in),optional,target::models5
      procedure(fvl_rblock_linearsystem_tfun)::tfun1
      procedure(fvl_rblock_linearsystem_afun)::afun1
      procedure(fvl_rblock_linearsystem_bfun)::bfun1
      procedure(fvl_rblock_linearsystem_pfun)::pfun1
      procedure(fvl_rblock_linearsystem_ufun)::ufun1
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun2
      procedure(fvl_rblock_linearsystem_afun),optional::afun2
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun2
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun2
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun2
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun3
      procedure(fvl_rblock_linearsystem_afun),optional::afun3
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun3
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun3
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun3
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun4
      procedure(fvl_rblock_linearsystem_afun),optional::afun4
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun4
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun4
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun4
      class(fvl_linearsolver),intent(in),target::linearsolver
      integer(kind=__fvl_integer_kind__),intent(in),optional::method
      integer(kind=__fvl_integer_kind__),intent(in),optional::predictor
      integer(kind=__fvl_integer_kind__),intent(in),optional::numlevels
      integer(kind=__fvl_integer_kind__),intent(in),optional::numsteps
      type(fvl_rblock_linearsystem)::res
      integer(kind=__fvl_integer_kind__)::i,j,aux
      class(fvl_model),dimension(:),pointer::models
#if(__fvl_debug__==1)
      if(present(method)) then
            if(method/=fvl_rblock_linearsystem_none .and. method/=fvl_rblock_linearsystem_coupled .and.&
                  method/=fvl_rblock_linearsystem_segregated) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: invalid argument value for <method>")
            end if
      end if
      if(present(predictor)) then
            if(predictor/=fvl_rblock_linearsystem_nopredictor .and. predictor/=fvl_rblock_linearsystem_taylorexpansion .and.&
                  predictor/=fvl_rblock_linearsystem_forwardeuler) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: invalid argument value for <predictor>")
            end if
      end if
      if(present(numlevels)) then
            if(numlevels<2) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: invalid argument value for <numlevels>")
            end if
      end if
      if(present(numsteps)) then
            if(numsteps<1) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: invalid argument value for <numsteps>")
            end if
      end if
      if(present(numlevels)) then
            if(numlevels==3 .and. (.not. present(models3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <models3>")
            end if
            if(numlevels==3 .and. (.not. present(tfun2))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <tfun2>")
            end if
            if(numlevels==3 .and. (.not. present(afun2))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <afun2>")
            end if
            if(numlevels==3 .and. (.not. present(bfun2))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <bfun2>")
            end if
            if(numlevels==3 .and. (.not. present(pfun2))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <pfun2>")
            end if
            if(numlevels==3 .and. (.not. present(ufun2))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <ufun2>")
            end if
            if(numlevels==4 .and. (.not. present(models4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <models4>")
            end if
            if(numlevels==4 .and. (.not. present(tfun3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <tfun3>")
            end if
            if(numlevels==4 .and. (.not. present(afun3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <afun3>")
            end if
            if(numlevels==4 .and. (.not. present(bfun3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <bfun3>")
            end if
            if(numlevels==4 .and. (.not. present(pfun3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <pfun3>")
            end if
            if(numlevels==4 .and. (.not. present(ufun3))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <ufun3>")
            end if
            if(numlevels==5 .and. (.not. present(models5))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <models5>")
            end if
            if(numlevels==5 .and. (.not. present(tfun4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <tfun4>")
            end if
            if(numlevels==5 .and. (.not. present(afun4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <afun4>")
            end if
            if(numlevels==5 .and. (.not. present(bfun4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <bfun4>")
            end if
            if(numlevels==5 .and. (.not. present(pfun4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <pfun4>")
            end if
            if(numlevels==5 .and. (.not. present(ufun4))) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init3: missing argument <ufun4>")
            end if
      end if
#endif
      ! initialize data
      res%numruntimes=1
      call fvl_allocate(res%runtimes,res%numruntimes)
      call res%runtimes(1)%setptr(runtime)
      res%tfun1=>tfun1
      res%afun1=>afun1
      res%bfun1=>bfun1
      res%pfun1=>pfun1
      res%ufun1=>ufun1
      if(present(tfun2)) then
            res%tfun2=>tfun2
      end if
      if(present(afun2)) then
            res%afun2=>afun2
      end if
      if(present(bfun2)) then
            res%bfun2=>bfun2
      end if
      if(present(pfun2)) then
            res%pfun2=>pfun2
      end if
      if(present(ufun2)) then
            res%ufun2=>ufun2
      end if
      if(present(tfun3)) then
            res%tfun3=>tfun3
      end if
      if(present(afun3)) then
            res%afun3=>afun3
      end if
      if(present(bfun3)) then
            res%bfun3=>bfun3
      end if
      if(present(pfun3)) then
            res%pfun3=>pfun3
      end if
      if(present(ufun3)) then
            res%ufun3=>ufun3
      end if
      if(present(tfun4)) then
            res%tfun4=>tfun4
      end if
      if(present(afun4)) then
            res%afun4=>afun4
      end if
      if(present(bfun4)) then
            res%bfun4=>bfun4
      end if
      if(present(pfun4)) then
            res%pfun4=>pfun4
      end if
      if(present(ufun4)) then
            res%ufun4=>ufun4
      end if
      res%linearsolver=>linearsolver
      if(present(method)) then
            res%method=method
      else
            res%method=fvl_rblock_linearsystem_method_default
      end if
      if(present(predictor)) then
            res%predictor=predictor
      else
            res%predictor=fvl_rblock_linearsystem_predictor_default
      end if
      if(present(numsteps)) then
            res%numsteps=numsteps
      else
            res%numsteps=fvl_rblock_linearsystem_numsteps_default
      end if
      if(present(numlevels)) then
            res%numlevels=numlevels
      else
            res%numlevels=fvl_rblock_linearsystem_numlevels_default
      end if
      ! initialise data
      res%nummodels=0
      do i=1,res%numsteps
            do j=1,res%numlevels
                  select case(j)
                        case(1)
                              models=>models1
                        case(2)
                              models=>models2
                        case(3)
                              models=>models3
                        case(4)
                              models=>models4
                        case(5)
                              models=>models5
                  end select
                  res%nummodels=res%nummodels+1
            end do
      end do
      call fvl_allocate(res%models,res%nummodels,res%numlevels,res%numsteps)
      res%pointsize=0
      res%datasize=0
      aux=0
      do i=1,res%numsteps
            do j=1,res%numlevels
                  select case(j)
                        case(1)
                              models=>models1
                        case(2)
                              models=>models2
                        case(3)
                              models=>models3
                        case(4)
                              models=>models4
                        case(5)
                              models=>models5
                  end select
                  call res%models(1,j,i)%setptr(models(i))
                  res%datasize=res%datasize+models(i)%getfieldssize()
                  if(i==1 .and. j==1) then
                        res%pointsize=res%pointsize+models(i)%getfieldssize()
                  else
                        aux=aux+models(i)%getfieldssize()
                  end if
#if(__fvl_debug__==1)
                  if(i/=1 .and. j/=1) then
                        if(aux/=res%pointsize) then
                              __fvl_log_error__("Internal error in call fvl_rblock_linearsystem_init3: inconsistent model sizes")
                        end if
                  end if
                  aux=0
#endif
            end do
      end do
      call res%make()
end function fvl_rblock_linearsystem_init3

function fvl_rblock_linearsystem_init4(runtimes,models1,models2,models3,models4,models5,&
      tfun1,afun1,bfun1,pfun1,ufun1,tfun2,afun2,bfun2,pfun2,ufun2,tfun3,afun3,bfun3,pfun3,ufun3,tfun4,afun4,bfun4,pfun4,ufun4,&
      linearsolver,dictfile,dictlabel) result(res)
      class(fvl_runtime),dimension(1:),contiguous,intent(in)::runtimes
      class(fvl_model),dimension(1:,1:),contiguous,intent(in),target::models1
      class(fvl_model),dimension(1:,1:),contiguous,intent(in),target::models2
      class(fvl_model),dimension(1:,1:),contiguous,intent(in),optional,target::models3
      class(fvl_model),dimension(1:,1:),contiguous,intent(in),optional,target::models4
      class(fvl_model),dimension(1:,1:),contiguous,intent(in),optional,target::models5
      procedure(fvl_rblock_linearsystem_tfun)::tfun1
      procedure(fvl_rblock_linearsystem_afun)::afun1
      procedure(fvl_rblock_linearsystem_bfun)::bfun1
      procedure(fvl_rblock_linearsystem_pfun)::pfun1
      procedure(fvl_rblock_linearsystem_ufun)::ufun1
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun2
      procedure(fvl_rblock_linearsystem_afun),optional::afun2
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun2
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun2
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun2
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun3
      procedure(fvl_rblock_linearsystem_afun),optional::afun3
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun3
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun3
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun3
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun4
      procedure(fvl_rblock_linearsystem_afun),optional::afun4
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun4
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun4
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun4
      class(fvl_linearsolver),intent(in),target::linearsolver
      class(fvl_dict_file),intent(in)::dictfile
      character(len=*),intent(in)::dictlabel
      type(fvl_rblock_linearsystem)::res
      integer(kind=__fvl_integer_kind__)::i,j,k,aux
      class(fvl_model),dimension(:,:),pointer::models
      ! initialize data
      res%numruntimes=size(runtimes,1)
      call fvl_allocate(res%runtimes,res%numruntimes)
      do i=1,res%numruntimes
            call res%runtimes(i)%setptr(runtimes(i))
      end do
      res%tfun1=>tfun1
      res%afun1=>afun1
      res%bfun1=>bfun1
      res%pfun1=>pfun1
      res%ufun1=>ufun1
      if(present(tfun2)) then
            res%tfun2=>tfun2
      end if
      if(present(afun2)) then
            res%afun2=>afun2
      end if
      if(present(bfun2)) then
            res%bfun2=>bfun2
      end if
      if(present(pfun2)) then
            res%pfun2=>pfun2
      end if
      if(present(ufun2)) then
            res%ufun2=>ufun2
      end if
      if(present(tfun3)) then
            res%tfun3=>tfun3
      end if
      if(present(afun3)) then
            res%afun3=>afun3
      end if
      if(present(bfun3)) then
            res%bfun3=>bfun3
      end if
      if(present(pfun3)) then
            res%pfun3=>pfun3
      end if
      if(present(ufun3)) then
            res%ufun3=>ufun3
      end if
      if(present(tfun4)) then
            res%tfun4=>tfun4
      end if
      if(present(afun4)) then
            res%afun4=>afun4
      end if
      if(present(bfun4)) then
            res%bfun4=>bfun4
      end if
      if(present(pfun4)) then
            res%pfun4=>pfun4
      end if
      if(present(ufun4)) then
            res%ufun4=>ufun4
      end if
      res%linearsolver=>linearsolver
      ! get parameters
      res%method=dictfile%getvalue(dictlabel,"method",fvl_rblock_linearsystem_method_default)
      res%predictor=dictfile%getvalue(dictlabel,"predictor",fvl_rblock_linearsystem_predictor_default)
      res%numlevels=dictfile%getvalue(dictlabel,"num_levels",fvl_rblock_linearsystem_numlevels_default)
      res%numsteps=dictfile%getvalue(dictlabel,"num_steps",fvl_rblock_linearsystem_numsteps_default)
      ! check parameters
      if(res%method/=fvl_rblock_linearsystem_none .and. res%method/=fvl_rblock_linearsystem_coupled .and.&
            res%method/=fvl_rblock_linearsystem_segregated) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: invalid argument value for <method>")
      end if
      if(res%predictor/=fvl_rblock_linearsystem_nopredictor .and. res%predictor/=fvl_rblock_linearsystem_taylorexpansion .and.&
            res%predictor/=fvl_rblock_linearsystem_forwardeuler) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: invalid argument value for <predictor>")
      end if
      if(res%numlevels<2) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: invalid argument value for <numlevels>")
      end if
      if(res%numsteps<1) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: invalid argument value for <numsteps>")
      end if
      if(res%numlevels==2 .and. (.not. present(models3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <models3>")
      end if
      if(res%numlevels==2 .and. (.not. present(tfun2))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <tfun2>")
      end if
      if(res%numlevels==2 .and. (.not. present(afun2))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <afun2>")
      end if
      if(res%numlevels==2 .and. (.not. present(bfun2))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <bfun2>")
      end if
      if(res%numlevels==2 .and. (.not. present(pfun2))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <pfun2>")
      end if
      if(res%numlevels==2 .and. (.not. present(ufun2))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <ufun2>")
      end if
      if(res%numlevels==3 .and. (.not. present(models4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <models4>")
      end if
      if(res%numlevels==3 .and. (.not. present(tfun3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <tfun3>")
      end if
      if(res%numlevels==3 .and. (.not. present(afun3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <afun3>")
      end if
      if(res%numlevels==3 .and. (.not. present(bfun3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <bfun3>")
      end if
      if(res%numlevels==3 .and. (.not. present(pfun3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <pfun3>")
      end if
      if(res%numlevels==3 .and. (.not. present(ufun3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <ufun3>")
      end if
      if(res%numlevels==4 .and. (.not. present(models5))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <models5>")
      end if
      if(res%numlevels==4 .and. (.not. present(tfun4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <tfun4>")
      end if
      if(res%numlevels==4 .and. (.not. present(afun4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <afun4>")
      end if
      if(res%numlevels==4 .and. (.not. present(bfun4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <bfun4>")
      end if
      if(res%numlevels==5 .and. (.not. present(pfun4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <pfun4>")
      end if
      if(res%numlevels==5 .and. (.not. present(ufun4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init4: missing argument <ufun4>")
      end if
      ! initialise data
      res%nummodels=0
      do i=1,res%numsteps
            do j=1,res%numlevels
                  select case(j)
                        case(1)
                              models=>models1
                        case(2)
                              models=>models2
                        case(3)
                              models=>models3
                        case(4)
                              models=>models4
                        case(5)
                              models=>models5
                  end select
                  res%nummodels=res%nummodels+size(models,1)
            end do
      end do
      call fvl_allocate(res%models,res%nummodels,res%numlevels,res%numsteps)
      res%pointsize=0
      res%datasize=0
      aux=0
      do i=1,res%numsteps
            do j=1,res%numlevels
                  select case(j)
                        case(1)
                              models=>models1
                        case(2)
                              models=>models2
                        case(3)
                              models=>models3
                        case(4)
                              models=>models4
                        case(5)
                              models=>models5
                  end select
                  do k=1,size(models,1)
                        call res%models(k,j,i)%setptr(models(k,i))
                        res%datasize=res%datasize+models(k,i)%getfieldssize()
                        if(i==1 .and. j==1) then
                              res%pointsize=res%pointsize+models(k,i)%getfieldssize()
                        else
                              aux=aux+models(k,i)%getfieldssize()
                        end if
                  end do
#if(__fvl_debug__==1)
                  if(i/=1 .and. j/=1) then
                        if(aux/=res%pointsize) then
                              __fvl_log_error__("Internal error in call fvl_rblock_linearsystem_init4: inconsistent model sizes")
                        end if
                  end if
                  aux=0
#endif
            end do
      end do
      call res%make()
end function fvl_rblock_linearsystem_init4

function fvl_rblock_linearsystem_init5(runtimes,models1,models2,models3,models4,models5,&
      tfun1,afun1,bfun1,pfun1,ufun1,tfun2,afun2,bfun2,pfun2,ufun2,tfun3,afun3,bfun3,pfun3,ufun3,tfun4,afun4,bfun4,pfun4,ufun4,&
      linearsolver,dictfile,dictlabel) result(res)
      class(fvl_ptr_runtime),dimension(1:),contiguous,intent(in)::runtimes
      class(fvl_ptr_model),dimension(1:,1:),contiguous,intent(in),target::models1
      class(fvl_ptr_model),dimension(1:,1:),contiguous,intent(in),target::models2
      class(fvl_ptr_model),dimension(1:,1:),contiguous,intent(in),optional,target::models3
      class(fvl_ptr_model),dimension(1:,1:),contiguous,intent(in),optional,target::models4
      class(fvl_ptr_model),dimension(1:,1:),contiguous,intent(in),optional,target::models5
      procedure(fvl_rblock_linearsystem_tfun)::tfun1
      procedure(fvl_rblock_linearsystem_afun)::afun1
      procedure(fvl_rblock_linearsystem_bfun)::bfun1
      procedure(fvl_rblock_linearsystem_pfun)::pfun1
      procedure(fvl_rblock_linearsystem_ufun)::ufun1
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun2
      procedure(fvl_rblock_linearsystem_afun),optional::afun2
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun2
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun2
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun2
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun3
      procedure(fvl_rblock_linearsystem_afun),optional::afun3
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun3
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun3
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun3
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun4
      procedure(fvl_rblock_linearsystem_afun),optional::afun4
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun4
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun4
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun4
      class(fvl_linearsolver),intent(in),target::linearsolver
      class(fvl_dict_file),intent(in)::dictfile
      character(len=*),intent(in)::dictlabel
      type(fvl_rblock_linearsystem)::res
      integer(kind=__fvl_integer_kind__)::i,j,k,aux
      class(fvl_runtime),pointer::runtime
      class(fvl_model),pointer::model
      class(fvl_ptr_model),dimension(:,:),pointer::models
#if(__fvl_debug__==1)
      if(size(runtimes,1)<1) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument size for <runtimes>")
      end if
      do i=1,size(runtimes,1)
            if(.not. associated(runtimes(i)%getptr())) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument value for <runtimes>")
            end if
      end do
      if(size(models1,1)<1) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument size for <models1>")
      end if
      if(size(models1,2)<1) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument size for <models1>")
      end if
      if(size(models2,1)<1) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument size for <models2>")
      end if
      if(size(models2,2)<1) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument size for <models2>")
      end if
      if(present(models3)) then
            if(size(models3,1)<1) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument size for <models3>")
            end if
            if(size(models3,2)<1) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument size for <models3>")
            end if
      end if
      if(present(models4)) then
            if(size(models4,1)<1) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument size for <models4>")
            end if
            if(size(models4,2)<1) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument size for <models4>")
            end if
      end if
      if(present(models5)) then
            if(size(models5,1)<1) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument size for <models5>")
            end if
            if(size(models5,2)<1) then
                  __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument size for <models5>")
            end if
      end if
      do i=1,size(models1,2)
            do j=1,size(models1,1)
                  if(.not. associated(models1(j,i)%getptr())) then
                        __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument value for <models1>")
                  end if
            end do
      end do
      do i=1,size(models2,2)
            do j=1,size(models2,1)
                  if(.not. associated(models2(j,i)%getptr())) then
                        __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument value for <models2>")
                  end if
            end do
      end do
      if(present(models3)) then
            do i=1,size(models3,2)
                  do j=1,size(models3,1)
                        if(.not. associated(models3(j,i)%getptr())) then
                              __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument value for <models3>")
                        end if
                  end do
            end do
      end if
      if(present(models4)) then
            do i=1,size(models4,2)
                  do j=1,size(models4,1)
                        if(.not. associated(models4(j,i)%getptr())) then
                              __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument value for <models4>")
                        end if
                  end do
            end do
      end if
      if(present(models5)) then
            do i=1,size(models5,2)
                  do j=1,size(models5,1)
                        if(.not. associated(models5(j,i)%getptr())) then
                              __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument value for <models5>")
                        end if
                  end do
            end do
      end if
#endif
      ! initialize data
      res%numruntimes=size(runtimes,1)
      call fvl_allocate(res%runtimes,res%numruntimes)
      do i=1,res%numruntimes
            runtime=>runtimes(i)%getptr()
            call res%runtimes(i)%setptr(runtime)
      end do
      res%tfun1=>tfun1
      res%afun1=>afun1
      res%bfun1=>bfun1
      res%pfun1=>pfun1
      res%ufun1=>ufun1
      if(present(tfun2)) then
            res%tfun2=>tfun2
      end if
      if(present(afun2)) then
            res%afun2=>afun2
      end if
      if(present(bfun2)) then
            res%bfun2=>bfun2
      end if
      if(present(pfun2)) then
            res%pfun2=>pfun2
      end if
      if(present(ufun2)) then
            res%ufun2=>ufun2
      end if
      if(present(tfun3)) then
            res%tfun3=>tfun3
      end if
      if(present(afun3)) then
            res%afun3=>afun3
      end if
      if(present(bfun3)) then
            res%bfun3=>bfun3
      end if
      if(present(pfun3)) then
            res%pfun3=>pfun3
      end if
      if(present(ufun3)) then
            res%ufun3=>ufun3
      end if
      if(present(tfun4)) then
            res%tfun4=>tfun4
      end if
      if(present(afun4)) then
            res%afun4=>afun4
      end if
      if(present(bfun4)) then
            res%bfun4=>bfun4
      end if
      if(present(pfun4)) then
            res%pfun4=>pfun4
      end if
      if(present(ufun4)) then
            res%ufun4=>ufun4
      end if
      res%linearsolver=>linearsolver
      ! get parameters
      res%method=dictfile%getvalue(dictlabel,"method",fvl_rblock_linearsystem_method_default)
      res%predictor=dictfile%getvalue(dictlabel,"predictor",fvl_rblock_linearsystem_predictor_default)
      res%numlevels=dictfile%getvalue(dictlabel,"num_levels",fvl_rblock_linearsystem_numlevels_default)
      res%numsteps=dictfile%getvalue(dictlabel,"num_steps",fvl_rblock_linearsystem_numsteps_default)
      ! check parameters
      if(res%method/=fvl_rblock_linearsystem_none .and. res%method/=fvl_rblock_linearsystem_coupled .and.&
            res%method/=fvl_rblock_linearsystem_segregated) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument value for <method>")
      end if
      if(res%predictor/=fvl_rblock_linearsystem_nopredictor .and. res%predictor/=fvl_rblock_linearsystem_taylorexpansion .and.&
            res%predictor/=fvl_rblock_linearsystem_forwardeuler) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument value for <predictor>")
      end if
      if(res%numlevels<2) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument value for <numlevels>")
      end if
      if(res%numsteps<1) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: invalid argument value for <numsteps>")
      end if
      if(res%numlevels==3 .and. (.not. present(models3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <models3>")
      end if
      if(res%numlevels==3 .and. (.not. present(tfun2))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <tfun2>")
      end if
      if(res%numlevels==3 .and. (.not. present(afun2))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <afun2>")
      end if
      if(res%numlevels==3 .and. (.not. present(bfun2))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <bfun2>")
      end if
      if(res%numlevels==3 .and. (.not. present(pfun2))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <pfun2>")
      end if
      if(res%numlevels==3 .and. (.not. present(ufun2))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <ufun2>")
      end if
      if(res%numlevels==4 .and. (.not. present(models4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <models4>")
      end if
      if(res%numlevels==4 .and. (.not. present(tfun3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <tfun3>")
      end if
      if(res%numlevels==4 .and. (.not. present(afun3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <afun3>")
      end if
      if(res%numlevels==4 .and. (.not. present(bfun3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <bfun3>")
      end if
      if(res%numlevels==4 .and. (.not. present(pfun3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <pfun3>")
      end if
      if(res%numlevels==4 .and. (.not. present(ufun3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <ufun3>")
      end if
      if(res%numlevels==5 .and. (.not. present(models5))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <models5>")
      end if
      if(res%numlevels==5 .and. (.not. present(tfun4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <tfun4>")
      end if
      if(res%numlevels==5 .and. (.not. present(afun4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <afun4>")
      end if
      if(res%numlevels==5 .and. (.not. present(bfun4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <bfun4>")
      end if
      if(res%numlevels==5 .and. (.not. present(pfun4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <pfun4>")
      end if
      if(res%numlevels==5 .and. (.not. present(ufun4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init5: missing argument <ufun4>")
      end if
      ! initialise data
      res%nummodels=0
      do i=1,res%numsteps
            do j=1,res%numlevels
                  select case(j)
                        case(1)
                              models=>models1
                        case(2)
                              models=>models2
                        case(3)
                              models=>models3
                        case(4)
                              models=>models4
                        case(5)
                              models=>models5
                  end select
                  res%nummodels=res%nummodels+size(models,1)
            end do
      end do
      call fvl_allocate(res%models,res%nummodels,res%numlevels,res%numsteps)
      res%pointsize=0
      res%datasize=0
      aux=0
      do i=1,res%numsteps
            do j=1,res%numlevels
                  select case(j)
                        case(1)
                              models=>models1
                        case(2)
                              models=>models2
                        case(3)
                              models=>models3
                        case(4)
                              models=>models4
                        case(5)
                              models=>models5
                  end select
                  do k=1,size(models,1)
                        model=>models(k,i)%getptr()
                        call res%models(k,j,i)%setptr(model)
                        res%datasize=res%datasize+model%getfieldssize()
                        if(i==1 .and. j==1) then
                              res%pointsize=res%pointsize+model%getfieldssize()
                        else
                              aux=aux+model%getfieldssize()
                        end if
                  end do
#if(__fvl_debug__==1)
                  if(i/=1 .and. j/=1) then
                        if(aux/=res%pointsize) then
                              __fvl_log_error__("Internal error in call fvl_rblock_linearsystem_init5: inconsistent model sizes")
                        end if
                  end if
                  aux=0
#endif
            end do
      end do
      call res%make()
end function fvl_rblock_linearsystem_init5

function fvl_rblock_linearsystem_init6(runtime,models1,models2,models3,models4,models5,&
      tfun1,afun1,bfun1,pfun1,ufun1,tfun2,afun2,bfun2,pfun2,ufun2,tfun3,afun3,bfun3,pfun3,ufun3,tfun4,afun4,bfun4,pfun4,ufun4,&
      linearsolver,dictfile,dictlabel) result(res)
      class(fvl_runtime),intent(in)::runtime
      class(fvl_model),dimension(1:),contiguous,intent(in),target::models1
      class(fvl_model),dimension(1:),contiguous,intent(in),target::models2
      class(fvl_model),dimension(1:),contiguous,intent(in),optional,target::models3
      class(fvl_model),dimension(1:),contiguous,intent(in),optional,target::models4
      class(fvl_model),dimension(1:),contiguous,intent(in),optional,target::models5
      procedure(fvl_rblock_linearsystem_tfun)::tfun1
      procedure(fvl_rblock_linearsystem_afun)::afun1
      procedure(fvl_rblock_linearsystem_bfun)::bfun1
      procedure(fvl_rblock_linearsystem_pfun)::pfun1
      procedure(fvl_rblock_linearsystem_ufun)::ufun1
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun2
      procedure(fvl_rblock_linearsystem_afun),optional::afun2
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun2
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun2
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun2
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun3
      procedure(fvl_rblock_linearsystem_afun),optional::afun3
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun3
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun3
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun3
      procedure(fvl_rblock_linearsystem_tfun),optional::tfun4
      procedure(fvl_rblock_linearsystem_afun),optional::afun4
      procedure(fvl_rblock_linearsystem_bfun),optional::bfun4
      procedure(fvl_rblock_linearsystem_pfun),optional::pfun4
      procedure(fvl_rblock_linearsystem_ufun),optional::ufun4
      class(fvl_linearsolver),intent(in),target::linearsolver
      class(fvl_dict_file),intent(in)::dictfile
      character(len=*),intent(in)::dictlabel
      type(fvl_rblock_linearsystem)::res
      integer(kind=__fvl_integer_kind__)::i,j,aux
      class(fvl_model),dimension(:),pointer::models
      ! initialize data
      res%numruntimes=1
      call fvl_allocate(res%runtimes,res%numruntimes)
      call res%runtimes(1)%setptr(runtime)
      res%tfun1=>tfun1
      res%afun1=>afun1
      res%bfun1=>bfun1
      res%pfun1=>pfun1
      res%ufun1=>ufun1
      if(present(tfun2)) then
            res%tfun2=>tfun2
      end if
      if(present(afun2)) then
            res%afun2=>afun2
      end if
      if(present(bfun2)) then
            res%bfun2=>bfun2
      end if
      if(present(pfun2)) then
            res%pfun2=>pfun2
      end if
      if(present(ufun2)) then
            res%ufun2=>ufun2
      end if
      if(present(tfun3)) then
            res%tfun3=>tfun3
      end if
      if(present(afun3)) then
            res%afun3=>afun3
      end if
      if(present(bfun3)) then
            res%bfun3=>bfun3
      end if
      if(present(pfun3)) then
            res%pfun3=>pfun3
      end if
      if(present(ufun3)) then
            res%ufun3=>ufun3
      end if
      if(present(tfun4)) then
            res%tfun4=>tfun4
      end if
      if(present(afun4)) then
            res%afun4=>afun4
      end if
      if(present(bfun4)) then
            res%bfun4=>bfun4
      end if
      if(present(pfun4)) then
            res%pfun4=>pfun4
      end if
      if(present(ufun4)) then
            res%ufun4=>ufun4
      end if
      res%linearsolver=>linearsolver
      ! get parameters
      res%method=dictfile%getvalue(dictlabel,"method",fvl_rblock_linearsystem_method_default)
      res%predictor=dictfile%getvalue(dictlabel,"predictor",fvl_rblock_linearsystem_predictor_default)
      res%numlevels=dictfile%getvalue(dictlabel,"num_levels",fvl_rblock_linearsystem_numlevels_default)
      res%numsteps=dictfile%getvalue(dictlabel,"num_steps",fvl_rblock_linearsystem_numsteps_default)
      ! check parameters
      if(res%method/=fvl_rblock_linearsystem_none .and. res%method/=fvl_rblock_linearsystem_coupled .and.&
            res%method/=fvl_rblock_linearsystem_segregated) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: invalid argument value for <method>")
      end if
      if(res%predictor/=fvl_rblock_linearsystem_nopredictor .and. res%predictor/=fvl_rblock_linearsystem_taylorexpansion .and.&
            res%predictor/=fvl_rblock_linearsystem_forwardeuler) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: invalid argument value for <predictor>")
      end if
      if(res%numlevels<2) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: invalid argument value for <numlevels>")
      end if
      if(res%numsteps<1) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: invalid argument value for <numsteps>")
      end if
      if(res%numlevels==3 .and. (.not. present(models3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <models3>")
      end if
      if(res%numlevels==3 .and. (.not. present(tfun2))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <tfun2>")
      end if
      if(res%numlevels==3 .and. (.not. present(afun2))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <afun2>")
      end if
      if(res%numlevels==3 .and. (.not. present(bfun2))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <bfun2>")
      end if
      if(res%numlevels==3 .and. (.not. present(pfun2))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <pfun2>")
      end if
      if(res%numlevels==3 .and. (.not. present(ufun2))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <ufun2>")
      end if
      if(res%numlevels==4 .and. (.not. present(models4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <models4>")
      end if
      if(res%numlevels==4 .and. (.not. present(tfun3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <tfun3>")
      end if
      if(res%numlevels==4 .and. (.not. present(afun3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <afun3>")
      end if
      if(res%numlevels==4 .and. (.not. present(bfun3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <bfun3>")
      end if
      if(res%numlevels==4 .and. (.not. present(pfun3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <pfun3>")
      end if
      if(res%numlevels==4 .and. (.not. present(ufun3))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <ufun3>")
      end if
      if(res%numlevels==5 .and. (.not. present(models5))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <models5>")
      end if
      if(res%numlevels==5 .and. (.not. present(tfun4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <tfun4>")
      end if
      if(res%numlevels==5 .and. (.not. present(afun4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <afun4>")
      end if
      if(res%numlevels==5 .and. (.not. present(bfun4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <bfun4>")
      end if
      if(res%numlevels==5 .and. (.not. present(pfun4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <pfun4>")
      end if
      if(res%numlevels==5 .and. (.not. present(ufun4))) then
            __fvl_log_error__("Input error in call fvl_rblock_linearsystem_init6: missing argument <ufun4>")
      end if
      ! initialise data
      res%nummodels=0
      do i=1,res%numsteps
            do j=1,res%numlevels
                  select case(j)
                        case(1)
                              models=>models1
                        case(2)
                              models=>models2
                        case(3)
                              models=>models3
                        case(4)
                              models=>models4
                        case(5)
                              models=>models5
                  end select
                  res%nummodels=res%nummodels+1
            end do
      end do
      call fvl_allocate(res%models,res%nummodels,res%numlevels,res%numsteps)
      res%pointsize=0
      res%datasize=0
      aux=0
      do i=1,res%numsteps
            do j=1,res%numlevels
                  select case(j)
                        case(1)
                              models=>models1
                        case(2)
                              models=>models2
                        case(3)
                              models=>models3
                        case(4)
                              models=>models4
                        case(5)
                              models=>models5
                  end select
                  call res%models(1,j,i)%setptr(models(i))
                  res%datasize=res%datasize+models(i)%getfieldssize()
                  if(i==1 .and. j==1) then
                        res%pointsize=res%pointsize+models(i)%getfieldssize()
                  else
                        aux=aux+models(i)%getfieldssize()
                  end if
#if(__fvl_debug__==1)
                  if(i/=1 .and. j/=1) then
                        if(aux/=res%pointsize) then
                              __fvl_log_error__("Internal error in call fvl_rblock_linearsystem_init6: inconsistent model sizes")
                        end if
                  end if
                  aux=0
#endif
            end do
      end do
      call res%make()
end function fvl_rblock_linearsystem_init6

subroutine fvl_rblock_linearsystem_delete(this)
      type(fvl_rblock_linearsystem),intent(inout)::this
      this%numruntimes=0
      this%nummodels=0
      this%numlevels=0
      this%numsteps=0
      this%pointsize=0
      this%datasize=0
      this%method=0
      this%predictor=0
      this%basesize=0
      this%currenttimestep=0
      this%tfun1=>null()
      this%afun1=>null()
      this%bfun1=>null()
      this%pfun1=>null()
      this%ufun1=>null()
      this%tfun2=>null()
      this%afun2=>null()
      this%bfun2=>null()
      this%pfun2=>null()
      this%ufun2=>null()
      this%tfun3=>null()
      this%afun3=>null()
      this%bfun3=>null()
      this%pfun3=>null()
      this%ufun3=>null()
      this%tfun4=>null()
      this%afun4=>null()
      this%bfun4=>null()
      this%pfun4=>null()
      this%ufun4=>null()
      call fvl_deallocate(this%runtimes)
      call fvl_deallocate(this%models)
      call fvl_deallocate(this%rhs)
      call fvl_deallocate(this%basecoefs)
      call fvl_deallocate(this%basecoefsrhs)
end subroutine fvl_rblock_linearsystem_delete

subroutine fvl_rblock_linearsystem_clean(this)
      type(fvl_rblock_linearsystem),intent(inout)::this
      this%numruntimes=0
      this%nummodels=0
      this%numlevels=0
      this%numsteps=0
      this%pointsize=0
      this%datasize=0
      this%method=0
      this%predictor=0
      this%basesize=0
      this%currenttimestep=0
      this%tfun1=>null()
      this%afun1=>null()
      this%bfun1=>null()
      this%pfun1=>null()
      this%ufun1=>null()
      this%tfun2=>null()
      this%afun2=>null()
      this%bfun2=>null()
      this%pfun2=>null()
      this%ufun2=>null()
      this%tfun3=>null()
      this%afun3=>null()
      this%bfun3=>null()
      this%pfun3=>null()
      this%ufun3=>null()
      this%tfun4=>null()
      this%afun4=>null()
      this%bfun4=>null()
      this%pfun4=>null()
      this%ufun4=>null()
      call fvl_clean(this%runtimes)
      call fvl_clean(this%models)
      call fvl_omp_clean(this%rhs)
      call fvl_omp_clean(this%basecoefs)
      call fvl_omp_clean(this%basecoefsrhs)
end subroutine fvl_rblock_linearsystem_clean

subroutine fvl_rblock_linearsystem_copy(this,that)
      type(fvl_rblock_linearsystem),intent(in)::this
      type(fvl_rblock_linearsystem),intent(out)::that
      that%numruntimes=this%numruntimes
      that%nummodels=this%nummodels
      that%numlevels=this%numlevels
      that%numsteps=this%numsteps
      that%pointsize=this%pointsize
      that%datasize=this%datasize
      that%method=this%method
      that%predictor=this%predictor
      that%basesize=this%basesize
      that%currenttimestep=this%currenttimestep
      that%tfun1=>this%tfun1
      that%afun1=>this%afun1
      that%bfun1=>this%bfun1
      that%pfun1=>this%pfun1
      that%ufun1=>this%ufun1
      that%tfun2=>this%tfun2
      that%afun2=>this%afun2
      that%bfun2=>this%bfun2
      that%pfun2=>this%pfun2
      that%ufun2=>this%ufun2
      that%tfun3=>this%tfun3
      that%afun3=>this%afun3
      that%bfun3=>this%bfun3
      that%pfun3=>this%pfun3
      that%ufun3=>this%ufun3
      that%tfun4=>this%tfun4
      that%afun4=>this%afun4
      that%bfun4=>this%bfun4
      that%pfun4=>this%pfun4
      that%ufun4=>this%ufun4
      call fvl_copy(this%runtimes,that%runtimes)
      call fvl_copy(this%models,that%models)
      call fvl_omp_copy(this%rhs,that%rhs)
      call fvl_omp_copy(this%basecoefs,that%basecoefs)
      call fvl_omp_copy(this%basecoefsrhs,that%basecoefsrhs)
end subroutine fvl_rblock_linearsystem_copy

function fvl_rblock_linearsystem_equal(this,that) result(res)
      type(fvl_rblock_linearsystem),intent(in)::this
      type(fvl_rblock_linearsystem),intent(in)::that
      logical(kind=__fvl_logical_kind__)::res
      if(this%numruntimes/=this%numruntimes .or.&
            this%nummodels/=this%nummodels .or.&
            that%numlevels/=this%numlevels .or.&
            that%numsteps/=this%numsteps .or.&
            that%pointsize/=this%pointsize .or.&
            that%datasize/=this%datasize .or.&
            that%method/=this%method .or.&
            that%predictor/=this%predictor .or.&
            that%basesize/=this%basesize .or.&
            .not. fvl_equal(that%currenttimestep,this%currenttimestep) .or.&
            .not. associated(that%tfun1,this%tfun1) .or.&
            .not. associated(that%afun1,this%afun1) .or.&
            .not. associated(that%bfun1,this%bfun1) .or.&
            .not. associated(that%pfun1,this%pfun1) .or.&
            .not. associated(that%ufun1,this%ufun1) .or.&
            .not. associated(that%tfun2,this%tfun2) .or.&
            .not. associated(that%afun2,this%afun2) .or.&
            .not. associated(that%bfun2,this%bfun2) .or.&
            .not. associated(that%pfun2,this%pfun2) .or.&
            .not. associated(that%ufun2,this%ufun2) .or.&
            .not. associated(that%tfun3,this%tfun3) .or.&
            .not. associated(that%afun3,this%afun3) .or.&
            .not. associated(that%bfun3,this%bfun3) .or.&
            .not. associated(that%pfun3,this%pfun3) .or.&
            .not. associated(that%ufun3,this%ufun3) .or.&
            .not. associated(that%tfun4,this%tfun4) .or.&
            .not. associated(that%afun4,this%afun4) .or.&
            .not. associated(that%bfun4,this%bfun4) .or.&
            .not. associated(that%pfun4,this%pfun4) .or.&
            .not. associated(that%ufun4,this%ufun4)) then
            res=.false.
            return
      else
            res=.true.
      end if
      if(.not. fvl_equal(this%runtimes,that%runtimes)) then
            res=.false.
            return
      end if
      if(.not. fvl_equal(this%models,that%models)) then
            res=.false.
            return
      end if
      if(.not. fvl_omp_equal(this%rhs,that%rhs)) then
            res=.false.
            return
      end if
      if(.not. fvl_omp_equal(this%basecoefs,that%basecoefs)) then
            res=.false.
            return
      end if
      if(.not. fvl_omp_equal(this%basecoefsrhs,that%basecoefsrhs)) then
            res=.false.
            return
      end if
end function fvl_rblock_linearsystem_equal

subroutine fvl_rblock_linearsystem_make(this)
      class(fvl_rblock_linearsystem),intent(inout)::this
      integer(kind=__fvl_integer_kind__)::k,l,r,m
      real(kind=__fvl_real_kind__),dimension(1:this%pointsize,1:this%numlevels,1:this%numsteps)::auxv1
      real(kind=__fvl_real_kind__),dimension(1:this%pointsize)::auxv2
      real(kind=__fvl_real_kind__),dimension(1:(this%numsteps+1)*this%numlevels,&
            1:(this%numsteps+1)*this%numlevels-this%numsteps)::auxm
      real(kind=__fvl_real_kind__),dimension(1:(this%numsteps+1)*this%numlevels,&
            1:(this%numsteps+1)*this%numlevels)::q
      procedure(fvl_rblock_linearsystem_afun),pointer::afun
      procedure(fvl_rblock_linearsystem_bfun),pointer::bfun
      procedure(fvl_rblock_linearsystem_ufun),pointer::ufun
      ! initialize data
      this%currenttimestep=0.0d0
      this%basesize=this%numsteps
      ! allocate memory
      call fvl_allocate(this%rhs,this%pointsize,this%basesize)
      call fvl_allocate(this%basecoefs,this%numlevels,this%numsteps,this%basesize)
      call fvl_allocate(this%basecoefsrhs,this%numlevels,this%basesize)
      ! compute physical equations initial conditions
      ! ! get solution
      ! call this%getsolution(auxv1)
      ! ! compute initial conditions
      ! do k=1,this%numlevels-1
      !       !select case(k)
      !             !case(1)
      !                   afun=>this%afun2
      !                   bfun=>this%bfun2
      !                   ufun=>this%ufun2
      !       !       case(2)
      !       !             afun=>this%afun2
      !       !             bfun=>this%bfun2
      !       !             ufun=>this%ufun2
      !       !       case(3)
      !       !             afun=>this%afun3
      !       !             bfun=>this%bfun3
      !       !             ufun=>this%ufun3
      !       !       case(4)
      !       !             afun=>this%afun4
      !       !             bfun=>this%bfun4
      !       !             ufun=>this%ufun4
      !       ! end select
      !
      !
      !
      !       call ufun(1)
      !       call afun(1,auxv1(:,k,1),auxv2)
      !       call bfun(1,auxv1(:,k+1,1))
      !       call fvl_omp_sub(auxv1(:,k+1,1),auxv2)
      !
      !       write(*,*) auxv1(:,k,1)
      !       write(*,*) " "
      !       write(*,*) auxv2
      !       write(*,*) " "
      !       write(*,*) auxv1(:,k+1,1)
      !       write(*,*) " "
      !       write(*,*) " "
      !       write(*,*) " "
      !       write(*,*) " "
      !       write(*,*) " "
      !       write(*,*) " "
      !       write(*,*) " "
      !       write(*,*) " "
      !
      ! end do
      ! ! set solution
      ! call this%setsolution(auxv1)
      ! compute structural equations coefficients
      m=0
      do while(m<this%numlevels*(this%numsteps+1)-this%basesize)
            m=m+1
            l=0
            do r=0,this%numsteps
                  do k=0,this%numlevels-1
                        l=l+1
                        if(m>k) then
                              auxm(l,m)=(fvl_real(fvl_fact(m-1))/fvl_real(fvl_fact(m-1-k)))*(fvl_real(r)**(m-1-k))
                        else
                              auxm(l,m)=0.0d0
                        end if
                  end do
            end do
      end do
      call fvl_qrfactq(auxm,q)
      ! extract coeffients
      do l=1,this%basesize
            m=1
            do k=1,this%numlevels
                  this%basecoefsrhs(k,l)=-q(this%numlevels*(this%numsteps+1)-this%basesize+l,m)
                  m=m+1
            end do
            do r=2,this%numsteps+1
                  do k=1,this%numlevels
                        this%basecoefs(k,this%numsteps-r+2,l)=q(this%numlevels*(this%numsteps+1)-this%basesize+l,m)
                        m=m+1
                  end do
            end do
      end do
end subroutine fvl_rblock_linearsystem_make

function fvl_rblock_linearsystem_gettfun1(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_tfun),pointer::res
      res=>this%tfun1
end function fvl_rblock_linearsystem_gettfun1

function fvl_rblock_linearsystem_getafun1(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_afun),pointer::res
      res=>this%afun1
end function fvl_rblock_linearsystem_getafun1

function fvl_rblock_linearsystem_getbfun1(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_bfun),pointer::res
      res=>this%bfun1
end function fvl_rblock_linearsystem_getbfun1

function fvl_rblock_linearsystem_getpfun1(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_pfun),pointer::res
      res=>this%pfun1
end function fvl_rblock_linearsystem_getpfun1

function fvl_rblock_linearsystem_getufun1(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_ufun),pointer::res
      res=>this%ufun1
end function fvl_rblock_linearsystem_getufun1

function fvl_rblock_linearsystem_gettfun2(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_tfun),pointer::res
      res=>this%tfun2
end function fvl_rblock_linearsystem_gettfun2

function fvl_rblock_linearsystem_getafun2(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_afun),pointer::res
      res=>this%afun2
end function fvl_rblock_linearsystem_getafun2

function fvl_rblock_linearsystem_getbfun2(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_bfun),pointer::res
      res=>this%bfun2
end function fvl_rblock_linearsystem_getbfun2

function fvl_rblock_linearsystem_getpfun2(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_pfun),pointer::res
      res=>this%pfun2
end function fvl_rblock_linearsystem_getpfun2

function fvl_rblock_linearsystem_getufun2(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_ufun),pointer::res
      res=>this%ufun2
end function fvl_rblock_linearsystem_getufun2

function fvl_rblock_linearsystem_gettfun3(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_tfun),pointer::res
      res=>this%tfun3
end function fvl_rblock_linearsystem_gettfun3

function fvl_rblock_linearsystem_getafun3(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_afun),pointer::res
      res=>this%afun3
end function fvl_rblock_linearsystem_getafun3

function fvl_rblock_linearsystem_getbfun3(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_bfun),pointer::res
      res=>this%bfun3
end function fvl_rblock_linearsystem_getbfun3

function fvl_rblock_linearsystem_getpfun3(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_pfun),pointer::res
      res=>this%pfun3
end function fvl_rblock_linearsystem_getpfun3

function fvl_rblock_linearsystem_getufun3(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_ufun),pointer::res
      res=>this%ufun3
end function fvl_rblock_linearsystem_getufun3

function fvl_rblock_linearsystem_gettfun4(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_tfun),pointer::res
      res=>this%tfun4
end function fvl_rblock_linearsystem_gettfun4

function fvl_rblock_linearsystem_getafun4(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_afun),pointer::res
      res=>this%afun4
end function fvl_rblock_linearsystem_getafun4

function fvl_rblock_linearsystem_getbfun4(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_bfun),pointer::res
      res=>this%bfun4
end function fvl_rblock_linearsystem_getbfun4

function fvl_rblock_linearsystem_getpfun4(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_pfun),pointer::res
      res=>this%pfun4
end function fvl_rblock_linearsystem_getpfun4

function fvl_rblock_linearsystem_getufun4(this) result(res)
      class(fvl_rblock_linearsystem),intent(in)::this
      procedure(fvl_rblock_linearsystem_ufun),pointer::res
      res=>this%ufun4
end function fvl_rblock_linearsystem_getufun4

subroutine fvl_rblock_linearsystem_evallhs(this,x,res)
      class(fvl_rblock_linearsystem),intent(in)::this
      real(kind=__fvl_real_kind__),dimension(1:this%datasize),intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:this%datasize),intent(out)::res
      integer(kind=__fvl_integer_kind__)::i,j,k,l,m
      real(kind=__fvl_real_kind__),dimension(1:this%pointsize)::auxv
      procedure(fvl_rblock_linearsystem_tfun),pointer::tfun
      procedure(fvl_rblock_linearsystem_afun),pointer::afun
      ! physical equations
      m=0
      k=0
      do i=1,this%numsteps
            do j=1,this%numlevels-1
                  select case(j)
                        case(1)
                              tfun=>this%tfun1
                              afun=>this%afun1
                        case(2)
                              tfun=>this%tfun2
                              afun=>this%afun2
                        case(3)
                              tfun=>this%tfun3
                              afun=>this%afun3
                        case(4)
                              tfun=>this%tfun4
                              afun=>this%afun4
                  end select
                  l=k+this%pointsize
                  call tfun(i,x(l+1:l+this%pointsize),res(m+1:m+this%pointsize))
                  call afun(i,x(k+1:k+this%pointsize),auxv)
                  call fvl_omp_add(res(m+1:m+this%pointsize),auxv)
                  k=k+this%pointsize
                  m=m+this%pointsize
            end do
            k=k+this%pointsize
      end do
      ! structural equations
      do i=1,this%basesize
            call fvl_omp_zeros(res(m+1:m+this%pointsize))
            l=0
            do j=1,this%numsteps
                  do k=1,this%numlevels
                        call fvl_omp_addmul(res(m+1:m+this%pointsize),&
                              x(l+1:l+this%pointsize),(this%currenttimestep**fvl_real(k-1))*this%basecoefs(k,j,i))
                        l=l+this%pointsize
                  end do
            end do
            m=m+this%pointsize
      end do
end subroutine fvl_rblock_linearsystem_evallhs

subroutine fvl_rblock_linearsystem_evalrhs(this,res)
      class(fvl_rblock_linearsystem),intent(in)::this
      real(kind=__fvl_real_kind__),dimension(1:this%datasize),intent(out)::res
      integer(kind=__fvl_integer_kind__)::i,j,m
      procedure(fvl_rblock_linearsystem_bfun),pointer::bfun
      ! physical equations
      m=0
      do i=1,this%numsteps
            do j=1,this%numlevels-1
                  select case(j)
                        case(1)
                              bfun=>this%bfun1
                        case(2)
                              bfun=>this%bfun2
                        case(3)
                              bfun=>this%bfun3
                        case(4)
                              bfun=>this%bfun4
                  end select
                  call bfun(i,res(m+1:m+this%pointsize))
                  m=m+this%pointsize
            end do
      end do
      ! structural equations
      do i=1,this%basesize
            call fvl_omp_assign(this%rhs(:,i),res(m+1:m+this%pointsize))
            m=m+this%pointsize
      end do
end subroutine fvl_rblock_linearsystem_evalrhs

subroutine fvl_rblock_linearsystem_precond(this,x,res)
      class(fvl_rblock_linearsystem),intent(in)::this
      real(kind=__fvl_real_kind__),dimension(1:this%datasize),intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:this%datasize),intent(out)::res
      integer(kind=__fvl_integer_kind__)::i,j,m
      procedure(fvl_rblock_linearsystem_pfun),pointer::pfun
      ! physical equations
      m=0
      do i=1,this%numsteps
            do j=1,this%numlevels-1
                  select case(j)
                        case(1)
                              pfun=>this%pfun1
                        case(2)
                              pfun=>this%pfun2
                        case(3)
                              pfun=>this%pfun3
                        case(4)
                              pfun=>this%pfun4
                  end select
                  call pfun(i,x(m+1:m+this%pointsize),res(m+1:m+this%pointsize))
                  m=m+this%pointsize
            end do
      end do
      ! structural equations
      call fvl_omp_assign(x(m+1:),res(m+1:))
end subroutine fvl_rblock_linearsystem_precond

subroutine fvl_rblock_linearsystem_update(this)
      class(fvl_rblock_linearsystem),intent(inout)::this
      integer(kind=__fvl_integer_kind__)::i,j
      real(kind=__fvl_real_kind__),dimension(1:this%pointsize,1:this%numlevels,1:this%numsteps)::auxv1
      real(kind=__fvl_real_kind__),dimension(1:this%pointsize)::auxv2,auxv3
      class(fvl_runtime),pointer::runtime
      procedure(fvl_rblock_linearsystem_afun),pointer::afun
      procedure(fvl_rblock_linearsystem_bfun),pointer::bfun
      procedure(fvl_rblock_linearsystem_ufun),pointer::ufun
      ! get timestep
      runtime=>this%runtimes(1)%getptr()
      this%currenttimestep=runtime%getcurrenttimestep()/fvl_real(this%numsteps)
      ! get solution
      call this%getsolution(auxv1)
      ! compute rhs
      ! structural equations
      do i=1,this%basesize
            call fvl_omp_mul(auxv1(:,1,1),this%basecoefsrhs(1,i),this%rhs(:,i))
            do j=2,this%numlevels
                  call fvl_omp_addmul(this%rhs(:,i),auxv1(:,j,1),(this%currenttimestep**fvl_real(j-1))*this%basecoefsrhs(j,i))
            end do
      end do
      ! extrapolate solution
      if(this%predictor==fvl_rblock_linearsystem_nopredictor) then
            return
      end if
      ! first step
      ! first level
      if(this%predictor==fvl_rblock_linearsystem_taylorexpansion) then
            call fvl_omp_assign(auxv1(:,1,1),auxv1(:,1,this%numsteps))
            do j=2,this%numlevels
                  call fvl_omp_addmul(auxv1(:,1,this%numsteps),auxv1(:,j,1),&
                        (this%currenttimestep**fvl_real(j-1))/fvl_fact(fvl_real(j-1)))
            end do
      else if(this%predictor==fvl_rblock_linearsystem_forwardeuler) then
            afun=>this%afun1
            bfun=>this%bfun1
            ufun=>this%ufun1
            call fvl_omp_assign(auxv1(:,1,1),auxv3)
            call ufun(1)
            call afun(1,auxv1(:,1,1),auxv2)
            call bfun(1,auxv1(:,1,this%numsteps))
            call fvl_omp_sub(auxv1(:,1,this%numsteps),auxv2)
            call fvl_omp_muladd(auxv1(:,1,this%numsteps),this%currenttimestep,auxv3)
      end if
      ! other levels
      do j=1,this%numlevels-1
            select case(j)
                  case(1)
                        afun=>this%afun1
                        bfun=>this%bfun1
                        ufun=>this%ufun1
                  case(2)
                        afun=>this%afun2
                        bfun=>this%bfun2
                        ufun=>this%ufun2
                  case(3)
                        afun=>this%afun3
                        bfun=>this%bfun3
                        ufun=>this%ufun3
                  case(4)
                        afun=>this%afun4
                        bfun=>this%bfun4
                        ufun=>this%ufun4
            end select
            call ufun(1)
            call afun(1,auxv1(:,j,this%numsteps),auxv2)
            call bfun(1,auxv1(:,j+1,this%numsteps))
            call fvl_omp_sub(auxv1(:,j+1,this%numsteps),auxv2)
      end do
      ! other steps
      do i=1,this%numsteps-1
            ! first level
            if(this%predictor==fvl_rblock_linearsystem_taylorexpansion) then
                  call fvl_omp_assign(auxv1(:,1,this%numsteps-i+1),auxv1(:,1,this%numsteps-i))
                  do j=2,this%numlevels
                        call fvl_omp_addmul(auxv1(:,1,this%numsteps-i),auxv1(:,j,this%numsteps-i+1),&
                              (this%currenttimestep**fvl_real(j-1))/fvl_fact(fvl_real(j-1)))
                  end do
            else if(this%predictor==fvl_rblock_linearsystem_forwardeuler) then
                  afun=>this%afun1
                  bfun=>this%bfun1
                  ufun=>this%ufun1
                  call fvl_omp_assign(auxv1(:,1,this%numsteps-i+1),auxv3)
                  call ufun(1)
                  call afun(1,auxv1(:,1,this%numsteps-i+1),auxv2)
                  call bfun(1,auxv1(:,1,this%numsteps-i))
                  call fvl_omp_sub(auxv1(:,1,this%numsteps-i),auxv2)
                  call fvl_omp_muladd(auxv1(:,1,this%numsteps-i),this%currenttimestep,auxv3)
            end if
            ! other levels
            do j=1,this%numlevels-1
                  select case(j)
                        case(1)
                              afun=>this%afun1
                              bfun=>this%bfun1
                              ufun=>this%ufun1
                        case(2)
                              afun=>this%afun2
                              bfun=>this%bfun2
                              ufun=>this%ufun2
                        case(3)
                              afun=>this%afun3
                              bfun=>this%bfun3
                              ufun=>this%ufun3
                        case(4)
                              afun=>this%afun4
                              bfun=>this%bfun4
                              ufun=>this%ufun4
                  end select
                  call ufun(i)
                  call afun(i,auxv1(:,j,this%numsteps-i),auxv2)
                  call bfun(i,auxv1(:,j+1,this%numsteps-i))
                  call fvl_omp_sub(auxv1(:,j+1,this%numsteps-i),auxv2)
            end do
      end do
      ! set solution
      call this%setsolution(auxv1)
end subroutine fvl_rblock_linearsystem_update

subroutine fvl_rblock_linearsystem_solve(this)
      class(fvl_rblock_linearsystem),intent(inout)::this
      ! compute solution
      if(this%method==fvl_rblock_linearsystem_coupled) then
            call this%solvecoupled()
      else if(this%method==fvl_rblock_linearsystem_segregated) then
            call this%solvesegregated()
      end if
end subroutine fvl_rblock_linearsystem_solve

subroutine fvl_rblock_linearsystem_solvecoupled(this)
      class(fvl_rblock_linearsystem),intent(inout)::this
      real(kind=__fvl_real_kind__),dimension(1:this%datasize)::auxv
      integer(kind=__fvl_integer_kind__)::i,j
      procedure(fvl_rblock_linearsystem_ufun),pointer::ufun
      ! update linear system
      do i=1,this%numsteps
            do j=1,this%numlevels-1
                  select case(j)
                        case(1)
                              call this%ufun1(i)
                        case(2)
                              call this%ufun2(i)
                        case(3)
                              call this%ufun3(i)
                        case(4)
                              call this%ufun4(i)
                  end select
            end do
      end do
      ! get solution
      call this%getsolution(auxv)
      ! compute solution
      call this%linearsolver%solve(this,auxv)
      ! set solution
      call this%setsolution(auxv)
end subroutine fvl_rblock_linearsystem_solvecoupled

subroutine fvl_rblock_linearsystem_solvesegregated(this)
      class(fvl_rblock_linearsystem),intent(inout)::this
      real(kind=__fvl_real_kind__),dimension(1:this%pointsize,1:this%numlevels,1:this%numsteps)::auxv1
      real(kind=__fvl_real_kind__),dimension(1:this%pointsize)::auxv2
      real(kind=__fvl_real_kind__),dimension(1:this%basesize)::auxv3,auxv4
      real(kind=__fvl_real_kind__),dimension(1:this%basesize,1:this%basesize)::auxm
      integer(kind=__fvl_integer_kind__)::i,j,k,l
      procedure(fvl_rblock_linearsystem_afun),pointer::afun
      procedure(fvl_rblock_linearsystem_bfun),pointer::bfun
      procedure(fvl_rblock_linearsystem_ufun),pointer::ufun
      ! update linear system
      do i=1,this%numsteps
            do j=1,this%numlevels-1
                  select case(j)
                        case(1)
                              call this%ufun1(i)
                        case(2)
                              call this%ufun2(i)
                        case(3)
                              call this%ufun3(i)
                        case(4)
                              call this%ufun4(i)
                  end select
            end do
      end do
      ! get solution
      call this%getsolution(auxv1)
      ! structural equations
      ! coeff matrix
      do i=1,this%basesize
            do j=1,this%numsteps
                  auxm(i,j)=this%basecoefs(1,j,i)
            end do
      end do
      ! right-hand side
      do l=1,this%pointsize
            do i=1,this%basesize
                  auxv3(i)=this%rhs(l,i)
                  do j=1,this%numsteps
                        do k=2,this%numlevels
                              auxv3(i)=auxv3(i)-auxv1(l,k,j)*(this%currenttimestep**fvl_real(k-1))*this%basecoefs(k,j,i)
                        end do
                  end do
            end do
            call fvl_solve(auxm,auxv3,auxv4)
            do i=1,this%basesize
                  auxv1(l,1,i)=auxv4(i)
            end do
      end do
      ! physical equations
      do i=1,this%numsteps
            do j=1,this%numlevels-1
                  select case(j)
                        case(1)
                              afun=>this%afun1
                              bfun=>this%bfun1
                              ufun=>this%ufun1
                        case(2)
                              afun=>this%afun2
                              bfun=>this%bfun2
                              ufun=>this%ufun2
                        case(3)
                              afun=>this%afun3
                              bfun=>this%bfun3
                              ufun=>this%ufun3
                        case(4)
                              afun=>this%afun4
                              bfun=>this%bfun4
                              ufun=>this%ufun4
                  end select
                  call ufun(i)
                  call afun(i,auxv1(:,j,i),auxv2)
                  call bfun(i,auxv1(:,j+1,i))
                  call fvl_omp_sub(auxv1(:,j+1,i),auxv2)
            end do
      end do
      ! set solution
      call this%setsolution(auxv1)
end subroutine fvl_rblock_linearsystem_solvesegregated

end module fvl_rblock_linearsystem_mod
! end of file
