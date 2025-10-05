! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author:  Ricardo Costa           |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0a                    |
!| |_|       \_/   |_|_|_| |  Release: June 2016               |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: 3D non-isothermal incompressible fluid flow post-processor
! Modification: October, 2023

#include "macros.f"

program fvl_post_thermoincompflow3d_main

use fvl_lib3d

implicit none

! ============================================================================
! DECLARE VARIABLES
! ============================================================================

! ----------------------------------------------------------------------------
! variables
! ----------------------------------------------------------------------------

integer(kind=__fvl_integer_kind__),parameter::maxnummeshes=10
character(len=__fvl_character_len__)::postdictpath,t_suffix,p_suffix,u_suffix,t_label,p_label,u_label,&
      intformat,realformat,primal_doflabel,diamond_doflabel,&
      meshesfilepaths(maxnummeshes),meshesfileforms(maxnummeshes),solutionsfilesdir(maxnummeshes),solutionsfilesform(maxnummeshes),&
      paramsfilesdir(maxnummeshes),paramsfilesform(maxnummeshes),errorsfilesfir(maxnummeshes),errorsfilesform(maxnummeshes),&
      refmeshesfilepaths(maxnummeshes),refmeshesfileforms(maxnummeshes),refsolutionsfilesdir(maxnummeshes),refsolutionsfilesform(maxnummeshes),&
      refparamsfilesdir(maxnummeshes),refparamsfilesform(maxnummeshes),tablesfilesdir
integer(kind=__fvl_integer_kind__)::i,j,nummeshes,convergence,p_fixconstant,tablesprecision,unit,element,&
      primal_numdofs(maxnummeshes),diamond_numdofs(maxnummeshes),&
      primal_numcells(maxnummeshes),diamond_numcells(maxnummeshes),steps(maxnummeshes),&
      linearsystemnumiters(maxnummeshes),fixpointnumiters(maxnummeshes)
real(kind=__fvl_real_kind__)::time,aux,cellvolume,primalmeshvolumes(maxnummeshes),diamondmeshvolumes(maxnummeshes),&
      times(maxnummeshes),fixpointresid(maxnummeshes),linearsystemresid(maxnummeshes),&
      meshestotaltime(maxnummeshes),modelstotaltime(maxnummeshes),schemestotaltime(maxnummeshes),solutiontotaltime(maxnummeshes),&
      postproctotaltime(maxnummeshes),executiontotaltime(maxnummeshes),&
      refmeshestotaltime(maxnummeshes),refmodelstotaltime(maxnummeshes),refschemestotaltime(maxnummeshes),refsolutiontotaltime(maxnummeshes),&
      refpostproctotaltime(maxnummeshes),refexecutiontotaltime(maxnummeshes),&
      meshestotaltimeratio(maxnummeshes),modelstotaltimeratio(maxnummeshes),schemestotaltimeratio(maxnummeshes),solutiontotaltimeratio(maxnummeshes),&
      postproctotaltimeratio(maxnummeshes),executiontotaltimeratio(maxnummeshes),&
      totalexecutionmemory(maxnummeshes),reftotalexecutionmemory(maxnummeshes),totalexecutionmemoryratio(maxnummeshes),&
      t_error1(maxnummeshes),t_error2(maxnummeshes),t_errorinf(maxnummeshes),&
      t_order1(maxnummeshes),t_order2(maxnummeshes),t_orderinf(maxnummeshes),&
      p_error1(maxnummeshes),p_error2(maxnummeshes),p_errorinf(maxnummeshes),&
      p_order1(maxnummeshes),p_order2(maxnummeshes),p_orderinf(maxnummeshes),&
      u_error1(maxnummeshes),u_error2(maxnummeshes),u_errorinf(maxnummeshes),&
      u_order1(maxnummeshes),u_order2(maxnummeshes),u_orderinf(maxnummeshes),&
      ux_error1(maxnummeshes),ux_error2(maxnummeshes),ux_errorinf(maxnummeshes),&
      ux_order1(maxnummeshes),ux_order2(maxnummeshes),ux_orderinf(maxnummeshes),&
      uy_error1(maxnummeshes),uy_error2(maxnummeshes),uy_errorinf(maxnummeshes),&
      uy_order1(maxnummeshes),uy_order2(maxnummeshes),uy_orderinf(maxnummeshes),&
      uz_error1(maxnummeshes),uz_error2(maxnummeshes),uz_errorinf(maxnummeshes),&
      uz_order1(maxnummeshes),uz_order2(maxnummeshes),uz_orderinf(maxnummeshes),&
      diamond_u_error1(maxnummeshes),diamond_u_error2(maxnummeshes),diamond_u_errorinf(maxnummeshes),&
      diamond_u_order1(maxnummeshes),diamond_u_order2(maxnummeshes),diamond_u_orderinf(maxnummeshes),&
      diamond_ux_error1(maxnummeshes),diamond_ux_error2(maxnummeshes),diamond_ux_errorinf(maxnummeshes),&
      diamond_ux_order1(maxnummeshes),diamond_ux_order2(maxnummeshes),diamond_ux_orderinf(maxnummeshes),&
      diamond_uy_error1(maxnummeshes),diamond_uy_error2(maxnummeshes),diamond_uy_errorinf(maxnummeshes),&
      diamond_uy_order1(maxnummeshes),diamond_uy_order2(maxnummeshes),diamond_uy_orderinf(maxnummeshes),&
      diamond_uz_error1(maxnummeshes),diamond_uz_error2(maxnummeshes),diamond_uz_errorinf(maxnummeshes),&
      diamond_uz_order1(maxnummeshes),diamond_uz_order2(maxnummeshes),diamond_uz_orderinf(maxnummeshes)
integer(kind=__fvl_integer_kind__),allocatable,dimension(:)::t_map,p_map,u_map,diamond_u_map,linearsystemiterhistogram
real(kind=__fvl_real_kind__),allocatable,dimension(:)::fixpointresidhistogram,linearsystemresidhistogram
real(kind=__fvl_real_kind__),allocatable,dimension(:)::t_refvals,t_exactvals,t_approxvals,t_errorvals,p_refvals,p_exactvals,p_approxvals,p_errorvals,&
      u_refvals2,diamond_u_refvals2,auxv
real(kind=__fvl_real_kind__),allocatable,dimension(:,:)::u_refvals1,u_exactvals,u_approxvals,u_errorvals,diamond_u_refvals1,diamond_u_exactvals,&
      diamond_u_approxvals,diamond_u_errorvals
logical::exist1,exist2
type(fvl_mesh3d)::primalmesh,diamondmesh,primalrefmesh,diamondrefmesh
type(fvl_file)::file
type(fvl_keys_file)::paramfile
type(fvl_data_file)::solutionsfile
type(fvl_dict_file)::controldictfile
type(fvl_patch3d)::primal_patch,diamond_patch
type(fvl_meshlocator3d)::primalmeshlocator,diamondmeshlocator
type(fvl_reconstmap3d)::t_reconstmap,p_reconstmap,ux_reconstmap,uy_reconstmap,uz_reconstmap,diamond_ux_reconstmap,diamond_uy_reconstmap,diamond_uz_reconstmap

! user-defined variables
#include "variables.f"

! reports variables
#include "fvl_writereports_hdr.f"

! dictionaries variables
#include "fvl_readdictionaries_hdr.f"

! mesh variables
#include "fvl_creatediamondmesh2d_hdr.f"

! fields variables
#include "fvl_writefields_hdr.f"

! model variables
#include "fvl_createmodels_hdr.f"

! field allocators
#include "fvl_setallocators_hdr.f"

! ============================================================================
! INITIALIZE LIBRARY
! ============================================================================

! initialize library
call fvl_init()

! check options
if(.not. fvl_checkopts("-h","--help","-m","--meshes","-o","--models","-s","--schemes",&
      "-l","--solution ","-c","--control","-p","--parallel")) then
      call fvl_logerror("Invalid option.")
end if

! check help option
if(fvl_findopt("-h") .or. fvl_findopt("--help")) then
      call fvl_loginfo("About: 3D non-isothermal non-Newtonian incompressible fluid flow post-processor.")
      call fvl_loginfo("Usage: fvl-postcell-thermoincompflow3d [options]")
      call fvl_loginfo("Options:")
      call fvl_loginfo("      -h|--help               Displays help message.")
      call fvl_loginfo("      -m|--meshes <file>      Reads meshes parameters from dictionary file <file> (default is setup/meshes.fvd).")
      call fvl_loginfo("      -o|--models <file>      Reads models parameters from dictionary file <file> (default is setup/models.fvd).")
      call fvl_loginfo("      -s|--schemes <file>     Reads schemes parameters from dictionary file <file> (default is setup/schemes.fvd).")
      call fvl_loginfo("      -l|--solution <file>    Reads solution parameters from dictionary file <file> (default is setup/solution.fvd).")
      call fvl_loginfo("      -c|--control <file>     Reads control parameters from dictionary file <file> (default is setup/control.fvd).")
      call fvl_loginfo("      -t|--post <file>        Reads post parameters from dictionary file <file> (default is setup/post.fvd).")
      call fvl_loginfo("      -p|--parallel           Runs in parallel with <OMP_NUM_THREADS> processes (default is serial).")
      stop
end if

! start total execution timer
call fvl_startwtimer(10,"execution_total")

! ============================================================================
! READ DICTIONARIES
! ============================================================================

! state message
call fvl_loginfo("Reading dictionaries...")

! templated code
#include "fvl_readdictionaries_src.f"

! read post parameters
convergence             = controldictfile%getvalue("post","convergence",1)
t_label                 = controldictfile%getvalue("post","t_label","")
p_label                 = controldictfile%getvalue("post","p_label","")
u_label                 = controldictfile%getvalue("post","u_label","")
p_fixconstant           = controldictfile%getvalue("post","p_fix_constant",0)
time                    = controldictfile%getvalue("mesh","current_time",0.0d0)
do i=1,maxnummeshes
      meshesfilepaths(i) = controldictfile%getvalue("meshes","file_path"//fvl_trim(fvl_char(i)),fvl_char_empty)
      meshesfileforms(i) = controldictfile%getvalue("meshes","file_form"//fvl_trim(fvl_char(i)),"binary")
end do
do i=1,maxnummeshes
      solutionsfilesdir(i) = controldictfile%getvalue("solutions","files_dir"//fvl_trim(fvl_char(i)),fvl_char_empty)
      solutionsfilesdir(i) = fvl_trim(solutionsfilesdir(i))//fvl_trim(fvl_shortchar(time))
      solutionsfilesform(i) = controldictfile%getvalue("solutions","files_form"//fvl_trim(fvl_char(i)),"binary")
end do
do i=1,maxnummeshes
      paramsfilesdir(i) = controldictfile%getvalue("params","files_dir"//fvl_trim(fvl_char(i)),fvl_char_empty)
      paramsfilesform(i) = controldictfile%getvalue("params","files_form"//fvl_trim(fvl_char(i)),"ascii")
end do
do i=1,maxnummeshes
      errorsfilesfir(i) = controldictfile%getvalue("errors","files_dir"//fvl_trim(fvl_char(i)),fvl_char_empty)
      errorsfilesfir(i) = fvl_trim(errorsfilesfir(i))//fvl_trim(fvl_shortchar(time))
      errorsfilesform(i) = controldictfile%getvalue("errors","files_form"//fvl_trim(fvl_char(i)),"binary")
end do
do i=1,maxnummeshes
      refmeshesfilepaths(i) = controldictfile%getvalue("ref_meshes","file_path"//fvl_trim(fvl_char(i)),fvl_char_empty)
      refmeshesfileforms(i) = controldictfile%getvalue("ref_meshes","file_form"//fvl_trim(fvl_char(i)),"binary")
end do
do i=1,maxnummeshes
      refsolutionsfilesdir(i) = controldictfile%getvalue("ref_solutions","files_dir"//fvl_trim(fvl_char(i)),fvl_char_empty)
      refsolutionsfilesdir(i) = fvl_trim(refsolutionsfilesdir(i))//fvl_trim(fvl_shortchar(time))
      refsolutionsfilesform(i) = controldictfile%getvalue("ref_solutions","files_form"//fvl_trim(fvl_char(i)),"binary")
end do
do i=1,maxnummeshes
      refparamsfilesdir(i) = controldictfile%getvalue("ref_params","files_dir"//fvl_trim(fvl_char(i)),fvl_char_empty)
      refparamsfilesform(i) = controldictfile%getvalue("ref_params","files_form"//fvl_trim(fvl_char(i)),"ascii")
end do
tablesfilesdir          = controldictfile%getvalue("tables","files_dir",fvl_char_empty)
tablesprecision         = controldictfile%getvalue("tables","precision",2)

! set time directories precision
call fvl_settimedirectoryprecision(solutiontimeprecision)

! ============================================================================
! INITIALIZE MESHES
! ============================================================================

do i=1,maxnummeshes

      ! ============================================================================
      ! INITIALIZE MESHES
      ! ============================================================================

      ! state message
      call fvl_loginfo("Initializing meshes...")

      ! ----------------------------------------------------------------------------
      ! mesh
      ! ----------------------------------------------------------------------------

      ! inspect file
      inquire(file=fvl_trim(meshesfilepaths(i)),exist=exist1)

      ! files exist
      if(exist1) then

            ! check mesh path
            if(i>1) then
                  if(fvl_trim(meshesfilepaths(i))==fvl_trim(meshesfilepaths(i-1))) then
                        go to 1000
                  end if
            end if

            ! state message
            call fvl_loginfo("  >>> Reading mesh")

            ! read mesh
            primalmesh = fvl_mesh3d(filepath=meshesfilepaths(i),fileform=meshesfileforms(i),dictfile=controldictfile,dictlabel="mesh")

            ! write mesh report
            call primalmesh % report("- mesh")

            ! state message
            call fvl_loginfo("  >>> Making diamond mesh")

            ! create diamond mesh
            diamondmesh = fvl_diamond_mesh3d(primalmesh,dictfile=controldictfile,dictlabel="mesh")

            ! write mesh report
            call diamondmesh % report("- diamond mesh")

            ! initialize variables
            primal_numcells(i) = primalmesh%getnumcells()
            diamond_numcells(i) = diamondmesh%getnumcells()
            primalmeshvolumes(i) = fvl_omp_sum(primalmesh%getcellvolumes())
            diamondmeshvolumes(i) = fvl_omp_sum(diamondmesh%getcellvolumes())

            ! continue
            1000 continue

      else

            ! initialize variables
            primal_numcells(i) = 0
            diamond_numcells(i) = 0
            t_error1(i) = -1.0d0
            t_error2(i) = -1.0d0
            t_errorinf(i) = -1.0d0
            t_order1(i) = -1.0d0
            t_order2(i) = -1.0d0
            t_orderinf(i) = -1.0d0
            p_error1(i) = -1.0d0
            p_error2(i) = -1.0d0
            p_errorinf(i) = -1.0d0
            p_order1(i) = -1.0d0
            p_order2(i) = -1.0d0
            p_orderinf(i) = -1.0d0
            ux_error1(i) = -1.0d0
            ux_error2(i) = -1.0d0
            ux_errorinf(i) = -1.0d0
            ux_order1(i) = -1.0d0
            ux_order2(i) = -1.0d0
            ux_orderinf(i) = -1.0d0
            uy_error1(i) = -1.0d0
            uy_error2(i) = -1.0d0
            uy_errorinf(i) = -1.0d0
            uy_order1(i) = -1.0d0
            uy_order2(i) = -1.0d0
            uy_orderinf(i) = -1.0d0
            uz_error1(i) = -1.0d0
            uz_error2(i) = -1.0d0
            uz_errorinf(i) = -1.0d0
            uz_order1(i) = -1.0d0
            uz_order2(i) = -1.0d0
            uz_orderinf(i) = -1.0d0
            u_error1(i) = -1.0d0
            u_error2(i) = -1.0d0
            u_errorinf(i) = -1.0d0
            u_order1(i) = -1.0d0
            u_order2(i) = -1.0d0
            u_orderinf(i) = -1.0d0
            diamond_ux_error1(i) = -1.0d0
            diamond_ux_error2(i) = -1.0d0
            diamond_ux_errorinf(i) = -1.0d0
            diamond_ux_order1(i) = -1.0d0
            diamond_ux_order2(i) = -1.0d0
            diamond_ux_orderinf(i) = -1.0d0
            diamond_uy_error1(i) = -1.0d0
            diamond_uy_error2(i) = -1.0d0
            diamond_uy_errorinf(i) = -1.0d0
            diamond_uy_order1(i) = -1.0d0
            diamond_uy_order2(i) = -1.0d0
            diamond_uy_orderinf(i) = -1.0d0
            diamond_uz_error1(i) = -1.0d0
            diamond_uz_error2(i) = -1.0d0
            diamond_uz_errorinf(i) = -1.0d0
            diamond_uz_order1(i) = -1.0d0
            diamond_uz_order2(i) = -1.0d0
            diamond_uz_orderinf(i) = -1.0d0
            diamond_u_error1(i) = -1.0d0
            diamond_u_error2(i) = -1.0d0
            diamond_u_errorinf(i) = -1.0d0
            diamond_u_order1(i) = -1.0d0
            diamond_u_order2(i) = -1.0d0
            diamond_u_orderinf(i) = -1.0d0
            meshestotaltime(i) = -1.0d0
            modelstotaltime(i) = -1.0d0
            schemestotaltime(i) = -1.0d0
            solutiontotaltime(i) = -1.0d0
            postproctotaltime(i) = -1.0d0
            executiontotaltime(i) = -1.0d0
            meshestotaltimeratio(i) = -1.0d0
            modelstotaltimeratio(i) = -1.0d0
            schemestotaltimeratio(i) = -1.0d0
            solutiontotaltimeratio(i) = -1.0d0
            postproctotaltimeratio(i) = -1.0d0
            executiontotaltimeratio(i) = -1.0d0
            reftotalexecutionmemory(i) = -1.0d0
            linearsystemnumiters(i) = -1
            linearsystemresid(i) = -1.0d0
            fixpointnumiters(i) = -1
            fixpointresid(i) = -1.0d0
            cycle

      end if

      ! ----------------------------------------------------------------------------
      ! reference mesh
      ! ----------------------------------------------------------------------------

      ! inspect file
      inquire(file=fvl_trim(refmeshesfilepaths(i)),exist=exist1)

      ! files exist
      if(exist1 .and. (.not. associated(t_valuefun) .or. .not. associated(p_valuefun) .or. .not. associated(u_valuefun))) then

            if(i>1) then
                  if(fvl_trim(refmeshesfilepaths(i))==fvl_trim(refmeshesfilepaths(i-1))) then
                        go to 2000
                  end if
            end if

            ! state message
            call fvl_loginfo("  >>> Reading reference mesh")

            ! read mesh
            primalrefmesh = fvl_mesh3d(filepath=refmeshesfilepaths(i),fileform=refmeshesfileforms(i),quadratorder=8)

            ! set current time
            call primalrefmesh % setcurrenttime(time)

            ! write mesh report
            call primalrefmesh % report("- reference mesh")

            ! state message
            call fvl_loginfo("  >>> Making reference diamond mesh")

            ! create diamond mesh
            diamondrefmesh = fvl_diamond_mesh3d(primalrefmesh,quadratorder=8)

            ! set current time
            call diamondrefmesh % setcurrenttime(time)

            ! write mesh report
            call diamondrefmesh % report("- reference diamond mesh")

            ! state message
            call fvl_loginfo("  >>> Computing mesh locators")

            ! mesh locator
            primalmeshlocator = fvl_meshlocator3d(primalrefmesh,controldictfile,"mesh_locator")
            diamondmeshlocator = fvl_meshlocator3d(diamondrefmesh,controldictfile,"mesh_locator")

            ! mesh patches
            primal_patch = fvl_patch3d(primalrefmesh)
            diamond_patch = fvl_patch3d(diamondrefmesh)

            ! ! state message
            ! call fvl_loginfo("  >>> Computing reconstruction mappers")
            !
            ! reconstruction map
            ! primal_reconstmap = fvl_reconstmap3d(primal_patch,primalrefmesh,primalmeshlocator,controldictfile,"reconst_map")
            ! diamond_reconstmap = fvl_reconstmap3d(diamond_patch,diamondrefmesh,diamondmeshlocator,controldictfile,"reconst_map")

            ! continue
            2000 continue

            ! reference mesh exists
            exist2 = .true.

      else

            ! reference mesh exists
            exist2 = .false.

      end if

      ! ----------------------------------------------------------------------------
      ! user-defined mesh
      ! ----------------------------------------------------------------------------

      ! state message
      call fvl_loginfo("  >>> Setting user-defined mesh")

! user-defined code
#include "meshes.f"

      ! ============================================================================
      ! INITIALIZE PHYSICAL MODEL
      ! ============================================================================

      ! state message
      call fvl_loginfo("Initializing physical models...")

! user-defined code
#include "models.f"

      ! ----------------------------------------------------------------------------
      ! exact temperature
      ! ----------------------------------------------------------------------------

      ! suffixes
      if(fvl_trim(t_label)/=fvl_char_empty) then
            if(fvl_trim(refsolutionsfilesform(i))=="ascii") then
                  t_suffix = fvl_trim(refsolutionsfilesdir(i))//"/"//fvl_trim(t_label)//".fvd3"
            else
                  t_suffix = fvl_trim(refsolutionsfilesdir(i))//"/"//fvl_trim(t_label)//".fvd3.bin"
            end if
      else
            go to 3000
      end if

      ! inspect files
      inquire(file=fvl_trim(t_suffix),exist=exist1)

      ! temperature file exists
      if(associated(t_valuefun)) then

            ! state message
            call fvl_loginfo("  >>> Computing exact temperature solution")

            ! allocate memory
            call fvl_allocate(t_exactvals,primal_numcells(i))

            ! compute exact values
            call primalmesh % evalcellmeanvalues(t_valuefun,t_exactvals)

            ! state message
            call fvl_loginfo("  >>> Writing exact temperature solution")

            ! suffixes
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  t_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(t_label)//"_exact.fvd3"
            else
                  t_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(t_label)//"_exact.fvd3.bin"
            end if

            ! write file
            solutionsfile = fvl_data_file(path=fvl_trim(t_suffix),form=errorsfilesform(i),action="write")
            call solutionsfile % open()
            call solutionsfile % writedata(t_exactvals,fvl_trim(t_label)//" Exact",fvl_data_file_cellmeanvalues,times(i),steps(i))
            call solutionsfile % close()

      else if(exist1 .and. exist2) then

            if(i>1) then
                  if(fvl_trim(refsolutionsfilesdir(i))==fvl_trim(refsolutionsfilesdir(i-1))) then
                        go to 3000
                  end if
            end if

            ! state message
            call fvl_loginfo("  >>> Reading reference temperature solution")

            ! read file
            solutionsfile = fvl_data_file(path=fvl_trim(t_suffix),form=refsolutionsfilesform(i),action="read")
            call solutionsfile % open()
            element = 0
            do while(element/=fvl_data_file_cellmeanvalues)
                  call solutionsfile % readdata(res=t_refvals,map=t_map,element=element,time=times(i),step=steps(i))
            end do
            call solutionsfile % close()

            ! check if data exists
            if(.not. allocated(t_refvals)) then
                  call fvl_logerror("Reference temperature solution not found")
            end if

            ! state message
            call fvl_loginfo("  >>> Computing temperature reconstruction map")

            ! reconstruction map
            t_reconstmap = fvl_reconstmap3d(primal_patch,primalrefmesh,primalmeshlocator,controldictfile,"reconst_map",t_refvals)

            ! state message
            3000 continue
            call fvl_loginfo("  >>> Computing exact temperature solution")

            ! allocate memory
            call fvl_allocate(t_exactvals,primal_numcells(i))

            ! map reference solution
            call t_reconstmap % mapcellmeanvalues(primalmesh,t_exactvals)

            ! state message
            call fvl_loginfo("  >>> Writing exact temperature solution")

            ! suffixes
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  t_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(t_label)//"_exact.fvd3"
            else
                  t_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(t_label)//"_exact.fvd3.bin"
            end if

            ! write file
            solutionsfile = fvl_data_file(path=fvl_trim(t_suffix),form=errorsfilesform(i),action="write")
            call solutionsfile % open()
            call solutionsfile % writedata(t_exactvals,fvl_trim(t_label)//" Exact",fvl_data_file_cellmeanvalues,times(i),steps(i))
            call solutionsfile % close()

      else

            ! allocate memory
            call fvl_allocate(t_exactvals,primal_numcells(i))

            ! clean memory
            call fvl_omp_clean(t_exactvals)

      end if

      ! continue
      3001 continue

      ! ----------------------------------------------------------------------------
      ! exact pressure
      ! ----------------------------------------------------------------------------

      ! suffixes
      if(fvl_trim(p_label)/=fvl_char_empty) then
            if(fvl_trim(refsolutionsfilesform(i))=="ascii") then
                  p_suffix = fvl_trim(refsolutionsfilesdir(i))//"/"//fvl_trim(p_label)//".fvd3"
            else
                  p_suffix = fvl_trim(refsolutionsfilesdir(i))//"/"//fvl_trim(p_label)//".fvd3.bin"
            end if
      else
            go to 4001
      end if

      ! inspect files
      inquire(file=fvl_trim(p_suffix),exist=exist1)

      ! pressure file exists
      if(associated(p_valuefun)) then

            ! state message
            call fvl_loginfo("  >>> Computing exact pressure solution")

            ! allocate memory
            call fvl_allocate(p_exactvals,primal_numcells(i))

            ! compute exact values
            call primalmesh % evalcellmeanvalues(p_valuefun,p_exactvals)

            ! state message
            call fvl_loginfo("  >>> Writing exact pressure solution")

            ! suffixes
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  p_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(p_label)//"_exact.fvd3"
            else
                  p_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(p_label)//"_exact.fvd3.bin"
            end if

            ! write file
            solutionsfile = fvl_data_file(path=fvl_trim(p_suffix),form=errorsfilesform(i),action="write")
            call solutionsfile % open()
            call solutionsfile % writedata(p_exactvals,fvl_trim(p_label)//" Exact",fvl_data_file_cellmeanvalues,times(i),steps(i))
            call solutionsfile % close()

      else if(exist1 .and. exist2) then

            if(i>1) then
                  if(fvl_trim(refsolutionsfilesdir(i))==fvl_trim(refsolutionsfilesdir(i-1))) then
                        go to 4000
                  end if
            end if

            ! state message
            call fvl_loginfo("  >>> Reading reference pressure solution")

            ! read file
            solutionsfile = fvl_data_file(path=fvl_trim(p_suffix),form=refsolutionsfilesform(i),action="read")
            call solutionsfile % open()
            element = 0
            do while(element/=fvl_data_file_cellmeanvalues)
                  call solutionsfile % readdata(res=p_refvals,map=p_map,element=element,time=times(i),step=steps(i))
            end do
            call solutionsfile % close()

            ! check if data exists
            if(.not. allocated(p_refvals)) then
                  call fvl_logerror("Reference pressure solution not found")
            end if

            ! state message
            call fvl_loginfo("  >>> Computing pressure reconstruction map")

            ! reconstruction map
            p_reconstmap = fvl_reconstmap3d(primal_patch,primalrefmesh,primalmeshlocator,controldictfile,"reconst_map",p_refvals)

            ! state message
            4000 continue
            call fvl_loginfo("  >>> Computing exact pressure solution")

            ! allocate memory
            call fvl_allocate(p_exactvals,primal_numcells(i))

            ! map reference solution
            call p_reconstmap % mapcellmeanvalues(primalmesh,p_exactvals)

            ! state message
            call fvl_loginfo("  >>> Writing exact pressure solution")

            ! suffixes
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  p_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(p_label)//"_exact.fvd3"
            else
                  p_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(p_label)//"_exact.fvd3.bin"
            end if

            ! write file
            solutionsfile = fvl_data_file(path=fvl_trim(p_suffix),form=errorsfilesform(i),action="write")
            call solutionsfile % open()
            call solutionsfile % writedata(p_exactvals,fvl_trim(p_label)//" Exact",fvl_data_file_cellmeanvalues,times(i),steps(i))
            call solutionsfile % close()

      else

            ! allocate memory
            call fvl_allocate(p_exactvals,primal_numcells(i))

            ! clean memory
            call fvl_omp_clean(p_exactvals)

      end if

      ! continue
      4001 continue

      ! ----------------------------------------------------------------------------
      ! exact velocity
      ! ----------------------------------------------------------------------------

      ! suffixes
      if(fvl_trim(u_label)/=fvl_char_empty) then
            if(refsolutionsfilesform(i)=="ascii") then
                  u_suffix = fvl_trim(refsolutionsfilesdir(i))//"/"//fvl_trim(u_label)//".fvd3"
            else
                  u_suffix = fvl_trim(refsolutionsfilesdir(i))//"/"//fvl_trim(u_label)//".fvd3.bin"
            end if
      else
            go to 5001
      end if

      ! inspect files
      inquire(file=fvl_trim(u_suffix),exist=exist1)

      ! velocity file exists
      if(associated(u_valuefun)) then

            ! state message
            call fvl_loginfo("  >>> Computing exact velocity solution")

            ! allocate memory
            call fvl_allocate(u_exactvals,3,primal_numcells(i))

            ! compute exact values
            call primalmesh % evalcellmeanvalues(u_valuefun,u_exactvals)

            ! state message
            call fvl_loginfo("  >>> Writing exact velocity solution")

            ! suffixes
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_exact.fvd3"
            else
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_exact.fvd3.bin"
            end if

            ! write file
            solutionsfile = fvl_data_file(path=fvl_trim(u_suffix),form=errorsfilesform(i),action="write")
            call solutionsfile % open()
            call solutionsfile % writedata(u_exactvals,fvl_trim(u_label)//" Exact",fvl_data_file_cellmeanvalues,times(i),steps(i))
            call solutionsfile % close()

      else if(exist1 .and. exist2) then

            if(i>1) then
                  if(fvl_trim(refsolutionsfilesdir(i))==fvl_trim(refsolutionsfilesdir(i-1))) then
                        go to 5000
                  end if
            end if

            ! state message
            call fvl_loginfo("  >>> Reading reference velocity solution")

            ! read file
            solutionsfile = fvl_data_file(path=fvl_trim(u_suffix),form=refsolutionsfilesform(i),action="read")
            call solutionsfile % open()
            element = 0
            do while(element/=fvl_data_file_cellmeanvalues)
                  call solutionsfile % readdata(res=u_refvals,map=u_map,element=element,time=times(i),step=steps(i))
            end do
            call solutionsfile % close()

            ! check if data exists
            if(.not. allocated(u_refvals1)) then
                  call fvl_logerror("Reference velocity solution not found")
            end if

            ! allocate memory
            call fvl_allocate(u_refvals2,3*primalrefmesh%getnumcells())

            ! reshape data
            u_refvals2 = reshape(transpose(u_refvals1),[3*primalrefmesh%getnumcells()],order=[1])

            ! state message
            call fvl_loginfo("  >>> Computing velocity reconstruction map")

            ! reconstruction map
            ux_reconstmap = fvl_reconstmap3d(primal_patch,primalrefmesh,primalmeshlocator,controldictfile,"reconst_map",&
                  u_refvals2(1:primalrefmesh%getnumcells()))
            uy_reconstmap = fvl_reconstmap3d(primal_patch,primalrefmesh,primalmeshlocator,controldictfile,"reconst_map",&
                  u_refvals2(primalrefmesh%getnumcells()+1:2*primalrefmesh%getnumcells()))
            uz_reconstmap = fvl_reconstmap3d(primal_patch,primalrefmesh,primalmeshlocator,controldictfile,"reconst_map",&
                  u_refvals2(2*primalrefmesh%getnumcells()+1:3*primalrefmesh%getnumcells()))

            ! state message
            5000 continue
            call fvl_loginfo("  >>> Computing exact velocity solution")

            ! allocate memory
            call fvl_allocate(u_exactvals,3,primal_numcells(i))
            call fvl_allocate(auxv,3*primal_numcells(i))

            ! map reference solution
            call ux_reconstmap % mapcellmeanvalues(primalmesh,auxv(1:primal_numcells(i)))
            call uy_reconstmap % mapcellmeanvalues(primalmesh,auxv(primal_numcells(i)+1:2*primal_numcells(i)))
            call uz_reconstmap % mapcellmeanvalues(primalmesh,auxv(2*primal_numcells(i)+1:3*primal_numcells(i)))
            u_exactvals = reshape(auxv,[3,primal_numcells(i)],order=[2,1])

            ! state message
            call fvl_loginfo("  >>> Writing exact velocity solution")

            ! suffixes
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_exact.fvd3"
            else
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_exact.fvd3.bin"
            end if

            ! write file
            solutionsfile = fvl_data_file(path=fvl_trim(u_suffix),form=errorsfilesform(i),action="write")
            call solutionsfile % open()
            call solutionsfile % writedata(u_exactvals,fvl_trim(u_label)//" Exact",fvl_data_file_cellmeanvalues,times(i),steps(i))
            call solutionsfile % close()

      else

            ! allocate memory
            call fvl_allocate(u_exactvals,3,primal_numcells(i))

            ! clean memory
            call fvl_omp_clean(u_exactvals)

      end if

      ! continue
      5001 continue

      ! ----------------------------------------------------------------------------
      ! exact diamond velocity
      ! ----------------------------------------------------------------------------

      ! suffixes
      if(fvl_trim(u_label)/=fvl_char_empty) then
            if(fvl_trim(refsolutionsfilesform(i))=="ascii") then
                  u_suffix = fvl_trim(refsolutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_diamond.fvd3"
            else
                  u_suffix = fvl_trim(refsolutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_diamond.fvd3.bin"
            end if
      else
            go to 6001
      end if

      ! inspect files
      inquire(file=fvl_trim(u_suffix),exist=exist1)

      ! velocity file exists
      if(associated(u_valuefun)) then

            ! state message
            call fvl_loginfo("  >>> Computing exact diamond velocity solution")

            ! allocate memory
            call fvl_allocate(diamond_u_exactvals,3,diamond_numcells(i))

            ! compute exact values
            call diamondmesh % evalcellmeanvalues(u_valuefun,diamond_u_exactvals)

            ! state message
            call fvl_loginfo("  >>> Writing exact diamond velocity solution")

            ! suffixes
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_diamond_exact.fvd3"
            else
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_diamond_exact.fvd3.bin"
            end if

            ! write file
            solutionsfile = fvl_data_file(path=fvl_trim(u_suffix),form=errorsfilesform(i),action="write")
            call solutionsfile % open()
            call solutionsfile % writedata(diamond_u_exactvals,fvl_trim(u_label)//" Exact",fvl_data_file_cellmeanvalues,times(i),steps(i))
            call solutionsfile % close()

      else if(exist1 .and. exist2) then

            if(i>1) then
                  if(fvl_trim(refsolutionsfilesdir(i))==fvl_trim(refsolutionsfilesdir(i-1))) then
                        go to 6000
                  end if
            end if

            ! state message
            call fvl_loginfo("  >>> Reading reference diamond velocity solution")

            ! read file
            solutionsfile = fvl_data_file(path=fvl_trim(u_suffix),form=refsolutionsfilesform(i),action="read")
            call solutionsfile % open()
            element = 0
            do while(element/=fvl_data_file_cellmeanvalues)
                  call solutionsfile % readdata(res=dimaond_u_refvals,map=diamond_u_map,element=element,time=times(i),step=steps(i))
            end do
            call solutionsfile % close()

            ! check if data exists
            if(.not. allocated(diamond_u_refvals1)) then
                  call fvl_logerror("Reference diamond velocity solution not found")
            end if

            ! allocate memory
            call fvl_allocate(diamond_u_refvals2,3*diamondrefmesh%getnumcells())

            ! reshape data
            diamond_u_refvals2 = reshape(transpose(diamond_u_refvals1),[3*diamondrefmesh%getnumcells()],order=[1])

            ! state message
            call fvl_loginfo("  >>> Computing diamond velocity reconstruction map")

            ! reconstruction map
            diamond_ux_reconstmap = fvl_reconstmap3d(diamond_patch,diamondrefmesh,diamondmeshlocator,controldictfile,"reconst_map",&
                  diamond_u_refvals2(1:diamondrefmesh%getnumcells()))
            diamond_uy_reconstmap = fvl_reconstmap3d(diamond_patch,diamondrefmesh,diamondmeshlocator,controldictfile,"reconst_map",&
                  diamond_u_refvals2(diamondrefmesh%getnumcells()+1:2*diamondrefmesh%getnumcells()))
            diamond_uz_reconstmap = fvl_reconstmap3d(diamond_patch,diamondrefmesh,diamondmeshlocator,controldictfile,"reconst_map",&
                  diamond_u_refvals2(2*diamondrefmesh%getnumcells()+1:3*diamondrefmesh%getnumcells()))

            ! state message
            6000 continue
            call fvl_loginfo("  >>> Computing exact diamond velocity solution")

            ! allocate memory
            call fvl_allocate(diamond_u_exactvals,3,diamond_numcells(i))
            call fvl_allocate(auxv,3*diamond_numcells(i))

            ! map reference solution
            call diamond_ux_reconstmap % mapcellmeanvalues(diamondmesh,auxv(1:diamond_numcells(i)))
            call diamond_uy_reconstmap % mapcellmeanvalues(diamondmesh,auxv(diamond_numcells(i)+1:2*diamond_numcells(i)))
            call diamond_uz_reconstmap % mapcellmeanvalues(diamondmesh,auxv(2*diamond_numcells(i)+1:3*diamond_numcells(i)))
            diamond_u_exactvals = reshape(auxv,[3,diamond_numcells(i)],order=[2,1])

            ! state message
            call fvl_loginfo("  >>> Writing exact diamond velocity solution")

            ! suffixes
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_diamond_exact.fvd3"
            else
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_diamond_exact.fvd3.bin"
            end if

            ! write file
            solutionsfile = fvl_data_file(path=fvl_trim(u_suffix),form=errorsfilesform(i),action="write")
            call solutionsfile % open()
            call solutionsfile % writedata(diamond_u_exactvals,fvl_trim(u_label)//" Exact",fvl_data_file_cellmeanvalues,times(i),steps(i))
            call solutionsfile % close()

      else

            ! allocate memory
            call fvl_allocate(diamond_u_exactvals,3,diamond_numcells(i))

            ! clean memory
            call fvl_omp_clean(diamond_u_exactvals)

      end if

      ! continue
      6001 continue

      ! ============================================================================
      ! POST-PROCESSING DATA
      ! ============================================================================

      ! state message
      call fvl_loginfo("Post-processing data...")

      ! ----------------------------------------------------------------------------
      ! approximate temperature
      ! ----------------------------------------------------------------------------

      ! suffixes
      if(fvl_trim(t_label)/=fvl_char_empty) then
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  t_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(t_label)//".fvd3"
            else
                  t_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(t_label)//".fvd3.bin"
            end if
      else
            go to 7000
      end if

      ! inspect files
      inquire(file=fvl_trim(t_suffix),exist=exist1)

      ! temperature file exists
      if(exist1) then

            ! state message
            call fvl_loginfo("  >>> Reading approximate temperature solution")

            ! read file
            solutionsfile = fvl_data_file(path=fvl_trim(t_suffix),form=solutionsfilesform(i),action="read")
            call solutionsfile % open()
            element = 0
            do while(element/=fvl_data_file_cellmeanvalues)
                  call solutionsfile % readdata(res=t_approxvals,map=t_map,element=element,time=times(i),step=steps(i))
            end do
            call solutionsfile % close()

            ! check if data exists
            if(.not. allocated(t_approxvals)) then
                  call fvl_logerror("Approximate temperature solution not found")
            end if

            ! state message
            call fvl_loginfo("  >>> Computing approximate temperature errors")

            ! allocate memory
            call fvl_allocate(t_errorvals,primal_numcells(i))
            call fvl_omp_clean(t_errorvals)

            ! compute errors
            do j=1,size(t_map,1)
                  t_errorvals(t_map(j)) = t_approxvals(j)-t_exactvals(t_map(j))
            end do

            ! compute error norms
            t_error1(i) = 0.0d0
            t_error2(i) = 0.0d0
            t_errorinf(i) = 0.0d0
            do j=1,primal_numcells(i)
                  cellvolume = primalmesh%getcellvolume(j)
                  t_error1(i) = t_error1(i)+abs(t_errorvals(j))*cellvolume
                  t_error2(i) = t_error2(i)+(t_errorvals(j)**2.0)*cellvolume
                  t_errorinf(i) = max(t_errorinf(i),abs(t_errorvals(j)))
            end do
            t_error1(i) = t_error1(i)/primalmeshvolumes(i)
            t_error2(i) = sqrt(t_error2(i)/primalmeshvolumes(i))

            ! compute convergence orders
            if(i>1) then
                  if(convergence==1) then
                        if(primal_numcells(i-1)>0 .and. primal_numcells(i-1)/=primal_numcells(i) .and. t_errorinf(i)>0.0d0) then
                              t_order1(i) = 3.0d0*log(t_error1(i-1)/t_error1(i))&
                                    /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                              t_order2(i) = 3.0d0*log(t_error2(i-1)/t_error2(i))&
                                    /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                              t_orderinf(i) = 3.0d0*log(t_errorinf(i-1)/t_errorinf(i))&
                                    /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                        else
                              t_order1(i) = -1.0d0
                              t_order2(i) = -1.0d0
                              t_orderinf(i) = -1.0d0
                        end if
                  else
                        if(steps(i-1)>0 .and. steps(i-1)/=steps(i) .and. t_errorinf(i)>0.0d0) then
                              t_order1(i) = log(t_error1(i-1)/t_error1(i))&
                                    /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              t_order2(i) = log(t_error2(i-1)/t_error2(i))&
                                    /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              t_orderinf(i) = log(t_errorinf(i-1)/t_errorinf(i))&
                                    /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                        else
                              t_order1(i) = -1.0d0
                              t_order2(i) = -1.0d0
                              t_orderinf(i) = -1.0d0
                        end if
                  end if
            else
                  t_order1(i) = -1.0d0
                  t_order2(i) = -1.0d0
                  t_orderinf(i) = -1.0d0
            end if

            ! state message
            call fvl_loginfo("  >>> Writing approximate temperature errors")

            ! suffixes
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  t_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(t_label)//"_error.fvd3"
            else
                  t_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(t_label)//"_error.fvd3.bin"
            end if

            ! write file
            solutionsfile = fvl_data_file(path=fvl_trim(t_suffix),form=errorsfilesform(i),action="write")
            call solutionsfile % open()
            call solutionsfile % writedata(t_errorvals,t_map,fvl_trim(t_label)//" Error",fvl_data_file_cellmeanvalues,times(i),steps(i))
            call solutionsfile % close()

      else

            t_error1(i) = -1.0d0
            t_error2(i) = -1.0d0
            t_errorinf(i) = -1.0d0
            t_order1(i) = -1.0d0
            t_order2(i) = -1.0d0
            t_orderinf(i) = -1.0d0

      end if

      ! continue
      7000 continue

      ! ----------------------------------------------------------------------------
      ! approximate pressure
      ! ----------------------------------------------------------------------------

      ! suffixes
      if(fvl_trim(p_label)/=fvl_char_empty) then
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  p_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(p_label)//".fvd3"
            else
                  p_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(p_label)//".fvd3.bin"
            end if
      else
            go to 8000
      end if

      ! inspect files
      inquire(file=fvl_trim(p_suffix),exist=exist1)

      ! pressure file exists
      if(exist1) then

            ! state message
            call fvl_loginfo("  >>> Reading approximate pressure solution")

            ! read file
            solutionsfile = fvl_data_file(path=fvl_trim(p_suffix),form=solutionsfilesform(i),action="read")
            call solutionsfile % open()
            element = 0
            do while(element/=fvl_data_file_cellmeanvalues)
                  call solutionsfile % readdata(res=p_approxvals,map=p_map,element=element,time=times(i),step=steps(i))
            end do
            call solutionsfile % close()

            ! check if data exists
            if(.not. allocated(p_approxvals)) then
                  call fvl_logerror("Approximate pressure solution not found")
            end if

            ! state message
            call fvl_loginfo("  >>> Computing approximate pressure errors")

            ! allocate memory
            call fvl_allocate(p_errorvals,primal_numcells(i))
            call fvl_omp_clean(p_errorvals)

            ! fix constant
            if(p_fixconstant==1) then
                  call fvl_omp_sub(p_exactvals,fvl_omp_dot(p_exactvals,primalmesh%getcellvolumes())/primalmeshvolumes(i))
                  call fvl_omp_sub(p_approxvals,fvl_omp_dot(p_approxvals,primalmesh%getcellvolumes())/primalmeshvolumes(i))
            end if

            ! compute errors
            do j=1,size(p_map,1)
                  p_errorvals(p_map(j)) = p_approxvals(j)-p_exactvals(p_map(j))
            end do

            ! compute error norms
            p_error1(i) = 0.0d0
            p_error2(i) = 0.0d0
            p_errorinf(i) = 0.0d0
            do j=1,primal_numcells(i)
                  cellvolume = primalmesh%getcellvolume(j)
                  p_error1(i) = p_error1(i)+abs(p_errorvals(j))*cellvolume
                  p_error2(i) = p_error2(i)+(p_errorvals(j)**2.0)*cellvolume
                  p_errorinf(i) = max(p_errorinf(i),abs(p_errorvals(j)))
            end do
            p_error1(i) = p_error1(i)/primalmeshvolumes(i)
            p_error2(i) = sqrt(p_error2(i)/primalmeshvolumes(i))

            ! compute convergence orders
            if(i>1) then
                  if(convergence==1) then
                        if(primal_numcells(i-1)>0 .and. primal_numcells(i-1)/=primal_numcells(i) .and. p_errorinf(i)>0.0d0) then
                              p_order1(i) = 3.0d0*log(p_error1(i-1)/p_error1(i))&
                                    /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                              p_order2(i) = 3.0d0*log(p_error2(i-1)/p_error2(i))&
                                    /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                              p_orderinf(i) = 3.0d0*log(p_errorinf(i-1)/p_errorinf(i))&
                                    /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                        else
                              p_order1(i) = -1.0d0
                              p_order2(i) = -1.0d0
                              p_orderinf(i) = -1.0d0
                        end if
                  else
                        if(steps(i-1)>0 .and. steps(i-1)/=steps(i) .and. p_errorinf(i)>0.0d0) then
                              p_order1(i) = log(p_error1(i-1)/p_error1(i))&
                                    /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              p_order2(i) = log(p_error2(i-1)/p_error2(i))&
                                    /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              p_orderinf(i) = log(p_errorinf(i-1)/p_errorinf(i))&
                                    /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                        else
                              p_order1(i) = -1.0d0
                              p_order2(i) = -1.0d0
                              p_orderinf(i) = -1.0d0
                        end if
                  end if
            else
                  p_order1(i) = -1.0d0
                  p_order2(i) = -1.0d0
                  p_orderinf(i) = -1.0d0
            end if

            ! state message
            call fvl_loginfo("  >>> Writing approximate pressure errors")

            ! suffixes
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  p_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(p_label)//"_error.fvd3"
            else
                  p_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(p_label)//"_error.fvd3.bin"
            end if

            ! write file
            solutionsfile = fvl_data_file(path=fvl_trim(p_suffix),form=errorsfilesform(i),action="write")
            call solutionsfile % open()
            call solutionsfile % writedata(p_errorvals,p_map,fvl_trim(p_label)//" Error",fvl_data_file_cellmeanvalues,times(i),steps(i))
            call solutionsfile % close()

      else

            p_error1(i) = -1.0d0
            p_error2(i) = -1.0d0
            p_errorinf(i) = -1.0d0
            p_order1(i) = -1.0d0
            p_order2(i) = -1.0d0
            p_orderinf(i) = -1.0d0

      end if

      ! continue
      8000 continue

      ! ----------------------------------------------------------------------------
      ! approximate velocity
      ! ----------------------------------------------------------------------------

      ! suffixes
      if(fvl_trim(u_label)/=fvl_char_empty) then
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//".fvd3"
            else
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//".fvd3.bin"
            end if
      else
            go to 9000
      end if

      ! inspect files
      inquire(file=fvl_trim(u_suffix),exist=exist1)

      ! velocity file exists
      if(exist1) then

            ! state message
            call fvl_loginfo("  >>> Reading approximate velocity solution")

            ! read file
            solutionsfile = fvl_data_file(path=fvl_trim(u_suffix),form=solutionsfilesform(i),action="read")
            call solutionsfile % open()
            element = 0
            do while(element/=fvl_data_file_cellmeanvalues)
                  call solutionsfile % readdata(res=u_approxvals,map=u_map,element=element,time=times(i),step=steps(i))
            end do
            call solutionsfile % close()

            ! check if data exists
            if(.not. allocated(u_approxvals)) then
                  call fvl_logerror("Approximate velocity solution not found")
            end if

            ! state message
            call fvl_loginfo("  >>> Computing approximate velocity errors")

            ! allocate memory
            call fvl_allocate(u_errorvals,2,primal_numcells(i))
            call fvl_omp_clean(u_errorvals)

            ! compute errors
            do j=1,size(u_map,1)
                  u_errorvals(:,u_map(j)) = u_approxvals(:,j)-u_exactvals(:,u_map(j))
            end do

            ! compute error norms
            ux_error1(i) = 0.0d0
            ux_error2(i) = 0.0d0
            ux_errorinf(i) = 0.0d0
            uy_error1(i) = 0.0d0
            uy_error2(i) = 0.0d0
            uy_errorinf(i) = 0.0d0
            uz_error1(i) = 0.0d0
            uz_error2(i) = 0.0d0
            uz_errorinf(i) = 0.0d0
            u_error1(i) = 0.0d0
            u_error2(i) = 0.0d0
            u_errorinf(i) = 0.0d0
            do j=1,primal_numcells(i)
                  cellvolume = primalmesh%getcellvolume(j)
                  ux_error1(i) = ux_error1(i)+abs(u_errorvals(1,j))*cellvolume
                  ux_error2(i) = ux_error2(i)+(u_errorvals(1,j)**2.0)*cellvolume
                  ux_errorinf(i) = max(ux_errorinf(i),abs(u_errorvals(1,j)))
                  uy_error1(i) = uy_error1(i)+abs(u_errorvals(2,j))*cellvolume
                  uy_error2(i) = uy_error2(i)+(u_errorvals(2,j)**2.0)*cellvolume
                  uy_errorinf(i) = max(uy_errorinf(i),abs(u_errorvals(2,j)))
                  uz_error1(i) = uz_error1(i)+abs(u_errorvals(3,j))*cellvolume
                  uz_error2(i) = uz_error2(i)+(u_errorvals(3,j)**2.0)*cellvolume
                  uz_errorinf(i) = max(uz_errorinf(i),abs(u_errorvals(3,j)))
                  aux=fvl_vectormag3(u_errorvals(:,j))
                  u_error1(i) = u_error1(i)+aux*cellvolume
                  u_error2(i) = u_error2(i)+(aux**2.0)*cellvolume
                  u_errorinf(i) = max(u_errorinf(i),aux)
            end do
            ux_error1(i) = ux_error1(i)/primalmeshvolumes(i)
            ux_error2(i) = sqrt(ux_error2(i)/primalmeshvolumes(i))
            uy_error1(i) = uy_error1(i)/primalmeshvolumes(i)
            uy_error2(i) = sqrt(uy_error2(i)/primalmeshvolumes(i))
            uz_error1(i) = uz_error1(i)/primalmeshvolumes(i)
            uz_error2(i) = sqrt(uz_error2(i)/primalmeshvolumes(i))
            u_error1(i) = u_error1(i)/primalmeshvolumes(i)
            u_error2(i) = sqrt(u_error2(i)/primalmeshvolumes(i))

            ! compute convergence orders
            if(i>1) then
                  if(convergence==1) then
                        if(primal_numcells(i-1)>0 .and. primal_numcells(i-1)/=primal_numcells(i)) then
                              if(ux_errorinf(i)>0.0d0) then
                                    ux_order1(i) = 3.0d0*log(ux_error1(i-1)/ux_error1(i))&
                                          /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                                    ux_order2(i) = 3.0d0*log(ux_error2(i-1)/ux_error2(i))&
                                          /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                                    ux_orderinf(i) = 3.0d0*log(ux_errorinf(i-1)/ux_errorinf(i))&
                                          /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                              else
                                    ux_order1(i) = -1.0d0
                                    ux_order2(i) = -1.0d0
                                    ux_orderinf(i) = -1.0d0
                              end if
                              if(uy_errorinf(i)>0.0d0) then
                                    uy_order1(i) = 3.0d0*log(uy_error1(i-1)/uy_error1(i))&
                                          /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                                    uy_order2(i) = 3.0d0*log(uy_error2(i-1)/uy_error2(i))&
                                          /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                                    uy_orderinf(i) = 3.0d0*log(uy_errorinf(i-1)/uy_errorinf(i))&
                                          /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                              else
                                    uy_order1(i) = -1.0d0
                                    uy_order2(i) = -1.0d0
                                    uy_orderinf(i) = -1.0d0
                              end if
                              if(uz_errorinf(i)>0.0d0) then
                                    uz_order1(i) = 3.0d0*log(uz_error1(i-1)/uz_error1(i))&
                                          /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                                    uz_order2(i) = 3.0d0*log(uz_error2(i-1)/uz_error2(i))&
                                          /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                                    uz_orderinf(i) = 3.0d0*log(uz_errorinf(i-1)/uz_errorinf(i))&
                                          /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                              else
                                    uz_order1(i) = -1.0d0
                                    uz_order2(i) = -1.0d0
                                    uz_orderinf(i) = -1.0d0
                              end if
                              if(u_errorinf(i)>0.0d0) then
                                    u_order1(i) = 3.0d0*log(u_error1(i-1)/u_error1(i))&
                                          /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                                    u_order2(i) = 3.0d0*log(u_error2(i-1)/u_error2(i))&
                                          /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                                    u_orderinf(i) = 3.0d0*log(u_errorinf(i-1)/u_errorinf(i))&
                                          /log(fvl_real(primal_numcells(i))/fvl_real(primal_numcells(i-1)))
                              else
                                    u_order1(i) = -1.0d0
                                    u_order2(i) = -1.0d0
                                    u_orderinf(i) = -1.0d0
                              end if
                        else
                              ux_order1(i) = -1.0d0
                              ux_order2(i) = -1.0d0
                              ux_orderinf(i) = -1.0d0
                              uy_order1(i) = -1.0d0
                              uy_order2(i) = -1.0d0
                              uy_orderinf(i) = -1.0d0
                              uz_order1(i) = -1.0d0
                              uz_order2(i) = -1.0d0
                              uz_orderinf(i) = -1.0d0
                              u_order1(i) = -1.0d0
                              u_order2(i) = -1.0d0
                              u_orderinf(i) = -1.0d0
                        end if
                  else
                        if(steps(i-1)>0 .and. steps(i-1)/=steps(i)) then
                              if(ux_errorinf(i)>0.0d0) then
                                    ux_order1(i) = log(ux_error1(i-1)/ux_error1(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    ux_order2(i) = log(ux_error2(i-1)/ux_error2(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    ux_orderinf(i) = log(ux_errorinf(i-1)/ux_errorinf(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              else
                                    ux_order1(i) = -1.0d0
                                    ux_order2(i) = -1.0d0
                                    ux_orderinf(i) = -1.0d0
                              end if
                              if(uy_errorinf(i)>0.0d0) then
                                    uy_order1(i) = log(uy_error1(i-1)/uy_error1(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    uy_order2(i) = log(uy_error2(i-1)/uy_error2(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    uy_orderinf(i) = log(uy_errorinf(i-1)/uy_errorinf(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              else
                                    uy_order1(i) = -1.0d0
                                    uy_order2(i) = -1.0d0
                                    uy_orderinf(i) = -1.0d0
                              end if
                              if(uz_errorinf(i)>0.0d0) then
                                    uz_order1(i) = log(uz_error1(i-1)/uz_error1(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    uz_order2(i) = log(uz_error2(i-1)/uz_error2(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    uz_orderinf(i) = log(uz_errorinf(i-1)/uz_errorinf(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              else
                                    uz_order1(i) = -1.0d0
                                    uz_order2(i) = -1.0d0
                                    uz_orderinf(i) = -1.0d0
                              end if
                              if(u_errorinf(i)>0.0d0) then
                                    u_order1(i) = log(u_error1(i-1)/u_error1(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    u_order2(i) = log(u_error2(i-1)/u_error2(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    u_orderinf(i) = log(u_errorinf(i-1)/u_errorinf(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              else
                                    u_order1(i) = -1.0d0
                                    u_order2(i) = -1.0d0
                                    u_orderinf(i) = -1.0d0
                              end if
                        else
                              ux_order1(i) = -1.0d0
                              ux_order2(i) = -1.0d0
                              ux_orderinf(i) = -1.0d0
                              uy_order1(i) = -1.0d0
                              uy_order2(i) = -1.0d0
                              uy_orderinf(i) = -1.0d0
                              uz_order1(i) = -1.0d0
                              uz_order2(i) = -1.0d0
                              uz_orderinf(i) = -1.0d0
                              u_order1(i) = -1.0d0
                              u_order2(i) = -1.0d0
                              u_orderinf(i) = -1.0d0
                        end if
                  end if
            else
                  ux_order1(i) = -1.0d0
                  ux_order2(i) = -1.0d0
                  ux_orderinf(i) = -1.0d0
                  uy_order1(i) = -1.0d0
                  uy_order2(i) = -1.0d0
                  uy_orderinf(i) = -1.0d0
                  uz_order1(i) = -1.0d0
                  uz_order2(i) = -1.0d0
                  uz_orderinf(i) = -1.0d0
                  u_order1(i) = -1.0d0
                  u_order2(i) = -1.0d0
                  u_orderinf(i) = -1.0d0
            end if

            ! state message
            call fvl_loginfo("  >>> Writing approximate velocity errors")

            ! suffixes
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_error.fvd3"
            else
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_error.fvd3.bin"
            end if

            ! write file
            solutionsfile = fvl_data_file(path=fvl_trim(u_suffix),form=errorsfilesform(i),action="write")
            call solutionsfile % open()
            call solutionsfile % writedata(u_errorvals,u_map,fvl_trim(u_label)//" Error",fvl_data_file_cellmeanvalues,times(i),steps(i))
            call solutionsfile % close()

      else

            ux_error1(i) = -1.0d0
            ux_error2(i) = -1.0d0
            ux_errorinf(i) = -1.0d0
            ux_order1(i) = -1.0d0
            ux_order2(i) = -1.0d0
            ux_orderinf(i) = -1.0d0
            uy_error1(i) = -1.0d0
            uy_error2(i) = -1.0d0
            uy_errorinf(i) = -1.0d0
            uy_order1(i) = -1.0d0
            uy_order2(i) = -1.0d0
            uy_orderinf(i) = -1.0d0
            uz_error1(i) = -1.0d0
            uz_error2(i) = -1.0d0
            uz_errorinf(i) = -1.0d0
            uz_order1(i) = -1.0d0
            uz_order2(i) = -1.0d0
            uz_orderinf(i) = -1.0d0
            u_error1(i) = -1.0d0
            u_error2(i) = -1.0d0
            u_errorinf(i) = -1.0d0
            u_order1(i) = -1.0d0
            u_order2(i) = -1.0d0
            u_orderinf(i) = -1.0d0

      end if

      ! continue
      9000 continue

      ! ----------------------------------------------------------------------------
      ! approximate diamond velocity
      ! ----------------------------------------------------------------------------

      ! suffixes
      if(fvl_trim(u_label)/=fvl_char_empty) then
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_diamond.fvd3"
            else
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_diamond.fvd3.bin"
            end if
      else
            go to 10000
      end if

      ! inspect files
      inquire(file=fvl_trim(u_suffix),exist=exist1)

      ! velocity file exists
      if(exist1) then

            ! state message
            call fvl_loginfo("  >>> Reading approximate diamond velocity solution")

            ! read file
            solutionsfile = fvl_data_file(path=fvl_trim(u_suffix),form=solutionsfilesform(i),action="read")
            call solutionsfile % open()
            element = 0
            do while(element/=fvl_data_file_cellmeanvalues)
                  call solutionsfile % readdata(res=diamond_u_approxvals,map=diamond_u_map,element=element,time=times(i),step=steps(i))
            end do
            call solutionsfile % close()

            ! check if data exists
            if(.not. allocated(diamond_u_approxvals)) then
                  call fvl_logerror("Approximate diamond velocity solution not found")
            end if

            ! state message
            call fvl_loginfo("  >>> Computing approximate diamond velocity errors")

            ! allocate memory
            call fvl_allocate(diamond_u_approxvals,2,diamond_numcells(i))
            call fvl_omp_clean(diamond_u_errorvals)

            ! compute errors
            do j=1,size(diamond_u_map,1)
                  diamond_u_errorvals(:,diamond_u_map(j)) = diamond_u_approxvals(:,j)-diamond_u_exactvals(:,diamond_u_map(j))
            end do

            ! compute error norms
            diamond_ux_error1(i) = 0.0d0
            diamond_ux_error2(i) = 0.0d0
            diamond_ux_errorinf(i) = 0.0d0
            diamond_uy_error1(i) = 0.0d0
            diamond_uy_error2(i) = 0.0d0
            diamond_uy_errorinf(i) = 0.0d0
            diamond_uz_error1(i) = 0.0d0
            diamond_uz_error2(i) = 0.0d0
            diamond_uz_errorinf(i) = 0.0d0
            diamond_u_error1(i) = 0.0d0
            diamond_u_error2(i) = 0.0d0
            diamond_u_errorinf(i) = 0.0d0
            do j=1,diamond_numcells(i)
                  cellvolume = diamondmesh%getcellvolume(j)
                  diamond_ux_error1(i) = diamond_ux_error1(i)+abs(diamond_u_errorvals(1,j))*cellvolume
                  diamond_ux_error2(i) = diamond_ux_error2(i)+(diamond_u_errorvals(1,j)**2.0)*cellvolume
                  diamond_ux_errorinf(i) = max(diamond_ux_errorinf(i),abs(diamond_u_errorvals(1,j)))
                  diamond_uy_error1(i) = diamond_uy_error1(i)+abs(diamond_u_errorvals(2,j))*cellvolume
                  diamond_uy_error2(i) = diamond_uy_error2(i)+(diamond_u_errorvals(2,j)**2.0)*cellvolume
                  diamond_uy_errorinf(i) = max(diamond_uy_errorinf(i),abs(diamond_u_errorvals(2,j)))
                  diamond_uz_error1(i) = diamond_uz_error1(i)+abs(diamond_u_errorvals(3,j))*cellvolume
                  diamond_uz_error2(i) = diamond_uz_error2(i)+(diamond_u_errorvals(3,j)**2.0)*cellvolume
                  diamond_uz_errorinf(i) = max(diamond_uz_errorinf(i),abs(diamond_u_errorvals(3,j)))
                  aux=fvl_vectormag3(diamond_u_errorvals(:,j))
                  diamond_u_error1(i) = diamond_u_error1(i)+aux*cellvolume
                  diamond_u_error2(i) = diamond_u_error2(i)+(aux**2.0)*cellvolume
                  diamond_u_errorinf(i) = max(diamond_u_errorinf(i),aux)
            end do
            diamond_ux_error1(i) = diamond_ux_error1(i)/diamondmeshvolumes(i)
            diamond_ux_error2(i) = sqrt(diamond_ux_error2(i)/diamondmeshvolumes(i))
            diamond_uy_error1(i) = diamond_uy_error1(i)/diamondmeshvolumes(i)
            diamond_uy_error2(i) = sqrt(diamond_uy_error2(i)/diamondmeshvolumes(i))
            diamond_uz_error1(i) = diamond_uz_error1(i)/diamondmeshvolumes(i)
            diamond_uz_error2(i) = sqrt(diamond_uz_error2(i)/diamondmeshvolumes(i))
            diamond_u_error1(i) = diamond_u_error1(i)/diamondmeshvolumes(i)
            diamond_u_error2(i) = sqrt(diamond_u_error2(i)/diamondmeshvolumes(i))

            ! compute convergence orders
            if(i>1) then
                  if(convergence==1) then
                        if(diamond_numcells(i-1)>0 .and. diamond_numcells(i-1)/=diamond_numcells(i)) then
                              if(diamond_ux_errorinf(i)>0.0d0) then
                                    diamond_ux_order1(i) = 3.0d0*log(diamond_ux_error1(i-1)/diamond_ux_error1(i))&
                                          /log(fvl_real(diamond_numcells(i))/fvl_real(diamond_numcells(i-1)))
                                    diamond_ux_order2(i) = 3.0d0*log(diamond_ux_error2(i-1)/diamond_ux_error2(i))&
                                          /log(fvl_real(diamond_numcells(i))/fvl_real(diamond_numcells(i-1)))
                                    diamond_ux_orderinf(i) = 3.0d0*log(diamond_ux_errorinf(i-1)/diamond_ux_errorinf(i))&
                                          /log(fvl_real(diamond_numcells(i))/fvl_real(diamond_numcells(i-1)))
                              else
                                    diamond_ux_order1(i) = -1.0d0
                                    diamond_ux_order2(i) = -1.0d0
                                    diamond_ux_orderinf(i) = -1.0d0
                              end if
                              if(diamond_uy_errorinf(i)>0.0d0) then
                                    diamond_uy_order1(i) = 3.0d0*log(diamond_uy_error1(i-1)/diamond_uy_error1(i))&
                                          /log(fvl_real(diamond_numcells(i))/fvl_real(diamond_numcells(i-1)))
                                    diamond_uy_order2(i) = 3.0d0*log(diamond_uy_error2(i-1)/diamond_uy_error2(i))&
                                          /log(fvl_real(diamond_numcells(i))/fvl_real(diamond_numcells(i-1)))
                                    diamond_uy_orderinf(i) = 3.0d0*log(diamond_uy_errorinf(i-1)/diamond_uy_errorinf(i))&
                                          /log(fvl_real(diamond_numcells(i))/fvl_real(diamond_numcells(i-1)))
                              else
                                    diamond_uy_order1(i) = -1.0d0
                                    diamond_uy_order2(i) = -1.0d0
                                    diamond_uy_orderinf(i) = -1.0d0
                              end if
                              if(diamond_uz_errorinf(i)>0.0d0) then
                                    diamond_uz_order1(i) = 3.0d0*log(diamond_uz_error1(i-1)/diamond_uz_error1(i))&
                                          /log(fvl_real(diamond_numcells(i))/fvl_real(diamond_numcells(i-1)))
                                    diamond_uz_order2(i) = 3.0d0*log(diamond_uz_error2(i-1)/diamond_uz_error2(i))&
                                          /log(fvl_real(diamond_numcells(i))/fvl_real(diamond_numcells(i-1)))
                                    diamond_uz_orderinf(i) = 3.0d0*log(diamond_uz_errorinf(i-1)/diamond_uz_errorinf(i))&
                                          /log(fvl_real(diamond_numcells(i))/fvl_real(diamond_numcells(i-1)))
                              else
                                    diamond_uz_order1(i) = -1.0d0
                                    diamond_uz_order2(i) = -1.0d0
                                    diamond_uz_orderinf(i) = -1.0d0
                              end if
                              if(diamond_u_errorinf(i)>0.0d0) then
                                    diamond_u_order1(i) = 3.0d0*log(diamond_u_error1(i-1)/diamond_u_error1(i))&
                                          /log(fvl_real(diamond_numcells(i))/fvl_real(diamond_numcells(i-1)))
                                    diamond_u_order2(i) = 3.0d0*log(diamond_u_error2(i-1)/diamond_u_error2(i))&
                                          /log(fvl_real(diamond_numcells(i))/fvl_real(diamond_numcells(i-1)))
                                    diamond_u_orderinf(i) = 3.0d0*log(diamond_u_errorinf(i-1)/diamond_u_errorinf(i))&
                                          /log(fvl_real(diamond_numcells(i))/fvl_real(diamond_numcells(i-1)))
                              else
                                    diamond_u_order1(i) = -1.0d0
                                    diamond_u_order2(i) = -1.0d0
                                    diamond_u_orderinf(i) = -1.0d0
                              end if
                        else
                              diamond_ux_order1(i) = -1.0d0
                              diamond_ux_order2(i) = -1.0d0
                              diamond_ux_orderinf(i) = -1.0d0
                              diamond_uy_order1(i) = -1.0d0
                              diamond_uy_order2(i) = -1.0d0
                              diamond_uy_orderinf(i) = -1.0d0
                              diamond_uz_order1(i) = -1.0d0
                              diamond_uz_order2(i) = -1.0d0
                              diamond_uz_orderinf(i) = -1.0d0
                              diamond_u_order1(i) = -1.0d0
                              diamond_u_order2(i) = -1.0d0
                              diamond_u_orderinf(i) = -1.0d0
                        end if
                  else
                        if(steps(i-1)>0 .and. steps(i-1)/=steps(i)) then
                              if(diamond_ux_errorinf(i)>0.0d0) then
                                    diamond_ux_order1(i) = log(diamond_ux_error1(i-1)/diamond_ux_error1(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    diamond_ux_order2(i) = log(diamond_ux_error2(i-1)/diamond_ux_error2(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    diamond_ux_orderinf(i) = log(diamond_ux_errorinf(i-1)/diamond_ux_errorinf(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              else
                                    diamond_ux_order1(i) = -1.0d0
                                    diamond_ux_order2(i) = -1.0d0
                                    diamond_ux_orderinf(i) = -1.0d0
                              end if
                              if(diamond_uy_errorinf(i)>0.0d0) then
                                    diamond_uy_order1(i) = log(diamond_uy_error1(i-1)/diamond_uy_error1(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    diamond_uy_order2(i) = log(diamond_uy_error2(i-1)/diamond_uy_error2(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    diamond_uy_orderinf(i) = log(diamond_uy_errorinf(i-1)/diamond_uy_errorinf(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              else
                                    diamond_uy_order1(i) = -1.0d0
                                    diamond_uy_order2(i) = -1.0d0
                                    diamond_uy_orderinf(i) = -1.0d0
                              end if
                              if(diamond_uz_errorinf(i)>0.0d0) then
                                    diamond_uz_order1(i) = log(diamond_uz_error1(i-1)/diamond_uz_error1(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    diamond_uz_order2(i) = log(diamond_uz_error2(i-1)/diamond_uz_error2(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    diamond_uz_orderinf(i) = log(diamond_uz_errorinf(i-1)/diamond_uz_errorinf(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              else
                                    diamond_uz_order1(i) = -1.0d0
                                    diamond_uz_order2(i) = -1.0d0
                                    diamond_uz_orderinf(i) = -1.0d0
                              end if
                              if(diamond_u_errorinf(i)>0.0d0) then
                                    diamond_u_order1(i) = log(diamond_u_error1(i-1)/diamond_u_error1(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    diamond_u_order2(i) = log(diamond_u_error2(i-1)/diamond_u_error2(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    diamond_u_orderinf(i) = log(diamond_u_errorinf(i-1)/diamond_u_errorinf(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              else
                                    diamond_u_order1(i) = -1.0d0
                                    diamond_u_order2(i) = -1.0d0
                                    diamond_u_orderinf(i) = -1.0d0
                              end if
                        else
                              diamond_ux_order1(i) = -1.0d0
                              diamond_ux_order2(i) = -1.0d0
                              diamond_ux_orderinf(i) = -1.0d0
                              diamond_uy_order1(i) = -1.0d0
                              diamond_uy_order2(i) = -1.0d0
                              diamond_uy_orderinf(i) = -1.0d0
                              diamond_uz_order1(i) = -1.0d0
                              diamond_uz_order2(i) = -1.0d0
                              diamond_uz_orderinf(i) = -1.0d0
                              diamond_u_order1(i) = -1.0d0
                              diamond_u_order2(i) = -1.0d0
                              diamond_u_orderinf(i) = -1.0d0
                        end if
                  end if
            else
                  diamond_ux_order1(i) = -1.0d0
                  diamond_ux_order2(i) = -1.0d0
                  diamond_ux_orderinf(i) = -1.0d0
                  diamond_uy_order1(i) = -1.0d0
                  diamond_uy_order2(i) = -1.0d0
                  diamond_uy_orderinf(i) = -1.0d0
                  diamond_uz_order1(i) = -1.0d0
                  diamond_uz_order2(i) = -1.0d0
                  diamond_uz_orderinf(i) = -1.0d0
                  diamond_u_order1(i) = -1.0d0
                  diamond_u_order2(i) = -1.0d0
                  diamond_u_orderinf(i) = -1.0d0
            end if

            ! state message
            call fvl_loginfo("  >>> Writing approximate diamond velocity errors")

            ! suffixes
            if(fvl_trim(solutionsfilesform(i))=="ascii") then
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_diamond_error.fvd3"
            else
                  u_suffix = fvl_trim(solutionsfilesdir(i))//"/"//fvl_trim(u_label)//"_diamond_error.fvd3.bin"
            end if

            ! write file
            solutionsfile = fvl_data_file(path=fvl_trim(u_suffix),form=errorsfilesform(i),action="write")
            call solutionsfile % open()
            call solutionsfile % writedata(diamond_u_errorvals,diamond_u_map,fvl_trim(u_label)//" Error",fvl_data_file_cellmeanvalues,times(i),steps(i))
            call solutionsfile % close()

      else

            diamond_ux_error1(i) = -1.0d0
            diamond_ux_error2(i) = -1.0d0
            diamond_ux_errorinf(i) = -1.0d0
            diamond_ux_order1(i) = -1.0d0
            diamond_ux_order2(i) = -1.0d0
            diamond_ux_orderinf(i) = -1.0d0
            diamond_uy_error1(i) = -1.0d0
            diamond_uy_error2(i) = -1.0d0
            diamond_uy_errorinf(i) = -1.0d0
            diamond_uy_order1(i) = -1.0d0
            diamond_uy_order2(i) = -1.0d0
            diamond_uy_orderinf(i) = -1.0d0
            diamond_uz_error1(i) = -1.0d0
            diamond_uz_error2(i) = -1.0d0
            diamond_uz_errorinf(i) = -1.0d0
            diamond_uz_order1(i) = -1.0d0
            diamond_uz_order2(i) = -1.0d0
            diamond_uz_orderinf(i) = -1.0d0
            diamond_u_error1(i) = -1.0d0
            diamond_u_error2(i) = -1.0d0
            diamond_u_errorinf(i) = -1.0d0
            diamond_u_order1(i) = -1.0d0
            diamond_u_order2(i) = -1.0d0
            diamond_u_orderinf(i) = -1.0d0

      end if

      ! continue
      10000 continue

      ! ----------------------------------------------------------------------------
      ! reference time parameters
      ! ----------------------------------------------------------------------------

      ! inspect file
      inquire(file=fvl_trim(refparamsfilesdir(i))//"/time.fvd",exist=exist1)

      ! file exists
      if(exist1) then

            ! state message
            call fvl_loginfo("  >>> Reading reference time parameters")

            ! read file
            paramfile = fvl_keys_file(path=fvl_trim(refparamsfilesdir(i))//"/time.fvd",form=refparamsfilesform(i),action="read")
            call paramfile % open()
            call paramfile % readkeys()
            refmeshestotaltime(i) = paramfile%getvalue("meshes_total",0.0d0)
            refmodelstotaltime(i) = paramfile%getvalue("models_total",0.0d0)
            refschemestotaltime(i) = paramfile%getvalue("schemes_total",0.0d0)
            refsolutiontotaltime(i) = paramfile%getvalue("solution_total",0.0d0)
            refpostproctotaltime(i) = paramfile%getvalue("postprocessing_total",0.0d0)
            refexecutiontotaltime(i) = paramfile%getvalue("execution_total",0.0d0)
            call paramfile % close()

      else

            refmeshestotaltime(i) = -1.0d0
            refmodelstotaltime(i) = -1.0d0
            refschemestotaltime(i) = -1.0d0
            refsolutiontotaltime(i) = -1.0d0
            refpostproctotaltime(i) = -1.0d0
            refexecutiontotaltime(i) = -1.0d0

      end if

      ! ----------------------------------------------------------------------------
      ! time parameters
      ! ----------------------------------------------------------------------------

      ! inspect file
      inquire(file=fvl_trim(paramsfilesdir(i))//"/time.fvd",exist=exist1)

      ! file exists
      if(exist1) then

            ! state message
            call fvl_loginfo("  >>> Reading time parameters")

            ! read file
            paramfile = fvl_keys_file(path=fvl_trim(paramsfilesdir(i))//"/time.fvd",form=paramsfilesform(i),action="read")
            call paramfile % open()
            call paramfile % readkeys()
            meshestotaltime(i) = paramfile%getvalue("meshes_total",0.0d0)
            modelstotaltime(i) = paramfile%getvalue("models_total",0.0d0)
            schemestotaltime(i) = paramfile%getvalue("schemes_total",0.0d0)
            solutiontotaltime(i) = paramfile%getvalue("solution_total",0.0d0)
            postproctotaltime(i) = paramfile%getvalue("postprocessing_total",0.0d0)
            executiontotaltime(i) = paramfile%getvalue("execution_total",0.0d0)
            call paramfile % close()

            ! compute time ratios
            if(refmeshestotaltime(i)>0.0d0) then
                  meshestotaltimeratio(i) = meshestotaltime(i)/refmeshestotaltime(i)
            else
                  meshestotaltimeratio(i) = -1.0d0
            end if
            if(refmodelstotaltime(i)>0.0d0) then
                  modelstotaltimeratio(i) = modelstotaltime(i)/refmodelstotaltime(i)
            else
                  modelstotaltimeratio(i) = -1.0d0
            end if
            if(refschemestotaltime(i)>0.0d0) then
                  schemestotaltimeratio(i) = schemestotaltime(i)/refschemestotaltime(i)
            else
                  schemestotaltimeratio(i) = -1.0d0
            end if
            if(refsolutiontotaltime(i)>0.0d0) then
                  solutiontotaltimeratio(i) = solutiontotaltime(i)/refsolutiontotaltime(i)
            else
                  solutiontotaltimeratio(i) = -1.0d0
            end if
            if(refpostproctotaltime(i)>0.0d0) then
                  postproctotaltimeratio(i) = postproctotaltime(i)/refpostproctotaltime(i)
            else
                  postproctotaltimeratio(i) = -1.0d0
            end if
            if(refexecutiontotaltime(i)>0.0d0) then
                  executiontotaltimeratio(i) = executiontotaltime(i)/refexecutiontotaltime(i)
            else
                  executiontotaltimeratio(i) = -1.0d0
            end if

      else

            meshestotaltime(i) = -1.0d0
            modelstotaltime(i) = -1.0d0
            schemestotaltime(i) = -1.0d0
            solutiontotaltime(i) = -1.0d0
            postproctotaltime(i) = -1.0d0
            executiontotaltime(i) = -1.0d0
            meshestotaltimeratio(i) = -1.0d0
            modelstotaltimeratio(i) = -1.0d0
            schemestotaltimeratio(i) = -1.0d0
            solutiontotaltimeratio(i) = -1.0d0
            postproctotaltimeratio(i) = -1.0d0
            executiontotaltimeratio(i) = -1.0d0

      end if

      ! ----------------------------------------------------------------------------
      ! reference memory parameters
      ! ----------------------------------------------------------------------------

      ! inspect file
      inquire(file=fvl_trim(refparamsfilesdir(i))//"/memory.fvd",exist=exist1)

      ! file exists
      if(exist1) then

            ! state message
            call fvl_loginfo("  >>> Reading reference memory parameters")

            ! read file
            paramfile = fvl_keys_file(path=fvl_trim(refparamsfilesdir(i))//"/memory.fvd",form=refparamsfilesform(i),action="read")
            call paramfile % open()
            call paramfile % readkeys()
            reftotalexecutionmemory(i) = paramfile%getvalue("execution_total",0.0d0)*1.0d-6
            call paramfile % close()

      else

            reftotalexecutionmemory(i) = -1.0d0

      end if

      ! ----------------------------------------------------------------------------
      ! memory parameters
      ! ----------------------------------------------------------------------------

      ! inspect file
      inquire(file=fvl_trim(paramsfilesdir(i))//"/memory.fvd",exist=exist1)

      ! file exists
      if(exist1) then

            ! state message
            call fvl_loginfo("  >>> Reading memory parameters")

            ! read file
            paramfile = fvl_keys_file(path=fvl_trim(paramsfilesdir(i))//"/memory.fvd",form=paramsfilesform(i),action="read")
            call paramfile % open()
            call paramfile % readkeys()
            totalexecutionmemory(i) = paramfile%getvalue("execution_total",0.0d0)*1.0d-6
            call paramfile % close()

            ! compute time ratios
            if(reftotalexecutionmemory(i)>0.0d0) then
                  totalexecutionmemoryratio(i) = totalexecutionmemory(i)/reftotalexecutionmemory(i)
            else
                  totalexecutionmemoryratio(i) = -1.0d0
            end if

      else

            totalexecutionmemory(i) = -1.0d0
            totalexecutionmemoryratio(i) = -1.0d0

      end if

      ! ----------------------------------------------------------------------------
      ! solution parameters
      ! ----------------------------------------------------------------------------

      ! inspect file
      inquire(file=fvl_trim(paramsfilesdir(i))//"/linear_solver.fvd",exist=exist1)

      if(exist1) then

            ! state message
            call fvl_loginfo("  >>> Reading solution parameters")

            ! read file
            file = fvl_file(path=fvl_trim(paramsfilesdir(i))//"/linear_solver.fvd",form=paramsfilesform(i),action="read")
            call file % open()
            call file % read(linearsystemiterhistogram)
            call file % read(linearsystemresidhistogram)
            call file % close()

            ! get parameters
            linearsystemnumiters(i) = fvl_ninteger(fvl_real(linearsystemiterhistogram(size(linearsystemiterhistogram,1)))/fvl_real(steps(i)))
            linearsystemresid(i) = linearsystemresidhistogram(size(linearsystemresidhistogram,1))

      else

            linearsystemnumiters(i) = -1
            linearsystemresid(i) = -1.0d0

      end if

      ! ----------------------------------------------------------------------------
      ! fixpoint parameters
      ! ----------------------------------------------------------------------------

      ! inspect file
      inquire(file=fvl_trim(paramsfilesdir(i))//"/fixpoint_solver.fvd",exist=exist1)

      if(exist1) then

            ! state message
            call fvl_loginfo("  >>> Reading fixpoint parameters")

            ! read file
            file = fvl_file(path=fvl_trim(paramsfilesdir(i))//"/fixpoint_solver.fvd",form=paramsfilesform(i),action="read")
            call file % open()
            call file % read(fixpointresidhistogram)
            call file % close()

            ! get parameters
            fixpointnumiters(i) = fvl_ninteger(fvl_real(size(fixpointresidhistogram,1))/fvl_real(steps(i)))
            fixpointresid(i) = fixpointresidhistogram(size(fixpointresidhistogram,1))

      else

            fixpointnumiters(i) = -1
            fixpointresid(i) = -1.0d0

      end if

end do

! ============================================================================
! WRITE DATA
! ============================================================================

! state message
call fvl_loginfo("Writing data...")

! write precision
intformat         = "(i0)"
realformat        = "(es20."//fvl_trim(fvl_char(tablesprecision))//")"

! number of meshes
nummeshes         = 0
do i=1,maxnummeshes
      if(primal_numcells(i)/=0) then
            nummeshes = i
      end if
end do

! dofs
if(convergence==1) then
      primal_numdofs(1:nummeshes) = primal_numcells(1:nummeshes)
      diamond_numdofs(1:nummeshes) = diamond_numcells(1:nummeshes)
      primal_doflabel = "N^{\text{P}}_{\text{C}}"
      diamond_doflabel = "N^{\text{D}}_{\text{C}}"
else
      primal_numdofs(1:nummeshes) = steps(1:nummeshes)
      diamond_numdofs(1:nummeshes) = steps(1:nummeshes)
      primal_doflabel = "N^{\text{P}}_{\text{S}}"
      diamond_doflabel = "N^{\text{D}}_{\text{S}}"
end if

! ----------------------------------------------------------------------------
! errors table
! ----------------------------------------------------------------------------

! write errors table
file              = fvl_file(path=fvl_trim(tablesfilesdir)//"errors_table.tex",form="ascii",action="write")
call file         % open()
unit              = file%getunit()

! temperature
if(fvl_trim(t_label)/=fvl_char_empty) then
      write(unit,"(a)") "\begin{table}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Temperature errors and convergence orders.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tabular}{@{}rrrrrrrrrr@{}}"
      write(unit,"(a)") "\toprule"
      write(unit,"(a)") "$"//fvl_trim(primal_doflabel)//"$ & & $E_{1}$ & $O_{1}$ & $E_{2}$ & $O_{2}$ & $E_{\infty}$ & $O_{\infty}$ \\"
      write(unit,"(a)") "\midrule"
      do i=1,nummeshes
            if(primal_numdofs(i)==0) then
                  write(unit,"(a)",advance="no") "--- & & "
                  write(unit,"(a)",advance="no") "--- & --- & "
                  ! write(unit,"(a)",advance="no") "--- & --- & "
                  write(unit,"(a)") "--- & --- \\"
            else
                  if(i==1) then
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(primal_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(t_error1(i),tablesprecision,0.0d0))//" & --- & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(t_error2(i),tablesprecision,0.0d0))//" & --- & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(t_errorinf(i),tablesprecision,0.0d0))//" & --- \\"
                  else
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(primal_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(t_error1(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(t_order1(i),tablesprecision,0.0d0))//" & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(t_error2(i),tablesprecision,0.0d0))//" & "&
                        !      //fvl_trim(fvl_latexchar(t_order2(i),tablesprecision,0.0d0))//" & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(t_errorinf(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(t_orderinf(i),tablesprecision,0.0d0))//" \\"
                  end if
            end if
      end do
      write(unit,"(a)") "\bottomrule"
      write(unit,"(a)") "\end{tabular}"
      write(unit,"(a)") "\end{table}"
      write(unit,"(a)") " "
end if

! pressure
if(fvl_trim(p_label)/=fvl_char_empty) then
      write(unit,"(a)") "\begin{table}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Pressure errors and convergence orders.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tabular}{@{}rrrrrrrrrr@{}}"
      write(unit,"(a)") "\toprule"
      write(unit,"(a)") "$"//fvl_trim(primal_doflabel)//"$ & & $E_{1}$ & $O_{1}$ & $E_{2}$ & $O_{2}$ & $E_{\infty}$ & $O_{\infty}$ \\"
      write(unit,"(a)") "\midrule"
      do i=1,nummeshes
            if(primal_numdofs(i)==0) then
                  write(unit,"(a)",advance="no") "--- & & "
                  write(unit,"(a)",advance="no") "--- & --- & "
                  ! write(unit,"(a)",advance="no") "--- & --- & "
                  write(unit,"(a)") "--- & --- \\"
            else
                  if(i==1) then
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(primal_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(p_error1(i),tablesprecision,0.0d0))//" & --- & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(p_error2(i),tablesprecision,0.0d0))//" & --- & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(p_errorinf(i),tablesprecision,0.0d0))//" & --- \\"
                  else
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(primal_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(p_error1(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(p_order1(i),tablesprecision,0.0d0))//" & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(p_error2(i),tablesprecision,0.0d0))//" & "&
                        !      //fvl_trim(fvl_latexchar(p_order2(i),tablesprecision,0.0d0))//" & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(p_errorinf(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(p_orderinf(i),tablesprecision,0.0d0))//" \\"
                  end if
            end if
      end do
      write(unit,"(a)") "\bottomrule"
      write(unit,"(a)") "\end{tabular}"
      write(unit,"(a)") "\end{table}"
      write(unit,"(a)") " "
end if

! velocity
if(fvl_trim(u_label)/=fvl_char_empty) then
      write(unit,"(a)") "\begin{table}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Velocity x-component errors and convergence orders.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tabular}{@{}rrrrrrrrrr@{}}"
      write(unit,"(a)") "\toprule"
      write(unit,"(a)") "$"//fvl_trim(primal_doflabel)//"$ & & $E_{1}$ & $O_{1}$ & $E_{2}$ & $O_{2}$ & $E_{\infty}$ & $O_{\infty}$ \\"
      write(unit,"(a)") "\midrule"
      do i=1,nummeshes
            if(primal_numdofs(i)==0) then
                  write(unit,"(a)",advance="no") "--- & & "
                  write(unit,"(a)",advance="no") "--- & --- & "
                  ! write(unit,"(a)",advance="no") "--- & --- & "
                  write(unit,"(a)") "--- & --- \\"
            else
                  if(i==1) then
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(primal_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(ux_error1(i),tablesprecision,0.0d0))//" & --- & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(ux_error2(i),tablesprecision,0.0d0))//" & --- & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(ux_errorinf(i),tablesprecision,0.0d0))//" & --- \\"
                  else
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(primal_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(ux_error1(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(ux_order1(i),tablesprecision,0.0d0))//" & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(ux_error2(i),tablesprecision,0.0d0))//" & "&
                        !       //fvl_trim(fvl_latexchar(ux_order2(i),tablesprecision,0.0d0))//" & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(ux_errorinf(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(ux_orderinf(i),tablesprecision,0.0d0))//" \\"
                  end if
            end if
      end do
      write(unit,"(a)") "\bottomrule"
      write(unit,"(a)") "\end{tabular}"
      write(unit,"(a)") "\end{table}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{table}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Velocity y-component errors and convergence orders.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tabular}{@{}rrrrrrrrrr@{}}"
      write(unit,"(a)") "\toprule"
      write(unit,"(a)") "$"//fvl_trim(primal_doflabel)//"$ & & $E_{1}$ & $O_{1}$ & $E_{2}$ & $O_{2}$ & $E_{\infty}$ & $O_{\infty}$ \\"
      write(unit,"(a)") "\midrule"
      do i=1,nummeshes
            if(primal_numdofs(i)==0) then
                  write(unit,"(a)",advance="no") "--- & & "
                  write(unit,"(a)",advance="no") "--- & --- & "
                  ! write(unit,"(a)",advance="no") "--- & --- & "
                  write(unit,"(a)") "--- & --- \\"
            else
                  if(i==1) then
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(primal_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(uy_error1(i),tablesprecision,0.0d0))//" & --- & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(uy_error2(i),tablesprecision,0.0d0))//" & --- & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(uy_errorinf(i),tablesprecision,0.0d0))//" & --- \\"
                  else
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(primal_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(uy_error1(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(uy_order1(i),tablesprecision,0.0d0))//" & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(uy_error2(i),tablesprecision,0.0d0))//" & "&
                        !       //fvl_trim(fvl_latexchar(uy_order2(i),tablesprecision,0.0d0))//" & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(uy_errorinf(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(uy_orderinf(i),tablesprecision,0.0d0))//" \\"
                  end if
            end if
      end do
      write(unit,"(a)") "\bottomrule"
      write(unit,"(a)") "\end{tabular}"
      write(unit,"(a)") "\end{table}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{table}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Velocity errors and convergence orders.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tabular}{@{}rrrrrrrrrr@{}}"
      write(unit,"(a)") "\toprule"
      write(unit,"(a)") "$"//fvl_trim(primal_doflabel)//"$ & & $E_{1}$ & $O_{1}$ & $E_{2}$ & $O_{2}$ & $E_{\infty}$ & $O_{\infty}$ \\"
      write(unit,"(a)") "\midrule"
      do i=1,nummeshes
            if(primal_numdofs(i)==0) then
                  write(unit,"(a)",advance="no") "--- & & "
                  write(unit,"(a)",advance="no") "--- & --- & "
                  ! write(unit,"(a)",advance="no") "--- & --- & "
                  write(unit,"(a)") "--- & --- \\"
            else
                  if(i==1) then
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(primal_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(u_error1(i),tablesprecision,0.0d0))//" & --- & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(u_error2(i),tablesprecision,0.0d0))//" & --- & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(u_errorinf(i),tablesprecision,0.0d0))//" & --- \\"
                  else
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(primal_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(u_error1(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(u_order1(i),tablesprecision,0.0d0))//" & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(u_error2(i),tablesprecision,0.0d0))//" & "&
                        !       //fvl_trim(fvl_latexchar(u_order2(i),tablesprecision,0.0d0))//" & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(u_errorinf(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(u_orderinf(i),tablesprecision,0.0d0))//" \\"
                  end if
            end if
      end do
      write(unit,"(a)") "\bottomrule"
      write(unit,"(a)") "\end{tabular}"
      write(unit,"(a)") "\end{table}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{table}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Diamond velocity x-component errors and convergence orders.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tabular}{@{}rrrrrrrrrr@{}}"
      write(unit,"(a)") "\toprule"
      write(unit,"(a)") "$"//fvl_trim(diamond_doflabel)//"$ & & $E_{1}$ & $O_{1}$ & $E_{2}$ & $O_{2}$ & $E_{\infty}$ & $O_{\infty}$ \\"
      write(unit,"(a)") "\midrule"
      do i=1,nummeshes
            if(diamond_numdofs(i)==0) then
                  write(unit,"(a)",advance="no") "--- & & "
                  write(unit,"(a)",advance="no") "--- & --- & "
                  ! write(unit,"(a)",advance="no") "--- & --- & "
                  write(unit,"(a)") "--- & --- \\"
            else
                  if(i==1) then
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(diamond_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(diamond_ux_error1(i),tablesprecision,0.0d0))//" & --- & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(diamond_ux_error2(i),tablesprecision,0.0d0))//" & --- & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(diamond_ux_errorinf(i),tablesprecision,0.0d0))//" & --- \\"
                  else
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(diamond_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(diamond_ux_error1(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(diamond_ux_order1(i),tablesprecision,0.0d0))//" & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(diamond_ux_error2(i),tablesprecision,0.0d0))//" & "&
                        !       //fvl_trim(fvl_latexchar(diamond_ux_order2(i),tablesprecision,0.0d0))//" & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(diamond_ux_errorinf(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(diamond_ux_orderinf(i),tablesprecision,0.0d0))//" \\"
                  end if
            end if
      end do
      write(unit,"(a)") "\bottomrule"
      write(unit,"(a)") "\end{tabular}"
      write(unit,"(a)") "\end{table}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{table}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Diamond velocity y-component errors and convergence orders.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tabular}{@{}rrrrrrrrrr@{}}"
      write(unit,"(a)") "\toprule"
      write(unit,"(a)") "$"//fvl_trim(diamond_doflabel)//"$ & & $E_{1}$ & $O_{1}$ & $E_{2}$ & $O_{2}$ & $E_{\infty}$ & $O_{\infty}$ \\"
      write(unit,"(a)") "\midrule"
      do i=1,nummeshes
            if(diamond_numdofs(i)==0) then
                  write(unit,"(a)",advance="no") "--- & & "
                  write(unit,"(a)",advance="no") "--- & --- & "
                  ! write(unit,"(a)",advance="no") "--- & --- & "
                  write(unit,"(a)") "--- & --- \\"
            else
                  if(i==1) then
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(diamond_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(diamond_uy_error1(i),tablesprecision,0.0d0))//" & --- & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(diamond_uy_error2(i),tablesprecision,0.0d0))//" & --- & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(diamond_uy_errorinf(i),tablesprecision,0.0d0))//" & --- \\"
                  else
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(diamond_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(diamond_uy_error1(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(diamond_uy_order1(i),tablesprecision,0.0d0))//" & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(diamond_uy_error2(i),tablesprecision,0.0d0))//" & "&
                        !       //fvl_trim(fvl_latexchar(diamond_uy_order2(i),tablesprecision,0.0d0))//" & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(diamond_uy_errorinf(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(diamond_uy_orderinf(i),tablesprecision,0.0d0))//" \\"
                  end if
            end if
      end do
      write(unit,"(a)") "\bottomrule"
      write(unit,"(a)") "\end{tabular}"
      write(unit,"(a)") "\end{table}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{table}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Diamond velocity errors and convergence orders.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tabular}{@{}rrrrrrrrrr@{}}"
      write(unit,"(a)") "\toprule"
      write(unit,"(a)") "$"//fvl_trim(diamond_doflabel)//"$ & & $E_{1}$ & $O_{1}$ & $E_{2}$ & $O_{2}$ & $E_{\infty}$ & $O_{\infty}$ \\"
      write(unit,"(a)") "\midrule"
      do i=1,nummeshes
            if(diamond_numdofs(i)==0) then
                  write(unit,"(a)",advance="no") "--- & & "
                  write(unit,"(a)",advance="no") "--- & --- & "
                  ! write(unit,"(a)",advance="no") "--- & --- & "
                  write(unit,"(a)") "--- & --- \\"
            else
                  if(i==1) then
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(diamond_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(diamond_u_error1(i),tablesprecision,0.0d0))//" & --- & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(diamond_u_error2(i),tablesprecision,0.0d0))//" & --- & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(diamond_u_errorinf(i),tablesprecision,0.0d0))//" & --- \\"
                  else
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(diamond_numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(diamond_u_error1(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(diamond_u_order1(i),tablesprecision,0.0d0))//" & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(diamond_u_error2(i),tablesprecision,0.0d0))//" & "&
                        !       //fvl_trim(fvl_latexchar(diamond_u_order2(i),tablesprecision,0.0d0))//" & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(diamond_u_errorinf(i),tablesprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(diamond_u_orderinf(i),tablesprecision,0.0d0))//" \\"
                  end if
            end if
      end do
      write(unit,"(a)") "\bottomrule"
      write(unit,"(a)") "\end{tabular}"
      write(unit,"(a)") "\end{table}"
      write(unit,"(a)") " "
end if

! close file
write(unit,"(a)") "% end of file"
call file          % close()

! ----------------------------------------------------------------------------
! solution table
! ----------------------------------------------------------------------------

! write solution table
file              = fvl_file(path=fvl_trim(tablesfilesdir)//"solution_table.tex",form="ascii",action="write")
call file         % open()
unit              = file%getunit()
write(unit,"(a)") "\begin{table}[!h]"
write(unit,"(a)") "\centering"
write(unit,"(a)") "\caption{Iterations and residuals.}"
write(unit,"(a)") "\label{}"
write(unit,"(a)") "\begin{tabular}{@{}rrrrrr@{}}"
write(unit,"(a)") "\toprule"
write(unit,"(a)") "$"//fvl_trim(primal_doflabel)//"$ & & $N_{\text{fixpoint}}$ & $R_{\text{fixpoint}}$ & $N_{\text{linearsystem}}$ & $R_{\text{linearsystem}}$ \\"
write(unit,"(a)") "\midrule"
do i=1,nummeshes
      if(primal_numdofs(i)==0) then
            write(unit,"(a)",advance="no") "--- & & "
            write(unit,"(a)",advance="no") "--- & --- & "
            write(unit,"(a)") "--- & --- \\"
      else
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(primal_numdofs(i)))//" & & "
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(fixpointnumiters(i)))//" & "&
                  //fvl_trim(fvl_latexchar(fixpointresid(i),tablesprecision,0.0d0))//" & "
            write(unit,"(a)") fvl_trim(fvl_latexchar(linearsystemnumiters(i)))//" & "&
                  //fvl_trim(fvl_latexchar(linearsystemresid(i),tablesprecision,0.0d0))//" \\"
      end if
end do
write(unit,"(a)") "\bottomrule"
write(unit,"(a)") "\end{tabular}"
write(unit,"(a)") "\end{table}"
write(unit,"(a)") " "
write(unit,"(a)") "% end of file"
call file         % close()

! ----------------------------------------------------------------------------
! time table
! ----------------------------------------------------------------------------

! write time table
file              = fvl_file(path=fvl_trim(tablesfilesdir)//"time_table.tex",form="ascii",action="write")
call file         % open()
unit              = file%getunit()
write(unit,"(a)") "\begin{table}[!h]"
write(unit,"(a)") "\centering"
write(unit,"(a)") "\caption{Time.}"
write(unit,"(a)") "\label{}"
write(unit,"(a)") "\begin{tabular}{@{}rrrrrrrrrrrrrrr@{}}"
write(unit,"(a)") "\toprule"
write(unit,"(a)",advance="no") "$"//fvl_trim(primal_doflabel)//"$ & & "
write(unit,"(a)",advance="no") "$T_{\text{meshes}}$ & $R_{\text{meshes}}$ & "
write(unit,"(a)",advance="no") "$T_{\text{models}}$ & $R_{\text{models}}$ & "
write(unit,"(a)",advance="no") "$T_{\text{schemes}}$ & $R_{\text{schemes}}$ & "
write(unit,"(a)",advance="no") "$T_{\text{solution}}$ & $R_{\text{solution}}$ & "
write(unit,"(a)",advance="no") "$T_{\text{postproc}}$ & $R_{\text{postproc}}$ & "
write(unit,"(a)") "$T_{\text{execution}}$ & $R_{\text{execution}}$ \\"
write(unit,"(a)") "\midrule"
do i=1,nummeshes
      if(primal_numdofs(i)==0) then
            write(unit,"(a)",advance="no") "--- & & "
            write(unit,"(a)",advance="no") "--- & --- & "
            write(unit,"(a)",advance="no") "--- & --- & "
            write(unit,"(a)",advance="no") "--- & --- & "
            write(unit,"(a)",advance="no") "--- & --- & "
            write(unit,"(a)",advance="no") "--- & --- & "
            write(unit,"(a)") "--- & --- \\"
      else
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(primal_numdofs(i)))//" & & "
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(meshestotaltime(i),tablesprecision,0.0d0))//" & "&
                  //fvl_trim(fvl_latexchar(meshestotaltimeratio(i),tablesprecision,0.0d0))//" & "
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(modelstotaltime(i),tablesprecision,0.0d0))//" & "&
                  //fvl_trim(fvl_latexchar(modelstotaltimeratio(i),tablesprecision,0.0d0))//" & "
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(schemestotaltime(i),tablesprecision,0.0d0))//" & "&
                  //fvl_trim(fvl_latexchar(schemestotaltimeratio(i),tablesprecision,0.0d0))//" & "
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(solutiontotaltime(i),tablesprecision,0.0d0))//" & "&
                  //fvl_trim(fvl_latexchar(solutiontotaltimeratio(i),tablesprecision,0.0d0))//" & "
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(postproctotaltime(i),tablesprecision,0.0d0))//" & "&
                  //fvl_trim(fvl_latexchar(postproctotaltimeratio(i),tablesprecision,0.0d0))//" & "
            write(unit,"(a)") fvl_trim(fvl_latexchar(executiontotaltime(i),tablesprecision,0.0d0))//" & "&
                  //fvl_trim(fvl_latexchar(executiontotaltimeratio(i),tablesprecision,0.0d0))//" \\"
      end if
end do
write(unit,"(a)") "\bottomrule"
write(unit,"(a)") "\end{tabular}"
write(unit,"(a)") "\end{table}"
write(unit,"(a)") " "
write(unit,"(a)") "% end of file"
call file         % close()

! ----------------------------------------------------------------------------
! memory table
! ----------------------------------------------------------------------------

! write iterations and residual table
file              = fvl_file(path=fvl_trim(tablesfilesdir)//"memory_table.tex",form="ascii",action="write")
call file         % open()
unit              = file%getunit()
write(unit,"(a)") "\begin{table}[!h]"
write(unit,"(a)") "\centering"
write(unit,"(a)") "\caption{Memory usage.}"
write(unit,"(a)") "\label{}"
write(unit,"(a)") "\begin{tabular}{@{}rrrr@{}}"
write(unit,"(a)") "\toprule"
write(unit,"(a)") "$"//fvl_trim(primal_doflabel)//"$ & & $M_{\text{execution}}$ & $R_{\text{execution}}$ \\"
write(unit,"(a)") "\midrule"
do i=1,nummeshes
      if(primal_numdofs(i)==0) then
            write(unit,"(a)") "--- & & --- & --- \\"
      else
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(primal_numdofs(i)))//" & & "
            write(unit,"(a)") fvl_trim(fvl_latexchar(totalexecutionmemory(i),tablesprecision,0.0d0))//" & "&
                  //fvl_trim(fvl_latexchar(totalexecutionmemoryratio(i),tablesprecision,0.0d0))//" \\"
      end if
end do
write(unit,"(a)") "\bottomrule"
write(unit,"(a)") "\end{tabular}"
write(unit,"(a)") "\end{table}"
write(unit,"(a)") " "
write(unit,"(a)") "% end of file"
call file         % close()

! ----------------------------------------------------------------------------
! errors plot
! ----------------------------------------------------------------------------

! write errors plot
file              = fvl_file(path=fvl_trim(tablesfilesdir)//"errors_plot.tex",form="ascii",action="write")
call file         % open()
unit              = file%getunit()

! temperature
if(fvl_trim(p_label)/=fvl_char_empty) then
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Temperature errors \textit{versus} mumber of degrees-of-freedom.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$"//fvl_trim(primal_doflabel)//"$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "% \printslope{3}{2,4,6,8}{800}{1E-11}{3500}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(primal_numdofs(i),intformat))//","//fvl_trim(fvl_char(t_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"

      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(primal_numdofs(i),intformat))//","//fvl_trim(fvl_char(t_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"

      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(primal_numdofs(i),intformat))//","//fvl_trim(fvl_char(t_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
end if

! pressure
if(fvl_trim(p_label)/=fvl_char_empty) then
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Pressure errors \textit{versus} mumber of degrees-of-freedom.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$"//fvl_trim(primal_doflabel)//"$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "% \printslope{3}{2,4,6,8}{800}{1E-11}{3500}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(primal_numdofs(i),intformat))//","//fvl_trim(fvl_char(p_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"

      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(primal_numdofs(i),intformat))//","//fvl_trim(fvl_char(p_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"

      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(primal_numdofs(i),intformat))//","//fvl_trim(fvl_char(p_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
end if

! velocity
if(fvl_trim(u_label)/=fvl_char_empty) then
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Velocity x-component errors \textit{versus} mumber of degrees-of-freedom.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$"//fvl_trim(primal_doflabel)//"$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "% \printslope{3}{2,4,6,8}{800}{1E-11}{3500}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(primal_numdofs(i),intformat))//","//fvl_trim(fvl_char(ux_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(primal_numdofs(i),intformat))//","//fvl_trim(fvl_char(ux_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(primal_numdofs(i),intformat))//","//fvl_trim(fvl_char(ux_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Velocity y-component errors \textit{versus} mumber of degrees-of-freedom.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$"//fvl_trim(primal_doflabel)//"$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "% \printslope{3}{2,4,6,8}{800}{1E-11}{3500}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(primal_numdofs(i),intformat))//","//fvl_trim(fvl_char(uy_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(primal_numdofs(i),intformat))//","//fvl_trim(fvl_char(uy_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(primal_numdofs(i),intformat))//","//fvl_trim(fvl_char(uy_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Velocity errors \textit{versus} mumber of degrees-of-freedom.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$"//fvl_trim(primal_doflabel)//"$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "% \printslope{3}{2,4,6,8}{800}{1E-11}{3500}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(primal_numdofs(i),intformat))//","//fvl_trim(fvl_char(u_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(primal_numdofs(i),intformat))//","//fvl_trim(fvl_char(u_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(primal_numdofs(i),intformat))//","//fvl_trim(fvl_char(u_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Diamond velocity x-component errors \textit{versus} mumber of degrees-of-freedom.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$"//fvl_trim(diamond_doflabel)//"$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "% \printslope{3}{2,4,6,8}{800}{1E-11}{3500}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(diamond_numdofs(i),intformat))//","//fvl_trim(fvl_char(diamond_ux_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(diamond_numdofs(i),intformat))//","//fvl_trim(fvl_char(diamond_ux_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(diamond_numdofs(i),intformat))//","//fvl_trim(fvl_char(diamond_ux_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Diamond velocity y-component errors \textit{versus} mumber of degrees-of-freedom.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$"//fvl_trim(diamond_doflabel)//"$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "% \printslope{3}{2,4,6,8}{800}{1E-11}{3500}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(diamond_numdofs(i),intformat))//","//fvl_trim(fvl_char(diamond_uy_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(diamond_numdofs(i),intformat))//","//fvl_trim(fvl_char(diamond_uy_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(diamond_numdofs(i),intformat))//","//fvl_trim(fvl_char(diamond_uy_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Diamond velocity errors \textit{versus} mumber of degrees-of-freedom.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$"//fvl_trim(diamond_doflabel)//"$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "% \printslope{3}{2,4,6,8}{800}{1E-11}{3500}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(diamond_numdofs(i),intformat))//","//fvl_trim(fvl_char(diamond_u_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(diamond_numdofs(i),intformat))//","//fvl_trim(fvl_char(diamond_u_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(diamond_numdofs(i),intformat))//","//fvl_trim(fvl_char(diamond_u_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
end if

! close file
write(unit,"(a)") "% end of file"
call file         % close()

! ----------------------------------------------------------------------------
! time plot
! ----------------------------------------------------------------------------

! write time plot
file              = fvl_file(path=fvl_trim(tablesfilesdir)//"time_plot.tex",form="ascii",action="write")
call file         % open()
unit              = file%getunit()

! temperature
if(fvl_trim(t_label)/=fvl_char_empty) then
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Temperature errors \textit{versus} total execution time.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$T_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(t_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(t_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(t_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
end if

! pressure
if(fvl_trim(p_label)/=fvl_char_empty) then
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Pressure errors \textit{versus} total execution time.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$T_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(p_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(p_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(p_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
end if

! velocity
if(fvl_trim(u_label)/=fvl_char_empty) then
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Velocity x-component errors \textit{versus} total execution time.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$T_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(ux_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(ux_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(ux_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Velocity y-component errors \textit{versus} total execution time.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$T_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(uy_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(uy_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(uy_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Velocity errors \textit{versus} total execution time.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$T_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(u_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(u_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(u_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Diamond velocity x-component errors \textit{versus} total execution time.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$T_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(diamond_ux_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(diamond_ux_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(diamond_ux_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Diamond velocity y-component errors \textit{versus} total execution time.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$T_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(diamond_uy_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(diamond_uy_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(diamond_uy_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Diamond velocity errors \textit{versus} total execution time.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$T_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(diamond_u_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(diamond_u_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(diamond_u_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
end if

! close file
write(unit,"(a)") "% end of file"
call file         % close()

! ----------------------------------------------------------------------------
! memory plot
! ----------------------------------------------------------------------------

! write memory plot
file              = fvl_file(path=fvl_trim(tablesfilesdir)//"memory_plot.tex",form="ascii",action="write")
call file         % open()
unit              = file%getunit()

! temperature
if(fvl_trim(t_label)/=fvl_char_empty) then
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Temperature errors \textit{versus} total execution memory.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$M_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(t_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(t_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(t_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
end if

! pressure
if(fvl_trim(p_label)/=fvl_char_empty) then
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Pressure errors \textit{versus} total execution memory.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$M_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(p_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(p_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(p_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
end if

! velocity
if(fvl_trim(u_label)/=fvl_char_empty) then
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Velocity x-component errors \textit{versus} total execution memory.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$M_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(ux_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(ux_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(ux_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Velocity y-component errors \textit{versus} total execution memory.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$M_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(uy_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(uy_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(uy_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Velocity errors \textit{versus} total execution memory.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$M_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(u_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(u_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(u_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Diamond velocity x-component errors \textit{versus} total execution memory.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$M_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(diamond_ux_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(diamond_ux_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(diamond_ux_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Diamond velocity y-component errors \textit{versus} total execution memory.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$M_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(diamond_uy_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(diamond_uy_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(diamond_uy_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Diamond velocity errors \textit{versus} total execution memory.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tikzpicture}[fixed point arithmetic]"
      write(unit,"(a)") "\begin{loglogaxis}["
      write(unit,"(a)") "% xmin=1,"
      write(unit,"(a)") "% xmax=1E7,"
      write(unit,"(a)") "% ymin=1E-12,"
      write(unit,"(a)") "% ymax=1E2,"
      write(unit,"(a)") "width={280},"
      write(unit,"(a)") "height={210},"
      write(unit,"(a)") "grid={major},"
      write(unit,"(a)") "major grid style={dashed},"
      write(unit,"(a)") "xlabel={$M_{\text{execution}}$},"
      write(unit,"(a)") "ylabel={$E$},"
      write(unit,"(a)") "xtick={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "ytick={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "xticklabels={1E-2,1E-1,1E0,1E1,1E2,1E3,1E4,1E5,1E6,1E7},"
      write(unit,"(a)") "yticklabels={1E2,1E0,1E-2,1E-4,1E-6,1E-8,1E-10,1E-12},"
      write(unit,"(a)") "every axis/.append style={font=\normalsize},"
      write(unit,"(a)") "legend style={at={(1.02,1)},anchor=north west,xshift=4pt},"
      write(unit,"(a)") "legend cell align=left,"
      write(unit,"(a)") "legend image post style={scale=1.0},"
      write(unit,"(a)") "]"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=o,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(diamond_u_error1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(diamond_u_error2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(diamond_u_errorinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
end if

! close file
write(unit,"(a)") "% end of file"
call file         % close()

! ----------------------------------------------------------------------------
! report
! ----------------------------------------------------------------------------

! write report
file              = fvl_file(path=fvl_trim(tablesfilesdir)//"report.tex",form="ascii",action="write")
call file         % open()
unit              = file%getunit()
write(unit,"(a)") "\documentclass[8pt,a4paper,landscape]{extreport}"
write(unit,"(a)") "\usepackage[utf8]{inputenc}"
write(unit,"(a)") "\usepackage{amsmath}"
write(unit,"(a)") "\usepackage{amsfonts}"
write(unit,"(a)") "\usepackage{amssymb}"
write(unit,"(a)") "\usepackage{booktabs}"
write(unit,"(a)") "\usepackage{pgfplots}"
write(unit,"(a)") "\usepackage{tikz}"
write(unit,"(a)") "\usetikzlibrary{fixedpointarithmetic}"
write(unit,"(a)") "\usepackage[nomessages]{fp}"
write(unit,"(a)") "\usepackage[left=2cm,right=2cm,top=2cm,bottom=2cm]{geometry}"
write(unit,"(a)") "\setlength{\tabcolsep}{5pt}"
write(unit,"(a)") "\begin{document}"
write(unit,"(a)") "\input{"//fvl_trim(fvl_basename(tablesfilesdir))//"errors_table.tex}"
write(unit,"(a)") "\input{"//fvl_trim(fvl_basename(tablesfilesdir))//"solution_table.tex}"
write(unit,"(a)") "\input{"//fvl_trim(fvl_basename(tablesfilesdir))//"time_table.tex}"
write(unit,"(a)") "\input{"//fvl_trim(fvl_basename(tablesfilesdir))//"memory_table.tex}"
write(unit,"(a)") "\input{"//fvl_trim(fvl_basename(tablesfilesdir))//"errors_plot.tex}"
write(unit,"(a)") "\input{"//fvl_trim(fvl_basename(tablesfilesdir))//"time_plot.tex}"
write(unit,"(a)") "\input{"//fvl_trim(fvl_basename(tablesfilesdir))//"memory_plot.tex}"
write(unit,"(a)") "\end{document}"
write(unit,"(a)") "% end of file"
call file         % close()

! ============================================================================
! STOP LIBRARY
! ============================================================================

call fvl_stop()

contains

! ============================================================================
! INCLUDE FUNCTIONS
! ============================================================================

! user-defined code
#include "functions.f"

! ----------------------------------------------------------------------------
! functions
! ----------------------------------------------------------------------------

subroutine fvl_set_t_exactfun(valuefun)
      procedure(fvl_scalarvaluefun)::valuefun
      t_valuefun=>valuefun
end subroutine fvl_set_t_exactfun

subroutine fvl_set_p_exactfun(valuefun)
      procedure(fvl_scalarvaluefun)::valuefun
      p_valuefun=>valuefun
end subroutine fvl_set_p_exactfun

subroutine fvl_set_u_exactfun(valuefun)
      procedure(fvl_vectorvaluefun)::valuefun
      u_valuefun=>valuefun
end subroutine fvl_set_u_exactfun

end program fvl_post_thermoincompflow3d_main
! end of file
