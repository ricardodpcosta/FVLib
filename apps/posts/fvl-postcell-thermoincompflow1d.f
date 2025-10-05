! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author:  Ricardo Costa           |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0a                    |
!| |_|       \_/   |_|_|_| |  Release: June 2016               |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: 2d non-isothermal incompressible fluid flow post-processor
! Modification: October, 2023

#include "macros.f"

program fvl_post_thermoincompflow2d_main

use fvl_lib2d

implicit none

! ============================================================================
! DECLARE VARIABLES
! ============================================================================

! ----------------------------------------------------------------------------
! variables
! ----------------------------------------------------------------------------

integer(kind=__fvl_integer_kind__),parameter::maxnumtests=10
character(len=__fvl_character_len__)::fieldlabel,tabledirectory,&
      intformat,realformat,doflabel,&
      testdirectories(maxnumtests),referencedirectories(maxnumtests),&
      testmeshfilepath,testmeshfileform,testsolutionfiledirectory,testsolutionsfilesform,&
      testreportfiledirectory,testreportsfilesform,testoutputfiledirectory,testoutputfileform,&
      referencemeshfilepath,referencemeshfileform,referencesolutionfiledirectory,referencesolutionsfilesform,&
      referencereportfiledirectory,referencereportsfilesform
integer(kind=__fvl_integer_kind__)::i,j,nummeshes,convergencetytpe,referencesolution,fixconstant,fieldtype,tableprecision,unit,&
      numdofs(maxnumtests),numcells(maxnumtests),refnumcells(maxnumtests),steps(maxnumtests),&
      linearsystemnumiters(maxnumtests),fixpointnumiters(maxnumtests)
real(kind=__fvl_real_kind__)::time,cellarea,meshareas(maxnumtests),&
      times(maxnumtests),fixpointresid(maxnumtests),linearsystemresid(maxnumtests),&
      meshestotaltime(maxnumtests),modelstotaltime(maxnumtests),schemestotaltime(maxnumtests),solutiontotaltime(maxnumtests),&
      postproctotaltime(maxnumtests),executiontotaltime(maxnumtests),&
      refmeshestotaltime(maxnumtests),refmodelstotaltime(maxnumtests),refschemestotaltime(maxnumtests),refsolutiontotaltime(maxnumtests),&
      refpostproctotaltime(maxnumtests),refexecutiontotaltime(maxnumtests),&
      meshestotaltimeratio(maxnumtests),modelstotaltimeratio(maxnumtests),schemestotaltimeratio(maxnumtests),solutiontotaltimeratio(maxnumtests),&
      postproctotaltimeratio(maxnumtests),executiontotaltimeratio(maxnumtests),&
      totalexecutionmemory(maxnumtests),reftotalexecutionmemory(maxnumtests),totalexecutionmemoryratio(maxnumtests),&
      errormag1(maxnumtests),errormag2(maxnumtests),errormaginf(maxnumtests),&
      ordermag1(maxnumtests),ordermag2(maxnumtests),ordermaginf(maxnumtests),&
      errorx1(maxnumtests),errorx2(maxnumtests),errorxinf(maxnumtests),&
      orderx1(maxnumtests),orderx2(maxnumtests),orderxinf(maxnumtests),&
      errory1(maxnumtests),errory2(maxnumtests),erroryinf(maxnumtests),&
      ordery1(maxnumtests),ordery2(maxnumtests),orderyinf(maxnumtests)
integer(kind=__fvl_integer_kind__),allocatable,dimension(:)::linearsystemiterhistogram
real(kind=__fvl_real_kind__),allocatable,dimension(:)::fixpointresidhistogram,linearsystemresidhistogram
real(kind=__fvl_real_kind__),allocatable,dimension(:)::scalarreferencesolution,scalartestsolution,scalartesterror
real(kind=__fvl_real_kind__),allocatable,dimension(:,:)::vectorreferencesolution,vectortestsolution,vectortesterror
logical(kind=__fvl_logical_kind__)::exist
type(fvl_mesh2d)::mesh,refmesh
type(fvl_file)::file
type(fvl_keys_file)::paramfile
type(fvl_data_file)::solutionsfile
type(fvl_patch2d)::domain_patch
type(fvl_scalarmodel2d)::scalarmodel
type(fvl_vectormodel2d)::vectormodel
type(fvl_scalarfield2d_selector)::scalarfieldselector
type(fvl_scalarboundcond2d_selector)::scalarboundcondselector
type(fvl_vectorfield2d_selector)::vectorfieldselector
type(fvl_vectorboundcond2d_selector)::vectorboundcondselector
! type(fvl_meshlocator2d)::meshlocator
! type(fvl_reconstmap2d)::reconstmap

! user-defined variables
#include "variables.f"

! reports variables
#include "fvl_writereports_hdr.f"

! dictionaries variables
#include "fvl_readdictionaries_hdr.f"

! function lists
fvl_scalarfield2d_functionlist=fvl_scalarfunctionlist2d()
fvl_vectorfield2d_functionlist=fvl_vectorfunctionlist2d()

! scalar model allocators
#undef __fvl_model__
#define __fvl_model__ scalarmodel
#include "fvl_setscalarallocator2d_hdr.f"

! vector model allocators
#undef __fvl_model__
#define __fvl_model__ vectormodel
#include "fvl_setvectorallocator2d_hdr.f"

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
      call fvl_loginfo("About: 2d non-isothermal non-Newtonian incompressible fluid flow post-processor.")
      call fvl_loginfo("Usage: fvl-postcell-thermoincompflow2d [options]")
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
convergencetype               = controldictfile%getvalue("post","convergence_type",1)
referencesolution             = controldictfile%getvalue("post","reference_solution",1)
fixconstant                   = controldictfile%getvalue("post","fix_constant",0)
fieldlabel                    = controldictfile%getvalue("post","field_label","t")
fieldtype                     = controldictfile%getvalue("post","field_type",1)
tabledirectory                = controldictfile%getvalue("post","table_directory","output/")
tableprecision                = controldictfile%getvalue("post","table_precision",2)
do i=1,maxnumtests
      testdirectories(i) = controldictfile%getvalue("post","test_directory_"//fvl_trim(fvl_char(i)))
end do
do i=1,maxnumtests
      referencedirectories(i) = controldictfile%getvalue("post","reference_directory_"//fvl_trim(fvl_char(i)))
end do

! read test case parameters
testmeshfilepath              = controldictfile%getvalue("test_case","mesh_file_path","mesh/mesh.fvm2.bin")
testmeshfileform              = controldictfile%getvalue("test_case","mesh_file_form","binary")
testsolutionfiledirectory     = controldictfile%getvalue("test_case","solution_file_directory","output/")
testsolutionfiledirectory     = fvl_trim(solutionfiledirectory)//fvl_trim(fvl_shortchar(time))
testsolutionsfilesform        = controldictfile%getvalue("test_case","solution_file_form","binary")
testreportfiledirectory       = controldictfile%getvalue("test_case","report_file_directory","output/")
testreportfiledirectory       = fvl_trim(reportfiledirectory)//fvl_trim(fvl_shortchar(time))
testreportsfilesform          = controldictfile%getvalue("test_case","report_file_form","ascii")
testoutputfiledirectory       = controldictfile%getvalue("test_case","output_file_directory","output/")
testoutputfiledirectory       = fvl_trim(outputfiledirectory)//fvl_trim(fvl_shortchar(time))
testoutputfileform            = controldictfile%getvalue("test_case","output_file_form","binary")

! read reference case parameters
refmeshfilepath               = controldictfile%getvalue("reference_case","ref_mesh_file_path","mesh/mesh.fvm2.bin")
refmeshfileform               = controldictfile%getvalue("reference_case","ref_mesh_file_form","binary")
refsolutionfiledirectory      = controldictfile%getvalue("reference_case","ref_solution_file_directory","output/")
refsolutionfiledirectory      = fvl_trim(refsolutionfiledirectory)//fvl_trim(fvl_shortchar(time))
refsolutionfileform           = controldictfile%getvalue("reference_case","ref_solution_file_form","binary")
refreportfiledirectory        = controldictfile%getvalue("reference_case","ref_report_file_directory","output/")
refreportfiledirectory        = fvl_trim(refreportfiledirectory)//fvl_trim(fvl_shortchar(time))
refreportfileform             = controldictfile%getvalue("reference_case","ref_report_file_form","ascii")

! ============================================================================
! LOOP TESTS
! ============================================================================

do i=1,maxnumtests

      ! ============================================================================
      ! INITIALIZE MESHES
      ! ============================================================================

      ! state message
      call fvl_loginfo("Initializing mesh...")

      ! inspect file
      inquire(file=fvl_trim(testdirectories(i))//fvl_trim(testmeshfilepath),exist=exist)

      ! file exists
      if(exist) then

            ! check mesh path
            if(i>1) then
                  if(fvl_trim(testdirectories(i))==fvl_trim(testdirectories(i))) then
                        go to 1000
                  end if
            end if

            ! state message
            call fvl_loginfo("  >>> Reading mesh")

            ! read mesh
            mesh = fvl_mesh2d(filepath=fvl_trim(testdirectories(i))//fvl_trim(meshfilepath),fileform=meshfileform,dictfile=controldictfile,dictlabel="mesh")

            ! write mesh report
            call mesh % report("- mesh")

            ! initialize variables
            numcells(i) = mesh%getnumcells()
            refnumcells(i) = mesh%getnumcells()
            meshareas(i) = fvl_omp_sum(mesh%getcellareas())

      else

            ! initialize variables
            numcells(i) = 0
            errormag1(i) = -1.0d0
            errormag2(i) = -1.0d0
            errormaginf(i) = -1.0d0
            ordermag1(i) = -1.0d0
            ordermag2(i) = -1.0d0
            ordermaginf(i) = -1.0d0
            errorx1(i) = -1.0d0
            errorx2(i) = -1.0d0
            errorxinf(i) = -1.0d0
            orderx1(i) = -1.0d0
            orderx2(i) = -1.0d0
            orderxinf(i) = -1.0d0
            errory1(i) = -1.0d0
            errory2(i) = -1.0d0
            erroryinf(i) = -1.0d0
            ordery1(i) = -1.0d0
            ordery2(i) = -1.0d0
            orderyinf(i) = -1.0d0
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

      1000 continue

      ! ============================================================================
      ! INITIALIZE REFERENCE MESH
      ! ============================================================================

      if(referencesolution==2) then

            ! state message
            call fvl_loginfo("Initializing reference mesh...")

            if(i>1) then
                  if(fvl_trim(refmeshesfilepath(i))==fvl_trim(refmeshesfilepath(i-1))) then
                        2000 continue
                  end if
            end if

            ! state message
            call fvl_loginfo("  >>> Reading reference mesh")

            ! read mesh
            mesh = fvl_mesh2d(filepath=fvl_trim(referencedirectories(i))//fvl_trim(refmeshfilepath),fileform=refmeshfileform,dictfile=controldictfile,dictlabel="mesh")

            ! write mesh report
            call refmesh % report("- reference mesh")

            ! state message
            call fvl_loginfo("  >>> Computing mesh locators")

            ! ! mesh locator
            ! meshlocator = fvl_meshlocator2d(refmesh,postdictfile,"mesh_locator")

            ! ! mesh patches
            ! domain_patch = fvl_patch2d(refmesh)

            ! ! state message
            ! call fvl_loginfo("  >>> Computing reconstruction map")
            !
            ! reconstruction map
            ! reconstmap = fvl_reconstmap2d(domain_patch,refmesh,meshlocator,postdictfile,"reconst_map")

            ! initialize variables
            refnumcells(i) = refmesh%getnumcells()

      end if

      2000 continue

! user-defined code
#include "meshes.f"

      ! ============================================================================
      ! READ REFERENCE SOLUTION
      ! ============================================================================

      ! state message
      call fvl_loginfo("Reading reference solution...")

! user-defined code
#include "models.f"

      if(fieldtype==1) then

            ! scalar model
            scalarmodel       = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_cellmeansfield,&
                                    modeldictfile=modelsdictfile,modeldictlabel"reference_model",&
                                    patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                                    fieldselector=scalarfieldselector,boundcondselector=scalarboundcondselector,&
                                    label=fieldlabel//"_post",filepath=refsolutionfiledirectory,fileform=refsolutionfileform)

            ! evaluate reference solution
            call fvl_allocate(scalarreferencesolution,refnumcells(i))
            call scalarmodel % evalcellmeanvalues(scalarreferencesolution)

            ! state message
            call fvl_loginfo("  >>> Writing reference solution")

            ! write file
            call scalarmodel % write()

      else

            ! vector model
            scalarmodel       = fvl_vectormodel2d(mesh=mesh,fieldtype=fvl_field2d_cellmeansfield,&
                                    modeldictfile=modelsdictfile,modeldictlabel="reference_model",&
                                    patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                                    fieldselector=scalarfieldselector,boundcondselector=scalarboundcondselector,&
                                    label=fieldlabel//"_post",filepath=refsolutionfiledirectory,fileform=refsolutionfileform)

            ! evaluate reference solution
            call fvl_allocate(vectorreferencesolution,2,refnumcells(i))
            call vectormodel % evalcellmeanvalues(vectorreferencesolution)

            ! state message
            call fvl_loginfo("  >>> Writing reference solution")

            ! write file
            call vectormodel % write()

      end if

      if(referencesolution==2) then

            if(i>1) then
                  if(fvl_trim(refsolutionfiledirectory(i))==fvl_trim(refsolutionfiledirectory(i-1))) then
                        go to 3000
                  end if
            end if

            if(fieldtype==1) then

                  ! state message
                  call fvl_loginfo("  >>> Computing reference solution")

                  ! allocate memory
                  call fvl_allocate(scalarreferencesolution,numcells(i))

                  ! ! map reference solution
                  ! call reconstmap % mapcellmeanvalues(mesh,scalarreferencesolution)

            else

                  ! state message
                  call fvl_loginfo("  >>> Computing reference solution")

                  ! allocate memory
                  call fvl_allocate(vectorreferencesolution,numcells(i))

                  ! ! map reference solution
                  ! call reconstmap % mapcellmeanvalues(mesh,vectorreferencesolution)

            end if

      end if

      3000 continue

      ! ============================================================================
      ! READ TEST SOLUTION
      ! ============================================================================

      ! state message
      call fvl_loginfo("Reading test solution...")

      if(fieldtype==1) then

            ! scalar model
            scalarmodel       = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_cellmeansfield,&
                                    modeldictfile=modelsdictfile,modeldictlabel"test_model",&
                                    patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                                    fieldselector=scalarfieldselector,boundcondselector=scalarboundcondselector,&
                                    label=fieldlabel//"_post",filepath=testsolutionfiledirectory,fileform=testsolutionfileform)

            ! allocate memory
            call fvl_allocate(scalartestsolution,numcells(i))

            ! evaluate test solution
            call scalarmodel % evalcellmeanvalues(scalartestsolution)

      else

            ! vector model
            scalarmodel       = fvl_vectormodel2d(mesh=mesh,fieldtype=fvl_field2d_cellmeansfield,&
                                    modeldictfile=modelsdictfile,modeldictlabel="test_model",&
                                    patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                                    fieldselector=scalarfieldselector,boundcondselector=scalarboundcondselector,&
                                    label=fieldlabel//"_post",filepath=testsolutionfiledirectory,fileform=testsolutionfileform)

            ! allocate memory
            call fvl_allocate(vectortestsolution,2,numcells(i))

            ! evaluate test solution
            call vectormodel % evalcellmeanvalues(vectortestsolution)

      end if

      ! ============================================================================
      ! COMPUTE TEST ERRORS
      ! ============================================================================

      ! state message
      call fvl_loginfo("Computing test errors...")

      ! ----------------------------------------------------------------------------
      ! scalar field
      ! ----------------------------------------------------------------------------

      if(fieldtype==1) then

            ! allocate memory
            call fvl_allocate(scalartesterror,numcells(i))
            call fvl_omp_clean(scalartesterror)

            ! fix constant
            if(fixconstant==1) then
                  call fvl_omp_sub(scalarreferencesolution,fvl_omp_dot(scalarreferencesolution,mesh%getcellareas())/meshareas(i))
                  call fvl_omp_sub(scalartestsolution,fvl_omp_dot(scalartestsolution,mesh%getcellareas())/meshareas(i))
            end if

            ! compute errors
            call fvl_omp_sub(scalarreferencesolution,scalartestsolution,scalartesterror)

            ! compute error norms
            errormag1(i) = 0.0d0
            errormag2(i) = 0.0d0
            errormaginf(i) = 0.0d0
            do j=1,numcells(i)
                  cellarea = mesh%getcellarea(j)
                  errormag1(i) = errormag1(i)+abs(scalartesterror(j))*cellarea
                  errormag2(i) = errormag2(i)+(scalartesterror(j)**2.0)*cellarea
                  errormaginf(i) = max(errormaginf(i),abs(scalartesterror(j)))
            end do
            errormag1(i) = errormag1(i)/meshareas(i)
            errormag2(i) = sqrt(errormag2(i)/meshareas(i))

            ! compute convergence orders
            if(i>1) then
                  if(convergence==1) then
                        if(numcells(i-1)>0 .and. numcells(i-1)/=numcells(i) .and. errormaginf(i)>0.0d0) then
                              ordermag1(i) = log(errormag1(i-1)/errormag1(i))&
                                    /log(fvl_real(numcells(i))/fvl_real(numcells(i-1)))
                              ordermag2(i) = log(errormag2(i-1)/errormag2(i))&
                                    /log(fvl_real(numcells(i))/fvl_real(numcells(i-1)))
                              ordermaginf(i) = log(errormaginf(i-1)/errormaginf(i))&
                                    /log(fvl_real(numcells(i))/fvl_real(numcells(i-1)))
                        else
                              ordermag1(i) = -1.0d0
                              ordermag2(i) = -1.0d0
                              ordermaginf(i) = -1.0d0
                        end if
                  else
                        if(steps(i-1)>0 .and. steps(i-1)/=steps(i) .and. errormaginf(i)>0.0d0) then
                              ordermag1(i) = log(errormag1(i-1)/errormag1(i))&
                                    /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              ordermag2(i) = log(errormag2(i-1)/errormag2(i))&
                                    /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              ordermaginf(i) = log(errormaginf(i-1)/errormaginf(i))&
                                    /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                        else
                              ordermag1(i) = -1.0d0
                              ordermag2(i) = -1.0d0
                              ordermaginf(i) = -1.0d0
                        end if
                  end if
            else
                  ordermag1(i) = -1.0d0
                  ordermag2(i) = -1.0d0
                  ordermaginf(i) = -1.0d0
            end if

            ! state message
            call fvl_loginfo("  >>> Writing test errors")

            !!!

      ! ----------------------------------------------------------------------------
      ! vector field
      ! ----------------------------------------------------------------------------

      else

            ! allocate memory
            call fvl_allocate(vectortesterror,numcells(i))
            call fvl_omp_clean(vectortesterror)

            ! compute errors
            call fvl_omp_sub(vectorreferencesolution,vectortestsolution,vectortesterror)

            ! compute error norms
            errorx1(i) = 0.0d0
            errorx2(i) = 0.0d0
            errorxinf(i) = 0.0d0
            errory1(i) = 0.0d0
            errory2(i) = 0.0d0
            erroryinf(i) = 0.0d0
            errormag1(i) = 0.0d0
            errormag2(i) = 0.0d0
            errormaginf(i) = 0.0d0
            do j=1,numcells(i)
                  cellarea = primalmesh%getcellarea(j)
                  errorx1(i) = errorx1(i)+abs(errormagvals(1,j))*cellarea
                  errorx2(i) = errorx2(i)+(errormagvals(1,j)**2.0)*cellarea
                  errorxinf(i) = max(errorxinf(i),abs(errormagvals(1,j)))
                  errory1(i) = errory1(i)+abs(errormagvals(2,j))*cellarea
                  errory2(i) = errory2(i)+(errormagvals(2,j)**2.0)*cellarea
                  erroryinf(i) = max(erroryinf(i),abs(errormagvals(2,j)))
                  aux=fvl_vectormag2(errormagvals(:,j))
                  errormag1(i) = errormag1(i)+aux*cellarea
                  errormag2(i) = errormag2(i)+(aux**2.0)*cellarea
                  errormaginf(i) = max(errormaginf(i),aux)
            end do
            errorx1(i) = errorx1(i)/primalmeshareas(i)
            errorx2(i) = sqrt(errorx2(i)/primalmeshareas(i))
            errory1(i) = errory1(i)/primalmeshareas(i)
            errory2(i) = sqrt(errory2(i)/primalmeshareas(i))
            errormag1(i) = errormag1(i)/primalmeshareas(i)
            errormag2(i) = sqrt(errormag2(i)/primalmeshareas(i))

            ! compute convergence orders
            if(i>1) then
                  if(convergence==1) then
                        if(numcells(i-1)>0 .and. numcells(i-1)/=numcells(i)) then
                              if(errorxinf(i)>0.0d0) then
                                    orderx1(i) = 2.0d0*log(errorx1(i-1)/errorx1(i))&
                                          /log(fvl_real(numcells(i))/fvl_real(numcells(i-1)))
                                    orderx2(i) = 2.0d0*log(errorx2(i-1)/errorx2(i))&
                                          /log(fvl_real(numcells(i))/fvl_real(numcells(i-1)))
                                    orderxinf(i) = 2.0d0*log(errorxinf(i-1)/errorxinf(i))&
                                          /log(fvl_real(numcells(i))/fvl_real(numcells(i-1)))
                              else
                                    orderx1(i) = -1.0d0
                                    orderx2(i) = -1.0d0
                                    orderxinf(i) = -1.0d0
                              end if
                              if(erroryinf(i)>0.0d0) then
                                    ordery1(i) = 2.0d0*log(errory1(i-1)/errory1(i))&
                                          /log(fvl_real(numcells(i))/fvl_real(numcells(i-1)))
                                    ordery2(i) = 2.0d0*log(errory2(i-1)/errory2(i))&
                                          /log(fvl_real(numcells(i))/fvl_real(numcells(i-1)))
                                    orderyinf(i) = 2.0d0*log(erroryinf(i-1)/erroryinf(i))&
                                          /log(fvl_real(numcells(i))/fvl_real(numcells(i-1)))
                              else
                                    ordery1(i) = -1.0d0
                                    ordery2(i) = -1.0d0
                                    orderyinf(i) = -1.0d0
                              end if
                              if(errormaginf(i)>0.0d0) then
                                    ordermag1(i) = 2.0d0*log(errormag1(i-1)/errormag1(i))&
                                          /log(fvl_real(numcells(i))/fvl_real(numcells(i-1)))
                                    ordermag2(i) = 2.0d0*log(errormag2(i-1)/errormag2(i))&
                                          /log(fvl_real(numcells(i))/fvl_real(numcells(i-1)))
                                    ordermaginf(i) = 2.0d0*log(errormaginf(i-1)/errormaginf(i))&
                                          /log(fvl_real(numcells(i))/fvl_real(numcells(i-1)))
                              else
                                    ordermag1(i) = -1.0d0
                                    ordermag2(i) = -1.0d0
                                    ordermaginf(i) = -1.0d0
                              end if
                        else
                              orderx1(i) = -1.0d0
                              orderx2(i) = -1.0d0
                              orderxinf(i) = -1.0d0
                              ordery1(i) = -1.0d0
                              ordery2(i) = -1.0d0
                              orderyinf(i) = -1.0d0
                              ordermag1(i) = -1.0d0
                              ordermag2(i) = -1.0d0
                              ordermaginf(i) = -1.0d0
                        end if
                  else
                        if(steps(i-1)>0 .and. steps(i-1)/=steps(i)) then
                              if(errorxinf(i)>0.0d0) then
                                    orderx1(i) = log(errorx1(i-1)/errorx1(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    orderx2(i) = log(errorx2(i-1)/errorx2(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    orderxinf(i) = log(errorxinf(i-1)/errorxinf(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              else
                                    orderx1(i) = -1.0d0
                                    orderx2(i) = -1.0d0
                                    orderxinf(i) = -1.0d0
                              end if
                              if(erroryinf(i)>0.0d0) then
                                    ordery1(i) = log(errory1(i-1)/errory1(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    ordery2(i) = log(errory2(i-1)/errory2(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    orderyinf(i) = log(erroryinf(i-1)/erroryinf(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              else
                                    ordery1(i) = -1.0d0
                                    ordery2(i) = -1.0d0
                                    orderyinf(i) = -1.0d0
                              end if
                              if(errormaginf(i)>0.0d0) then
                                    ordermag1(i) = log(errormag1(i-1)/errormag1(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    ordermag2(i) = log(errormag2(i-1)/errormag2(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                                    ordermaginf(i) = log(errormaginf(i-1)/errormaginf(i))&
                                          /log(fvl_real(steps(i))/fvl_real(steps(i-1)))
                              else
                                    ordermag1(i) = -1.0d0
                                    ordermag2(i) = -1.0d0
                                    ordermaginf(i) = -1.0d0
                              end if
                        else
                              orderx1(i) = -1.0d0
                              orderx2(i) = -1.0d0
                              orderxinf(i) = -1.0d0
                              ordery1(i) = -1.0d0
                              ordery2(i) = -1.0d0
                              orderyinf(i) = -1.0d0
                              ordermag1(i) = -1.0d0
                              ordermag2(i) = -1.0d0
                              ordermaginf(i) = -1.0d0
                        end if
                  end if
            else
                  orderx1(i) = -1.0d0
                  orderx2(i) = -1.0d0
                  orderxinf(i) = -1.0d0
                  ordery1(i) = -1.0d0
                  ordery2(i) = -1.0d0
                  orderyinf(i) = -1.0d0
                  ordermag1(i) = -1.0d0
                  ordermag2(i) = -1.0d0
                  ordermaginf(i) = -1.0d0
            end if

            ! state message
            call fvl_loginfo("  >>> Writing approximate velocity errors")

            !!!

      end if

      ! ============================================================================
      ! COMPUTE TEST METRICS
      ! ============================================================================

      ! state message
      call fvl_loginfo("Computing test metrics...")

      ! ----------------------------------------------------------------------------
      ! reference time parameters
      ! ----------------------------------------------------------------------------

      ! inspect file
      inquire(file=fvl_trim(refreportfiledirectory(i))//"/time.fvd",exist=exist)

      ! file exists
      if(exist) then

            ! state message
            call fvl_loginfo("  >>> Reading reference time parameters")

            ! read file
            paramfile = fvl_keys_file(path=fvl_trim(refreportfiledirectory(i))//"/time.fvd",form=refreportfileform(i),action="read")
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
      inquire(file=fvl_trim(reportfiledirectory(i))//"/time.fvd",exist=exist)

      ! file exists
      if(exist) then

            ! state message
            call fvl_loginfo("  >>> Reading time parameters")

            ! read file
            paramfile = fvl_keys_file(path=fvl_trim(reportfiledirectory(i))//"/time.fvd",form=reportfileform(i),action="read")
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
      inquire(file=fvl_trim(refreportfiledirectory(i))//"/memory.fvd",exist=exist)

      ! file exists
      if(exist) then

            ! state message
            call fvl_loginfo("  >>> Reading reference memory parameters")

            ! read file
            paramfile = fvl_keys_file(path=fvl_trim(refreportfiledirectory(i))//"/memory.fvd",form=refreportfileform(i),action="read")
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
      inquire(file=fvl_trim(reportfiledirectory(i))//"/memory.fvd",exist=exist)

      ! file exists
      if(exist) then

            ! state message
            call fvl_loginfo("  >>> Reading memory parameters")

            ! read file
            paramfile = fvl_keys_file(path=fvl_trim(reportfiledirectory(i))//"/memory.fvd",form=reportfileform(i),action="read")
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
      inquire(file=fvl_trim(reportfiledirectory(i))//"/linear_solver.fvd",exist=exist)

      if(exist) then

            ! state message
            call fvl_loginfo("  >>> Reading solution parameters")

            ! read file
            file = fvl_file(path=fvl_trim(reportfiledirectory(i))//"/linear_solver.fvd",form=reportfileform(i),action="read")
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
      inquire(file=fvl_trim(reportfiledirectory(i))//"/fixpoint_solver.fvd",exist=exist)

      if(exist) then

            ! state message
            call fvl_loginfo("  >>> Reading fixpoint parameters")

            ! read file
            file = fvl_file(path=fvl_trim(reportfiledirectory(i))//"/fixpoint_solver.fvd",form=reportfileform(i),action="read")
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
realformat        = "(es20."//fvl_trim(fvl_char(tableprecision))//")"

! number of meshes
nummeshes         = 0
do i=1,maxnumtests
      if(numcells(i)/=0) then
            nummeshes = i
      end if
end do

! dofs
if(convergencetype==1) then
      numdofs(1:nummeshes) = numcells(1:nummeshes)
      doflabel = "N^{\text{P}}_{\text{C}}"
else
      numdofs(1:nummeshes) = steps(1:nummeshes)
      doflabel = "N^{\text{P}}_{\text{S}}"
end if

! ----------------------------------------------------------------------------
! errors table
! ----------------------------------------------------------------------------

! write errors table
file              = fvl_file(path=fvl_trim(tabledirectory)//"errors_table.tex",form="ascii",action="write")
call file         % open()
unit              = file%getunit()

write(unit,"(a)") "\begin{table}[!h]"
write(unit,"(a)") "\centering"
write(unit,"(a)") "\caption{Field magnitude errors and convergence orders.}"
write(unit,"(a)") "\label{}"
write(unit,"(a)") "\begin{tabular}{@{}rrrrrrrrrr@{}}"
write(unit,"(a)") "\toprule"
write(unit,"(a)") "$"//fvl_trim(doflabel)//"$ & & $E_{1}$ & $O_{1}$ & $E_{2}$ & $O_{2}$ & $E_{\infty}$ & $O_{\infty}$ \\"
write(unit,"(a)") "\midrule"
do i=1,nummeshes
      if(numdofs(i)==0) then
            write(unit,"(a)",advance="no") "--- & & "
            write(unit,"(a)",advance="no") "--- & --- & "
            ! write(unit,"(a)",advance="no") "--- & --- & "
            write(unit,"(a)") "--- & --- \\"
      else
            if(i==1) then
                  write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(numdofs(i)))//" & & "
                  write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(errormag1(i),tableprecision,0.0d0))//" & --- & "
                  ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(errormag2(i),tableprecision,0.0d0))//" & --- & "
                  write(unit,"(a)") fvl_trim(fvl_latexschar(errormaginf(i),tableprecision,0.0d0))//" & --- \\"
            else
                  write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(numdofs(i)))//" & & "
                  write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(errormag1(i),tableprecision,0.0d0))//" & "&
                        //fvl_trim(fvl_latexchar(ordermag1(i),tableprecision,0.0d0))//" & "
                  ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(errormag2(i),tableprecision,0.0d0))//" & "&
                  !      //fvl_trim(fvl_latexchar(ordermag2(i),tableprecision,0.0d0))//" & "
                  write(unit,"(a)") fvl_trim(fvl_latexschar(errormaginf(i),tableprecision,0.0d0))//" & "&
                        //fvl_trim(fvl_latexchar(ordermaginf(i),tableprecision,0.0d0))//" \\"
            end if
      end if
end do
write(unit,"(a)") "\bottomrule"
write(unit,"(a)") "\end{tabular}"
write(unit,"(a)") "\end{table}"
write(unit,"(a)") " "

if(fieldtype==2) then
      write(unit,"(a)") "\begin{table}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Field x-component errors and convergence orders.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tabular}{@{}rrrrrrrrrr@{}}"
      write(unit,"(a)") "\toprule"
      write(unit,"(a)") "$"//fvl_trim(doflabel)//"$ & & $E_{1}$ & $O_{1}$ & $E_{2}$ & $O_{2}$ & $E_{\infty}$ & $O_{\infty}$ \\"
      write(unit,"(a)") "\midrule"
      do i=1,nummeshes
            if(numdofs(i)==0) then
                  write(unit,"(a)",advance="no") "--- & & "
                  write(unit,"(a)",advance="no") "--- & --- & "
                  ! write(unit,"(a)",advance="no") "--- & --- & "
                  write(unit,"(a)") "--- & --- \\"
            else
                  if(i==1) then
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(errorx1(i),tableprecision,0.0d0))//" & --- & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(errorx2(i),tableprecision,0.0d0))//" & --- & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(errorxinf(i),tableprecision,0.0d0))//" & --- \\"
                  else
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(errorx1(i),tableprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(orderx1(i),tableprecision,0.0d0))//" & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(errorx2(i),tableprecision,0.0d0))//" & "&
                        !       //fvl_trim(fvl_latexchar(orderx2(i),tableprecision,0.0d0))//" & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(errorxinf(i),tableprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(orderxinf(i),tableprecision,0.0d0))//" \\"
                  end if
            end if
      end do
      write(unit,"(a)") "\bottomrule"
      write(unit,"(a)") "\end{tabular}"
      write(unit,"(a)") "\end{table}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{table}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Field y-component errors and convergence orders.}"
      write(unit,"(a)") "\label{}"
      write(unit,"(a)") "\begin{tabular}{@{}rrrrrrrrrr@{}}"
      write(unit,"(a)") "\toprule"
      write(unit,"(a)") "$"//fvl_trim(doflabel)//"$ & & $E_{1}$ & $O_{1}$ & $E_{2}$ & $O_{2}$ & $E_{\infty}$ & $O_{\infty}$ \\"
      write(unit,"(a)") "\midrule"
      do i=1,nummeshes
            if(numdofs(i)==0) then
                  write(unit,"(a)",advance="no") "--- & & "
                  write(unit,"(a)",advance="no") "--- & --- & "
                  ! write(unit,"(a)",advance="no") "--- & --- & "
                  write(unit,"(a)") "--- & --- \\"
            else
                  if(i==1) then
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(errory1(i),tableprecision,0.0d0))//" & --- & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(errory2(i),tableprecision,0.0d0))//" & --- & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(erroryinf(i),tableprecision,0.0d0))//" & --- \\"
                  else
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(numdofs(i)))//" & & "
                        write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(errory1(i),tableprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(ordery1(i),tableprecision,0.0d0))//" & "
                        ! write(unit,"(a)",advance="no") fvl_trim(fvl_latexschar(errory2(i),tableprecision,0.0d0))//" & "&
                        !       //fvl_trim(fvl_latexchar(ordery2(i),tableprecision,0.0d0))//" & "
                        write(unit,"(a)") fvl_trim(fvl_latexschar(erroryinf(i),tableprecision,0.0d0))//" & "&
                              //fvl_trim(fvl_latexchar(orderyinf(i),tableprecision,0.0d0))//" \\"
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
file              = fvl_file(path=fvl_trim(tabledirectory)//"solution_table.tex",form="ascii",action="write")
call file         % open()
unit              = file%getunit()
write(unit,"(a)") "\begin{table}[!h]"
write(unit,"(a)") "\centering"
write(unit,"(a)") "\caption{Iterations and residuals.}"
write(unit,"(a)") "\label{}"
write(unit,"(a)") "\begin{tabular}{@{}rrrrrr@{}}"
write(unit,"(a)") "\toprule"
write(unit,"(a)") "$"//fvl_trim(doflabel)//"$ & & $N_{\text{fixpoint}}$ & $R_{\text{fixpoint}}$ & $N_{\text{linearsystem}}$ & $R_{\text{linearsystem}}$ \\"
write(unit,"(a)") "\midrule"
do i=1,nummeshes
      if(numdofs(i)==0) then
            write(unit,"(a)",advance="no") "--- & & "
            write(unit,"(a)",advance="no") "--- & --- & "
            write(unit,"(a)") "--- & --- \\"
      else
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(numdofs(i)))//" & & "
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(fixpointnumiters(i)))//" & "&
                  //fvl_trim(fvl_latexchar(fixpointresid(i),tableprecision,0.0d0))//" & "
            write(unit,"(a)") fvl_trim(fvl_latexchar(linearsystemnumiters(i)))//" & "&
                  //fvl_trim(fvl_latexchar(linearsystemresid(i),tableprecision,0.0d0))//" \\"
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
file              = fvl_file(path=fvl_trim(tabledirectory)//"time_table.tex",form="ascii",action="write")
call file         % open()
unit              = file%getunit()
write(unit,"(a)") "\begin{table}[!h]"
write(unit,"(a)") "\centering"
write(unit,"(a)") "\caption{Time.}"
write(unit,"(a)") "\label{}"
write(unit,"(a)") "\begin{tabular}{@{}rrrrrrrrrrrrrrr@{}}"
write(unit,"(a)") "\toprule"
write(unit,"(a)",advance="no") "$"//fvl_trim(doflabel)//"$ & & "
write(unit,"(a)",advance="no") "$T_{\text{meshes}}$ & $R_{\text{meshes}}$ & "
write(unit,"(a)",advance="no") "$T_{\text{models}}$ & $R_{\text{models}}$ & "
write(unit,"(a)",advance="no") "$T_{\text{schemes}}$ & $R_{\text{schemes}}$ & "
write(unit,"(a)",advance="no") "$T_{\text{solution}}$ & $R_{\text{solution}}$ & "
write(unit,"(a)",advance="no") "$T_{\text{postproc}}$ & $R_{\text{postproc}}$ & "
write(unit,"(a)") "$T_{\text{execution}}$ & $R_{\text{execution}}$ \\"
write(unit,"(a)") "\midrule"
do i=1,nummeshes
      if(numdofs(i)==0) then
            write(unit,"(a)",advance="no") "--- & & "
            write(unit,"(a)",advance="no") "--- & --- & "
            write(unit,"(a)",advance="no") "--- & --- & "
            write(unit,"(a)",advance="no") "--- & --- & "
            write(unit,"(a)",advance="no") "--- & --- & "
            write(unit,"(a)",advance="no") "--- & --- & "
            write(unit,"(a)") "--- & --- \\"
      else
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(numdofs(i)))//" & & "
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(meshestotaltime(i),tableprecision,0.0d0))//" & "&
                  //fvl_trim(fvl_latexchar(meshestotaltimeratio(i),tableprecision,0.0d0))//" & "
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(modelstotaltime(i),tableprecision,0.0d0))//" & "&
                  //fvl_trim(fvl_latexchar(modelstotaltimeratio(i),tableprecision,0.0d0))//" & "
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(schemestotaltime(i),tableprecision,0.0d0))//" & "&
                  //fvl_trim(fvl_latexchar(schemestotaltimeratio(i),tableprecision,0.0d0))//" & "
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(solutiontotaltime(i),tableprecision,0.0d0))//" & "&
                  //fvl_trim(fvl_latexchar(solutiontotaltimeratio(i),tableprecision,0.0d0))//" & "
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(postproctotaltime(i),tableprecision,0.0d0))//" & "&
                  //fvl_trim(fvl_latexchar(postproctotaltimeratio(i),tableprecision,0.0d0))//" & "
            write(unit,"(a)") fvl_trim(fvl_latexchar(executiontotaltime(i),tableprecision,0.0d0))//" & "&
                  //fvl_trim(fvl_latexchar(executiontotaltimeratio(i),tableprecision,0.0d0))//" \\"
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
file              = fvl_file(path=fvl_trim(tabledirectory)//"memory_table.tex",form="ascii",action="write")
call file         % open()
unit              = file%getunit()
write(unit,"(a)") "\begin{table}[!h]"
write(unit,"(a)") "\centering"
write(unit,"(a)") "\caption{Memory usage.}"
write(unit,"(a)") "\label{}"
write(unit,"(a)") "\begin{tabular}{@{}rrrr@{}}"
write(unit,"(a)") "\toprule"
write(unit,"(a)") "$"//fvl_trim(doflabel)//"$ & & $M_{\text{execution}}$ & $R_{\text{execution}}$ \\"
write(unit,"(a)") "\midrule"
do i=1,nummeshes
      if(numdofs(i)==0) then
            write(unit,"(a)") "--- & & --- & --- \\"
      else
            write(unit,"(a)",advance="no") fvl_trim(fvl_latexchar(numdofs(i)))//" & & "
            write(unit,"(a)") fvl_trim(fvl_latexchar(totalexecutionmemory(i),tableprecision,0.0d0))//" & "&
                  //fvl_trim(fvl_latexchar(totalexecutionmemoryratio(i),tableprecision,0.0d0))//" \\"
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
file              = fvl_file(path=fvl_trim(tabledirectory)//"errors_plot.tex",form="ascii",action="write")
call file         % open()
unit              = file%getunit()

write(unit,"(a)") "\begin{figure}[!h]"
write(unit,"(a)") "\centering"
write(unit,"(a)") "\caption{Field magnitude errors \textit{versus} mumber of degrees-of-freedom.}"
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
write(unit,"(a)") "xlabel={$"//fvl_trim(doflabel)//"$},"
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
      write(unit,"(a)") "("//fvl_trim(fvl_char(numdofs(i),intformat))//","//fvl_trim(fvl_char(errormag1(i),realformat))//")"
end do
write(unit,"(a)") "};"
write(unit,"(a)") "\addlegendentry{$E_{1}$}"

write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
do i=1,nummeshes
      write(unit,"(a)") "("//fvl_trim(fvl_char(numdofs(i),intformat))//","//fvl_trim(fvl_char(errormag2(i),realformat))//")"
end do
write(unit,"(a)") "};"
write(unit,"(a)") "\addlegendentry{$E_{2}$}"

write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
do i=1,nummeshes
      write(unit,"(a)") "("//fvl_trim(fvl_char(numdofs(i),intformat))//","//fvl_trim(fvl_char(errormaginf(i),realformat))//")"
end do
write(unit,"(a)") "};"
write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
write(unit,"(a)") "\end{loglogaxis}"
write(unit,"(a)") "\end{tikzpicture}"
write(unit,"(a)") "\end{figure}"
write(unit,"(a)") " "

if(fieldtype==2) then
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Field x-component errors \textit{versus} mumber of degrees-of-freedom.}"
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
      write(unit,"(a)") "xlabel={$"//fvl_trim(doflabel)//"$},"
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
            write(unit,"(a)") "("//fvl_trim(fvl_char(numdofs(i),intformat))//","//fvl_trim(fvl_char(errorx1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(numdofs(i),intformat))//","//fvl_trim(fvl_char(errorx2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(numdofs(i),intformat))//","//fvl_trim(fvl_char(errorxinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Field y-component errors \textit{versus} mumber of degrees-of-freedom.}"
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
      write(unit,"(a)") "xlabel={$"//fvl_trim(doflabel)//"$},"
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
            write(unit,"(a)") "("//fvl_trim(fvl_char(numdofs(i),intformat))//","//fvl_trim(fvl_char(errory1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(numdofs(i),intformat))//","//fvl_trim(fvl_char(errory2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(numdofs(i),intformat))//","//fvl_trim(fvl_char(erroryinf(i),realformat))//")"
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
file              = fvl_file(path=fvl_trim(tabledirectory)//"time_plot.tex",form="ascii",action="write")
call file         % open()
unit              = file%getunit()

write(unit,"(a)") "\begin{figure}[!h]"
write(unit,"(a)") "\centering"
write(unit,"(a)") "\caption{Field magnitude errors \textit{versus} total execution time.}"
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
      write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(errormag1(i),realformat))//")"
end do
write(unit,"(a)") "};"
write(unit,"(a)") "\addlegendentry{$E_{1}$}"
write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
do i=1,nummeshes
      write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(errormag2(i),realformat))//")"
end do
write(unit,"(a)") "};"
write(unit,"(a)") "\addlegendentry{$E_{2}$}"
write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
do i=1,nummeshes
      write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(errormaginf(i),realformat))//")"
end do
write(unit,"(a)") "};"
write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
write(unit,"(a)") "\end{loglogaxis}"
write(unit,"(a)") "\end{tikzpicture}"
write(unit,"(a)") "\end{figure}"
write(unit,"(a)") " "

if(fieldtype==2) then
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Field x-component errors \textit{versus} total execution time.}"
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
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(errorx1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(errorx2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(errorxinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Field y-component errors \textit{versus} total execution time.}"
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
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(errory1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(errory2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(executiontotaltime(i),realformat))//","//fvl_trim(fvl_char(erroryinf(i),realformat))//")"
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
file              = fvl_file(path=fvl_trim(tabledirectory)//"memory_plot.tex",form="ascii",action="write")
call file         % open()
unit              = file%getunit()

write(unit,"(a)") "\begin{figure}[!h]"
write(unit,"(a)") "\centering"
write(unit,"(a)") "\caption{Field magnitude errors \textit{versus} total execution memory.}"
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
      write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(errormag1(i),realformat))//")"
end do
write(unit,"(a)") "};"
write(unit,"(a)") "\addlegendentry{$E_{1}$}"
write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
do i=1,nummeshes
      write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(errormag2(i),realformat))//")"
end do
write(unit,"(a)") "};"
write(unit,"(a)") "\addlegendentry{$E_{2}$}"
write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
do i=1,nummeshes
      write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(errormaginf(i),realformat))//")"
end do
write(unit,"(a)") "};"
write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
write(unit,"(a)") "\end{loglogaxis}"
write(unit,"(a)") "\end{tikzpicture}"
write(unit,"(a)") "\end{figure}"
write(unit,"(a)") " "

if(fieldtype==2) then
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Field x-component errors \textit{versus} total execution memory.}"
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
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(errorx1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(errorx2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(errorxinf(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{\infty}$}"
      write(unit,"(a)") "\end{loglogaxis}"
      write(unit,"(a)") "\end{tikzpicture}"
      write(unit,"(a)") "\end{figure}"
      write(unit,"(a)") " "
      write(unit,"(a)") "\begin{figure}[!h]"
      write(unit,"(a)") "\centering"
      write(unit,"(a)") "\caption{Field y-component errors \textit{versus} total execution memory.}"
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
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(errory1(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{1}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=square,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(errory2(i),realformat))//")"
      end do
      write(unit,"(a)") "};"
      write(unit,"(a)") "\addlegendentry{$E_{2}$}"
      write(unit,"(a)") "\addplot[style=thick,solid,mark=triangle,mark size=2pt,mark options={solid},color=black] coordinates{"
      do i=1,nummeshes
            write(unit,"(a)") "("//fvl_trim(fvl_char(totalexecutionmemory(i),realformat))//","//fvl_trim(fvl_char(erroryinf(i),realformat))//")"
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
file              = fvl_file(path=fvl_trim(tabledirectory)//"report.tex",form="ascii",action="write")
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
write(unit,"(a)") "\input{"//fvl_trim(fvl_basename(tabledirectory))//"errors_table.tex}"
write(unit,"(a)") "\input{"//fvl_trim(fvl_basename(tabledirectory))//"solution_table.tex}"
write(unit,"(a)") "\input{"//fvl_trim(fvl_basename(tabledirectory))//"time_table.tex}"
write(unit,"(a)") "\input{"//fvl_trim(fvl_basename(tabledirectory))//"memory_table.tex}"
write(unit,"(a)") "\input{"//fvl_trim(fvl_basename(tabledirectory))//"errors_plot.tex}"
write(unit,"(a)") "\input{"//fvl_trim(fvl_basename(tabledirectory))//"time_plot.tex}"
write(unit,"(a)") "\input{"//fvl_trim(fvl_basename(tabledirectory))//"memory_plot.tex}"
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

end program fvl_post_thermoincompflow2d_main
! end of file
