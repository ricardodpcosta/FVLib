! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: 3D non-isothermal Newtonian incompressible fluid flow solver
! Modification: February, 2025

#include "macros.f"

program fvl_pro_thermoincompflow3d_main

use fvl_lib3d

implicit none

! ============================================================================
! DECLARE VARIABLES
! ============================================================================

! ----------------------------------------------------------------------------
! variables
! ----------------------------------------------------------------------------

integer(kind=__fvl_integer_kind__)::i,p_dofs,u_dofs,&
      staticmumodel,staticrhomodel,staticnumodel,staticpmodel,staticumodel,staticvmodel,staticfupmodel,staticgupmodel
real(kind=__fvl_real_kind__),allocatable,dimension(:)::p_staticsource,ux_staticsource
real(kind=__fvl_real_kind__),allocatable,dimension(:)::p_globalsource,ux_globalsource
real(kind=__fvl_real_kind__),allocatable,dimension(:)::p_staticrhs,ux_staticrhs
real(kind=__fvl_real_kind__),allocatable,dimension(:)::p_globalrhs,ux_globalrhs
type(fvl_calc_scalarfield3d)::p_innerfield
type(fvl_calc_vectorfield3d)::u_innerfield,v_innerfield
type(fvl_const_scalarfield3d)::xi_innerfield
type(fvl_diamondpro_gradientscheme3d)::p_gradscheme
type(fvl_centredpro_interpolate3d)::ux_interpolate
type(fvl_centredpro_integratescheme3d)::ux_integscheme
type(fvl_centredpro_vectorlaplacianscheme3d)::ux_lapscheme
type(fvl_upwindpro_vectordivergencescheme3d)::ux_divscheme
type(fvl_primalpro_divergencescheme3d)::u_divscheme
type(fvl_primalpro_interpolate3d)::u_interpolate
type(fvl_lil_rspmat)::ux_staticmat1,ux_integmat,ux_globalmat1,p_precondmat,ux_precondmat
type(fvl_lil_brspmat)::p_globalmat,ux_staticmat2,ux_globalmat2,u_globalmat
type(fvl_csr_rspmat)::ux_staticmat_comp1,ux_integmat_comp,ux_globalmat_comp1
type(fvl_csr_brspmat)::p_globalmat_comp,ux_staticmat_comp2,ux_globalmat_comp2,u_globalmat_comp
type(fvl_linearsolver)::linearsolver
type(fvl_schur_preconditioner)::preconditioner
type(fvl_fixpointsolver)::fixpointsolver
type(fvl_basic_timesolver)::timesolver

! ----------------------------------------------------------------------------
! user-defined variables
! ----------------------------------------------------------------------------

#include "variables.f"

! ----------------------------------------------------------------------------
! reports variables
! ----------------------------------------------------------------------------

#include "fvl_writereports_hdr.f"

! ----------------------------------------------------------------------------
! dictionaries variables
! ----------------------------------------------------------------------------

#include "fvl_readdictionaries_hdr.f"

! ----------------------------------------------------------------------------
! mesh variables
! ----------------------------------------------------------------------------

#include "fvl_creatediamondmesh3d_hdr.f"

! ----------------------------------------------------------------------------
! fields variables
! ----------------------------------------------------------------------------

#include "fvl_writefields_hdr.f"

! ----------------------------------------------------------------------------
! generic model variables
! ----------------------------------------------------------------------------

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ k
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ cp
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ mu
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ rho
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ t
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field3d_cellmeansfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ p
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field3d_cellmeansfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ u
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field3d_cellmeansfield
#include "fvl_vectorfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ v
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_vectorfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ fup
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field3d_cellmeansfield
#include "fvl_vectorfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ gup
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field3d_cellmeansfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ alpha
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ nu
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ vp
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_scalarfield_models3d_hdr.f"

! ----------------------------------------------------------------------------
! specific model functions
! ----------------------------------------------------------------------------

! #include "fvl_p_models3d_hdr.f"
! #include "fvl_u_models3d_hdr.f"

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
      call fvl_loginfo("About: 3D non-isothermal Newtonian incompressible fluid flow solver.")
      call fvl_loginfo("Usage: fvl-pro-thermoincompflow3d [options]")
      call fvl_loginfo("Options:")
      call fvl_loginfo("      -h|--help               Displays help message.")
      call fvl_loginfo("      -m|--meshes <file>      Reads meshes parameters from dictionary file <file> (default is setup/meshes.fvd).")
      call fvl_loginfo("      -o|--models <file>      Reads models parameters from dictionary file <file> (default is setup/models.fvd).")
      call fvl_loginfo("      -s|--schemes <file>     Reads schemes parameters from dictionary file <file> (default is setup/schemes.fvd).")
      call fvl_loginfo("      -l|--solution <file>    Reads solution parameters from dictionary file <file> (default is setup/solution.fvd).")
      call fvl_loginfo("      -c|--control <file>     Reads control parameters from dictionary file <file> (default is setup/control.fvd).")
      call fvl_loginfo("      -p|--parallel           Runs in parallel with <OMP_NUM_THREADS> processes (default is serial).")
      stop
end if

! start total execution timer
call fvl_startwtimer(10,"execution_total")

! ============================================================================
! READ CASE PARAMETERS
! ============================================================================

! state message
call fvl_loginfo("Reading case parameters...")

#include "fvl_readdictionaries_src.f"

staticmumodel           = controldictfile%getvalue("solver","static_mu_model",0)
staticrhomodel          = controldictfile%getvalue("solver","static_rho_model",0)
staticpmodel            = controldictfile%getvalue("solver","static_p_model",0)
staticumodel            = controldictfile%getvalue("solver","static_u_model",0)
staticvmodel            = controldictfile%getvalue("solver","static_v_model",0)
staticfupmodel          = controldictfile%getvalue("solver","static_fup_model",0)
staticgupmodel          = controldictfile%getvalue("solver","static_gup_model",0)
staticnumodel           = staticmumodel*staticrhomodel

! set time directories precision
call fvl_settimedirectoryprecision(solutiontimeprecision)

! ============================================================================
! INITIALIZE MESHES
! ============================================================================

! state message
call fvl_loginfo("Initializing meshes...")

! start meshes initialization timer
call fvl_startwtimer(20,"meshes_total")

! ----------------------------------------------------------------------------
! mesh
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Initializing mesh")

#include "fvl_creatediamondmesh3d_src.f"

! stop meshes initialization timer
call fvl_stopwtimer(20)

! optimize memory
call fvl_optimize_memory1()

! save memory usage
call fvl_setvmhwm(1,"meshes_total")

! ----------------------------------------------------------------------------
! user-defined mesh
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Setting user-defined mesh")

#include "mesh.f"

! ----------------------------------------------------------------------------
! allocate vectors
! ----------------------------------------------------------------------------

! initialize variables
p_dofs                  = mesh%getnumcells()
u_dofs                  = 3*diamondmesh%getnumcells()

! allocate memory
call fvl_allocate(p_staticsource,p_dofs)
call fvl_allocate(p_globalsource,p_dofs)
call fvl_allocate(ux_staticsource,2,u_dofs/2)
call fvl_allocate(ux_globalsource,2,u_dofs/2)
call fvl_allocate(p_staticrhs,p_dofs)
call fvl_allocate(p_globalrhs,p_dofs)
call fvl_allocate(ux_staticrhs,2,u_dofs/2)
call fvl_allocate(ux_globalrhs,2,u_dofs/2)

! clean memory
call fvl_omp_clean(p_staticsource)
call fvl_omp_clean(p_globalsource)
call fvl_omp_clean(ux_staticsource)
call fvl_omp_clean(ux_globalsource)
call fvl_omp_clean(p_staticrhs)
call fvl_omp_clean(p_globalrhs)
call fvl_omp_clean(ux_staticrhs)
call fvl_omp_clean(ux_globalrhs)
call fvl_clean(p_globalmat)
call fvl_clean(ux_staticmat1)
call fvl_clean(ux_staticmat2)
call fvl_clean(ux_integmat)
call fvl_clean(ux_globalmat1)
call fvl_clean(ux_globalmat2)
call fvl_clean(u_globalmat)
call fvl_clean(p_globalmat_comp)
call fvl_clean(ux_staticmat_comp1)
call fvl_clean(ux_staticmat_comp2)
call fvl_clean(ux_integmat_comp)
call fvl_clean(ux_globalmat_comp1)
call fvl_clean(ux_globalmat_comp2)
call fvl_clean(u_globalmat_comp)
call fvl_clean(p_precondmat)
call fvl_clean(ux_precondmat)

! ============================================================================
! INITIALIZE PHYSICAL MODEL
! ============================================================================

! state message
call fvl_loginfo("Initializing physical models...")

! start models initialization timer
call fvl_startwtimer(30,"models_total")

! ----------------------------------------------------------------------------
! physical models
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Initializing physical models")

! dynamic viscosity model
mu_model                = fvl_scalarmodel3d(mesh=diamondmesh,fieldtype=fvl_field3d_edgequadratsfield,&
                              dictfile=modelsdictfile,dictlabel="mu")

! density model
rho_model               = fvl_scalarmodel3d(mesh=diamondmesh,fieldtype=fvl_field3d_edgequadratsfield,&
                              dictfile=modelsdictfile,dictlabel="rho")

! pressure model
p_model                 = fvl_scalarmodel3d(mesh=mesh,fieldtype=fvl_field3d_cellmeansfield,&
                              dictfile=modelsdictfile,dictlabel="p")

! velocity model
u_model                 = fvl_vectormodel3d(mesh=diamondmesh,fieldtype=fvl_field3d_cellmeansfield,&
                              dictfile=modelsdictfile,dictlabel="u")

! velocity model
v_model                 = fvl_vectormodel3d(mesh=diamondmesh,fieldtype=fvl_field3d_edgequadratsfield,&
                              dictfile=modelsdictfile,dictlabel="v")

! forces model
fup_model               = fvl_vectormodel3d(mesh=diamondmesh,fieldtype=fvl_field3d_cellmeansfield,&
                              dictfile=modelsdictfile,dictlabel="fup")

! mass model
gup_model               = fvl_scalarmodel3d(mesh=mesh,fieldtype=fvl_field3d_cellmeansfield,&
                              dictfile=modelsdictfile,dictlabel="gup")

! ----------------------------------------------------------------------------
! user-defined models
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Setting user-defined models")

#include "model.f"

! ! initialize pressure pseudo patches
! call p_initialize_pseudo_patches()

! ! initialize velocity pseudo patches
! call u_initialize_pseudo_patches()

! ----------------------------------------------------------------------------
! computed fields
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Initializing computed fields")

! kinematic viscosity model
nu_model                = -mu_model/rho_model

! specific volume model
vp_model                = 1.0d0/rho_model

! pressure field
p_innerfield            = fvl_calc_scalarfield3d(patch=cellspatch,type=fvl_field3d_cellmeansfield,&
                              label="p",filepath=solutionfilesdir,fileform=solutionfilesform)

! velocity field
u_innerfield            = fvl_calc_vectorfield3d(patch=diamond_cellspatch,type=fvl_field3d_cellmeansfield,&
                              label="u_diamond",filepath=solutionfilesdir,fileform=solutionfilesform)

! xi field
xi_innerfield           = fvl_const_scalarfield3d(patch=cellspatch,type=fvl_field3d_cellmeansfield,&
                              label="xi",filepath=solutionfilesdir,fileform=solutionfilesform)

! velocity field
v_innerfield            = fvl_calc_vectorfield3d(patch=diamond_edgespatch,type=fvl_field3d_edgequadratsfield,&
                              label="v",filepath=solutionfilesdir,fileform=solutionfilesform)

! set fields
call solutionfields(1)  % setptr(u_innerfield)
call solutionfields(2)  % setptr(p_innerfield)
call solutionfields(3)  % setptr(xi_innerfield)

! intialize pressure field
call p_model            % evalcellmeanvalues(p_innerfield)

! intialize velocity field
call u_model            % evalcellmeanvalues(u_innerfield)

! velocity model
call v_model            % setfield(v_innerfield)

! ----------------------------------------------------------------------------
! write fields
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Writing fields")

#include "fvl_writefields_src.f"

! stop models initialization timer
call fvl_stopwtimer(30)

! save memory usage
call fvl_setvmhwm(2,"models_total")

! ============================================================================
! COMPUTE NUMERICAL SCHEMES
! ============================================================================

! state message
call fvl_loginfo("Computing numerical schemes...")

! start schemes initialization timer
call fvl_startwtimer(40,"schemes_total")

! ----------------------------------------------------------------------------
! pressure gradient scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Computing gradient(1/rho,p) scheme")

! start timer
call fvl_startwtimer(41,"gradient(1/rho,p)_scheme")

! initialize scheme
p_gradscheme            = fvl_diamondpro_gradientscheme3d(coefsmodel=vp_model,condsmodel=p_model,&
                              mesh=diamondmesh,primalmesh=mesh,dictfile=schemesdictfile,dictlabel="gradient(1/rho,p)",&
                              staticcoeffields=(staticrhomodel==1),staticboundconds=.true.)

! evaluate coefficients matrix
if(staticrhomodel==1) then
      call p_gradscheme       % evalmatrix(p_globalmat)
end if

! ! optimize memory
! call p_gradscheme%deallocatelhs()

! compress coefficients matrix
if(staticrhomodel==1) then
      call p_globalmat_comp%compress(p_globalmat)
      call fvl_delete(p_globalmat)
end if

! stop timer
call fvl_stopwtimer(41)

! ----------------------------------------------------------------------------
! velocity interpolation scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Computing interpolate(u) scheme")

! start timer
call fvl_startwtimer(42,"interpolate(u)_scheme")

! initialize scheme
ux_interpolate          = fvl_centredpro_interpolate3d(patch=diamond_facespatch,mesh=diamondmesh,&
                              type=fvl_centredpro_interpolate3d_edgequadratvalues,&
                              dictfile=schemesdictfile,dictlabel="interpolate(u)")

! deallocate memory
! call fvl_delete(ux_interpolate)

! stop timer
call fvl_stopwtimer(42)

! ----------------------------------------------------------------------------
! velocity integrate scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Computing integrate(u) scheme")

! start timer
call fvl_startwtimer(43,"integrate(u)_scheme")

! initialize scheme
ux_integscheme          = fvl_centredpro_integratescheme3d(mesh=diamondmesh,&
                              dictfile=schemesdictfile,dictlabel="integrate(u)",&
                              staticcoeffields=.true.)

! evaluate coefficients matrix
call ux_integscheme    % evalmatrix(ux_integmat)

! compress coefficients matrix
call ux_integmat_comp   % compress(ux_integmat)

! deallocate memory
call fvl_delete(ux_integscheme)
call fvl_delete(ux_integmat)

! stop timer
call fvl_stopwtimer(43)

! ----------------------------------------------------------------------------
! velocity laplacian scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Computing laplacian(nu,u) scheme")

! start timer
call fvl_startwtimer(44,"laplacian(nu,u)_scheme")

! initialize scheme
ux_lapscheme            = fvl_centredpro_vectorlaplacianscheme3d(coefsmodel=nu_model,condsmodel=u_model,&
                              mesh=diamondmesh,dictfile=schemesdictfile,dictlabel="laplacian(nu,u)",&
                              staticcoeffields=(staticnumodel==1),staticboundconds=.true.)

! evaluate coefficients matrix
if(staticnumodel==1) then
      call ux_lapscheme       % evalmatrix(ux_staticmat1,ux_staticmat2)
end if

! ! optimize memory
! call u_lapscheme%deallocatelhs()

! stop timer
call fvl_stopwtimer(44)

! ----------------------------------------------------------------------------
! velocity divergence scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Computing divergence(v,u) scheme")

! start timer
call fvl_startwtimer(45,"divergence(v,u)_scheme")

! initialize scheme
ux_divscheme            = fvl_upwindpro_vectordivergencescheme3d(coefsmodel=v_model,condsmodel=u_model,&
                              mesh=diamondmesh,dictfile=schemesdictfile,dictlabel="divergence(v,u)",&
                              staticcoeffields=(staticvmodel==1),staticboundconds=.true.)

! evaluate coefficients matrix
if(staticvmodel==1) then
      call ux_divscheme       % evaladdmatrix(ux_staticmat1,ux_staticmat2)
end if

! ! optimize memory
! call u_divscheme%deallocatelhs()

! compress coefficients matrix
if(staticnumodel==1 .or. staticvmodel==1) then
      call ux_staticmat_comp1%compress(ux_staticmat1)
      call fvl_delete(ux_staticmat1)
      call ux_staticmat_comp2%compress(ux_staticmat2)
      call fvl_delete(ux_staticmat2)
end if

! stop timer
call fvl_stopwtimer(45)

! ----------------------------------------------------------------------------
! velocity divergence scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Computing divergence(u) scheme")

! start timer
call fvl_startwtimer(46,"divergence(u)_scheme")

! initialize scheme
u_divscheme         = fvl_primalpro_divergencescheme3d(condsmodel=u_model,&
                              mesh=mesh,diamondmesh=diamondmesh,dictfile=schemesdictfile,dictlabel="divergence(u)",&
                              staticcoeffields=.true.,staticboundconds=.true.)

! evaluate coefficients matrix
call u_divscheme    % evaladdmatrix(u_globalmat)

! ! optimize memory
! call u_divscheme%deallocatelhs()

! compress coefficients matrix
call u_globalmat_comp%compress(u_globalmat)
call fvl_delete(u_globalmat)

! stop timer
call fvl_stopwtimer(46)

! stop timer
call fvl_stopwtimer(40)

! save memory usage
call fvl_setvmhwm(3,"schemes_total")

! ============================================================================
! COMPUTE SOLUTION
! ============================================================================

! state message
call fvl_loginfo("Computing solution...")

! start solution timer
call fvl_startwtimer(50,"solution_total")

! linear solver
linearsolver            = fvl_linearsolver(dictfile=solutiondictfile,dictlabel="linear_solver")

! preconditioner
preconditioner          = fvl_schur_preconditioner(dictfile=solutiondictfile,dictlabel="preconditioner")

! fixpoint solver
fixpointsolver          = fvl_fixpointsolver(dictfile=solutiondictfile,dictlabel="fixpoint_solver")

! time solver
timesolver              = fvl_basic_timesolver(runtimes=runtimes,fields=solutionfields(1:3),&
                              tfun=tfun,afun=afun,bfun=bfun,pfun=pfun,updatefun=updatefun,&
                              dictfile=solutiondictfile,dictlabel="time_solver")

! timesolver              = fvl_multistep_timesolver(runtime=mesh,field=tz_innerfield,tfun=tfun,afun=afun,bfun=bfun,pfun=pfun,&
!                               updatefun=updatefun,dictfile=solutiondictfile,dictlabel="time_solver")

! timesolver              = fvl_compact_timesolver(runtime=mesh,zfield=tz_innerfield,dfield=td_innerfield,sfield=ts_innerfield,&
!                               matrices1=t_globalmat_comp,matrices2=t_globalmat_comp,&
!                               tfun1=tfun,afun1=afun,bfun1=bfun,pfun1=pfun,tfun2=tfun,afun2=afun,bfun2=bfun2,pfun2=pfun,&
!                               updatefun=updatefun,staticcoefsmatrix=(staticnumodel==1 .and. staticpsimodel==1),&
!                               dictfile=solutiondictfile,dictlabel="time_solver")

! compute solution
call timesolver         % solve(linearsolver,fixpointsolver)

! stop solution timer
call fvl_stopwtimer(50)

! save memory usage
call fvl_setvmhwm(4,"solution_total")

! ----------------------------------------------------------------------------
! write histograms
! ----------------------------------------------------------------------------

! linear system
call linearsolver       % writehistograms(fvl_trim(solutionfilesdir)//"/linear_solver.fvd","ascii")

! fixpoint solver
call fixpointsolver     % writehistograms(fvl_trim(solutionfilesdir)//"/fixpoint_solver.fvd","ascii")

! optimize memory
call fvl_optimize_memory4()

! ============================================================================
! POST-PROCESSING DATA
! ============================================================================

! state message
call fvl_loginfo("Post-processing data...")

! start post-processing timer
call fvl_startwtimer(60,"postprocessing_total")

! ----------------------------------------------------------------------------
! interpolate velocity from diamond to primal mesh
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Interpolating velocity field")

! velocity field
v_innerfield            = fvl_calc_vectorfield3d(patch=cellspatch,type=fvl_field3d_cellmeansfield,&
                              label="u",filepath=solutionfilesdir,fileform=solutionfilesform)

! velocity interpolation scheme
u_interpolate           = fvl_primalpro_interpolate3d(patch=cellspatch,mesh=mesh,diamondmesh=diamondmesh,&
                              type=fvl_primalpro_interpolate3d_cellmeanvalues,&
                              dictfile=schemesdictfile,dictlabel="interpolate(u)")

! loop over time directories
do i=1,numtimedirectories

      ! set step and time
      call mesh               % setcurrentstep(steps(i))
      call mesh               % setcurrenttime(times(i))
      call diamondmesh        % setcurrentstep(steps(i))
      call diamondmesh        % setcurrenttime(times(i))

      ! read velocity field
      call u_innerfield       % read()

      ! interpolate velocity field
      call u_interpolate      % evalvalues(u_innerfield,v_innerfield)

      ! write velocity field
      call v_innerfield       % write()

end do

! deallocate memory
call fvl_delete(u_interpolate)

! ----------------------------------------------------------------------------
! user-defined post-processing
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Running user-defined post-processing")

#include "post.f"

! stop post-processing timer
call fvl_stopwtimer(60)

! save memory usage
call fvl_setvmhwm(5,"postprocessing_total")

! ============================================================================
! STOP LIBRARY
! ============================================================================

! stop total execution timer
call fvl_stopwtimer(10)

! save memory usage
call fvl_setvmhwm(6,"execution_total")

! ----------------------------------------------------------------------------
! write reports
! ----------------------------------------------------------------------------

#include "fvl_writereports_src.f"

! stop library
call fvl_stop()

contains

! ============================================================================
! INCLUDE FUNCTIONS
! ============================================================================

! ----------------------------------------------------------------------------
! generic model functions
! ----------------------------------------------------------------------------

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ k
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ cp
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ mu
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ rho
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ t
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field3d_cellmeansfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ p
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field3d_cellmeansfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ u
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field3d_cellmeansfield
#include "fvl_vectorfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ v
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_vectorfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ fup
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field3d_cellmeansfield
#include "fvl_vectorfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ gup
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field3d_cellmeansfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ alpha
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ nu
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_scalarfield_models3d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ vp
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field3d_edgequadratsfield
#include "fvl_scalarfield_models3d_hdr.f"

! ----------------------------------------------------------------------------
! specific model functions
! ----------------------------------------------------------------------------

! #include "fvl_u_models3d_src.f"
! #include "fvl_p_models3d_src.f"

! ----------------------------------------------------------------------------
! user-defined functions
! ----------------------------------------------------------------------------

#include "functions.f"

! ----------------------------------------------------------------------------
! built-in functions
! ----------------------------------------------------------------------------

subroutine tfun(x,res)
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:size(x,1)),intent(out)::res

      ! start timer
      call fvl_startwtimer(51,"solution_tfun")

      ! momentum equation
      if(ux_integmat_comp%getmaxnumrowvalues()/=0) then
            call ux_integmat_comp%blockmul3(x(1:u_dofs),res(1:u_dofs))
      else
            call fvl_omp_assign(x(1:u_dofs),res(1:u_dofs))
      end if

      ! continuity equation
      call fvl_omp_assign(0.0d0,res(u_dofs+1:p_dofs+u_dofs))

      ! compatibility equation
      res(p_dofs+u_dofs+1)=0.0d0

      ! stop timer
      call fvl_stopwtimer(51)
end subroutine tfun

subroutine afun(x,res)
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:size(x,1)),intent(out)::res

      ! start timer
      call fvl_startwtimer(52,"solution_afun")

      ! momentum equation
      call ux_globalmat_comp1%blockmul3(x(1:u_dofs),res(1:u_dofs))
      call ux_globalmat_comp2%addrowblockmul3(x(1:u_dofs),res(1:u_dofs))
      call p_globalmat_comp%addcolblockmul3(x(u_dofs+1:u_dofs+p_dofs),res(1:u_dofs))

      ! continuity equation
      call u_globalmat_comp%rowblockmul3(x(1:u_dofs),res(u_dofs+1:u_dofs+p_dofs))
      call fvl_omp_add(res(u_dofs+1:u_dofs+p_dofs),-x(p_dofs+u_dofs+1))

      ! compatibility equation
      res(p_dofs+u_dofs+1)=fvl_omp_dot(x(u_dofs+1:u_dofs+p_dofs),mesh%getcellvolumes())

      ! stop timer
      call fvl_stopwtimer(52)
end subroutine afun

subroutine bfun(res)
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(inout)::res

      ! start timer
      call fvl_startwtimer(53,"solution_afun")

      ! momentum equation
      call fvl_omp_add(reshape(ux_globalrhs,[u_dofs]),reshape(ux_globalsource,[u_dofs]),res(1:u_dofs))

      ! continuity equation
      call fvl_omp_add(p_globalrhs,p_globalsource,res(u_dofs+1:p_dofs+u_dofs))

      ! compatibility equation
      res(p_dofs+u_dofs+1)=0.0d0

      ! stop timer
      call fvl_stopwtimer(53)
end subroutine bfun

subroutine pfun(x,res)
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:size(x,1)),intent(out)::res

      ! start timer
      call fvl_startwtimer(54,"solution_pfun")

      ! momentum equation
      if(ux_precondmat%getmaxnumrowvalues()/=0) then
            call ux_precondmat%blockmul3(x(1:u_dofs),res(1:u_dofs))
      else
            call fvl_omp_assign(x(1:u_dofs),res(1:u_dofs))
      end if

      ! continuity equation
      if(p_precondmat%getmaxnumrowvalues()/=0) then
            call p_precondmat%vecmul(x(u_dofs+1:p_dofs+u_dofs),res(u_dofs+1:p_dofs+u_dofs))
      else
            call fvl_omp_assign(x(u_dofs+1:p_dofs+u_dofs),res(u_dofs+1:p_dofs+u_dofs))
      end if

      ! compatibility equation
      res(p_dofs+u_dofs+1)=x(p_dofs+u_dofs+1)

      ! stop timer
      call fvl_stopwtimer(54)
end subroutine pfun

subroutine updatefun(updatefields,updatetime,writefields)
      logical(kind=__fvl_logical_kind__),intent(in)::updatefields
      logical(kind=__fvl_logical_kind__),intent(in)::updatetime
      logical(kind=__fvl_logical_kind__),intent(in)::writefields
      logical(kind=__fvl_logical_kind__),save::restart=.false.

      real(kind=__fvl_real_kind__),dimension(1:p_dofs+u_dofs+1)::auxv
      real(kind=__fvl_real_kind__),dimension(1:p_dofs+u_dofs+1)::auxvv

      ! start timer
      call fvl_startwtimer(55,"solution_updatefun")

      ! update fields
      if(updatefields) then

            ! update velocity field
            if(staticvmodel/=1) then

                  ! interpolate velocity field
                  call ux_interpolate%evalvalues(u_innerfield,v_innerfield)

                  ! update velocity divergence scheme
                  call ux_divscheme%update()

            end if

            ! first iteration
            if(.not. restart) then

                  ! clean memory
                  call fvl_omp_clean(ux_staticrhs)
                  call fvl_omp_clean(p_staticrhs)

                  ! momentum source mean-values
                  if(staticfupmodel==1) then
                        call fup_model%evalcellmeanvalues(ux_staticsource)
                  end if

                  ! continuity source mean-values
                  if(staticgupmodel==1) then
                        call gup_model%evalcellmeanvalues(p_staticsource)
                  end if

                  ! evaluate pressure gradient right-hand side
                  if(staticrhomodel==1) then
                        call p_gradscheme%evaladdrhs(ux_staticrhs)
                  end if

                  ! evaluate velocity laplacian right-hand side
                  if(staticnumodel==1) then
                        call ux_lapscheme%evaladdrhs(ux_staticrhs)
                  end if

                  ! evaluate velocity divergence right-hand side
                  if(staticvmodel==1) then
                        call ux_divscheme%evaladdrhs(ux_staticrhs)
                  end if

                  ! evaluate velocity divergence right-hand side
                  call u_divscheme%evaladdrhs(p_staticrhs)

            end if

            ! clean memory
            call fvl_omp_clean(ux_globalsource)
            call fvl_omp_clean(p_globalsource)
            call fvl_omp_clean(ux_globalrhs)
            call fvl_omp_clean(p_globalrhs)

            ! momentum source mean-values
            if(staticfupmodel/=1) then
                  call fup_model%evalcellmeanvalues(ux_globalsource)
            end if

            ! continuity source mean-values
            if(staticgupmodel/=1) then
                  call gup_model%evalcellmeanvalues(p_globalsource)
            end if

            ! evaluate velocity laplacian right-hand side
            if(staticnumodel/=1) then
                  call ux_lapscheme%evaladdrhs(ux_globalrhs)
            end if

            ! evaluate velocity divergence right-hand side
            if(staticvmodel/=1) then
                  call ux_divscheme%evaladdrhs(ux_globalrhs)
            end if

            ! evaluate pressure gradient right-hand side
            if(staticrhomodel/=1) then
                  call p_gradscheme%evaladdrhs(ux_globalrhs)
            end if

            ! evaluate velocity divergence right-hand side
            ! call u_divscheme%evaladdrhs(p_globalrhs)

            ! evaluate pressure gradient coefficients matrix
            if(staticrhomodel/=1) then
                  call p_gradscheme%evaladdmatrix(p_globalmat)
            end if

            ! evaluate velocity laplacian coefficients matrix
            if(staticnumodel/=1) then
                  call ux_lapscheme%evalmatrix(ux_globalmat1,ux_globalmat2)
            end if

            ! evaluate velocity divergence coefficients matrix
            if(staticvmodel/=1) then
                  call ux_divscheme%evaladdmatrix(ux_globalmat1,ux_globalmat2)
            end if

            ! add source terms
            call fvl_omp_add(ux_globalsource,ux_staticsource)
            call fvl_omp_add(p_globalsource,p_staticsource)

            ! add right-hand sides
            call fvl_omp_add(ux_globalrhs,ux_staticrhs)
            call fvl_omp_add(p_globalrhs,p_staticrhs)

            ! add coefficients matrices
            if(staticrhomodel/=1) then
                  call p_globalmat_comp%compress(p_globalmat)
                  call fvl_delete(p_globalmat)
            end if

            ! add coefficients matrices
            if(staticrhomodel/=1 .or. staticvmodel/=1) then
                  call ux_globalmat1%addmatrix(ux_staticmat_comp1)
                  call ux_globalmat_comp1%compress(ux_globalmat1)
                  call fvl_delete(ux_globalmat1)
                  call ux_globalmat2%addmatrix(ux_staticmat_comp2)
                  call ux_globalmat_comp2%compress(ux_globalmat2)
                  call fvl_delete(ux_globalmat2)
            end if

            ! compute preconditioner
            if(.not. restart .or. staticnumodel/=1 .or. staticrhomodel/=1 .or. staticvmodel/=1) then
                  if(preconditioner%updatecounter()) then
                        call preconditioner%evalprecond(ux_globalmat_comp1,p_globalmat_comp,u_globalmat_comp,ux_precondmat,p_precondmat)
                  end if
            end if

            ! save restart
            restart=.true.


            !
            ! auxv=0.0d0
            ! auxv(1:u_dofs)=1.0d0
            !
            ! ! momentum equation
            ! call ux_globalmat_comp1%blockmul3(auxv(1:u_dofs),auxvv(1:u_dofs))
            ! call ux_globalmat_comp2%addrowblockmul3(auxv(1:u_dofs),auxvv(1:u_dofs))
            !
            ! call p_globalmat_comp%addcolblockmul3(auxv(u_dofs+1:u_dofs+p_dofs),auxvv(1:u_dofs))
            !
            ! ! continuity equation
            ! call u_globalmat_comp%rowblockmul3(auxv(1:u_dofs),auxvv(u_dofs+1:u_dofs+p_dofs))
            !
            !
            ! ! momentum equation
            ! call fvl_omp_add(ux_globalrhs,ux_globalsource,auxv(1:u_dofs))
            !
            ! ! continuity equation
            ! call fvl_omp_add(p_globalrhs,p_globalsource,auxv(u_dofs+1:p_dofs+u_dofs))
            !
            !
            !
            ! write(*,*) auxv(1:u_dofs/3)-auxvv(1:u_dofs/3)
            ! write(*,*) " "
            ! write(*,*) auxv(1:2*u_dofs/3)-auxvv(1:2*u_dofs/3)
            ! write(*,*) " "
            ! write(*,*) auxv(1:u_dofs)-auxvv(1:u_dofs)
            ! write(*,*) " "
            ! write(*,*) auxv(u_dofs+1:p_dofs+u_dofs)-auxvv(u_dofs+1:p_dofs+u_dofs)
            ! write(*,*) " "
            ! !
            ! ! STOP


      end if

      ! update time
      if(updatetime) then

            ! reset restart
            restart=.false.

      end if

      ! write fields
      if(writefields) then

#include "fvl_writefields_src.f"

      end if

      ! stop timer
      call fvl_stopwtimer(55)
end subroutine updatefun

subroutine fvl_optimize_memory1()
      ! call fvl_deallocate(mesh%celltypes)
      ! call fvl_deallocate(mesh%vertcodes)
      ! ! call fvl_deallocate(mesh%edgecodes)
      ! ! call fvl_deallocate(mesh%cellcodes)
      ! call fvl_deallocate(mesh%boundverts)
      ! call fvl_deallocate(mesh%boundedges)
      ! call fvl_deallocate(mesh%boundcells)
      ! call fvl_deallocate(mesh%verttoedgemaps)
      ! call fvl_deallocate(mesh%verttoedgesizes)
      ! call fvl_deallocate(mesh%verttocellmaps)
      ! call fvl_deallocate(mesh%verttocellsizes)
      ! call fvl_deallocate(mesh%edgetovertmaps)
      ! call fvl_deallocate(mesh%edgetovertsizes)
      ! ! call fvl_deallocate(mesh%edgetocellmaps)
      ! ! call fvl_deallocate(mesh%edgetocellsizes)
      ! call fvl_deallocate(mesh%celltovertmaps)
      ! call fvl_deallocate(mesh%celltovertsizes)
      ! ! call fvl_deallocate(mesh%celltoedgemaps)
      ! ! call fvl_deallocate(mesh%celltoedgesizes)
      ! ! call fvl_deallocate(mesh%celltocellmaps)
      ! ! call fvl_deallocate(mesh%celltocellsizes)
      ! call fvl_deallocate(mesh%primcelltoedgemaps)
      ! call fvl_deallocate(mesh%primcelltoedgesizes)
      ! call fvl_deallocate(mesh%edgetoprimcellmaps)
      ! call fvl_deallocate(mesh%edgetoprimcellsizes)
      ! call fvl_deallocate(mesh%primedgetocellmaps)
      ! call fvl_deallocate(mesh%primedgetocellsizes)
      ! call fvl_deallocate(mesh%celltoprimedgemaps)
      ! call fvl_deallocate(mesh%celltoprimedgesizes)
      ! call fvl_deallocate(mesh%vertindices)
      ! call fvl_deallocate(mesh%edgeindices)
      ! call fvl_deallocate(mesh%cellindices)
      ! call fvl_deallocate(mesh%vertcentres)
      ! ! call fvl_deallocate(mesh%edgecentres)
      ! ! call fvl_deallocate(mesh%cellcentres)
      ! ! call fvl_deallocate(mesh%edgenormals)
      ! ! call fvl_deallocate(mesh%edgevolumes)
      ! ! call fvl_deallocate(mesh%cellvolumes)
      ! ! call fvl_deallocate(mesh%numedgequadratpoints)
      ! ! call fvl_deallocate(mesh%edgequadratpoints)
      ! ! call fvl_deallocate(mesh%edgequadratweights)
      ! ! call fvl_deallocate(mesh%numcellquadratpoints)
      ! ! call fvl_deallocate(mesh%cellquadratpoints)
      ! ! call fvl_deallocate(mesh%cellquadratweights)
      ! ! call fvl_deallocate(mesh%numedgeboundpoints)
      ! ! call fvl_deallocate(mesh%edgeboundpoints)
      ! ! call fvl_deallocate(mesh%edgeboundnormals)
      ! ! call fvl_deallocate(mesh%edgeboundtangents)
      ! ! call fvl_deallocate(mesh%edgeboundbitangents)
end subroutine fvl_optimize_memory1

subroutine fvl_optimize_memory2()
      ! if(staticalphamodel==1 .and. statictmodel==1 .and. statictdmodel==1) then
      !       call fvl_delete(u_lapscheme)
      ! else if(staticalphamodel==1) then
      !       ! call fvl_deallocate(u_lapscheme%celltoedgemaps)
      !       call fvl_deallocate(u_lapscheme%edgesalias)
      !       ! call fvl_deallocate(u_lapscheme%edgecoeffields)
      !       call fvl_deallocate(u_lapscheme%edgedegrees)
      !       call fvl_deallocate(u_lapscheme%edgefluxcoefs)
      !       ! call fvl_deallocate(u_lapscheme%edgestencils)
      !       call fvl_deallocate(u_lapscheme%inneredgefastrecs)
      !       ! call fvl_deallocate(u_lapscheme%boundedgefastcrecs)
      !       call fvl_deallocate(u_lapscheme%boundedgefastrecs)
      !       ! call fvl_deallocate(u_lapscheme%boundfastcrecs)
      ! end if
end subroutine fvl_optimize_memory2

subroutine fvl_optimize_memory3()
      ! if(staticpsimodel==1 .and. statictmodel==1 .and. statictdmodel==1) then
      !       call fvl_delete(u_divscheme)
      ! else if(staticpsimodel==1) then
      !       ! call fvl_deallocate(u_divscheme%celltoedgemaps)
      !       call fvl_deallocate(u_divscheme%edgesalias)
      !       ! call fvl_deallocate(u_divscheme%edgecoeffields)
      !       ! call fvl_deallocate(u_divscheme%upwindcells)
      !       call fvl_deallocate(u_divscheme%boundedgedegrees)
      !       call fvl_deallocate(u_divscheme%celldegrees)
      !       ! call fvl_deallocate(u_divscheme%edgefluxcoefs)
      !       ! call fvl_deallocate(u_divscheme%boundedgestencils)
      !       call fvl_deallocate(u_divscheme%cellstencils)
      !       ! call fvl_deallocate(u_divscheme%edgefastcrecs)
      !       call fvl_deallocate(u_divscheme%edgefastrecs)
      !       ! call fvl_deallocate(u_divscheme%boundfastcrecs)
      ! end if
end subroutine fvl_optimize_memory3

subroutine fvl_optimize_memory4()
      ! call fvl_delete(p_lapscheme)
      ! call fvl_delete(u_integscheme)
      ! call fvl_delete(u_lapscheme)
      ! call fvl_delete(u_divscheme)

      ! call fvl_delete(u_integmat)
      ! call fvl_delete(u_staticmat)

      ! call fvl_deallocate(u_globalmat)
      ! call fvl_delete(ux_precondmat)
      ! call fvl_delete(u_integmat_comp)
      ! call fvl_delete(u_staticmat_comp)
      ! call fvl_deallocate(u_globalmatcomp)
end subroutine fvl_optimize_memory4

end program fvl_pro_thermoincompflow3d_main
! end of file
