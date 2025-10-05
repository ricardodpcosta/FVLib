! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: 2D non-isothermal Newtonian incompressible fluid flow solver
! Modification: February, 2025

#include "fvl_psi_pseudo_boundcond2d_mod.f"
#include "fvl_omega_pseudo_boundcond2d_mod.f"

#include "macros.f"

program fvl_pro_thermopsiomega2d_main

use fvl_lib2d
use fvl_psi_pseudo_boundcond2d_mod
use fvl_omega_pseudo_boundcond2d_mod

implicit none

! ============================================================================
! DECLARE VARIABLES
! ============================================================================

! ----------------------------------------------------------------------------
! variables
! ----------------------------------------------------------------------------

integer(kind=__fvl_integer_kind__)::i,psi_dofs,omega_dofs,total_dofs,&
      staticmumodel,staticrhomodel,staticnumodel,staticpsimodel,staticomegamodel,staticumodel,staticfupmodel,&
      staticfpsiomegamodel,staticgpsiomegamodel
real(kind=__fvl_real_kind__),allocatable,dimension(:)::psi_staticsource,psi_globalsource
real(kind=__fvl_real_kind__),allocatable,dimension(:)::omega_staticsource,omega_globalsource
real(kind=__fvl_real_kind__),allocatable,dimension(:)::psi_staticrhs,psi_globalrhs
real(kind=__fvl_real_kind__),allocatable,dimension(:)::omega_staticrhs,omega_globalrhs,xi_globalrhs
type(fvl_calc_scalarfield2d)::psi_innerfield,omega_innerfield,omegaconv_innerfield,omegadiff_innerfield,phi_innerfield
type(fvl_calc_vectorfield2d)::u_innerfield,normal_innerfield
type(fvl_calc_scalarfield2d)::psinn_boundfield
type(fvl_centredpro_interpolate2d)::psi_interpolate
type(fvl_centredpro_scalarboundscheme2d)::psi_boundscheme
type(fvl_centredpro_laplacianscheme2d)::psi_lapscheme
type(fvl_centredpro_integratescheme2d)::omega_integscheme
type(fvl_centredpro_scalarinterscheme2d)::omega_interscheme1,omega_interscheme2
type(fvl_centredpro_laplacianscheme2d)::omega_lapscheme
type(fvl_upwindpro_divergencescheme2d)::omega_divscheme
type(fvl_lil_rspmat)::omega_staticmat,omega_integmat,omega_globalmat,psi_globalmat,psi_precondmat,omega_precondmat
type(fvl_csr_rspmat)::omega_staticmat_comp,psi_globalmat_comp,omega_integmat_comp,omega_globalmat_comp
type(fvl_linearsolver)::linearsolver
type(fvl_basic_preconditioner)::preconditioner
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

#include "fvl_createmesh2d_hdr.f"

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
#define __fvl_field__ cp
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_scalarfield_models2d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ k
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_scalarfield_models2d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ mu
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_scalarfield_models2d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ rho
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_scalarfield_models2d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ t
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_scalarfield_models2d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ omega
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_scalarfield_models2d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ psi
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_scalarfield_models2d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ u
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_vectorfield_models2d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ fup
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_vectorfield_models2d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ fpsiomega
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_scalarfield_models2d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ gpsiomega
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_scalarfield_models2d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ alpha
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_scalarfield_models2d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ nu
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_scalarfield_models2d_hdr.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ phi
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_scalarfield_models2d_hdr.f"

! ----------------------------------------------------------------------------
! specific model functions
! ----------------------------------------------------------------------------

#include "fvl_omega_models2d_hdr.f"
#include "fvl_psi_models2d_hdr.f"

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
      call fvl_loginfo("About: 2D non-isothermal non-Newtonian incompressible fluid flow solver.")
      call fvl_loginfo("Usage: fvl-pro-thermopsiomega2d [options]")
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
staticpsimodel          = controldictfile%getvalue("solver","static_psi_model",0)
staticomegamodel        = controldictfile%getvalue("solver","static_omega_model",0)
staticumodel            = controldictfile%getvalue("solver","static_u_model",0)
staticfupmodel          = controldictfile%getvalue("solver","static_fup_model",0)
staticfpsiomegamodel    = controldictfile%getvalue("solver","static_fpsiomega_model",0)
staticgpsiomegamodel    = controldictfile%getvalue("solver","static_gpsiomega_model",0)
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

#include "fvl_createmesh2d_src.f"

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
psi_dofs                = mesh%getnumcells()
omega_dofs              = mesh%getnumcells()
total_dofs              = psi_dofs+omega_dofs

! allocate memory
call fvl_allocate(psi_staticsource,psi_dofs)
call fvl_allocate(psi_globalsource,psi_dofs)
call fvl_allocate(omega_staticsource,omega_dofs)
call fvl_allocate(omega_globalsource,omega_dofs)
call fvl_allocate(psi_staticrhs,psi_dofs)
call fvl_allocate(psi_globalrhs,psi_dofs)
call fvl_allocate(omega_staticrhs,omega_dofs)
call fvl_allocate(omega_globalrhs,omega_dofs)
call fvl_allocate(xi_globalrhs,100)

! clean memory
call fvl_omp_clean(psi_staticsource)
call fvl_omp_clean(omega_staticsource)
call fvl_omp_clean(psi_globalsource)
call fvl_omp_clean(omega_globalsource)
call fvl_omp_clean(psi_staticrhs)
call fvl_omp_clean(omega_staticrhs)
call fvl_omp_clean(psi_globalrhs)
call fvl_omp_clean(omega_globalrhs)
call fvl_clean(psi_globalmat)
call fvl_clean(omega_staticmat)
call fvl_clean(omega_globalmat)
call fvl_clean(psi_globalmat_comp)
call fvl_clean(omega_staticmat_comp)
call fvl_clean(omega_integmat_comp)
call fvl_clean(omega_globalmat_comp)
call fvl_clean(psi_precondmat)
call fvl_clean(omega_precondmat)

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
mu_model                = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_edgequadratsfield,&
                              dictfile=modelsdictfile,dictlabel="mu")

! density model
rho_model               = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_edgequadratsfield,&
                              dictfile=modelsdictfile,dictlabel="rho")

! streamfunction model
psi_model               = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_cellmeansfield,&
                              dictfile=modelsdictfile,dictlabel="psi")

! vorticity model
omega_model             = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_cellmeansfield,&
                              dictfile=modelsdictfile,dictlabel="omega")

! velocity model
u_model                 = fvl_vectormodel2d(mesh=mesh,fieldtype=fvl_field2d_edgequadratsfield,&
                              dictfile=modelsdictfile,dictlabel="u")

! forces model
fup_model               = fvl_vectormodel2d(mesh=mesh,fieldtype=fvl_field2d_edgemeansfield,&
                              dictfile=modelsdictfile,dictlabel="fup")

! forces model
fpsiomega_model         = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_cellmeansfield,&
                              dictfile=modelsdictfile,dictlabel="fpsiomega")

! mass source model
gpsiomega_model         = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_cellmeansfield,&
                              dictfile=modelsdictfile,dictlabel="gpsiomega")

! normal velocity model
phi_model               = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_edgequadratsfield,&
                              dictfile=modelsdictfile,dictlabel="phi")

! ----------------------------------------------------------------------------
! user-defined models
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Setting user-defined models")

#include "model.f"

! initialize streamfunction pseudo patches
call psi_initialize_pseudo_patches()

! initialize vorticity pseudo patches
call omega_initialize_pseudo_patches()

! ----------------------------------------------------------------------------
! computed fields
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Initializing computed fields")

! kinematic viscosity model
nu_model                = -mu_model/rho_model

! streamfunction field
psi_innerfield          = fvl_calc_scalarfield2d(patch=domainpatch,type=fvl_field2d_cellmeansfield,&
                              label="psi",filepath=solutionfilesdir,fileform=solutionfilesform)

! vorticity field
omega_innerfield        = fvl_calc_scalarfield2d(patch=domainpatch,type=fvl_field2d_cellmeansfield,&
                              label="omega",filepath=solutionfilesdir,fileform=solutionfilesform)

! velocity field
u_innerfield            = fvl_calc_vectorfield2d(patch=domainpatch,type=fvl_field2d_edgequadratsfield,&
                              label="u",filepath=solutionfilesdir,fileform=solutionfilesform)

! psinn field
psinn_boundfield        = fvl_calc_scalarfield2d(patch=omega_pseudo_patch,type=fvl_field2d_edgeboundsfield,&
                              label="psinn",filepath=solutionfilesdir,fileform=solutionfilesform)

! omegaconv field
omegaconv_innerfield     = fvl_calc_scalarfield2d(patch=psi_pseudo_patch,type=fvl_field2d_edgemeansfield,&
                              label="omegaconv",filepath=solutionfilesdir,fileform=solutionfilesform)

! omegadiff field
omegadiff_innerfield     = fvl_calc_scalarfield2d(patch=psi_pseudo_patch,type=fvl_field2d_edgemeansfield,&
                              label="omegadiff",filepath=solutionfilesdir,fileform=solutionfilesform)

! normal field
normal_innerfield       = fvl_calc_vectorfield2d(patch=psi_pseudo_patch,type=fvl_field2d_edgequadratsfield,&
                              label="normal",filepath=solutionfilesdir,fileform=solutionfilesform)

! phi field
phi_innerfield          = fvl_calc_scalarfield2d(patch=psi_pseudo_patch,type=fvl_field2d_edgequadratsfield,&
                              label="phi",filepath=solutionfilesdir,fileform=solutionfilesform)

! set fields
call solutionfields(1)  % setptr(omega_innerfield)
call solutionfields(2)  % setptr(psi_innerfield)
do i=1,psi_numpseudoboundconds
      call solutionfields(i+2)      % setptr(psi_singlefields(i))
end do

! intialize streamfunction field
call psi_model          % evalcellmeanvalues(psi_innerfield)

! intialize vorticity field
call omega_model        % evalcellmeanvalues(omega_innerfield)

! velocity model
call u_model            % setfield(u_innerfield)

! phi model
call phi_model          % setfield(phi_innerfield)

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
! streamfunction interpolation scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Computing interpolate(psi) scheme")

! start timer
call fvl_startwtimer(41,"interpolate(psi)_scheme")

! initialize scheme
psi_interpolate         = fvl_centredpro_interpolate2d(patch=domainpatch,mesh=mesh,&
                              type=fvl_centredpro_interpolate2d_edgequadratgrads,&
                              dictfile=schemesdictfile,dictlabel="interpolate(psi)")

if(omega_numpseudoboundconds>0) then

! initialize scheme
psi_boundscheme         = fvl_centredpro_scalarboundscheme2d(condsmodel=psi_model,mesh=mesh,&
                              patch=omega_pseudo_patch,evaltype=fvl_centredpro_scalarboundscheme2d_pointderivnns,&
                              dictfile=schemesdictfile,dictlabel="interpolate(psi)",&
                              staticcoeffields=.true.,staticboundconds=.true.)

end if

! deallocate memory
! call fvl_delete(psi_interpolate)
! call fvl_delete(psi_boundscheme)

! stop timer
call fvl_stopwtimer(41)

! ----------------------------------------------------------------------------
! streamfunction laplacian scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Computing laplacian(psi) scheme")

! start timer
call fvl_startwtimer(42,"laplacian(psi)_scheme")

! initialize scheme
psi_lapscheme           = fvl_centredpro_laplacianscheme2d(condsmodel=psi_model,&
                              mesh=mesh,dictfile=schemesdictfile,dictlabel="laplacian(psi)",&
                              staticcoeffields=.true.,staticboundconds=.true.)

! evaluate coefficients matrix
call psi_lapscheme      % evalmatrix(psi_globalmat)

! ! optimize memory
! call psi_lapscheme%deallocatelhs()

! compress coefficients matrix
call psi_globalmat_comp%compress(psi_globalmat)
call fvl_delete(psi_globalmat)

! stop timer
call fvl_stopwtimer(42)

! ----------------------------------------------------------------------------
! vorticity integrate scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Computing integrate(omega) scheme")

! start timer
call fvl_startwtimer(43,"integrate(omega)_scheme")

! initialize scheme
omega_integscheme       = fvl_centredpro_integratescheme2d(mesh=mesh,&
                              dictfile=schemesdictfile,dictlabel="integrate(omega)",&
                              staticcoeffields=.true.)

! evaluate coefficients matrix
call omega_integscheme  % evalmatrix(omega_integmat)

! compress coefficients matrix
call omega_integmat_comp      % compress(omega_integmat)

! deallocate memory
call fvl_delete(omega_integscheme)
call fvl_delete(omega_integmat)

! stop timer
call fvl_stopwtimer(43)

! ----------------------------------------------------------------------------
! vorticity interpolation scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Computing interpolate(omega) scheme")

! start timer
call fvl_startwtimer(44,"interpolate(omega)_scheme")

if(psi_numpseudoboundconds>0) then

! initialize scheme
omega_interscheme1      = fvl_centredpro_scalarinterscheme2d(coefsmodel=phi_model,condsmodel=omega_model,mesh=mesh,&
                              patch=psi_pseudo_patch,evaltype=fvl_centredpro_scalarinterscheme2d_meanvalues,&
                              dictfile=schemesdictfile,dictlabel="interpolate(omega)",&
                              staticcoeffields=.true.,staticboundconds=.true.)

! initialize scheme
omega_interscheme2      = fvl_centredpro_scalarinterscheme2d(coefsmodel=nu_model,condsmodel=omega_model,mesh=mesh,&
                              patch=psi_pseudo_patch,evaltype=fvl_centredpro_scalarinterscheme2d_meanderivns,&
                              dictfile=schemesdictfile,dictlabel="interpolate(omega)",&
                              staticcoeffields=.true.,staticboundconds=.true.)

end if

! stop timer
call fvl_stopwtimer(44)

! ----------------------------------------------------------------------------
! vorticity laplacian scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Computing laplacian(nu,omega) scheme")

! start timer
call fvl_startwtimer(45,"laplacian(nu,omega)_scheme")

! initialize scheme
omega_lapscheme         = fvl_centredpro_laplacianscheme2d(coefsmodel=nu_model,condsmodel=omega_model,&
                              mesh=mesh,dictfile=schemesdictfile,dictlabel="laplacian(nu,omega)",&
                              staticcoeffields=(staticnumodel==1),staticboundconds=.true.)

! evaluate coefficients matrix
if(staticnumodel==1) then
      call omega_lapscheme    % evalmatrix(omega_staticmat)
end if

! ! optimize memory
! call omega_lapscheme%deallocatelhs()

! stop timer
call fvl_stopwtimer(45)

! ----------------------------------------------------------------------------
! vorticity divergence scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Computing divergence(u,omega) scheme")

! start timer
call fvl_startwtimer(46,"divergence(u,omega)_scheme")

! initialize scheme
omega_divscheme         = fvl_upwindpro_divergencescheme2d(coefsmodel=u_model,condsmodel=omega_model,&
                              mesh=mesh,dictfile=schemesdictfile,dictlabel="divergence(u,omega)",&
                              staticcoeffields=(staticumodel==1),staticboundconds=.true.)

! evaluate coefficients matrix
if(staticumodel==1) then
      call omega_divscheme    % evaladdmatrix(omega_staticmat)
end if

! ! optimize memory
! call omega_divscheme%deallocatelhs()

! compress coefficients matrix
if(staticnumodel==1 .or. staticumodel==1) then
      call omega_staticmat_comp%compress(omega_staticmat)
      call fvl_delete(omega_staticmat)
end if

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
preconditioner          = fvl_basic_preconditioner(dictfile=solutiondictfile,dictlabel="preconditioner")

! fixpoint solver
fixpointsolver          = fvl_fixpointsolver(dictfile=solutiondictfile,dictlabel="fixpoint_solver")

! time solver
timesolver              = fvl_basic_timesolver(runtime=mesh,fields=solutionfields(1:psi_numpseudoboundconds+2),&
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
! compute velocity field
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo("  -> Computing velocity field")

! velocity field
u_innerfield            = fvl_calc_vectorfield2d(patch=domainpatch,type=fvl_field2d_edgemeansfield,&
                              label="u",filepath=solutionfilesdir,fileform=solutionfilesform)

! streamfunction interpolation scheme
psi_interpolate         = fvl_centredpro_interpolate2d(patch=domainpatch,mesh=mesh,&
                              type=fvl_centredpro_interpolate2d_edgemeansgrads,&
                              dictfile=schemesdictfile,dictlabel="interpolate(psi)")

! loop over time directories
do i=1,numtimedirectories

      ! set step and time
      call mesh               % setcurrentstep(steps(i))
      call mesh               % setcurrenttime(times(i))

      ! read streamfunction field
      call psi_innerfield     % read()

      ! interpolate streamfunction field
      call psi_interpolate    % evalgrads(psi_innerfield,u_innerfield)

      ! compute velocity field
      call u_innerfield       % perpfield()

      ! write fields
      call u_innerfield       % write()

end do

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
#define __fvl_field__ cp
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_scalarfield_models2d_src.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ k
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_scalarfield_models2d_src.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ mu
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_scalarfield_models2d_src.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ rho
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_scalarfield_models2d_src.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ t
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_scalarfield_models2d_src.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ omega
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_scalarfield_models2d_src.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ psi
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_scalarfield_models2d_src.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ u
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_vectorfield_models2d_src.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ fup
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_vectorfield_models2d_src.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ fpsiomega
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_scalarfield_models2d_src.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ gpsiomega
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_scalarfield_models2d_src.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ alpha
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_scalarfield_models2d_src.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ nu
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_scalarfield_models2d_src.f"

#undef __fvl_field__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_field__ phi
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_scalarfield_models2d_src.f"

! ----------------------------------------------------------------------------
! specific model functions
! ----------------------------------------------------------------------------

#include "fvl_omega_models2d_src.f"
#include "fvl_psi_models2d_src.f"

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

      ! vorticity equation
      if(omega_integmat_comp%getmaxnumrowvalues()/=0) then
            call omega_integmat_comp%vecmul(x(1:omega_dofs),res(1:omega_dofs))
      else
            call fvl_omp_assign(x(1:omega_dofs),res(1:omega_dofs))
      end if

      ! streamfuntion equation
      call fvl_omp_assign(0.0d0,res(omega_dofs+1:total_dofs))

      ! compatibility equation
      if(psi_numpseudoboundconds>0) then
            ! compute left-hand side
            do i=1,psi_numpseudoboundconds
                  res(total_dofs+i)=0.0d0
            end do
      end if

      ! stop timer
      call fvl_stopwtimer(51)
end subroutine tfun

subroutine afun(x,res)
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:size(x,1)),intent(out)::res

      ! start timer
      call fvl_startwtimer(52,"solution_afun")

      ! streamfunction pseudo boundary condition
      if(psi_numpseudoboundconds>0) then

            ! update boundary conditions
            do i=1,psi_numpseudoboundconds
                  call psi_singlefields(i)%setvalue(x(total_dofs+i))
                  call psi_pseudoboundconds(i)%updatefields()
            end do

      end if

      ! vorticity pseudo boundary condition
      if(omega_numpseudoboundconds>0) then

            ! interpolate velocity field
            call psi_boundscheme%evalboundpointderivnns(x(omega_dofs+1:total_dofs),psinn_boundfield)

            ! update boundary conditions
            do i=1,omega_numpseudoboundconds
                  call omega_pseudoboundconds(i)%updatefields()
            end do

      end if

      ! vorticity equation
      if(omega_numpseudoboundconds>0) then

            ! compute left-hand side
            call omega_lapscheme%evalrhs(res(1:omega_dofs))
            call omega_divscheme%evaladdrhs(res(1:omega_dofs))
            call fvl_omp_sub(omega_globalrhs,res(1:omega_dofs),res(1:omega_dofs))
            call omega_globalmat_comp%addvecmul(x(1:omega_dofs),res(1:omega_dofs))

      else

            ! compute left-hand side
            call omega_globalmat_comp%vecmul(x(1:omega_dofs),res(1:omega_dofs))

      end if

      ! streamfunction equation
      if(psi_numpseudoboundconds>0) then

            ! compute left-hand side
            call psi_lapscheme%evalrhs(res(omega_dofs+1:total_dofs))
            call fvl_omp_sub(psi_globalrhs,res(omega_dofs+1:total_dofs),res(omega_dofs+1:total_dofs))
            if(omega_integmat_comp%getmaxnumrowvalues()/=0) then
                  call omega_integmat_comp%addvecmul(x(1:omega_dofs),res(omega_dofs+1:total_dofs))
            else
                  call fvl_omp_add(res(omega_dofs+1:total_dofs),x(1:omega_dofs))
            end if
            call psi_globalmat_comp%addvecmul(x(omega_dofs+1:total_dofs),res(omega_dofs+1:total_dofs))

      else

            ! compute left-hand side
            if(omega_integmat_comp%getmaxnumrowvalues()/=0) then
                  call omega_integmat_comp%vecmul(x(1:omega_dofs),res(omega_dofs+1:total_dofs))
            else
                  call fvl_omp_assign(x(1:omega_dofs),res(omega_dofs+1:total_dofs))
            end if
            call psi_globalmat_comp%addvecmul(x(omega_dofs+1:total_dofs),res(omega_dofs+1:total_dofs))

      end if

      ! compatibility equation
      if(psi_numpseudoboundconds>0) then

            ! interpolate vorticity field
            call omega_interscheme1%evaledgemeanvalues(x(1:omega_dofs),omegaconv_innerfield)
            call omega_interscheme2%evaledgemeanderivns(x(1:omega_dofs),omegadiff_innerfield)

            ! compute left-hand side
            do i=1,psi_numpseudoboundconds
                  res(total_dofs+i)=psi_pseudoboundconds(i)%evalcirculation()+xi_globalrhs(i)
            end do

      end if

      ! stop timer
      call fvl_stopwtimer(52)
end subroutine afun

subroutine bfun(res)
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(inout)::res

      ! start timer
      call fvl_startwtimer(53,"solution_afun")

      ! vorticity equation
      call fvl_omp_add(omega_globalrhs,omega_globalsource,res(1:omega_dofs))

      ! streamfunction equation
      call fvl_omp_add(psi_globalrhs,psi_globalsource,res(omega_dofs+1:total_dofs))

      ! compatibility equation
      if(psi_numpseudoboundconds>0) then
            res(total_dofs+1:total_dofs+psi_numpseudoboundconds)=xi_globalrhs(1:psi_numpseudoboundconds)
      end if

      ! stop timer
      call fvl_stopwtimer(53)
end subroutine bfun

subroutine pfun(x,res)
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:size(x,1)),intent(out)::res

      ! start timer
      call fvl_startwtimer(54,"solution_pfun")

      ! vorticity equation
      if(omega_precondmat%getmaxnumrowvalues()/=0) then
            call omega_precondmat%vecmul(x(1:omega_dofs),res(1:omega_dofs))
      else
            call fvl_omp_assign(x(1:omega_dofs),res(1:omega_dofs))
      end if

      ! streamfuntion equation
      if(psi_precondmat%getmaxnumrowvalues()/=0) then
            call psi_precondmat%vecmul(x(omega_dofs+1:total_dofs),res(omega_dofs+1:total_dofs))
      else
            call fvl_omp_assign(x(omega_dofs+1:total_dofs),res(omega_dofs+1:total_dofs))
      end if

      ! compatibility equation
      if(psi_numpseudoboundconds>0) then
            res(total_dofs+1:total_dofs+psi_numpseudoboundconds)=x(total_dofs+1:total_dofs+psi_numpseudoboundconds)
      end if

      ! stop timer
      call fvl_stopwtimer(54)
end subroutine pfun

subroutine updatefun(updatefields,updatetime,writefields)
      logical(kind=__fvl_logical_kind__),intent(in)::updatefields
      logical(kind=__fvl_logical_kind__),intent(in)::updatetime
      logical(kind=__fvl_logical_kind__),intent(in)::writefields
      logical(kind=__fvl_logical_kind__),save::restart=.false.

      ! start timer
      call fvl_startwtimer(55,"solution_updatefun")

      ! update fields
      if(updatefields) then

            ! update boundary conditions
            if(psi_numpseudoboundconds>0) then

                  ! update fields
                  do i=1,psi_numpseudoboundconds
                        call psi_singlefields(i)%setvalue(0.0d0)
                        call psi_pseudoboundconds(i)%updatefields()
                  end do

            end if

            ! update boundary conditions
            if(omega_numpseudoboundconds>0) then

                  ! interpolate streamfunction field
                  call psi_boundscheme%evalboundpointderivnnsrhs(psinn_boundfield)

                  ! update fields
                  do i=1,omega_numpseudoboundconds
                        call omega_pseudoboundconds(i)%updatefields()
                  end do

            end if

            ! first iteration
            if(.not. restart) then

                  ! initialize normal field
                  call normal_innerfield%setedgequadratvalues(mesh%getedgenormals())

            end if

            ! update velocity field
            if(staticumodel/=1) then

                  ! interpolate streamfunction field
                  call psi_interpolate%evalgrads(psi_innerfield,u_innerfield)

                  ! compute velocity field
                  call u_innerfield%perpfield()

                  ! compute normal velocity field
                  call phi_innerfield%dotfield(u_innerfield,normal_innerfield)

                  ! update vorticity divergence scheme
                  call omega_divscheme%update()

            end if

            ! first iteration
            if(.not. restart) then

                  ! clean memory
                  call fvl_omp_clean(omega_staticrhs)
                  call fvl_omp_clean(psi_staticrhs)

                  ! momentum source mean-values
                  if(staticfpsiomegamodel==1) then
                        call fpsiomega_model%evalcellmeanvalues(omega_staticsource)
                  end if

                  ! continuity source mean-values
                  if(staticgpsiomegamodel==1) then
                        call gpsiomega_model%evalcellmeanvalues(psi_staticsource)
                  end if

                  ! evaluate streamfunction laplacian right-hand side
                  call psi_lapscheme%evaladdrhs(psi_staticrhs)

                  ! evaluate vorticity laplacian right-hand side
                  if(staticnumodel==1) then
                        call omega_lapscheme%evaladdrhs(omega_staticrhs)
                  end if

                  ! evaluate vorticity divergence right-hand side
                  if(staticumodel==1) then
                        call omega_divscheme%evaladdrhs(omega_staticrhs)
                  end if

            end if

            ! clean memory
            call fvl_omp_clean(omega_globalsource)
            call fvl_omp_clean(psi_globalsource)
            call fvl_omp_clean(omega_globalrhs)
            call fvl_omp_clean(psi_globalrhs)

            ! momentum source mean-values
            if(staticfpsiomegamodel/=1) then
                  call fpsiomega_model%evalcellmeanvalues(omega_globalsource)
            end if

            ! continuity source mean-values
            if(staticgpsiomegamodel/=1) then
                  call gpsiomega_model%evalcellmeanvalues(psi_globalsource)
            end if

            ! evaluate vorticity laplacian right-hand side
            if(staticnumodel/=1) then
                  call omega_lapscheme%evaladdrhs(omega_globalrhs)
            end if

            ! evaluate vorticity divergence right-hand side
            if(staticumodel/=1) then
                  call omega_divscheme%evaladdrhs(omega_globalrhs)
            end if

            ! evaluate vorticity laplacian coefficients matrix
            if(staticnumodel/=1) then
                  call omega_lapscheme%evalmatrix(omega_globalmat)
            end if

            ! evaluate vorticity divergence coefficients matrix
            if(staticumodel/=1) then
                  call omega_divscheme%evaladdmatrix(omega_globalmat)
            end if

            ! interpolate vorticity field
            if(psi_numpseudoboundconds>0) then

                  ! recompute interpolation scheme
                  call omega_interscheme1%makereconst()

                  ! compute convective circulation field
                  call omega_interscheme1%evaledgemeanvaluesrhs(omegaconv_innerfield)

                  ! compute diffusive circulation field
                  call omega_interscheme2%evaledgemeanderivnsrhs(omegadiff_innerfield)

                  ! compute right-hand side
                  do i=1,psi_numpseudoboundconds
                        xi_globalrhs(i)=-psi_pseudoboundconds(i)%evalcirculation()
                  end do

            end if

            ! add source terms
            call fvl_omp_add(omega_globalsource,omega_staticsource)
            call fvl_omp_add(psi_globalsource,psi_staticsource)

            ! add right-hand sides
            call fvl_omp_add(omega_globalrhs,omega_staticrhs)
            call fvl_omp_add(psi_globalrhs,psi_staticrhs)

            ! add coefficients matrices
            if(staticnumodel/=1 .or. staticumodel/=1) then
                  call omega_globalmat%addmatrix(omega_staticmat_comp)
                  call omega_globalmat_comp%compress(omega_globalmat)
                  call fvl_delete(omega_globalmat)
            end if

            ! first iteration
            if(.not. restart) then

                  ! compute preconditioner
                  call preconditioner%evalprecond(psi_globalmat_comp,psi_precondmat)

                  ! reset counter
                  call preconditioner%resetcounter()

            end if

            ! compute preconditioner
            if(.not. restart .or. staticnumodel/=1 .or. staticumodel/=1) then
                  if(preconditioner%updatecounter()) then
                        call preconditioner%evalprecond(omega_globalmat_comp,omega_precondmat)
                  end if
            end if

            ! save restart
            restart=.true.

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
      call fvl_deallocate(mesh%celltypes)
      call fvl_deallocate(mesh%vertcodes)
      ! call fvl_deallocate(mesh%edgecodes)
      ! call fvl_deallocate(mesh%cellcodes)
      call fvl_deallocate(mesh%boundverts)
      call fvl_deallocate(mesh%boundedges)
      call fvl_deallocate(mesh%boundcells)
      call fvl_deallocate(mesh%verttoedgemaps)
      call fvl_deallocate(mesh%verttoedgesizes)
      call fvl_deallocate(mesh%verttocellmaps)
      call fvl_deallocate(mesh%verttocellsizes)
      call fvl_deallocate(mesh%edgetovertmaps)
      call fvl_deallocate(mesh%edgetovertsizes)
      ! call fvl_deallocate(mesh%edgetocellmaps)
      ! call fvl_deallocate(mesh%edgetocellsizes)
      call fvl_deallocate(mesh%celltovertmaps)
      call fvl_deallocate(mesh%celltovertsizes)
      ! call fvl_deallocate(mesh%celltoedgemaps)
      ! call fvl_deallocate(mesh%celltoedgesizes)
      ! call fvl_deallocate(mesh%celltocellmaps)
      ! call fvl_deallocate(mesh%celltocellsizes)
      call fvl_deallocate(mesh%primcelltoedgemaps)
      call fvl_deallocate(mesh%primcelltoedgesizes)
      call fvl_deallocate(mesh%edgetoprimcellmaps)
      call fvl_deallocate(mesh%edgetoprimcellsizes)
      call fvl_deallocate(mesh%primedgetocellmaps)
      call fvl_deallocate(mesh%primedgetocellsizes)
      call fvl_deallocate(mesh%celltoprimedgemaps)
      call fvl_deallocate(mesh%celltoprimedgesizes)
      call fvl_deallocate(mesh%vertindices)
      call fvl_deallocate(mesh%edgeindices)
      call fvl_deallocate(mesh%cellindices)
      call fvl_deallocate(mesh%vertcentres)
      ! call fvl_deallocate(mesh%edgecentres)
      ! call fvl_deallocate(mesh%cellcentres)
      ! call fvl_deallocate(mesh%edgenormals)
      ! call fvl_deallocate(mesh%edgeareas)
      ! call fvl_deallocate(mesh%cellvolumes)
      ! call fvl_deallocate(mesh%numedgequadratpoints)
      ! call fvl_deallocate(mesh%edgequadratpoints)
      ! call fvl_deallocate(mesh%edgequadratweights)
      ! call fvl_deallocate(mesh%numcellquadratpoints)
      ! call fvl_deallocate(mesh%cellquadratpoints)
      ! call fvl_deallocate(mesh%cellquadratweights)
      ! call fvl_deallocate(mesh%numedgeboundpoints)
      ! call fvl_deallocate(mesh%edgeboundpoints)
      ! call fvl_deallocate(mesh%edgeboundnormals)
      ! call fvl_deallocate(mesh%edgeboundtangents)
      ! call fvl_deallocate(mesh%edgeboundbitangents)
end subroutine fvl_optimize_memory1

subroutine fvl_optimize_memory2()
      ! if(staticalphamodel==1 .and. statictmodel==1 .and. statictdmodel==1) then
      !       call fvl_delete(omega_lapscheme)
      ! else if(staticalphamodel==1) then
      !       ! call fvl_deallocate(omega_lapscheme%celltoedgemaps)
      !       call fvl_deallocate(omega_lapscheme%edgesalias)
      !       ! call fvl_deallocate(omega_lapscheme%edgecoeffields)
      !       call fvl_deallocate(omega_lapscheme%edgedegrees)
      !       call fvl_deallocate(omega_lapscheme%edgefluxcoefs)
      !       ! call fvl_deallocate(omega_lapscheme%edgestencils)
      !       call fvl_deallocate(omega_lapscheme%inneredgefastrecs)
      !       ! call fvl_deallocate(omega_lapscheme%boundedgefastcrecs)
      !       call fvl_deallocate(omega_lapscheme%boundedgefastrecs)
      !       ! call fvl_deallocate(omega_lapscheme%boundfastcrecs)
      ! end if
end subroutine fvl_optimize_memory2

subroutine fvl_optimize_memory3()
      ! if(staticpsimodel==1 .and. statictmodel==1 .and. statictdmodel==1) then
      !       call fvl_delete(omega_divscheme)
      ! else if(staticpsimodel==1) then
      !       ! call fvl_deallocate(omega_divscheme%celltoedgemaps)
      !       call fvl_deallocate(omega_divscheme%edgesalias)
      !       ! call fvl_deallocate(omega_divscheme%edgecoeffields)
      !       ! call fvl_deallocate(omega_divscheme%upwindcells)
      !       call fvl_deallocate(omega_divscheme%boundedgedegrees)
      !       call fvl_deallocate(omega_divscheme%celldegrees)
      !       ! call fvl_deallocate(omega_divscheme%edgefluxcoefs)
      !       ! call fvl_deallocate(omega_divscheme%boundedgestencils)
      !       call fvl_deallocate(omega_divscheme%cellstencils)
      !       ! call fvl_deallocate(omega_divscheme%edgefastcrecs)
      !       call fvl_deallocate(omega_divscheme%edgefastrecs)
      !       ! call fvl_deallocate(omega_divscheme%boundfastcrecs)
      ! end if
end subroutine fvl_optimize_memory3

subroutine fvl_optimize_memory4()
      ! call fvl_delete(psi_lapscheme)
      ! call fvl_delete(omega_integscheme)
      ! call fvl_delete(omega_lapscheme)
      ! call fvl_delete(omega_divscheme)

      ! call fvl_delete(omega_integmat)
      ! call fvl_delete(omega_staticmat)

      ! call fvl_deallocate(omega_globalmat)
      ! call fvl_delete(omega_precondmat)
      ! call fvl_delete(omega_integmat_comp)
      ! call fvl_delete(omega_staticmat_comp)
      ! call fvl_deallocate(omega_globalmatcomp)
end subroutine fvl_optimize_memory4

end program fvl_pro_thermopsiomega2d_main
! end of file
