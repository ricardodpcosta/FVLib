! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: 2D non-isothermal Newtonian incompressible fluid flow solver
! Modification: February, 2025

#include "fvl_boussinesq_forcesource2d_mod.f"
#include "fvl_rrblock_linearsystem_mod.f"

#include "macros.f"

program fvl_pro_thermoincompflow2d_main

use fvl_lib2d
use fvl_boussinesq_forcesource2d_mod
use fvl_rrblock_linearsystem_mod

#define __fvl_max_num_steps__ 5

implicit none

! ============================================================================
! DECLARE VARIABLES
! ============================================================================

! ----------------------------------------------------------------------------
! variables
! ----------------------------------------------------------------------------

integer(kind=__fvl_integer_kind__),parameter::numsteps=__fvl_max_num_steps__
integer(kind=__fvl_integer_kind__)::step,t_dofs,p_dofs,u_dofs,& !!!!!
      statickmodel,staticcpmodel,staticmumodel,staticrhomodel,staticrhodmodel,staticalphamodel,staticnumodel,&
      statictmodel,statictcmodel,staticpmodel,staticumodel,staticucmodel,staticvcmodel,&
      staticftmodel,staticfpumodel,staticgpumodel,staticgmodel,staticbetamodel
real(kind=__fvl_real_kind__),allocatable,dimension(:,:)::t_globalsource,t_staticrhs,t_globalrhs!,t_staticsource
real(kind=__fvl_real_kind__),allocatable,dimension(:,:)::td_globalsource,td_staticrhs,td_globalrhs!,td_staticsource
real(kind=__fvl_real_kind__),allocatable,dimension(:,:)::tdd_globalsource,tdd_staticrhs,tdd_globalrhs!,tdd_staticsource
real(kind=__fvl_real_kind__),allocatable,dimension(:,:)::pu_globalsource
real(kind=__fvl_real_kind__),allocatable,dimension(:)::pu_staticrhs,pu_globalrhs!,pu_staticsource
real(kind=__fvl_real_kind__),allocatable,dimension(:)::u_globalsource,u_staticrhs,u_globalrhs!,u_staticsource
type(fvl_const_scalarfield2d)::xi_field
type(fvl_ptr_model)::tx_models(4*__fvl_max_num_steps__)
type(fvl_ptr_model)::pu_models(3)
type(fvl_ptr_model)::tpu_models(7)
type(fvl_diamondpro_interpolate2d)::t_interpolate
type(fvl_centredpro_integratescheme2d)::t_integscheme
type(fvl_centredpro_laplacianscheme2d)::t_lapscheme
type(fvl_upwindpro_divergencescheme2d)::t_divscheme
type(fvl_diamondpro_gradientscheme2d)::p_gradscheme
type(fvl_centredpro_interpolate2d)::ux_interpolate
type(fvl_centredpro_integratescheme2d)::ux_integscheme
type(fvl_centredpro_vectorlaplacianscheme2d)::ux_lapscheme
type(fvl_upwindpro_vectordivergencescheme2d)::ux_divscheme
type(fvl_primalpro_divergencescheme2d)::u_divscheme
type(fvl_primalpro_interpolate2d)::u_interpolate
type(fvl_lil_rspmat)::t_precondmat,p_precondmat,ux_precondmat
type(fvl_lil_rspmat)::t_staticmat,t_integmat,t_globalmat
type(fvl_csr_rspmat)::t_staticmat_comp,t_integmat_comp,t_globalmat_comp
type(fvl_lil_brspmat)::p_globalmat
type(fvl_csr_brspmat)::p_globalmat_comp
type(fvl_lil_rspmat)::ux_staticmat1,ux_integmat,ux_globalmat1
type(fvl_lil_brspmat)::ux_staticmat2,ux_globalmat2
type(fvl_csr_rspmat)::ux_staticmat_comp1,ux_integmat_comp,ux_globalmat_comp1
type(fvl_csr_brspmat)::ux_staticmat_comp2,ux_globalmat_comp2
type(fvl_lil_brspmat)::u_globalmat
type(fvl_csr_brspmat)::u_globalmat_comp
type(fvl_stationary_linearsystem)::t_linearsystem
type(fvl_thetarule_linearsystem)::pu_linearsystem
type(fvl_linearsolver)::t_linearsolver,pu_linearsolver
type(fvl_preconditioner)::t_preconditioner,pu_preconditioner
type(fvl_fixpointiterator)::t_fixpointiterator,p_fixpointiterator,u_fixpointiterator
type(fvl_timeiterator)::timeiterator

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

! model allocators
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
      call fvl_loginfo("About: 2D non-isothermal non-Newtonian incompressible fluid flow solver.")
      call fvl_loginfo("Usage: fvl-pro-thermoincompflow1d [options]")
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

! read control parameters
statickmodel            = controldictfile%getvalue("solver","static_k_model",0)
staticcpmodel           = controldictfile%getvalue("solver","static_cp_model",0)
staticmumodel           = controldictfile%getvalue("solver","static_mu_model",0)
staticrhomodel          = controldictfile%getvalue("solver","static_rho_model",0)
staticrhodmodel         = controldictfile%getvalue("solver","static_rhod_model",0)
statictmodel            = controldictfile%getvalue("solver","static_t_model",0)
statictcmodel           = controldictfile%getvalue("solver","static_tc_model",0)
staticpmodel            = controldictfile%getvalue("solver","static_p_model",0)
staticumodel            = controldictfile%getvalue("solver","static_u_model",0)
staticucmodel           = controldictfile%getvalue("solver","static_uc_model",0)
staticvcmodel           = controldictfile%getvalue("solver","static_vc_model",0)
staticftmodel           = controldictfile%getvalue("solver","static_ft_model",0)
staticfpumodel          = controldictfile%getvalue("solver","static_fpu_model",0)
staticgpumodel          = controldictfile%getvalue("solver","static_gpu_model",0)
staticgmodel            = controldictfile%getvalue("solver","static_g_model",0)
staticbetamodel         = controldictfile%getvalue("solver","static_beta_model",0)
staticalphamodel        = statickmodel*staticrhomodel*staticcpmodel
staticnumodel           = staticmumodel*staticrhodmodel

! set time directories precision
call fvl_settimedirectoryprecision(solutiontimeprecision)

! ============================================================================
! INITIALIZE MESHES
! ============================================================================

! state message
call fvl_loginfo("Initializing meshes...")

! start meshes initialization timer
call fvl_startwtimer(20,"meshes_total")

! templated code
#include "fvl_creatediamondmesh2d_src.f"

! user-defined code
#include "meshes.f"

! stop meshes initialization timer
call fvl_stopwtimer(20)

! save memory usage
call fvl_setvmhwm(1,"meshes_total")

! ============================================================================
! INITIALIZE PHYSICAL MODELS
! ============================================================================

! state message
call fvl_loginfo("Initializing physical models...")

! start models initialization timer
call fvl_startwtimer(30,"models_total")

! user-defined code
#include "models.f"

! ----------------------------------------------------------------------------
! temperature
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing temperature model")

do step=1,numsteps

! temperature model
t_models(step)          = fvl_scalarmodel2d_init44(mesh=mesh,fieldtype=fvl_field2d_cellmeansfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="t_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=t_fieldselector,boundcondselector=t_boundcondselector,&
                              label="t",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.true.)

! temperature model
td_models(step)         = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_cellmeansfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="td_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=t_fieldselector,boundcondselector=t_boundcondselector,&
                              label="td",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.true.)

! temperature model
tdd_models(step)        = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_cellmeansfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="tdd_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=t_fieldselector,boundcondselector=t_boundcondselector,&
                              label="tdd",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.true.)

! temperature model
tddd_models(step)       = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_cellmeansfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="tddd_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=t_fieldselector,boundcondselector=t_boundcondselector,&
                              label="tddd",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.true.)

end do

! ----------------------------------------------------------------------------
! pressure
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing pressure model")

! pressure model
p_model                 = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_cellmeansfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="p_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=p_fieldselector,boundcondselector=p_boundcondselector,&
                              label="p",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.true.)

! ----------------------------------------------------------------------------
! velocity
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing velocity model")

! velocity model
u_model                 = fvl_vectormodel2d(mesh=diamondmesh,fieldtype=fvl_field2d_cellmeansfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="u_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=u_fieldselector,boundcondselector=u_boundcondselector,&
                              label="u_diamond",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.true.)

! velocity model
uc_model                 = fvl_vectormodel2d(mesh=mesh,fieldtype=fvl_field2d_edgequadratsfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="u_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=u_fieldselector,boundcondselector=u_boundcondselector,&
                              label="uc",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.true.)

! velocity model
vc_model                 = fvl_vectormodel2d(mesh=diamondmesh,fieldtype=fvl_field2d_edgequadratsfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="u_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=u_fieldselector,boundcondselector=u_boundcondselector,&
                              label="vc",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.true.)

! ----------------------------------------------------------------------------
! xi
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing xi model")

! xi field
xi_field                = fvl_const_scalarfield2d(patch=cellspatch,type=fvl_field2d_cellmeansfield,&
                              label="xi",filepath=solutionfiledirectory,fileform=solutionfileform)

! xi model
xi_model                = fvl_scalarmodel2d(mesh=mesh,field=xi_field,&
                              label="xi",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.false.)

! ----------------------------------------------------------------------------
! thermal conductivity
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing thermal conductivity model")

! thermal conductivity model
k_model                 = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_edgequadratsfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="k_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=k_fieldselector,boundcondselector=k_boundcondselector,&
                              label="k",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.false.)

! ----------------------------------------------------------------------------
! heat capacity
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing heat capacity model")

! heat capacity model
cp_model                = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_edgequadratsfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="cp_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=cp_fieldselector,boundcondselector=cp_boundcondselector,&
                              label="cp",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.false.)

! ----------------------------------------------------------------------------
! dynamic viscosity
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing dynamic viscosity model")

! dynamic viscosity model
mu_model                = fvl_scalarmodel2d(mesh=diamondmesh,fieldtype=fvl_field2d_edgequadratsfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="mu_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=mu_fieldselector,boundcondselector=mu_boundcondselector,&
                              label="mu",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.false.)

! ----------------------------------------------------------------------------
! density
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing density model")

! density model
rho_model               = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_edgequadratsfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="rho_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=rho_fieldselector,boundcondselector=rho_boundcondselector,&
                              label="rho",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.false.)

! density model
rhod_model               = fvl_scalarmodel2d(mesh=diamondmesh,fieldtype=fvl_field2d_edgequadratsfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="rho_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=rho_fieldselector,boundcondselector=rho_boundcondselector,&
                              label="rhod",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.false.)

! ----------------------------------------------------------------------------
! thermal diffusivity
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing thermal diffusivity model")

! thermal diffusivity model
alpha_model             = -k_model/(rho_model*cp_model)

! ----------------------------------------------------------------------------
! kinematic viscosity
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing kinematic viscosity model")

! kinematic viscosity model
nu_model                = -mu_model/rhod_model

! ----------------------------------------------------------------------------
! specific volume
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing specific volume model")

! specific volume model
vp_model                = 1.0d0/rhod_model

! ----------------------------------------------------------------------------
! heat source
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing heat source model")

! temperature model
tc_model                = fvl_scalarmodel2d(mesh=diamondmesh,fieldtype=fvl_field2d_cellquadratsfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="t_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=t_fieldselector,boundcondselector=t_boundcondselector,&
                              label="tc",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.true.)

! heat source model
ft_model                = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_cellquadratsfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="ft_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=ft_fieldselector,boundcondselector=ft_boundcondselector,&
                              label="ft",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.false.)

! heat source model
ftd_model                = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_cellquadratsfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="ftd_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=ft_fieldselector,boundcondselector=ft_boundcondselector,&
                              label="ftd",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.false.)

! heat source model
ftdd_model              = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_cellquadratsfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="ftdd_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=ft_fieldselector,boundcondselector=ft_boundcondselector,&
                              label="ftdd",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.false.)

! ----------------------------------------------------------------------------
! force source
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing force source model")

! force source model
fpu_model               = fvl_vectormodel2d(mesh=diamondmesh,fieldtype=fvl_field2d_cellquadratsfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="fpu_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=fpu_fieldselector,boundcondselector=fpu_boundcondselector,&
                              label="fpu",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.false.)

! ----------------------------------------------------------------------------
! mass source
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing mass source model")

! mass source model
gpu_model               = fvl_scalarmodel2d(mesh=mesh,fieldtype=fvl_field2d_cellquadratsfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="gpu_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=gpu_fieldselector,boundcondselector=gpu_boundcondselector,&
                              label="gpu",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.false.)

! ----------------------------------------------------------------------------
! gravitational acceleration
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing gravitational acceleration model")

! gravitational acceleration model
g_model               = fvl_vectormodel2d(mesh=diamondmesh,fieldtype=fvl_field2d_cellquadratsfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="g_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=g_fieldselector,boundcondselector=g_boundcondselector,&
                              label="g",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.false.)

! ----------------------------------------------------------------------------
! thermal expansion
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Initializing thermal expansion model")

! thermal expansion model
beta_model               = fvl_scalarmodel2d(mesh=diamondmesh,fieldtype=fvl_field2d_cellquadratsfield,&
                              modeldictfile=modelsdictfile,modeldictlabel="beta_model",&
                              patchesdictfile=meshesdictfile,patchesdictlabel="patches",&
                              fieldselector=beta_fieldselector,boundcondselector=beta_boundcondselector,&
                              label="beta",filepath=solutionfiledirectory,fileform=solutionfileform,calculatedfields=.false.)

! ----------------------------------------------------------------------------
! set models
! ----------------------------------------------------------------------------

! set models
do step=1,__fvl_max_num_steps__
      call tx_models((step-1)*4+1)%setptr(t_models(step))
      call tx_models((step-1)*4+2)%setptr(td_models(step))
      call tx_models((step-1)*4+3)%setptr(tdd_models(step))
      call tx_models((step-1)*4+4)%setptr(tddd_models(step))
end do

call pu_models(1)%setptr(u_model)
call pu_models(2)%setptr(p_model)
call pu_models(3)%setptr(xi_model)

call tpu_models(1)%setptr(t_model)
call tpu_models(2)%setptr(td_model)
call tpu_models(3)%setptr(tdd_model)
call tpu_models(4)%setptr(tddd_model)
call tpu_models(5)%setptr(u_model)
call tpu_models(6)%setptr(p_model)
call tpu_models(7)%setptr(xi_model)

call fvl_allocate(solutionmodels,7)
call solutionmodels(1)%setptr(t_model)
call solutionmodels(2)%setptr(td_model)
call solutionmodels(3)%setptr(tdd_model)
call solutionmodels(4)%setptr(tddd_model)
call solutionmodels(5)%setptr(u_model)
call solutionmodels(6)%setptr(p_model)
call solutionmodels(7)%setptr(xi_model)

! ----------------------------------------------------------------------------
! allocate vectors
! ----------------------------------------------------------------------------

! initialize variables
t_dofs                  = mesh%getnumcells()
p_dofs                  = mesh%getnumcells()
u_dofs                  = 2*diamondmesh%getnumcells()

! allocate memory
! call fvl_allocate(t_staticsource,t_dofs)
call fvl_allocate(t_globalsource,t_dofs,numsteps)
call fvl_allocate(t_staticrhs,t_dofs,numsteps)
call fvl_allocate(t_globalrhs,t_dofs,numsteps)
call fvl_allocate(td_globalsource,t_dofs,numsteps)
call fvl_allocate(td_staticrhs,t_dofs,numsteps)
call fvl_allocate(td_globalrhs,t_dofs,numsteps)
call fvl_allocate(tdd_globalsource,t_dofs,numsteps)
call fvl_allocate(tdd_staticrhs,t_dofs,numsteps)
call fvl_allocate(tdd_globalrhs,t_dofs,numsteps)
! call fvl_allocate(u_staticsource,p_dofs)
call fvl_allocate(u_globalsource,p_dofs)
call fvl_allocate(u_staticrhs,p_dofs)
call fvl_allocate(u_globalrhs,p_dofs)
! call fvl_allocate(pu_staticsource,2,u_dofs/2)
call fvl_allocate(pu_globalsource,2,u_dofs/2)
call fvl_allocate(pu_staticrhs,u_dofs)
call fvl_allocate(pu_globalrhs,u_dofs)

! clean memory
! call fvl_omp_clean(t_staticsource)
call fvl_omp_clean(t_globalsource)
call fvl_omp_clean(t_staticrhs)
call fvl_omp_clean(t_globalrhs)
call fvl_omp_clean(td_globalsource)
call fvl_omp_clean(td_staticrhs)
call fvl_omp_clean(td_globalrhs)
call fvl_omp_clean(tdd_globalsource)
call fvl_omp_clean(tdd_staticrhs)
call fvl_omp_clean(tdd_globalrhs)
! call fvl_omp_clean(pu_staticsource)
call fvl_omp_clean(pu_globalsource)
call fvl_omp_clean(pu_staticrhs)
call fvl_omp_clean(pu_globalrhs)
! call fvl_omp_clean(u_staticsource)
call fvl_omp_clean(u_globalsource)
call fvl_omp_clean(u_staticrhs)
call fvl_omp_clean(u_globalrhs)
call fvl_clean(t_precondmat)
call fvl_clean(p_precondmat)
call fvl_clean(ux_precondmat)
call fvl_clean(t_staticmat)
call fvl_clean(t_integmat)
call fvl_clean(t_globalmat)
call fvl_clean(t_staticmat_comp)
call fvl_clean(t_integmat_comp)
call fvl_clean(t_globalmat_comp)
call fvl_clean(p_globalmat)
call fvl_clean(p_globalmat_comp)
call fvl_clean(ux_staticmat1)
call fvl_clean(ux_integmat)
call fvl_clean(ux_globalmat1)
call fvl_clean(ux_staticmat2)
call fvl_clean(ux_globalmat2)
call fvl_clean(ux_staticmat_comp1)
call fvl_clean(ux_integmat_comp)
call fvl_clean(ux_globalmat_comp1)
call fvl_clean(ux_staticmat_comp2)
call fvl_clean(ux_globalmat_comp2)
call fvl_clean(u_globalmat)
call fvl_clean(u_globalmat_comp)

! ! initialize pressure pseudo patches
! call t_initialize_pseudo_patches()

! ! initialize pressure pseudo patches
! call p_initialize_pseudo_patches()

! ! initialize velocity pseudo patches
! call u_initialize_pseudo_patches()

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

! ============================================================================
! ENERGY BALANCE EQUATION
! ============================================================================

! ----------------------------------------------------------------------------
! temperature integrate scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Computing temperature integrate scheme")

! start timer
call fvl_startwtimer(43,"t_integrate_scheme")

! initialize scheme
t_integscheme           = fvl_centredpro_integratescheme2d(mesh=mesh,&
                              dictfile=schemesdictfile,dictlabel="t_integrate_scheme",&
                              staticcoeffields=.true.)

! evaluate coefficients matrix
call t_integscheme%evalmatrix(t_integmat)

! compress coefficients matrix
call t_integmat_comp%compress(t_integmat)

! deallocate memory
call fvl_delete(t_integscheme)
call fvl_delete(t_integmat)

! stop timer
call fvl_stopwtimer(43)

! ----------------------------------------------------------------------------
! temperature laplacian scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Computing temperature laplacian scheme")

! start timer
call fvl_startwtimer(44,"t_laplacian_scheme")

! initialize scheme
t_lapscheme             = fvl_centredpro_laplacianscheme2d(coefsmodel=alpha_model,condsmodel=t_model,&
                              mesh=mesh,dictfile=schemesdictfile,dictlabel="t_laplacian_scheme",&
                              staticcoeffields=(staticalphamodel==1),staticboundconds=.true.)

! static terms
if(staticalphamodel==1) then
      ! evaluate coefficients matrix
      call t_lapscheme%evalmatrix(t_staticmat)
      ! evaluate right-hand side
      call t_lapscheme%updatemodels(alpha_model,t_model)
      call t_lapscheme%evalrhs(t_staticrhs(:,1))
      call t_lapscheme%updatemodels(alpha_model,td_model)
      call t_lapscheme%evalrhs(td_staticrhs(:,1))
      call t_lapscheme%updatemodels(alpha_model,tdd_model)
      call t_lapscheme%evalrhs(tdd_staticrhs(:,1))
      ! deallocate memory
      call fvl_delete(t_lapscheme)
end if

! stop timer
call fvl_stopwtimer(44)

! ----------------------------------------------------------------------------
! temperature divergence scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Computing temperature divergence scheme")

! start timer
call fvl_startwtimer(45,"t_divergence_scheme")

! initialize scheme
t_divscheme             = fvl_upwindpro_divergencescheme2d(coefsmodel=uc_model,condsmodel=t_model,&
                              mesh=mesh,dictfile=schemesdictfile,dictlabel="t_divergence_scheme",&
                              staticcoeffields=(staticucmodel==1),staticboundconds=.true.)

! static terms
if(staticucmodel==1) then
      ! evaluate coefficients matrix
      call t_divscheme%evaladdmatrix(t_staticmat)
      ! evaluate right-hand side
      call t_divscheme%updatemodels(uc_model,t_model)
      call t_divscheme%evaladdrhs(t_staticrhs(:,1))
      call t_divscheme%updatemodels(uc_model,td_model)
      call t_divscheme%evaladdrhs(td_staticrhs(:,1))
      call t_divscheme%updatemodels(uc_model,tdd_model)
      call t_divscheme%evaladdrhs(tdd_staticrhs(:,1))
      ! deallocate memory
      call fvl_delete(t_divscheme)
end if

! stop timer
call fvl_stopwtimer(45)

! ----------------------------------------------------------------------------
! velocity interpolation scheme
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Computing velocity interpolate scheme")

! start timer
call fvl_startwtimer(42,"u_interpolate_scheme")

! initialize scheme
u_interpolate           = fvl_primalpro_interpolate2d(patch=edgespatch,mesh=mesh,diamondmesh=diamondmesh,& !!!!
                              type=fvl_primalpro_interpolate2d_edgequadratvalues,&
                              dictfile=schemesdictfile,dictlabel="u_interpolate_scheme")

! stop timer
call fvl_stopwtimer(42)

! ----------------------------------------------------------------------------
! static terms
! ----------------------------------------------------------------------------

! static terms
if(staticalphamodel==1 .and. staticucmodel==1) then
      ! compress coefficients matrices
      call t_globalmat_comp%compress(t_staticmat)
      ! deallocate memory
      call fvl_delete(t_staticmat)
else if(staticalphamodel==1 .or. staticucmodel==1) then
      ! compress coefficients matrices
      call t_staticmat_comp%compress(t_staticmat)
      ! deallocate memory
      call fvl_delete(t_staticmat)
end if

! static terms
if(staticalphamodel==1 .and. staticucmodel==1) then
      ! set global right-hand side
      call fvl_omp_assign(t_staticrhs(:,1),t_globalrhs(:,1))
      call fvl_omp_assign(td_staticrhs(:,1),td_globalrhs(:,1))
      call fvl_omp_assign(tdd_staticrhs(:,1),tdd_globalrhs(:,1))
      ! deallocate memory
      call fvl_deallocate(t_staticrhs)
      call fvl_deallocate(td_staticrhs)
      call fvl_deallocate(tdd_staticrhs)
end if

! static terms
if(staticftmodel==1) then
      ! evaluate heat source term
      call ft_model%evalcellquadrattomeanvalues(t_globalsource(:,1))
      call ftd_model%evalcellquadrattomeanvalues(td_globalsource(:,1))
      call ftdd_model%evalcellquadrattomeanvalues(tdd_globalsource(:,1))
end if

! ! static terms
! if(staticftmodel==1) then
!       ! set global source term
!       call fvl_omp_assign(t_staticsource,t_globalsource)
!       ! deallocate memory
!       call fvl_deallocate(t_staticsource)
! end if

! ! ============================================================================
! ! MOMENTUM BALANCE EQUATION
! ! ============================================================================
!
! ! ----------------------------------------------------------------------------
! ! pressure gradient scheme
! ! ----------------------------------------------------------------------------
!
! ! state message
! call fvl_loginfo(">>> Computing pressure gradient scheme")
!
! ! start timer
! call fvl_startwtimer(41,"p_gradient_scheme")
!
! ! initialize scheme
! p_gradscheme            = fvl_diamondpro_gradientscheme2d(coefsmodel=vp_model,condsmodel=p_model,&
!                               mesh=diamondmesh,primalmesh=mesh,dictfile=schemesdictfile,dictlabel="p_gradient_scheme",&
!                               staticcoeffields=(staticrhodmodel==1),staticboundconds=.true.)
!
! ! static terms
! if(staticrhodmodel==1) then
!       ! evaluate coefficients matrix
!       call p_gradscheme%evalmatrix(p_globalmat)
!       ! compress coefficients matrix
!       call p_globalmat_comp%compress(p_globalmat)
!       ! deallocate memory
!       call fvl_delete(p_globalmat)
!       ! evaluate right-hand side
!       call p_gradscheme%evalrhs(pu_staticrhs)
!       ! deallocate memory
!       call fvl_delete(p_gradscheme)
! end if
!
! ! stop timer
! call fvl_stopwtimer(41)
!
! ! ----------------------------------------------------------------------------
! ! velocity integrate scheme
! ! ----------------------------------------------------------------------------
!
! ! state message
! call fvl_loginfo(">>> Computing velocity integrate scheme")
!
! ! start timer
! call fvl_startwtimer(43,"u_integrate_scheme")
!
! ! initialize scheme
! ux_integscheme          = fvl_centredpro_integratescheme2d(mesh=diamondmesh,&
!                               dictfile=schemesdictfile,dictlabel="u_integrate_scheme",&
!                               staticcoeffields=.true.)
!
! ! evaluate coefficients matrix
! call ux_integscheme%evalmatrix(ux_integmat)
!
! ! compress coefficients matrix
! call ux_integmat_comp%compress(ux_integmat)
!
! ! deallocate memory
! call fvl_delete(ux_integscheme)
! call fvl_delete(ux_integmat)
!
! ! stop timer
! call fvl_stopwtimer(43)
!
! ! ----------------------------------------------------------------------------
! ! velocity laplacian scheme
! ! ----------------------------------------------------------------------------
!
! ! state message
! call fvl_loginfo(">>> Computing velocity vector laplacian scheme")
!
! ! start timer
! call fvl_startwtimer(44,"u_vector_laplacian_scheme")
!
! ! initialize scheme
! ux_lapscheme            = fvl_centredpro_vectorlaplacianscheme2d(coefsmodel=nu_model,condsmodel=u_model,&
!                               mesh=diamondmesh,dictfile=schemesdictfile,dictlabel="u_vector_laplacian_scheme",&
!                               staticcoeffields=(staticnumodel==1),staticboundconds=.true.)
!
! ! static terms
! if(staticnumodel==1) then
!       ! evaluate coefficients matrix
!       call ux_lapscheme%evalmatrix(ux_staticmat1,ux_staticmat2)
!       ! evaluate right-hand side
!       call ux_lapscheme%evaladdrhs(pu_staticrhs)
!       ! deallocate memory
!       call fvl_delete(ux_lapscheme)
! end if
!
! ! stop timer
! call fvl_stopwtimer(44)
!
! ! ----------------------------------------------------------------------------
! ! velocity divergence scheme
! ! ----------------------------------------------------------------------------
!
! ! state message
! call fvl_loginfo(">>> Computing velocity vector divergence scheme")
!
! ! start timer
! call fvl_startwtimer(45,"u_vector_divergence_scheme")
!
! ! initialize scheme
! ux_divscheme            = fvl_upwindpro_vectordivergencescheme2d(coefsmodel=vc_model,condsmodel=u_model,&
!                               mesh=diamondmesh,dictfile=schemesdictfile,dictlabel="u_vector_divergence_scheme",&
!                               staticcoeffields=(staticvcmodel==1),staticboundconds=.true.)
!
! ! static terms
! if(staticvcmodel==1) then
!       ! evaluate coefficients matrix
!       call ux_divscheme%evaladdmatrix(ux_staticmat1,ux_staticmat2)
!       ! evaluate right-hand side
!       call ux_divscheme%evaladdrhs(pu_staticrhs)
!       ! deallocate memory
!       call fvl_delete(ux_divscheme)
! end if
!
! ! stop timer
! call fvl_stopwtimer(45)
!
! ! ----------------------------------------------------------------------------
! ! temperature interpolation scheme
! ! ----------------------------------------------------------------------------
!
! ! state message
! call fvl_loginfo(">>> Computing temperature interpolate scheme")
!
! ! start timer
! call fvl_startwtimer(42,"t_interpolate_scheme")
!
! ! initialize scheme
! t_interpolate           = fvl_diamondpro_interpolate2d(patch=diamond_cellspatch,mesh=diamondmesh,primalmesh=mesh,&
!                               type=fvl_diamondpro_interpolate2d_cellquadratvalues,&
!                               dictfile=schemesdictfile,dictlabel="t_interpolate_scheme")
!
! ! stop timer
! call fvl_stopwtimer(42)
!
! ! ----------------------------------------------------------------------------
! ! velocity interpolation scheme
! ! ----------------------------------------------------------------------------
!
! ! state message
! call fvl_loginfo(">>> Computing velocity interpolate scheme")
!
! ! start timer
! call fvl_startwtimer(42,"u_interpolate_scheme")
!
! ! initialize scheme
! ux_interpolate          = fvl_centredpro_interpolate2d(patch=diamond_edgespatch,mesh=diamondmesh,&
!                               type=fvl_centredpro_interpolate2d_edgequadratvalues,&
!                               dictfile=schemesdictfile,dictlabel="u_interpolate_scheme")
!
! ! stop timer
! call fvl_stopwtimer(42)
!
! ! ----------------------------------------------------------------------------
! ! static terms
! ! ----------------------------------------------------------------------------
!
! ! static terms
! if(staticnumodel==1 .and. staticvcmodel==1) then
!       ! compress coefficients matrices
!       call ux_globalmat_comp1%compress(ux_staticmat1)
!       call ux_globalmat_comp2%compress(ux_staticmat2)
!       ! deallocate memory
!       call fvl_delete(ux_staticmat1)
!       call fvl_delete(ux_staticmat2)
! else if(staticnumodel==1 .or. staticvcmodel==1) then
!       ! compress coefficients matrices
!       call ux_staticmat_comp1%compress(ux_staticmat1)
!       call ux_staticmat_comp2%compress(ux_staticmat2)
!       ! deallocate memory
!       call fvl_delete(ux_staticmat1)
!       call fvl_delete(ux_staticmat2)
! end if
!
! ! static terms
! if(staticrhodmodel==1 .and. staticnumodel==1 .and. staticvcmodel==1) then
!       ! set global right-hand side
!       call fvl_omp_assign(pu_staticrhs,pu_globalrhs)
!       ! deallocate memory
!       call fvl_deallocate(pu_staticrhs)
! end if
!
! ! static terms
! if(staticfpumodel==1) then
!       ! evaluate momentum source term
!       call fpu_model%evalcellquadrattomeanvalues(pu_globalsource)
! end if
!
! ! ! static terms
! ! if(staticfpumodel==1) then
! !       ! set global source term
! !       call fvl_omp_assign(pu_staticsource,pu_globalsource)
! !       ! deallocate memory
! !       call fvl_deallocate(pu_staticsource)
! ! end if
!
! ! ============================================================================
! ! MASS CONSERVATION EQUATION
! ! ============================================================================
!
! ! ----------------------------------------------------------------------------
! ! velocity divergence scheme
! ! ----------------------------------------------------------------------------
!
! ! state message
! call fvl_loginfo(">>> Computing velocity divergence scheme")
!
! ! start timer
! call fvl_startwtimer(46,"u_divergence_scheme")
!
! ! initialize scheme
! u_divscheme             = fvl_primalpro_divergencescheme2d(condsmodel=u_model,&
!                               mesh=mesh,diamondmesh=diamondmesh,dictfile=schemesdictfile,dictlabel="u_divergence_scheme",&
!                               staticcoeffields=.true.,staticboundconds=.true.)
!
! ! evaluate coefficients matrix
! call u_divscheme%evalmatrix(u_globalmat)
!
! ! compress coefficients matrix
! call u_globalmat_comp%compress(u_globalmat)
! call fvl_delete(u_globalmat)
!
! ! evaluate right-hand side
! call u_divscheme%evalrhs(u_globalrhs)
!
! ! deallocate memory
! call fvl_delete(u_divscheme)
!
! ! stop timer
! call fvl_stopwtimer(46)
!
! ! stop timer
! call fvl_stopwtimer(40)
!
! ! save memory usage
! call fvl_setvmhwm(3,"schemes_total")
!
! ! ----------------------------------------------------------------------------
! ! static terms
! ! ----------------------------------------------------------------------------
!
! ! static terms
! if(staticgpumodel==1) then
!       ! evaluate mass source term
!       call gpu_model%evalcellquadrattomeanvalues(u_globalsource)
! end if
!
! ! ! static terms
! ! if(staticgpumodel==1) then
! !       ! set global source term
! !       call fvl_omp_assign(u_staticsource,u_globalsource)
! !       ! deallocate memory
! !       call fvl_deallocate(u_staticsource)
! ! end if

! ============================================================================
! COMPUTE SOLUTION
! ============================================================================

! state message
call fvl_loginfo("Computing solution...")

! start solution timer
call fvl_startwtimer(50,"solution_total")

! ----------------------------------------------------------------------------
! linear systems
! ----------------------------------------------------------------------------

! temperature linear system
t_linearsystem          = fvl_stationary_linearsystem(runtime=mesh,model=t_model,&
                              tfun=t_tfun,afun=t_afun,bfun=t_bfun1,pfun=t_pfun,ufun=t_ufun1,&
                              linearsolver=t_linearsolver,dictfile=solutiondictfile,dictlabel="t_linear_system")

! ! temperature linear system
! t_linearsystem          = fvl_thetarule_linearsystem(runtime=mesh,model=t_model,&
!                               tfun=t_tfun,afun=t_afun,bfun=t_bfun1,pfun=t_pfun,ufun=t_ufun1,&
!                               linearsolver=t_linearsolver,dictfile=solutiondictfile,dictlabel="t_linear_system")

! ! temperature linear system
! t_linearsystem          = fvl_rrblock_linearsystem(runtime=mesh,&
!                               models1=t_models,models2=td_models,models3=tdd_models,models4=tddd_models,&
!                               tfun1=t_tfun,afun1=t_afun,bfun1=t_bfun1,pfun1=t_pfun,ufun1=t_ufun1,&
!                               tfun2=t_tfun,afun2=t_afun,bfun2=t_bfun2,pfun2=t_pfun,ufun2=t_ufun2,&
!                               tfun3=t_tfun,afun3=t_afun,bfun3=t_bfun3,pfun3=t_pfun,ufun3=t_ufun3,&
!                               linearsolver=t_linearsolver,dictfile=solutiondictfile,dictlabel="t_linear_system")

! ! pressure-velocity linear system
! pu_linearsystem         = fvl_thetarule_linearsystem(runtimes=runtimes,models=pu_models,&
!                               tfun=pu_tfun,afun=pu_afun,bfun=pu_bfun,pfun=pu_pfun,ufun=pu_ufun,&
!                               linearsolver=pu_linearsolver,dictfile=solutiondictfile,dictlabel="pu_linear_system")

! ----------------------------------------------------------------------------
! preconditioners
! ----------------------------------------------------------------------------

! temperature preconditioner
t_preconditioner        = fvl_preconditioner(dictfile=solutiondictfile,dictlabel="t_preconditioner")

! ! pressure-velocity preconditioner
! pu_preconditioner       = fvl_preconditioner(dictfile=solutiondictfile,dictlabel="pu_preconditioner")

! ----------------------------------------------------------------------------
! linear solvers
! ----------------------------------------------------------------------------

! temperature linear solver
t_linearsolver          = fvl_linearsolver(dictfile=solutiondictfile,dictlabel="t_linear_solver")

! ! pressure-velocity linear solver
! pu_linearsolver         = fvl_linearsolver(dictfile=solutiondictfile,dictlabel="pu_linear_solver")

! ----------------------------------------------------------------------------
! fixpoint iterators
! ----------------------------------------------------------------------------

! temperature fixpoint iterator
t_fixpointiterator      = fvl_fixpointiterator(runtime=mesh,model=t_model,& !!!!!!!!!!!!!!!!!!!!!
                              dictfile=solutiondictfile,dictlabel="t_fixpoint_iterator")

! ! pressure fixpoint iterator
! p_fixpointiterator      = fvl_fixpointiterator(runtimes=runtimes,models=p_models,&
!                               dictfile=solutiondictfile,dictlabel="p_fixpoint_iterator")
!
! ! velocity fixpoint iterator
! u_fixpointiterator      = fvl_fixpointiterator(runtimes=runtimes,models=u_models,&
!                               dictfile=solutiondictfile,dictlabel="u_fixpoint_iterator")

! ----------------------------------------------------------------------------
! time iterator
! ----------------------------------------------------------------------------

! time iterator
timeiterator            = fvl_timeiterator(runtimes=runtimes,models=solutionmodels,writetime=writetime,&
                              dictfile=solutiondictfile,dictlabel="time_iterator")

! ----------------------------------------------------------------------------
! solution algorithm
! ----------------------------------------------------------------------------

! restart time interator
call timeiterator%restart()

! compute solution
do while(timeiterator%iterate())

      ! state message
      call fvl_loginfo("")

      ! update time
      call timeiterator%updatetime()

      ! state message
      call fvl_loginfo("")

      ! restart linearsolvers
      call t_linearsolver%restart()
      ! call pu_linearsolver%restart()

      ! restart fixpoint iterators
      call t_fixpointiterator%restart()
      ! call p_fixpointiterator%restart()
      ! call u_fixpointiterator%restart()

      ! restart preconditioners
      call t_preconditioner%restart() !!!
      ! call pu_preconditioner%restart() !!!

      ! update time
      !call t_updatetime()
      !call pu_updatetime()

      ! update linear systems
      call t_linearsystem%update()
      ! call pu_linearsystem%update()

      do while(t_fixpointiterator%iterate())! .or. p_fixpointiterator%iterate() .or. u_fixpointiterator%iterate())

            ! state message
            call fvl_loginfo(">>> Solving for temperature")

            ! solve linear system
            call t_linearsystem%solve()
            ! call fvl_rblock_linearsystem_solvesegregatedd(t_linearsystem)

            ! state message
            call fvl_loginfo("")

            ! check convergence and relax solution
            call t_fixpointiterator%updatefields()

            ! ! state message
            ! call fvl_loginfo("")
            !
            ! ! state message
            ! call fvl_loginfo(">>> Solving for pressure and velocity")

            ! ! solve linear system
            ! call pu_linearsystem%solve()
            !
            ! ! check convergence and relax solution
            ! call p_fixpointiterator%updatefields()
            ! call u_fixpointiterator%updatefields()

            ! state message
            call fvl_loginfo("")

      end do

      ! check convergence
      call timeiterator%updatefields()

end do

! stop solution timer
call fvl_stopwtimer(50)

! save memory usage
call fvl_setvmhwm(4,"solution_total")

! ----------------------------------------------------------------------------
! write histograms
! ----------------------------------------------------------------------------

! time iterator
call timeiterator%writehistogram(fvl_trim(solutionfiledirectory)//"/timeiterator.fvd","ascii")

! ============================================================================
! POST-PROCESSING DATA
! ============================================================================

! state message
call fvl_loginfo("Post-processing data...")

! start post-processing timer
call fvl_startwtimer(60,"postprocessing_total")

! deallocate memory
call fvl_delete(u_interpolate)

! ----------------------------------------------------------------------------
! user-defined post-processing
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Running user-defined post-processing")

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
! model functions
! ----------------------------------------------------------------------------

#include "fvl_createmodels_src.f"

! ----------------------------------------------------------------------------
! user-defined functions
! ----------------------------------------------------------------------------

#include "functions.f"

! ----------------------------------------------------------------------------
! built-in functions
! ----------------------------------------------------------------------------

! ============================================================================
! ENERGY BALANCE EQUATION
! ============================================================================

subroutine t_tfun(step,x,res)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:size(x,1)),intent(out)::res
      ! start timer
      call fvl_startwtimer(51,"solution_tfun")
      ! energy equation
      if(t_integmat_comp%getmaxnumrowvalues()/=0) then
            call t_integmat_comp%vecmul(x(1:t_dofs),res(1:t_dofs))
      else
            call fvl_omp_assign(x(1:t_dofs),res(1:t_dofs))
      end if
      ! stop timer
      call fvl_stopwtimer(51)
end subroutine t_tfun

subroutine t_afun(step,x,res)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:size(x,1)),intent(out)::res
      ! start timer
      call fvl_startwtimer(52,"solution_afun")
      ! energy equation
      call t_globalmat_comp%vecmul(x(1:t_dofs),res(1:t_dofs))
      ! stop timer
      call fvl_stopwtimer(52)
end subroutine t_afun

subroutine t_bfun1(step,res)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(inout)::res
      ! start timer
      call fvl_startwtimer(53,"solution_afun")
      ! energy equation
      call fvl_omp_add(t_globalrhs(:,step),t_globalsource(:,step),res(1:t_dofs))
      ! stop timer
      call fvl_stopwtimer(53)
end subroutine t_bfun1

subroutine t_bfun2(step,res)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(inout)::res
      ! start timer
      call fvl_startwtimer(53,"solution_afun")
      ! energy equation
      call fvl_omp_add(td_globalrhs(:,step),td_globalsource(:,step),res(1:t_dofs))
      ! stop timer
      call fvl_stopwtimer(53)
end subroutine t_bfun2

subroutine t_bfun3(step,res)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(inout)::res
      ! start timer
      call fvl_startwtimer(53,"solution_afun")
      ! energy equation
      call fvl_omp_add(tdd_globalrhs(:,step),tdd_globalsource(:,step),res(1:t_dofs))
      ! stop timer
      call fvl_stopwtimer(53)
end subroutine t_bfun3

subroutine t_pfun(step,x,res)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:size(x,1)),intent(out)::res
      ! start timer
      call fvl_startwtimer(54,"solution_pfun")
      ! momentum equation
      if(t_precondmat%getmaxnumrowvalues()/=0) then
            call t_precondmat%vecmul(x(1:t_dofs),res(1:t_dofs))
      else
            call fvl_omp_assign(x(1:t_dofs),res(1:t_dofs))
      end if
      ! stop timer
      call fvl_stopwtimer(54)
end subroutine t_pfun

subroutine t_updatetime()
      ! start timer
      call fvl_startwtimer(55,"solution_ufun")

      ! stop timer
      call fvl_stopwtimer(55)
end subroutine t_updatetime

subroutine t_ufun1(step)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      ! start timer
      call fvl_startwtimer(55,"solution_ufun")

      ! update models
      call t_lapscheme%updatemodels(alpha_model,t_models(step))
      call t_divscheme%updatemodels(uc_model,t_models(step))




      ! update velocity field
      if(staticucmodel/=1) then
            ! interpolate velocity field
            call u_interpolate%evalvalues(u_model,uc_model)
            ! update temperature divergence scheme
            call t_divscheme%update()
      end if

      ! evaluate temperature laplacian coefficients matrix
      if(staticalphamodel/=1) then
            call t_lapscheme%evalmatrix(t_globalmat)
      end if
      ! evaluate temperature divergence coefficients matrix
      if(staticucmodel/=1) then
            call t_divscheme%evaladdmatrix(t_globalmat)
      end if
      ! add and compress coefficients matrices
      if(staticalphamodel==1 .and. staticucmodel==1) then
            ! do nothing
      else if(staticalphamodel/=1 .and. staticucmodel/=1) then
            call t_globalmat_comp%compress(t_globalmat)
            call fvl_delete(t_globalmat)
      else
            call t_globalmat%addmatrix(t_staticmat_comp)
            call t_globalmat_comp%compress(t_globalmat)
            call fvl_delete(t_globalmat)
      end if




      ! clean memory
      if(staticalphamodel/=1 .or. staticucmodel/=1) then
            call fvl_omp_zeros(t_globalrhs(:,step))
      end if
      ! evaluate temperature laplacian right-hand side
      if(staticalphamodel/=1) then
            call t_lapscheme%evalrhs(t_globalrhs(:,step))
      end if
      ! evaluate temperature divergence right-hand side
      if(staticucmodel/=1) then
            call t_divscheme%evaladdrhs(t_globalrhs(:,step))
      end if
      ! add right-hand sides
      if(staticalphamodel==1 .and. staticucmodel==1) then
            ! do nothing
      else if(staticalphamodel/=1 .and. staticucmodel/=1) then
            ! do nothing
      else
            call fvl_omp_add(t_globalrhs(:,step),t_staticrhs(:,step))
      end if

      ! clean memory
      if(staticftmodel/=1) then
            call fvl_omp_zeros(t_globalsource(:,step))
      end if
      ! evaluate heat source term
      if(staticftmodel/=1) then
            call ft_model%evalcellquadrattomeanvalues(t_globalsource(:,step))
      end if
      ! ! add source terms
      ! if(.not. (staticftmodel/=1)) then
      !       call fvl_omp_add(t_globalsource(:,step),t_staticsource(:,step))
      ! end if





      ! compute preconditioner
      ! if(staticalphamodel/=1 .or. staticucmodel/=1) then
            if(t_preconditioner%iterate()) then
                  call t_preconditioner%evalprecond(t_globalmat_comp,t_precondmat)
            end if
      ! end if




      ! stop timer
      call fvl_stopwtimer(55)
end subroutine t_ufun1

subroutine t_ufun2(step)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      ! start timer
      call fvl_startwtimer(55,"solution_ufun")

      ! update models
      call t_lapscheme%updatemodels(alpha_model,td_models(step))
      call t_divscheme%updatemodels(uc_model,td_models(step))

      ! ! update velocity field
      ! if(staticucmodel/=1) then
      !       ! interpolate velocity field
      !       call u_interpolate%evalvalues(u_model,uc_model)
      !       ! update temperature divergence scheme
      !       call t_divscheme%update()
      ! end if
      !
      ! ! evaluate temperature laplacian coefficients matrix
      ! if(staticalphamodel/=1) then
      !       call t_lapscheme%evalmatrix(t_globalmat)
      ! end if
      ! ! evaluate temperature divergence coefficients matrix
      ! if(staticucmodel/=1) then
      !       call t_divscheme%evaladdmatrix(t_globalmat)
      ! end if
      ! ! add and compress coefficients matrices
      ! if(staticalphamodel==1 .and. staticucmodel==1) then
      !       ! do nothing
      ! else if(staticalphamodel/=1 .and. staticucmodel/=1) then
      !       call t_globalmat_comp%compress(t_globalmat)
      !       call fvl_delete(t_globalmat)
      ! else
      !       call t_globalmat%addmatrix(t_staticmat_comp)
      !       call t_globalmat_comp%compress(t_globalmat)
      !       call fvl_delete(t_globalmat)
      ! end if

      ! clean memory
      if(staticalphamodel/=1 .or. staticucmodel/=1) then
            call fvl_omp_zeros(td_globalrhs(:,step))
      end if
      ! evaluate temperature laplacian right-hand side
      if(staticalphamodel/=1) then
            call t_lapscheme%evalrhs(td_globalrhs(:,step))
      end if
      ! evaluate temperature divergence right-hand side
      if(staticucmodel/=1) then
            call t_divscheme%evaladdrhs(td_globalrhs(:,step))
      end if
      ! add right-hand sides
      if(staticalphamodel==1 .and. staticucmodel==1) then
            ! do nothing
      else if(staticalphamodel/=1 .and. staticucmodel/=1) then
            ! do nothing
      else
            call fvl_omp_add(td_globalrhs(:,step),td_staticrhs(:,step))
      end if

      ! clean memory
      if(staticftmodel/=1) then
            call fvl_omp_zeros(td_globalsource(:,step))
      end if
      ! evaluate heat source term
      if(staticftmodel/=1) then
            call ftd_model%evalcellquadrattomeanvalues(td_globalsource(:,step))
      end if
      ! ! add source terms
      ! if(.not. (staticftmodel/=1)) then
      !       call fvl_omp_add(td_globalsource(:,step),td_staticsource(:,step))
      ! end if

      ! ! compute preconditioner
      ! ! if(staticalphamodel/=1 .or. staticucmodel/=1) then
      !       if(t_preconditioner%iterate()) then
      !             call t_preconditioner%evalprecond(t_globalmat_comp,t_precondmat)
      !       end if
      ! ! end if

      ! stop timer
      call fvl_stopwtimer(55)
end subroutine t_ufun2

subroutine t_ufun3(step)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      ! start timer
      call fvl_startwtimer(55,"solution_ufun")

      ! update models
      call t_lapscheme%updatemodels(alpha_model,tdd_models(step))
      call t_divscheme%updatemodels(uc_model,tdd_models(step))

      ! ! update velocity field
      ! if(staticucmodel/=1) then
      !       ! interpolate velocity field
      !       call u_interpolate%evalvalues(u_model,uc_model)
      !       ! update temperature divergence scheme
      !       call t_divscheme%update()
      ! end if
      !
      ! ! evaluate temperature laplacian coefficients matrix
      ! if(staticalphamodel/=1) then
      !       call t_lapscheme%evalmatrix(t_globalmat)
      ! end if
      ! ! evaluate temperature divergence coefficients matrix
      ! if(staticucmodel/=1) then
      !       call t_divscheme%evaladdmatrix(t_globalmat)
      ! end if
      ! ! add and compress coefficients matrices
      ! if(staticalphamodel==1 .and. staticucmodel==1) then
      !       ! do nothing
      ! else if(staticalphamodel/=1 .and. staticucmodel/=1) then
      !       call t_globalmat_comp%compress(t_globalmat)
      !       call fvl_delete(t_globalmat)
      ! else
      !       call t_globalmat%addmatrix(t_staticmat_comp)
      !       call t_globalmat_comp%compress(t_globalmat)
      !       call fvl_delete(t_globalmat)
      ! end if

      ! clean memory
      if(staticalphamodel/=1 .or. staticucmodel/=1) then
            call fvl_omp_zeros(tdd_globalrhs(:,step))
      end if
      ! evaluate temperature laplacian right-hand side
      if(staticalphamodel/=1) then
            call t_lapscheme%evalrhs(tdd_globalrhs(:,step))
      end if
      ! evaluate temperature divergence right-hand side
      if(staticucmodel/=1) then
            call t_divscheme%evaladdrhs(tdd_globalrhs(:,step))
      end if
      ! add right-hand sides
      if(staticalphamodel==1 .and. staticucmodel==1) then
            ! do nothing
      else if(staticalphamodel/=1 .and. staticucmodel/=1) then
            ! do nothing
      else
            call fvl_omp_add(tdd_globalrhs(:,step),tdd_staticrhs(:,step))
      end if

      ! clean memory
      if(staticftmodel/=1) then
            call fvl_omp_zeros(tdd_globalsource(:,step))
      end if
      ! evaluate heat source term
      if(staticftmodel/=1) then
            call ftdd_model%evalcellquadrattomeanvalues(tdd_globalsource(:,step))
      end if
      ! ! add source terms
      ! if(.not. (staticftmodel/=1)) then
      !       call fvl_omp_add(tdd_globalsource(:,step),tdd_staticsource(:,step))
      ! end if

      ! ! compute preconditioner
      ! ! if(staticalphamodel/=1 .or. staticucmodel/=1) then
      !       if(t_preconditioner%iterate()) then
      !             call t_preconditioner%evalprecond(t_globalmat_comp,t_precondmat)
      !       end if
      ! ! end if

      ! stop timer
      call fvl_stopwtimer(55)
end subroutine t_ufun3

! ============================================================================
! MOMENTUM BALANCE AND MASS CONSERVATION EQUATION
! ============================================================================

subroutine pu_tfun(step,x,res)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:size(x,1)),intent(out)::res
      ! start timer
      call fvl_startwtimer(51,"solution_tfun")
      ! momentum equation
      if(ux_integmat_comp%getmaxnumrowvalues()/=0) then
            call ux_integmat_comp%blockmul2(x(1:u_dofs),res(1:u_dofs))
      else
            call fvl_omp_assign(x(1:u_dofs),res(1:u_dofs))
      end if
      ! continuity equation
      call fvl_omp_assign(0.0d0,res(u_dofs+1:u_dofs+p_dofs))
      ! compatibility equation
      res(u_dofs+p_dofs+1)=0.0d0
      ! stop timer
      call fvl_stopwtimer(51)
end subroutine pu_tfun

subroutine pu_afun(step,x,res)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:size(x,1)),intent(out)::res
      ! start timer
      call fvl_startwtimer(52,"solution_afun")
      ! momentum equation
      call ux_globalmat_comp1%blockmul2(x(1:u_dofs),res(1:u_dofs))
      call ux_globalmat_comp2%addrowblockmul2(x(1:u_dofs),res(1:u_dofs))
      call p_globalmat_comp%addcolblockmul2(x(u_dofs+1:u_dofs+p_dofs),res(1:u_dofs))
      ! continuity equation
      call u_globalmat_comp%rowblockmul2(x(1:u_dofs),res(u_dofs+1:u_dofs+p_dofs))
      call fvl_omp_add(res(u_dofs+1:u_dofs+p_dofs),-x(u_dofs+p_dofs+1))
      ! compatibility equation
      res(u_dofs+p_dofs+1)=fvl_omp_dot(x(u_dofs+1:u_dofs+p_dofs),mesh%getcellareas())
      ! stop timer
      call fvl_stopwtimer(52)
end subroutine pu_afun

subroutine pu_bfun(step,res)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(inout)::res
      ! start timer
      call fvl_startwtimer(53,"solution_afun")
      ! momentum equation
      call fvl_omp_add(pu_globalrhs,reshape(pu_globalsource,[u_dofs]),res(1:u_dofs))
      ! continuity equation
      call fvl_omp_add(u_globalrhs,u_globalsource,res(u_dofs+1:u_dofs+p_dofs))
      ! compatibility equation
      res(u_dofs+p_dofs+1)=0.0d0
      ! stop timer
      call fvl_stopwtimer(53)
end subroutine pu_bfun

subroutine pu_pfun(step,x,res)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      real(kind=__fvl_real_kind__),dimension(1:),contiguous,intent(in)::x
      real(kind=__fvl_real_kind__),dimension(1:size(x,1)),intent(out)::res
      ! start timer
      call fvl_startwtimer(54,"solution_pfun")
      ! momentum equation
      if(ux_precondmat%getmaxnumrowvalues()/=0) then
            call ux_precondmat%blockmul2(x(1:u_dofs),res(1:u_dofs))
      else
            call fvl_omp_assign(x(1:u_dofs),res(1:u_dofs))
      end if
      ! continuity equation
      if(p_precondmat%getmaxnumrowvalues()/=0) then
            call p_precondmat%vecmul(x(u_dofs+1:u_dofs+p_dofs),res(u_dofs+1:u_dofs+p_dofs))
      else
            call fvl_omp_assign(x(u_dofs+1:u_dofs+p_dofs),res(u_dofs+1:u_dofs+p_dofs))
      end if
      ! compatibility equation
      res(u_dofs+p_dofs+1)=x(u_dofs+p_dofs+1)
      ! stop timer
      call fvl_stopwtimer(54)
end subroutine pu_pfun

subroutine pu_updatetime()
      ! start timer
      call fvl_startwtimer(55,"solution_ufun")

      ! stop timer
      call fvl_stopwtimer(55)
end subroutine pu_updatetime

subroutine pu_ufun(step)
      integer(kind=__fvl_integer_kind__),intent(in)::step
      integer(kind=__fvl_integer_kind__)::i
      class(fvl_vectorfield2d),pointer::field

      ! start timer
      call fvl_startwtimer(55,"solution_ufun")

      ! evaluate pressure gradient coefficients matrix
      if(staticrhodmodel/=1) then
            call p_gradscheme%evalmatrix(p_globalmat)
      end if
      ! compress coefficients matrices
      if(staticrhodmodel/=1) then
            call p_globalmat_comp%compress(p_globalmat)
            call fvl_delete(p_globalmat)
      end if

      ! update velocity field
      if(staticvcmodel/=1) then
            ! interpolate velocity field
            call ux_interpolate%evalvalues(u_model,vc_model)
            ! update velocity divergence scheme
            call ux_divscheme%makeupwind() !!!!!
      end if

      ! update temperature field
      if(staticfpumodel/=1) then
            ! interpolate temperature field
            call t_interpolate%evalvalues(t_model,tc_model)
            ! update source term
            do i=1,fpu_model%getnumfields()
                  field=>fpu_model%getfield(i)
                  select type(field)
                        class is(fvl_boussinesq_forcesource2d)
                              call field%update(tc_model,g_model,beta_model)
                  end select
            end do
      end if

      ! evaluate velocity laplacian coefficients matrix
      if(staticnumodel/=1) then
            call ux_lapscheme%evalmatrix(ux_globalmat1,ux_globalmat2)
      end if
      ! evaluate velocity divergence coefficients matrix
      if(staticvcmodel/=1) then
            call ux_divscheme%evaladdmatrix(ux_globalmat1,ux_globalmat2)
      end if
      ! add and compress coefficients matrices
      if(staticnumodel==1 .and. staticvcmodel==1) then
            ! do nothing
      else if(staticnumodel/=1 .and. staticvcmodel/=1) then
            call ux_globalmat_comp1%compress(ux_globalmat1)
            call fvl_delete(ux_globalmat1)
            call ux_globalmat_comp2%compress(ux_globalmat2)
            call fvl_delete(ux_globalmat2)
      else
            call ux_globalmat1%addmatrix(ux_staticmat_comp1)
            call ux_globalmat2%addmatrix(ux_staticmat_comp2)
            call ux_globalmat_comp1%compress(ux_globalmat1)
            call fvl_delete(ux_globalmat1)
            call ux_globalmat_comp2%compress(ux_globalmat2)
            call fvl_delete(ux_globalmat2)
      end if

      ! clean memory
      if(staticnumodel/=1 .or. staticvcmodel/=1) then
            call fvl_omp_zeros(pu_globalrhs)
      end if
      ! evaluate pressure gradient right-hand side
      if(staticrhodmodel/=1) then
            call p_gradscheme%evalrhs(pu_globalrhs)
      end if
      ! evaluate velocity laplacian right-hand side
      if(staticnumodel/=1) then
            call ux_lapscheme%evaladdrhs(pu_globalrhs)
      end if
      ! evaluate velocity divergence right-hand side
      if(staticvcmodel/=1) then
            call ux_divscheme%evaladdrhs(pu_globalrhs)
      end if
      ! add right-hand sides
      if(staticrhodmodel==1 .and. staticnumodel==1 .and. staticvcmodel==1) then
            ! do nothing
      else if(staticrhodmodel/=1 .and. staticnumodel/=1 .and. staticvcmodel/=1) then
            ! do nothing
      else
            call fvl_omp_add(pu_globalrhs,pu_staticrhs)
      end if

      ! clean memory
      if(staticfpumodel/=1) then
            call fvl_omp_zeros(pu_globalsource)
      end if
      ! evaluate momentum source term
      if(staticfpumodel/=1) then
            call fpu_model%evalcellquadrattomeanvalues(pu_globalsource)
      end if
      ! ! add source terms
      ! if(.not. (staticfpumodel/=1)) then
      !       call fvl_omp_add(pu_globalsource,pu_staticsource)
      ! end if

      ! clean memory
      if(staticgpumodel/=1) then
            call fvl_omp_zeros(u_globalsource)
      end if
      ! evaluate mass source term
      if(staticgpumodel/=1) then
            call gpu_model%evalcellquadrattomeanvalues(u_globalsource)
      end if
      ! ! add source terms
      ! if(.not. (staticgpumodel/=1)) then
      !       call fvl_omp_add(pu_globalsource,pu_staticsource)
      ! end if

      ! compute preconditioner
      ! if(staticrhodmodel/=1 .or. staticnumodel/=1 .or. staticvcmodel/=1) then
            if(pu_preconditioner%iterate()) then
                  call pu_preconditioner%evalschurprecond(ux_globalmat_comp1,p_globalmat_comp,u_globalmat_comp,ux_precondmat,p_precondmat)
            end if
      ! end if

      ! stop timer
      call fvl_stopwtimer(55)
end subroutine pu_ufun

subroutine writetime()

      ! state message
      call fvl_loginfo(">>> Writing fields")

      ! write fields
#include "fvl_writefields_src.f"

      ! write reports
#include "fvl_writereports_src.f"

      ! write histograms
      call t_linearsolver%writehistogram(fvl_trim(currenttimedirectory)//"/t_linearsolver.fvd","ascii")
      !call pu_linearsolver%writehistogram(fvl_trim(currenttimedirectory)//"/pu_linearsolver.fvd","ascii")
      call t_fixpointiterator%writehistogram(fvl_trim(currenttimedirectory)//"/t_fixpointiterator.fvd","ascii")
      !call p_fixpointiterator%writehistogram(fvl_trim(currenttimedirectory)//"/p_fixpointiterator.fvd","ascii") !!!!!
      !call u_fixpointiterator%writehistogram(fvl_trim(currenttimedirectory)//"/u_fixpointiterator.fvd","ascii") !!!!!

end subroutine writetime


function fvl_scalarmodel2d_init44(mesh,fieldtype,modeldictfile,modeldictlabel,patchesdictfile,patchesdictlabel,&
      fieldselector,boundcondselector,label,filepath,fileform,calculatedfields) result(res)
      class(fvl_mesh2d),intent(in),target::mesh
      integer(kind=__fvl_integer_kind__),intent(in)::fieldtype
      class(fvl_dict_file),intent(in)::modeldictfile
      character(len=*),intent(in)::modeldictlabel
      class(fvl_dict_file),intent(in)::patchesdictfile
      character(len=*),intent(in)::patchesdictlabel
      class(fvl_scalarfield2d_selector)::fieldselector
      class(fvl_scalarboundcond2d_selector)::boundcondselector
      character(len=*),intent(in),optional::label
      character(len=*),intent(in),optional::filepath
      character(len=*),intent(in),optional::fileform
      logical(kind=__fvl_logical_kind__),optional::calculatedfields
      type(fvl_scalarmodel2d)::res
      character(len=__fvl_character_len__)::type
      integer(kind=__fvl_integer_kind__)::i,numdicts,numcodes,patchtype
      character(len=__fvl_character_len__),dimension(1:fvl_scalarmodel2d_maxsize)::dictlabels
      integer(kind=__fvl_integer_kind__),dimension(1:fvl_scalarmodel2d_maxnumcodes)::codes
      type(fvl_dict_file)::modeldict,patchesdict,fieldsdict,boundcondsdict
      class(fvl_patch2d),pointer::patch
      class(fvl_scalarfield2d),pointer::field1,field2
      class(fvl_scalarboundcond2d),pointer::boundcond
#if(__fvl_debug__==1)
      if(fieldtype/=fvl_field2d_vertcentresfield .and. fieldtype/=fvl_field2d_vertboundsfield &
            .and. fieldtype/=fvl_field2d_edgecentresfield .and. fieldtype/=fvl_field2d_edgemeansfield&
            .and. fieldtype/=fvl_field2d_edgequadratsfield .and. fieldtype/=fvl_field2d_edgeboundsfield&
            .and. fieldtype/=fvl_field2d_cellcentresfield .and. fieldtype/=fvl_field2d_cellmeansfield&
            .and. fieldtype/=fvl_field2d_cellquadratsfield .and. fieldtype/=fvl_field2d_cellboundsfield&
            .and. fieldtype/=fvl_field2d_cellboundsfield) then
            __fvl_log_error__("Input error in call fvl_scalarmodel2d_init4: invalid argument value for <fieldtype>")
      end if
      if(present(fileform)) then
            if(fvl_trim(fileform)/="ascii" .and. fvl_trim(fileform)/="binary") then
                  __fvl_log_error__("Input error in call fvl_scalarmodel2d_init4: invalid argument value for <fileform>")
            end if
      end if
#endif
      ! initialize data
      if(present(label)) then
            res%label=label
      else
            res%label=fvl_scalarmodel2d_label_default
      end if
      if(present(filepath)) then
            res%filepath=filepath
      else
            res%filepath=fvl_scalarmodel2d_filepath_default
      end if
      if(present(fileform)) then
            res%fileform=fileform
      else
            res%fileform=fvl_scalarmodel2d_fileform_default
      end if
      if(present(calculatedfields)) then
            res%calculatedfields=calculatedfields
      else
            res%calculatedfields=fvl_scalarmodel2d_calculatedfields_default
      end if
      res%mesh=>mesh
      if(fieldtype==fvl_field2d_vertcentresfield) then
            patchtype=fvl_patch2d_vertspatch
      else if(fieldtype==fvl_field2d_vertboundsfield) then
            patchtype=fvl_patch2d_boundvertspatch
      else if(fieldtype==fvl_field2d_edgecentresfield .or. fieldtype==fvl_field2d_edgemeansfield&
            .or. fieldtype==fvl_field2d_edgequadratsfield) then
            patchtype=fvl_patch2d_edgespatch
      else if(fieldtype==fvl_field2d_edgeboundsfield) then
            patchtype=fvl_patch2d_boundedgespatch
      else if(fieldtype==fvl_field2d_cellcentresfield .or. fieldtype==fvl_field2d_cellmeansfield&
            .or. fieldtype==fvl_field2d_cellquadratsfield) then
            patchtype=fvl_patch2d_cellspatch
      else if(fieldtype==fvl_field2d_cellboundsfield) then
            patchtype=fvl_patch2d_boundcellspatch
      else if(fieldtype==fvl_field2d_boundsfield) then
            patchtype=fvl_patch2d_boundspatch
      end if
      res%numfieldpatches=0
      res%numboundcondpatches=0
      res%numfields=0
      res%numboundconds=0
      ! allocate memory
      call fvl_allocate(res%fieldpatches,fvl_scalarmodel2d_maxsize)
      call fvl_allocate(res%boundcondpatches,fvl_scalarmodel2d_maxsize)
      call fvl_allocate(res%fields,fvl_scalarmodel2d_maxsize)
      call fvl_allocate(res%boundconds,fvl_scalarmodel2d_maxsize)
      ! extract model dictionary
      call modeldictfile%extractdict(modeldictlabel,modeldict)
      ! extract patch dictionary
      call patchesdictfile%extractdict(patchesdictlabel,patchesdict)
      ! extract field dictionaries
      call modeldict%extractdict("fields",fieldsdict)
      numdicts=fieldsdict%extractdictlabels(dictlabels)
      ! initialize fields



      do i=1,numdicts
            ! get patch codes
            numcodes=patchesdict%getarray(dictlabels(i),"codes",0,codes)
            ! create patch
            call fvl_setnumthreads(4)

            res%numfieldpatches=res%numfieldpatches+1
            allocate(patch,source=fvl_patch2d(mesh,patchtype,codes(1:numcodes),dictlabels(i)))
            call res%fieldpatches(res%numfieldpatches)%setptr(patch)

            write(*,*) patch%indices

            call fvl_setnumthreads(4)

            ! get field type
            type=fieldsdict%getvalue(dictlabels(i),"type","")
            ! create field
            res%numfields=res%numfields+1
            allocate(field1,source=fieldselector%allocate(type,patch,fieldtype,res%label,fieldsdict,dictlabels(i)))
            call res%fields(res%numfields)%setptr(field1)
      end do


      ! extract boundconds dictionary
      call modeldict%extractdict("bound_conds",boundcondsdict)
      numdicts=boundcondsdict%extractdictlabels(dictlabels)
      ! initialize boundconds
      do i=1,numdicts
            ! get patch codes
            numcodes=patchesdict%getarray(dictlabels(i),"codes",0,codes)
            ! create patch
            res%numboundcondpatches=res%numboundcondpatches+1
            allocate(patch,source=fvl_patch2d(mesh,fvl_patch2d_boundedgespatch,codes(1:numcodes),dictlabels(i)))
            call res%boundcondpatches(res%numboundcondpatches)%setptr(patch)
            ! get boundcond type
            type=boundcondsdict%getvalue(dictlabels(i),"type","")
            ! create boundcond
            res%numboundconds=res%numboundconds+1
            allocate(boundcond,source=boundcondselector%allocate(type,patch,res%label,boundcondsdict,dictlabels(i)))
            call res%boundconds(res%numboundconds)%setptr(boundcond)
      end do
      ! read calculated fields
      call res%read()
      ! convert to calculated
      if(res%calculatedfields) then
            do i=1,res%numfields
                  field1=>res%fields(i)%getptr()
                  select type(field1)
                        class is(fvl_calc_scalarfield2d)
                              ! do nothing
                        class default
                              allocate(field2,source=fvl_calc_scalarfield2d(field1))
                              call res%fields(i)%setptr(field2)
                              deallocate(field1)
                  end select
            end do
      end if
end function fvl_scalarmodel2d_init44

end program fvl_pro_thermoincompflow2d_main
! end of file
