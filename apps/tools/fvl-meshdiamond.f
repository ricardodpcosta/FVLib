! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: generate staggered dimamond meshes
! Modification: March, 2025

#include "macros.f"

program fvl_meshdiamond_main

use fvl_lib

implicit none

! ============================================================================
! VARIABLES
! ============================================================================

character(len=__fvl_character_len__)::meshinpath,meshoutpath,meshinform,meshoutform,meshinformat,meshoutformat
type(fvl_file)::meshinfile,meshoutfile

! ============================================================================
! INITIALIZE LIBRARY
! ============================================================================

! initialize global parameters
call fvl_init()

! check arguments
if(.not. fvl_checkopts("-h","--help","-p","--parallel"," ")) then
      call fvl_logerror("Invalid option.")
end if

! check help option
if(fvl_findopt("-h") .or. fvl_findopt("--help")) then
      call fvl_loginfo("About: Utility tool for creation of dimamond meshes.")
      call fvl_loginfo("Usage: fvl-meshdiamond <input_mesh_file> <output_mesh_file> [options]")
      call fvl_loginfo("Arguments:")
      call fvl_loginfo("      <input_mesh_file>       Input mesh file in msh/fvm2[.bin]/fvm3[.bin] format.")
      call fvl_loginfo("      <output_mesh_file>      Output mesh file (same format as the input mesh).")
      call fvl_loginfo("Options:")
      call fvl_loginfo("      -h|--help               Displays help message.")
      call fvl_loginfo("      -p|--parallel           Runs in parallel with <OMP_NUM_THREADS> processes (default: serial).")
      stop
end if

! check arguments
if(fvl_getnumopts()<2) then
      call fvl_logerror("Missing arguments")
end if

! get arguments
meshinpath=adjustr(fvl_getarg(1))
meshoutpath=adjustr(fvl_getarg(2))

! ============================================================================
! RUN UTILITY TOOL
! ============================================================================

! check input mesh form
if(meshinpath(__fvl_character_len__-3:__fvl_character_len__)==".bin") then
      meshinform="binary"
      meshinpath=meshinpath(1:__fvl_character_len__-4)
      meshinpath=adjustr(meshinpath)
else
      meshinform="ascii"
end if

! check input mesh format
if(meshinpath(__fvl_character_len__-3:__fvl_character_len__)==".msh") then
      meshinformat="msh"
else if(meshinpath(__fvl_character_len__-4:__fvl_character_len__)==".fvm2") then
      meshinformat="fvm2"
else if(meshinpath(__fvl_character_len__-4:__fvl_character_len__)==".fvm3") then
      meshinformat="fvm3"
else
      call fvl_logerror("Invalid input mesh format")
end if

! check output mesh form
if(meshoutpath(__fvl_character_len__-3:__fvl_character_len__)==".bin") then
      meshoutform="binary"
      meshoutpath=meshoutpath(1:__fvl_character_len__-4)
      meshoutpath=adjustr(meshoutpath)
else
      meshoutform="ascii"
end if

! check output mesh format
if(meshoutpath(__fvl_character_len__-3:__fvl_character_len__)==".msh") then
      meshoutformat="msh"
else if(meshoutpath(__fvl_character_len__-4:__fvl_character_len__)==".fvm2") then
      meshoutformat="fvm2"
else if(meshoutpath(__fvl_character_len__-4:__fvl_character_len__)==".fvm3") then
      meshoutformat="fvm3"
else
      call fvl_logerror("Invalid output mesh format")
end if

! initialize files
meshinpath=adjustl(fvl_getarg(1))
meshoutpath=adjustl(fvl_getarg(2))
meshinfile=fvl_file(path=meshinpath,form=meshinform,action="read")
meshoutfile=fvl_file(path=meshoutpath,form=meshoutform,action="write")

! create diamond mesh
if(fvl_trim(meshinformat)=="fvm2") then
      if(fvl_trim(meshoutformat)=="fvm2") then
            call fvl_loginfo("Generating diamond mesh from "//fvl_trim(meshinpath)//" to "//fvl_trim(meshoutpath)//"...")
            call fvl_meshdiamond2d(meshinfile,meshoutfile)
      else
            call fvl_logerror("Invalid input mesh format")
      end if
else if(fvl_trim(meshinformat)=="fvm3") then
      if(fvl_trim(meshoutformat)=="fvm3") then
            call fvl_loginfo("Creating diamond mesh from "//fvl_trim(meshinpath)//" to "//fvl_trim(meshoutpath)//"...")
            call fvl_meshdiamond3d(meshinfile,meshoutfile)
      else
            call fvl_logerror("Invalid input mesh format")
      end if
else if(fvl_trim(meshinformat)=="msh") then
      if(fvl_trim(meshoutformat)=="msh") then
            call fvl_loginfo("Creating diamond mesh from "//fvl_trim(meshinpath)//" to "//fvl_trim(meshoutpath)//"...")
            call fvl_gmshdiamond3d(meshinfile,meshoutfile)
      else
            call fvl_logerror("Invalid input mesh format")
      end if
else
      call fvl_logerror("Invalid output mesh format")
end if

! ============================================================================
! STOP LIBRARY
! ============================================================================

call fvl_stop()

end program fvl_meshdiamond_main
! end of file
