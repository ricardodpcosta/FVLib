! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: check and analyse meshes
! Modification: March, 2025

#include "macros.f"

program fvl_meshreport_main

use fvl_lib

implicit none

! ============================================================================
! VARIABLES
! ============================================================================

character(len=__fvl_character_len__)::meshinpath,meshinform,meshinformat
type(fvl_file)::meshinfile

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
      call fvl_loginfo("About: Utility tool for mesh inspection and reporting.")
      call fvl_loginfo("Usage: fvl-meshreport <mesh_file> [options]")
      call fvl_loginfo("Arguments:")
      call fvl_loginfo("      <input_mesh_file>       Input mesh file in msh/fvm1[.bin]/fvm2[.bin]/fvm3[.bin] format.")
      call fvl_loginfo("Options:")
      call fvl_loginfo("      -h|--help               Displays help message.")
      call fvl_loginfo("      -p|--parallel           Runs in parallel with <OMP_NUM_THREADS> processes (default: serial).")
      stop
end if

! check arguments
if(fvl_getnumopts()<1) then
      call fvl_logerror("Missing arguments")
end if

! get arguments
meshinpath=adjustr(fvl_getarg(1))

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
else if(meshinpath(__fvl_character_len__-4:__fvl_character_len__)==".fvm1") then
      meshinformat="fvm1"
else if(meshinpath(__fvl_character_len__-4:__fvl_character_len__)==".fvm2") then
      meshinformat="fvm2"
else if(meshinpath(__fvl_character_len__-4:__fvl_character_len__)==".fvm3") then
      meshinformat="fvm3"
else
      call fvl_logerror("Invalid input mesh format")
end if

! initialize files
meshinpath=adjustl(fvl_getarg(1))
meshinfile=fvl_file(path=meshinpath,form=meshinform,action="read")

! report mesh
if(fvl_trim(meshinformat)=="msh") then
      call fvl_gmshreport3d(meshinfile)
else if(fvl_trim(meshinformat)=="fvm1") then
      call fvl_meshreport1d(meshinfile)
else if(fvl_trim(meshinformat)=="fvm2") then
      call fvl_meshreport2d(meshinfile)
else if(fvl_trim(meshinformat)=="fvm3") then
      call fvl_meshreport3d(meshinfile)
else
      call fvl_logerror("Invalid input mesh format")
end if

! ============================================================================
! STOP LIBRARY
! ============================================================================

call fvl_stop()

end program fvl_meshreport_main
! end of file
