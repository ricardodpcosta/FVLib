! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: deform mesh
! Modification: March, 2025

#include "macros.f"

program fvl_meshdeform_main

use fvl_lib

implicit none

! ============================================================================
! VARIABLES
! ============================================================================

character(len=__fvl_character_len__)::deformation,meshinpath,meshoutpath,meshinform,meshoutform,meshinformat,meshoutformat,dictfilepath
type(fvl_file)::meshinfile,meshoutfile
type(fvl_dict_file)::dictfile

! ============================================================================
! INITIALIZE LIBRARY
! ============================================================================

! initialize global parameters
call fvl_init()

! check arguments
if(.not. fvl_checkopts("-h","--help","-d","--dict","-p","--parallel"," ")) then
      call fvl_logerror("Invalid option.")
end if

! check help option
if(fvl_findopt("-h") .or. fvl_findopt("--help")) then
      call fvl_loginfo("About: Utility tool for mesh deformation.")
      call fvl_loginfo("Usage: fvl-meshdeform <deformation> <input_mesh_file> <output_mesh_file> [options]")
      call fvl_loginfo("Arguments:")
      call fvl_loginfo("      <deformation>           Deformation to perform (available: perturb, ellipsoid, popcorn, asteroid).")
      call fvl_loginfo("      <input_mesh_file>       Input mesh file in msh/fvm2[.bin]/fvm3[.bin] format.")
      call fvl_loginfo("      <output_mesh_file>      Output mesh file in msh/fvm2[.bin]/fvm3[.bin] format.")
      call fvl_loginfo("Options:")
      call fvl_loginfo("      -h|--help               Displays help message.")
      call fvl_loginfo("      -d|--dict <file>        Reads parameters from dictionary file <file> (default: setup/tools.fvd).")
      call fvl_loginfo("      -p|--parallel           Runs in parallel with <OMP_NUM_THREADS> processes (default: serial).")
      stop
end if

! check arguments
if(fvl_getnumopts()<3) then
      call fvl_logerror("Missing arguments")
end if

! get arguments
deformation=adjustr(fvl_getarg(1))
meshinpath=adjustr(fvl_getarg(2))
meshoutpath=adjustr(fvl_getarg(3))

! check param file option
if(fvl_findopt("-d",dictfilepath) .or. fvl_findopt("--dict",dictfilepath)) then
      if(dictfilepath==" ") then
            call fvl_logerror("Missing argument for option.")
      end if
else
      dictfilepath = "setup/tools.fvd"
end if

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
meshinpath=adjustl(fvl_getarg(2))
meshoutpath=adjustl(fvl_getarg(3))
meshinfile=fvl_file(path=meshinpath,form=meshinform,action="read")
meshoutfile=fvl_file(path=meshoutpath,form=meshoutform,action="write")
dictfile=fvl_dict_file(path=dictfilepath,form="ascii",action="read")

! read dictionary
if(dictfile%exist()) then
      call dictfile%open()
      call dictfile%readdicts()
      call dictfile%close()
else
      call fvl_logerror("Dictionary file not found")
end if

! create deformed mesh
call fvl_loginfo("Creating deformed mesh from "//fvl_trim(meshinpath)//" to "//fvl_trim(meshoutpath)//"...")
if(fvl_trim(meshinformat)=="fvm2") then
      if(fvl_trim(meshoutformat)=="fvm2") then
            if(fvl_trim(deformation)=="perturb") then
                  call fvl_meshperturb2d(meshinfile,meshoutfile,dictfile)
            else
                  call fvl_logerror("Invalid mesh deformation")
            end if
      else
            call fvl_logerror("Invalid input mesh format")
      end if
else if(fvl_trim(meshinformat)=="fvm3") then
      if(fvl_trim(meshoutformat)=="fvm3") then
            if(fvl_trim(deformation)=="perturb") then
                  call fvl_meshperturb3d(meshinfile,meshoutfile,dictfile)
            else if(fvl_trim(deformation)=="ellipsoid") then
                  call fvl_meshellipsoid3d(meshinfile,meshoutfile,dictfile)
            else if(fvl_trim(deformation)=="popcorn") then
                  call fvl_meshpopcorn3d(meshinfile,meshoutfile,dictfile)
            else if(fvl_trim(deformation)=="asteroid") then
                  call fvl_meshasteroid3d(meshinfile,meshoutfile,dictfile)
            else
                  call fvl_logerror("Invalid mesh deformation")
            end if
      else
            call fvl_logerror("Invalid input mesh format")
      end if
else if(fvl_trim(meshinformat)=="msh") then
      if(fvl_trim(meshoutformat)=="msh") then
            if(fvl_trim(deformation)=="perturb") then
                  call fvl_gmshperturb3d(meshinfile,meshoutfile,dictfile)
            else if(fvl_trim(deformation)=="ellipsoid") then
                  call fvl_gmshellipsoid3d(meshinfile,meshoutfile,dictfile)
            else if(fvl_trim(deformation)=="popcorn") then
                  call fvl_gmshpopcorn3d(meshinfile,meshoutfile,dictfile)
            else if(fvl_trim(deformation)=="asteroid") then
                  call fvl_gmshasteroid3d(meshinfile,meshoutfile,dictfile)
            else
                  call fvl_logerror("Invalid mesh deformation")
            end if
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

end program fvl_meshdeform_main
! end of file
