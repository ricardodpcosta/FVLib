! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: search vortex centre
! Modification: March, 2025

#include "macros.f"

program fvl_vortexcentre_main

use fvl_lib

implicit none

! ============================================================================
! VARIABLES
! ============================================================================

character(len=__fvl_character_len__)::meshinpath,meshoutpath,datainpath,dataoutpath,meshinform,meshoutform,datainform,dataoutform,&
      meshinformat,meshoutformat,datainformat,dataoutformat,dictfilepath
type(fvl_file)::meshinfile,meshoutfile
type(fvl_data_file)::datainfile,dataoutfile
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
      call fvl_loginfo("About: Utility tool for searching vortex centre.")
      call fvl_loginfo("Usage: fvl-vortexcentre <input_mesh_file> <input_data_file> <output_mesh_file> <output_data_file> [options]")
      call fvl_loginfo("Arguments:")
      call fvl_loginfo("      <input_mesh_file>       Input mesh file in msh/fvm2[.bin]/fvm3[.bin] format.")
      call fvl_loginfo("      <input_data_file>       Input data file in msh/fvd2[.bin]/fvd3[.bin] format.")
      call fvl_loginfo("      <output_mesh_file>      Output mesh file in msh/fvm2[.bin]/fvm3[.bin] format.")
      call fvl_loginfo("      <output_data_file>      Output data file in msh/fvd2[.bin]/fvd3[.bin] format.")
      call fvl_loginfo("Options:")
      call fvl_loginfo("      -h|--help               Displays help message.")
      call fvl_loginfo("      -d|--dict <file>        Reads parameters from dictionary file <file> (default: setup/tools.fvd).")
      call fvl_loginfo("      -p|--parallel           Runs in parallel with <OMP_NUM_THREADS> processes (default: serial).")
      stop
end if

! check arguments
if(fvl_getnumopts()<4) then
      call fvl_logerror("Missing arguments")
end if

! get arguments
meshinpath=adjustr(fvl_getarg(1))
datainpath=adjustr(fvl_getarg(2))
meshoutpath=adjustr(fvl_getarg(3))
dataoutpath=adjustr(fvl_getarg(4))

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

! check input data form
if(datainpath(__fvl_character_len__-3:__fvl_character_len__)==".bin") then
      datainform="binary"
      datainpath=datainpath(1:__fvl_character_len__-4)
      datainpath=adjustr(datainpath)
else
      datainform="ascii"
end if

! check input data format
! if(datainpath(__fvl_character_len__-3:__fvl_character_len__)==".msh") then
!       datainformat="msh"
! else
if(datainpath(__fvl_character_len__-4:__fvl_character_len__)==".fvd2") then
      datainformat="fvd2"
else if(datainpath(__fvl_character_len__-4:__fvl_character_len__)==".fvd3") then
      datainformat="fvd3"
else
      call fvl_logerror("Invalid input data format")
end if

! check output data form
if(dataoutpath(__fvl_character_len__-3:__fvl_character_len__)==".bin") then
      dataoutform="binary"
      dataoutpath=dataoutpath(1:__fvl_character_len__-4)
      dataoutpath=adjustr(dataoutpath)
else
      dataoutform="ascii"
end if

! check output data format
! if(dataoutpath(__fvl_character_len__-3:__fvl_character_len__)==".msh") then
!       dataoutformat="msh"
! else
if(dataoutpath(__fvl_character_len__-4:__fvl_character_len__)==".fvd2") then
      dataoutformat="fvd2"
else if(dataoutpath(__fvl_character_len__-4:__fvl_character_len__)==".fvd3") then
      dataoutformat="fvd3"
else
      call fvl_logerror("Invalid output data format")
end if

! initialize files
meshinpath=adjustl(fvl_getarg(1))
datainpath=adjustl(fvl_getarg(2))
meshoutpath=adjustl(fvl_getarg(3))
dataoutpath=adjustl(fvl_getarg(4))
meshinfile=fvl_file(path=meshinpath,form=meshinform,action="read")
datainfile=fvl_data_file(path=datainpath,form=datainform,action="read")
meshoutfile=fvl_file(path=meshoutpath,form=meshoutform,action="write")
dataoutfile=fvl_data_file(path=dataoutpath,form=dataoutform,action="write")
dictfile=fvl_dict_file(path=dictfilepath,form="ascii",action="read")

! read dictionary
if(dictfile%exist()) then
      call dictfile%open()
      call dictfile%readdicts()
      call dictfile%close()
else
      call fvl_logerror("Dictionary file not found")
end if

! compute vortexcentre
call fvl_loginfo("Computing vortexcentre from "//fvl_trim(datainpath)//" to "//fvl_trim(dataoutpath)//"...")
if(fvl_trim(meshinformat)=="fvm2") then
      if(fvl_trim(meshoutformat)=="fvm2") then
            if(fvl_trim(datainformat)=="fvd2") then
                  if(fvl_trim(dataoutformat)=="fvd2") then
                        call fvl_vortexcentre2d(meshinfile,datainfile,meshoutfile,dataoutfile,dictfile)
                  else
                        call fvl_logerror("Invalid output data format")
                  end if
            else
                  call fvl_logerror("Invalid input data format")
            end if
      else
            call fvl_logerror("Invalid output mesh format")
      end if
else if(fvl_trim(meshinformat)=="fvm3") then
      if(fvl_trim(meshoutformat)=="fvm3") then
            if(fvl_trim(datainformat)=="fvd3") then
                  if(fvl_trim(dataoutformat)=="fvd3") then
                        call fvl_vortexcentre3d(meshinfile,datainfile,meshoutfile,dataoutfile,dictfile)
                  else
                        call fvl_logerror("Invalid output data format")
                  end if
            else
                  call fvl_logerror("Invalid input data format")
            end if
      else
            call fvl_logerror("Invalid output mesh format")
      end if
else
      call fvl_logerror("Invalid input mesh format")
end if

! ============================================================================
! STOP LIBRARY
! ============================================================================

call fvl_stop()

end program fvl_vortexcentre_main
! end of file
