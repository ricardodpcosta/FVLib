! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: operate data
! Modification: March, 2025

#include "macros.f"

program fvl_dataoperate_main

use fvl_lib

implicit none

! ============================================================================
! VARIABLES
! ============================================================================

character(len=__fvl_character_len__)::operation,datainpath,dataoutpath,datainform,dataoutform,datainformat,dataoutformat
type(fvl_data_file)::datainfile,dataoutfile

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
      call fvl_loginfo("About: Utility tool for data operations.")
      call fvl_loginfo("Usage: fvl-dataoperate <operation> <input_data_file> <output_data_file> [options]")
      call fvl_loginfo("Arguments:")
      call fvl_loginfo("      <operation>             Operation to perform (available: magfield, sqrfield, perpfield).")
      call fvl_loginfo("      <input_data_file>       Input data file in msh/fvd2[.bin]/fvd3[.bin] format.")
      call fvl_loginfo("      <output_data_file>      Output data file in msh/fvd2[.bin]/fvd3[.bin] format.")
      call fvl_loginfo("Options:")
      call fvl_loginfo("      -h|--help               Displays help message.")
      call fvl_loginfo("      -p|--parallel           Runs in parallel with <OMP_NUM_THREADS> processes (default: serial).")
      stop
end if

! check arguments
if(fvl_getnumopts()<3) then
      call fvl_logerror("Missing arguments")
end if

! get arguments
operation=adjustr(fvl_getarg(1))
datainpath=adjustr(fvl_getarg(2))
dataoutpath=adjustr(fvl_getarg(3))

! ============================================================================
! RUN UTILITY TOOL
! ============================================================================

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
datainpath=adjustl(fvl_getarg(2))
dataoutpath=adjustl(fvl_getarg(3))
datainfile=fvl_data_file(path=datainpath,form=datainform,action="read")
dataoutfile=fvl_data_file(path=dataoutpath,form=dataoutform,action="write")

! compute dataoperate
if(fvl_trim(datainformat)=="fvd2") then
      if(fvl_trim(dataoutformat)=="fvd2") then
            if(fvl_trim(operation)=="magfield") then
                  call fvl_loginfo("Computing magfield from "//fvl_trim(datainpath)//" to "//fvl_trim(dataoutpath)//"...")
                  call fvl_datamagfield2d(datainfile,dataoutfile)
            else if(fvl_trim(operation)=="sqrfield") then
                  call fvl_loginfo("Computing sqrfield from "//fvl_trim(datainpath)//" to "//fvl_trim(dataoutpath)//"...")
                  call fvl_datasqrfield2d(datainfile,dataoutfile)
            else if(fvl_trim(operation)=="perpfield") then
                  call fvl_loginfo("Computing perpfield from "//fvl_trim(datainpath)//" to "//fvl_trim(dataoutpath)//"...")
                  call fvl_dataperpfield2d(datainfile,dataoutfile)
            else
                  call fvl_logerror("Invalid data operation")
            end if
      else
            call fvl_logerror("Invalid output data format")
      end if
else if(fvl_trim(datainformat)=="fvd3") then
      if(fvl_trim(dataoutformat)=="fvd3") then
            if(fvl_trim(operation)=="magfield") then
                  call fvl_loginfo("Computing magfield from "//fvl_trim(datainpath)//" to "//fvl_trim(dataoutpath)//"...")
                  call fvl_datamagfield3d(datainfile,dataoutfile)
            else if(fvl_trim(operation)=="sqrfield") then
                  call fvl_loginfo("Computing sqrfield from "//fvl_trim(datainpath)//" to "//fvl_trim(dataoutpath)//"...")
                  call fvl_datasqrfield3d(datainfile,dataoutfile)
            ! else if(fvl_trim(operation)=="perpfield") then
            !       call fvl_loginfo("Computing perpfield from "//fvl_trim(datainpath)//" to "//fvl_trim(dataoutpath)//"...")
            !       call fvl_dataperpfield3d(datainfile,dataoutfile)
            else
                  call fvl_logerror("Invalid data operation")
            end if
      else
            call fvl_logerror("Invalid output data format")
      end if
else
      call fvl_logerror("Invalid input data format")
end if

! ============================================================================
! STOP LIBRARY
! ============================================================================

call fvl_stop()

end program fvl_dataoperate_main
! end of file
