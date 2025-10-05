! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: check and analyse dictionaries
! Modification: March, 2025

#include "macros.f"

program fvl_dictreport_main

use fvl_lib

implicit none

! ============================================================================
! VARIABLES
! ============================================================================

character(len=__fvl_character_len__)::dictpath,dictform,dictformat
type(fvl_dict_file)::dictfile

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
      call fvl_loginfo("About: Utility tool for dict inspection and reporting.")
      call fvl_loginfo("Usage: fvl-dictreport <dict_file> [options]")
      call fvl_loginfo("Arguments:")
      call fvl_loginfo("      <dict_file>             Input dictionary file in .fvd format.")
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
dictpath=adjustr(fvl_getarg(1))

! ============================================================================
! RUN UTILITY TOOL
! ============================================================================

! check input dict form
if(dictpath(__fvl_character_len__-3:__fvl_character_len__)==".bin") then
      ! dictform="binary"
      ! dictpath=dictpath(__fvl_character_len__-3:__fvl_character_len__)
      ! dictpath=adjustr(dictpath)
      call fvl_logerror("Invalid input dict format")
else
      dictform="ascii"
end if

! check input dict format
if(dictpath(__fvl_character_len__-3:__fvl_character_len__)==".fvd") then
      dictformat="fvd"
else
      call fvl_logerror("Invalid input dict format")
end if

! initialize files
dictpath=adjustl(dictpath)
dictfile=fvl_dict_file(path=dictpath,form=dictform,action="read")

! report dictionary
call dictfile%open()
call dictfile%readdicts()
call dictfile%report()
call dictfile%close()

! ============================================================================
! STOP LIBRARY
! ============================================================================

call fvl_stop()

end program fvl_dictreport_main
! end of file
