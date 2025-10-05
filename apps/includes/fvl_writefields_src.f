! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: dictionaries
! Modification: April, 2024

! ----------------------------------------------------------------------------
! write models
! ----------------------------------------------------------------------------

! create directory
call fvl_createtimedirectory(solutionfiledirectory,mesh%getcurrenttime(),currenttimedirectory)

! write models
do solutionmodeli=1,size(solutionmodels,1)
      solutionmodel=>solutionmodels(solutionmodeli)%getptr()
      if(associated(solutionmodel)) then
            call solutionmodel%write()
      end if
end do

! write time report
call fvl_writewtimereport(fvl_trim(currenttimedirectory)//"/time.fvd","ascii")
call fvl_writewtimereport(fvl_trim(solutionfiledirectory)//"/time.fvd","ascii")

! write memory report
call fvl_writevmhwmreport(fvl_trim(currenttimedirectory)//"/memory.fvd","ascii")
call fvl_writevmhwmreport(fvl_trim(solutionfiledirectory)//"/memory.fvd","ascii")

! end of file
