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
! write reports
! ----------------------------------------------------------------------------

! write time report
call fvl_writewtimereport(fvl_trim(solutionfiledirectory)//"/time.fvd","ascii")

! write memory report
call fvl_writevmhwmreport(fvl_trim(solutionfiledirectory)//"/memory.fvd","ascii")

! end of file
