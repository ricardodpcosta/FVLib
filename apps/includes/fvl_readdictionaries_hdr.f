! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: dictionaries
! Modification: February, 2025

! ============================================================================
! VARIABLES
! ============================================================================

character(len=__fvl_character_len__)::meshesdictpath,modelsdictpath,schemesdictpath,solutiondictpath,controldictpath,&
      postdictpath,solutionfiledirectory,solutionfileform
integer(kind=__fvl_integer_kind__)::solutiontimeprecision
type(fvl_dict_file)::meshesdictfile,modelsdictfile,schemesdictfile,solutiondictfile,controldictfile,postdictfile

! end of file
