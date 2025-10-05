! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: dictionaries
! Modification: April, 2024

! ============================================================================
! VARIABLES
! ============================================================================

integer(kind=__fvl_integer_kind__)::solutionmodeli
character(len=__fvl_character_len__)::currenttimedirectory
class(fvl_model),pointer::solutionmodel
type(fvl_ptr_model),dimension(:),allocatable::solutionmodels

! end of file
