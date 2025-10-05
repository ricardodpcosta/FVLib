! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: scalar model
! Modification: February, 2025

! ============================================================================
! CREATE MODEL
! ============================================================================

type(fvl_scalarmodel3d),target::__fvl_print_model__()_models(__fvl_num_steps__)
type(fvl_scalarmodel3d),pointer::__fvl_print_model__()_model=>__fvl_print_model__()_models(1)
type(fvl_scalarfield3d_selector)::__fvl_print_model__()_fieldselector
type(fvl_scalarboundcond3d_selector)::__fvl_print_model__()_boundcondselector

! end of file
