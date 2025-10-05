! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: scalar allocator
! Modification: March, 2025

! ============================================================================
! SET ALLOCATORS
! ============================================================================

! initialise field selector
__fvl_print_model__()_fieldselector=fvl_scalarfield3d_selector()

! initialise boundcond selector
__fvl_print_model__()_boundcondselector=fvl_scalarboundcond3d_selector()

! set field allocators
call __fvl_print_model__()_fieldselector%setallocator("constant",fvl_const_scalarfield3d_allocate)
call __fvl_print_model__()_fieldselector%setallocator("function",fvl_func_scalarfield3d_allocate)
call __fvl_print_model__()_fieldselector%setallocator("calculated",fvl_calc_scalarfield3d_allocate)

! set boundcond allocators
call __fvl_print_model__()_boundcondselector%setallocator("constant_value",fvl_constvalue_scalarboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("constant_deriv",fvl_constderiv_scalarboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("constant_mixed",fvl_constmixed_scalarboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("constant_combined",fvl_constcombined_scalarboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("function_value",fvl_constvalue_scalarboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("function_deriv",fvl_constderiv_scalarboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("function_mixed",fvl_constmixed_scalarboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("function_combined",fvl_constcombined_scalarboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("calculated_value",fvl_calcvalue_scalarboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("calculated_deriv",fvl_calcderiv_scalarboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("calculated_mixed",fvl_calcmixed_scalarboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("calculated_combined",fvl_calccombined_scalarboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("symmetry",fvl_symmetry_scalarboundcond3d_allocate)

! end of file
