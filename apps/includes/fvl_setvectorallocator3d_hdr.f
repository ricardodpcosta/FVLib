! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: vector allocator
! Modification: March, 2025

! ============================================================================
! SET ALLOCATORS
! ============================================================================

! initialise field selector
__fvl_print_model__()_fieldselector=fvl_vectorfield3d_selector()

! initialise boundcond selector
__fvl_print_model__()_boundcondselector=fvl_vectorboundcond3d_selector()

! set field allocators
call __fvl_print_model__()_fieldselector%setallocator("constant",fvl_const_vectorfield3d_allocate)
call __fvl_print_model__()_fieldselector%setallocator("function",fvl_func_vectorfield3d_allocate)
call __fvl_print_model__()_fieldselector%setallocator("calculated",fvl_calc_vectorfield3d_allocate)

! set boundcond allocators
call __fvl_print_model__()_boundcondselector%setallocator("constant_value",fvl_constvalue_vectorboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("constant_deriv",fvl_constderiv_vectorboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("constant_mixed",fvl_constmixed_vectorboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("constant_dirmixed",fvl_constdirmixed_vectorboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("constant_combined",fvl_constcombined_vectorboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("function_value",fvl_constvalue_vectorboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("function_deriv",fvl_constderiv_vectorboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("function_mixed",fvl_constmixed_vectorboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("function_dirmixed",fvl_constdirmixed_vectorboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("function_combined",fvl_constcombined_vectorboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("calculated_value",fvl_calcvalue_vectorboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("calculated_deriv",fvl_calcderiv_vectorboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("calculated_mixed",fvl_calcmixed_vectorboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("calculated_dirmixed",fvl_calcdirmixed_vectorboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("calculated_combined",fvl_calccombined_vectorboundcond3d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("symmetry",fvl_symmetry_vectorboundcond3d_allocate)

! end of file
