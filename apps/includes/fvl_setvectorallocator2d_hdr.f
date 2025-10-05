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
__fvl_print_model__()_fieldselector=fvl_vectorfield2d_selector()

! initialise boundcond selector
__fvl_print_model__()_boundcondselector=fvl_vectorboundcond2d_selector()

! set field allocators
call __fvl_print_model__()_fieldselector%setallocator("constant",fvl_const_vectorfield2d_allocate)
call __fvl_print_model__()_fieldselector%setallocator("function",fvl_func_vectorfield2d_allocate)
call __fvl_print_model__()_fieldselector%setallocator("calculated",fvl_calc_vectorfield2d_allocate)

! set boundcond allocators
call __fvl_print_model__()_boundcondselector%setallocator("constant_value",fvl_constvalue_vectorboundcond2d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("constant_deriv",fvl_constderiv_vectorboundcond2d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("constant_mixed",fvl_constmixed_vectorboundcond2d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("constant_dirmixed",fvl_constdirmixed_vectorboundcond2d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("constant_combined",fvl_constcombined_vectorboundcond2d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("function_value",fvl_constvalue_vectorboundcond2d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("function_deriv",fvl_constderiv_vectorboundcond2d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("function_mixed",fvl_constmixed_vectorboundcond2d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("function_dirmixed",fvl_constdirmixed_vectorboundcond2d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("function_combined",fvl_constcombined_vectorboundcond2d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("calculated_value",fvl_calcvalue_vectorboundcond2d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("calculated_deriv",fvl_calcderiv_vectorboundcond2d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("calculated_mixed",fvl_calcmixed_vectorboundcond2d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("calculated_dirmixed",fvl_calcdirmixed_vectorboundcond2d_allocate)
call __fvl_print_model__()_boundcondselector%setallocator("calculated_combined",fvl_calccombined_vectorboundcond2d_allocate)
!!!!call __fvl_print_model__()_boundcondselector%setallocator("symmetry",fvl_symmetry_vectorboundcond2d_allocate)

! end of file
