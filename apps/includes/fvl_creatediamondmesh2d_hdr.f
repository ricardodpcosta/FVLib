! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: create meshes
! Modification: February, 2025

! ============================================================================
! VARIABLES
! ============================================================================

type(fvl_ptr_runtime)::runtimes(2)
type(fvl_mesh2d)::mesh,diamondmesh
type(fvl_patch2d)::cellspatch,diamond_cellspatch,edgespatch,diamond_edgespatch,boundarypatch,diamond_boundarypatch

! end of file
