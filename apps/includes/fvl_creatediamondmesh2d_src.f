! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: create meshes
! Modification: February, 2025

! ----------------------------------------------------------------------------
! create mesh
! ----------------------------------------------------------------------------

! state message
call fvl_loginfo(">>> Reading mesh")

! read mesh
mesh=fvl_mesh2d(dictfile1=controldictfile,dictlabel1="mesh",dictfile2=meshesdictfile,dictlabel2="mesh")

! write mesh report
call mesh%report("- mesh")

! state message
call fvl_loginfo(">>> Making diamond mesh")

! create diamond mesh
diamondmesh=fvl_diamond_mesh2d(mesh=mesh,dictfile=meshesdictfile,dictlabel="diamond_mesh")

! write mesh report
call diamondmesh%report("- diamond mesh")

! runtimes
call runtimes(1)%setptr(mesh)
call runtimes(2)%setptr(diamondmesh)

! domain patches
edgespatch=fvl_patch2d(mesh=mesh,type=fvl_patch2d_edgespatch,codes=fvl_integer_emptyarray,label="fluid")
cellspatch=fvl_patch2d(mesh=mesh,type=fvl_patch2d_cellspatch,codes=fvl_integer_emptyarray,label="fluid")
diamond_edgespatch=fvl_patch2d(mesh=diamondmesh,type=fvl_patch2d_edgespatch,codes=fvl_integer_emptyarray,label="fluid")
diamond_cellspatch=fvl_patch2d(mesh=diamondmesh,type=fvl_patch2d_cellspatch,codes=fvl_integer_emptyarray,label="fluid")

! boundary patches
boundarypatch=fvl_patch2d(mesh=mesh,type=fvl_patch2d_boundedgespatch,codes=fvl_integer_emptyarray)
diamond_boundarypatch=fvl_patch2d(mesh=diamondmesh,type=fvl_patch2d_boundedgespatch,codes=fvl_integer_emptyarray)

! end of file
