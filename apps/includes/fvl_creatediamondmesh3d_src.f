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
mesh=fvl_mesh3d(dictfile1=controldictfile,dictlabel1="mesh",dictfile2=meshesdictfile,dictlabel2="mesh")

! write mesh report
call mesh%report("- mesh")

! state message
call fvl_loginfo(">>> Making diamond mesh")

! create diamond mesh
diamondmesh=fvl_diamond_mesh3d(mesh=mesh,dictfile=meshesdictfile,dictlabel="diamond_mesh")

! write mesh report
call diamondmesh%report("- diamond mesh")

! runtimes
call runtimes(1)%setptr(mesh)
call runtimes(2)%setptr(diamondmesh)

! domain patches
facespatch=fvl_patch3d(mesh=mesh,type=fvl_patch3d_facespatch,codes=fvl_integer_emptyarray)
cellspatch=fvl_patch3d(mesh=mesh,type=fvl_patch3d_cellspatch,codes=fvl_integer_emptyarray)
diamond_facespatch=fvl_patch3d(mesh=diamondmesh,type=fvl_patch3d_facespatch,codes=fvl_integer_emptyarray)
diamond_cellspatch=fvl_patch3d(mesh=diamondmesh,type=fvl_patch3d_cellspatch,codes=fvl_integer_emptyarray)

! boundary patches
boundarypatch=fvl_patch3d(mesh=mesh,type=fvl_patch3d_boundfacespatch,codes=fvl_integer_emptyarray)
diamond_boundarypatch=fvl_patch3d(mesh=diamondmesh,type=fvl_patch3d_boundfacespatch,codes=fvl_integer_emptyarray)

! end of file
