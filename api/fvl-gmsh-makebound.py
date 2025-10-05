#!/usr/local/bin/python3
# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#|  _ _ _   _   _   _      |                                   |
#| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
#| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
#| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
#| |_|       \_/   |_|_|_| |  Release: January, 2022           |
#|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
# About: generate a boundary discrete representation
# Modification: March, 2022

# import modules
import gmsh
import sys
import math

# initialize the Gmsh API
gmsh.initialize()

# check arguments
for opt in sys.argv:
    if(opt[0]=="-"):
        if(opt!="-h" and opt!="--help" and opt!="-p" and opt!="--popup"):
            print("[FVL] ERROR: Invalid option.")
            quit()

# check help option
for opt in sys.argv:
    if(opt=="-h" or opt=="--help"):
          print("[FVL] About: Utility tool to generate a boundary discrete representation.")
          print("[FVL] Usage: fvl-gmsh-makebound <model-file> <data-file> [options]")
          print("[FVL] Arguments:")
          print("[FVL]      <model-file>             Input model file in geo format.")
          print("[FVL]      <data-file>              Output data file in msh format.")
          print("[FVL] Options:")
          print("[FVL]      -h|--help                Displays help message.")
          print("[FVL]      -p|--popup               Lauch the Gmsh graphical interface.")
          quit()

# check arguments
if(len(sys.argv)<3):
    print("[FVL] ERROR: Missing arguments")
    quit()

# get arguments
modelpath = sys.argv[1]
datapath = sys.argv[2]

# open a file with model data to create a new model
gmsh.open(modelpath)

# generate a surface mesh of the current model
gmsh.model.mesh.generate(dim=2)

# retrieve the exact normals and the curvature at all the mesh nodes,
# i.e. not normals and curvatures computed from the mesh, but directly
# evaluated on the geometry)
normals = []
curvatures = []

# for each surface in the model
for entity in gmsh.model.getEntities(dim=2):
    # retrieve the surface tag
    surfacetag = entity[1]

    # get the mesh nodes on the surface, including those on the boundary
    # (contrary to internal nodes, which store their parametric coordinates,
    # boundary nodes will be reparametrized on the surface in order to compute
    # their parametric coordinates, the result being different when
    # reparametrized on another adjacent surface)
    tags,coord,param = gmsh.model.mesh.getNodes(dim=2,tag=surfacetag,\
        includeBoundary=False,returnParametricCoord=True)

    # get the surface normals on all the points on the surface corresponding to
    # the parametric coordinates of the nodes
    norm = gmsh.model.getNormal(tag=surfacetag,parametricCoord=param)

    # In the same way, get the curvature
    curv = gmsh.model.getCurvature(dim=2,tag=surfacetag,parametricCoord=param)

    # store the normals and the curvatures so that they can be written on disk
    # and displed as list-based post-processing views
    for i in range(0,len(coord),3):
        normals.append(coord[i])
        normals.append(coord[i+1])
        normals.append(coord[i+2])
        normals.append(norm[i])
        normals.append(norm[i+1])
        normals.append(norm[i+2])
        curvatures.append(coord[i])
        curvatures.append(coord[i+1])
        curvatures.append(coord[i+2])
        curvatures.append(curv[i//3])

# create a list-based vector view on points to display the normals, and a scalar
# view on points to display the curvatures
vn = gmsh.view.add(name="Normals")
gmsh.view.addListData(tag=vn,dataType="VP",numEle=len(normals)//6,data=normals)
vc = gmsh.view.add(name="Curvatures")
gmsh.view.addListData(tag=vc,dataType="SP",numEle=len(curvatures)//4,data=curvatures)

# launch the GUI to see the results
if "--popup" in sys.argv:
    gmsh.fltk.run()

# write normals and curvatures to a file
gmsh.option.setNumber("Mesh.MshFileVersion",2.0)
gmsh.view.write(tag=vn,fileName=datapath,append=False)
gmsh.view.write(tag=vc,fileName=datapath,append=True)

# finalize the Gmsh API
gmsh.finalize()

# end of file
