#!/usr/local/bin/python3
# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#|  _ _ _   _   _   _      |                                   |
#| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
#| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
#| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
#| |_|       \_/   |_|_|_| |  Release: January, 2022           |
#|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
# About: generate a mesh
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
          print("[FVL] About: Utility tool to generate a mesh.")
          print("[FVL] Usage: fvl-gmsh-makemesh <model-file> <mesh-file> [options]")
          print("[FVL] Arguments:")
          print("[FVL]      <model-file>             Input model file in geo format.")
          print("[FVL]      <mesh-file>              Output mesh file in msh format.")
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
meshpath = sys.argv[2]

# open a file with model data to create a new model
gmsh.open(modelpath)

# generate a mesh of the current model
gmsh.model.mesh.generate(dim=3)

# launch the GUI to see the results
if "--popup" in sys.argv:
    gmsh.fltk.run()

# write mesh to a file
gmsh.option.setNumber("Mesh.MshFileVersion",2.0)
gmsh.write(fileName=meshpath)

# finalize the Gmsh API
gmsh.finalize()

# end of file
