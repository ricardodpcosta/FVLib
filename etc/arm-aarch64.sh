# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#|  _ _ _   _   _   _      |                                   |
#| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
#| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
#| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
#| |_|       \_/   |_|_|_| |  Release: June 2021               |
#|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
# About: configuration file for fvl-make
# Modification: January, 2025
# compiler command
FVL_CMP="armflang";
# create archive command
FVL_AR="armllvm-ar";
# --plugin=/opt/FJSVxos/devkit/aarch64/libexec/gcc/aarch64-linux-gnu/8/liblto_plugin.so";
# update archive command
FVL_RANLIB="armllvm-ranlib";
# --plugin=/opt/FJSVxos/devkit/aarch64/libexec/gcc/aarch64-linux-gnu/8/liblto_plugin.so";
# openmp compatibility
FVL_OMP_VERSION="5";
# openmp simd directives
FVL_OMP_SIMD="0";
# gnu instrinsics compatibility
FVL_GNU_INTRINSICS="0";
# gnu directives compatibility
FVL_GNU_DIRECTIVES="0";
# gnu instrinsics compatibility
FVL_INTEL_INTRINSICS="0";
# gnu directives compatibility
FVL_INTEL_DIRECTIVES="0";
# fuji instrinsics compatibility
FVL_FUJI_INTRINSICS="0";
# fuji directives compatibility
FVL_FUJI_DIRECTIVES="0";
# arm instrinsics compatibility
FVL_ARM_INTRINSICS="1";
# arm directives compatibility
FVL_ARM_DIRECTIVES="1";
# production mode
FVL_LANG_FLAGS_OPT="
      -std=f2008
      -ffree-form
      -ffixed-line-length-none
      ";
FVL_PREPROC_FLAGS_OPT="
      -cpp
      -r8
      -fopenmp
      ";
FVL_WARN_FLAGS_OPT="
      -Wall
      -Wno-unused-command-line-argument
      ";
FVL_OPTIM_FLAGS_OPT="
      -march=armv8.2-a+sve
      -mtune=native
      -mcpu=native
      -O3
      -fvectorize
      -fno-signed-zeros
      -finline-functions
      -flto
      -armpl
      -fsimdmath
      ";
FVL_DEBUG_FLAGS_OPT="
      -g0
      -fno-trapping-math
      ";
FVL_DIR_FLAGS_OPT="
      -I./
      -I$FVL_INC_PATH_OPT
      -I$FVL_SRC_PATH/includes
      -I$FVL_APPS_PATH/includes
      -J$FVL_INC_PATH_OPT
      ";
FVL_LIB_FLAGS_OPT="
      -L$FVL_LIB_PATH_OPT
      -lfvl
      -lm
      ";
FVL_DEF_FLAGS_OPT="
      -Dfvlmpi=1
      -Dfvlcharacterlen=250
      -Dfvlshortkind=2
      -Dfvlintegerkind=4
      -Dfvlfloatkind=4
      -Dfvlrealkind=8
      -Dfvlquadkind=16
      -Dfvlcomplexkind=8
      -Dfvllogicalkind=4
      -Dfvlfloatprecision=1.0d-25
      -Dfvlrealprecision=1.0d-250
      -Dfvlquadprecision=1.0d-250
      -Dfvlalignbytes=8
      -Dfvldebug=0
      -Dfvlverbose=2
      -Dfvlmaxnumprocesses=100
      -Dfvlmaxnumthreads=100
      -Dfvlstacksize=1024000
      -Dfvlompversion=$FVL_OMP_VERSION
      -Dfvlompsimd=$FVL_OMP_SIMD
      -Dfvlgnuintrinsics=$FVL_GNU_INTRINSICS
      -Dfvlgnudirectives=$FVL_GNU_DIRECTIVES
      -Dfvlintelintrinsics=$FVL_INTEL_INTRINSICS
      -Dfvlinteldirectives=$FVL_INTEL_DIRECTIVES
      -Dfvlfujiintrinsics=$FVL_FUJI_INTRINSICS
      -Dfvlfujidirectives=$FVL_FUJI_DIRECTIVES
      ";
# debuggging mode
FVL_LANG_FLAGS_DBG="
      -std=f2008
      -ffree-form
      -ffixed-line-length-none
      ";
FVL_PREPROC_FLAGS_DBG="
      -cpp
      -r8
      -fopenmp
      ";
FVL_WARN_FLAGS_DBG="
      -Wall
      -Wno-unused-command-line-argument
      ";
FVL_OPTIM_FLAGS_DBG="
      -march=armv8.2-a+sve
      -mtune=native
      -mcpu=native
      -O0
      -fno-vectorize
      -fsigned-zeros
      -fno-inline-functions
      -fno-lto
      -fno-simdmath
      ";
FVL_DEBUG_FLAGS_DBG="
      -g
      -ftrapping-math
      ";
FVL_DIR_FLAGS_DBG="
      -I./
      -I$FVL_INC_PATH_DBG
      -I$FVL_SRC_PATH/includes
      -I$FVL_APPS_PATH/includes
      -J$FVL_INC_PATH_DBG
      ";
FVL_LIB_FLAGS_DBG="
      -L$FVL_LIB_PATH_DBG
      -lfvl
      -lm
      ";
FVL_DEF_FLAGS_DBG="
      -Dfvlmpi=1
      -Dfvlcharacterlen=250
      -Dfvlshortkind=2
      -Dfvlintegerkind=4
      -Dfvlfloatkind=4
      -Dfvlrealkind=8
      -Dfvlquadkind=16
      -Dfvlcomplexkind=8
      -Dfvllogicalkind=4
      -Dfvlfloatprecision=1.0d-25
      -Dfvlrealprecision=1.0d-250
      -Dfvlquadprecision=1.0d-250
      -Dfvlalignbytes=8
      -Dfvldebug=0
      -Dfvlverbose=2
      -Dfvlmaxnumprocesses=100
      -Dfvlmaxnumthreads=100
      -Dfvlstacksize=1024000
      -Dfvlompversion=$FVL_OMP_VERSION
      -Dfvlompsimd=$FVL_OMP_SIMD
      -Dfvlgnuintrinsics=$FVL_GNU_INTRINSICS
      -Dfvlgnudirectives=$FVL_GNU_DIRECTIVES
      -Dfvlintelintrinsics=$FVL_INTEL_INTRINSICS
      -Dfvlinteldirectives=$FVL_INTEL_DIRECTIVES
      -Dfvlfujiintrinsics=$FVL_FUJI_INTRINSICS
      -Dfvlfujidirectives=$FVL_FUJI_DIRECTIVES
      ";
# profiling mode
FVL_LANG_FLAGS_PRF="
      -std=f2008
      -ffree-form
      -ffixed-line-length-none
      ";
FVL_PREPROC_FLAGS_PRF="
      -cpp
      -r8
      -fopenmp
      ";
FVL_WARN_FLAGS_PRF="
      -Wall
      -Wno-unused-command-line-argument
      ";
FVL_OPTIM_FLAGS_PRF="
      -march=armv8.2-a+sve
      -mtune=native
      -mcpu=native
      -O2
      -fvectorize
      -fno-signed-zeros
      -finline-functions
      -flto
      -armpl
      -fsimdmath
      ";
FVL_DEBUG_FLAGS_PRF="
      -g0
      -fno-trapping-math
      ";
FVL_DIR_FLAGS_PRF="
      -I./
      -I$FVL_INC_PATH_PRF
      -I$FVL_SRC_PATH/includes
      -I$FVL_APPS_PATH/includes
      -J$FVL_INC_PATH_PRF
      ";
FVL_LIB_FLAGS_PRF="
      -L$FVL_LIB_PATH_PRF
      -lfvl
      -lm
      ";
FVL_DEF_FLAGS_PRF="
      -Dfvlmpi=1
      -Dfvlcharacterlen=250
      -Dfvlshortkind=2
      -Dfvlintegerkind=4
      -Dfvlfloatkind=4
      -Dfvlrealkind=8
      -Dfvlquadkind=16
      -Dfvlcomplexkind=8
      -Dfvllogicalkind=4
      -Dfvlfloatprecision=1.0d-25
      -Dfvlrealprecision=1.0d-250
      -Dfvlquadprecision=1.0d-250
      -Dfvlalignbytes=8
      -Dfvldebug=0
      -Dfvlverbose=2
      -Dfvlmaxnumprocesses=100
      -Dfvlmaxnumthreads=100
      -Dfvlstacksize=1024000
      -Dfvlompversion=$FVL_OMP_VERSION
      -Dfvlompsimd=$FVL_OMP_SIMD
      -Dfvlgnuintrinsics=$FVL_GNU_INTRINSICS
      -Dfvlgnudirectives=$FVL_GNU_DIRECTIVES
      -Dfvlintelintrinsics=$FVL_INTEL_INTRINSICS
      -Dfvlinteldirectives=$FVL_INTEL_DIRECTIVES
      -Dfvlfujiintrinsics=$FVL_FUJI_INTRINSICS
      -Dfvlfujidirectives=$FVL_FUJI_DIRECTIVES
      ";
# end of file
