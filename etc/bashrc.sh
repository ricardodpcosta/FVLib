# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#|  _ _ _   _   _   _      |                                   |
#| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
#| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
#| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
#| |_|       \_/   |_|_|_| |  Release: June 2021               |
#|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
# About: configuration file for fvl-make
# Modification: March, 2022
# absolute path to sourced script
export FVL_SCRIPT_PATH=$(cd $(dirname "${BASH_SOURCE}");pwd);
# absolute path to installation
export FVL_PREFIX=$(dirname "$FVL_SCRIPT_PATH");
# installation name
export FVL_NAME=$(basename "$FVL_PREFIX");
# release version
export FVL_VERSION="1.0";
# remove environment variables from path
export PATH=${PATH#$FVL_ETC_PATH:};
export PATH=${PATH#$FVL_API_PATH:};
export PATH=${PATH#$FVL_BIN_PATH_OPT:};
export PATH=${PATH#$FVL_BIN_PATH_DBG:};
export PATH=${PATH#$FVL_BIN_PATH_PRF:};
# general paths
export FVL_ETC_PATH=$FVL_PREFIX/etc;
export FVL_API_PATH=$FVL_PREFIX/api;
export FVL_SRC_PATH=$FVL_PREFIX/src;
export FVL_APPS_PATH=$FVL_PREFIX/apps;
export FVL_TOOLS_PATH=$FVL_PREFIX/apps/tools;
export FVL_SOLVERS_PATH=$FVL_PREFIX/apps/solvers;
export FVL_POSTS_PATH=$FVL_PREFIX/apps/posts;
export FVL_BUILD_PATH=$FVL_PREFIX/build;
# production mode paths
export FVL_INC_PATH_OPT=$FVL_BUILD_PATH/inc_opt;
export FVL_LIB_PATH_OPT=$FVL_BUILD_PATH/lib_opt;
export FVL_BIN_PATH_OPT=$FVL_BUILD_PATH/bin_opt;
# debugging mode paths
export FVL_INC_PATH_DBG=$FVL_BUILD_PATH/inc_dbg;
export FVL_LIB_PATH_DBG=$FVL_BUILD_PATH/lib_dbg;
export FVL_BIN_PATH_DBG=$FVL_BUILD_PATH/bin_dbg;
# profiling mode paths
export FVL_INC_PATH_PRF=$FVL_BUILD_PATH/inc_prf;
export FVL_LIB_PATH_PRF=$FVL_BUILD_PATH/lib_prf;
export FVL_BIN_PATH_PRF=$FVL_BUILD_PATH/bin_prf;
# export environment variables to path
export PATH=$FVL_BIN_PATH_PRF:$PATH;
export PATH=$FVL_BIN_PATH_DBG:$PATH;
export PATH=$FVL_BIN_PATH_OPT:$PATH;
export PATH=$FVL_API_PATH:$PATH;
export PATH=$FVL_ETC_PATH:$PATH;
# end of file
