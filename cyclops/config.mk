### INSTALLATION TARGET DIRECTORY (for make install)
INSTALL_DIR = /usr/local

### LINK TIME LIBRARIES AND FLAGS
#libraries and flags for link time (irrelevant if only building CTF lib and not examples/tests)
LIB_PATH     = 
LIB_FILES    =  -llapack   -lblas -lctf
LINKFLAGS    = 
LD_LIB_PATH  = 
SO_LIB_PATH  = 
SO_LIB_FILES =   -lblas  -llapack -lctf
LDFLAGS      = 


### COMPILE TIME INCLUDES AND FLAGS
#C++ compiler 
CXX         = mpicxx
#includes for compile time
INCLUDES    = 
#optimization flags, some intel compiler versions may run into errors when using -fast or -ipo
CXXFLAGS    = -O3 -std=c++0x -DOMP_OFF -Wall 
#command to make library out of object files
AR          = ar

#macros to be defined throughout the code, use the below in combination with appropriate external libraries
#Include in DEFS -DUSE_LAPACK to build with LAPACK functionality, 
#                -DUSE_SCALAPACK to build with SCALAPACK functionality
#                -DUSE_BATCH_GEMM to build {without, with} batched BLAS routines
#                -DUSE_MKL to build with MKL sparse matrix kernels
#                -DUSE_HPTT to build with optimized tensor transposition routines from HPTT library
DEFS        = -D_POSIX_C_SOURCE=200112L -D__STDC_LIMIT_MACROS -D_DARWIN_C_SOURCE -DFTN_UNDERSCORE=1 -DUSE_LAPACK 


### Optional: PROFILING AND TUNING
#uncomment below to enable performance profiling
#DEFS       += -DPROFILE -DPMPI
#also uncomment below to enable performance profiling without epochs (started when first CTF World created and stopped when destroyed)
#DEFS       += -DAUTO_PROFILE
#uncomment below to enable automatic performance tuning (loses reproducibility of results)
#recommended usage is to run model_trainer with -DTUNE at scale for suitable duration to obtain suitable architecutral parameters, 
#then recompile with parameters without -DTUNE
#Note: -DTUNE requires lapack (include -mkl or -llapack in LIBS) and the inclusion of above performance profiling flags
#DEFS       += -DTUNE

### Optional: DEBUGGING AND VERBOSITY
#uncomment below to enable CTF execution output (1 for basic contraction information on start-up and contractions)
#DEFS       += -DVERBOSE=1
#uncomment to set debug level to dump information about mapping and internal CTF actions and activate asserts
#DEFS       += -DDEBUG=1

### FULL COMPILE COMMAND AND LIBRARIES
#used to compile all plain C++ files
FCXX        = $(CXX) $(CXXFLAGS) $(DEFS) $(INCLUDES)
#link-line for all executables
LIBS        = $(LIB_PATH) $(LIB_FILES) $(LINKFLAGS)
#compiler for CUDA files (used to compile CUDA code only when -DOFFLOAD and -DUSE_CUDA are in DEFS, otherwise should be same as FCXX with -x c++)
OFFLOAD_CXX = $(CXX) -x c++ $(CXXFLAGS) $(DEFS) $(INCLUDES)
