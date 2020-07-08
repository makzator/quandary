# Set location of PETSC
#PETSC_DIR=/usr/workspace/wsa/wave/petsc-3.13
#PETSC_ARCH=arch-linux-c-debug

# Set location of SLEPC
SLEPC_DIR=${HOME}/Software/slepc-3.13.3

# Set Braid location, or comment out
#BRAID_DIR = ${HOME}/Numerics/xbraid_solveadjointwithxbraid

# Choose to run sanity tests
SANITY_CHECK = false

#######################################################
# Typically no need to change anything below
#######################################################

# Add optional Slepc
ifdef SLEPC_DIR
CXX_OPT = -DWITH_SLEPC
LDFLAGS_OPT = -L${SLEPC_DIR}/lib -L${SLEPC_DIR}/${PETSC_ARCH}/lib -lslepc 
INCDIR_OPT = -I${SLEPC_DIR}/${PETSC_ARCH}/include -I${SLEPC_DIR}/include
endif

# Add optional Braid include and library location
ifdef BRAID_DIR
BRAID_INC_DIR = $(BRAID_DIR)/braid
BRAID_LIB_FILE = $(BRAID_DIR)/braid/libbraid.a
CXX_OPT += -DWITH_BRAID
LDFLAGS_OPT += ${BRAID_LIB_FILE}
INCDIR_OPT += -I$(BRAID_INC_DIR) 
endif

# Add optional sanity check
ifeq ($(SANITY_CHECK), true)
CXX_OPT += -DSANITY_CHECK
endif

# Set direction of source and header files, and build direction
SRC_DIR   = src
INC_DIR   = include
BUILD_DIR = build

# list all source and object files
SRC_FILES  = $(wildcard $(SRC_DIR)/*.cpp)
SRC_FILES += $(wildcard $(SRC_DIR)/*/*.cpp)
OBJ_FILES  = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC_FILES))

# set include directory
INC = -I$(INC_DIR) -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include ${INCDIR_OPT} 

# Set Library paths and flags
LDPATH  = ${PETSC_DIR}/${PETSC_ARCH}/lib
LDFLAGS = -lm -lblas -llapack -L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc ${LDFLAGS_OPT} 

# Set compiler and flags 
CXX=mpicxx
CXXFLAGS= -O3 -std=c++11 -lstdc++ $(CXX_OPT)


# Rule for linking main
main: $(OBJ_FILES)
	$(CXX) -o $@ $(OBJ_FILES) $(LDFLAGS) -L$(LDPATH)

# Rule for building all src files
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(INC) 
	@$(CXX) -MM $< -MP -MT $@ -MF $(@:.o=.d) $(INC) 


.PHONY: all cleanup

# use 'make cleanup' to remove object files and executable
cleanup:
	rm -fr $(BUILD_DIR) 
	rm -f  main 


# include the dependency files
-include $(OBJ_FILES:.o=.d)
