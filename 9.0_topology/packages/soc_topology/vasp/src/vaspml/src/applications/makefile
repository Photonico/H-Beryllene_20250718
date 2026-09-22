.PHONY: \
	all \
	build \
	clean

# Names of all applications.
APP = mlab2json \
      vaspml-predict \
      vaspml-predict_nvidia \
      test_cublas \
      test_linalg \
      test_linalg2 \
      test_mpiSHMEM \
      test_MultiType \
      test_predict_mpi \
      test_scalapack

# Location of source files.
BASE_DIR = ../..
SRC_DIR = $(BASE_DIR)/src/applications
RES_DIR = $(BASE_DIR)/res
LIB_DIR = $(BASE_DIR)/lib
INC_DIR = $(BASE_DIR)/include
BIN_DIR = $(BASE_DIR)/bin

# Include compiler/preprocessor/linker information.
include $(BASE_DIR)/makefile.include

# Common variables (e.g. for colored output)
include $(RES_DIR)/makefile.common

# Derive source files from application names.
SRC = $(addsuffix .cpp,$(APP))

# Static pattern rule to fetch all changed source files.
$(SRC): %: $(SRC_DIR)/%
	@printf $(LOG_INOUT) "$(C_VGREY)Collecting source $< " $(LOG_ARROW) "$(C_VPURPLE)$(@)$(C_OFF)\n"
	rsync -u $< .

LIB_PREREQ = $(LIB_DIR)/libvaspml.a
LIB_LINK = -lvaspml
fetch: headers

# Fetch all headers from lib directories (may be required for coverage analysis).
headers:
	@printf $(LOG_INOUT) "$(C_VGREY)Collecting headers from libvaspml " $(LOG_ARROW) "$(C_VPURPLE)*.hpp$(C_OFF)\n"
	rsync -u ../libvaspml/*.hpp .

# Fetch all headers (Fortran) from lib directories (may be required for coverage analysis).
headers_f: headers
	@printf $(LOG_INOUT) "$(C_VGREY)Collecting headers from libvaspml_f " $(LOG_ARROW) "$(C_VPURPLE)*.hpp$(C_OFF)\n"
	rsync -u ../libvaspml_f/*.hpp .

# Default compile rule.
$(APP): %: %.cpp $(LIB_PREREQ)
	@printf $(LOG_INOUT) "$(C_VGREY)Compiling $< " $(LOG_ARROW) "$(C_VCYAN)$(@)$(C_OFF)\n"
	$(CXX) $(CXXFLAGS) $(INCLUDE) -I$(INC_DIR) $(CPP_OPTIONS) -c $<
	$(LD) $(CXXFLAGS) $(INCLUDE) -I$(INC_DIR) $(CPP_OPTIONS) -o $@ $@.o -L$(LIB_DIR) $(LIB_LINK) $(LDFLAGS)
	@printf $(LOG_INOUT) "$(C_VGREY)Copying executable $(@)" $(LOG_ARROW) "$(C_VRED)$(BIN_DIR)/$(@)$(C_OFF)\n";
	cp $@ $(BIN_DIR)

# Build all applications.
build: fetch $(APP)

# Remove all sources, objects, temporary and resulting files.
ifneq ($(realpath $(SRC_DIR)),$(CURDIR))
clean:
	$(RM) *.cpp *.o $(APP)
else
clean:
	echo "ERROR: Invoking clean target in source folder, this is dangerous, aborting..."
	exit 1
endif

# Silence compiler output (re-enable with "VERBOSE=1" make argument).
$(VERBOSE).SILENT:
