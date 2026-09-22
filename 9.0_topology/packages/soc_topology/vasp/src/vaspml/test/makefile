# Location of source files.
BASE_DIR = ..
RES_DIR = $(BASE_DIR)/res
BUILD_DIR = build
SRC_DIR = src
COV_DIR = $(BASE_DIR)/build

# Include compiler/preprocessor/linker information.
include $(BASE_DIR)/makefile.include

# Common variables (e.g. for colored output)
include $(RES_DIR)/makefile.common

TARGETS = libvaspml
CLEAN_TARGETS = $(patsubst %, clean-%, $(TARGETS))

.PHONY: all \
	$(TARGETS) \
	coverage \
	clean \
	$(CLEAN_TARGETS) \
	clean-coverage

all: $(TARGETS)

$(TARGETS): common
	@printf "Building and running test part: \"$(@)\" in directory $(BUILD_DIR)/$@...\n"
	if [ ! -d $(BUILD_DIR)/$@ ] ; then mkdir -p $(BUILD_DIR)/$@ ; fi
	cp $(SRC_DIR)/$@/makefile $(BUILD_DIR)/$@
	$(MAKE) -C $(BUILD_DIR)/$@ depend
	$(MAKE) -C $(BUILD_DIR)/$@ run
	@printf -- "$(LOG_SEPARATOR)\n"

common:
	@printf "Building common test infrastructure in directory $(BUILD_DIR)/$@...\n"
	if [ ! -d $(BUILD_DIR)/$@ ] ; then mkdir -p $(BUILD_DIR)/$@ ; fi
	cp $(SRC_DIR)/$@/makefile $(BUILD_DIR)/$@
	$(MAKE) -C $(BUILD_DIR)/$@ depend
	$(MAKE) -C $(BUILD_DIR)/$@ build
	@printf -- "$(LOG_SEPARATOR)\n"

coverage:
	@printf "$(C_OFF)Generating coverage report in directory cov-html...\n$(C_VGREY)"
	gcovr --xml-pretty --exclude-unreachable-branches -e $(COV_DIR)/benchmarks -e $(COV_DIR)/applications --print-summary -o cov-xml/coverage.xml --html-details cov-html/index.html -r $(COV_DIR)
	@printf -- "$(C_OFF)$(LOG_SEPARATOR)\n"

clean: $(CLEAN_TARGETS) clean-common clean-coverage

$(CLEAN_TARGETS):
	@printf "Removing build directory for test part:   \"$(subst clean-,,$@)\".\n"
	$(RM) -r $(BUILD_DIR)/$(subst clean-,,$@)
	@printf -- "$(LOG_SEPARATOR)\n"

clean-common:
	@printf "Removing build directory for common test infrastructure.\n"
	$(RM) -r $(BUILD_DIR)/$(subst clean-,,$@)
	$(RM) include/*
	$(RM) lib/*
	@printf -- "$(LOG_SEPARATOR)\n"

clean-coverage:
	@printf "Removing coverage data.\n"
	$(RM) cov-html/*.html cov-xml/*.xml
	@printf -- "$(LOG_SEPARATOR)\n"

# Silence compiler output (re-enable with "VERBOSE=1" make argument).
$(VERBOSE).SILENT:
