VASPML_BUILD_DIR = build
VASPML_SRC_DIR = src
VASPML_RES_DIR = res

# Include compiler/preprocessor/linker information.
include makefile.include

# Common variables (e.g. for colored output)
include $(VASPML_RES_DIR)/makefile.common

TARGETS = libvaspml applications benchmarks
CLEAN_TARGETS = $(patsubst %, clean-%, $(TARGETS))

.PHONY: all \
        pre \
        $(TARGETS) \
        $(CLEAN_TARGETS) \
        clean

all: $(TARGETS)

ifneq (,$(findstring no-logo,$(MAKE_OPTIONS)))
logo: ;
else
logo:
	@printf "$(C_VPURPLE) _    _____   _____ ____  __  _____  \n"
	@printf "$(C_VBLUE)| |  / /   | / ___// __ \/  |/  / /  \n"
	@printf "$(C_VCYAN)| | / / /| | \__ \/ /_/ / /|_/ / /   \n"
	@printf "$(C_VGREEN)| |/ / ___ |___/ / ____/ /  / / /___ \n"
	@printf "$(C_VRED)|___/_/  |_/____/_/   /_/  /_/_____/ $(C_OFF)\n\n"
endif

pre: logo

libvaspml: pre
	@printf "Building VASPml C++ library in directory $(VASPML_BUILD_DIR)/$@...\n"
	if [ ! -d $(VASPML_BUILD_DIR)/$@ ] ; then mkdir -p $(VASPML_BUILD_DIR)/$@ ; fi
	cp $(VASPML_SRC_DIR)/$@/makefile $(VASPML_BUILD_DIR)/$@
	$(MAKE) -C $(VASPML_BUILD_DIR)/$@ depend
	$(MAKE) -C $(VASPML_BUILD_DIR)/$@ build
	@printf -- "$(LOG_SEPARATOR)\n"

applications: libvaspml
	@printf "Building VASPml applications in directory $(VASPML_BUILD_DIR)/$@...\n"
	if [ ! -d $(VASPML_BUILD_DIR)/$@ ] ; then mkdir -p $(VASPML_BUILD_DIR)/$@ ; fi
	cp $(VASPML_SRC_DIR)/$@/makefile $(VASPML_BUILD_DIR)/$@
	$(MAKE) -C $(VASPML_BUILD_DIR)/$@ build
	@printf -- "$(LOG_SEPARATOR)\n"

benchmarks: libvaspml applications
	@printf "Building VASPml benchmarks in directory $(VASPML_BUILD_DIR)/$@...\n"
	if [ ! -d $(VASPML_BUILD_DIR)/$@ ] ; then mkdir -p $(VASPML_BUILD_DIR)/$@ ; fi
	cp $(VASPML_SRC_DIR)/$@/makefile $(VASPML_BUILD_DIR)/$@
	$(MAKE) -C $(VASPML_BUILD_DIR)/$@ build
	@printf -- "$(LOG_SEPARATOR)\n"

code-format-check:
	find $(patsubst %, src/%, $(TARGETS)) -regex '.*\.\(cpp\|hpp\|c\|h\)' -exec $(VASPML_RES_DIR)/run-clang-format.py {} \;

code-format:
	find $(patsubst %, src/%, $(TARGETS)) -regex '.*\.\(cpp\|hpp\|c\|h\)' -exec clang-format -style=file -i {} \;

clean: $(CLEAN_TARGETS)

$(CLEAN_TARGETS):
	@printf "Removing build directory for target:   \"$(subst clean-,,$@)\".\n"
	$(RM) -r $(VASPML_BUILD_DIR)/$(subst clean-,,$@)
	$(RM) bin/*
	$(RM) lib/*
	$(RM) include/*
	@printf -- "$(LOG_SEPARATOR)\n"

$(VERBOSE).SILENT:
