# --------
# Compiler
# --------
CXX             = mpic++

# --------------
# Compiler flags
# --------------
CXXFLAGS       := -O3 -fopenmp
CXXFLAGS       += -g
CXXFLAGS       += -std=c++17 -pedantic-errors
CXXFLAGS       += -Wall -Wextra
# Extra flags required for testing and coverage.
CXXFLAGS       += --coverage -fno-default-inline -fno-inline -fno-inline-small-functions -fno-elide-constructors

# -------------------
# Include directories
# -------------------
INCLUDE        := -I$(OPENBLAS_ROOT)/include

# -------
# Linking
# -------
LD              = $(CXX)
LDFLAGS        := -L$(OPENBLAS_ROOT)/lib -lopenblas

# ------------
# Preprocessor
# ------------
CPP_OPTIONS    += -DVASPML_USE_CBLAS
#CPP_OPTIONS    += -DVASPML_DEBUG_LEVEL=3

# --------------------
# Dependency generator
# --------------------
CPP_DEP         = $(CXX)

# --------------
# Archiving tool
# --------------
AR              = ar
ARFLAGS         = rcs

# -----------------------
# Additional make options
# -----------------------
#MAKE_OPTIONS   += --no-color
#MAKE_OPTIONS   += --no-logo

# -----------------------------------
# Boost library
# (optional, only needed for testing)
# -----------------------------------
BOOST_ROOT     ?=
BOOST_INCLUDE   = -isystem${BOOST_ROOT}/include
BOOST_LIB       = -L${BOOST_ROOT}/lib -lboost_unit_test_framework
