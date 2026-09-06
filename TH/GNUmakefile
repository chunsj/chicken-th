# GNUmakefile for TH -- builds libTHTensor as a shared library.
#
#   make                     build with default BLAS (openblas on Linux, accelerate on macOS)
#   make BLAS=mkl            explicit BLAS backend (openblas|mkl|accelerate|netlib|auto)
#   make BLAS=mkl OMP=1      enable OpenMP threading
#   make MARCH=znver1        override -march/-mtune (x86_64 only)
#   make EXTRA_CFLAGS="-g"   append extra compile flags
#   make CFLAGS="-O2 -g"     fully override compile flags
#   make info                print the resolved configuration
#   make clean && make       full rebuild after switching BLAS/OMP
#   make install PREFIX=...  install library and headers (default /usr/local)
#   make BUILD_DIR=...       output directory for objects and library (default build)

VERSION    := 0
OS         := $(shell uname -s)
ARCH       := $(shell uname -m)

# --- Platform (output filename and shared-library linker flags) ------------
ifeq ($(OS),Darwin)
  DLEXT  := dylib
  SHARED := -dynamiclib -install_name libTHTensor.$(VERSION).$(DLEXT)
else
  DLEXT  := so
  SHARED := -shared -Wl,-soname,libTHTensor.$(VERSION).$(DLEXT)
endif

LIBNAME    := libTHTensor.$(VERSION).$(DLEXT)
SOFTLINK   := libTHTensor.$(DLEXT)

# --- Compiler ---------------------------------------------------------------
CC ?= gcc

OPT      ?= -O3
MARCH    ?= native
MTUNE    ?= $(MARCH)

ifeq ($(ARCH),x86_64)
  SIMD_FLAGS ?= -msse4.2 -mfma -mavx -mavx2 -march=$(MARCH) -mtune=$(MTUNE)
endif

ifeq ($(OS),Darwin)
  PIC  :=
else
  PIC  := -fPIC
endif

EXTRA_CFLAGS ?=
CFLAGS ?= $(OPT) $(SIMD_FLAGS) $(OMP_CFLAGS) $(PIC) $(EXTRA_CFLAGS)
CPPFLAGS ?=
DEPFLAGS := -MMD -MP

# --- OpenMP ------------------------------------------------------------------
OMP ?= 0
CC_IS_CLANG := $(shell $(CC) --version 2>/dev/null | grep -q clang && echo 1)

ifeq ($(OMP),1)
  ifeq ($(OS),Darwin)
    ifeq ($(CC_IS_CLANG),1)
      OMP_CFLAGS := -Xpreprocessor -fopenmp
      OMP_LIBS   := -lomp
    else
      OMP_CFLAGS := -fopenmp
      OMP_LIBS   := -fopenmp
    endif
  else
    OMP_CFLAGS := -fopenmp
    OMP_LIBS   := -fopenmp
  endif
else
  OMP_CFLAGS :=
  OMP_LIBS   :=
endif

# --- BLAS backend ------------------------------------------------------------
# Default per platform: openblas on Linux, accelerate on macOS.
ifeq ($(OS),Darwin)
  DEFAULT_BLAS := accelerate
else
  DEFAULT_BLAS := openblas
endif

OPENBLAS_SONAME := $(shell ldconfig -p 2>/dev/null | grep -m1 'libopenblas\.so\.0' | awk '{print $$1}')

BLAS ?= $(DEFAULT_BLAS)

ifeq ($(BLAS),auto)
  ifeq ($(OS),Darwin)
    override BLAS := accelerate
  else
    override BLAS := $(shell { (ldconfig -p 2>/dev/null | grep -q 'libopenblas\.so' && echo openblas) || \
                                (ldconfig -p 2>/dev/null | grep -q libmkl_intel_lp64 && echo mkl) || \
                                echo netlib; })
  endif
  BLAS_AUTO := 1
endif

ifeq ($(BLAS),mkl)
  BLAS_LIBS := -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
else ifeq ($(BLAS),openblas)
  ifneq ($(shell ldconfig -p 2>/dev/null | grep -q 'libopenblas\.so ' && echo 1),)
    BLAS_LIBS := -lopenblas -lpthread
  else ifneq ($(OPENBLAS_SONAME),)
    BLAS_LIBS := -l:$(OPENBLAS_SONAME) -lpthread
  else
    BLAS_LIBS := -lopenblas -lpthread
  endif
else ifeq ($(BLAS),accelerate)
  BLAS_LIBS := -framework Accelerate
else ifeq ($(BLAS),netlib)
  BLAS_LIBS := -lcblas -llapack -lblas -lpthread
else
  $(error unsupported BLAS '$(BLAS)' -- use openblas, mkl, accelerate, netlib or auto)
endif

ifeq ($(OS),Darwin)
  LDLIBS ?=
else
  LDLIBS ?= -lm -lpthread
endif
LDFLAGS ?=

# cblas.h is required for THGeneral.c's `#include <cblas.h>`
ifneq ($(shell test -f /usr/include/x86_64-linux-gnu/cblas.h && echo 1),)
  CPPFLAGS += -I/usr/include/x86_64-linux-gnu
else ifneq ($(shell test -f /usr/include/openblas/cblas.h && echo 1),)
  CPPFLAGS += -I/usr/include/openblas
else ifneq ($(shell test -f /usr/local/include/cblas.h && echo 1),)
  CPPFLAGS += -I/usr/local/include
endif

# --- THConfig.h selection ----------------------------------------------------
# Shipped variants: linux (openblas default, mkl), macos (accelerate default,
# openblas). The active THConfig.h is left untouched unless it matches one of the
# shipped variants (yours may be customized).
ifeq ($(BLAS),mkl)
  THCONFIG ?= THConfig.h.linux.mkl
else ifeq ($(BLAS),openblas)
  ifeq ($(OS),Darwin)
    THCONFIG ?= THConfig.h.macos.openblas
  else
    THCONFIG ?= THConfig.h.linux.openblas
  endif
else ifeq ($(BLAS),accelerate)
  THCONFIG ?= THConfig.h.macos.accelerate
else
  THCONFIG ?= THConfig.h.macos.accelerate
endif

KNOWN_CONFIGS := THConfig.h.linux.mkl THConfig.h.linux.openblas \
                 THConfig.h.macos.accelerate THConfig.h.macos.openblas

# --- Sources ------------------------------------------------------------------
BUILD_DIR ?= build

SRCS      := $(wildcard *.c)
ifeq ($(ARCH),x86_64)
  SIMD_SRCS := vector/AVX.c vector/AVX2.c
else ifneq ($(filter aarch64 armv7l armv8l,$(ARCH)),)
  SIMD_SRCS := vector/NEON.c
else ifneq ($(filter ppc64el ppc64le ppc64,$(ARCH)),)
  SIMD_SRCS := vector/VSX.c
else
  SIMD_SRCS :=
endif

ROOT_OBJ := $(addprefix $(BUILD_DIR)/,$(SRCS:.c=.o))
SIMD_OBJ := $(addprefix $(BUILD_DIR)/,$(notdir $(SIMD_SRCS:.c=.o)))
OBJS     := $(ROOT_OBJ) $(SIMD_OBJ)
DEPS     := $(OBJS:.o=.d)

# --- Rules --------------------------------------------------------------------
.PHONY: all clean config info help install uninstall .thconfig static

all: $(BUILD_DIR)/$(LIBNAME)

$(BUILD_DIR):
	mkdir -p $@

.thconfig:
	@if ! cmp -s THConfig.h $(THCONFIG); then \
	  managed=1; \
	  for v in $(KNOWN_CONFIGS); do cmp -s THConfig.h $$v && managed=0 && break; done; \
	  if [ $$managed = 0 ]; then \
	    cp $(THCONFIG) THConfig.h; \
	    echo "  [config] THConfig.h  <- $(THCONFIG)  (BLAS=$(BLAS))"; \
	  else \
	    echo "  [config] THConfig.h customized -- leaving unchanged"; \
	    echo "           (BLAS=$(BLAS) expects THConfig.h to match $(THCONFIG))"; \
	  fi; \
	fi
	@touch $@

$(ROOT_OBJ): $(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	$(CC) $(DEPFLAGS) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

ifneq ($(SIMD_OBJ),)
$(SIMD_OBJ): $(BUILD_DIR)/%.o: vector/%.c | $(BUILD_DIR)
	$(CC) $(DEPFLAGS) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
endif

$(BUILD_DIR)/$(LIBNAME): .thconfig $(OBJS) | $(BUILD_DIR)
	$(CC) $(SHARED) -o $@ $(OBJS) $(BLAS_LIBS) $(OMP_LIBS) $(LDLIBS) $(LDFLAGS)
	@ln -sf $(LIBNAME) $(BUILD_DIR)/$(SOFTLINK)
	@echo "  built $(BUILD_DIR)/$(LIBNAME)  (CC=$(CC) BLAS=$(BLAS) OMP=$(OMP))"

static: $(BUILD_DIR)/libTHTensor.a

$(BUILD_DIR)/libTHTensor.a: .thconfig $(OBJS) | $(BUILD_DIR)
	rm -f $@
	ar rcs $@ $(OBJS)
	@echo "  built $(BUILD_DIR)/libTHTensor.a  ($(words $(OBJS)) objects)"

config: .thconfig

clean:
	rm -rf $(BUILD_DIR)

# --- Install ------------------------------------------------------------------
PREFIX      ?= /usr/local
INCDIR      := $(DESTDIR)$(PREFIX)/include/TH
LIBDIR      := $(DESTDIR)$(PREFIX)/lib
HDRS        := $(wildcard *.h)

install: $(BUILD_DIR)/$(LIBNAME)
	install -d $(LIBDIR) $(INCDIR)/generic $(INCDIR)/vector
	install -m 644 $(HDRS) $(INCDIR)
	install -m 644 generic/*.h $(INCDIR)/generic/
	install -m 644 vector/*.h $(INCDIR)/vector/
	install -m 755 $(BUILD_DIR)/$(LIBNAME) $(LIBDIR)
	ln -sf $(LIBNAME) $(LIBDIR)/$(SOFTLINK)

uninstall:
	rm -f $(LIBDIR)/$(LIBNAME) $(LIBDIR)/$(SOFTLINK)
	rm -rf $(INCDIR)

# --- Info ----------------------------------------------------------------------
info:
	@echo "OS        : $(OS)"
	@echo "ARCH      : $(ARCH)"
	@echo "CC        : $(CC)"
	@echo "OPT       : $(OPT)"
	@echo "SIMD      : $(SIMD_FLAGS)"
	@echo "DEFAULT   : $(DEFAULT_BLAS)"
	@echo "BLAS      : $(BLAS)$(if $(BLAS_AUTO), (auto-detected))"
	@echo "BLAS_LIBS : $(BLAS_LIBS)"
	@echo "THCONFIG  : $(THCONFIG)"
	@echo "OMP       : $(OMP)"
	@echo "CFLAGS    : $(CFLAGS)"
	@echo "LIBS      : $(BLAS_LIBS) $(OMP_LIBS) $(LDLIBS)"
	@echo "OUTPUT    : $(BUILD_DIR)/$(LIBNAME) -> $(BUILD_DIR)/$(SOFTLINK)"

help:
	@echo "Targets: all, config, info, clean, install, uninstall"
	@echo
	@echo "Variables (all settable via make VAR=...):"
	@echo "  BLAS       openblas|mkl|accelerate|netlib|auto   (default: openblas Linux, accelerate macOS)"
	@echo "  OMP        0|1                                   enable OpenMP"
	@echo "  MARCH      -march value, e.g. native, znver1     (x86_64 only)"
	@echo "  CC, OPT, CFLAGS, CPPFLAGS, LDFLAGS, EXTRA_CFLAGS"
	@echo "  BUILD_DIR  object and library output dir (default build)"
	@echo "  PREFIX     install prefix (default /usr/local)"

-include $(DEPS)