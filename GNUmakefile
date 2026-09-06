# Make targets for building, installing and cleaning the th egg.
#
# Requires: CHICKEN 5, make, cc/gcc, ar and an OpenBLAS runtime library
# (see README). Building delegates to chicken-install, so the .egg file
# stays the single source of build truth.
#
#   make            build in place (th.so, import lib, static object, ...)
#   make install    install into a CHICKEN repository
#   make check      run the test suite against an installation
#   make clean      remove all generated files
#   make info       show where things will go

EGG             := th
DEFAULT_PREFIX  := /usr/local/lib/chicken
DEFAULT_REPO    := $(DEFAULT_PREFIX)/11

# chicken-install installs into /usr/local by default. If that is not
# writable, fall back to a per-user repository. CHICKEN_INSTALL_PREFIX /
# CHICKEN_INSTALL_REPOSITORY given on the command line take precedence.
ROOT_WRITABLE   := $(shell test -w /usr/local && echo 1)
ifeq ($(ROOT_WRITABLE),1)
INSTALL_ENV     :=
REPO            := $(or $(CHICKEN_INSTALL_REPOSITORY),/usr/local/lib/chicken/11)
else
INSTALL_ENV     := CHICKEN_INSTALL_PREFIX="$(or $(CHICKEN_INSTALL_PREFIX),$(DEFAULT_PREFIX))" \
                   CHICKEN_INSTALL_REPOSITORY="$(or $(CHICKEN_INSTALL_REPOSITORY),$(DEFAULT_REPO))"
REPO            := $(or $(CHICKEN_INSTALL_REPOSITORY),$(DEFAULT_REPO))
endif

.PHONY: all build install check test clean info

all: build

build:
	chicken-install -s -no-install

install:
	$(INSTALL_ENV) chicken-install -s

check test: install
	cd $(HOME) && \
	CHICKEN_REPOSITORY=$(REPO) CHICKEN_INCLUDE_PATH=$(REPO) \
	csi -q -s $(CURDIR)/tests/run.scm

clean:
	rm -f th.so th.import.so th.o th.static.o th.link th.import.scm \
	      th.build.sh th.install.sh
	rm -rf TH/build

info:
	@echo "chicken   : $(shell chicken-install -version 2>/dev/null || csc -version 2>/dev/null | head -1)"
	@echo "repo      : $(REPO)"
	@echo "install   : $(shell echo '$(INSTALL_ENV)' | sed -n 's/.*CHICKEN_INSTALL_PREFIX="\([^"]*\)".*/\1/p')"
	@echo "root mode : $(ROOT_WRITABLE)"
