# Darwin - specific options

# The generic container templates (HashMapTemplate.h, AssociativeArrayTemplate.h)
# build module names by token pasting:
#     #define CONCAT(A,B) IDENTITY(A)IDENTITY(B)
#     #define MODULE_NAME DEFAULT_MODULE_NAME(HASH_TYPE)
# which relies on the preprocessor emitting adjacent macro expansions with
# nothing between them. Apple's /usr/bin/cpp (clang) inserts a space, so
# AttributeHashMap.F90 preprocesses to "module AttributeHashMap _mod" instead
# of "module AttributeHashMap_mod". The dependency scan then records the wrong
# module name and the build fails to order AttributeHashMap.o before its users.
# Prefer a real GNU cpp when one is installed (Homebrew gcc ships cpp-NN).
GNU_CPP := $(firstword $(foreach c,cpp-15 cpp-14 cpp-13 cpp-12,$(shell command -v $(c) 2>/dev/null)))
ifneq ($(GNU_CPP),)
  CPP = $(GNU_CPP) -P -traditional
else
  CPP = cpp -P -traditional
endif

# /usr/bin/cpp on macOS is traditional-by-default, so base.mk's
# "filter-out -traditional" cannot produce a normal-mode preprocessor for C.
# Name one explicitly; used only for scanning .c sources, where traditional
# mode cannot parse the current SDK headers.
CPP_C = clang -E -P

# The Fortran and C compilers here can target different architectures: a
# Homebrew gfortran under /usr/local is an x86_64 (Rosetta) build while
# /usr/bin/gcc is native arm64. Objects of different architectures cannot be
# combined into one static archive -
#   ranlib: archive member cputype (16777228) does not match previous archive
#           members cputype (16777223) (all members must match)
# - and -m64 does not switch a native arm64 clang over. Pin the C compiler to
# whatever the Fortran compiler targets. CFLAGS_MACHINE rather than CFLAGS:
# CFLAGS is a name the environment may already define (module systems set it),
# and it must not leak into the build on platforms that never set it here.
ifeq ($(COMPILER),gfortran)
  FC_ARCH := $(shell gfortran -dumpmachine 2>/dev/null | sed 's/-.*//;s/aarch64/arm64/')
  ifneq ($(FC_ARCH),)
    CFLAGS_MACHINE += -arch $(FC_ARCH)
  endif
endif

CPPFLAGS = -DMACHINE_MAC

# this is a hack to work around Xcode/MacPorts bug with Xcode 11.x.x
# (the compiler can't find the system include files)
### For now commenting it out - use explicit flags in ~/.modelErc
# XCODE_VERSION = $(word 2,$(shell xcodebuild -version))
#XCODE_VERSION_MAJOR = $(word 1,$(subst ., ,$(XCODE_VERSION)))
# ifeq ($(XCODE_VERSION_MAJOR),11)
#   CPATH_HACK=CPATH=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include 
# endif

# Homebrew's gcc is built against the SDK current at its release; when the
# installed Command Line Tools are newer, its linker cannot find libSystem
# and every link fails with "ld: library 'System' not found" - even for a
# hello-world program. Point the toolchain at the SDK actually installed.
# (Equivalent to exporting SDKROOT in the shell before building.)
SDKROOT ?= $(shell xcrun --show-sdk-path 2>/dev/null)
export SDKROOT
