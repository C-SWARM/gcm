AC_PREREQ([2.63])
m4_define([LIBGCM_MAJOR], [0])
m4_define([LIBGCM_MINOR], [0])
m4_define([LIBGCM_PATCH], [0])
AC_INIT([gcm], LIBGCM_MAJOR.LIBGCM_MINOR.LIBGCM_PATCH, [gcm-bugs@nd.edu])

: ${CXXFLAGS="-O3 -g"}

AC_LANG([C++])

# Standard preliminary autotools stuff.
AC_CONFIG_SRCDIR([utils/include/math_help.h])
#AC_CONFIG_HEADERS([include/config.h include/hpx/types.h])
AC_CONFIG_AUX_DIR([m4])
AC_CONFIG_MACRO_DIR([m4])

AX_PROG_CXX_MPI([test "x" == x], [have_mpi=yes])
AX_CXX_COMPILE_STDCXX_14([ext], [mandatory])

AM_INIT_AUTOMAKE([1.11 nostdinc dist-bzip2 subdir-objects foreign tar-ustar -Wall
                  -Werror parallel-tests color-tests])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
LT_PREREQ([2.4.2])
LT_INIT([disable-shared])

AC_ARG_WITH(ttl,
 [AS_HELP_STRING([--with-ttl=<path>], [Specify path to TTL installation.])],
 [CPPFLAGS="$CPPFLAGS -I$withval/include"])

AC_CHECK_HEADER([ttl/ttl.h], [], [AC_MSG_ERROR(Failed to find ttl.h)])

AC_ARG_WITH(mkl,
 [AS_HELP_STRING([--with-mkl], [Path to Intel MKL @<:@default=$MKLROOT@:>@])],
 [], [with_mkl=$MKLROOT])

CPPFLAGS="$CPPFLAGS -m64 -I$with_mkl/include"
AC_CHECK_HEADERS([$cblas_h $lapack_h], [], [AC_MSG_ERROR(Could not find cblas headers)])

LDFLAGS="$LDFLAGS -L$with_mkl/lib/intel64 -Wl,-rpath,$with_mkl/lib/intel64"
LIBS="$LIBS -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
AC_SEARCH_LIBS([dgetri], [], [], [AC_MSG_ERROR(Could not link MKL at $with_mkl)])
AC_SEARCH_LIBS([dgetrf], [], [], [AC_MSG_ERROR(Could not link MKL at $with_mkl)])

AC_CONFIG_FILES(
 [Makefile
  J2_plasticity/Makefile
  constitutive_model_handle/Makefile
  crystal_plasticity/Makefile
  damage/Makefile
  data_structure/Makefile
  elasticity/Makefile
  lib/Makefile
  material/Makefile
  utils/Makefile
  J2_plasticity/src/Makefile
  constitutive_model_handle/src/Makefile
  crystal_plasticity/src/Makefile
  damage/src/Makefile
  data_structure/src/Makefile
  elasticity/src/Makefile
  material/src/Makefile
  utils/src/Makefile
  test/Makefile
  test/J2_plasticity/Makefile
  test/crystal_plasticity_grains/Makefile
  test/crystal_plasticity_single_crystal/Makefile
  test/crystal_plasticity_with_thermal/Makefile
  test/damage/Makefile
  test/damage_split/Makefile
  test/data_structure/Makefile
  test/elasticity/Makefile
  test/slip_system/Makefile])

AC_CONFIG_FILES([test/run.sh], [chmod +x test/run.sh])
AC_CONFIG_FILES([test/crystal_plasticity_with_thermal/test_crystal_plasticity_with_thermal.sh], [chmod +x test/crystal_plasticity_with_thermal/test_crystal_plasticity_with_thermal.sh])

AM_CONDITIONAL([HAVE_MPI], [test "x$have_mpi" == xyes])

AC_OUTPUT
