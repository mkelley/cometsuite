AC_INIT([CometSuite],
        [0.9.5],
        [Michael S. Kelley <msk@astro.umd.edu>],
        [cometsuite])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_MACRO_DIR(m4)
AM_INIT_AUTOMAKE

AC_PROG_LIBTOOL
AC_PROG_CC
AC_PROG_CXX

AC_CHECK_LIB([m], [abs])
if test x"${ac_cv_lib_m_abs}" = xno; then
   AC_MSG_ERROR(libm not found.)
fi

AC_CHECK_LIB([CCfits], [exit], [], AC_MSG_ERROR([libCCfits not found.]))
if test x"${ac_cv_lib_CCfits_exit}" = xno; then
   AC_MSG_ERROR(libCCfits not found.)
fi

AC_CHECK_LIB([cfitsio], [exit], [], AC_MSG_ERROR([libcfitsio not found.]))
if test x"${ac_cv_lib_cftitsio_exit}" = xno; then
   AC_MSG_ERROR(libcfitsio not found.)
fi

AC_CHECK_LIB([getxyz], [getxyz])
if test x"${ac_cv_lib_getxyz_getxyz}" = xno; then
   AC_MSG_ERROR(libgetxyz not found.)
fi

AC_CHECK_HEADER([getxyz.h], [])
if test x"${ac_cv_header_getxyz_h}" = xno; then
   AC_MSG_ERROR([getxyz.h not found.])
fi

AM_MISSING_PROG(HELP2MAN, help2man, $missing_dir)

AC_CONFIG_FILES([
   Makefile
   doc/Makefile
   api/Doxyfile
   src/Makefile
   data/Makefile
   test/Makefile
])

AC_OUTPUT
