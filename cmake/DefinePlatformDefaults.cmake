# Set system vars

if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(LINUX TRUE)
endif(CMAKE_SYSTEM_NAME MATCHES "Linux")

if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
	set (OSX TRUE)
endif (CMAKE_SYSTEM_NAME MATCHES "Darwin")

if (CMAKE_SYSTEM_NAME MATCHES "FreeBSD")
    set(FREEBSD TRUE)
    set(BSD TRUE)
endif (CMAKE_SYSTEM_NAME MATCHES "FreeBSD")

if (CMAKE_SYSTEM_NAME MATCHES "OpenBSD")
    set(OPENBSD TRUE)
    set(BSD TRUE)
endif (CMAKE_SYSTEM_NAME MATCHES "OpenBSD")

if (CMAKE_SYSTEM_NAME MATCHES "NetBSD")
    set(NETBSD TRUE)
    set(BSD TRUE)
endif (CMAKE_SYSTEM_NAME MATCHES "NetBSD")

if (CMAKE_SYSTEM_NAME MATCHES "(Solaris|SunOS)")
    set(SOLARIS TRUE)
endif (CMAKE_SYSTEM_NAME MATCHES "(Solaris|SunOS)")

if (CMAKE_SYSTEM_NAME MATCHES "OS2")
    set(OS2 TRUE)
endif (CMAKE_SYSTEM_NAME MATCHES "OS2")

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    # require at least ifort 14.0
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.0)
        message("\n ---> IFORT COMPILER VERSION: " ${CMAKE_CXX_COMPILER_VERSION})
        message(FATAL_ERROR "---> IFORT version must be at least 14.0!")
    endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    # require at least gcc 4.9.1
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9.1)
        message("\n  ---> GNU COMPILER VERSION: " ${CMAKE_CXX_COMPILER_VERSION})
        message(FATAL_ERROR "---> GNU version must be at least 4.9.1!")
    endif()
else()
    message( FATAL_ERROR "Unrecognized compiler. Please use ifort or gfortran" )
endif()

