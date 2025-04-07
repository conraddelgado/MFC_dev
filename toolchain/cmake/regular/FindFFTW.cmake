# Attempt to find FFTW
# URL:  https://www.fftw.org/
# DOCS: https://cmake.org/cmake/help/latest/command/find_library.html
#       https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html

include(FindPackageHandleStandardArgs)


find_path(FFTW_INCLUDE_DIR
    NAMES         fftw3.h
    PATH_SUFFIXES fftw fftw3
)

find_library(FFTW_LIBRARY
    NAMES         fftw3
    PATH_SUFFIXES fftw fftw3
    NAMES_PER_DIR
)


FIND_PACKAGE_HANDLE_STANDARD_ARGS(
    FFTW
    REQUIRED_VARS
        FFTW_LIBRARY
        FFTW_INCLUDE_DIR
)


if (FFTW_FOUND AND NOT TARGET FFTW::FFTW)
    set(FFTW_LIBRARIES    "${FFTW_LIBRARY}")
    set(FFTW_INCLUDE_DIRS "${FFTW_INCLUDE_DIR}")

    add_library(FFTW::FFTW UNKNOWN IMPORTED)

    set_target_properties(FFTW::FFTW PROPERTIES
        IMPORTED_LOCATION             "${FFTW_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
    )
endif()

# # find FFTW MPI
# find_path(FFTW_MPI_INCLUDE_DIR
#     NAMES         fftw3.h
#     PATH_SUFFIXES fftw fftw3
# )

# find_library(FFTW_MPI_LIBRARY
#     NAMES         fftw3_mpi
#     PATH_SUFFIXES fftw fftw3
#     NAMES_PER_DIR
# )

# FIND_PACKAGE_HANDLE_STANDARD_ARGS(
#     FFTW_MPI
#     REQUIRED_VARS
#         FFTW_MPI_LIBRARY
#         FFTW_MPI_INCLUDE_DIR
# )

# # create the imported target for FFTW MPI if found
# if (FFTW_MPI_FOUND AND NOT TARGET FFTW_MPI::FFTW_MPI)
#     set(FFTW_MPI_LIBRARIES    "${FFTW_MPI_LIBRARY}")
#     set(FFTW_MPI_INCLUDE_DIRS "${FFTW_MPI_INCLUDE_DIR}")

#     add_library(FFTW_MPI::FFTW_MPI UNKNOWN IMPORTED)

#     set_target_properties(FFTW_MPI::FFTW_MPI PROPERTIES
#         IMPORTED_LOCATION             "${FFTW_MPI_LIBRARIES}"
#         INTERFACE_INCLUDE_DIRECTORIES "${FFTW_MPI_INCLUDE_DIRS}"
#     )
# endif()
