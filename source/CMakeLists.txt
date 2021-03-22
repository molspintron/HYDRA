cmake_minimum_required (VERSION 3.8)

project("HYDRA")

# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

find_package(MPI REQUIRED)
include_directories(_PATH)
find_package(OpenMP)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)





#message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

add_executable(HYDRA main.cpp InputOutput.cpp)

target_link_libraries(HYDRA ${MPI_LIBRARIES})

target_link_libraries (HYDRA Eigen3::Eigen)

if(MPI_COMPILE_FLAGS)
    set_target_properties(HYDRA PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(HYDRA PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

#target_link_directories(HYDRA PUBLIC MPI:MPI_CXX)

if(OpenMP_CXX_FOUND)
    target_link_libraries(HYDRA OpenMP::OpenMP_CXX)
endif()



#[[
#The following block is for automatically generated documentation and requires installation of
#Doxygen, Sphinx, Breath and sphinx_rtd_theme


# Doxygen

# look for Doxygen package
find_package(Doxygen)

if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs_doxygen/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.out)
    set(DOXYGEN_OUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs_doxygen/xml)
    set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUT_DIR}/index.xml)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # Note: do not put "ALL" - this builds docs together with application EVERY TIME!
    add_custom_command( OUTPUT ${DOXYGEN_INDEX_FILE}
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            DEPENDS main.cpp InputOutput.cpp
            MAIN_DEPENDENCY ${DOXYGEN_OUT} ${DOXYGEN_IN}
            VERBATIM )
    add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})
else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

#SPHINX
find_package(Sphinx REQUIRED)

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/docs_sphinx)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs_sphinx)
set(SPHINX_INDEX ${SPHINX_BUILD}/index.html)
set(SPHINX_CONF ${SPHINX_SOURCE}/conf.py)

add_custom_command(OUTPUT ${SPHINX_INDEX}
        COMMAND ${SPHINX_EXECUTABLE} -b html
        -Dbreathe_projects.HYDRA=${DOXYGEN_OUT_DIR}
        ${SPHINX_SOURCE} ${SPHINX_BUILD}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS ${SPHINX_SOURCE}/index.rst ${DOXYGEN_INDEX_FILE}
        MAIN_DEPENDENCY ${SPHINX_CONF}
        MAIN_DEPENDENCY ${SPHINX_CONF}
        COMMENT "Generating documentation with Sphinx")
add_custom_target(Sphinx ALL DEPENDS ${SPHINX_INDEX})
]]