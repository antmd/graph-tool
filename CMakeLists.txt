cmake_minimum_required( VERSION 2.8 )

project(graph-tool)

set( CMAKE_INSTALL_PREFIX /usr/local )
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")

include(CheckCXXCompilerFlag)

check_cxx_compiler_flag("-std=c++11" SUPPORTS_CXX11)
if( SUPPORTS_CXX11 )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    check_cxx_compiler_flag("-std=c++0x" SUPPORTS_CXX0X)
    if( SUPPORTS_CXX0X )
        set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    else()
        message("C++11 is not supported")
    endif()
endif()

find_package(Boost COMPONENTS graph math_c99 python QUIET REQUIRED)
find_package(PythonLibs)
find_package(NumPy)

set(GRAPH-TOOL_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
set(GRAPH-TOOL_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src )

set( CMAKE_INCLUDE_CURRENT_DIR ON )

# Globbing is not recommended -- CMake cannot tell when new files are added
# But for a generic skeleton...
file( GLOB_RECURSE 
    GRAPH-TOOL_INCLUDE_FILES 
    ${GRAPH-TOOL_SRC_DIR}/*.h
    ${GRAPH-TOOL_SRC_DIR}/graph/*.hh
    ${GRAPH-TOOL_SRC_DIR}/boost-workaround/*.hpp
    )

file( GLOB_RECURSE 
    GRAPH-TOOL_SOURCE_FILES
    ${GRAPH-TOOL_SRC_DIR}/boost-workaround/*.cpp
    ${GRAPH-TOOL_SRC_DIR}/graph/*.cc
    )

file(
    COPY ${GRAPH-TOOL_INCLUDE_FILES} 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include/GRAPH-TOOL
    )

include_directories( ${Boost_INCLUDE_DIR} )
include_directories(${GRAPH-TOOL_SRC_DIR}/graph)
include_directories(${GRAPH-TOOL_SRC_DIR}/boost-workaround)
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${NUMPY_INCLUDES})

set( ALL_SOURCES ${GRAPH-TOOL_INCLUDE_FILES} ${GRAPH-TOOL_SOURCE_FILES} )

add_compile_options(-Wall -Wextra -ftemplate-backtrace-limit=0 -std=gnu11 -std=c++1y -stdlib=libc++ -I/usr/local/include -DNDEBUG -std=gnu++11 -ftemplate-depth-250 -Wno-deprecated -Wno-unknown-pragmas -O3 -fvisibility=default -fvisibility-inlines-hidden -Wno-unknown-pragmas -I/usr/local/Cellar/google-sparsehash/2.0.2/include  -DSPARSEHASH_PREFIX=sparsehash)

add_executable( GRAPH-TOOL ${ALL_SOURCES} )

target_link_libraries( GRAPH-TOOL 
    ${PYTHON_LIBRARIES} 
    ${Boost_LIBRARIES} 
    )

# INSTALL
# -------

install (
    TARGETS GRAPH-TOOL
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    )

