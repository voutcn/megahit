cmake_minimum_required(VERSION 3.10)
project(megahit)

set(CMAKE_CXX_STANDARD 11)

include_directories(src)
include_directories(src/sparsepp)
FILE(GLOB_RECURSE ASSEMBLY_SOURCE "src/assembly/*.cpp")
FILE(GLOB_RECURSE IDBA_SOURCE "src/lib_idba/*.cpp")
FILE(GLOB_RECURSE SDBG_SOURCE "src/sdbg/*.cpp")
LIST(APPEND OTHER_SOURCE
        src/asm_core.cpp
        src/assembler.cpp
        src/build_read_lib.cpp
        src/iterate_edges.cpp
        src/local_assemble.cpp
        src/local_assembler.cpp
        src/options_description.cpp
        src/sequence_manager.cpp
        )

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXXH_INLINE_ALL -D__STDC_FORMAT_MACROS")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function -funroll-loops -fprefetch-loop-arrays")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcnt -mbmi2")

set(CMAKE_CXX_FLAGS_DEBUG "-g -O1")

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}: ${CMAKE_CXX_FLAGS}")
set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(megahit_asm_core ${OTHER_SOURCE} ${ASSEMBLY_SOURCE} ${IDBA_SOURCE} ${SDBG_SOURCE})
target_link_libraries(megahit_asm_core ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})