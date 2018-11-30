cmake_minimum_required(VERSION 3.10)
project(megahit)

set(CMAKE_CXX_STANDARD 11)

include_directories(src)
include_directories(src/sparsepp)
FILE(GLOB_RECURSE ASMBL_SOURCE "src/assembly/*.cpp")
FILE(GLOB_RECURSE LCASM_SOURCE "src/local_asm/*.cpp")
FILE(GLOB_RECURSE IDBA_SOURCE "src/lib_idba/*.cpp")
FILE(GLOB_RECURSE SDBG_SOURCE "src/sdbg/*.cpp")
FILE(GLOB_RECURSE CX1_SOURCE "src/cx1/*.cpp")
FILE(GLOB_RECURSE TOOLKIT_SOURCE "src/tools/*.cpp")

LIST(APPEND OTHER_SOURCE
        src/main_asm_core_entry.cpp
        src/main_assemble.cpp
        src/build_read_lib.cpp
        src/main_iterate.cpp
        src/main_local_assemble.cpp
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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__XFILE__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

set(CMAKE_CXX_FLAGS_DEBUG "-g -O1")

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}: ${CMAKE_CXX_FLAGS}")

add_executable(megahit_asm_core ${OTHER_SOURCE} ${ASMBL_SOURCE} ${IDBA_SOURCE} ${SDBG_SOURCE} ${LCASM_SOURCE})
target_link_libraries(megahit_asm_core ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_executable(megahit_sdbg_build ${SDBG_SOURCE} ${CX1_SOURCE}
        src/main_sdbg_build.cpp src/main_sdbg_build.cpp src/kthread.cpp src/options_description.cpp
        src/sequence_manager.cpp)
target_link_libraries(megahit_sdbg_build ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_executable(megahit_toolkit ${TOOLKIT_SOURCE})
target_link_libraries(megahit_toolkit ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_custom_target(
        simple_test
        COMMAND cp ${CMAKE_SOURCE_DIR}/src/megahit .
        COMMAND cp -r ${CMAKE_SOURCE_DIR}/test_data .
        COMMAND rm -fr test_data/megahit_out
        COMMAND ./megahit --12 test_data/readsInterleaved1.fa.gz,test_data/readsInterleaved2.fa.bz2,test_data/readsInterleaved3.fa -o test_data/megahit_out -t 2
)
add_dependencies(simple_test megahit_asm_core megahit_toolkit megahit_sdbg_build)