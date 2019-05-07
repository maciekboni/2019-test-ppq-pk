cmake_minimum_required(VERSION 3.1)
project(ppq-pk)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

if (APPLE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -O0 -fno-omit-frame-pointer")
endif ()

if (WIN32)
    string(REGEX REPLACE "/MD" "/MT" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
    string(REGEX REPLACE "/MD" "/MT" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
endif ()
#find_package(Boost REQUIRED)
#include_directories(${Boost_INCLUDE_DIRS})

find_package(GSL REQUIRED)


set(SRC_FILES main.cpp pkdynamics_ppq.cpp pkdynamics_ppq.h pkpd_dha.cpp pkpd_dha.h)

add_executable(run_ppq_pk ${SRC_FILES})
target_link_libraries(run_ppq_pk PRIVATE GSL::gsl GSL::gslcblas)


install(TARGETS run_ppq_pk DESTINATION ${PROJECT_SOURCE_DIR}/bin)

