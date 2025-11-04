cmake_minimum_required(VERSION 3.1)
project(ppq-pk)

#set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Using the latest C++
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

#Add all .cpp files and .h files here
# Add pkpd_mfq.cpp pkpd_mfq.h when needed
set(SRC_FILES main.cpp globals.cpp globals.h general_functions.cpp general_functions.h pkpd_simulation_functions.cpp pkpd_simulation_functions.h pkpd_artemether.cpp pkpd_artemether.h pkpd_dha.cpp pkpd_dha.h pkpd_artesunate.cpp pkpd_artesunate.h pkpd_lum.cpp pkpd_lum.h pkpd_ppq.cpp pkpd_ppq.h pkpd_adq.cpp pkpd_adq.h)


#set(SRC_FILES main.cpp pkpd_adq.cpp pkpd_adq.h pkpd_dha.cpp pkpd_dha.h pkpd_lum.cpp pkpd_lum.h pkpd_ppq.cpp pkpd_ppq.h)
add_executable(run_ppq_pk ${SRC_FILES})
target_link_libraries(run_ppq_pk PRIVATE GSL::gsl GSL::gslcblas)


install(TARGETS run_ppq_pk DESTINATION ${PROJECT_SOURCE_DIR}/bin)

