cmake_minimum_required(VERSION 3.27)
project(greedy_cross_cpp)

set(CMAKE_CXX_STANDARD 17)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")
FIND_PACKAGE(Boost REQUIRED COMPONENTS system)
FIND_PACKAGE(Armadillo REQUIRED COMPONENTS system)
FIND_PACKAGE(Eigen3 REQUIRED COMPONENTS system)

if(Boost_FOUND AND Armadillo_FOUND AND Eigen3_FOUND)
    include_directories(${Boost_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS} ${Eigen3_INCLUDE_DIRS})
    add_executable(greedy_cross_cpp
            ../include/maxvol.h
            ../src/main.cpp
            ../utility/vector_concatenation.h
            ../include/greedy_cross.h
            ../utility/info.h
            ../utility/timing.h
            ../utility/change_vector_type.h
            ../tests/test_functions.h
    )
    target_link_libraries(greedy_cross_cpp ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES} Eigen3::Eigen)
endif()
