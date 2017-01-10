cmake_minimum_required(VERSION 2.8)
PROJECT(ParticleFilter)

SET(CMAKE_BUILD_TYPE "RelWithDebInfo")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")

find_package(VTK REQUIRED)


include(${VTK_USE_FILE})

add_subdirectory(dataLoader)
set(LIBS LaserDataLoader)


add_subdirectory(gui)
set(LIBS ${LIBS} Gui)

add_subdirectory(laserSimulator)
set(LIBS ${LIBS} lasersimulator)

include_directories ("${PROJECT_SOURCE_DIR}")

add_executable(ParticleFilter pf_main.cpp)

target_link_libraries(ParticleFilter ${LIBS} ${VTK_LIBRARIES})
