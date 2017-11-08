cmake_minimum_required(VERSION 2.8.3)
project(ras_particle_filter)

add_compile_options(-std=c++11)

find_package(catkin REQUIRED COMPONENTS roscpp geometry_msgs sensor_msgs)
find_package(Armadillo REQUIRED)

catkin_package(
  INCLUDE_DIRS include
  CATKIN_DEPENDS roscpp geometry_msgs sensor_msgs
#  DEPENDS system_lib
)

include_directories(
 include
 ${catkin_INCLUDE_DIRS}
 ${ARMADILLO_INCLUDE_DIRS}
)

add_executable(particle_filter src/particle_filter.cpp)

target_link_libraries(particle_filter
  ${catkin_LIBRARIES}
  -lopenblas
  -llapack
  ${ARMADILLO_LIBRARIES}
)

add_dependencies(particle_filter ras_particle_filter_generate_messages_cpp)
