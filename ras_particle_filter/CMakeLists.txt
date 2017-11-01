cmake_minimum_required(VERSION 2.8.3)
project(ras_particle_filter)

add_compile_options(-std=c++11)

find_package(catkin REQUIRED COMPONENTS roscpp geometry_msgs message_generation ras_line_detector)
find_package(Armadillo REQUIRED)

add_message_files(
  FILES
  Particles.msg
)

generate_messages(DEPENDENCIES geometry_msgs)

catkin_package(
  INCLUDE_DIRS include
  CATKIN_DEPENDS roscpp geometry_msgs message_runtime ras_line_detector
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
