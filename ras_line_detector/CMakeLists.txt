cmake_minimum_required(VERSION 2.8.3)
project(ras_line_detector)

find_package(catkin REQUIRED COMPONENTS rospy sensor_msgs message_generation geometry_msgs)

add_message_files(
  FILES
  LineSegment.msg
  LineSegmentList.msg
)

generate_messages(
  DEPENDENCIES
  sensor_msgs
)

catkin_package(CATKIN_DEPENDS rospy sensor_msgs message_runtime geometry_msgs)

catkin_install_python(
  PROGRAMS scripts/line_detector.py scripts/draw_lines.py
  DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
