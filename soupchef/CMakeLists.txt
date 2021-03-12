cmake_minimum_required(VERSION 3.17)
project(geo1004.soupchef)

add_executable(soupchef 
  main.cpp
)
set_target_properties(soupchef PROPERTIES 
  CXX_STANDARD 17
)

# Make the executable a default target to build & run in Visual Studio
set_property(DIRECTORY ${PROJECT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT soupchef)

