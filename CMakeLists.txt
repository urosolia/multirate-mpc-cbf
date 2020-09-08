cmake_minimum_required (VERSION 3.10)
project(multirate-mpc-cbf VERSION 1.0 LANGUAGES CXX)

## Compile as C++17, supported in ROS Kinetic and newer
add_compile_options(-std=c++17)

# Include utilities
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Set various flags
set(default_build_type "Release")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3 -DNDEBUG")

# OSQP version option
option(OSQP_EMBEDDED "Use embedded version of osqp." OFF)

# Chose osqp version
if(OSQP_EMBEDDED)
  message(STATUS "Using embedded version of osqp.")
  set(OSQP_MODE "osqp_embedded")
  set(OSQP_MODE_LIB "osqp_embedded")
else()
  set(OSQP_MODE "osqp")
  set(OSQP_MODE_LIB "osqpstatic")
endif()

# Release type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Find OSQP library and headers
find_package(${OSQP_MODE} REQUIRED CONFIG)

# Add sources
set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/cbf_test.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/mpc_test.cpp")
file(GLOB "${PROJECT_NAME}_HEADERS" "${CMAKE_CURRENT_SOURCE_DIR}/include/mpc/*.h*")

# Make library
add_library(${PROJECT_NAME} STATIC ${SRC} "${${PROJECT_NAME}_HEADERS}")
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

# Set include directory
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                                  "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>")

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${${PROJECT_NAME}_HEADERS}")

# Link with other libraries
target_link_libraries(${PROJECT_NAME} PUBLIC ${OSQP_MODE}::${OSQP_MODE_LIB}
                                             -ldl)

# Set installation options
# ----------------------------------------------

install(TARGETS ${PROJECT_NAME} 
        EXPORT  ${PROJECT_NAME}Targets
        ARCHIVE DESTINATION "lib/${PROJECT_NAME}"
        LIBRARY DESTINATION "lib/${PROJECT_NAME}"
        RUNTIME DESTINATION "bin/${PROJECT_NAME}"
        PUBLIC_HEADER DESTINATION "include/${PROJECT_NAME}")

install(EXPORT "${PROJECT_NAME}Targets"
        FILE "${PROJECT_NAME}Targets.cmake"
        NAMESPACE "${PROJECT_NAME}::"
        DESTINATION "share/${PROJECT_NAME}")

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION share/${PROJECT_NAME}/cmake
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION "share/${PROJECT_NAME}")

# uninstall target
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake)