cmake_minimum_required(VERSION 2.8.11)
project(main)

add_executable(${PROJECT_NAME} main.cpp)
add_definitions(-DPROJECT_NAME="${PROJECT_NAME}")
target_include_directories(${PROJECT_NAME} PRIVATE ${LIBIGL_INCLUDE_DIRS})
target_compile_definitions(${PROJECT_NAME} PRIVATE ${LIBIGL_DEFINITIONS})
target_link_libraries(${PROJECT_NAME} ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES})
