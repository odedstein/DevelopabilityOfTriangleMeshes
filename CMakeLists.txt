cmake_minimum_required(VERSION 2.6)
project(developableflow)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(LIBIGL QUIET)

if (NOT LIBIGL_FOUND)
message(FATAL_ERROR "libigl not found --- You can download it using: \n git clone --recursive https://github.com/libigl/libigl.git ${PROJECT_SOURCE_DIR}/../libigl")
endif()

# Compilation flags: adapt to your needs
if(MSVC)
# Enable parallel compilation
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
else()
# Libigl requires a modern C++ compiler that supports c++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "." )
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")

# libigl options: choose between header only and compiled static library
# Header-only is preferred for small projects. For larger projects the static build
# considerably reduces the compilation times
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" OFF)

# add a customizable menu bar
option(LIBIGL_WITH_NANOGUI     "Use Nanogui menu"   ON)

# libigl options: choose your dependencies (by default everything is OFF except opengl)
option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_OPENGL      "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW"           ON)
option(LIBIGL_WITH_BBW         "Use BBW"            OFF)
option(LIBIGL_WITH_EMBREE      "Use Embree"         OFF)
option(LIBIGL_WITH_PNG         "Use PNG"            OFF)
option(LIBIGL_WITH_TETGEN      "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE    "Use Triangle"       OFF)
option(LIBIGL_WITH_XML         "Use XML"            OFF)
option(LIBIGL_WITH_LIM         "Use LIM"            OFF)
option(LIBIGL_WITH_COMISO      "Use CoMiso"         OFF)
option(LIBIGL_WITH_MATLAB      "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK       "Use MOSEK"          OFF) # This option is not supported yet
option(LIBIGL_WITH_CGAL        "Use CGAL"           OFF)
if(LIBIGL_WITH_CGAL) # Do not remove or move this block, the cgal build system fails without it
find_package(CGAL REQUIRED)
set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "CGAL's CMAKE Setup is super annoying ")
include(${CGAL_USE_FILE})
endif()

# Adding libigl: choose the path to your local copy libigl
# This is going to compile everything you requested
#message(FATAL_ERROR "${PROJECT_SOURCE_DIR}/../libigl/cmake")
add_subdirectory("${LIBIGL_INCLUDE_DIR}/../shared/cmake" "libigl")

#Add MPFR
#find_package(MPFR)
#if(MPFR_FOUND)
#string(STRIP ${MPFR_INCLUDES} MPFR_INCLUDES)
#string(STRIP ${MPFR_LIBRARIES} MPFR_LIBRARIES)
#message("ALIB: ${MPFR_LIBRARIES}")
#set(LIBIGL_INCLUDE_DIRS "${LIBIGL_INCLUDE_DIRS};${MPFR_INCLUDES};${LIBIGL_INCLUDE_DIR}/../external/nanogui/ext/eigen/unsupported/test/mpreal/")
#set(LIBIGL_LIBRARIES "${LIBIGL_EXTRA_LIBRARIES};${MPFR_LIBRARIES}")
#message("-- MPFR found")
#endif()

#Add BFF (Thank you, Rohan)
#find_package(SuiteSparse REQUIRED)
#set(LIBIGL_INCLUDE_DIRS "${LIBIGL_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}/bff/mesh/include;${PROJECT_SOURCE_DIR}/bff/linear-algebra/include;${PROJECT_SOURCE_DIR}/bff/project/include;${PROJECT_SOURCE_DIR}/bff/deps/nanogui/include;${PROJECT_SOURCE_DIR}/bff/viewer/include")
#set(LIBIGL_LIBRARIES "${LIBIGL_LIBRARIES};${PROJECT_SOURCE_DIR}/bff/build/libbff.a;${CHOLMOD_LIBRARY};${LAPACK_LIBRARIES}")

# If static, make sure to compile everything in include
if(LIBIGL_USE_STATIC_LIBRARY)
file(GLOB SOURCES_IGL
"${PROJECT_SOURCE_DIR}/include/developableflow/*.cpp"
"${PROJECT_SOURCE_DIR}/include/tools/*.cpp"
"${PROJECT_SOURCE_DIR}/include/viewer/*.cpp")
add_library(${PROJECT_NAME} STATIC ${SOURCES_IGL})
target_include_directories(${PROJECT_NAME} PRIVATE "${LIBIGL_INCLUDE_DIRS}")
target_compile_definitions(${PROJECT_NAME} PRIVATE "-DIGL_STATIC_LIBRARY -DPROJECT_NAME=\"${PROJECT_NAME}\"")
list(APPEND LIBIGL_LIBRARIES "${PROJECT_NAME}")
endif()

# libigl information
message("libigl includes: ${LIBIGL_INCLUDE_DIRS}")
message("libigl libraries: ${LIBIGL_LIBRARIES}")
message("libigl extra sources: ${LIBIGL_EXTRA_SOURCES}")
message("libigl extra libraries: ${LIBIGL_EXTRA_LIBRARIES}")
message("libigl definitions: ${LIBIGL_DEFINITIONS}")

# Prepare the build environment
include_directories(${LIBIGL_INCLUDE_DIRS})
add_definitions(${LIBIGL_DEFINITIONS})

# Add your project files and includes
include_directories("include")


# Add all subdirectories that should be compiled
add_subdirectory("applications/main")
