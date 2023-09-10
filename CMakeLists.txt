project(donut)

cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE Debug)

include(ExternalProject)

ExternalProject_Add(glm
    URL https://github.com/g-truc/glm/releases/download/0.9.9.8/glm-0.9.9.8.zip
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    CMAKE_ARGS -DCMAKE_CXX_FLAGS=-w -DCMAKE_C_FLAGS=-w
)

ExternalProject_Get_Property(glm SOURCE_DIR)
set(GLM_INCLUDE_DIR ${SOURCE_DIR})
include_directories(${GLM_INCLUDE_DIR})

add_executable(
    donut
    donut.cpp
    ${SOURCE_DIR}
)

target_compile_options(donut PRIVATE -fconstexpr-ops-limit=173554432)

ExternalProject_Add(glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG 3.3.4
    CMAKE_ARGS -DGLFW_BUILD_DOCS=OFF -DGLFW_BUILD_TESTS=OFF -DGLFW_BUILD_EXAMPLES=OFF
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(glfw source_dir binary_dir)
set(GLFW_INCLUDE_DIR ${source_dir}/include)
set(GLFW_LIBRARY ${binary_dir}/src/libglfw3.a)
include_directories(${GLFW_INCLUDE_DIR})
add_dependencies(donut glfw)
target_link_libraries(donut ${GLFW_LIBRARY})

find_package(X11 REQUIRED)
find_package(OpenGL REQUIRED)
include_directories(${X11_INCLUDE_DIR} ${OPENGL_INCLUDE_DIRS})
target_link_libraries(donut ${X11_LIBRARIES} ${OPENGL_LIBRARIES})
