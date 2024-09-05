cmake_minimum_required(VERSION 3.4...3.18)
project(pyDTFE)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_PREFIX_PATH "/home/shuren/Code/CGAL-5.6.1/MyInstall")
# include_directories("/home/shuren/.local/lib/python3.9/site-packages/pybind11/include")
# include_directories("/home/shuren/Code/CGAL-5.6.1/MyInstall")
# set(Python3_EXECUTABLE "/home/shuren/.conda/envs/shuren_env3.9/bin")
find_package(Python REQUIRED COMPONENTS Interpreter Development)
# set(BOOST_ROOT "/path/to/Boost")
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_C_COMPILER gcc)

find_package(pybind11 REQUIRED)
find_package(CGAL REQUIRED)
find_package(MPI REQUIRED)

add_library(src_dtfe MODULE ./src/main.cpp)

target_link_libraries(src_dtfe PRIVATE MPI::MPI_CXX pybind11::module)
target_link_libraries(src_dtfe PRIVATE CGAL::CGAL pybind11::module)
target_include_directories(src_dtfe PRIVATE
    ${Python3_INCLUDE_DIRS}
    ${CGAL_INCLUDE_DIRS}
    ${pybind11_INCLUDE_DIRS}
    ${MPI_INCLUDE_PATH}
)

# Specify the C++ standard
# target_compile_features(src_dtfe PRIVATE cxx_std_11)

#Set the output to be a Python module
set_target_properties(src_dtfe PROPERTIES
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
)



