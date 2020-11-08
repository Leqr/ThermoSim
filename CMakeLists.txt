project(CFDSim)
set(CMAKE_CXX_FLAGS "-O3 -W -Wall -pedantic -std=c++11")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

#headers are added as executables for displaying purpose in the file manager while using xcode IDE.
SET(SOURCE src/main.cpp src/Exporter.cpp src/Simulator.cpp)
SET(HEADERS src/Exporter.hpp src/Simulator.hpp)
add_executable(CFDSim ${SOURCE} ${HEADERS})

target_link_libraries (CFDSim Eigen3::Eigen)
