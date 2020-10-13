project(CFDSim)
set(CMAKE_CXX_FLAGS "-W -Wall -pedantic -std=c++11")

#headers are added as executables for displaying purpose in the file manager while using xcode IDE.
SET(SOURCE src/main.cpp src/Exporter.cpp src/Simulator.cpp)
SET(HEADERS src/Exporter.hpp src/Simulator.hpp)
add_executable(CFDSim ${SOURCE} ${HEADERS})
