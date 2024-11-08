cmake_minimum_required(VERSION 3.1)
project(Project_3)
set(CMAKE_CXX_STANDARD 11)

include_directories(Core)

set(SOURCES Core/Edge.cpp Core/Face.cpp Core/FException.cpp Core/HalfEdge.cpp 
			Core/OBJFileReader.cpp Core/Parser.cpp Core/Point.cpp Core/Solid.cpp Core/SolidDelegate.cpp 
			Core/StringTokenizer.cpp Core/TopologyException.cpp Core/Trait.cpp Core/Vertex.cpp)

add_executable(main main.cpp ${SOURCES})
