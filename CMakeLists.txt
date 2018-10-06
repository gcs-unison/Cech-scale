cmake_minimum_required (VERSION 2.6)
project (Cech-scale)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -lm")
#set(CMAKE_BINARY_DIR ${PROJECT_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

#sets optimization for release build
#compile: cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
if(CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
    set(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS})
endif()

#sets debug flag for debug build
#compile: cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Debug
if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
    set(CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS})
endif()

#when cross compiling to windows, add flags to run correctly
#cross compile windows: cmake -H. -Bwindows -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=toolchains/Toolchain-mingw32.cmake
if(CMAKE_TOOLCHAIN_FILE MATCHES Toolchain-mingw32.cmake)
    set(LINKER_FLAGS "-static -static-libgcc -static-libstdc++")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LINKER_FLAGS}")
endif()

#adds header folders
include_directories(header)

#adds all source files
file(GLOB SOURCES "src/*.cpp")

#indicates resulting file
add_executable(Cech-scale ${SOURCES})

#install(TARGETS cech DESTINATION ${PROJECT_SOURCE_DIR})
