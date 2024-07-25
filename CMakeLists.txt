cmake_minimum_required(VERSION 3.18.0)
project(DetailedCableRouting VERSION 1.0.0)
set(PROJECT_NAMESPACE dcrp)

option(USE_GUROBI "use gurobi" OFF)

# unifiy all system output path 
set(OUTPUT_ITEMS RUNTIME LIBRARY ARCHIVE)
set(OUTPUT_MODES RELEASE DEBUG)
foreach(OUTPUT_ITEM ${OUTPUT_ITEMS})
    set(CMAKE_${OUTPUT_ITEM}_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
    foreach(OUTPUT_MODE ${OUTPUT_MODES})
        set(CMAKE_${OUTPUT_ITEM}_OUTPUT_DIRECTORY_${OUTPUT_MODE} ${CMAKE_CURRENT_SOURCE_DIR}/bin)
    endforeach()
endforeach()

if(MSVC)
    set(USE_GUROBI OFF)
    list(APPEND CMAKE_C_FLAGS "/utf-8 /EHsc")
    list(APPEND CMAKE_CXX_FLAGS "/utf-8 /EHsc")
endif()

if(USE_GUROBI)
    list(APPEND CMAKE_MODULE_PATH  "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
    include(FindGurobi)
    include_directories(${GUROBI_INCLUDE_DIRS})
endif()

add_subdirectory(src)