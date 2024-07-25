include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB_RECURSE SOURCES "*.c*")
file(GLOB_RECURSE HEADERS "*.h*")

list(REMOVE_ITEM SOURCES "main*.cc")

add_library(DCRP STATIC ${SOURCES} ${HEADERS})
add_executable(DCRPVDLNS main_vdlns.cc)
target_link_libraries(DCRPVDLNS DCRP)