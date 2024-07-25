include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
    cyclecovering.cc
    dcrp.cc
    point.cc
    random.cc
    shortestpath.cc
    simplegraph.cc
    tsp.cc
)
set(HEADERS
    argheap.hpp
    cyclecovering.hh
    dcrp.hh
    point.hh
    random.hh
    shortestpath.hh
    simplegraph.hh
    tsp.hh
)

add_library(DCRP STATIC ${SOURCES} ${HEADERS})

add_executable(DCRPVDLNS main_vdlns.cc)
target_link_libraries(DCRPVDLNS DCRP)

if(USE_GUROBI)
    add_library(DCRPGRB_LIB STATIC cyclecovering_grb.hh cyclecovering_grb.cc)
    target_link_libraries(DCRPGRB_LIB PUBLIC DCRP ${GUROBI_LIBRARIES})
    target_include_directories(DCRPGRB_LIB PUBLIC "${GUROBI_INCLUDE_DIR}")

    add_executable(DCRPGRB main_grb.cc)
    target_link_libraries(DCRPGRB DCRPGRB_LIB ${GUROBI_LIBRARIES})
    target_include_directories(DCRPGRB PUBLIC "${GUROBI_INCLUDE_DIR}")
endif(USE_GUROBI)
