SET(GraphMiningSrc
graph.cpp
labeled_graph.cpp
dataloader.cpp
pattern.cpp
prefix.cpp
#schedule.cpp
schedule_IEP.cpp
vertex_set.cpp
motif_generator.cpp
graphmpi.cpp
common.cpp
disjoint_set_union.cpp
set_operation.cpp
)


link_directories(${CMAKE_SOURCE_DIR}/../libs)

ADD_LIBRARY(graph_mining SHARED ${GraphMiningSrc}) 

# ADD_EXECUTABLE(fsm_test fsm_test.cpp)
# TARGET_LINK_LIBRARIES(fsm_test graph_mining)

# ADD_EXECUTABLE(fsm_vertex_test fsm_vertex_test.cpp)
# TARGET_LINK_LIBRARIES(fsm_vertex_test graph_mining)

# ADD_EXECUTABLE(clique_test clique_test.cpp)
# TARGET_LINK_LIBRARIES(clique_test graph_mining)

ADD_EXECUTABLE(motif_counting_test motif_counting_test.cpp)
TARGET_LINK_LIBRARIES(motif_counting_test PRIVATE -Wall ${OpenMP_CXX_FLAGS} graph_mining)

# ADD_EXECUTABLE(pm_test pm_test.cpp)
# TARGET_LINK_LIBRARIES(pm_test graph_mining)

ADD_EXECUTABLE(approx approx.cpp)
TARGET_LINK_LIBRARIES(approx graph_mining)

