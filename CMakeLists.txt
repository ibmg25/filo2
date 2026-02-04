cmake_minimum_required(VERSION 3.22)
project(filo2)

set(CMAKE_CXX_STANDARD 17)

get_filename_component(CONCORDE_DIR "${CMAKE_SOURCE_DIR}/../concorde" ABSOLUTE)
message(STATUS "Looking for Concorde in: ${CONCORDE_DIR}")

# Definimos TODAS las librerías necesarias.
# Concorde separa sus módulos en archivos .a distintos.
set(CONCORDE_LIBS
    "${CONCORDE_DIR}/concorde.a"
    "${CONCORDE_DIR}/EDGEGEN/edgegen.a"
    "${CONCORDE_DIR}/LINKERN/linkern.a"
    "${CONCORDE_DIR}/KDTREE/kdtree.a"
    "${CONCORDE_DIR}/UTIL/util.a"
    "${CONCORDE_DIR}/FMATCH/fmatch.a"
)

# Verificamos que existan (al menos la principal y edgegen que falló)
if(EXISTS "${CONCORDE_DIR}/concorde.a" AND EXISTS "${CONCORDE_DIR}/EDGEGEN/edgegen.a")
    message(STATUS "Found Concorde Libraries.")
else()
    message(FATAL_ERROR "Missing Concorde libraries.\n Expected: ${CONCORDE_DIR}/EDGEGEN/edgegen.a\n Make sure you ran 'make' inside the concorde folder correctly.")
endif()

set(WARNING_FLAGS "-Wall -Wextra -Wpedantic -Wuninitialized")
set(OPT_FLAGS "-O2 -fno-tree-vectorize")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fsanitize=undefined -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OPT_FLAGS}")

set(CMAKE_VERBOSE_MAKEFILE ON)

option(ENABLE_VERBOSE "Enable verbose output" ON)
option(ENABLE_GUI "Enable graphical interface" OFF)

message("-- Build options")

if(ENABLE_VERBOSE)
    message("--- Verbose output ENABLED")
    add_definitions(-DVERBOSE)
else()
    message("--- Verbose output DISABLED")
endif()

if(ENABLE_GUI)
    message("--- Graphical interface ENABLED")
    add_definitions(-DGUI)
    set(LIBRARIES ${LIBRARIES} glfw GL)
    set(SOURCE ${SOURCE} Renderer.hpp)
else()
    message("--- Graphical interface DISABLED")
endif()

add_executable(filo2

    base/BinaryHeap.hpp
    base/BitMatrix.hpp
    base/FixedSizeValueStack.hpp
    base/Flat2DVector.hpp
    base/functor.hpp
    base/KDTree.cpp
    base/KDTree.hpp
    base/LRUCache.hpp
    base/macro.hpp
    base/NonCopyable.hpp
    base/PrettyPrinter.hpp
    base/SmallFlatMap.hpp
    base/SmallFlatSet.hpp
    base/SparseIntSet.hpp
    base/Timer.hpp
    base/VectorView.hpp
    base/Welford.hpp

    instance/Instance.cpp
    instance/Instance.hpp
    instance/Parser.cpp
    instance/Parser.hpp

    localsearch/AbstractOperator.hpp
    localsearch/EjectionChain.hpp
    localsearch/LocalSearch.hpp
    localsearch/OneOneExchange.hpp
    localsearch/OneZeroExchange.hpp
    localsearch/RevThreeOneExchange.hpp
    localsearch/RevThreeThreeExchange.hpp
    localsearch/RevThreeTwoExchange.hpp
    localsearch/RevThreeZeroExchange.hpp
    localsearch/RevTwoOneExchange.hpp
    localsearch/RevTwoTwoExchange.hpp
    localsearch/RevTwoZeroExchange.hpp
    localsearch/SplitExchange.hpp
    localsearch/TailsExchange.hpp
    localsearch/ThreeOneExchange.hpp
    localsearch/ThreeThreeExchange.hpp
    localsearch/ThreeTwoExchange.hpp
    localsearch/ThreeZeroExchange.hpp
    localsearch/TwoOneExchange.hpp
    localsearch/TwoOptExchange.hpp
    localsearch/TwoTwoExchange.hpp
    localsearch/TwoZeroExchange.hpp

    movegen/MoveGenerators.hpp

    opt/routemin.hpp
    opt/RuinAndRecreate.hpp
    opt/SimulatedAnnealing.hpp

    solution/construction.hpp
    solution/construction.cpp
    solution/savings.hpp
    solution/Solution.hpp
    solution/Solution.cpp

    main.cpp
    Parameters.hpp
    Renderer.hpp
)

target_include_directories(filo2 PRIVATE
    ${CONCORDE_DIR}
    ${CONCORDE_DIR}/INCLUDE
    ${CONCORDE_DIR}/LINKERN
)

target_link_libraries(filo2 PUBLIC ${LIBRARIES} ${CONCORDE_LIBS} m)
