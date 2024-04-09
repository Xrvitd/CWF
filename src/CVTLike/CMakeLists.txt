cmake_minimum_required(VERSION 3.17)

set(BGAL_CVTLike_SRC        
		CPD.cpp
		CVT.cpp
		)
# Get static lib
find_package(OpenMP REQUIRED)
add_library(CVTLike STATIC ${BGAL_CVTLike_SRC})
target_link_libraries(CVTLike Algorithm BaseShape Model Tessellation2D Tessellation3D Optimization OpenMP::OpenMP_CXX ${Boost_LIBRARIES})
set_target_properties(CVTLike PROPERTIES VERSION ${VERSION})
set_target_properties(CVTLike PROPERTIES CLEAN_DIRECT_OUTPUT 1)

add_definitions(-DNOMINMAX)

target_include_directories(CVTLike PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)
