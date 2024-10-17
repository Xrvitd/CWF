cmake_minimum_required(VERSION 3.17)

add_executable(MAIN main.cpp)
target_link_libraries(MAIN Algorithm BaseShape Draw Geodesic Integral Model Optimization PointCloudProcessing PQP Reconstruction Tessellation2D Tessellation3D CVTLike)
target_include_directories(MAIN PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)


add_executable(cwf cwf.cpp)
target_link_libraries(cwf Algorithm BaseShape Draw Geodesic Integral Model Optimization PointCloudProcessing PQP Reconstruction Tessellation2D Tessellation3D CVTLike)
target_include_directories(cwf PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)
