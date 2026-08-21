cmake_minimum_required(VERSION 3.20)
project(triang_repro CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CGAL_DIR  "D:/github/CGAL-6.2/CGAL-6.2"   CACHE PATH "CGAL 6.2")
set(BOOST_DIR "D:/github/boost_1_91_0"        CACHE PATH "Boost 源码树")
set(GMP_DIR   "${CGAL_DIR}/auxiliary/gmp"     CACHE PATH "GMP/MPFR")

add_executable(triang_repro triang_repro.cpp)
add_executable(gate_probe gate_probe.cpp)
target_include_directories(triang_repro PRIVATE
  "${CGAL_DIR}/include"
  "${BOOST_DIR}"
  "${GMP_DIR}/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/../src"
)
target_compile_definitions(triang_repro PRIVATE CGAL_HEADER_ONLY=1)
target_link_libraries(triang_repro PRIVATE
  "${GMP_DIR}/lib/gmp.lib"
  "${GMP_DIR}/lib/mpfr.lib"
)
if(MSVC)
  target_compile_options(triang_repro PRIVATE /bigobj /Zc:__cplusplus /permissive- /utf-8)
endif()
add_custom_command(TARGET triang_repro POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${GMP_DIR}/bin/gmp-10.dll"
    "${GMP_DIR}/bin/mpfr-6.dll"
    $<TARGET_FILE_DIR:triang_repro>)

foreach(t gate_probe)
  target_include_directories(${t} PRIVATE
    "${CGAL_DIR}/include" "${BOOST_DIR}" "${GMP_DIR}/include" "${CMAKE_CURRENT_SOURCE_DIR}/../src")
  target_compile_definitions(${t} PRIVATE CGAL_HEADER_ONLY=1)
  target_link_libraries(${t} PRIVATE "${GMP_DIR}/lib/gmp.lib" "${GMP_DIR}/lib/mpfr.lib")
  if(MSVC)
    target_compile_options(${t} PRIVATE /bigobj /Zc:__cplusplus /permissive- /utf-8)
  endif()
endforeach()
