# polygon3split —— 验证测试工程(文档 docs/polygon-split-by-folded-surface.md §8)
#
# 手动 header-only 模式(不使用发行版 CGALConfig.cmake —— 它是打包机上生成的,
# 且依赖 config 模式的 Boost 查找;本机 Boost 是源码树+b2 构建产物,无 BoostConfig)。
# 直接给定三个依赖的 include/lib 路径,定义 CGAL_HEADER_ONLY。
#
# 用户环境(2026-08):
#   CGAL 6.2:  D:/github/CGAL-6.2/CGAL-6.2(include/ + auxiliary/gmp: 头文件/.lib/.dll)
#   Boost 1.91: D:/github/boost_1_91_0(纯头文件使用,PMP 用到的模块均为 header-only)
#   工具链: CMake 4.3 + Visual Studio 18 2026(MSVC v145)
#
# 路径可用 -D 覆盖。

cmake_minimum_required(VERSION 3.20)
project(polygon3split_tests CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CGAL_DIR  "D:/github/CGAL-6.2/CGAL-6.2"   CACHE PATH "CGAL 6.2 发行版根目录")
set(BOOST_DIR "D:/github/boost_1_91_0"        CACHE PATH "Boost 源码树(含 boost/ 头文件目录)")
set(GMP_DIR   "${CGAL_DIR}/auxiliary/gmp"     CACHE PATH "GMP/MPFR(发行版自带)")

add_executable(split_test tests/split_test.cpp)
target_include_directories(split_test PRIVATE
  src
  "${CGAL_DIR}/include"
  "${BOOST_DIR}"
  "${GMP_DIR}/include"
)
target_compile_definitions(split_test PRIVATE
  CGAL_HEADER_ONLY=1          # 无预编译 CGAL 库,全部头文件方式
)
target_link_libraries(split_test PRIVATE
  "${GMP_DIR}/lib/gmp.lib"
  "${GMP_DIR}/lib/mpfr.lib"
)

if(MSVC)
  target_compile_options(split_test PRIVATE /bigobj /Zc:__cplusplus /permissive- /utf-8)
endif()

# 运行时需要 auxiliary/gmp/bin 下的 gmp-10.dll / mpfr-6.dll,复制到输出目录
add_custom_command(TARGET split_test POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${GMP_DIR}/bin/gmp-10.dll"
    "${GMP_DIR}/bin/mpfr-6.dll"
    $<TARGET_FILE_DIR:split_test>)
