#-------------------- CMAKE Configuration
cmake_minimum_required(VERSION 3.21)
project(asdf2fits VERSION 0.1.0 LANGUAGES C CXX)

if(NOT WIN32)
    string(ASCII 27 Esc)
    set(ColourReset "${Esc}[m")
    set(ColourBold  "${Esc}[1m")
    set(Red         "${Esc}[31m")
    set(Green       "${Esc}[32m")
    set(Yellow      "${Esc}[33m")
    set(Blue        "${Esc}[34m")
    set(Magenta     "${Esc}[35m")
    set(Cyan        "${Esc}[36m")
    set(White       "${Esc}[37m")
    set(BoldRed     "${Esc}[1;31m")
    set(BoldGreen   "${Esc}[1;32m")
    set(BoldYellow  "${Esc}[1;33m")
    set(BoldBlue    "${Esc}[1;34m")
    set(BoldMagenta "${Esc}[1;35m")
    set(BoldCyan    "${Esc}[1;36m")
    set(BoldWhite   "${Esc}[1;37m")
endif()

message("${BoldCyan}----------- asdf2fits ${PROJECT_VERSION} -----------${ColourReset}")
message("${Green}-- Project source directory :${ColourReset} ${CMAKE_SOURCE_DIR}")

if(DEFINED ENV{DeepSkyROOT})
    set(DeepSkyROOT $ENV{DeepSkyROOT})
else()
    set(DeepSkyROOT ${CMAKE_SOURCE_DIR}/..)
endif()
get_filename_component(DeepSkyROOT "${DeepSkyROOT}" ABSOLUTE)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${DeepSkyROOT} CACHE PATH "default install path" FORCE)
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Same language level as the DeepSky libraries whose headers are compiled here.
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED OFF)
set(CMAKE_CXX_EXTENSIONS OFF)

option(WITH_ASAN "Build asdf2fits and its tests with AddressSanitizer" OFF)
option(BUILD_TESTS "Build the GoogleTest suite" ON)

message("${BoldYellow}-- CMAKE INSTALL PREFIX : ${ColourReset}${CMAKE_INSTALL_PREFIX}")
message("${Green}-- DeepSkyROOT :${ColourReset} ${DeepSkyROOT}")
message("${Green}-- BUILD TYPE :${ColourReset} ${CMAKE_BUILD_TYPE}")
message("${Green}-- CXX COMPILER :${ColourReset} ${CMAKE_CXX_COMPILER}")
message("${Green}-- USE CXX STANDARD :${ColourReset} ${CMAKE_CXX_STANDARD}")
if(WITH_ASAN)
    message("${Blue}-- AddressSanitizer ACTIVATED${ColourReset}")
endif()

foreach(sub libasdf/CMakeLists.txt libfyaml/CMakeLists.txt libmd/configure.ac)
    if(NOT EXISTS ${CMAKE_SOURCE_DIR}/${sub})
        message(FATAL_ERROR "${BoldRed}${sub} is missing: run 'git submodule update --init'${ColourReset}")
    endif()
endforeach()

#-------------------- Dependencies built from the submodules
message("${Red}----------- SUBMODULES : libmd, libfyaml, libasdf -----------${ColourReset}")
include(ExternalProject)

# libmd and libfyaml are linked statically into libasdf and stay in this private stage directory;
# only libasdf itself and its headers reach the install prefix.
set(STAGE_DIR ${CMAKE_BINARY_DIR}/stage)
file(MAKE_DIRECTORY ${STAGE_DIR}/include ${STAGE_DIR}/lib)

find_program(GLIBTOOLIZE glibtoolize)
if(GLIBTOOLIZE)
    set(AUTORECONF_ENV LIBTOOLIZE=${GLIBTOOLIZE})
    message("-- ${Green}Found${ColourReset} glibtoolize: ${GLIBTOOLIZE}")
endif()

message("${Blue}-- libmd    ${ColourReset} static, kept in ${STAGE_DIR}")
message("${Blue}-- libfyaml ${ColourReset} static, kept in ${STAGE_DIR}")
message("${Blue}-- libasdf  ${ColourReset} shared, installed in ${CMAKE_INSTALL_PREFIX}/lib")

# autogen writes into its source tree: work on a copy so that the submodule stays clean.
# The copy has no usable .git, so the version is given through .dist-version.
set(LIBMD_COPY ${CMAKE_BINARY_DIR}/external/libmd-src)
ExternalProject_Add(libmd_external
    PREFIX            ${CMAKE_BINARY_DIR}/external/libmd
    DOWNLOAD_COMMAND  ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/libmd ${LIBMD_COPY}
    COMMAND           ${CMAKE_COMMAND} -E rm -f ${LIBMD_COPY}/.git
    SOURCE_DIR        ${LIBMD_COPY}
    PATCH_COMMAND     ${CMAKE_COMMAND} -DOUT=${LIBMD_COPY}/.dist-version -P ${CMAKE_SOURCE_DIR}/cmake/WriteLibmdVersion.cmake
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E chdir ${LIBMD_COPY} ${CMAKE_COMMAND} -E env ${AUTORECONF_ENV} ./autogen
    COMMAND           ${CMAKE_COMMAND} -E chdir ${LIBMD_COPY} ./configure --prefix=${STAGE_DIR} --libdir=${STAGE_DIR}/lib --enable-static --disable-shared --with-pic
    BUILD_COMMAND     ${CMAKE_COMMAND} -E chdir ${LIBMD_COPY} make
    INSTALL_COMMAND   ${CMAKE_COMMAND} -E chdir ${LIBMD_COPY} make install
    BUILD_BYPRODUCTS  ${STAGE_DIR}/lib/libmd.a
)

ExternalProject_Add(libfyaml_external
    PREFIX     ${CMAKE_BINARY_DIR}/external/libfyaml
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/libfyaml
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_INSTALL_PREFIX=${STAGE_DIR}
        -DCMAKE_INSTALL_LIBDIR=lib
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5
        -DBUILD_SHARED_LIBS=OFF
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DBUILD_TESTING=OFF
        -DENABLE_NETWORK=OFF
        -DENABLE_REFLECTION=OFF
        -DENABLE_PYTHON_BINDINGS=OFF
    BUILD_BYPRODUCTS ${STAGE_DIR}/lib/libfyaml.a
)

set(LIBASDF_NAME ${CMAKE_SHARED_LIBRARY_PREFIX}asdf${CMAKE_SHARED_LIBRARY_SUFFIX})

# libasdf is built from a patched copy, until the fix is available upstream. The hash in the command
# makes the copy start again whenever the patch changes.
find_program(PATCH_PROGRAM patch REQUIRED)
set(LIBASDF_COPY ${CMAKE_BINARY_DIR}/external/libasdf-src)
set(LIBASDF_PATCH ${CMAKE_SOURCE_DIR}/cmake/patches/libasdf-document-end-marker.patch)
file(SHA256 ${LIBASDF_PATCH} LIBASDF_PATCH_HASH)
message("${Blue}-- libasdf patch ${ColourReset} cmake/patches/libasdf-document-end-marker.patch")
ExternalProject_Add(libasdf_external
    DEPENDS    libmd_external libfyaml_external
    PREFIX     ${CMAKE_BINARY_DIR}/external/libasdf
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -DSRC=${CMAKE_SOURCE_DIR}/libasdf -DDEST=${LIBASDF_COPY}
        -DPATCH_FILE=${LIBASDF_PATCH} -DPATCH_PROGRAM=${PATCH_PROGRAM} -DPATCH_HASH=${LIBASDF_PATCH_HASH}
        -P ${CMAKE_SOURCE_DIR}/cmake/PrepareLibasdf.cmake
    SOURCE_DIR ${LIBASDF_COPY}
    BINARY_DIR ${CMAKE_BINARY_DIR}/external/libasdf-build
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_INSTALL_PREFIX=${STAGE_DIR}
        -DCMAKE_INSTALL_LIBDIR=lib
        -DCMAKE_INSTALL_NAME_DIR=@rpath
        -DENABLE_TOOL=OFF
        -DUSE_STATGRAB=OFF
        -DENABLE_TESTING=OFF
        -DFYAML_NO_PKGCONFIG=YES
        -DFYAML_LIBDIR=${STAGE_DIR}/lib
        -DFYAML_INCLUDEDIR=${STAGE_DIR}/include
        -DCMAKE_C_FLAGS=-I${STAGE_DIR}/include
        -DCMAKE_EXE_LINKER_FLAGS=-L${STAGE_DIR}/lib
        -DCMAKE_SHARED_LINKER_FLAGS=-L${STAGE_DIR}/lib
    BUILD_BYPRODUCTS ${STAGE_DIR}/lib/${LIBASDF_NAME}
)

add_library(asdf_c SHARED IMPORTED)
set_target_properties(asdf_c PROPERTIES IMPORTED_LOCATION ${STAGE_DIR}/lib/${LIBASDF_NAME})
add_dependencies(asdf_c libasdf_external)

#-------------------- DeepSkyFits
message("${Red}----------- SEARCH DeepSkyFits -----------${ColourReset}")
set(DSTfits_MIN_VERSION 1.5)
# Only look in DeepSkyROOT: the RPATH points there.
find_package(DSTfits ${DSTfits_MIN_VERSION} REQUIRED PATHS ${DeepSkyROOT} NO_DEFAULT_PATH)
message("-- ${Green}Found${ColourReset} DSTfits ${DSTfits_VERSION}: ${DSTfits_DIR}")

set(DEEPSKY_INCLUDES ${DeepSkyROOT}/include ${DeepSkyROOT}/include/cfitsio ${DeepSkyROOT}/include/Minuit2)
set(DEEPSKY_RPATH    "${DeepSkyROOT}/lib")
# This build comes first: an installed copy in DeepSkyROOT must not replace what the tests are testing.
set(CMAKE_BUILD_RPATH "${CMAKE_BINARY_DIR};${STAGE_DIR}/lib;${DEEPSKY_RPATH}")
set(CMAKE_INSTALL_RPATH "${DEEPSKY_RPATH}")

#-------------------- libasdf2fits
message("${Red}----------- libasdf2fits -----------${ColourReset}")
set(ASDF2FITS_SOURCES
    src/ASDFmanager.cxx
    src/ASDFheader.cxx
    src/ASDFimage.cxx
    src/ASDFtable.cxx
    src/fits2asdf.cxx
)
set(ASDF2FITS_PUBLIC_HEADERS src/ASDFmanager.h src/fits2asdf.h)

if(EXISTS ${CMAKE_SOURCE_DIR}/src/ASDFmanager.cxx)
    add_library(asdf2fits SHARED ${ASDF2FITS_SOURCES})
    add_dependencies(asdf2fits libasdf_external)
    set_target_properties(asdf2fits PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
    )
    target_include_directories(asdf2fits
        PUBLIC  ${CMAKE_SOURCE_DIR}/src ${DEEPSKY_INCLUDES}
        PRIVATE ${STAGE_DIR}/include
    )
    target_link_libraries(asdf2fits PUBLIC DSTfits::DSTfits PRIVATE asdf_c)
    message("${Blue}-- Public headers ${ColourReset} ${ASDF2FITS_PUBLIC_HEADERS}")
    if(WITH_ASAN)
        target_compile_options(asdf2fits PUBLIC -fsanitize=address -fno-omit-frame-pointer -g)
        target_link_options(asdf2fits PUBLIC -fsanitize=address)
    endif()

    install(TARGETS asdf2fits LIBRARY DESTINATION lib)
    install(FILES ${ASDF2FITS_PUBLIC_HEADERS} DESTINATION include/asdf2fits)

    include(CMakePackageConfigHelpers)
    configure_package_config_file(
        ${CMAKE_SOURCE_DIR}/cmake/asdf2fitsConfig.cmake.in
        ${CMAKE_BINARY_DIR}/generated/asdf2fitsConfig.cmake
        INSTALL_DESTINATION lib/cmake/asdf2fits
    )
    write_basic_package_version_file(
        ${CMAKE_BINARY_DIR}/generated/asdf2fitsConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
    )
    install(FILES
        ${CMAKE_BINARY_DIR}/generated/asdf2fitsConfig.cmake
        ${CMAKE_BINARY_DIR}/generated/asdf2fitsConfigVersion.cmake
        DESTINATION lib/cmake/asdf2fits
    )
else()
    message("${Yellow}-- src/ASDFmanager.cxx not found: libasdf2fits is skipped${ColourReset}")
endif()

# libasdf and its public headers (the libmd and libfyaml archives are already inside the library)
install(DIRECTORY ${STAGE_DIR}/lib/ DESTINATION lib FILES_MATCHING PATTERN "libasdf*${CMAKE_SHARED_LIBRARY_SUFFIX}*" PATTERN "pkgconfig" EXCLUDE PATTERN "cmake" EXCLUDE)
install(FILES ${STAGE_DIR}/include/asdf.h DESTINATION include)
install(DIRECTORY ${STAGE_DIR}/include/asdf DESTINATION include)

#-------------------- Tests
if(BUILD_TESTS AND TARGET asdf2fits AND EXISTS ${CMAKE_SOURCE_DIR}/test/gTest_asdf2fits.cxx)
    message("${Red}----------- GoogleTest -----------${ColourReset}")
    set(ASDF2FITS_TESTS ON)
    enable_testing()
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
    add_subdirectory(googletest EXCLUDE_FROM_ALL)
    include(GoogleTest)

    add_executable(gTest_asdf2fits test/gTest_asdf2fits.cxx)
    target_link_libraries(gTest_asdf2fits asdf2fits asdf_c gtest gtest_main)
    target_include_directories(gTest_asdf2fits PRIVATE ${STAGE_DIR}/include)
    target_compile_definitions(gTest_asdf2fits PRIVATE
        ASDF_FIXTURES="${CMAKE_SOURCE_DIR}/libasdf/tests/fixtures")
    # A timeout, as a parser bug in libasdf shows up as a test that never ends
    gtest_discover_tests(gTest_asdf2fits WORKING_DIRECTORY ${CMAKE_BINARY_DIR} PROPERTIES TIMEOUT 120)
    message("${Blue}-- Test executable ${ColourReset} gTest_asdf2fits (run with ctest)")
else()
    set(ASDF2FITS_TESTS OFF)
endif()

#-------------------- Summary
message("${BoldGreen}----------- asdf2fits ${PROJECT_VERSION} configured -----------${ColourReset}")
message("${Green}-- Build type       :${ColourReset} ${CMAKE_BUILD_TYPE}")
message("${Green}-- AddressSanitizer :${ColourReset} ${WITH_ASAN}")
message("${Green}-- Tests            :${ColourReset} ${ASDF2FITS_TESTS}")
message("${Green}-- DSTfits          :${ColourReset} ${DSTfits_VERSION}")
message("${Green}-- Install prefix   :${ColourReset} ${CMAKE_INSTALL_PREFIX}")
message("${BoldGreen}-------------------------------------------------${ColourReset}")
