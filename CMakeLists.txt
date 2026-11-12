#-------------------- CMAKE Configuration
cmake_minimum_required(VERSION 3.21)

# Everything is generated in a build directory, never in the sources: from the project directory,
# 'cmake --preset default' (or 'cmake -B build'); from a build directory, 'cmake <path to project>'.
get_filename_component(_src "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component(_bin "${CMAKE_BINARY_DIR}" REALPATH)
if(_src STREQUAL _bin)
    message(FATAL_ERROR
        "In-source builds are not allowed. Remove what this run created:\n"
        "    rm -rf CMakeCache.txt CMakeFiles\n"
        "then, from the project directory:\n"
        "    cmake --preset default        (or: cmake -B build)\n"
        "    cmake --build build -j\n"
        "or from a build directory: mkdir build && cd build && cmake .. && make")
endif()

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

# Products of this build: executables in bin/, libraries in lib/ (objects stay in CMakeFiles/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

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

# libasdf is still young and changes fast: say which version is built, with which local patches and
# which known bugs, so that updating the submodule isn't forgotten.
# The version is the one libasdf declares in its own CMakeLists.txt.
file(STRINGS ${CMAKE_SOURCE_DIR}/libasdf/CMakeLists.txt _libasdf_version_line REGEX "^set\\(PACKAGE_VERSION ")
string(REGEX REPLACE "^set\\(PACKAGE_VERSION +([^ )]+).*" "\\1" LIBASDF_VERSION "${_libasdf_version_line}")
if(NOT LIBASDF_VERSION)
    set(LIBASDF_VERSION "unknown")
endif()
set(LIBASDF_COMMIT "")
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/libasdf/.git)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} -C ${CMAKE_SOURCE_DIR}/libasdf log -1 "--format=%h of %cs"
        OUTPUT_VARIABLE LIBASDF_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
endif()
if(LIBASDF_COMMIT)
    set(LIBASDF_DESCRIPTION "${LIBASDF_VERSION} (commit ${LIBASDF_COMMIT})")
else()
    set(LIBASDF_DESCRIPTION "${LIBASDF_VERSION}")
endif()

# Bugs of this libasdf that have no patch, only a workaround; each has a test that fails once it is fixed.
set(LIBASDF_KNOWN_BUGS
    "asdf_set_*(file, \"a/b\") fails (ASDF_VALUE_ERR_OOM) when mapping a exists; test LibasdfKnownBugs.SetUnderAnExistingMapping"
)

message("${BoldMagenta}-- libasdf ${LIBASDF_DESCRIPTION}${ColourReset}  https://github.com/asdf-format/libasdf")

# libasdf is built from a patched copy, until the fixes are available upstream: every
# cmake/patches/libasdf-*.patch is applied, in name order. The hash in the command makes the copy
# start again whenever a patch is added, changed or removed.
find_program(PATCH_PROGRAM patch REQUIRED)
set(LIBASDF_COPY ${CMAKE_BINARY_DIR}/external/libasdf-src)
file(GLOB LIBASDF_PATCHES ${CMAKE_SOURCE_DIR}/cmake/patches/libasdf-*.patch)
list(SORT LIBASDF_PATCHES)
list(LENGTH LIBASDF_PATCHES LIBASDF_PATCH_COUNT)
set(LIBASDF_PATCH_HASH "")
set(LIBASDF_PATCH_NAMES "")
foreach(patch IN LISTS LIBASDF_PATCHES)
    file(SHA256 ${patch} hash)
    string(APPEND LIBASDF_PATCH_HASH ${hash})
    get_filename_component(name ${patch} NAME)
    list(APPEND LIBASDF_PATCH_NAMES ${name})
    # The first line of a patch says what it fixes
    file(STRINGS ${patch} _summary LIMIT_COUNT 1)
    message("${Blue}--   patch ${ColourReset} cmake/patches/${name}: ${_summary}")
endforeach()
foreach(bug IN LISTS LIBASDF_KNOWN_BUGS)
    message("${Yellow}--   known bug, not patched ${ColourReset} ${bug}")
endforeach()
string(SHA256 LIBASDF_PATCH_HASH "${LIBASDF_PATCH_HASH}")
# Configure again when a patch changes (its hash) or is added or removed (the directory), or when the
# submodule is updated (its version)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${LIBASDF_PATCHES} ${CMAKE_SOURCE_DIR}/cmake/patches
    ${CMAKE_SOURCE_DIR}/libasdf/CMakeLists.txt)
ExternalProject_Add(libasdf_external
    DEPENDS    libmd_external libfyaml_external
    PREFIX     ${CMAKE_BINARY_DIR}/external/libasdf
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -DSRC=${CMAKE_SOURCE_DIR}/libasdf -DDEST=${LIBASDF_COPY}
        -DPATCH_DIR=${CMAKE_SOURCE_DIR}/cmake/patches -DPATCH_PROGRAM=${PATCH_PROGRAM} -DPATCH_HASH=${LIBASDF_PATCH_HASH}
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
set(CMAKE_BUILD_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY};${STAGE_DIR}/lib;${DEEPSKY_RPATH}")
set(CMAKE_INSTALL_RPATH "${DEEPSKY_RPATH}")

#-------------------- libasdf2fits
message("${Red}----------- libasdf2fits -----------${ColourReset}")
set(ASDF2FITS_SOURCES
    src/ASDFmanager.cxx
    src/ASDFheader.cxx
    src/ASDFimage.cxx
    src/ASDFtable.cxx
    src/ASDFnode.cxx
    src/ASDFarray.cxx
    src/ASDFcube.cxx
    src/ASDFtransform.cxx
    src/ASDFgwcs.cxx
    src/ASDFtree.cxx
    src/ASDFwrite.cxx
    src/fits2asdf.cxx
)
set(ASDF2FITS_PUBLIC_HEADERS src/ASDFmanager.h src/fits2asdf.h src/ASDFnode.h src/ASDFarray.h src/ASDFcube.h src/ASDFtransform.h src/ASDFgwcs.h)

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
    # Written in the history entries of the files it saves
    target_compile_definitions(asdf2fits PRIVATE ASDF2FITS_VERSION="${PROJECT_VERSION}")
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
        ASDF_FIXTURES="${CMAKE_SOURCE_DIR}/libasdf/tests/fixtures"
        LIBASDF_VERSION="${LIBASDF_VERSION}"
        ASDF2FITS_TEST_DIR="${CMAKE_BINARY_DIR}/Testing/Temporary"
        ASDF2FITS_SAMPLES="${CMAKE_SOURCE_DIR}/test/data")
    # Real Roman and JWST files, too large for git: downloaded on request only (see cmake/FetchSamples.cmake)
    add_custom_target(fetch_samples
        COMMAND ${CMAKE_COMMAND} -DDEST=${CMAKE_SOURCE_DIR}/test/data -P ${CMAKE_SOURCE_DIR}/cmake/FetchSamples.cmake
        COMMENT "Downloading the Roman and JWST sample files"
        USES_TERMINAL)
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
message("${Green}-- libasdf          :${ColourReset} ${LIBASDF_DESCRIPTION}, ${LIBASDF_PATCH_COUNT} local patch(es)")
foreach(name IN LISTS LIBASDF_PATCH_NAMES)
    message("${Green}--                    ${ColourReset} patch ${name}")
endforeach()
list(LENGTH LIBASDF_KNOWN_BUGS _bug_count)
if(_bug_count GREATER 0)
    message("${Yellow}--                     ${_bug_count} known bug(s) without a patch: see above${ColourReset}")
endif()
message("${Yellow}--                     libasdf is in early development: check for a newer version, and drop the${ColourReset}")
message("${Yellow}--                     patches and workarounds it makes useless (a patch that no longer applies stops the build)${ColourReset}")
message("${Green}-- Install prefix   :${ColourReset} ${CMAKE_INSTALL_PREFIX}")
message("${BoldGreen}-------------------------------------------------${ColourReset}")
