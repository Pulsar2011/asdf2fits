# Downloads the real ASDF files used to test the GWCS support into test/data (ignored by git: *.asdf, *.fits).
#
#   cmake -P cmake/FetchSamples.cmake                  # all files, ~450 MB
#   cmake -DSAMPLES=small -P cmake/FetchSamples.cmake  # only the GWCS objects, ~35 kB
#   cmake --build build --target fetch_samples         # same as the first, from a configured build
#
# A file already there with the right SHA-256 isn't downloaded again. Tests that need a missing file are skipped.
#
# Samples:
#  - gwcs_nircamwcs.asdf, gwcs_miriwcs.asdf: JWST NIRCam and MIRI imaging GWCS objects alone (no pixels), from the
#    test data of the gwcs package, pinned to one commit.
#  - roman_wfi01_f106_cal.asdf: a simulated Roman WFI L2 image (romancal 0.18), 4088x4088, with its GWCS in
#    meta/wcs; meta/wcsinfo/s_region gives the sky position of the corners of the bounding box.
#  - jwst_nrcb1_cal.fits: a JWST NIRCam L2 image (SMACS 0723, program 2736). The GWCS is in the ASDF extension of
#    the FITS file; the SCI header has a TAN-SIP approximation of it (A_ORDER 5).

if(NOT DEFINED SAMPLES)
    set(SAMPLES all)
endif()
if(NOT SAMPLES MATCHES "^(all|small)$")
    message(FATAL_ERROR "SAMPLES must be all or small, not '${SAMPLES}'")
endif()
if(NOT DEFINED DEST)
    get_filename_component(DEST "${CMAKE_CURRENT_LIST_DIR}/../test/data" ABSOLUTE)
endif()
file(MAKE_DIRECTORY "${DEST}")

set(GWCS_COMMIT 78d104ed54879e6f594c5b20881910413df99f8d)
set(GWCS_DATA https://raw.githubusercontent.com/spacetelescope/gwcs/${GWCS_COMMIT}/gwcs/tests/data)

# name | size class | SHA-256 | URL
set(sample_list
    "gwcs_nircamwcs.asdf|small|2b79d06269bdbfa14e4656dbc4628e407a12414efb4a4335e04a5a8b4bffecda|${GWCS_DATA}/nircamwcs.asdf"
    "gwcs_miriwcs.asdf|small|f3f4e6a6c02fd6a73f5c88f7b3e4c3f821047cac47656194e227a1c486bb2c07|${GWCS_DATA}/miriwcs.asdf"
    "roman_wfi01_f106_cal.asdf|large|920c202e912003173c8105a7c9f8f8654ff0f145b9753968a7d54ff55412c824|https://stpubdata.s3.amazonaws.com/roman/nexus/soc_simulations/tutorial_data/nexus_early_access/nexus_early_access/r0003201001001001004_0001_wfi01_f106_cal.asdf"
    "jwst_nrcb1_cal.fits|large|6e69c8d7cc745a6bdd432d60c69d2a09391db043e305435fa4c9a6f610d9762d|https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/jw02736001001_02105_00001_nrcb1_cal.fits"
)

set(failed "")
foreach(entry IN LISTS sample_list)
    string(REPLACE "|" ";" fields "${entry}")
    list(GET fields 0 name)
    list(GET fields 1 class)
    list(GET fields 2 sha)
    list(GET fields 3 url)

    if(SAMPLES STREQUAL "small" AND NOT class STREQUAL "small")
        continue()
    endif()

    set(path "${DEST}/${name}")
    if(EXISTS "${path}")
        file(SHA256 "${path}" have)
        if(have STREQUAL sha)
            message(STATUS "${name}: already there")
            continue()
        endif()
        message(STATUS "${name}: wrong checksum, downloading again")
    endif()

    message(STATUS "${name}: downloading")
    # To a temporary name, so that an interrupted download never leaves a file that looks complete
    file(DOWNLOAD "${url}" "${path}.part" STATUS status TLS_VERIFY ON SHOW_PROGRESS)
    list(GET status 0 code)
    if(code EQUAL 0)
        file(SHA256 "${path}.part" have)
        if(have STREQUAL sha)
            file(RENAME "${path}.part" "${path}")
            continue()
        endif()
        set(status "checksum ${have} instead of ${sha}")
    endif()
    file(REMOVE "${path}.part")
    message(WARNING "${name}: ${status}")
    list(APPEND failed ${name})
endforeach()

if(failed)
    message(FATAL_ERROR "Not downloaded: ${failed}")
endif()
message(STATUS "Samples in ${DEST}")
