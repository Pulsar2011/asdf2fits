# Copies the libasdf submodule to DEST and applies PATCH there, so that the submodule stays clean.
# Run with: cmake -DSRC=... -DDEST=... -DPATCH_FILE=... -DPATCH_PROGRAM=... -P PrepareLibasdf.cmake

# libasdf initialises its own submodules from its .git, which the copy doesn't have.
if(EXISTS "${SRC}/.git")
    execute_process(
        COMMAND git submodule update --init asdf-standard third_party/STC tests/munit
        WORKING_DIRECTORY "${SRC}"
        RESULT_VARIABLE _git_result
    )
    if(_git_result)
        message(FATAL_ERROR "Could not initialise the submodules of libasdf")
    endif()
endif()

file(REMOVE_RECURSE "${DEST}")
file(COPY "${SRC}/" DESTINATION "${DEST}" PATTERN ".git" EXCLUDE)

execute_process(
    COMMAND "${PATCH_PROGRAM}" -p1 --forward --batch -i "${PATCH_FILE}"
    WORKING_DIRECTORY "${DEST}"
    RESULT_VARIABLE _patch_result
)
if(_patch_result)
    message(FATAL_ERROR "${PATCH_FILE} doesn't apply to this libasdf. "
        "If the new version contains the fix, remove the patch from asdf2fits/CMakeLists.txt.")
endif()
