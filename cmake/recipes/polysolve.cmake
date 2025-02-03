# PolySolve (https://github.com/polyfem/polysolve)
# License: MIT

if(TARGET polysolve)
    return()
endif()

message(STATUS "Third-party: creating target 'polysolve'")

include(CPM)
CPMAddPackage("gh:maxpaik16/polysolve#40bc4bb0bde08128b908d40a20626101b229382b")
#target_link_libraries(polyfem "/Users/mpaik/Documents/Documents-CIMS-PHD-AP18/NYU/Daniele/Development/polysolve/build/libpolysolve.a")