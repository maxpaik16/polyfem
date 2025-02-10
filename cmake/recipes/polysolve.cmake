# PolySolve (https://github.com/polyfem/polysolve)
# License: MIT

if(TARGET polysolve)
    return()
endif()

message(STATUS "Third-party: creating target 'polysolve'")

include(CPM)
CPMAddPackage("gh:maxpaik16/polysolve#b174162505b340dab8c32537c5211418de6b01bd")
#target_link_libraries(polyfem "/Users/mpaik/Documents/Documents-CIMS-PHD-AP18/NYU/Daniele/Development/polysolve/build/libpolysolve.a")