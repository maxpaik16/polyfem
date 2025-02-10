# PolySolve (https://github.com/polyfem/polysolve)
# License: MIT

if(TARGET polysolve)
    return()
endif()

message(STATUS "Third-party: creating target 'polysolve'")

include(CPM)
CPMAddPackage("gh:maxpaik16/polysolve#af0a111563709b12e301a77fec123512b4cd0f1e")
#target_link_libraries(polyfem "/Users/mpaik/Documents/Documents-CIMS-PHD-AP18/NYU/Daniele/Development/polysolve/build/libpolysolve.a")