# PolySolve (https://github.com/polyfem/polysolve)
# License: MIT

if(TARGET polysolve)
    return()
endif()

message(STATUS "Third-party: creating target 'polysolve'")

include(CPM)
CPMAddPackage("gh:maxpaik16/polysolve#9cdee7185e87a77efb4635e441973c63f60aebc0")
#target_link_libraries(polyfem "/Users/mpaik/Documents/Documents-CIMS-PHD-AP18/NYU/Daniele/Development/polysolve/build/libpolysolve.a")