# PolySolve (https://github.com/polyfem/polysolve)
# License: MIT

if(TARGET polysolve)
    return()
endif()

message(STATUS "Third-party: creating target 'polysolve'")

include(CPM)
CPMAddPackage("gh:maxpaik16/polysolve#11ff15746768f9f0d650c96d76ed29dfc42487d5")
#target_link_libraries(polyfem "/Users/mpaik/Documents/Documents-CIMS-PHD-AP18/NYU/Daniele/Development/polysolve/build/libpolysolve.a")