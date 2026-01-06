# PolySolve (https://github.com/polyfem/polysolve)
# License: MIT

if(TARGET polysolve)
    return()
endif()

message(STATUS "Third-party: creating target 'polysolve'")

include(CPM)
CPMAddPackage("gh:maxpaik16/polysolve#f006c435fe7c13fe8ba1094b732ad5277d74f5bf")
