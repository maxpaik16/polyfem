# PolySolve (https://github.com/polyfem/polysolve)
# License: MIT

if(TARGET polysolve)
    return()
endif()

message(STATUS "Third-party: creating target 'polysolve'")

include(CPM)
CPMAddPackage("gh:maxpaik16/polysolve#6bbb31ac556f34c393ecfca638fa7b5ce6ad4ae9")
