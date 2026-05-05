# IPC Toolkit (https://github.com/ipc-sim/ipc-toolkit)
# License: MIT

if(TARGET ipc::toolkit)
    return()
endif()

message(STATUS "Third-party: creating target 'ipc::toolkit'")

set(IPC_TOOLKIT_WITH_CUDA OFF CACHE INTERNAL "" FORCE)

include(CPM)
CPMAddPackage("gh:maxpaik16/ipc-toolkit#68a87c116e448f8d22537ab5995d91affad78b25")
