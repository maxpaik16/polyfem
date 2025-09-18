# IPC Toolkit (https://github.com/ipc-sim/ipc-toolkit)
# License: MIT

if(TARGET ipc::toolkit)
    return()
endif()

message(STATUS "Third-party: creating target 'ipc::toolkit'")

set(IPC_TOOLKIT_WITH_CUDA OFF CACHE INTERNAL "" FORCE)

include(CPM)
CPMAddPackage("gh:geometryprocessing/GCP-toolkit#7185ecd8b7b46aadae78db688c6029a82a21dcca")
