add_subdirectory(zupply)

if (NOT shlibpp_FOUND)
    add_subdirectory(sharedlibpp)
endif()
