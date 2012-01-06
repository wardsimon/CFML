# Add the CFML directory.
add_subdirectory(CFML)

# In case of a GUI build, add the CFGL directory.
if(${GUI})
    add_subdirectory(CFGL)
endif()
