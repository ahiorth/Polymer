SET ( MAINSRC
          solution.cpp
)


SET( executable   polymerX )

INCLUDE_DIRECTORIES   ("${PROJECT_SOURCE_DIR}/polymer")

ADD_EXECUTABLE(        ${executable} ${MAINSRC} )
if(OpenMP_CXX_FOUND)
    target_link_libraries(${executable} PUBLIC polymer OpenMP::OpenMP_CXX)
    message("CMAKE: COMPILING WITH OPENMP")
else()
    target_link_libraries(${executable} polymer) 
endif()




#install (TARGETS ${executable}
#         RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin)
