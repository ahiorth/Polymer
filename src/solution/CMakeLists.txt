SET ( MAINSRC
          solution.cpp
)


SET( executable   polymerX )

INCLUDE_DIRECTORIES   ("${PROJECT_SOURCE_DIR}/polymer")

ADD_EXECUTABLE(        ${executable} ${MAINSRC} )
target_link_libraries(${executable} polymer)


