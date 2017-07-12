cmake_minimum_required(VERSION 2.8.9)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

set (OE_DIR "~/openeye/toolkits")

set (OE_LIBS
    oespicoli
    oezap
    oegrid
    oefizzchem
    oedepict
    oebio
    oechem
    oesystem
    oeplatform
)

include_directories ( ${OE_DIR}/include )
link_directories ( ${OE_DIR}/lib )

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/)
add_executable(text2hex ${OE_DIR}/examples/text2hex.cpp)

set (LIBS "${OMPLIBS};${OE_LIBS};${SYSRASTERLIBS};${SYSLIBS};${SYSZLIB}")

set (SRCS
    "comfa comfa.cpp"
)
add_definitions("${OMPCFLAGS}")

include (${OE_DIR}/examples/ExamplesTemplate.cmake)
