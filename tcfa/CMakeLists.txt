cmake_minimum_required(VERSION 2.8.9)

set (OE_DIR "~/openeye/toolkits")

set (OE_LIBS
    oedepict
    oebio
    oechem
    oesystem
    oeplatform
)

include_directories ( ${OE_DIR}/include )
link_directories ( ${OE_DIR}/lib )

set (LIBS "${OE_LIBS};${SYSLIBS};${SYSRASTERLIBS};${SYSZLIB}")

set (SRCS
    "tcfa tcfa.cpp"
)
set (CMAKE_BUILD_TYPE Debug)

include (${OE_DIR}/examples/ExamplesTemplate.cmake)
