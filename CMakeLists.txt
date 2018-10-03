#Inspired by https://github.com/dmonopoly/gtest-cmake-example
#           https://github.com/misaka-10032/boost-mpi-examples

cmake_minimum_required(VERSION 2.8)

SET(CMAKE_CXX_COMPILER g++)


set(PROJECT_NAME CT-INT)
project("${PROJECT_NAME}")

set(CMAKE_CXX_FLAGS " -Wall -Wextra -std=c++14 -O3 -Ofast -fomit-frame-pointer")

set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/src/Includes)

find_package(LAPACK REQUIRED)
find_package(Boost REQUIRED mpi serialization filesystem system)


set(LIBRARIES_EXEC ${LAPACK_LIBRARIES} ${Boost_LIBRARIES} armadillo)
set(LIBRARIES_TEST ${LIBRARIES_EXEC} gtest gtest_main)




# Executables

################################
set(EXECECUTABLES 
    dmft
)

foreach(executable ${EXECECUTABLES})
    list(APPEND ${executable}_src ${PROJECT_SOURCE_DIR}/src/${executable}.cpp)
    add_executable(${executable}
                ${${executable}_src}
                  )
endforeach()

foreach(executable ${EXECECUTABLES})
    target_link_libraries(${executable}  ${LIBRARIES_EXEC})
    install(TARGETS ${executable} DESTINATION $ENV{HOME}/bin)
endforeach()


#--------------------end Executables ----------------------------------------------------------



