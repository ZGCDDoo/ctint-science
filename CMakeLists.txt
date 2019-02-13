#Inspired by https://github.com/dmonopoly/gtest-cmake-example

cmake_minimum_required(VERSION 2.8)

SET(CMAKE_CXX_COMPILER g++)
option(SNAP OFF)


set(PROJECT_NAME CT-INT)
project("${PROJECT_NAME}")

set(CMAKE_CXX_FLAGS " -Wall -Wextra -std=c++14 -O3 -Ofast -fomit-frame-pointer")


find_package(LAPACK REQUIRED)
find_package(Boost REQUIRED filesystem system)


set(LIBRARIES_EXEC ${LAPACK_LIBRARIES} ${Boost_LIBRARIES} armadillo)




# Executables

################################
set(EXECECUTABLES 
    ctint-science
)

foreach(executable ${EXECECUTABLES})
    list(APPEND ${executable}_src ${PROJECT_SOURCE_DIR}/src/${executable}.cpp)
    add_executable(${executable}
                ${${executable}_src}
                ${PROJECT_SOURCE_DIR}/src/Utilities.cpp
                  )
endforeach()

foreach(executable ${EXECECUTABLES})
    target_include_directories(${executable} PRIVATE ${PROJECT_SOURCE_DIR}/include)
    target_link_libraries(${executable}  PRIVATE ${LIBRARIES_EXEC})
    if(SNAP)
        install(TARGETS ${executable} DESTINATION /bin)
    else()
        install(TARGETS ${executable} DESTINATION $ENV{HOME}/bin)
    endif()
endforeach()


#--------------------end Executables ----------------------------------------------------------



