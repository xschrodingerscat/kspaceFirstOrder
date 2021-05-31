
project(KSpaceSolver)


# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(HDF5_DIR ${EBROOTHDF5})
set(FFT_DIR ${EBROOTFFTW})
set(ZLIB_DIR ${EBROOTZLIB})
set(SZIP_DIR ${EBROOTSZIP})


# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)


# Now simply link against gtest or gtest_main as needed. Eg
include_directories(
        lib
        /usr/local/include
        ${gtest_SOURCE_DIR}/include
)
link_directories(
        /usr/local/lib
)

set(APP_SOURCE
        lib/Containers/MatrixContainer.cpp
        lib/Containers/OutputStreamContainer.cpp
        lib/Hdf5/Hdf5File.cpp
        lib/Hdf5/Hdf5FileHeader.cpp
        lib/KSpaceSolver/KSpaceFirstOrderSolver.cpp
        lib/Logger/Logger.cpp
        lib/MatrixClasses/BaseFloatMatrix.cpp
        lib/MatrixClasses/BaseIndexMatrix.cpp
        lib/MatrixClasses/ComplexMatrix.cpp
        lib/MatrixClasses/FftwComplexMatrix.cpp
        lib/MatrixClasses/FftwRealMatrix.cpp
        lib/MatrixClasses/IndexMatrix.cpp
        lib/MatrixClasses/RealMatrix.cpp
        lib/OutputStreams/BaseOutputStream.cpp
        lib/OutputStreams/IndexOutputStream.cpp
        lib/OutputStreams/CuboidOutputStream.cpp
        lib/OutputStreams/WholeDomainOutputStream.cpp
        lib/Parameters/CommandLineParameters.cpp
        lib/KSpace/KConfig.cpp
        lib/Parameters/Parameters.cpp
        lib/KSpace/KSpaceSolverFluid.cpp
        lib/KSpace/KSpaceSolverElastic.cpp
        lib/KSpace/KMatrixUnitTest.cpp
        src/KSpaceUnitTest.cpp
        lib/KSpace/KVersion.h lib/KSpace/KOutput.cpp lib/KSpace/KOutput.h lib/KSpace/KIndexOutputStream.cpp lib/KSpace/KIndexOutputStream.h)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -ggdb -ffast-math -fassociative-math")

#set(CMAKE_BUILD_TYPE "Debug")
#set(CMAKE_BUILD_TYPE "Release")
#set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
#set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3")

add_executable(${PROJECT_NAME} ${APP_SOURCE})
target_link_libraries(${PROJECT_NAME}
        gtest_main
        hdf5
        hdf5_hl
        fftw3
        )

#add_test(NAME path_finder_test COMMAND path_finder)

add_custom_command(TARGET ${PROJECT_NAME}
        POST_BUILD
        COMMAND echo "copy compile_commands.json to for YCM."
        COMMAND ${CMAKE_COMMAND} -E copy compile_commands.json ../
        )



