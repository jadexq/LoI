cmake_minimum_required(VERSION 3.3.1)
project(code)

set(CMAKE_CXX_STANDARD 14)

include_directories(/users/student/include/eigen3 .)

add_executable(exe
    main.cpp Constants.cpp Generate.cpp Initial.cpp Mcmc.cpp
    rnorm.cpp mvnorm.cpp rgamma.cpp normalPdf.cpp normalPdfVec.cpp runif.cpp
    runifVec.cpp rbeta.cpp
)
