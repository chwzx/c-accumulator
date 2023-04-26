cmake_minimum_required(VERSION 3.24)
project(NTL_test)

set(CMAKE_CXX_STANDARD 17)

LINK_DIRECTORIES(/usr/local/lib)

link_libraries(/usr/local/lib/libntl.a)
#被链接的库相对路径
link_directories(/usr/local/lib /lib)
#被链接库的头文件相对路径
include_directories(/usr/local/include)


#add_executable(NTL_test main.cpp classgroup.cpp classgroup.h)
add_executable(NTL_test main.cpp util.cpp util.h hash/Hash.cpp hash/Hash.h proof/poe.cpp proof/poe.h group/classgroup.cpp proof/poke2.cpp proof/poke2.h accumulator.cpp accumulator.h vector_commitment.cpp vector_commitment.h)
#target_link_libraries(NTL_test pbc)s

target_link_libraries(NTL_test PRIVATE ${PROJECT_SOURCE_DIR}/lib/libblake3.so -lm -lgmp -lntl -lpthread -lpbc -lgmpxx -lb2s)