cmake_minimum_required(VERSION 3.16)
project(hdfs-cpp VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(BUILD_TESTS "Build unit tests" ON)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Find required packages
find_package(Threads REQUIRED)

# Try to find gRPC first (it will bring in protobuf)
find_package(gRPC CONFIG QUIET)

if(gRPC_FOUND)
    message(STATUS "Found gRPC via CMake config")
    # gRPC brings its own protobuf, use it
    if(NOT TARGET protobuf::libprotobuf)
        find_package(Protobuf REQUIRED)
    endif()
else()
    # gRPC not found via config, find protobuf separately
    find_package(Protobuf REQUIRED)
    
    # Try pkg-config for gRPC
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(GRPC QUIET grpc++ grpc)
        if(GRPC_FOUND)
            set(gRPC_FOUND TRUE)
        endif()
    endif()
    
    if(NOT gRPC_FOUND)
        # Try Homebrew paths on macOS
        if(APPLE)
            set(HOMEBREW_PREFIX "/opt/homebrew" CACHE PATH "Homebrew install prefix")
            if(EXISTS "${HOMEBREW_PREFIX}/lib/libgrpc++.dylib")
                set(gRPC_FOUND TRUE)
                set(GRPC_INCLUDE_DIRS "${HOMEBREW_PREFIX}/include")
                set(GRPC_LIBRARIES 
                    "${HOMEBREW_PREFIX}/lib/libgrpc++.dylib"
                    "${HOMEBREW_PREFIX}/lib/libgrpc.dylib"
                    "${HOMEBREW_PREFIX}/lib/libgpr.dylib"
                    "${HOMEBREW_PREFIX}/lib/libaddress_sorting.dylib"
                )
                # Find grpc_cpp_plugin
                find_program(GRPC_CPP_PLUGIN grpc_cpp_plugin 
                    HINTS "${HOMEBREW_PREFIX}/bin" "/usr/local/bin")
            endif()
        endif()
    endif()
    
    if(NOT gRPC_FOUND)
        message(STATUS "gRPC not found, will use stub implementations")
        set(USE_GRPC_STUBS TRUE)
    endif()
endif()

# Find or download spdlog
find_package(spdlog QUIET)
if(NOT spdlog_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.12.0
    )
    FetchContent_MakeAvailable(spdlog)
endif()

# Find yaml-cpp for configuration (optional)
find_package(yaml-cpp QUIET)
if(yaml-cpp_FOUND)
    set(HAVE_YAML_CPP TRUE)
else()
    message(STATUS "yaml-cpp not found, using simple config parser")
    set(HAVE_YAML_CPP FALSE)
endif()

# Find or download abseil
find_package(absl QUIET)
if(NOT absl_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        abseil-cpp
        GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
        GIT_TAG 20230802.1
    )
    set(ABSL_PROPAGATE_CXX_STD ON)
    FetchContent_MakeAvailable(abseil-cpp)
endif()

# Proto files
set(PROTO_DIR ${CMAKE_CURRENT_SOURCE_DIR}/proto)
set(PROTO_OUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/proto)
file(MAKE_DIRECTORY ${PROTO_OUT_DIR})

# Find protoc executable
if(TARGET protobuf::protoc)
    set(PROTOC_EXECUTABLE $<TARGET_FILE:protobuf::protoc>)
elseif(Protobuf_PROTOC_EXECUTABLE)
    set(PROTOC_EXECUTABLE ${Protobuf_PROTOC_EXECUTABLE})
else()
    find_program(PROTOC_EXECUTABLE protoc REQUIRED)
endif()

# Function to generate protobuf files (gRPC optional)
function(generate_proto PROTO_FILE)
    get_filename_component(PROTO_NAME ${PROTO_FILE} NAME_WE)
    set(PROTO_SRC ${PROTO_OUT_DIR}/${PROTO_NAME}.pb.cc)
    set(PROTO_HDR ${PROTO_OUT_DIR}/${PROTO_NAME}.pb.h)
    
    if(TARGET gRPC::grpc_cpp_plugin)
        set(GRPC_SRC ${PROTO_OUT_DIR}/${PROTO_NAME}.grpc.pb.cc)
        set(GRPC_HDR ${PROTO_OUT_DIR}/${PROTO_NAME}.grpc.pb.h)
        add_custom_command(
            OUTPUT ${PROTO_SRC} ${PROTO_HDR} ${GRPC_SRC} ${GRPC_HDR}
            COMMAND ${PROTOC_EXECUTABLE}
                --cpp_out=${PROTO_OUT_DIR}
                --grpc_out=${PROTO_OUT_DIR}
                --plugin=protoc-gen-grpc=$<TARGET_FILE:gRPC::grpc_cpp_plugin>
                -I${PROTO_DIR}
                ${PROTO_FILE}
            DEPENDS ${PROTO_FILE}
            COMMENT "Generating protobuf/gRPC for ${PROTO_NAME}"
        )
        set(PROTO_SRCS ${PROTO_SRC} ${GRPC_SRC} PARENT_SCOPE)
    elseif(GRPC_CPP_PLUGIN)
        set(GRPC_SRC ${PROTO_OUT_DIR}/${PROTO_NAME}.grpc.pb.cc)
        set(GRPC_HDR ${PROTO_OUT_DIR}/${PROTO_NAME}.grpc.pb.h)
        add_custom_command(
            OUTPUT ${PROTO_SRC} ${PROTO_HDR} ${GRPC_SRC} ${GRPC_HDR}
            COMMAND ${PROTOC_EXECUTABLE}
                --cpp_out=${PROTO_OUT_DIR}
                --grpc_out=${PROTO_OUT_DIR}
                --plugin=protoc-gen-grpc=${GRPC_CPP_PLUGIN}
                -I${PROTO_DIR}
                ${PROTO_FILE}
            DEPENDS ${PROTO_FILE}
            COMMENT "Generating protobuf/gRPC for ${PROTO_NAME}"
        )
        set(PROTO_SRCS ${PROTO_SRC} ${GRPC_SRC} PARENT_SCOPE)
    else()
        # Protobuf only, no gRPC
        add_custom_command(
            OUTPUT ${PROTO_SRC} ${PROTO_HDR}
            COMMAND ${PROTOC_EXECUTABLE}
                --cpp_out=${PROTO_OUT_DIR}
                -I${PROTO_DIR}
                ${PROTO_FILE}
            DEPENDS ${PROTO_FILE}
            COMMENT "Generating protobuf for ${PROTO_NAME}"
        )
        set(PROTO_SRCS ${PROTO_SRC} PARENT_SCOPE)
    endif()
endfunction()

# Generate proto files
file(GLOB PROTO_FILES ${PROTO_DIR}/*.proto)
set(ALL_PROTO_SRCS "")
foreach(PROTO_FILE ${PROTO_FILES})
    generate_proto(${PROTO_FILE})
    list(APPEND ALL_PROTO_SRCS ${PROTO_SRCS})
endforeach()

# Proto library
add_library(hdfs_proto STATIC ${ALL_PROTO_SRCS})
target_include_directories(hdfs_proto PUBLIC ${PROTO_OUT_DIR} ${Protobuf_INCLUDE_DIRS})
if(TARGET gRPC::grpc++)
    target_link_libraries(hdfs_proto PUBLIC protobuf::libprotobuf gRPC::grpc++)
elseif(GRPC_LIBRARIES)
    target_include_directories(hdfs_proto PUBLIC ${GRPC_INCLUDE_DIRS})
    target_link_libraries(hdfs_proto PUBLIC protobuf::libprotobuf ${GRPC_LIBRARIES})
else()
    target_link_libraries(hdfs_proto PUBLIC protobuf::libprotobuf)
    target_compile_definitions(hdfs_proto PUBLIC HDFS_NO_GRPC)
endif()

# Common library
add_library(hdfs_common STATIC
    src/common/config.cpp
    src/common/logging.cpp
    src/common/utils/checksum.cpp
    src/common/utils/thread_pool.cpp
    src/common/utils/uuid.cpp
    src/common/rpc/rpc_server.cpp
    src/common/rpc/rpc_client.cpp
    src/common/ec/erasure_coding.cpp
)
target_include_directories(hdfs_common PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(hdfs_common PUBLIC
    hdfs_proto
    spdlog::spdlog
    absl::strings
    absl::time
    absl::random_random
    Threads::Threads
)
if(HAVE_YAML_CPP)
    target_link_libraries(hdfs_common PUBLIC yaml-cpp)
    target_compile_definitions(hdfs_common PUBLIC HAVE_YAML_CPP)
endif()

# NameNode library
add_library(hdfs_namenode STATIC
    src/namenode/namenode.cpp
    src/namenode/namespace.cpp
    src/namenode/block_manager.cpp
    src/namenode/datanode_manager.cpp
    src/namenode/edit_log.cpp
    src/namenode/fsimage.cpp
    src/namenode/snapshot.cpp
    src/namenode/quota.cpp
    src/namenode/ha/ha_context.cpp
    src/namenode/ha/standby_checkpoint.cpp
    src/namenode/ha/failover.cpp
    src/namenode/federation/block_pool.cpp
    src/namenode/federation/namespace_info.cpp
    src/namenode/federation/router.cpp
)
target_include_directories(hdfs_namenode PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(hdfs_namenode PUBLIC hdfs_common)

# DataNode library
add_library(hdfs_datanode STATIC
    src/datanode/datanode.cpp
    src/datanode/block_pool.cpp
    src/datanode/block_receiver.cpp
    src/datanode/block_sender.cpp
    src/datanode/heartbeat.cpp
    src/datanode/replication.cpp
)
target_include_directories(hdfs_datanode PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(hdfs_datanode PUBLIC hdfs_common)

# JournalNode library
add_library(hdfs_journalnode STATIC
    src/journalnode/journal.cpp
)
target_include_directories(hdfs_journalnode PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(hdfs_journalnode PUBLIC hdfs_common)

# Client library
add_library(hdfs_client STATIC
    src/client/hdfs_client.cpp
    src/client/input_stream.cpp
    src/client/output_stream.cpp
    src/client/filesystem.cpp
)
target_include_directories(hdfs_client PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(hdfs_client PUBLIC hdfs_common)

# NameNode executable
add_executable(namenode src/namenode/main.cpp)
target_link_libraries(namenode PRIVATE hdfs_namenode)

# DataNode executable
add_executable(datanode src/datanode/main.cpp)
target_link_libraries(datanode PRIVATE hdfs_datanode)

# JournalNode executable
add_executable(journalnode src/journalnode/main.cpp)
target_link_libraries(journalnode PRIVATE hdfs_journalnode hdfs_namenode)

# CLI executable
add_executable(hdfs-cli src/cli/main.cpp)
target_link_libraries(hdfs-cli PRIVATE hdfs_client)

# Router executable
add_executable(router src/router/main.cpp)
target_link_libraries(router PRIVATE hdfs_namenode)

# Install targets
install(TARGETS namenode datanode journalnode hdfs-cli router
    RUNTIME DESTINATION bin
)
install(TARGETS hdfs_client hdfs_common hdfs_proto
    ARCHIVE DESTINATION lib
)
install(DIRECTORY include/hdfs DESTINATION include)

# Tests
if(BUILD_TESTS)
    enable_testing()
    find_package(GTest QUIET)
    if(NOT GTest_FOUND)
        include(FetchContent)
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG v1.14.0
        )
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
    endif()
    
    # Unit tests for common utilities only (minimal dependencies)
    add_executable(hdfs_common_tests
        tests/unit/test_checksum.cpp
    )
    target_link_libraries(hdfs_common_tests PRIVATE
        hdfs_common
        GTest::gtest_main
    )
    include(GoogleTest)
    gtest_discover_tests(hdfs_common_tests)
endif()

