cmake_minimum_required(VERSION 3.14)
set(CMAKE_C_STANDARD 11)

set(ENVSDK $ENV{PLAYDATE_SDK_PATH})

if (NOT ${ENVSDK} STREQUAL "")
	# Convert path from Windows
	file(TO_CMAKE_PATH ${ENVSDK} SDK)
else()
	execute_process(
			COMMAND bash -c "egrep '^\\s*SDKRoot' $HOME/.Playdate/config"
			COMMAND head -n 1
			COMMAND cut -c9-
			OUTPUT_VARIABLE SDK
			OUTPUT_STRIP_TRAILING_WHITESPACE
	)
endif()

if (NOT EXISTS ${SDK})
	message(FATAL_ERROR "SDK Path not found; set ENV value PLAYDATE_SDK_PATH")
	return()
endif()

set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
set(CMAKE_XCODE_GENERATE_SCHEME TRUE)

# Game Name Customization
set(PLAYDATE_GAME_NAME isometric)
set(PLAYDATE_GAME_DEVICE isometric_DEVICE)

project(${PLAYDATE_GAME_NAME} C ASM)

# Add the display folder to the include path
include_directories(${CMAKE_SOURCE_DIR}/display)

if (TOOLCHAIN STREQUAL "armgcc")
	add_executable(${PLAYDATE_GAME_DEVICE} playdate/main.c display/display.c ${IMAGES})
else()
	add_library(${PLAYDATE_GAME_NAME} SHARED playdate/main.c display/display.c ${IMAGES})
endif()

include(${SDK}/C_API/buildsupport/playdate_game.cmake)

