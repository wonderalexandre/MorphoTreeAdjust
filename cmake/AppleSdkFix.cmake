if(APPLE AND CMAKE_OSX_SYSROOT AND NOT EXISTS "${CMAKE_OSX_SYSROOT}")
  set(MTA_STALE_MACOS_SDK "${CMAKE_OSX_SYSROOT}")
  execute_process(
    COMMAND xcrun --sdk macosx --show-sdk-path
    OUTPUT_VARIABLE MTA_ACTIVE_MACOS_SDK
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )

  if(EXISTS "${MTA_ACTIVE_MACOS_SDK}")
    message(STATUS "Resetting stale CMAKE_OSX_SYSROOT from ${MTA_STALE_MACOS_SDK} to ${MTA_ACTIVE_MACOS_SDK}")
    set(CMAKE_OSX_SYSROOT "${MTA_ACTIVE_MACOS_SDK}" CACHE PATH "The product will be built against the headers and libraries located inside the indicated SDK." FORCE)
  else()
    message(STATUS "Clearing stale CMAKE_OSX_SYSROOT ${MTA_STALE_MACOS_SDK}")
    unset(CMAKE_OSX_SYSROOT CACHE)
  endif()

  unset(MTA_ACTIVE_MACOS_SDK)
  unset(MTA_STALE_MACOS_SDK)
endif()
