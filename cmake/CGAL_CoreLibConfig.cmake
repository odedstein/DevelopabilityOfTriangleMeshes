set(WITH_CGAL_Core "ON")

# The else condition of this code is never used in an installed
# version since it cannot happen there. Note also that for
# CMake<=2.8.11 (detected by the absence of CMP0024), the else()
# condition is never used.
if(NOT POLICY CMP0024 OR NOT CGAL_BUILDING_LIBS)
  if(NOT MSVC AND NOT CGAL_HEADER_ONLY)
    get_property(CGAL_Core_LIBRARY TARGET CGAL::CGAL_Core PROPERTY LOCATION)
  else()
    set(CGAL_Core_LIBRARY "")
  endif()
else()
  # We are currently in a CGAL Build and CGALExports.cmake has not
  # necessarily been created yet. Just alias the targets. Also don't
  # access the LOCATION property here to set lib_LIBRARY, since those
  # targets are not imported and this is disallowed by CMP0026. Just
  # set it to the target name.
  if(TARGET CGAL_Core AND NOT TARGET CGAL::CGAL_Core AND NOT CGAL_HEADER_ONLY)
    add_library(CGAL::CGAL_Core ALIAS CGAL_Core)
    set(CGAL_Core_LIBRARY CGAL::CGAL_Core)
  else()
    set(CGAL_Core_LIBRARY "")
  endif()
endif()


# 3RD_PARTY variables.
set(CGAL_Core_3RD_PARTY_INCLUDE_DIRS   "")
set(CGAL_Core_3RD_PARTY_DEFINITIONS    "")
set(CGAL_Core_3RD_PARTY_LIBRARIES_DIRS "")
set(CGAL_Core_3RD_PARTY_LIBRARIES      "")
