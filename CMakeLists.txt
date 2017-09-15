cmake_minimum_required (VERSION 2.6)
project (phonon)

set(VERSION "0.1.0a4")
set(AUTHOR "Jiao Lin")
set(AUTHOR_EMAILE "jiao.lin@gmail.com")
set(VENDOR "MCViNE")
set(SUMMARY "mcvine.phonon python subpackage")
set(LICENSE "BSD")
set(RPM_APP_GROUP "Development/Tools")
set(GITHUB_URL "https://github.com/mcvine/phonon")
set (HOMEPAGE ${GITHUB_URL})

# dependencies
set (CONDA_PKGING_BUILD_DEPS "mcvine")
set (CONDA_PKGING_RUN_DEPS "mcvine")

# variables
set (PARENT_PATH "mcvine")
set (PARENT_PY_PKG "mcvine")
set (PARENT_GITHUB_URL "https://github.com/mcvine/mcvine")
set (PARENT_VERSION "1.3")
set (PARENT_RELEASE "1")
set (PKG_NAME "${PARENT_PY_PKG}.${PROJECT_NAME}")


# needs python
find_package ( PythonLibs REQUIRED )
find_package ( PythonInterp REQUIRED )
# python installation path
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_PY_PATH "${INSTALL_LIB_DIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages")

#
enable_testing()


# copy all python source over to build directory
file(COPY ${PROJECT_NAME} DESTINATION "modules/mcvine/")

# install
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/
  DESTINATION ${INSTALL_PY_PATH}
  )
