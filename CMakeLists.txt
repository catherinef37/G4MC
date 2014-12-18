#----------------------------------------------------------------------------

# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(g4mc)

enable_language (Fortran)


#ROOT path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
find_package(ROOT REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/HRSTransport/
    ${PROJECT_SOURCE_DIR}/HRSTransport/G2PTrans400016/
    ${PROJECT_SOURCE_DIR}/HRSTransport/G2PTrans484816/
    ${PROJECT_SOURCE_DIR}/HRSTransport/G2PTrans484816R00/
    ${PROJECT_SOURCE_DIR}/HRSTransport/GDHTransLargeX0/
    ${PROJECT_SOURCE_DIR}/HRSTransport/GDHTransSTD/
    ${PROJECT_SOURCE_DIR}/HRSTransport/hamcPREXTrans/
    ${PROJECT_SOURCE_DIR}/HRSTransport/HRSTransSTD/

    ${PROJECT_SOURCE_DIR}/XSModel/Bosted/
    ${PROJECT_SOURCE_DIR}/XSModel/Compton/
    ${PROJECT_SOURCE_DIR}/XSModel/Elas/
    ${PROJECT_SOURCE_DIR}/XSModel/QFS_N_EPC/
    ${PROJECT_SOURCE_DIR}/XSModel/Wiser
    )


#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})


#debugging symbols, complains
set(CMAKE_CXX_FLAGS "-g -Wall -pg")

#----------------------------------------------------------------------------
# Generate dictionaries as well (so we can write objects to ROOT files

include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

file(GLOB transsources ${PROJECT_SOURCE_DIR}/HRSTransport/*.cc 
                       ${PROJECT_SOURCE_DIR}/HRSTransport/G2PTrans400016/*cc  

                       ${PROJECT_SOURCE_DIR}/HRSTransport/G2PTrans484816/*.cc 
                       ${PROJECT_SOURCE_DIR}/HRSTransport/G2PTrans484816/*.f 

                       ${PROJECT_SOURCE_DIR}/HRSTransport/G2PTrans484816R00/*.cc

                       ${PROJECT_SOURCE_DIR}/HRSTransport/GDHTransLargeX0/*.cc 
                       ${PROJECT_SOURCE_DIR}/HRSTransport/GDHTransLargeX0/*.f 
                       ${PROJECT_SOURCE_DIR}/HRSTransport/GDHTransSTD/*.cc
                       ${PROJECT_SOURCE_DIR}/HRSTransport/GDHTransSTD/*.f

                       ${PROJECT_SOURCE_DIR}/HRSTransport/hamcPREXTrans/*.cc
                       ${PROJECT_SOURCE_DIR}/HRSTransport/hamcPREXTrans/*.f

                       ${PROJECT_SOURCE_DIR}/HRSTransport/HRSTransSTD/*.cc
                       ${PROJECT_SOURCE_DIR}/HRSTransport/HRSTransSTD/*.f
                       )

file(GLOB transheaders ${PROJECT_SOURCE_DIR}/HRSTransport/*.hh
                       ${PROJECT_SOURCE_DIR}/HRSTransport/G2PTrans400016/*.hh 

                       ${PROJECT_SOURCE_DIR}/HRSTransport/G2PTrans484816/*.hh 

                       ${PROJECT_SOURCE_DIR}/HRSTransport/G2PTrans484816R00/*.hh

                       ${PROJECT_SOURCE_DIR}/HRSTransport/GDHTransLargeX0/*.hh 

                       ${PROJECT_SOURCE_DIR}/HRSTransport/GDHTransSTD/*.hh

                       ${PROJECT_SOURCE_DIR}/HRSTransport/hamcPREXTrans/*.hh

                       ${PROJECT_SOURCE_DIR}/HRSTransport/HRSTransSTD/*.hh
    )

file(GLOB xssources ${PROJECT_SOURCE_DIR}/XSModel/*.cc 
    ${PROJECT_SOURCE_DIR}/XSModel/Bosted/*.cc
    ${PROJECT_SOURCE_DIR}/XSModel/Bosted/*.f

    ${PROJECT_SOURCE_DIR}/XSModel/Compton/*.cc
    ${PROJECT_SOURCE_DIR}/XSModel/Compton/*.f

    ${PROJECT_SOURCE_DIR}/XSModel/Elas/*.cc
    ${PROJECT_SOURCE_DIR}/XSModel/Elas/*.f

    ${PROJECT_SOURCE_DIR}/XSModel/QFS_N_EPC/*.cc
    ${PROJECT_SOURCE_DIR}/XSModel/QFS_N_EPC/*.f

    ${PROJECT_SOURCE_DIR}/XSModel/Wiser/*.cc
    ${PROJECT_SOURCE_DIR}/XSModel/Wiser/*.f
    )

file(GLOB xsheaders ${PROJECT_SOURCE_DIR}/XSModel/*.hh 
    ${PROJECT_SOURCE_DIR}/XSModel/Bosted/*.hh
    ${PROJECT_SOURCE_DIR}/XSModel/Compton/*.hh
    ${PROJECT_SOURCE_DIR}/XSModel/Elas/*.hh
    ${PROJECT_SOURCE_DIR}/XSModel/QFS_N_EPC/*.hh
    ${PROJECT_SOURCE_DIR}/XSModel/Wiser/*.hh
    )



#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(G4MC G4MC.cc ${sources} ${headers})
add_library(HRSTransport ${transsources} ${transheaders})
add_library(XSModel ${xssources} ${xsheaders})
target_link_libraries(G4MC ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} HRSTransport XSModel)

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(g4mc DEPENDS G4MC)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS G4MC DESTINATION bin)


