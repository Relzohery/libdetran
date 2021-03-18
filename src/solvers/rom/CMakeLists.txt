# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(ROM_SRC
  ${SRC_DIR}/ROMBasis.cc
  ${SRC_DIR}/OperatorProjection.cc
  ${SRC_DIR}/ROMSolver.cc
  ${SRC_DIR}/KineticMatrices.cc
  ${SRC_DIR}/TransientSolver.cc
  ${SRC_DIR}/DEIM.cc
  ${SRC_DIR}/offline_stage.cc
  ${SRC_DIR}/ParametricRom.cc
  PARENT_SCOPE
)
