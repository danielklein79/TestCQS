ordinary_structures <- c(
  "S", "D", "CS", "CT"
)

doubly_structures <- c(
  "BS_D", "BS_CS", "BS_CT",
  "BD_S", "BD_CS", "BD_CT",
  "BCS_S", "BCS_D", "BCS_CS", "BCS_CT",
  "BCT_S", "BCT_D", "BCT_CS", "BCT_CT",
  "DCS"
)

UBS_structure <- "UBS_CS"

structure_labels <- c(
  S = "Sphericity",
  D = "Diagonality",
  CS = "Compound symmetry",
  CT = "Circular Toeplitz",
  BS_D = "Block spherical with diagonal blocks",
  BS_CS = "Block spherical with compound symmetry blocks",
  BS_CT = "Block spherical with circular Toeplitz blocks",
  BD_S = "Block diagonal with spherical blocks",
  BD_CS = "Block diagonal with compound symmetry blocks",
  BD_CT = "Block diagonal with circular Toeplitz blocks",
  BCS_S = "Block compound symmetry with spherical blocks",
  BCS_D = "Block compound symmetry with diagonal blocks",
  BCS_CS = "Block compound symmetry with compound symmetry blocks",
  BCS_CT = "Block compound symmetry with circular Toeplitz blocks",
  BCT_S = "Block circular Toeplitz with spherical blocks",
  BCT_D = "Block circular Toeplitz with diagonal blocks",
  BCT_CS = "Block circular Toeplitz with compound symmetry blocks",
  BCT_CT = "Block circular Toeplitz with circular Toeplitz blocks",
  UBS_CS = "Unifrom block structure with compound symmetry block eigenvalue",
  DCS = "Double complete symmetry"
)
