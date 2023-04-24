######################################################################
#
# Compilers has to set following variables:

# - AGRESSIVE_FINAL_L         :: list of objects that are agresivelly
#                                optimized
# - PASIVE_FINAL_L            :: list of objects that are
#                                non-agressivelly optimized
# - AGRESSIVE_FLAG_S          :: Agresive flag
# - PASIVE_FLAG_S             :: Non-agresive flag
#
# For each compilers optimization flags for EXT_OPT=ON/OFF
# has to be set.
#

add_subdirectory(nvhpc)
add_subdirectory(GNU)
add_subdirectory(intel)
add_subdirectory(oneapi)
add_subdirectory(ibm)
add_subdirectory(appleclang)

#
#####################################################################
