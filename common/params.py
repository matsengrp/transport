import os

CSV_OUTPUT = "csv"
DEFAULT_LAMBDA = 0.1
DEFAULT_NEIGHBOR_RADIUS = 50.5
DIST_MATRICES = "dist_matrices"
DMAX = 200
JSON_OUTPUT = "json_output"
OUTPUT = "output"
TMP_OUTPUT = "tmp_output"
TRB_MOUSE_HMM = "data/hmms/TRB_mouse.hmm"

OUTPUT_DIRNAME = "output"

CSV_OUTPUT_DIRNAME = os.path.join(OUTPUT_DIRNAME, CSV_OUTPUT)
DIST_MATRICES_DIRNAME = os.path.join(OUTPUT_DIRNAME, "dist_matrices")
JSON_OUTPUT_DIRNAME = os.path.join(OUTPUT_DIRNAME, "json")
TMP_OUTPUT_DIRNAME = "tmp_output"

DIRECTORIES = {
    OUTPUT: OUTPUT_DIRNAME,
    TMP_OUTPUT: TMP_OUTPUT_DIRNAME,
    DIST_MATRICES: DIST_MATRICES_DIRNAME,
    JSON_OUTPUT: JSON_OUTPUT_DIRNAME,
}

IEL_DATA_DIR = '/loc/no-backup/pbradley/share/pot_data/iels_tcrs_by_mouse/'

