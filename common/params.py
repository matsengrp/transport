import os

CSV_OUTPUT = "csv"
DEFAULT_LAMBDA = 0.01
DEFAULT_NEIGHBOR_RADIUS = 48.5
DIST_MATRICES = "dist_matrices"
DMAX = 200
HMM_OUTPUT = "hmm_output"
JSON_OUTPUT = "json_output"
OUTPUT = "output"
TMP_OUTPUT = "tmp_output"
TRB_HUMAN_CDR3_HMM = "data/hmms/trb_human_cdr3.hmm"
TRB_HUMAN_CDR3_STO = "data/hmms/trb_human_cdr3.sto"
TRB_MOUSE_HMM = "data/hmms/TRB_mouse.hmm"
TRB_MOUSE_CDR3_HMM = "data/hmms/trb_mouse_cdr3.hmm"
TRB_MOUSE_CDR3_STO = "data/hmms/trb_mouse_cdr3.sto"

OUTPUT_DIRNAME = "output"

CSV_OUTPUT_DIRNAME = os.path.join(OUTPUT_DIRNAME, CSV_OUTPUT)
DIST_MATRICES_DIRNAME = os.path.join(OUTPUT_DIRNAME, "dist_matrices")
JSON_OUTPUT_DIRNAME = os.path.join(OUTPUT_DIRNAME, "json")
HMM_OUTPUT_DIRNAME = os.path.join(OUTPUT_DIRNAME, "hmm")
TMP_OUTPUT_DIRNAME = "tmp_output"

DIRECTORIES = {
    OUTPUT: OUTPUT_DIRNAME,
    TMP_OUTPUT: TMP_OUTPUT_DIRNAME,
    DIST_MATRICES: DIST_MATRICES_DIRNAME,
    JSON_OUTPUT: JSON_OUTPUT_DIRNAME,
    HMM_OUTPUT: HMM_OUTPUT_DIRNAME,
}

IEL_DATA_DIR = '/loc/no-backup/pbradley/share/pot_data/iels_tcrs_by_mouse/'

