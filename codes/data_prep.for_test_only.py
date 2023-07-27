import sys
import tensorflow as tf
import numpy as np
from library_utils import detect_normalization_params, import_labels, import_data_from_text_export_to_tfrecord, import_data_from_text_export_to_tfrecord_testonly


test_filename = sys.argv[1]
test_output_filename = sys.argv[2]
norm_params_filename = sys.argv[3]

norm_params = np.load(norm_params_filename, allow_pickle = True)[()] 
import_data_from_text_export_to_tfrecord_testonly(test_filename, test_output_filename, norm_params, "log_min_max")

