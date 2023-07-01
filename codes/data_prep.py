import sys
import tensorflow as tf
import numpy as np
from library_utils import detect_normalization_params, import_labels, import_data_from_text_export_to_tfrecord, import_data_from_text_export_to_tfrecord_testonly


label_filename = sys.argv[1]
filename = sys.argv[2]
test_filename = sys.argv[3]
train_output_filename = sys.argv[4]
valid_output_filename = sys.argv[5]
test_output_filename = sys.argv[6]

label_dict, train_test_split = import_labels(label_filename)
norm_params = detect_normalization_params(filename, train_test_split)
import_data_from_text_export_to_tfrecord(filename, train_output_filename, valid_output_filename, norm_params, label_dict, train_test_split, "log_min_max")
import_data_from_text_export_to_tfrecord_testonly(test_filename, test_output_filename, norm_params, label_dict, "log_min_max")

