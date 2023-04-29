import sys
import tensorflow as tf
import numpy as np
from library_utils import predict_on_test_set_with_trained_model


#DEBUG = True
DEBUG = False
if DEBUG:
	test_filename = "/u/project/xjzhou/shuoli/tissue_deconvolution/supervised_deconvolution/alpha_value_data/alpha_matrix/test.loose.hypo.top100.tfrecords"
	model_filename = "kl_Adam_64_loose.hypo.100_1024,512,128,64_0.3,0.2,0.2,0.1_0.001.model.h5"
	pred_filename = "kl_Adam_64_loose.hypo.100_1024,512,128,64_0.3,0.2,0.2,0.1_0.001.y_pred.txt"
	num_features = 19209
	has_label = True
else:
	test_filename = sys.argv[1]
	model_filename = sys.argv[2]
	pred_filename = sys.argv[3]
	num_features = int(sys.argv[4])
	has_label = (sys.argv[5] == "True")


feature_shape = (num_features, )
predict_on_test_set_with_trained_model(model_filename, test_filename, pred_filename, feature_shape, has_label)

print("python done")
