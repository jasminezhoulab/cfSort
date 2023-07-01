import sys
import tensorflow as tf
import numpy as np
from library_utils import predict_on_test_set_with_trained_model


num_features = 10183
has_label = False
test_filename = sys.argv[1]
model_filename = sys.argv[2]
pred_filename = sys.argv[3]


feature_shape = (num_features, )
predict_on_test_set_with_trained_model(model_filename, test_filename, pred_filename, feature_shape, has_label)

print("python done")
