import sys, os
import numpy as np
from datetime import datetime
import logging
import gc
import tensorflow as tf
import collections


train_file = sys.argv[1]
valid_file = sys.argv[2]
samples_list_file = sys.argv[3]
num_nodes = [int(i) for i in sys.argv[4].split(',')]
rate_dropout = [float(i) for i in sys.argv[5].split(',')]
learning_rate = float(sys.argv[6])
batch_size = int(sys.argv[7])
loss = sys.argv[8]
optimizer = sys.argv[9]
output_prefix = sys.argv[10]
epoch = int(sys.argv[11])
global_seed = int(sys.argv[12])

tf.random.set_seed(global_seed)

from library_utils import NN_multiple_response_batch_variedLR_startMiddle

num_features = 10183
TOTAL_TRAINING_SAMPLES = 259484
step = round(float(TOTAL_TRAINING_SAMPLES)/batch_size) -2
feature_shape = (num_features,)

model_parameter = {"num_nodes": num_nodes, "rate_dropout": rate_dropout, "step": step, "learning_rate": learning_rate, "batch_size": batch_size, "loss": loss, "optimizer": optimizer, "epoch": epoch, "epoch_start" : 0, "existing_model": "None", "new_model" : True, "initializer" : "default", }

NN_multiple_response_batch_variedLR_startMiddle(model_parameter, train_file, valid_file, output_prefix, feature_shape, True)




