import sys, os
import numpy as np
from datetime import datetime
import logging
import gc
import tensorflow as tf
import collections


train_file = sys.argv[1]
test_file = sys.argv[2]
samples_list_file = sys.argv[3]
num_nodes = [int(i) for i in sys.argv[4].split(',')]
rate_dropout = [float(i) for i in sys.argv[5].split(',')]
step = int(sys.argv[6]) + 1
learning_rate = float(sys.argv[7])
batch_size = int(sys.argv[8])
loss = sys.argv[9]
optimizer = sys.argv[10]
output_prefix = sys.argv[11]
epoch = int(sys.argv[12])
num_features = int(sys.argv[13])
fixed_learning_rate = (sys.argv[14] == "True")
epoch_start = int(sys.argv[15])
model_file = sys.argv[16]
new_model = (sys.argv[17] == "True")
is_cpu = (sys.argv[18] == "True")
initializer = sys.argv[19]
global_seed = int(sys.argv[20])
train_size = int(sys.argv[21])

tf.random.set_seed(global_seed)

from library_utils import NN_multiple_response_batch_variedLR_startMiddle


TOTAL_TRAINING_SAMPLES = train_size
step = round(float(TOTAL_TRAINING_SAMPLES)/batch_size) -2
feature_shape = (num_features,)

model_parameter = {"num_nodes": num_nodes, "rate_dropout": rate_dropout, "step": step, "learning_rate": learning_rate, "batch_size": batch_size, "loss": loss, "optimizer": optimizer, "epoch": epoch, "epoch_start" : epoch_start, "existing_model": model_file, "new_model" : new_model, "initializer" : initializer, }

NN_multiple_response_batch_variedLR_startMiddle(model_parameter, train_file, test_file, output_prefix, feature_shape, fixed_learning_rate)




