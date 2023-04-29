import sys, gzip, re
import time
from datetime import datetime
import numpy as np
from scipy import sparse
from collections import defaultdict
from scipy.special import expit, logit
import logging
import gc
import tensorflow as tf
import random
import collections



####################################################################################
####################################################################################
## Data normalization and loading
####################################################################################
####################################################################################


NUM_TISSUES = 29
LABEL_SHAPE = (NUM_TISSUES,)


def detect_normalization_params(filename, train_test_split, norm_option = "log_min_max"):
	params = {}
	with open(filename, 'r') as f:
		for line in f:
			items = line.rstrip().split('\t')
			sample_id = items[0]
			num_features = len(items) - 1
			sample_train_test_split = train_test_split.get(sample_id)
			if sample_train_test_split == '1':
				features = np.array([float(items[i]) if items[i] != "NA" else 0.0 for i in range(1, len(items))])
				features = np.log2(features + 1)
				if len(params) == 0:
					params["min"] = features
					params["max"] = features
				else:
					params["min"] = np.minimum(features, params["min"])
					params["max"] = np.maximum(features, params["max"])
	return params


def import_labels(filename):
	labels = {}
	train_test_split = {}
	with open(filename, 'r') as f:
		for line in f:
			items = line.rstrip().split('\t')
			sample_id = items[0]
			train_test_split[sample_id] = items[1]
			labels[sample_id] = np.array([float(items[i]) for i in range(2, 2 + NUM_TISSUES)])
	return labels, train_test_split


def normalize_features(features, norm_params, norm_option = "log_min_max"):
	if norm_option == "log_min_max":
		diff = norm_params["max"] - norm_params["min"]
		diff = np.where(diff == 0.0, 1.0, diff)
		norm_features = (np.log2(features + 1) - norm_params["min"])/diff
	return norm_features

def _float_feature(value):
	return tf.train.Feature(float_list=tf.train.FloatList(value=value))

def _bytes_feature(value):
	if isinstance(value, type(tf.constant(0))):
		value = value.numpy() # BytesList won't unpack a string from an EagerTensor.
	return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))

def import_data_from_text_export_to_tfrecord(filename, train_output_filename, test_output_filename, norm_params, label_dict, train_test_split, norm_option = "log_min_max"):
	with open(filename, 'r') as f:
		with tf.io.TFRecordWriter(train_output_filename) as train_writer:
			with tf.io.TFRecordWriter(test_output_filename) as test_writer:
				for line in f:
					items = line.rstrip().split('\t')
					sample_id = items[0]
					x_features = np.array([float(items[i]) if items[i] != "NA" else 0.0 for i in range(1, len(items))])
					x_norm_features = normalize_features(x_features, norm_params, norm_option)
					y_labels = label_dict.get(sample_id)
					if not type(y_labels) is np.ndarray:
						print(sample_id)
						continue
					feature_dict = {"features": _float_feature(x_norm_features), "labels": _float_feature(y_labels)}
					sample = tf.train.Example(features = tf.train.Features(feature = feature_dict))
					sample_train_test_split = train_test_split.get(sample_id)
					if sample_train_test_split == '1':
						train_writer.write(sample.SerializeToString())
					elif sample_train_test_split == '2':
						test_writer.write(sample.SerializeToString())
					elif not sample_train_test_split:
						continue
	return




def import_data_from_text_export_to_tfrecord_testonly(filename, test_output_filename, norm_params, label_dict, norm_option = "log_min_max"):
	with open(filename, 'r') as f:
		with tf.io.TFRecordWriter(test_output_filename) as test_writer:
			for line in f:
				items = line.rstrip().split('\t')
				sample_id = items[0]
				x_features = np.array([float(items[i]) if items[i] != "NA" and items[i] != "" else 0.0 for i in range(1, len(items))])
				x_norm_features = normalize_features(x_features, norm_params, norm_option)
#				y_labels = label_dict[sample_id]
#				feature_dict = {"features": _float_feature(x_norm_features), "labels": _float_feature(y_labels)}
				feature_dict = {"features": _float_feature(x_norm_features)}
				sample = tf.train.Example(features = tf.train.Features(feature = feature_dict))
				test_writer.write(sample.SerializeToString())
	return






def _parse_function(example_proto, feature_shape):
	features = {'features': tf.io.FixedLenFeature(shape=feature_shape, dtype=tf.float32),'labels': tf.io.FixedLenFeature(shape=LABEL_SHAPE, dtype=tf.float32)}
	return tf.io.parse_single_example(example_proto, features)


def _parse_function_no_label(example_proto, feature_shape):
	features = {'features': tf.io.FixedLenFeature(shape=feature_shape, dtype=tf.float32)}
	return tf.io.parse_single_example(example_proto, features)


def load_test_tfrecorddataset(test_filename, feature_shape, has_label = True):
	if has_label:
		test_dataset = tf.data.TFRecordDataset(test_filename).map(lambda x: _parse_function(x, feature_shape) )
		X_test = []
		y_test = []
		for record in test_dataset:
			X_test.append(record["features"].numpy())
			y_test.append(record["labels"].numpy())
		X_test = tf.cast(np.array(X_test), tf.float32)
		y_test = tf.cast(np.array(y_test), tf.float32)
		return X_test, y_test
	else:
		test_dataset = tf.data.TFRecordDataset(test_filename).map(lambda x: _parse_function_no_label(x, feature_shape) )
		X_test = []
		for record in test_dataset:
			X_test.append(record["features"].numpy())
		X_test = tf.cast(np.array(X_test), tf.float32)
		return X_test





####################################################################################
####################################################################################
## Multiple response variable neural network
####################################################################################
####################################################################################


def compute_mae_loss(logits, targets):
	loss = tf.reduce_mean(input_tensor=tf.math.abs(logits - targets))
	return loss


def sample_scaling(x, scaling_option = "log_min_max"):
	if True: #scaling_option == "log_min_max":
		x = np.log2(x + 1)
		mms = pp.MinMaxScaler(feature_range=(0, 1), copy=True)
		x = mms.fit_transform(x.T).T
	return x, mms


def sample_scaling_test(x, scaler):
	x = np.log2(x + 1)
	x = scaler.transform(x.T).T
	return x



def output_during_training_results(model, loss_test, y_pred, step, output_prefix):
	np.savetxt(output_prefix + '_' + str(step) + ".y_pred.txt", y_pred, delimiter = "\t", fmt = "%f")
	#model.save(output_prefix + '_' + str(step) + ".model.h5")
	g = open(output_prefix + '_' + str(step) + ".loss_test.txt", 'w')
	g.write(str(float(loss_test)) + '\n')
	g.close()
	return


def output_after_training_results(model, loss_test_list, loss_train_list, y_pred, output_prefix):
	g = open(output_prefix + ".loss_test_train.txt", 'w')
	for i in range(len(loss_test_list)):
		g.write(str(loss_test_list[i]) + '\t' + str(loss_train_list[i]) + '\n')
	g.close()
	#model.save(output_prefix + ".model.h5")
	np.savetxt(output_prefix + ".y_pred.txt", y_pred, delimiter = "\t", fmt = "%f")
	return



def compute_loss(logits, targets, loss_function):
	if loss_function == "mae":
		loss = compute_mae_loss(logits, targets)
	elif loss_function == "mse":
		loss = compute_mse_loss(logits, targets)
	elif loss_function == "mape":
		loss = compute_mape_loss(logits, targets)
	elif loss_function == "kl":
		loss = compute_kl_loss(logits, targets)
	return loss





def predict_on_test_set_with_trained_model(model_filename, test_filename, pred_filename, feature_shape, has_label = True):
	if has_label:
		X_test, y_test = load_test_tfrecorddataset(test_filename, feature_shape)
		model = tf.keras.models.load_model(model_filename)
		y_pred = model.predict(X_test)
		np.savetxt(pred_filename, y_pred, delimiter = "\t", fmt = "%f")
	else:
		X_test = load_test_tfrecorddataset(test_filename, feature_shape, has_label)
		model = tf.keras.models.load_model(model_filename)
		y_pred = model.predict(X_test)
		np.savetxt(pred_filename, y_pred, delimiter = "\t", fmt = "%f")
	return


def load_initializer(initializer_id):
	if initializer_id == "random_normal":
		return tf.keras.initializers.RandomNormal(mean=0., stddev=1.)
	elif initializer_id == "random_uniform":
		return tf.keras.initializers.RandomUniform(minval=-0.05, maxval=0.05)
	elif initializer_id == "glorot_normal":
		return tf.keras.initializers.GlorotNormal()
	elif initializer_id == "glorot_uniform":
		return tf.keras.initializers.GlorotUniform()



def NN_multiple_response_batch_variedLR_startMiddle(model_parameter, train_filename, test_filename, output_prefix, feature_shape, fixed_learning_rate = False):
	n_classes = NUM_TISSUES
	if model_parameter["new_model"]:
		model = tf.keras.Sequential()
		n_layers = len(model_parameter["num_nodes"])
		initializer = load_initializer(model_parameter["initializer"])
		if initializer == "default":
			for i_layer in range(n_layers):
				model.add(tf.keras.layers.Dense(model_parameter["num_nodes"][i_layer]))
				model.add(tf.keras.layers.BatchNormalization())
				model.add(tf.keras.layers.Activation(tf.nn.relu))
				model.add(tf.keras.layers.Dropout(model_parameter["rate_dropout"][i_layer]))
			model.add(tf.keras.layers.Dense(n_classes, activation=tf.nn.softmax))
		else:
			for i_layer in range(n_layers):
				model.add(tf.keras.layers.Dense(model_parameter["num_nodes"][i_layer], kernel_initializer=initializer))
				model.add(tf.keras.layers.BatchNormalization())
				model.add(tf.keras.layers.Activation(tf.nn.relu))
				model.add(tf.keras.layers.Dropout(model_parameter["rate_dropout"][i_layer]))
			model.add(tf.keras.layers.Dense(n_classes, activation=tf.nn.softmax, kernel_initializer=initializer))
	else:
		EXISTING_MODEL = model_parameter["existing_model"]
		model = tf.keras.models.load_model(EXISTING_MODEL)
		print("LOAD EXISTING MODEL:", EXISTING_MODEL)
	### LOAD DATA
	train_dataset = tf.data.TFRecordDataset(train_filename).map(lambda x: _parse_function(x, feature_shape))

	### TRAIN MODEL
	patience = 10
	wait = 0
	best = 100000
	EPOCH_START = model_parameter["epoch_start"]
	if fixed_learning_rate:
		optimizer = tf.keras.optimizers.Adam(learning_rate = model_parameter["learning_rate"])
	else:
		optimizer = tf.keras.optimizers.Adam(tf.keras.optimizers.schedules.ExponentialDecay(initial_learning_rate=model_parameter["learning_rate"], decay_steps = model_parameter["step"], decay_rate = 0.95))
	X_test, y_test = load_test_tfrecorddataset(test_filename, feature_shape)
	loss_train_list = []
	loss_test_list = []
	np.savetxt(output_prefix + ".y_true.txt", y_test, delimiter = "\t", fmt = "%f")
	for epoch in range(EPOCH_START, model_parameter["epoch"]):
		print("\nStart of epoch %d" % (epoch,))
		start_time = time.time()
		data = train_dataset.shuffle(500000).repeat(1).batch(batch_size=model_parameter["batch_size"])
		data_iter = iter(data)
		tmp_loss_train_list = []
		tmp_loss_test_list = []
		for step in range(model_parameter["step"]):
			batch = data_iter.get_next()
			x = batch["features"]
			y = batch["labels"]
			with tf.GradientTape() as tape:
				logits = model(x, training=True)
				loss = compute_loss(logits, y, model_parameter["loss"])
			grads = tape.gradient(loss, model.trainable_weights)
			optimizer.apply_gradients(zip(grads, model.trainable_weights))
			y_pred = model.predict(X_test)
			loss_test = compute_loss(y_pred, y_test, model_parameter["loss"])
			tmp_loss_train_list.append(float(loss))
			tmp_loss_test_list.append(float(loss_test))
			# Log every 200 batches.
			if step % 200 == 0:
				print("Training loss (for one batch) at step %d: %.8f" % (step, float(loss)))
				print("Validation loss (for one batch) at step %d: %.8f" % (step, float(loss_test)))
		avg_training_loss = float(np.mean(tmp_loss_train_list))
		avg_validation_loss = float(np.mean(tmp_loss_test_list))
		print("Average training loss over epoch: %.8f" % (avg_training_loss))
		print("Average validation loss over epoch: %.8f" % (avg_validation_loss))
		loss_train_list.append(avg_training_loss)
		loss_test_list.append(avg_validation_loss)
		gc.collect()
		wait += 1
		output_during_training_results(model, loss_test, y_pred, epoch, output_prefix)
		if avg_validation_loss < best:
			best = avg_validation_loss
			model.save(output_prefix + ".model.h5")
			wait = 0
		if wait >= patience:
			break

	output_after_training_results(model, loss_test_list, loss_train_list, y_pred, output_prefix)
	return





