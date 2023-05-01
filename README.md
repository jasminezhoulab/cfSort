# cfSort

## Overview

`cfSort` is a supervised approach for sensitive and accurate tissue deconvolution in cfDNA. It is built upon comprehensive tissue markers derived from 521 tissue samples collected from GTEx project, covering 29 major human tissue types. `cfSort` is implemented in python using tensorflow. Here we provide the core codes and the pre-trained DNN models.

`cfSort` can be freely used for educational and research purposes by non-profit institutions and U.S. government agencies only under the UCLA Academic Software License. For information on the use for a commercial purpose or by a commercial or for-profit entity, please contact Prof. Xianghong Jasmine Zhou (XJZhou@mednet.ucla.edu).

## Pre-requisite packages

cfSort was developed on UCLA Hoffman2 cluster. All the analyses in this study were done there. The system information of UCLA Hoffman2 cluster is:
Linux login1 3.10.0-1160.80.1.el7.x86_64 #1 SMP Tue Nov 8 15:48:59 UTC 2022 x86_64 x86_64 x86_64 GNU/Linux
There is no non-standard hardware required.

1) python3 version 3.8.5
2) numpy 1.23.4
3) scipy 1.8.0
4) tensorflow 2.8.0


## Usage

There are two core scripts: 
1) model_train.py, to train the DNN models
2) model_pred.py, to estimate tissue composition of testing samples using the pre-trained DNN models.

### Train DNN models

model_train.py requires the following inputs:
1) train_file: training data in the tfrecord format, where each sample is formatted as \{"features": \<numerical features at markers\>, label: \<ground truth tissue composition\>\}.
2) test_file: validation data in the tfrecord format, where each sample is formatted as \{"features": \<numerical features at markers\>, label: \<ground truth tissue composition\>\}.
3) samples_list_file: a list of sample ids used in the training.
4) num_nodes: number of nodes in each hidden layer, comma delimited, e.g. "1024,512,32".
5) rate_dropout: dropout rates in each hidden layer, comma delimited, e.g. "0.05,0.01,0". Note that the numbers in rate_dropout has to be consistent with the number in num_nodes.
6) step: number of steps per epoch. This number can be automatically learned at the beginning of the training.
7) learning_rate: starting learning rate of the DNN, e.g. 0.001
8) batch_size: number of samplers per minibatch, e.g. 32
9) loss: name of the loss function, e.g. mae
10) optimizer: optimizer, e.g. Adam.
11) output_prefix: prefix of output files
12) epoch: max epochs in the training process
13) num_features: number of numerical features in the training and testing data
14) fixed_learning_rate: False
15) epoch_start: 0
16) model_file: path to a model ".h5" file.
17) new_model: boolean, if training starts from a new model. If False, a model_file is required.
18) is_cpu: boolean, if training is done on CPU
19) initializer: initializer used in the dense layers
20) global_seed: random seed
21) train_size: number of training samples

model_train.py will output logging files (\<output_prefix\>.\*.txt) and a final model file (\<output_prefix\>.model.h5).


### Predict with pre-trained DNNs

model_pred.py requires the following inputs:
1) test_filename: testing data in the tfrecord format, where each sample is formatted as \{"features": \<numerical features at markers\>\}.
2) model_filename: pre-trained DNN model, e.g. DNN1.h5
3) pred_filename: output file name
4) num_features: number of numerical features in the testing data
5) has_label: False

model_pred.py saves the predictions in the pred_filename. The predictions is in tab-delimited format, containing the tissue compositions of the testing samples in the same order as in the testing data tfrecord file. 

We provide two pre-trained DNNs which are the component DNNs of the cfSort: DNN1.h5 and DNN2.h5. Due to the limit of file size on GitHub, we deposit these two models on https://zenodo.org/record/7884243. cfSort ensembles the two DNNs by averaging the tissue compositions.


## Citation
