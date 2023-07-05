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

There are three core scripts: 
1) data_prep.py, to prepare training, validation, and testing data matrices to tensorflow tfrecords.  
2) model_train.py, to train the DNN models
3) model_pred.py, to estimate tissue composition of testing samples using the pre-trained DNN models.


### Prepare data

data_prep.py requires the following inputs:
1) label_filename: the file containing labels of the training and validation data, e.g. demo/input.label.txt. This file contains sample id (1st column), train/validation group id (2nd column, 0 = testing, 1 = training, 2 = validation), and ground-truth tissue fraction (3-31th columns).
2) filename: the data matrix file containing training and validation samples, e.g. demo/input.train_valid.txt. This file contains sample id (1st column), and features (2+ columns). The features were ordered as in the marker file (marker/cfSort_markers.txt.gz, column marker_cluster_index).
3) test_filename: the data matrix file containing testing samples, e.g. demo/input.test.txt. This file has the same format as the training/validation data. 
4) train_output_filename: the output filename for the training data, e.g. demo/output.train.tfrecords.
5) valid_output_filename: the output filename for the validation data, e.g. demo/output.valid.tfrecords.
6) test_output_filename: the output filename for the testing data, e.g. demo/output.test.tfrecords.

usage: python data_prep.py label_filename filename test_filename train_output_filename valid_output_filename test_output_filename

The demo output was in the demo folder.

### Train DNN models

model_train.py requires the following inputs:
1) train_file: training data in the tfrecord format, where each sample is formatted as \{"features": \<numerical features at markers\>, label: \<ground truth tissue composition\>\}.
2) test_file: validation data in the tfrecord format, where each sample is formatted as \{"features": \<numerical features at markers\>, label: \<ground truth tissue composition\>\}.
3) samples_list_file: a list of sample ids used in the training.
4) num_nodes: number of nodes in each hidden layer, comma delimited, e.g. "1024,512,32".
5) rate_dropout: dropout rates in each hidden layer, comma delimited, e.g. "0.05,0.01,0". Note that the numbers in rate_dropout has to be consistent with the number in num_nodes.
6) learning_rate: starting learning rate of the DNN, e.g. 0.001
7) batch_size: number of samplers per minibatch, e.g. 32
8) loss: name of the loss function, e.g. mae
9) optimizer: optimizer, e.g. Adam.
10) output_prefix: prefix of output files
11) epoch: max epochs in the training process
12) global_seed: random seed

model_train.py will output logging files (\<output_prefix\>.\*.txt) and a final model file (\<output_prefix\>.model.h5).

usage: python model_train.py train_file test_file samples_list_file num_nodes rate_dropout learning_rate batch_size loss optimizer output_prefix epoch global_seed

### Predict with pre-trained DNNs

model_pred.py requires the following inputs:
1) test_filename: testing data in the tfrecord format, where each sample is formatted as \{"features": \<numerical features at markers\>\}.
2) model_filename: pre-trained DNN model, e.g. DNN1.h5
3) pred_filename: output file name

model_pred.py saves the predictions in the pred_filename. The predictions is in tab-delimited format, containing the tissue compositions of the testing samples in the same order as in the testing data tfrecord file. 

usage: python model_pred.py test_filename model_filename pred_filename

We provide two pre-trained DNNs which are the component DNNs of the cfSort: DNN1.h5 and DNN2.h5. Due to the limit of file size on GitHub, we deposit these two models on https://zenodo.org/record/7884243. cfSort ensembles the two DNNs by averaging the tissue compositions.




## Citation
Li S, Zeng W, Ni X, Liu Q, Li W, Stackpole ML, Zhou Y, Gower A, Krysan K, Ahuja P, Lu DS, Raman SS, Hsu W, Aberle DR, Magyar CE, French SW, Han SB, Garon EB, Agopian VG, Wong WH, Dubinett SM, Zhou XJ. Comprehensive tissue deconvolution of cell-free DNA by deep learning for disease diagnosis and monitoring. Proc Natl Acad Sci U S A. 2023 Jul 11;120(28):e2305236120. doi: 10.1073/pnas.2305236120. Epub 2023 Jul 3. PMID: 37399400. [link(https://www.pnas.org/doi/10.1073/pnas.2305236120)]
