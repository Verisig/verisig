#!/usr/bin/python

from keras.models import Sequential
from keras import models
from keras import optimizers
from keras.layers.core import Dense, Dropout, Activation
from keras.layers.normalization import BatchNormalization
from keras.layers.advanced_activations import LeakyReLU
from keras.regularizers import l2
import numpy as np
import sys
import yaml

def main(argv):
    input_filename = argv[0]
    output_filename = argv[1]

    model = models.load_model(input_filename)

    dnn_dict = {}
    dnn_dict['weights'] = {}
    dnn_dict['offsets'] = {}
    dnn_dict['activations'] = {}

    layer_count = 1
    for layer in model.layers:
        if len(layer.get_weights()) > 0:
            dnn_dict['weights'][layer_count] = []
            for row in layer.get_weights()[0].T:
                a = []
                for column in row:
                    a.append(float(column))
                dnn_dict['weights'][layer_count].append(a)
            
            dnn_dict['offsets'][layer_count] = []
            for row in layer.get_weights()[1].T:
                dnn_dict['offsets'][layer_count].append(float(row))

            if 'Sigmoid' in str(layer.output):
                dnn_dict['activations'][layer_count] = 'Sigmoid'
            elif 'Tanh' in str(layer.output):
                dnn_dict['activations'][layer_count] = 'Tanh'
            elif 'Relu' in str(layer.output):
                dnn_dict['activations'][layer_count] = 'Relu'
            else:
                dnn_dict['activations'][layer_count] = 'Linear'

            layer_count += 1

    with open(output_filename, 'w') as f:
        yaml.dump(dnn_dict, f)

if __name__ == '__main__':
    main(sys.argv[1:])