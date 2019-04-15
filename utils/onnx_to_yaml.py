#!/usr/bin/python

import onnx
from onnx import numpy_helper
import sys
import yaml

def main(argv):
    input_filename = argv[0]
    output_filename = argv[1]

    model = onnx.load(input_filename)
    weights = model.graph.initializer

    dnn_dict = {}
    dnn_dict['weights'] = {}
    dnn_dict['offsets'] = {}
    dnn_dict['activations'] = {}

    for layer_index in range(len(weights) // 2):
        dnn_dict['weights'][layer_index + 1] = []
        for row in numpy_helper.to_array(weights[layer_index * 2 + 1]):
            a = []
            for column in row:
                a.append(float(column))
            dnn_dict['weights'][layer_index + 1].append(a)
        
        dnn_dict['offsets'][layer_index + 1] = []
        for row in numpy_helper.to_array(weights[(layer_index + 1) * 2]):
            dnn_dict['offsets'][layer_index + 1].append(float(row))

        if layer_index <= (len(weights) // 2) - 2:
            #Assuming Tanh activations for hidden layers
            dnn_dict['activations'][layer_index + 1] = 'Tanh'
        else:
            #Assuming a linear last layer
            dnn_dict['activations'][layer_index + 1] = 'Linear'

    with open(output_filename, 'w') as f:
        yaml.dump(dnn_dict, f)

if __name__ == '__main__':
    main(sys.argv[1:])