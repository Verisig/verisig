#!/usr/bin/python3

import subprocess

verisig_path = '../../verisig'
flowstar_path = '../../flowstar/flowstar'

legend = ['X1_LOWER', 'X1_UPPER']
test_set = [
    [-0.59, -0.585],
    [-0.585, -0.58],
    [-0.58, -0.57],
    [-0.57, -0.55],
    [-0.55, -0.53],
    [-0.53, -0.5],
    [-0.50, -0.48],
    [-0.48, -0.46],
    [-0.46, -0.45],
    [-0.45, -0.44],
    [-0.44, -0.43],
    [-0.43, -0.42],
    [-0.42, -0.415],
    [-0.415, -0.41],
    [-0.41, -0.4]
]

print("Building the base model...")
subprocess.run([verisig_path, '-vc=MC_multi.yml', '-o' ,'-nf', 'MC.xml', 'sig16x16.yml'])

with open('MC.model', 'r') as f:
    model = f.read()

for test in test_set:
    print("=========================================")
    print("Running test with initial conditions of: ")
    test_model = model
    for i in range(len(legend)):
        print(legend[i] + '=' + str(test[i]), end=', ')
        test_model = test_model.replace(legend[i], str(test[i]))
    print()
    print("=========================================")

    subprocess.run(flowstar_path, input=test_model, shell=True, universal_newlines=True)
    print()
    print()
