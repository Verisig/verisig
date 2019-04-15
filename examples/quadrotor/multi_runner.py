#!/usr/bin/python3

import subprocess

verisig_path = '../../verisig'
flowstar_path = '../../flowstar/flowstar'

legend = ['X1_LOWER_BOUND', 'X1_UPPER_BOUND', 'X2_LOWER_BOUND', 'X2_UPPER_BOUND']
test_set = [
    [-0.05, -0.025, -0.05, -0.025],
    [-0.025, 0, -0.05, -0.025],
    [0, 0.025, -0.05, -0.025],
    [0.025, 0.05, -0.05, -0.025],
    [-0.05, -0.025, -0.025, 0],
    [-0.025, 0, -0.025, 0],
    [0, 0.025, -0.025, 0],
    [0.025, 0.05, -0.025, 0],
    [-0.05, -0.025, 0, 0.025],
    [-0.025, 0, 0, 0.025],
    [0, 0.025, 0, 0.025],
    [0.025, 0.05, 0, 0.025],
    [-0.05, -0.025, 0.025, 0.05],
    [-0.025, 0, 0.025, 0.05],
    [0, 0.025, 0.025, 0.05],
    [0.025, 0.05, 0.025, 0.05],
]

print("Building the base model...")
subprocess.run([verisig_path, '-vc=quadrotor_MPC_multi.yml', '-o' ,'-nf', 'quadrotor_MPC.xml', 'tanh20x20.yml'])

with open('quadrotor_MPC.model', 'r') as f:
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
