#!/usr/bin/python3

import os
import subprocess
import multiprocessing

verisig_path = '../../verisig'
flowstar_path = '../../flowstar/flowstar'
output_path = 'output'
dnn_yaml = 'tanh.yml'

if not os.path.exists(output_path):
        os.mkdir(output_path)

legend = ['X1_LOWER', 'X1_UPPER']

test_set = []

lb = 90
ub = 110
step = 1

# this is a workaround for a precision issue in python that generates an extra instance
num_instances = (ub - lb) / step

count = 0
while count < num_instances:
    test_set.append([lb, lb + step])
    
    lb += step
    count += 1

print("Building the base model...")
subprocess.run([verisig_path, '-vc=ACC_multi.yml', '-o' ,'-nf', 'ACC.xml', dnn_yaml])

with open('ACC.model', 'r') as f:
    model = f.read()


#===========================================================================================
# Begin Parallel Function
#===========================================================================================
def evaluate_conditions(conditions):
    test_model = model
    for i in range(len(legend)):
        test_model = test_model.replace(legend[i], str(conditions[i]))

    with open(output_path + '/ACC_' + str(conditions[0]) + '.txt', 'w') as f:
        subprocess.run(flowstar_path + ' ' + dnn_yaml , input=test_model, shell=True, universal_newlines=True, stdout=f)
#===========================================================================================
# End Parallel Function
#===========================================================================================

print("Starting parallel verification")
num_parallel = multiprocessing.cpu_count() // 2
with multiprocessing.Pool(processes=num_parallel) as pool:
    pool.map(evaluate_conditions, test_set)
