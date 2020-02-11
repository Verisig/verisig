from UUV import World
import numpy as np
import random
from keras.models import Sequential
from keras import models
from keras import optimizers
from keras.layers.core import Dense, Dropout, Activation
from keras.layers.normalization import BatchNormalization
from keras.layers.advanced_activations import LeakyReLU
from keras.regularizers import l2
from six.moves import cPickle as pickle
import matplotlib.pyplot as plt
import sys

def relu(x):
    relu = np.maximum(0, x)

    return relu

#this is just for testing purposes
def relu_predict(model, inputs):
    weights = {}
    offsets = {}

    layerCount = 1

    for layer in model.layers:
        if len(layer.get_weights()) > 0:
            weights[layerCount] = layer.get_weights()[0]
            offsets[layerCount] = layer.get_weights()[1]

            layerCount += 1

    curNeurons = inputs

    for layer in range(layerCount-1):
        curNeurons = curNeurons.dot(weights[layer + 1]) + offsets[layer + 1]

        if layer <= layerCount - 3:
            curNeurons = relu(curNeurons)

    return curNeurons

#this is just for testing purposes
def tanh_predict(model, inputs):
    weights = {}
    offsets = {}

    layerCount = 1

    for layer in model.layers:
        
        if len(layer.get_weights()) > 0:
            weights[layerCount] = layer.get_weights()[0]
            offsets[layerCount] = layer.get_weights()[1]

            layerCount += 1

    curNeurons = inputs

    for layer in range(layerCount-1):
        curNeurons = curNeurons.dot(weights[layer + 1]) + offsets[layer + 1]

        curNeurons = np.tanh(curNeurons)

    return curNeurons

def sigmoid(x):
    sigm = 1. / (1. + np.exp(-x))

    return sigm

def swish_predict(model, inputs):
    weights = {}
    offsets = {}

    layerCount = 1

    for layer in model.layers:
        if len(layer.get_weights()) > 0:
            weights[layerCount] = layer.get_weights()[0]
            offsets[layerCount] = layer.get_weights()[1]

            layerCount += 1

    curNeurons = inputs

    for layer in range(layerCount-1):
        curNeurons = curNeurons.dot(weights[layer + 1]) + offsets[layer + 1]

        if layer <= layerCount - 3:
            curNeurons = curNeurons * sigmoid(curNeurons)
            #curNeurons = relu(curNeurons)

    return curNeurons

def normalize(s):
    mean = [2.5]
    spread = [5.0]
    return (s - mean) / spread

def main(argv):

    input_filename = argv[0]
    
    model = models.load_model(input_filename)

    y_pos = 25

    episode_length = 500
    
    env = World(0.0, y_pos, 0.0, episode_length)

    rew = 0

    prev_pos_y = env.pos_y

    #observation = env.reset()
    observation = np.array([0, -y_pos, y_pos])

    u = np.array([0, 0.48556, 45])    

    allX = []
    allY = []
    allR = []

    numTrajectories = 100

    for step in range(numTrajectories):

        observation = env.reset()

        rew = 0

        for e in range(episode_length):

            observation, reward, done, info = env.step(u)
        
            u = np.radians(5) * model.predict(observation.reshape(1,len(observation)))[0]

            if done:
                
                break

            rew += reward

        allX.append(env.allX)
        allY.append(env.allY)
        allR.append(rew)


    #print(np.mean(allR))
    #print('number of crashes: ' + str(num_unsafe))
    
    fig = plt.figure(figsize=(12,10))
    
    #plt.ylim((-1,11))
    #plt.xlim((-1.75,10.25))
    #plt.suptitle('Simulated trajectories of the F1/10 Car', fontsize=30)
    #plt.tick_params(labelsize=20)

    for i in range(numTrajectories):
        plt.plot(allX[i], allY[i], 'r-')

    #plt.savefig('simulations.pdf', format='pdf', bbox_inches = 'tight', pad_inches = 0)
    plt.show()


    #w.plot_lidar()
    
if __name__ == '__main__':
    main(sys.argv[1:])
