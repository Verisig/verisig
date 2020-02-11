import gym
from gym import spaces
import numpy as np
import scipy.io as sio
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math

# car parameters
DEPTH = 45.0
MIN_HEADING_ACTION = math.radians(-5.0)
MAX_HEADING_ACTION = math.radians(5.0)
MIN_SPEED_ACTION = 0.0
MAX_SPEED_ACTION = 1.5433

MIN_SPEED = 0.51444
MAX_SPEED = 2.50
FIXED_SPEED_ACTION = 0.48556

# training parameters
STEP_REWARD_GAIN = 0.5
HEADING_REWARD_GAIN = 5
INPUT_REWARD_GAIN = -0.5
RANGE_REWARD_PENALTY = -0.1
CRASH_PENALTY = -100

GOAL_RANGE = 30.0

MAX_DISTANCE = 60
MIN_DISTANCE = 10

PIPE_LENGTH = 400

class World:

    def __init__(self, pos_x, pos_y, heading, episode_length):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.heading = heading

        model = sio.loadmat('model.mat')

        self.A = model['A']
        self.B = model['B']
        self.C = model['C']
        self.D = model['D']
        self.K = model['K']

        self.x = np.array([[0.0], [0.0], [0.0], [0.0]])
        #u = np.array([[0.0], [2.0], [45.0]])
        self.e = np.array([[0.0], [0.0], [0.0]])

        # step parameters
        self.cur_step = 0
        self.episode_length = episode_length

        # storage
        self.allX = []
        self.allY = []
        self.allX.append(self.pos_x)
        self.allY.append(self.pos_y)


        # parameters needed for consistency with gym environments
        self.obs_low = np.array( [math.radians(-180.0), -1.0 * MAX_DISTANCE, MIN_DISTANCE])
        self.obs_high = np.array( [math.radians(180.0), -1.0 * MIN_DISTANCE, MAX_DISTANCE])

        self.action_space = spaces.Box(low=MIN_HEADING_ACTION, high=MAX_HEADING_ACTION, shape=(1,))
        self.observation_space = spaces.Box(low=self.obs_low, high=self.obs_high)

        self._max_episode_steps = episode_length

    def reset(self):
        self.x = np.array([[0.0], [0.0], [0.0], [0.0]])
        self.e = np.array([[0.0], [0.0], [0.0]])
        
        self.cur_step = 0
        self.pos_x = (np.random.random() * 5.0)
        self.pos_y = 25.0 + (np.random.random() * 10.0)
        self.heading = 0.0 #math.radians(-5.0) + (np.random.random() * math.radians(10.0))

        self.allX = []
        self.allY = []
        self.allX.append(self.pos_x)
        self.allY.append(self.pos_y)

        pipe_heading = -1.0 * self.heading
        stbd_range = self.pos_y / math.cos(self.heading)
        port_range = -1.0 * stbd_range
        measurements = np.array([pipe_heading, port_range, stbd_range])
        return measurements


    def step(self, action):
        self.cur_step += 1

        heading_delta = action[0]
        speed = FIXED_SPEED_ACTION

        # Constrain turning input
        if heading_delta > MAX_HEADING_ACTION:
            heading_delta = MAX_HEADING_ACTION

        if heading_delta < MIN_HEADING_ACTION:
            heading_delta = MIN_HEADING_ACTION

        if speed < MIN_SPEED_ACTION:
            speed = MIN_SPEED_ACTION
        
        if speed > MAX_SPEED_ACTION:
            speed = MAX_SPEED_ACTION
        
        abs_heading = self.heading + heading_delta
        abs_heading = abs_heading if abs_heading < math.pi else abs_heading - (2*math.pi)

        #print(abs_heading)

        u = np.array([[abs_heading], [MIN_SPEED + speed], [45.0]])
        
        y = np.dot(self.C,self.x) + np.dot(self.D,u) + self.e
        self.x = np.dot(self.A,self.x) + np.dot(self.B,u) + np.dot(self.K,self.e)

        self.heading = y[0][0]
        self.heading = self.heading if self.heading < math.pi else self.heading - (2*math.pi)
        self.pos_x += y[1][0] * math.cos(self.heading)
        self.pos_y += y[1][0] * -1.0 * math.sin(self.heading)

        #print(self.heading)
        #print(self.cur_step)

        terminal = False
        reward = 0.0
        
        if self.pos_x > PIPE_LENGTH:
            #print("off the end")
            terminal = True
        elif self.pos_x < -10.0:
            #print("off the beginning")
            terminal = True

        if self.pos_y > MAX_DISTANCE or self.pos_y < MIN_DISTANCE:
            #print("too far")
            terminal = True
            reward += CRASH_PENALTY

        if self.cur_step == self.episode_length:
            terminal = True

        self.allX.append(self.pos_x)
        self.allY.append(self.pos_y)
        

        # Measurements
        pipe_heading = -1.0 * self.heading
        stbd_range = self.pos_y / math.cos( self.heading )
        port_range = -1.0 * stbd_range
        measurements = np.array([pipe_heading, port_range, stbd_range])
        
        # Compute reward
        reward += STEP_REWARD_GAIN
        reward += INPUT_REWARD_GAIN * abs(heading_delta)
        
        if( abs(pipe_heading) < math.radians(5.0) ):
            pass#reward += HEADING_REWARD_GAIN - abs(pipe_heading)

        reward += RANGE_REWARD_PENALTY * abs(stbd_range - GOAL_RANGE)

        return measurements, reward, terminal, -1

    def plot_trajectory(self):
        fig = plt.figure()

        plt.plot(np.array([0.0, PIPE_LENGTH]), np.array([0.0, 0.0]), 'b', linewidth=3)

        plt.plot(self.allX, self.allY, 'r--')

        plt.show()
