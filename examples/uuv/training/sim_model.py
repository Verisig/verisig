from UUV import World
import numpy as np
import random
from keras import models

MAX_EP_STEPS = 100

if __name__ == '__main__':
    model = models.load_model('models/tanh_3_32x32.h5')

    y_pos = 28.01
    
    env = World(0.0, y_pos, 0.0, MAX_EP_STEPS)

    rew = 0

    prev_pos_y = env.pos_y

    observation = np.array([0, -y_pos, y_pos])

    u = np.array([0, 0.48556, 45])

    for e in range(MAX_EP_STEPS):

        observation, reward, done, info = env.step(u)
                    
        u = np.radians(5) * model.predict(observation.reshape(1,len(observation)))[0]
        
        rew += reward

        prev_pos_y = env.pos_y

        if done:
            break
    

    print(model.predict(observation.reshape(1,len(observation)))[0])
    print('heading: ' + str(env.heading))
    print('x pos: ' + str(env.pos_x))
    print('y pos: ' + str(env.pos_y))
    #env.plot_trajectory()
    
