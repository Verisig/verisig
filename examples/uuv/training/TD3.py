import copy
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from keras.models import Sequential
from keras.layers import Dense
from keras import optimizers
from keras import models

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Implementation of Twin Delayed Deep Deterministic Policy Gradients (TD3)
# Paper: https://arxiv.org/abs/1802.09477

class UnitNormClipper(object):

        def __init__(self, frequency=5):
                self.frequency = frequency

        def __call__(self, module):
                # don't normalize the last layer
                if hasattr(module, 'weight') and module.weight.data.shape[0] > 1:
                        w = module.weight.data
                        w.div_(torch.norm(w, 2).expand_as(w))

                if hasattr(module, 'bias') and module.bias.data.shape[0] > 1:
                        w = module.bias.data
                        w.div_(torch.norm(w, 1).expand_as(w))                        

class Actor(nn.Module):
        def __init__(self, state_dim, action_dim, max_action):
                super(Actor, self).__init__()

                layer_size = 32
                
                self.l1 = nn.Linear(state_dim, layer_size)
                self.l2 = nn.Linear(layer_size, layer_size)
                self.l3 = nn.Linear(layer_size, action_dim)

                # self.l1 = nn.Conv1d(1, 16, 4)
                # self.l2 = nn.Conv1d(16, 16, 4)
                # self.l3 = nn.Linear(240, 64)
                # self.l4 = nn.Linear(64, 1)
                
                self.max_action = max_action
                

        def forward(self, state):

                a = F.tanh(self.l1(state))
                a = F.tanh(self.l2(a))
                
                return self.max_action * torch.tanh(self.l3(a))                


class Critic(nn.Module):
        def __init__(self, state_dim, action_dim):
                super(Critic, self).__init__()

                layer_size = 256
                
                # Q1 architecture
                self.l1 = nn.Linear(state_dim + action_dim, layer_size)
                self.l2 = nn.Linear(layer_size, layer_size)
                self.l3 = nn.Linear(layer_size, 1)

                # self.l1 = nn.Conv1d(1, 16, 4)
                # self.l2 = nn.Conv1d(16, 16, 4)
                # self.l3 = nn.Linear(240 + action_dim, 64)
                # self.l4 = nn.Linear(64, 1)

                # Q2 architecture
                self.l4 = nn.Linear(state_dim + action_dim, layer_size)
                self.l5 = nn.Linear(layer_size, layer_size)
                self.l6 = nn.Linear(layer_size, 1)

                # self.l5 = nn.Conv1d(1, 16, 4)
                # self.l6 = nn.Conv1d(16, 16, 4)
                # self.l7 = nn.Linear(240 + action_dim, 64)
                # self.l8 = nn.Linear(64, 1)


        def forward(self, state, action):

                sa = torch.cat([state, action], 1)
                sa_ext = sa.unsqueeze(1)
                state_ext = state.unsqueeze(1)

                # q1 = F.relu(self.l1(state_ext))
                # q1 = F.relu(self.l2(q1))
                # q1 = q1.view(-1, 240)
                # q1 = torch.cat([q1, action], 1)
                # q1 = F.relu(self.l3(q1))
                # q1 = self.l4(q1)

                q1 = F.relu(self.l1(sa))
                q1 = F.relu(self.l2(q1))
                q1 = self.l3(q1)
                

                # q2 = F.relu(self.l5(state_ext))
                # q2 = F.relu(self.l6(q2))
                # q2 = q2.view(-1, 240)
                # q2 = torch.cat([q2, action], 1)
                # q2 = F.relu(self.l7(q2))
                # q2 = self.l8(q2)

                q2 = F.relu(self.l4(sa))
                q2 = F.relu(self.l5(q2))
                q2 = self.l6(q2)

                
                return q1, q2


        def Q1(self, state, action):
                sa = torch.cat([state, action], 1)
                sa_ext = sa.unsqueeze(1)
                state_ext = state.unsqueeze(1)

                # q1 = F.relu(self.l1(state_ext))
                # q1 = F.relu(self.l2(q1))
                # q1 = q1.view(-1, 240)
                # q1 = torch.cat([q1, action], 1)
                # q1 = F.relu(self.l3(q1))
                # q1 = self.l4(q1)

                q1 = F.relu(self.l1(sa))
                q1 = F.relu(self.l2(q1))
                q1 = self.l3(q1)                
                return q1


class TD3(object):
        def __init__(
                self,
                state_dim,
                action_dim,
                max_action,
                discount=0.99,
                tau=0.005,
                policy_noise=0.2,
                noise_clip=0.5,
                policy_freq=2
        ):

                self.actor = Actor(state_dim, action_dim, max_action).to(device)
                self.actor_target = copy.deepcopy(self.actor)
                self.actor_optimizer = torch.optim.Adam(self.actor.parameters(), lr=3e-4)

                self.critic = Critic(state_dim, action_dim).to(device)
                self.critic_target = copy.deepcopy(self.critic)
                self.critic_optimizer = torch.optim.Adam(self.critic.parameters(), lr=3e-4)

                self.max_action = max_action
                self.discount = discount
                self.tau = tau
                self.policy_noise = policy_noise
                self.noise_clip = noise_clip
                self.policy_freq = policy_freq

                self.total_it = 0

                self.clipper = UnitNormClipper()


        def select_action(self, state):
                
                state = torch.FloatTensor(state.reshape(1, -1)).to(device)
                return self.actor(state).cpu().data.numpy().flatten()


        def train(self, replay_buffer, batch_size=100):
                self.total_it += 1

                # Sample replay buffer 
                state, action, next_state, reward, not_done = replay_buffer.sample(batch_size)

                with torch.no_grad():
                        # Select action according to policy and add clipped noise
                        noise = (
                                torch.randn_like(action) * self.policy_noise
                        ).clamp(-self.noise_clip, self.noise_clip)
                        
                        next_action = (
                                self.actor_target(next_state) + noise
                        ).clamp(-self.max_action, self.max_action)

                        # Compute the target Q value
                        target_Q1, target_Q2 = self.critic_target(next_state, next_action)
                        target_Q = torch.min(target_Q1, target_Q2)
                        target_Q = reward + not_done * self.discount * target_Q

                # Get current Q estimates
                current_Q1, current_Q2 = self.critic(state, action)

                # Compute critic loss
                critic_loss = F.mse_loss(current_Q1, target_Q) + F.mse_loss(current_Q2, target_Q)

                # Optimize the critic
                self.critic_optimizer.zero_grad()
                critic_loss.backward()
                self.critic_optimizer.step()

                # Delayed policy updates
                if self.total_it % self.policy_freq == 0:

                        # Compute actor losse
                        actor_loss = -self.critic.Q1(state, self.actor(state)).mean()
                        
                        # Optimize the actor 
                        self.actor_optimizer.zero_grad()
                        actor_loss.backward()
                        self.actor_optimizer.step()

                        # normalize actor weights and biases
                        self.actor.apply(self.clipper)

                        # Update the frozen target models
                        for param, target_param in zip(self.critic.parameters(), self.critic_target.parameters()):
                                target_param.data.copy_(self.tau * param.data + (1 - self.tau) * target_param.data)

                        for param, target_param in zip(self.actor.parameters(), self.actor_target.parameters()):
                                target_param.data.copy_(self.tau * param.data + (1 - self.tau) * target_param.data)


        def save(self, filename):

                kmodel = Sequential()

                kmodel.add(Dense(units=self.actor.l1.weight.shape[0], activation='tanh', weights = [self.actor.l1.weight.data.numpy().T, self.actor.l1.bias.data.numpy()], input_dim=self.actor.l1.weight.shape[1]))
                kmodel.add(Dense(units=self.actor.l2.weight.shape[0], activation='tanh', weights = [self.actor.l2.weight.data.numpy().T, self.actor.l2.bias.data.numpy()]))
                kmodel.add(Dense(units=1, activation = 'tanh', weights = [self.actor.l3.weight.data.numpy().T, self.actor.l3.bias.data.numpy()]))
                optimizer = optimizers.RMSprop(lr=0.00025, rho=0.9, epsilon=1e-06)
                kmodel.compile(loss="mse", optimizer=optimizer)
                
                kmodel.save(filename + str(self.actor.l1.weight.shape[0]) +\
                                           'x' + str(self.actor.l2.weight.shape[0]) + '.h5')

                      
                #torch.save(self.critic.state_dict(), filename + "_critic")
                #torch.save(self.critic_optimizer.state_dict(), filename + "_critic_optimizer")
                #torch.save(self.actor.state_dict(), filename + "_actor")
                #torch.save(self.actor_optimizer.state_dict(), filename + "_actor_optimizer")


        def load(self, filename):
                self.critic.load_state_dict(torch.load(filename + "_critic"))
                self.critic_optimizer.load_state_dict(torch.load(filename + "_critic_optimizer"))
                self.actor.load_state_dict(torch.load(filename + "_actor"))
                self.actor_optimizer.load_state_dict(torch.load(filename + "_actor_optimizer"))
