# UUV Controller Training Environment

## Setup
The training environment requires python3 and pip3. Use `pip3 install -r requirements.txt` to install the required python3 libraries.

The __model.mat__ file contains the matrix entries for the identified model.

## Training
The command `python3 train_td3.py` will run the training environment, creating __results__ and __models__ folders. The trained controller will be saved under the __models__ folder at regular intervals. A suitable controller should be trained within 10 minutes.

## Visualization
The command `python3 plotTrajectories.py <controller>` can be used to visualize the trajectories created by the controller.
