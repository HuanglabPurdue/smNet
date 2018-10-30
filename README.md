# smNet
This software is distributed as accompanying software for manuscript:*'Analysing complex single molecule emission patterns with deep learning'* by Peiyi Zhang, Sheng, Liu, Abhishek Chaurasia, Donghan Ma, Michael J. Mlodzianoski, Eugenio Culurciello and Fang Huang.
## Instructions for Single Molecule Deep Neural Network ('smNet')
**The code has been developed and tested on the following system and packages:**\
Ubuntu16.04, Torch7, CUDA8.0, cuDNN-8.0, NCCL1, MatlabR2015a, mattorch1.0-0
### Packages Installation Instruction
#### 1. Install Torch
Run the following commands in a terminal:
```
git clone https://github.com/torch/distro.git ~/torch --recursive
cd ~/torch; bash install-deps;
./install.sh
source ~/.bashrc
luarocks install image
luarocks list
```
(reference to http://torch.ch/docs/getting-started.html#)
