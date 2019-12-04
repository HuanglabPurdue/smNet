# smNet (Single Molecule Deep Neural Network)
This software is distributed as accompanying software for manuscript: P. Zhang, S. Liu, A. Chaurasia, D. Ma, M. J. Mlodzianoski, E. Culurciello and F. Huang, "Analyzing complex single molecule emission patterns with deep learning"(2018) **Nature Methods**, Advanced Online Publication, doi: https://doi.org/10.1038/s41592-018-0153-5


## Files included in this package
### Content of smNet software (Pytorch)
**Sample Data:**
* train/data.mat: A small training dataset containing simulated PSFs.
* train/label.mat: Labels of the training dataset (ground truth of the parameters)
* train/CRLB.mat: Calculated CRLB used in training
* test/data.mat: A small test dataset
* test/label.mat: Underlying true position to compare with smNet results
* test/CRLB.mat: Underlying true position to compare with smNet results

**smNet(Pytorch) Source Code:** coded in Python
* main.py: Main script for training smNet
* model.py: Script for smNet architecture definition
* opts.py: Definitions of user adjustable variables
* test.py: Script for testing smNet
* /setup/loss.py: Script that includes CRLB weighted loss function definition


### Content of smNet software (torch)
**Sample Data:**
* train/data.mat: A small training dataset containing simulated PSFs.
* train/label.mat: Labels of the training dataset (ground truth of the parameters)
* train/CRLB.mat: Calculated CRLB used in training
* test/testdata.mat: A small test dataset
* test/testlabel.mat: Underlying true positions to compare with smNet results
* test/CRLB.mat: CRLB to evaluate precision performance

**smNet(torch) Source Code:** coded in Lua (http://www.lua.org/)
* main.lua: Main script for training smNet
* opts.lua: Definitions of user adjustable variables
* test/test.lua: Script for testing smNet
* load/loaddata.lua: Script for loading training and validation data
* load/loadtestdata.lua: Script for loading test data
* load/loadweight.lua: Script for loading CRLB for training

**Note:**: smNet can be trained with >1M images. Additional training and validation datasets are available upon request. We are working on a complete version of smNet in python. We will update shortly.

### Content of PSF toolbox
**Matlab classes:**
* OptimPR_Ast -- for phase retrieval
* PRPSF -- for phase retrieval
* PSF_zernike -- simulation of PSFs with optical aberration
* PSF_interp -- simulation of PSFs from interpolation
* PSF_DH -- simulation of double-helix PSFs
* DipoleField -- simulation of dipole PSFs
* CalCRLB -- calculation of CRLB for non-dipole PSF model
* CalCRLB_dipole -- calculation of CRLB for dipole PSF model
* Zernike_Polynomials -- generation of Zernike polynomials
* OTFrescale -- OTF rescaling of the simulated PSFs

**Example codes:**
* PR_example.m
* PSF_astigmatism_example.m
* PSF_complex_simulate_example.m
* PSF_complex_interp_example.m
* PSF_doublehelix_example.m
* PSF_dipole_example.m

**Test data:**
* complexPSF_recorded.mat
* astigmatimPSF_PR_result_optimAst.mat
* complexPSF_samplepsf.mat

## Instruction on Using smNet software (Pytorch)
**The code has been developed and tested on the following system and packages:**\
Ubuntu16.04LTS, Python3.6.9, Pytorch0.4.0, CUDA10.1, MatlabR2015a
(Detailed package installation instructions are provided in the following text)

### 1. Setup
Copy loss.py and __init__.py from 'setup/' to '.../nn/modules/', which is usually under following directory:
/home/username/anaconda3/lib/python3.6/site-packages/torch/nn/modules/

### 2. Training
Run the following command in a terminal:
```
python main.py --datapath '/SampleData directory/SampleData/xyz/train/' --save 'result directory' --imHeight 32 --imWidth 32
--batchsize 20 --datasize 10000 --mode z --crange 8 --foldernum 1 --weighting 1
```
**Note:**
1. Each iteration will save a model named by the iteration number.
2. The user could open errorplot.png from the result directory to observe the evolution of training
and validation errors.
3. imHeight, imWidth, datasize, mode and save are user adjustable variables. Mode here includes xy, z, aberration, alpha and beta. See opts.py for the
definitions of other user adjustable variables.
4. The user adjustable variables for training will be saved in /result directory/opt.txt
5. The training and validation errors for each iteration will be saved in /result directory/error.log (The 1st column is training error and the 2nd column is validation error)

### 3. Testing
Run the following commands line by line in a terminal:
```
python test.py --datapath '/SampleData directory/SampleData/xyz/test' --save 'result directory'  --checkpoint 100 --mode z
```
**Note:**
1. Result directory is the same for both training and testing
2. The ‘checkpoint’ is the iteration number at the stop criterion from the training step. See **Supplementary Note 4.1** in smNet Manuscript (https://www.nature.com/articles/s41592-018-0153-5#Sec13)

## Instruction on Using smNet software (torch)
**The code has been developed and tested on the following system and packages:**\
Ubuntu16.04, Torch7, CUDA8.0, cuDNN-8.0, NCCL1, MatlabR2015a, mattorch1.0-0\
(Detailed package installation instructions are provided in the following text)

### 1. Training
Run the following command in a terminal:
```
th main.lua --datapath 'SampleData directory/SampleData/train/' --imHeight 32 --imWidth 32 --datasize
10000 --mode z --save ‘result directory’
```
**Note:**
1. Each iteration will save a model named by the iteration number.
2. The user could open error.log.eps from the result directory to observe the evolution of training
and validation errors.
3. imHeight, imWidth, datasize, mode and save are user adjustable variables. See opts.lua for the
definitions of other user adjustable variables.

### 2. Testing
Run the following commands line by line in a terminal:
```
cd test
th test.lua --testdatapath 'SampleData directory/SampleData/test/' --imHeight 32 --imWidth 32 --
numSubregion 410 --mode z --save ‘result directory’ --modelNum 200 --batchSize 200
```
**Note:**
1. Result directory is the same for both training and testing
2. The ‘modelNum’ is the iteration number at the stop criterion from the training step. See **Supplementary Note 4.1** in smNet Manuscript (https://www.nature.com/articles/s41592-018-0153-5#Sec13)
3. Increasing batch size (‘batchSize’) during testing increases speed of forward propagation.

## Instruction on Using PSF toolbox
**The code has been tested in following system and packages:**\
Windows7, Windows10, MatlabR2015a, DIPimage(http://www.diplib.org/)

1. Change current folder in Matlab to PSF toolbox.
2. Run each example code in PSF toolbox/examples/.\
**Note: ** for PR_example.m, it requires user to select the center of the PSF in the pop up window.
3. Type ‘help classname’ in Matlab command window for detailed help on each Matlab class.

## Packages Installation Instruction
### 1. Install Torch
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
### 2. Install Nvidia Driver
Download "NVIDIA-Linux-x86_64-[version#].run" and run the following commands:
```
cd ~/Downloads
sudo service lightdm stop
./NVIDIA-Linux-x86_64-[version#].run
```
### 3. Install CUDA
Download cuda "runfile(local)". And run the following commands in the terminal:
```
cd ~/Downloads
sudo sh cuda_[filename].run
```
Follow the installation instructions (do not install Nvidia driver). After installation finished:
```
vim ~/.bashrc
```
And then append the following lines to the opened file to define the CUDA Toolkit PATH variables:
```
#CUDA Toolkit
export CUDA_HOME=/usr/local/cuda-8.0
export LD_LIBRARY_PATH=${CUDA_HOME}/lib64:$LD_LIBRARY_PATH
export PATH=${CUDA_HOME}/bin:${PATH}
```
Save and close the file. Then run the following command in terminal:
```
source ~/.bashrc
```
### 4. Install CUDNN
Download cudnn-8.0-linux-x64-v6.0.tgz. And run the following commands line by line in the terminal:
```
cd Downloads
tar -zxf cudnn-8.0-linux-x64-v6.0.tgz
cd cuda
sudo cp lib64/* /usr/local/cuda/lib64/
sudo cp include/* /usr/local/cuda/include/
luarocks install cudnn
```
### 5. Install nccl
Run the following commands line by line in the terminal:
```
git clone https://github.com/NVIDIA/nccl
cd nccl
make install
sudo ldconfig
luarocks install nccl
```
### 6. Install Matlab
Download Matlab R2015a, unzip the package and run the following commands line by line in the terminal:
```
sudo ./install
sudo apt-get install matlab-support
```
(reference to https://help.ubuntu.com/community/MATLAB)
### 7. Install mattorch (After Matlab is installed)
Run the following commands line by line in the terminal:
```
git clone https://github.com/clementfarabet/lua---mattorch.git
cd lua---mattorch
mkdir build
cd build
cmake -DMATLAB_ROOT=/usr/local/MATLAB/R2015a/ ..
cd ..
luarocks make
```
### 8. Install gnuplot
Run the following commands line by line in the terminal:
```
sudo apt-get install gnuplot
sudo apt-get install gnuplot-qt
```
## Acknowledgement:
We thank K. A. Lidke and M. J. Wester from University of New Mexico for initial contribution to the
PSF toolbox software. We thank E. B. Kromann from Technical University of Denmark for sharing the
phase unwrapping code.

