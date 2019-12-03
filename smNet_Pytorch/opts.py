
####################################################################################################
##Script for user adjustable variables
##
##(C) Copyright 2019                The Huang Lab
##
##    All rights reserved           Weldon School of Biomedical Engineering
##                                  Purdue University
##                                  West Lafayette, Indiana
##                                  USA
##
##    Peiyi Zhang, Dec 2019
####################################################################################################

import argparse

#Command line options
def options():

    parser = argparse.ArgumentParser()

    #Training Related:
    parser.add_argument('--learningRate', type=float, default=1e-5) #learning rate
    parser.add_argument('--maxepoch', type=int, default=1000)       #maximum number of training iterations
    parser.add_argument('--savenum', type=int, default=1)           #save model per # iterations
    parser.add_argument('--plot', type=int, default=1)              #1: plot training/testing error; 0: no plot
    parser.add_argument('--modelnum', type=int, default=1)          #start to train from xx model (this is in case training stops in the middle)
    parser.add_argument('--trainSplit', type=float, default=0.8)    #percentage of training data
    parser.add_argument('--batchsize', type=int, default=64)       #batch size
    parser.add_argument('--datasize', type=int, default=100000)     #training and validation data size
    parser.add_argument('--foldernum', type=int, default=1)         #data number that you want to concatenate together (this is in case mat file is too large to save)
    parser.add_argument('--crange', type=int, default=20)            #2: The range of HardTanh for training xy (unit pixel); 8: The range of HardTanh for training z (unit 10 x micron)

    #Testing Related
    parser.add_argument('--checkpoint', type=int, default=0)        #iteration number at the stop criterion
    parser.add_argument('--testbatchsize', type=int, default=1000)  #batch for testing

    #Training and Testing Related
    parser.add_argument('--mode', type=str, default='z')            #choose 'xy', 'z', 'alpha', 'beta', 'aberration' for xy, z, polar angle, azimuthal angle and aberration estimation respectively
    parser.add_argument('--weighting', type=int, default=1)         #0: without CRLB weighting; 1: with CRLB weighting
    parser.add_argument('--datapath', type=str)                     #training and validation dataset location
    parser.add_argument('--labelsize', type=int, default=12)        #number of zernike modes
    parser.add_argument('--save', type=str)                         #save trained model, error plot and result here
    parser.add_argument('--imHeight', type=int, default=16)         #image height
    parser.add_argument('--imWidth', type=int, default=16)          #image width
    parser.add_argument('--channel', type=int, default=2)           #2:biplane; 1:single plane
    parser.add_argument('--nGPU', type=int, default=4)              #number of GPUs to use

    opt = parser.parse_args()

    return opt
