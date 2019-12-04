####################################################################################################
##Script for testing smNet
##
##(C) Copyright 2019                The Huang Lab
##
##    All rights reserved           Weldon School of Biomedical Engineering
##                                  Purdue University
##                                  West Lafayette, Indiana
##                                  USA
##
##    Peiyi Zhang, Feb 2019
####################################################################################################
#import tapackages
import torch
import torch.nn as nn
import torch.nn.parallel
import os
import scipy.io as sio
import numpy as np
from torch.autograd import Variable
import argparse
from model import Net
from opts import options
from tqdm import tqdm
import h5py

#########################################
# define functions
def loaddata(opt, item):
    datadir = os.path.join(opt.datapath, item)
    print('\nloading data from: ', datadir)

    mat_contents = h5py.File(os.path.join(datadir, 'data.mat'),'r')
    dataset = torch.from_numpy(np.array(mat_contents.get('imsp'))).float()
    if dataset.dim() < 4:
        print('(Error) require 4D input. Plese modify your input to be: imsz x imsz x nchannel x datasize')
        exit()

    opt.datasize = dataset.size(0)
    opt.imWidth = dataset.size(2)
    opt.imHeight = dataset.size(3)

    N = opt.datasize
    if N<opt.batchsize:
        opt.batchsize = N

    print('test data size:', dataset.size(0), 'x', dataset.size(1), 'x', dataset.size(2), 'x', dataset.size(3))
    ##normalizing
    for i in range(N):
        dataset[i]=dataset[i]/torch.max(dataset[i])

    return dataset

def test(model, opt, data, item):
    count, totalErr, err = 0, 0, 0
    x = torch.Tensor(opt.batchsize, opt.channel, opt.imWidth, opt.imHeight).cuda()
    model.eval()
    output = torch.Tensor()
    with torch.no_grad():
        for i in tqdm(range(0, opt.datasize, opt.batchsize), ncols = 100, desc = "Testing (# of batches)"):
            for ii in range(opt.batchsize):
                if (i+ii) > opt.datasize - 1:
                    break
                else:
                    x[ii] = data[i+ii]
            xin = Variable(x)
            out = model.forward(xin)
            output = torch.cat([output, out.cpu()], 0)

    savefolder = 'result_from_model' + str(opt.checkpoint)
    savepath = os.path.join(opt.save, savefolder, item)
    if not os.path.exists(savepath):
        os.makedirs(savepath)
    sio.savemat(os.path.join(savepath, 'result.mat'), {opt.mode: output[0:opt.datasize].numpy()})
    print('result saved at: ', savepath)


def main():
    torch.manual_seed(12)

    #load arguments
    opt = options()

    for item in os.listdir(opt.datapath):

        #load data
        data = loaddata(opt, item)

        #model
        model = torch.load(os.path.join(opt.save, 'model{}.pth'.format(opt.checkpoint)))

        #test
        test(model, opt, data, item)

if __name__ == '__main__':
    main()
