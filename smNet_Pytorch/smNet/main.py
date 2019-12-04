####################################################################################################
##Script for training smNet
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
#import packages
import torch
import torch.nn as nn
import torch.nn.parallel
import os
import h5py
import numpy as np
from torch.autograd import Variable
import argparse
from opts import options
from model import Net
from tqdm import tqdm
import scipy.io as sio
import matplotlib.pyplot as plt
##########################################
# define functions
def loaddata(opt, trsize, valsize):
    dataset = torch.Tensor()
    label = torch.Tensor()
    CRLB_std = torch.Tensor()
    for ff in range(1, opt.foldernum+1):
        if opt.foldernum == 1 and ~os.path.exists(os.path.join(opt.datapath, str(ff), 'data.mat')):
            datadir = opt.datapath
        else:
            datadir = os.path.join(opt.datapath, str(ff))
        print('\nloading from: ', datadir)

        # load PSFs
        mat_contents = h5py.File(os.path.join(datadir, 'data.mat'),'r')
        temp = torch.from_numpy(np.array(mat_contents.get('imsp'))).float()
        dataset = torch.cat([dataset, temp],0)
        if dataset.dim() < 4:
            print('(Error) require 4D input. Plese modify your input to be: imsz x imsz x nchannel x datasize')
            exit()

        # load labels
        mat_contents_2 = sio.loadmat(os.path.join(datadir, 'label.mat'))
        if opt.mode == 'xy':
            temp2 = torch.from_numpy(mat_contents_2['xy_in_pixel']).float()
            opt.labelsize = 2
            if not temp2.size(1) == 2:
                print('(Error) require label to be Nx2 for training (x,y) position')
                exit()
        if opt.mode == 'z':
            opt.labelsize = 1
            temp2 = torch.from_numpy(mat_contents_2['z_in_10xmicron']).float()
            if not temp2.size(1) == 1:
                print('(Error) require label to be Nx1 for training z position')
                exit()
        if opt.mode == 'alpha':
            opt.labelsize = 1
            temp2 = torch.from_numpy(mat_contents_2['polar_in_degree']).float()
            if not temp2.size(1) == 1:
                print('(Error) require label to be Nx1 for training polar angle')
                exit()
        if opt.mode == 'beta':
            opt.labelsize = 1
            temp2 = torch.from_numpy(mat_contents_2['azimuthal_in_degree']).float()
            if not temp2.size(1) == 1:
                print('(Error) require label to be Nx1 for training azimuthal angle')
                exit()
        if opt.mode == 'aberration':
            temp2 = torch.from_numpy(mat_contents_2['label']).float()


        label = torch.cat([label, temp2], 0)

        # load CRLB weighting
        if opt.weighting == 1:
            mat_contents_3 = sio.loadmat(os.path.join(datadir, 'CRLB.mat'))
            if opt.mode == 'aberration':
                varname = 'CRLB'
            else:
                varname = 'CRLB' + opt.mode
            temp3 = torch.from_numpy(mat_contents_3[varname]).float()
            CRLB_std = torch.cat([CRLB_std, temp3], 0)

    N = opt.datasize

    print('training for estimating ', opt.mode)
    # normalizing
    for i in range(N):
        dataset[i]=dataset[i]/torch.max(dataset[i])

    # split into train and validation data
    trdata = dataset[0:trsize]
    valdata = dataset[trsize:N]
    trlabel = label[0:trsize, :]
    vallabel = label[trsize:N, :]
    if opt.weighting == 1:
        trw = CRLB_std[0:trsize, :]
        valw = CRLB_std[trsize:N, :]
    else:
        trw = None
        valw = None
    print('training data size:', trdata.size(0), 'x', trdata.size(1), 'x', trdata.size(2), 'x', trdata.size(3))
    print('training label size:', trlabel.size(0), 'x', trlabel.size(1))
    if opt.weighting == 1:
        print('training CRLB size:', trw.size(0), 'x', trw.size(1))
    return {'trainData': trdata, 'trainLabel':trlabel, 'trainWeight': trw}, {'valData': valdata, 'valLabel':vallabel, 'valWeight':valw}

def train(epoch, model, optimizer, opt, data, size):
    trainData = data['trainData']
    trainLabel = data['trainLabel']
    tCRLB_std = data['trainWeight']

    shuffle = torch.randperm(size)
    count, totalErr, err = 0, 0, 0
    x = torch.Tensor(opt.batchsize, opt.channel, opt.imWidth, opt.imHeight).cuda()
    yt = torch.Tensor(opt.batchsize, opt.labelsize).cuda()
    std = torch.Tensor(opt.batchsize, opt.labelsize).cuda()
    model.train()

    for i in tqdm(range(0, size, opt.batchsize), ncols = 100, desc = "Training (# of batches):"):
        for ii in range(opt.batchsize):
            if (i+ii) > size-1:
                break
            else:
                x[ii] = trainData[shuffle[i+ii]]
                yt[ii] = trainLabel[shuffle[i+ii]]
                if opt.weighting == 1:
                    std[ii] = tCRLB_std[shuffle[i+ii]]

        optimizer.zero_grad()
        xin = Variable(x, requires_grad = True)
        ytout = Variable(yt)
        crlb = Variable(std)
        out = model.forward(xin)

        if opt.weighting == 1:
            loss = nn.CRLBMSELoss()
            loss.cuda()
            err = loss(out, ytout, std)
        else:
            loss = nn.MSELoss()
            loss.cuda()
            err = loss(out, ytout)

        err.backward()
        optimizer.step()
        totalErr = totalErr + err.item()
        count += 1
    trainErr = totalErr/count
    print('training error:', trainErr, ' (A.U.)')
    return trainErr


def validate(epoch, model, opt, data, size):
    valData = data['valData']
    valLabel = data['valLabel']
    vCRLB_std = data['valWeight']

    shuffle_val = torch.randperm(size)
    count, totalErr, err = 0, 0, 0
    x = torch.Tensor(opt.batchsize, opt.channel, opt.imWidth, opt.imHeight).cuda()
    yt = torch.Tensor(opt.batchsize, opt.labelsize).cuda()
    std = torch.Tensor(opt.batchsize, opt.labelsize).cuda()
    model.eval()

    with torch.no_grad():
    	for i in tqdm(range(0, size, opt.batchsize), ncols = 100, desc = "Validating (# of batches):"):
            for ii in range(opt.batchsize):
                if (i+ii) > size-1:
                    break
                else:
                    x[ii] = valData[shuffle_val[i+ii]]
                    yt[ii] = valLabel[shuffle_val[i+ii]]
                    if opt.weighting == 1:
                        std[ii] = vCRLB_std[shuffle_val[i+ii]]
            xin = Variable(x)
            ytout = Variable(yt)
            out = model.forward(xin)
            crlb = Variable(std)
            if opt.weighting == 1:
                loss = nn.CRLBMSELoss()
                loss.cuda()
                err = loss(out, ytout, std)
            else:
                loss = nn.MSELoss()
                loss.cuda()
                err = loss(out, ytout)
            totalErr = totalErr + err.item()
            count += 1
    valErr = totalErr/count
    print('validation error:', valErr, ' (A.U.)')
    return valErr

def plot(opt, ep, tErr, vErr):
    fig, ax = plt.subplots(1,1)
    plt.plot(ep, tErr, color = 'black', label='Train error')
    plt.plot(ep, vErr, color = 'red', label='validation error')
    plt.xlabel('Epoch')
    plt.ylabel('Error')
    ax.set_yscale('log')
    plt.legend()
    plotname = os.path.join(opt.save, 'errorplot.png')
    plt.savefig(plotname)

def main():
    torch.manual_seed(12)

    #load arguments
    opt = options()
    if not os.path.exists(opt.save):
        os.mkdir(opt.save)

    #load data
    trsize = int(opt.trainSplit * opt.datasize)
    valsize = opt.datasize - trsize
    trainset, valset = loaddata(opt, trsize, valsize)

    #model
    files = [f for f in os.listdir(opt.save) if f.endswith(".pth")]#
    if files == []:
        model = Net(opt.channel,opt.labelsize)
        model = torch.nn.DataParallel(model).cuda()
    else:
        model = torch.load(os.path.join(opt.save, 'model{}.pth'.format(opt.modelnum)))
        opt.save = os.path.join(opt.save, str(opt.modelnum))
        os.mkdir(opt.save)

    #define optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr = opt.learningRate)

    #train
    tErr, vErr, ep = [], [], []
    for epoch in range(int(opt.modelnum), opt.maxepoch):
        print('\nepoch:', epoch)
        trainErr = train(epoch, model, optimizer, opt, trainset, trsize)
        valErr = validate(epoch, model, opt, valset, valsize)

        #plot error
        ep.append(epoch)
        tErr.append(trainErr)
        vErr.append(valErr)
        plot(opt, ep, tErr, vErr)

        #save error
        errname = os.path.join(opt.save, 'error.log')
        np.savetxt(errname, np.column_stack((tErr, vErr)))

        #save model
        if epoch % opt.savenum == 0:
            modelname =opt.save + '/model{}.pth'.format(epoch)
            torch.save(model, modelname)

        #save variable
        if epoch == 1:
            arg = [['arguments' + ': ', 'values']]
            for key in opt.__dict__:
                arg.append([key + ': ', opt.__dict__.get(key)])
            with open(os.path.join(opt.save, 'opt.txt'), 'w') as f:
                for item in arg:
                    f.write('%s\n' % item)

if __name__ == '__main__':
    main()
