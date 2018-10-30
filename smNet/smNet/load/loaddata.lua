----------------------------------------------------------------------------------------------------
-- Script for loading training and validation data
--
--(C) Copyright 2017                The Huang Lab
--
--    All rights reserved           Weldon School of Biomedical Engineering
--                                  Purdue University
--                                  West Lafayette, Indiana
--                                  USA
--
--    Author: Peiyi Zhang, December 2017
----------------------------------------------------------------------------------------------------
require 'image'
require 'xlua'
require 'cudnn'
require 'cunn'
require 'cutorch'
require 'mattorch'
----------------------------------------------------------------------
local N = opt.datasize
local trsize  = opt.trsize * N
local valsize = N - trsize
local dirRoot  = opt.datapath

local red      = '\27[31m'
local green    = '\27[32m'
local resetCol = '\27[0m'
local dataCache={
                 trainData={data,labels},
                 valData={data,labels}
                }
----------------------------------------------------------------------
print '\n\27[31m\27[4mPreparing data for training\27[0m'

-- define loader
local function forwardseq()

   -- load images
   local imageName = dirRoot ..'/data.mat'
   local inputImg = mattorch.load(imageName).imsp:transpose(2,3):float()
   inputImg:resize(N,opt.channels,opt.imHeight,opt.imWidth)


   -- preprocess input images
   for i=1,N do
   	inputImg[i][1]=inputImg[i][1]/torch.max(inputImg[i][1])
   end


   dataCache.trainData.data=inputImg[{{1,trsize}}]:float()
   dataCache.valData.data=inputImg[{{trsize+1,N}}]:float()

   dataCache.trainData[1]=nil
   dataCache.valData[1]=nil


   --load labels
   local labelName=dirRoot ..'/label.mat'
   local label=mattorch.load(labelName)
   if opt.mode == 'z' then
      dataCache.trainData.label=label.label[{{3}, {1,trsize}}]:t():float()
      dataCache.valData.label=label.label[{{3}, {trsize+1,N}}]:t():float()
   else if opt.mode == 'xy' then
      dataCache.trainData.label=label.label[{{1, 2}, {1,trsize}}]:t():float()
      dataCache.valData.label=label.label[{{1, 2}, {trsize+1,N}}]:t():float()
   else if opt.mode == 'Alpha' then
      dataCache.trainData.label=label.label[{{4},{1,trsize}}]:t():float()
      dataCache.valData.label=label.label[{{4},{trsize+1,N}}]:t():float()
   else if opt.mode == 'Beta' then
      dataCache.trainData.label=label.label[{{5},{1,trsize}}]:t():float()
      dataCache.valData.label=label.label[{{5},{trsize+1,N}}]:t():float()
   end
   end

end
end

end

-----------------------------------------------------------------------------------
-- Main section
-----------------------------------------------------------------------------------
print(red .. "\nGetting input images and labels" .. resetCol)

-- load training and validation data
forwardseq()

trainData=dataCache.trainData
valData=dataCache.valData

print(green .. "Loaded input images and labels!!!" .. resetCol)

collectgarbage()

print(string.format("%s# of training   data :%s %d", green, resetCol, trsize))
print(string.format("%s# of validation data :%s %d", green, resetCol, valsize))
print(string.format("Frame res is: %dx%d",opt.imHeight, opt.imWidth))
---------------------------------------------------------------------------------

return {trainData, valData}
