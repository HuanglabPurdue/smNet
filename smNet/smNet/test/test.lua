----------------------------------------------------------------------------------------------------
-- Script for testing smNet
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
require 'nn'
require 'image'
require 'xlua'
require 'pl'
require 'cudnn'
require 'cunn'
require 'cutorch'
require 'gnuplot'
require 'mattorch'
----------------------------------------------------------------------------------------------------
-- initial settings
----------------------------------------------------------------------------------------------------
-- load options
local opts = require '../opts'
opt = opts.parse(arg)

-- fixed seed (for repeatable experiments)
torch.manualSeed(12)

local testy=torch.Tensor():cuda()
local x=torch.Tensor(opt.batchSize, opt.channels, opt.imHeight, opt.imWidth):cuda()

----------------------------------------------------------------------------------------------------
-- load test data
----------------------------------------------------------------------------------------------------
local testData=dofile('../load/loadtestdata.lua') --separate test data

----------------------------------------------------------------------------------------------------
--load model
----------------------------------------------------------------------------------------------------
local filename = paths.concat(opt.save, 'model' .. opt.modelNum .. '.net')
local model=torch.load(filename)
model:evaluate()

-- create blank tensor for saving result
local posx=torch.DoubleTensor(opt.numSubregion)
local posy=torch.DoubleTensor(opt.numSubregion)
local pos=torch.DoubleTensor(opt.numSubregion)
----------------------------------------------------------------------------------------------------
-- Main Section
----------------------------------------------------------------------------------------------------
print('\n\27[32m==> Localization using model:' .. filename .. ':\27[0m')


for k=1,testData.data:size(1),opt.batchSize do
   for kk=0,opt.batchSize-1 do
      xlua.progress(k+kk, testData.data:size(1))
      if(k+kk > testData.data:size(1)) then break end
      x[kk+1] = testData.data[k+kk]:cuda()
   end
   local testoutput=model:forward(x)
   testy=torch.cat(testy,testoutput,1)
   if(k+opt.batchSize-1 > testData.data:size(1)) then break end
end

----------------------------------------------------------------------------------------------------
-- save result
----------------------------------------------------------------------------------------------------
local savefolder = paths.concat(opt.save, 'result')
if not paths.dirp(savefolder) then paths.mkdir(savefolder) end
if opt.mode == 'z' then
   pos=testy[{{1,opt.numSubregion}}]:double()
   mattorch.save(savefolder.. '/posz.mat',pos)
   pos = nil
else if opt.mode == 'xy' then
   posx=testy[{{1,opt.numSubregion}}]:t()[1]:double()
   posy=testy[{{1,opt.numSubregion}}]:t()[2]:double()
   mattorch.save(savefolder .. '/posx.mat',posx)
   mattorch.save(savefolder .. '/posy.mat',posy)
   posx = nil
   posy = nil
else if opt.mode == 'Alpha' then
   pos=testy[{{1,opt.numSubregion}}]:double()
   mattorch.save(savefolder.. '/posa.mat',pos)
   pos = nil
else if opt.mode == 'Beta' then
   pos=testy[{{1,opt.numSubregion}}]:double()
   mattorch.save(savefolder.. '/posb.mat',pos)
   pos = nil
else
   print("\nPlease select mode between xy and z")
end
end
end
end

x = nil
testy = nil
model = nil
collectgarbage()

