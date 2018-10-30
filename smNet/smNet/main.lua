----------------------------------------------------------------------------------------------------
-- Main script for training smNet
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

----------------------------------------------------------------------------------------------------
-- include packages
----------------------------------------------------------------------------------------------------
require 'nn'
require 'image'
require 'xlua'
require 'pl'
require 'cudnn'
require 'cunn'
require 'cutorch'
require 'gnuplot'
require 'optim'
require 'nccl'
----------------------------------------------------------------------------------------------------
-- initial settings
----------------------------------------------------------------------------------------------------
-- load options
local opts = require 'opts'
opt = opts.parse(arg)

-- fixed seed for repeatable experiments
torch.manualSeed(12)

torch.setdefaulttensortype('torch.DoubleTensor')

-- Save options used for training smNet
local savepath
if opt.startNum > 1 then
   savepath = paths.concat(opt.save, opt.startNum)
else savepath = opt.save
end

if not paths.dirp(savepath) then paths.mkdir(savepath) end
local filename = paths.concat(savepath, 'opt.txt')
local file = io.open(filename, 'w')
for i,v in pairs(opt) do
   file:write(tostring(i)..' : '..tostring(v)..'\n')
end

-- save training and validation error
local errorLogger = optim.Logger(paths.concat(savepath, 'error.log'))

-- define learning rate
local config={
      learningRate = opt.learningRate
}

local timeRecord = torch.Tensor(opt.maxepoch, 2)


----------------------------------------------------------------------------------------------------
--setup multi-GPU
----------------------------------------------------------------------------------------------------
local function makeDataParallelTable(model, nGPU)
   if nGPU > 1 then
      local gpus = torch.range(1, nGPU):totable()
      local fastest, benchmark = cudnn.fastest, cudnn.benchmark

      local dpt = nn.DataParallelTable(1, true, true)
      :add(model, gpus)
      :threads(function()
         local cudnn = require 'cudnn'
         cudnn.fastest, cudnn.benchmark = fastest, benchmark
      end)
      dpt.gradInput = nil
      model = dpt:cuda()
   end
   return model
end
----------------------------------------------------------------------------------------------------
-- load data and weights(CRLB)
----------------------------------------------------------------------------------------------------
local data = require 'load/loaddata'
local trainData=data[1]
local valData=data[2]

if opt.weighting == 1 then
    local weight =require 'load/loadweight'
    local trainWeight=weight[1]
    local valWeight=weight[2]
end
----------------------------------------------------------------------------------------------------
-- define architecture and cost function
----------------------------------------------------------------------------------------------------
local Convolution = nn.SpatialConvolution
local SBatchNorm = nn.SpatialBatchNormalization

-- identity mapping definition
local function shortcut(nInputPlane, nOutputPlane, stride)
   if stride > 1 or nInputPlane ~= nOutputPlane then
      -- 1x1 convolution
      return nn.Sequential()
      :add(Convolution(nInputPlane, nOutputPlane, 1, 1, stride, stride,0,0):noBias())
      :add(SBatchNorm(nOutputPlane))
   else
      return nn.Identity()
   end
end

-- residual block definition
local function block(iChannels, oChannels, stride)
   local n = oChannels
   local s = nn.Sequential()
   s:add(Convolution(iChannels, n/4, 3, 3, 1, 1, 1, 1))
   s:add(SBatchNorm(n/4, 1e-3))
   s:add(nn.PReLU())
   s:add(Convolution(n/4, n/2, 3, 3, stride, stride, 1, 1))
   s:add(SBatchNorm(n/2, 1e-3))
   s:add(nn.PReLU())
   s:add(Convolution(n/2, n, 3, 3, 1, 1, 1, 1))
   s:add(SBatchNorm(n, 1e-3))

   return nn.Sequential()
         :add(nn.ConcatTable()
         :add(s)
         :add(shortcut(iChannels, n, stride)))
         :add(nn.CAddTable(true))
         :add(nn.PReLU())
end


--Start to train from certain model
local function start_from_model()
   local checkmodels = paths.dir(opt.save)
   for i, file in pairs(checkmodels) do
      if string.find(file, opt.startNum .. ".net") then print (file) end
   end
   local modelname = paths.concat(opt.save, 'model' .. opt.startNum .. '.net')
   local start_model = torch.load(modelname)
   return start_model
end


-- architecture definition
local model

if opt.startNum > 1 then
   model = start_from_model()
else
   model = nn.Sequential()

   model:add(nn.SpatialConvolution(opt.channels, 64, 7, 7, 1, 1, 3, 3))
   model:add(SBatchNorm(64, 1e-3))
   model:add(nn.PReLU())

   model:add(nn.SpatialConvolution(64, 128, 5, 5, 1, 1, 2, 2))
   model:add(SBatchNorm(128, 1e-3))
   model:add(nn.PReLU())

   model:add(block(128,128,1))
   model:add(block(128,128,1))
   model:add(block(128,128,1))

   model:add(block(128,256,1))

   model:add(block(256,256,1))
   model:add(block(256,256,1))
   model:add(block(256,256,1))


   model:add(nn.SpatialConvolution(256, 128, 1, 1, 1, 1, 0, 0))
   model:add(SBatchNorm(128, 1e-3))
   model:add(nn.PReLU())

   model:add(nn.SpatialConvolution(128, 64, 1, 1, 1, 1, 0, 0))
   model:add(SBatchNorm(64, 1e-3))
   model:add(nn.PReLU())


   model:add(nn.SpatialConvolution(64, 1, 1, 1, 1, 1, 0, 0))
   model:add(SBatchNorm(1, 1e-3))
   model:add(nn.PReLU())


   model:add(nn.View(opt.imHeight*opt.imWidth))
   model:add(nn.Linear(opt.imHeight*opt.imWidth,10))
   model:add(nn.PReLU())

   if opt.mode == 'z' then
      model:add(nn.Linear(10,1))
      model:add(nn.HardTanh(-opt.zrange,opt.zrange))
   else if opt.mode == 'xy' then
      model:add(nn.Linear(10,2))
      model:add(nn.HardTanh(-opt.xyrange,opt.xyrange))
   else if opt.mode == 'Alpha' then
      model:add(nn.Linear(10,1))
      model:add(nn.HardTanh(0, 90))
   else if opt.mode == 'Beta' then
      model:add(nn.Linear(10,1))
      else
      print("\27[31m \nError: \n\tPlease select training mode between xy and z")
   end
   end
   end
   end

   if opt.nGPU > 1 then
      model=makeDataParallelTable(model, opt.nGPU)
   else
      model:cuda()
   end
end

local labelsize
if opt.mode == 'xy' then
   labelsize = 2
else
   labelsize = 1
end


-- cost function definition
local loss = nn.MSECriterion()
loss:cuda()

----------------------------------------------------------------------------------------------------
-- setups for training
----------------------------------------------------------------------------------------------------

-- created blank tensor for store batch
local x=torch.Tensor(opt.batchSize, opt.channels, opt.imHeight, opt.imWidth):cuda()
local yt=torch.Tensor(opt.batchSize, labelsize):cuda()
local crlb=torch.Tensor(opt.batchSize, labelsize):cuda()

-- get parameters from smNet
local w, dE_dw = model:getParameters()


-- define trainer
local eval_E = function(w)
   model:zeroGradParameters()
   local y = model:forward(x)
   local dE_dy

   if opt.weighting == 1 then
      -- calculate CRLB weighted error
      err = loss:forward(torch.cdiv(y, crlb), torch.cdiv(yt, crlb))
      -- calculate CRLB weighted derivative
      dE_dy = loss:backward(torch.cdiv(y, torch.cmul(crlb,crlb)), torch.cdiv(yt, torch.cmul(crlb,crlb)))
   else
      err = loss:forward(y,yt)
      dE_dy = loss:backward(y,yt)
   end
   model:backward(x, dE_dy)
   return err, dE_dw
end

----------------------------------------------------------------------------------------------------
-- Main Section
---------------------------------------------------------------------------------------------------

for epoch = opt.startNum, opt.maxepoch do
   print("\n\27[31m\27[4mIteration # \27[0m" .. epoch)

   local loop_start = sys.clock()
   -- shuffle data order
   local shuffle = torch.randperm(trainData.data:size(1))
   local shuffle_val = torch.randperm(valData.data:size(1))

   -- training
   model:training()
   local totalErr = 0
   local err
   print('\n\27[32m==> Training:\27[0m')

   local trtime=sys.clock()

   for i=1, trainData.data:size(1), opt.batchSize do
       for ii=0, opt.batchSize-1 do
           if (i+ii)>trainData.data:size(1) then break end
           xlua.progress(i+ii, trainData.data:size(1))
           x[ii+1] = trainData.data[shuffle[i+ii]]
           yt[ii+1] = trainData.label[shuffle[i+ii]]

           if opt.weighting == 1 then
               crlb[ii+1] = trainWeight[shuffle[i+ii]]
           end

        end
       local _, errt = optim.adam(eval_E, w, config)
       totalErr=totalErr+errt[1]
   end

   trtime=sys.clock()-trtime
   local trainErr = totalErr/(trainData.data:size(1)/opt.batchSize)
   print("Training Error: " .. trainErr .. " (a.u.)   Training Time: " .. trtime .. ' (s)')

   -- validating
   model:evaluate()
   totalErr=0
   print('\n\27[32m==> Validating:\27[0m')
   for j=1, valData.data:size(1), opt.batchSize do
       for jj=0, opt.batchSize-1 do
           if (j+jj)>valData.data:size(1) then break end
           xlua.progress(j+jj, valData.data:size(1))
           x[jj+1] = valData.data[shuffle_val[j+jj]]:cuda()
           yt[jj+1] = valData.label[shuffle_val[j+jj]]:cuda()
           if opt.weighting == 1 then
               crlb[jj+1] = valWeight[shuffle_val[j+jj]]
           end
       end
       local valy=model:forward(x)
       local err_val
       if opt.weighting == 1 then
           err_val=loss:updateOutput(torch.cdiv(valy,crlb),torch.cdiv(yt,crlb));
       else
           err_val=loss:updateOutput(valy,yt);
       end
       totalErr=totalErr+err_val
   end
   local valErr=totalErr/(valData.data:size(1)/opt.batchSize)
   print("Validating Error: " .. valErr .. " (a.u.)")

   --save model, error and plot
   if epoch % opt.savenum == 0 then
       local filename = paths.concat(savepath, 'model' .. epoch .. '.net')
       torch.save(filename, model)
   end

   timeRecord[epoch][1]=sys.clock()-loop_start
   timeRecord[epoch][2]=trtime

   mattorch.save(savepath .. '/time.mat', timeRecord)

   errorLogger:add{
                   ['Training error'] = trainErr,
                   ['Validating error'] = valErr
                   }

   if opt.plot then
      errorLogger:style{['Training error'] = '-',
      ['Validating error'] = '-'}
      errorLogger:display(false)
      errorLogger:plot()
      errorLogger:setlogscale(true)
   end
end
