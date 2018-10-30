----------------------------------------------------------------------------------------------------
-- Script for loading CRLB for training
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
local N =  opt.datasize
local trsize  = opt.trsize * N
local valsize = N - trsize
local dirRoot  = opt.datapath

local red      = '\27[31m'
local green    = '\27[32m'
local resetCol = '\27[0m'
local dataCache={
                 trainWeight,
                 valWeight
                }
----------------------------------------------------------------------
print '\n\27[31m\27[4mPreparing CRLB weight for training\27[0m'

-- define loader
local function forwardseq()
   --load weight
   local weightName=dirRoot ..'/CRLB.mat'
   local weight=mattorch.load(weightName)
   if opt.mode == 'z' then
      dataCache.trainWeight=(weight.CRLB:t()[{{1,trsize}}]:t()[{{3}}]):t():float()
      dataCache.valWeight=(weight.CRLB:t()[{{trsize+1,N}}]:t()[{{3}}]):t():float()
   else if opt.mode == 'xy' then
      dataCache.trainWeight=(weight.CRLB:t()[{{1,trsize}}]:t()[{{1,2}}]):t():float()
      dataCache.valWeight=(weight.CRLB:t()[{{trsize+1,N}}]:t()[{{1,2}}]):t():float()
   else if opt.mode == 'Alpha' then
      dataCache.trainWeight=(weight.CRLB:t()[{{1,trsize}}]:t()[{{4}}]):t():float()
      dataCache.valWeight=(weight.CRLB:t()[{{trsize+1,N}}]:t()[{{4}}]):t():float()
   else if opt.mode == 'Beta' then
      dataCache.trainWeight=(weight.CRLB:t()[{{1,trsize}}]:t()[{{5}}]):t():float()
      dataCache.valWeight=(weight.CRLB:t()[{{trsize+1,N}}]:t()[{{5}}]):t():float()
   end
   end
end
end
end

--------------------------------------------------------------
-- Main section
-----------------------------------------------------------------------------------
local loadedFromCache = false
print(red .. "\nGetting CRLB weightings" .. resetCol)

-- load CRLB for training and validation
forwardseq()
trainWeight=dataCache.trainWeight
valWeight=dataCache.valWeight

print(green .. "Loaded CRLB weightings!!!" .. resetCol)

collectgarbage()
---------------------------------------------------------------------------------

return {trainWeight, valWeight}
