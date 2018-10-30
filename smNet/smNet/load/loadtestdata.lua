----------------------------------------------------------------------------------------------------
-- Script for loading test data
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
local N = opt.numSubregion

local dirRoot  = opt.testdatapath
local inputH = opt.imHeight
local inputW = opt.imWidth

local red      = '\27[31m'
local green    = '\27[32m'
local resetCol = '\27[0m'
--------------------------------------------------------------------------------
-- Initialize data structures:
--------------------------------------------------------------------------------
local testData   = {data}

print '\n\27[31m\27[4mPreparing testdata\27[0m'

-- define loader
local function forwardseq()
   -- load image
   local imageName = dirRoot ..'/testdata.mat'
   print(imageName)
   print(mattorch.load(imageName))
   local inputImg = mattorch.load(imageName).imsp[{{1,N}}]:transpose(2,3):float()

   inputImg:resize(N,1,opt.imHeight,opt.imWidth)

   -- preprocess images
   for i=1,N do
	inputImg[i][1]=inputImg[i][1]/torch.max(inputImg[i][1])
   end

   testData.data=inputImg

   testData[1]=nil
end

-----------------------------------------------------------------------------------
-- Main section
-----------------------------------------------------------------------------------
print(red .. "\nGetting test images" .. resetCol)

-- load test data
forwardseq()

print(green .. "Loaded test images!!!" .. resetCol)

collectgarbage()

print(string.format("%s# of test data :%s %d", green, resetCol, N))
print(string.format("Frame res is: %dx%d",inputH, inputW))

return testData
