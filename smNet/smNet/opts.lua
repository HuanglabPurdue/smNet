----------------------------------------------------------------------------------------------------
-- Definitions of user adjustable variables
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

local opts = {}

lapp = require 'pl.lapp'
function opts.parse(arg)
   local opt = lapp [[
   Command line options:
   Training Related:
   --learningRate          (default 1e-5)        learning rate
   --maxepoch              (default 10000)       maximum number of training iterations
   --weighting             (default 1)           0: without CRLB weighting; 1: with CRLB weighting
   --savenum               (default 1)           save model per # interations
   --plot                  (default 1)           plot training/testing error in error.log.epc
   --datapath              (default /)           training and validation dataset location
   --datasize              (default 10000)       training data size
   --xyrange               (default 2)           The range of HardTanh for training x, y (unit pixel)
   --zrange                (default 20)          The range of HardTanh for training z (unit 10micron)
   --trsize                (default 0.75)        The percentage of training data
   --startNum              (default 1)           The model number that you want to load and start to train from

   Testing Related:
   --modelNum              (default 0)           iteration number at the stop criterion
   --testdatapath          (default /)           test data location
   --numSubregion          (default 0)           test data size

   Training and Testing Related:
   --batchSize             (default 128)         batch size
   --mode                  (default z)           choose from: xy, z, Alpha, Beta
   --save                  (default /)           save trained model, error plot and test result here
   --imHeight              (default 0)           image height
   --imWidth               (default 0)           image width
   --channels              (default 1)           grayscale or RGB image
   --nGPU                  (default 4)           number of GPUs to use (the same for training and testing)
 ]]

   return opt
end


return opts
