----------------------------------------------------------------------
-- This script makes predictions using a DeepSEA model
-- for a minimal example run 'luajit 2_DeepSEA.lua -test_file_h5 example.h5'
----------------------------------------------------------------------


require 'paths'
require 'nn'
require 'math'
require 'torch' 
require 'hdf5'
require 'string'
require 'sys'
----------------------------------------------------------------------


print 'Processing options'

cmd = torch.CmdLine()
cmd:text()
cmd:text('ENCODE Loss Function')
cmd:text()
cmd:text('Options:')
cmd:option('-seed', 1, 'fixed input seed for repeatable experiments')
cmd:option('-threads', 8, 'number of threads')
cmd:option('-type', 'float', 'type: double | float | cuda')
cmd:option('-netPath','','overide which model to use for prediction')
cmd:option('-test_file_h5','','test file in h5 format')
cmd:option('-setDevice', 1, 'specify which gpu to use. only effective when type is cuda.')

cmd:text()
opt = cmd:parse(arg or {})

-- nb of threads and fixed seed (for repeatable experiments)
if opt.type == 'float' then
   print('Switch to float')
   torch.setdefaulttensortype('torch.FloatTensor')
elseif opt.type == 'cuda' then
   print('Switch to CUDA')
   require 'cutorch'
   require 'cunn'
   torch.setdefaulttensortype('torch.FloatTensor')
   cutorch.setDevice(opt.setDevice)
   --print(  cutorch.getDeviceProperties(cutorch.getDevice()) )
end
torch.setnumthreads(opt.threads)
torch.manualSeed(opt.seed)


batchSize = 1024

--batchSize = 400
--Read Input Data File

if opt.test_file_h5 ~= '' then
   test_file = opt.test_file_h5
   require 'hdf5'
    loaded = hdf5.open(test_file,'r')
    testData = {
       data = loaded:read('/testxdata'):all(),
    }

end


--Read DeepSEA model

if opt.netPath~='' then
   filename = opt.netPath
else
   filename =  './deepsea.cpu'
end

model = torch.load(filename)
if opt.type == 'cuda' then
   model:cuda()
else
   model:float()
end
model:evaluate() -- turn off dropout for production mode 





nfeats = 4
width = 1000
noutputs = 919
height = 1   
predSize = testData.data:size(1)

torch.setdefaulttensortype('torch.FloatTensor')


--Make predictions with DeepSEA model

alloutputs = torch.DoubleTensor(predSize, noutputs)

for t = 1,predSize,batchSize do
  print(t)
  collectgarbage()
  local inputs = torch.FloatTensor(math.min(batchSize, predSize-t+1), nfeats, width, height)


  -- create mini batch
  k = 1
  for i = t,math.min(t+batchSize-1,predSize) do
     -- load new sample
     input = testData.data[i]:float()
     inputs[k]= input
     k = k + 1
  end

  if opt.type == 'cuda' then
     inputs = inputs:cuda()
  end


  local output = model:forward(inputs):double()
  k=1
  for ii = t,math.min(t+batchSize-1,predSize) do
      alloutputs[ii] = output[k]
      k=k+1
  end
end





--Save predictions

if opt.test_file_h5~='' then
   f=hdf5.open( opt.test_file_h5 .. '.pred' .. '.h5', 'w')
   f:write('/pred',alloutputs)
   f:close()
end



