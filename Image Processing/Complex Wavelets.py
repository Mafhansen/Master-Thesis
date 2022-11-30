from ij import IJ
from edfgui import Parameters, ExtendedDepthOfField 

# Create and setup the Parameters object

#string = getArguments();
#arg = split(string, "\t");

param = Parameters()

# List its attributes or see https://github.com/fiji-BIG/Extended_Depth_Field/blob/6213d53959cef11516a567b3b1120ace4b3e8e07/src/main/java/edfgui/Parameters.java

# Set projection method
param.edfMethod = 3    # 0=Sobel, 1=Variance, 2=RealWavelet (default), 3=ComplexWavelet

#param.varWindowSize = 5

#param.splineOrder = 5
#param.nScales = 11

#param.qualitySettings = 5
#param.topologySettings = 1
param.subBandCC = False
param.majCC = False
param.nScales = 11
param.daubechielength = 6

#print(dir(param))

# Create and set up the EDF object
EDF = ExtendedDepthOfField(IJ.getImage(),param)

EDF.process() # returns void, only display the image