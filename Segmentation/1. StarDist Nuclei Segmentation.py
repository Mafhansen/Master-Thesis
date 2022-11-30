#@ DatasetIOService io
#@ CommandService command

""" This example runs stardist on all tif files in a folder
Full list of Parameters: 
res = command.run(StarDist2D, False,
			 "input", imp, "modelChoice", "Versatile (fluorescent nuclei)",
			 "modelFile","/path/to/TF_SavedModel.zip",
			 "normalizeInput",True, "percentileBottom",1, "percentileTop",99.8,
			 "probThresh",0.5, "nmsThresh", 0.3, "outputType","Label Image",
			 "nTiles",1, "excludeBoundary",2, "verbose",1, "showCsbdeepProgress",1, "showProbAndDist",0).get();			
"""


from de.csbdresden.stardist import StarDist2D 
from glob import glob
from ij import IJ
from ij.plugin.frame import RoiManager
from ij.gui import Roi
import os

# run stardist on all tiff files in <indir> and save the label image to <outdir>
indir   = os.path.expanduser("C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220322_CODEX_fetalheart_7pcw/Image Processing/background_subtracted/New folder")
outdir  = os.path.expanduser("C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220322_CODEX_fetalheart_7pcw/Image Processing/background_subtracted/New folder")

for f in sorted(glob(os.path.join(indir,"*.tif"))):
	print "processing ", f
  
	imp = io.open(f)
  
	res = command.run(StarDist2D, False,
			"input", imp,
			"modelChoice", "Versatile (fluorescent nuclei)",
			"normalizeInput",True, "percentileBottom",1, "percentileTop",99.8,
			"probThresh",0.4, "nmsThresh", 0.3, "outputType","Both",
		    "nTiles",20, "excludeBoundary",2, "verbose",False, "showCsbdeepProgress",False, "showProbAndDist",True).get()
	label = res.getOutput("label")
	prob = res.getOutput("prob")
	rm = RoiManager.getInstance()
	
	rm.runCommand("Save", os.path.join(outdir, "ROISet.zip"))
	rm.close();
	
	io.save(prob, os.path.join(outdir,"prob."+os.path.basename(f)))
	io.save(label, os.path.join(outdir,"label."+os.path.basename(f)))
	
	
	
	
	
