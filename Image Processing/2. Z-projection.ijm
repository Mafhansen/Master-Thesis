
macro "projection" {
	
	//input_dir = getDirectory("Choose input directory.");
	output_dir = getDirectory("Choose output directory.");
	
	// Define expeiment parameters.
	nCycles = getNumber("Enter number of cycles: ", 1);
	nChannels = getNumber("Enter number of channels: ", 1);
	nTiles = getNumber("Enter number of tiles: ", 1);
	nZ = getNumber("Enter number of focal planes: ", 1);
	
	
	setBatchMode(true);
	
	sep = File.separator;
	
	File.makeDirectory(output_dir+"z_projected"+sep);

	for (cyc=1; cyc<=nCycles; cyc++) {
		if (cyc<10) {
			cyc_idx = "cyc00"+cyc+"_reg001";
		} else {
			cyc_idx = "cyc0"+cyc+"_reg001";
		}
		
		for (ch=1; ch<=nChannels; ch++) {
			for (tile=1; tile<=nTiles; tile++) {
				if (tile<10) {
					tile_idx = "1_0000"+tile;
				} else {
					tile_idx = "1_000"+tile;
				}
				
				// Retrieve and duplicate tiles.		
				open(output_dir+ sep+ "deconvoluted_tiles"+ sep+ cyc_idx+ "_"+ tile_idx+ "_CH"+ ch+ ".tif");
				run("Duplicate...", "title="+ cyc_idx+ "_"+ tile_idx+ "_CH"+ ch+ "_duplicate.tif duplicate");
				close(cyc_idx+ "_"+ tile_idx+ "_CH"+ ch+ ".tif");
				
				// Remove bottom z layer due to observed intensity skew in acquired images.
				selectWindow(cyc_idx+ "_"+ tile_idx+ "_CH"+ ch+ "_duplicate.tif");
				setSlice(1);
				run("Delete Slice");
				
				// Run the complex wavelet algorithm.
				runMacro("Complex_Wavelets.py");
				
				rename(cyc_idx+ "_"+ tile_idx+ "_CH"+ ch+ "_focus.tif");
				close(cyc_idx+ "_"+ tile_idx+ "_CH"+ ch+ "_duplicate.tif");
				selectWindow(cyc_idx+ "_"+ tile_idx+ "_CH"+ ch+ "_focus.tif");
				
				// Save projected image.
				saveAs("Tiff", output_dir+ sep+ "z_projected"+ sep+ cyc_idx+ "_"+ tile_idx+ "_CH"+ ch);
				
				close("*");				
				run("Collect Garbage");
			}
		}
	}
	setBatchMode(false);
}
