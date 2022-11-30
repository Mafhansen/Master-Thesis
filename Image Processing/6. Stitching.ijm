// The macro assumes that images are saved in the format:
//
// output_dir > tiles > cycle_idx > channel > tile_idx
//
// Images in the tiles folder are all deconvoluted, background corrected, drift corrected, and had their z-stack projected into one z-plane.
// In the tiles folder, there is one folder for each cycle denoted by cycle_idx. Within each such folder are folders dividing the images into separate channels.
// Tile_idx are of the format 1_{iiiii} where {iiiii} is between 00000 to 99999.


macro "stitching" {
	//arg_str = split(getArgument(), "\t");
	//input_dir = arg_str[0];	output_dir = arg_str[1];
	//nCycles = arg_str[2];	nChannels = arg_str[3];	nTiles = arg_str[4];	nZ = arg_str[5];	grid_x = arg_str[6];	grid_y = arg_str[7];
	
	//input_dir = getDirectory("Choose input directory.");
	output_dir = getDirectory("Choose output directory.");
	
	nCycles = 14; //getNumber("Enter number of cycles: ", 1);
	nChannels = 4; //getNumber("Enter number of channels: ", 1);
	nTiles = 20; //getNumber("Enter number of tiles: ", 1);
	//nZ = getNumber("Enter number of z-planes: ", 1);
	grid_x = 5; //getNumber("Enter number of tiles in x direction: ", 1);
	grid_y = 4; //getNumber("Enter number of tiles in y direction: ", 1);
	
	
	setBatchMode(true);
	
	sep = File.separator;

	
	//Generate tile configuration file.

	parDAPI = "type=[Grid: snake by rows] order=[Right & Down                ] grid_size_x="+ grid_x+ " grid_size_y="+ grid_y+ " tile_overlap=20.0 first_file_index_i=1 directory=["+ output_dir+ "tiles"+ sep+ "cyc001_reg001"+ sep+ "CH1"+ sep+ "] file_names=1_{iiiii}.tif output_textfile_name=TileConfiguration.txt fusion_method=[Linear Blending] regression_threshold=0.70 max/avd_displacement_threshold=2.5 absolute_displacement_threshold=5 compute_overlap computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]";
	run("Grid/Collection stitching", parDAPI);
	rename("cyc001_reg001_CH1");	
	
	close("*");	
	
	
	//Stitch rest of images based on generated tile configuration file.
	
	File.makeDirectory(output_dir+ sep+ "stitched_images"+ sep);
	
	file_temp = output_dir+ sep+ "tiles"+ sep+ "cyc001_reg001"+ sep+ "CH1"+ sep+ "TileConfiguration.registered.txt";
	file_source = output_dir+ sep+ "stitched_images"+ sep+ "TileConfiguration.registered.txt";
	File.copy(file_temp, file_source);
	
	for (cyc=1; cyc<=nCycles; cyc++) {
		if (cyc<10) {
			cyc_idx = "cyc00"+cyc+"_reg001";
		} else if (cyc >= 10 && cyc < 100) {
			cyc_idx = "cyc0"+cyc+"_reg001";
		} else {
			cyc_idx = "cyc"+cyc+"_reg001";
		}
		for (ch=1; ch<=nChannels; ch++) {
			file_copy = output_dir+ sep+ "tiles"+ sep+ cyc_idx+ sep+ "CH"+ ch+ sep+ "TileConfiguration.registered.txt";
			File.copy(file_source, file_copy);
			
			parameters = "type=[Positions from file] order=[Defined by TileConfiguration] directory=["+ output_dir+ "tiles"+ sep+ cyc_idx+ sep+ "CH"+ ch+ sep+ "] layout_file=[TileConfiguration.registered.txt] fusion_method=[Linear Blending] computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]";
			run("Grid/Collection stitching", parameters);
			rename(cyc_idx+ "_CH"+ ch);
			
			saveAs("Tiff", output_dir+ "stitched_images"+ sep+ cyc_idx+ "_CH"+ ch);
			
			close(cyc_idx+ "_CH"+ ch+ ".tif");
			
		}
	}
	setBatchMode(false);
}
