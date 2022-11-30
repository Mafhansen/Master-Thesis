import numpy as np, tifffile as tf, pandas as pd, os, scipy.ndimage

if __name__ == '__main__':
    
    nChannels = 4;  nCycles = 14;   nTiles = 20;    nZ = 9
    
    data_dir = 'C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220322_CODEX_fetalheart_7pcw/Image Processing/' 
    data_folder = data_dir+ 'z_projected/'
    
    align_folder = data_dir+ 'CH1_background_subtracted/'
    
    # Load shading profiles for each channel.
    img_adj_vignette_CH2 = tf.imread('C:/users/matti/Desktop/avg_CH2_res.tif')
    img_adj_vignette_CH2 = (img_adj_vignette_CH2 - 0) / (np.max(img_adj_vignette_CH2) - 0)
    img_adj_vignette_CH3 = tf.imread('C:/users/matti/Desktop/avg_CH3_res.tif')
    img_adj_vignette_CH3 = (img_adj_vignette_CH3 - 0) / (np.max(img_adj_vignette_CH3) - 0)
    img_adj_vignette_CH4 = tf.imread('C:/users/matti/Desktop/avg_CH4_res.tif')
    img_adj_vignette_CH4 = (img_adj_vignette_CH4 - 0) / (np.max(img_adj_vignette_CH4) - 0)
    
    # Load .csv file with drift corrected overlaps.
    df = pd.read_csv(data_dir+ 'drift_correction.csv', delimiter=',', header=None)
    
    max_val = dict()
    
    if not os.path.isdir(data_dir + 'background_subtracted/'):
        os.mkdir(data_dir + 'background_subtracted/')
    
    
    for tile in range(nTiles):
        
        print(tile)
        
        tile_idx = f'1_{(tile+1):05d}_'
        
        df_shift_x = list();    df_shift_y = list();
        
        # Retrieve drift corrected overlaps.
        for cyc in range(nCycles):
            df_shift_x.append(int(df[tile][cyc].strip("[]").split(', ')[0]))
            df_shift_y.append(int(df[tile][cyc].strip("[]").split(', ')[1]))

        # Define most extreme shifts in x and y directions across all cycles for a tile.
        shift_x_max = max(df_shift_x);  shift_y_max = max(df_shift_y)
        shift_x_min = min(df_shift_x);  shift_y_min = min(df_shift_y)
                 
        x_slice = list();  y_slice = list()
        
        # Define indeces to use for a tile across all cycles.
        for cyc in range(nCycles):
            x_slice.append([df_shift_x[cyc]- shift_x_min, 2048+ shift_x_min- shift_x_max+ df_shift_x[cyc]])
            y_slice.append([df_shift_y[cyc]- shift_y_min, 2048+ shift_y_min- shift_y_max+ df_shift_y[cyc]])
        
        
        for ch in range(nChannels):
            
            ch_idx = f'CH{ch+1}.tif'
            
            # Load blank cycle.
            img_blank = tf.imread(data_folder+ 'cyc001_reg001_'+ tile_idx+ ch_idx)
            img_blank = img_blank[x_slice[0][0]:x_slice[0][1],y_slice[0][0]:y_slice[0][1]]
            img_blank = scipy.ndimage.median_filter(input=img_blank, size=3, mode='mirror')
            
            if ch+1 == 2:
                img_adj = img_adj_vignette_CH2[x_slice[0][0]:x_slice[0][1],y_slice[0][0]:y_slice[0][1]]
            elif ch+1 == 3:
                img_adj = img_adj_vignette_CH3[x_slice[0][0]:x_slice[0][1],y_slice[0][0]:y_slice[0][1]]
            elif ch+1 == 4:
                img_adj = img_adj_vignette_CH4[x_slice[0][0]:x_slice[0][1],y_slice[0][0]:y_slice[0][1]]
            
            for cyc in range(nCycles):
                
                cyc_idx = f'cyc{(cyc+1):03d}_reg001_'
                
                # Perform background substraction and shading correction.
                if (ch+1) == 1: # DAPI channel.
                    
                    img = tf.imread(align_folder+ cyc_idx+ tile_idx+ ch_idx)
                    img = img[x_slice[cyc][0]:x_slice[cyc][1],y_slice[cyc][0]:y_slice[cyc][1]]
                    tf.imwrite(data_dir+ 'background_subtracted/'+ cyc_idx+ tile_idx+ f'CH{ch+1}_background.tif', img)
                    
                    img = tf.imread(data_folder+ cyc_idx+ tile_idx+ ch_idx)
                    img = img[x_slice[cyc][0]:x_slice[cyc][1],y_slice[cyc][0]:y_slice[cyc][1]]
                    
                elif cyc+1 == 1 or cyc+1 == nCycles: # Blank cycles.
                    
                    img = tf.imread(data_folder+ cyc_idx+ tile_idx+ ch_idx)
                    img = img[x_slice[cyc][0]:x_slice[cyc][1],y_slice[cyc][0]:y_slice[cyc][1]]
                
                else:
                    
                    img = tf.imread(data_folder+ cyc_idx+ tile_idx+ ch_idx)
                    img = img[x_slice[cyc][0]:x_slice[cyc][1],y_slice[cyc][0]:y_slice[cyc][1]]
                    
                    img -= img_blank
                    img = img / img_adj
                
                # Clip intensities exceeding the 16-bit image depth.
                img[img < 0] = 0
                img[img > 65535] = 65535
                
                # Save background substracted and shading corrected images.
                tf.imwrite(data_dir+ 'background_subtracted/'+ cyc_idx+ tile_idx+ ch_idx,img[15:img.shape[0]-15,15:img.shape[1]-15])
                
            
            
            
            
            
            
            
            
            
    