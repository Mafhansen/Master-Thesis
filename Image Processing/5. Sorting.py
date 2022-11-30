import numpy as np, tifffile as tf, os

if __name__ == '__main__':
    
    ### Experiment specifications. nChannels : number of channels, nCycles : number of cycles, nTiles : number of tiles.
    nChannels = 4;  nCycles = 14;   nTiles = 20
    
    data_dir = 'C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220322_CODEX_fetalheart_7pcw/Image Processing/' 
    data_folder = data_dir+ 'background_subtracted/'
    
    
    if not os.path.isdir(data_dir + 'tiles/'):
        os.mkdir(data_dir + 'tiles/')
        
    
    for cyc in range(nCycles):
        print(cyc)
        cyc_idx = f'cyc{(cyc+1):03d}_reg001'
        
        if not os.path.isdir(data_dir+ 'tiles/'+ cyc_idx):
            os.mkdir(data_dir + 'tiles/'+ cyc_idx)
        
        for ch in range(nChannels):
            ch_idx = f'CH{ch+1}'
            if ch+1 == 1:
                if not os.path.isdir(data_dir+ 'tiles/'+ cyc_idx+ '/'+ ch_idx):
                    os.mkdir(data_dir + 'tiles/'+ cyc_idx+ '/'+ ch_idx)
                if not os.path.isdir(data_dir+ 'tiles/'+ cyc_idx+ '/'+ ch_idx+ '_background'):
                    os.mkdir(data_dir + 'tiles/'+ cyc_idx+ '/'+ ch_idx+ '_background')
            else:    
                if not os.path.isdir(data_dir+ 'tiles/'+ cyc_idx+ '/'+ ch_idx):
                    os.mkdir(data_dir + 'tiles/'+ cyc_idx+ '/'+ ch_idx)
            
            for tile in range(nTiles):
                tile_idx = f'1_{(tile+1):05d}'
                
                if ch+1 == 1:
                    img = tf.imread(data_folder+ cyc_idx+ '_'+ tile_idx+ '_'+ ch_idx+ '_background.tif')
                    tf.imwrite(data_dir + 'tiles/'+ cyc_idx+ '/'+ ch_idx+ '_background'+ '/'+ tile_idx+ '.tif',img)
                
                img = tf.imread(data_folder+ cyc_idx+ '_'+ tile_idx+ '_'+ ch_idx+ '.tif')
                tf.imwrite(data_dir + 'tiles/'+ cyc_idx+ '/'+ ch_idx+ '/'+ tile_idx+ '.tif',img)
                
        