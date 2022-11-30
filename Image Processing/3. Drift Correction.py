'''
Aligns a stack of images that have been projected into the focal plane of best focus as estimated by complex wavelets.
The code was designed for a NIVIDIA GeForce GTX 1650 with Max-Q design GPU.
'''


import numpy as np, tifffile as tf, csv
from numba import cuda

# Define max drift.
x_drift = 30;    y_drift = 10

@cuda.jit
def Drift_Correction(A, B, C):
    
    tx = cuda.threadIdx.x;      ty = cuda.threadIdx.y
    bx = cuda.blockIdx.x;       by = cuda.blockIdx.y
    dim_x = cuda.blockDim.x;    dim_y = cuda.blockDim.y
        
    if tx + bx*dim_x >= A.shape[0] or ty + by*dim_y >= A.shape[1]:
        return
        
    for x in range(-x_drift,x_drift+1):
        for y in range(-y_drift,y_drift+1):
            
            if tx+bx*dim_x + x < 0 or tx+bx*dim_x + x >= A.shape[0] or ty + by*dim_y + y < 0 or ty + by*dim_y + y >= A.shape[1]:
                continue
            
            cuda.syncthreads()
            
            cuda.atomic.add(C,(x+x_drift,y+y_drift),((B[tx+bx*dim_x+x,ty+by*dim_y+y]-A[tx+bx*dim_x,ty+by*dim_y])**2)/((A.shape[0]-abs(x))*(A.shape[1]-abs(y))))
                        
            cuda.syncthreads()
            
            
    
if __name__ == '__main__':
    
    nCycles = 14;   nTiles = 20
    
    data_dir = 'C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220322_CODEX_fetalheart_7pcw/Image Processing/'
    data_folder = data_dir+ 'CH1_background_subtracted/'
    
    C = np.zeros(shape=(2*x_drift+1,2*y_drift+1),dtype=np.float64) # Temporary storage for best match.
    C_store = np.zeros(shape=(nTiles,nCycles,2*x_drift+1,2*y_drift+1),dtype=np.float64) # Storage for best match.
    
    for tile in range(nTiles):
        
        tile_idx = f'1_{(tile+1):05d}'
        
        # Tile from first cycle and DAPI channel.
        A = tf.imread(data_folder+ 'cyc001_reg001_'+ tile_idx+ '_CH1.tif').astype(np.float64)
               
        for i in drift_layers:            
            for cyc in range(nCycles):
                
                cyc_idx = f'cyc{(cyc+1):03d}_reg001_'
                
                # Tile from subsequent cycles to be corrected.
                B = tf.imread(data_folder+ cyc_idx+ tile_idx+ '_CH1.tif').astype(np.float64)
                
                # Reset temporary storage for best match.
                C *= 0
                
                # Define parameters for GPU.
                threadsperblock = (32,32)
                blockspergrid = ((A.shape[0] + (threadsperblock[0] - 1)) // threadsperblock[0], 
                    (A.shape[1] + (threadsperblock[1] - 1)) // threadsperblock[1])
                
                # Send matrices to GPU. Corners are exluded due non-consistent intensity variations resulting from processing.
                A_device = cuda.to_device(np.copy(A[10:A.shape[0]-10,10:A.shape[1]-10],order='C'))
                B_device = cuda.to_device(np.copy(B[10:B.shape[0]-10,10:B.shape[1]-10],order='C'))
                C_device = cuda.to_device(C)
                
                Drift_Correction[blockspergrid, threadsperblock](A_device,B_device,C_device)
                    
                C_device.copy_to_host(C)
                
                C_store[tile,cyc,:,:] = np.add(C_store[tile,cyc,:,:],C)

                print(f'tile {tile+1} cycle {cyc+1}, layer {i} : ',np.where(C_store[tile,cyc,:,:] == np.min(C_store[tile,cyc,:,:])), np.min(C_store[tile,cyc,:,:]))
        
        # Save parsing producing best overlap.
        if tile == 0:
            lst = list()
            for cyc in range(nCycles):
                lst.append([[np.where(C_store[tile,cyc,:,:] == np.min(C_store[tile,cyc,:,:]))[0][0]-x_drift,np.where(C_store[tile,cyc,:,:] == np.min(C_store[tile,cyc,:,:]))[1][0]-y_drift]])
        else:
            for cyc in range(nCycles):
                lst[cyc].append([np.where(C_store[tile,cyc,:,:] == np.min(C_store[tile,cyc,:,:]))[0][0]-x_drift,np.where(C_store[tile,cyc,:,:] == np.min(C_store[tile,cyc,:,:]))[1][0]-y_drift])

    # Write parsing producing best overlap to .csv file.
    with open(data_dir+ 'drift_correction.csv', 'w', newline='') as file:
        wr = csv.writer(file)
        
        wr.writerows(lst)
    
    
    
    
    
    
    
    
    
    
    
    
    