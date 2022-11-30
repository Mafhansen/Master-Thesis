from numba import cuda
import numpy as np, tifffile as tf, os
from glob import glob

@cuda.jit
def convolve(O, K, R):
    '''
    The function calculates the denominator in the iterative Richard-Lucy deconvolution method.
    
    The function applies a convolution kernel to an image. Boundaries are reflective when considering indicies exceeding the boundary of the image array.

    O :     float32 device array. Image array to be convolved. 
    K :     float32 device array. Kernel used to convolve the image.
    R :     N/A device array. Temporary array used to save convolution results
    '''

    ## Define array index threads.
    tx = cuda.threadIdx.x; ty = cuda.threadIdx.y
    bx = cuda.blockIdx.x; by = cuda.blockIdx.y
    dim_x = cuda.blockDim.x; dim_y = cuda.blockDim.y

    ## Stop any out of bounds threads.
    if tx + bx*dim_x >= R.shape[0] and ty + by*dim_y >= R.shape[1]: 
        return

    ## Iterate over the kernel dimensions.
    for i in range(K.shape[0]):
        for j in range(K.shape[1]):
            cuda.syncthreads()
            R[tx+bx*dim_x,ty+by*dim_y] += K[i,j]*O[round((-1)*abs(O.shape[0]-abs(tx+bx*dim_x+i-(K.shape[0]-1)/2)))-2*round((tx+bx*dim_x+i+1) / (2*(O.shape[0]+(K.shape[0]-1)/2))),round((-1)*abs(O.shape[0]-abs(ty+by*dim_y+j-(K.shape[0]-1)/2)))-2*round((ty+by*dim_y+j+1) / (2*(O.shape[0]+(K.shape[0]-1)/2)))]
            cuda.syncthreads()

@cuda.jit
def GPU_RL(O, IMG, KT, R):
    '''
    The function compiles the inputted entries for the Richardson-Lucy method.

    The function applies a convolution kernel to an image devided element-wise with convolved iteration array. Boundaries are reflective when considering indicies exceeding the boundary of the image array.

    O :     float32 device array. Image array to be convolved. 
    IMG :   float32 device array. Original image.
    KT :    float32 device array. Transposed kernel used to convolve the image.
    R :     N/A device array. Temporary array carrying results from the convolve() function.
    '''

    ## Define array index threads.
    tx = cuda.threadIdx.x; ty = cuda.threadIdx.y
    bx = cuda.blockIdx.x; by = cuda.blockIdx.y
    dim_x = cuda.blockDim.x; dim_y = cuda.blockDim.y

    ## Stop any out of bounds threads.
    if tx + bx*dim_x >= R.shape[0] and ty + by*dim_y >= R.shape[1]: 
        return

    val = 0
    
    ## Iterate over kernel dimensions.
    cuda.syncthreads()
    for i in range(KT.shape[0]):
        for j in range(KT.shape[1]):
            cuda.syncthreads()
            if R[round((-1)*abs(O.shape[0]-abs(tx+bx*dim_x+i-(KT.shape[0]-1)/2)))-2*round((tx+bx*dim_x+i+1) / (2*(O.shape[0]+(KT.shape[0]-1)/2))),round((-1)*abs(O.shape[0]-abs(ty+by*dim_y+j-(KT.shape[0]-1)/2)))-2*round((ty+by*dim_y+j+1) / (2*(O.shape[0]+(KT.shape[0]-1)/2)))] == 0:
                val += 0
            else:
                val += KT[i,j]*IMG[round((-1)*abs(O.shape[0]-abs(tx+bx*dim_x+i-(KT.shape[0]-1)/2)))-2*round((tx+bx*dim_x+i+1) / (2*(O.shape[0]+(KT.shape[0]-1)/2))),round((-1)*abs(O.shape[0]-abs(ty+by*dim_y+j-(KT.shape[0]-1)/2)))-2*round((ty+by*dim_y+j+1) / (2*(O.shape[0]+(KT.shape[0]-1)/2)))] / R[round((-1)*abs(O.shape[0]-abs(tx+bx*dim_x+i-(KT.shape[0]-1)/2)))-2*round((tx+bx*dim_x+i+1) / (2*(O.shape[0]+(KT.shape[0]-1)/2))),round((-1)*abs(O.shape[0]-abs(ty+by*dim_y+j-(KT.shape[0]-1)/2)))-2*round((ty+by*dim_y+j+1) / (2*(O.shape[0]+(KT.shape[0]-1)/2)))]
            cuda.syncthreads()
    cuda.syncthreads()
    
    ## Update pixel values.
    O[tx+bx*dim_x,ty+by*dim_y] *= val


def RL_deconvolve(O,IMG,PSF,PSFT,file_name):

    ## Defining number of threads and blocks to use in the device (GPU).
    threadsperblock = (32,32)
    blockspergrid = ((O.shape[1] + (threadsperblock[0] - 1)) // threadsperblock[0], 
        (O.shape[2] + (threadsperblock[1] - 1)) // threadsperblock[1])

    ## Initialize iteration dependent variables.
    O_old = np.ones(shape=O.shape)*2
    tol = 1e-3
    
    ## Continue looping when the stopping criterion is not fulfilled.
    while tol < np.sum(abs(np.subtract(O,O_old,dtype=np.float64)),dtype=np.float64)/np.sum(O_old,dtype=np.float64):
        
        O_old = np.copy(O)

        ## Iterate over the number of z planes in the image.
        for i in range(O.shape[0]):
            
            ## Send image arrays to device (GPU).
            O_device = cuda.to_device(O[i,:,:])
            IMG_device = cuda.to_device(IMG[i,:,:])
            PSF_device = cuda.to_device(PSF[i,:,:])
            PSFT_device = cuda.to_device(PSFT[i,:,:])
            temp_device = cuda.device_array(shape=[O.shape[1],O.shape[2]])

            ## Execute the Richardson-Lucy method.
            convolve[blockspergrid, threadsperblock] (O_device, PSF_device, temp_device)
            GPU_RL[blockspergrid, threadsperblock] (O_device, IMG_device, PSFT_device, temp_device)

            ## Transfer updated image back to host (CPU).
            O_device.copy_to_host(O[i,:,:])

    ## Save image as .tif file.
    tf.imwrite(file_name,O)


if __name__ == '__main__':

    nCycles = 14; nChannels = 4; nTiles = 25; nZ = 9

    psf_folder = 'C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220420_CODEX_fetalheart_5pcw/Image Processing/'
    experiment_folder = 'C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220420_CODEX_fetalheart_5pcw/'

    if not os.path.isdir(psf_folder + 'deconvoluted_tiles/'):
        os.mkdir(psf_folder + 'deconvoluted_tiles/')

    for cyc in range(11,nCycles):
        if cyc + 1 < 10:
            cyc_idx = 'cyc00' + str(cyc+1) + '_reg001'
        elif cyc + 1 >= 10:
            cyc_idx = 'cyc0' + str(cyc+1) + '_reg001'

        img_folder = experiment_folder + cyc_idx+ '/'

        for ch in range(nChannels):
            psf_file = 'PSF_CH' + str(ch+1) + '.tif'
            psf_path = psf_folder + psf_file

            for tile in range(nTiles):
                if tile + 1 < 10:
                    tile_idx = '1_0000' + str(tile+1)
                elif tile + 1 >= 10:
                    tile_idx = '1_000' + str(tile+1)

                img_tile = glob(img_folder + tile_idx + '_Z00?_CH' + str(ch+1) + '.tif')

                print(img_folder + tile_idx + '_Z00?_CH' + str(ch+1) + '.tif')
                
                img_tile0 = tf.imread(img_tile[0]); img_tile1 = tf.imread(img_tile[1]); img_tile2 = tf.imread(img_tile[2])
                img_tile3 = tf.imread(img_tile[3]); img_tile4 = tf.imread(img_tile[4]); img_tile5 = tf.imread(img_tile[5])
                img_tile6 = tf.imread(img_tile[6]); img_tile7 = tf.imread(img_tile[7]); img_tile8 = tf.imread(img_tile[8])

                img = np.array([img_tile0,img_tile1,img_tile2,img_tile3,img_tile4,img_tile5,img_tile6,img_tile7,img_tile8])
                O = np.copy(img)
                
                psf = tf.imread(psf_path).astype(np.float16)
                psf = np.stack((psf,psf,psf,psf,psf,psf,psf,psf,psf))
                psfT = np.copy(psf.swapaxes(1,2), order='C')

                save_file_path = psf_folder + 'deconvoluted_tiles/' + cyc_idx + '_' + tile_idx + '_CH' + str(ch+1) + '.tif'

                RL_deconvolve(O,img,psf,psfT,save_file_path)