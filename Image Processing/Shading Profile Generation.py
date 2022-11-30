import tifffile as tf, numpy as np
from numba import cuda
from glob import glob
import scipy.ndimage

@cuda.jit
def convolve(img, imgRes, kernel):
    
    ## Define array index threads.
    tx = cuda.threadIdx.x; ty = cuda.threadIdx.y
    bx = cuda.blockIdx.x; by = cuda.blockIdx.y
    dim_x = cuda.blockDim.x; dim_y = cuda.blockDim.y

    ## Stop any out of bounds threads.
    if tx + bx*dim_x >= imgRes.shape[0] or ty + by*dim_y >= imgRes.shape[1]: 
        return

    ## Iterate over the kernel dimensions.
    for i in range(kernel.shape[0]):
        for j in range(kernel.shape[1]):
            cuda.syncthreads()
            imgRes[tx+bx*dim_x,ty+by*dim_y] += kernel[i,j]*img[round((-1)*abs(img.shape[0]-abs(tx+bx*dim_x+i-(kernel.shape[0]-1)/2)))-2*round((tx+bx*dim_x+i+1) / (2*(img.shape[0]+(kernel.shape[0]-1)/2))),round((-1)*abs(img.shape[0]-abs(ty+by*dim_y+j-(kernel.shape[0]-1)/2)))-2*round((ty+by*dim_y+j+1) / (2*(img.shape[0]+(kernel.shape[0]-1)/2)))]
            cuda.syncthreads()


## Specification of folders with z projected tiles
data_folder_5pcw = "C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220420_CODEX_fetalheart_5pcw/Image Processing/z_projected/"
data_folder_7pcw = "C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220322_CODEX_fetalheart_7pcw/Image Processing/z_projected/"


## Retrieve file paths for background images for a specific channel
img_cyc001_5pcw = glob(data_folder_5pcw+ 'cyc001_reg001_1_000??_CH2'+ '.tif')
img_cyc001_7pcw = glob(data_folder_5pcw+ 'cyc001_reg001_1_000??_CH2'+ '.tif')
img_lst = img_cyc001_5pcw + img_cyc001_7pcw


## Average over all retrieved tiles for background
img_avg = np.empty(shape=(2048,2048))
for i in img_lst:
    img_avg += tf.imread(i) / len(img_lst)

img_blank = scipy.ndimage.median_filter(input=img_avg, size=3, mode='mirror')

## Display averaged background image
tf.imshow(img_blank)


## Retrieve file paths for marker images for a specific channel for all experiments
img_cyc002 = glob(data_folder_5pcw+ 'cyc003_reg001_1_000??_CH2'+ '.tif')
img_cyc003 = glob(data_folder_5pcw+ 'cyc004_reg001_1_000??_CH2'+ '.tif')
img_cyc004 = glob(data_folder_5pcw+ 'cyc005_reg001_1_000??_CH2'+ '.tif')
img_cyc006 = glob(data_folder_5pcw+ 'cyc006_reg001_1_000??_CH2'+ '.tif')
img_cyc007 = glob(data_folder_5pcw+ 'cyc008_reg001_1_000??_CH2'+ '.tif')
img_cyc008 = glob(data_folder_5pcw+ 'cyc009_reg001_1_000??_CH2'+ '.tif')
img_cyc009 = glob(data_folder_5pcw+ 'cyc011_reg001_1_000??_CH2'+ '.tif')
img_cyc011 = glob(data_folder_5pcw+ 'cyc012_reg001_1_000??_CH2'+ '.tif')
#img_cyc012 = glob(data_folder_5pcw+ 'cyc012_reg001_1_000??_CH3'+ '.tif')
#img_cyc013 = glob(data_folder_5pcw+ 'cyc013_reg001_1_000??_CH3'+ '.tif')
img_5pcw = img_cyc002 + img_cyc003 + img_cyc004 +img_cyc006 +img_cyc007 +img_cyc008 + img_cyc009 +img_cyc011 #+img_cyc012 +img_cyc013

img_cyc002 = glob(data_folder_7pcw+ 'cyc003_reg001_1_000??_CH2'+ '.tif')
img_cyc003 = glob(data_folder_7pcw+ 'cyc004_reg001_1_000??_CH2'+ '.tif')
img_cyc004 = glob(data_folder_7pcw+ 'cyc005_reg001_1_000??_CH2'+ '.tif')
img_cyc006 = glob(data_folder_7pcw+ 'cyc006_reg001_1_000??_CH2'+ '.tif')
img_cyc007 = glob(data_folder_7pcw+ 'cyc008_reg001_1_000??_CH2'+ '.tif')
img_cyc008 = glob(data_folder_7pcw+ 'cyc009_reg001_1_000??_CH2'+ '.tif')
img_cyc009 = glob(data_folder_7pcw+ 'cyc011_reg001_1_000??_CH2'+ '.tif')
img_cyc011 = glob(data_folder_7pcw+ 'cyc012_reg001_1_000??_CH2'+ '.tif')
#img_cyc012 = glob(data_folder_5pcw+ 'cyc012_reg001_1_000??_CH3'+ '.tif')
#img_cyc013 = glob(data_folder_5pcw+ 'cyc013_reg001_1_000??_CH3'+ '.tif')
img_7pcw = img_cyc002 + img_cyc003 + img_cyc004 +img_cyc006 +img_cyc007 +img_cyc008 + img_cyc009 +img_cyc011 #+img_cyc012 +img_cyc013

img_lst = img_5pcw + img_7pcw


## Average over all retrieved images
img_avg = np.empty(shape=(2048,2048))
for i in img_lst:
    img_avg += tf.imread(i) / len(img_lst)


tf.imshow(img_avg)      # Display averaged image
img_avg -= img_blank    # Subtract estimated background
tf.imshow(img_avg)      # Display background subtracted averaged image


## Smooth averaged image using a 501x501 pixel^2 mean kernel
N = 501
mean_kernel = np.ones(shape=(N,N))/(N*N)
img_avg_res = np.empty_like(img_avg) 

img_avg_device = cuda.to_device(img_avg)
img_avg_res_device = cuda.to_device(img_avg_res)
mean_kernel_device = cuda.to_device(mean_kernel)

threadsperblock = (32,32)
blockspergrid = ((img_avg_res.shape[0] + (threadsperblock[0] - 1)) // threadsperblock[0], 
                 (img_avg_res.shape[1] + (threadsperblock[1] - 1)) // threadsperblock[1])

convolve[blockspergrid, threadsperblock] (img_avg_device, img_avg_res_device, mean_kernel_device)
img_avg_res_device.copy_to_host(img_avg_res)

tf.imshow(img_avg_res)  # Display smoothed average image

## Save results
tf.imwrite('C:/Users/matti/Desktop/avg_CH2_res.tif', img_avg_res)
tf.imwrite('C:/Users/matti/Desktop/avg_CH2.tif', img_avg)





