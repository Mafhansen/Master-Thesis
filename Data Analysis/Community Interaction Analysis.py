import pandas as pd, numpy as np, tifffile as tf, matplotlib.pyplot as plt, seaborn as sns
import networkx as nx
import copy
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

## Load images.
obj_mask_path = "C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220322_CODEX_fetalheart_7pcw/Image Processing/Membrane Objects.tiff"
obj_mask = tf.imread(obj_mask_path)

img_background_path = "C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220322_CODEX_fetalheart_7pcw/Image Processing/stitched_images/cyc013_reg001_CH4.tif"
img_background = tf.imread(img_background_path)

img_seg_path = "C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220322_CODEX_fetalheart_7pcw/Image Processing/Membrane Segmentation.tiff"
img_seg = tf.imread(img_seg_path)
img_seg[img_seg == 255] = 1

## Load exported data from file "BB200X Seurat.R".
data_path = "C:/Users/matti/Desktop/utag_7pcw_test.csv"
dat=pd.read_table(data_path,sep=',')
dat.index = dat['ObjNo']

## Remove filtered out objects from image.
obj_unique = np.unique(obj_mask)
for i in obj_unique:
    if i in dat['ObjNo']:
        continue
    elif i != 0:        
        obj_mask[obj_mask == i] = 0

## Optional: Crop images
obj_mask_crop = obj_mask#[5500:6500,6000:7000]
img_seg_crop = img_seg#[5500:6500,6000:7000]
img_background_crop = img_background#[5500:6500,6000:7000]


# dat_crop = dat[dat['x']<7000]
# dat_crop = dat_crop[dat_crop['x']>=6000]
# dat_crop = dat_crop[dat_crop['y']<6500]
# dat_crop = dat_crop[dat_crop['y']>=5500]
dat_crop = copy.deepcopy(dat)

## Create graph
G = nx.Graph()
for i in dat_crop['ObjNo']:
    G.add_node(i, pos=(dat_crop.loc[i]['y'],dat.loc[i]['x']))
wd_size = 75
for i in G.nodes.keys():
    obj_temp = obj_mask_crop[max(int(G.nodes[i]['pos'][0])-wd_size,0):min(int(G.nodes[i]['pos'][0])+wd_size, obj_mask_crop.shape[0]),
                             max(int(G.nodes[i]['pos'][1])-wd_size,0):min(int(G.nodes[i]['pos'][1])+wd_size,obj_mask_crop.shape[0])]
    seg_temp = img_seg_crop[max(int(G.nodes[i]['pos'][0])-wd_size,0):min(int(G.nodes[i]['pos'][0])+wd_size,obj_mask_crop.shape[0]),
                            max(int(G.nodes[i]['pos'][1])-wd_size,0):min(int(G.nodes[i]['pos'][1])+wd_size,obj_mask_crop.shape[0])]
    obj_seg = obj_temp * seg_temp
    indices = np.where(obj_seg == i)
    weight = {}
    for j in range(len(indices[1])):
        
        for k in range(3):
            for l in range(3):
                if indices[0][j]-1+k < 0 or indices[0][j]-1+k >= obj_temp.shape[0] or indices[1][j]-1+l < 0 or indices[1][j]-1+l >= obj_temp.shape[1]:
                    continue
                else:
                    if (obj_temp[indices[0][j]-1+k,indices[1][j]-1+l] != i) and obj_temp[indices[0][j]-1+k,indices[1][j]-1+l] != 0:
                        if obj_temp[indices[0][j]-1+k,indices[1][j]-1+l] not in weight.keys():
                            weight[obj_temp[indices[0][j]-1+k,indices[1][j]-1+l]] = 1
                        else:
                            weight[obj_temp[indices[0][j]-1+k,indices[1][j]-1+l]] += 1
                    else:
                        continue

    if weight == {}:
        continue
    else:
        for j in weight.keys():
            if j in dat_crop['ObjNo']:
                G.add_edge(i, j, weight=weight[j])
            else:
                continue
       
Clusters = list(np.unique(dat_crop['cluster']))

## Plot node degree districbutions.
deg_dict = {}
for i in Clusters:
    deg_dict[i] = []
for i in G.nodes.keys():
    deg_dict[dat_crop['cluster'][i]].append(G.degree[i])
for i in deg_dict.keys():
    sns.histplot(deg_dict[i])
fig, ax = plt.subplots(); sns.histplot(deg_dict[16],binwidth=1); ax.set_xlim(1,12)


weights_out = {}
for i in G.nodes.keys():
    weights_out[i] = []
    weight = {}
    for j in Clusters:
        weight[j] = {}
        for k in Clusters:
            weight[j][k] = 0
    for j in G._adj[i].keys():
        weight[dat_crop['cluster'][i]][dat_crop['cluster'][j]] += G._adj[i][j]['weight']
    node_neighbors = G.neighbors(i)
    for j in node_neighbors:
        for k in G._adj[j].keys():
            weight[dat_crop['cluster'][j]][dat_crop['cluster'][k]] += G._adj[j][k]['weight']
    for j in weight.keys():
        for k in weight[j].values():
            weights_out[i].append(k)

## Community analysis.
from networkx.algorithms import community
communities = community.girvan_newman(G)
communities = community.louvain_communities(G)

community_dict = {}
counter = 0
for i in communities:
    if len(i) <= 3:
        continue
    else:
        community_dict[counter] = i
        counter += 1

community_dict_out = {}
for i in community_dict.keys():
    community_dict_out[i] = {}
    for j in Clusters:
        community_dict_out[i][j] = 0
    for j in community_dict[i]:
        community_dict_out[i][dat_crop['cluster'].loc[j]] += 1

temp = np.empty(shape=(16,len(community_dict_out)))
for i in community_dict_out.keys():
    temp_max = np.array(list(community_dict_out[i].values()))
    temp_max = temp_max / np.max(temp_max)
    
    temp[:,i] = temp_max

temp = pd.DataFrame(temp)
temp.index = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]

g = sns.clustermap(temp, figsize=(9,6),cmap='cividis')
ax = g.ax_heatmap
ax.set_xlabel("")
