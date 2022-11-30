import pandas as pd, numpy as np, tifffile as tf, seaborn as sns
import networkx as nx
from griottes import plot_2D
import copy
import random
from random import choices

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
obj_mask_crop = obj_mask
img_seg_crop = img_seg
img_background_crop = img_background

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

## Create Fig. 6A.
node_dict = {}
for i in G._adj:
    if dat_crop['cluster'][i] == 8:
        col = 'red'
    elif dat_crop['cluster'][i] == 1:
        col = 'blue'
    elif dat_crop['cluster'][i] == 2:
        col = 'green'
    elif dat_crop['cluster'][i] == 3:
        col = 'lime'
    elif dat_crop['cluster'][i] == 7:
        col = 'wheat'
    elif dat_crop['cluster'][i] == 9:
        col = 'orange'
    elif dat_crop['cluster'][i] == 11:
        col = 'steelblue'
    elif dat_crop['cluster'][i] == 14:
        col = 'slateblue'
    elif dat_crop['cluster'][i] == 17:
        col = 'turquoise'
    elif dat_crop['cluster'][i] == 20:
        col = 'yellow'

    node_dict[i] = {'label' : dat_crop['ObjNo'][i], 'cell_properties' : dat_crop['cluster'][i], 'color' : col}

nx.set_node_attributes(G, node_dict)
plot_2D(G, background_image = img_background_crop, include_weights = True, figsize = (5,5), alpha_line = 0.5, 
        scatterpoint_size = 4, legend = True, edge_color = 'w', line_factor = 0.02)

## Initialize heatmaps for Fig. 6B.
weight_heatmap = np.zeros(shape=(23,23))
count_heatmap = np.zeros(shape=(23,23))
count_dict = dict()
Clusters = np.unique(dat_crop['cluster']).tolist()
for i in Clusters:
    count_dict[i] = {}
    for j in Clusters:
        count_dict[i][j] = 0

count_dict = {}
for i in G._adj.keys():
    count_dict[i] = {}
    for j in Clusters:
        count_dict[i][j] = 0

for i in G._adj.keys():
    for j in G._adj[i].keys():
        
        count_dict[i][dat_crop['cluster'].loc[j]] += G._adj[i][j]['weight']
        
count_dict2 = {}
for i in Clusters:
    count_dict2[i] = {}
    for j in Clusters:
        count_dict2[i][j] = []

for i in count_dict.keys():
    for j in count_dict2.keys():
        count_dict2[dat_crop['cluster'].loc[i]][j].append(count_dict[i][j])

## Welch Analysis of Variance with permutation.
permutations = 1000
out = {}

for i in Clusters:
    
    W_lst = [];     obs = []
    dfn = 0;        w_k = [];       group_mean = [];        group_size = [];    cluster_idents = []
    for j in count_dict2[i].keys():
        if sum(count_dict2[i][j]) == 0:
            continue
        dfn += 1
        
        group_mean.append(np.mean(count_dict2[i][j]))
        s = sum([(x - group_mean[-1])**2 for x in count_dict2[i][j]]) / (len(count_dict2[i][j]) - 1)
        w_k.append(len(count_dict2[i][j]) / s)
        group_size.append(len(count_dict2[i][j]))
        obs += count_dict2[i][j]
        
        cluster_idents.append(j)
    
    if dfn == 1:
        out[i] = float('nan')
        continue
    
    w = sum(w_k)
    mean_total = sum([w_k[i] * group_mean[i] / w for i in range(len(w_k))])
    dfd = (dfn**2 - 1) / (3*sum([(1 - w_k[i]/w)**2 / (group_size[i] - 1) for i in range(len(w_k))]))
    SS_within = sum([1 + (2*(dfn-2) / (dfn**2 - 1)) * (1/(group_size[i] - 1)) * (1 - w_k[i]/w) for i in range(len(w_k))])
    dfn -= 1
    SS_among = sum([(1 / dfn) * w_k[i] * (group_mean[i] - mean_total)**2 for i in range(len(w_k))])
    W = SS_among / SS_within
    
    
    rnd_idx = cluster_idents * len(count_dict2[i][j])
    random.shuffle(rnd_idx)
    
    counter = 0
    while counter < permutations:
        per_dict = {};      redo = False        
        
        for j in cluster_idents:
            per_dict[j] = []
        
        for j in range(len(obs)):
            per_dict[rnd_idx[j]].append(obs[j])
    
        dfn = 0;        w_k = [];       group_mean = [];        group_size = []
        for j in per_dict.keys():
            if len(per_dict[j]) <= 1:
                redo = True
                break
            dfn += 1
            group_mean.append(np.mean(per_dict[j]))
            s = sum([(x - group_mean[-1])**2 for x in per_dict[j]]) / (len(per_dict[j]) - 1)
            w_k.append(len(per_dict[j]) / (s + 1e-30))
            group_size.append(len(per_dict[j]))
        
        if redo == True or dfn == 1:
            continue
        
        w = sum(w_k)
        
        mean_total = sum([w_k[i] * group_mean[i] / w for i in range(len(w_k))])
        dfd = (dfn**2 - 1) / (3*sum([(1 - w_k[i]/w)**2 / (group_size[i] - 1) for i in range(len(w_k))]))
        SS_within = sum([1 + (2*(dfn-2) / (dfn**2 - 1)) * (1/(group_size[i] - 1)) * (1 - w_k[i]/w) for i in range(len(w_k))])
        dfn -= 1
        SS_among = sum([(1 / dfn) * w_k[i] * (group_mean[i] - mean_total)**2 for i in range(len(w_k))])
        W_lst.append(SS_among / SS_within)
        counter += 1
    
    Statistic = [k for k in W_lst if k > W]
    p = len(Statistic) / permutations
    out[i] = p 
    print('Cell-type : ',i,', Number of interactions : ',len(obs),', F-val : ',W,', p-val : ',p)

## Non-parametric boostrap method with t-statistic computed according to Welch t-test.
permutations = 10000

p_dict1 = {}
p_dict2 = {}
per_dict = {}
obs_dict = {}
for i in count_dict2.keys():
    obs_dict[i] = {}
    per_dict[i] = {}
    p_dict1[i] = {}
    p_dict2[i] = {}
    for j in count_dict2[i].keys():
        obs_dict[i][j] = []
        per_dict[i][j] = []
        obs = []
        for k in count_dict2[i].keys():
            if j == k:
                continue
            else:
                obs += count_dict2[i][k]
         
        t_mean = np.mean(obs)
        t_std = np.std(obs)
        
        obs_mean = np.mean(count_dict2[i][j])
        obs_std = np.std(count_dict2[i][j])
        
        t = (obs_mean - t_mean) / (t_std / len(count_dict2[i][j]) + t_std / len(obs))
        
        x = [i - t_mean + (t_mean + obs_mean)/2 for i in obs]
        y = [i - obs_mean + (t_mean + obs_mean)/2 for i in count_dict2[i][j]]
        
        for per in range(permutations):
            per_obs = choices(x,k=int(len(x) * 0.2))
            obs_obs = choices(y,k=int(len(y) * 0.2))
                        
            per_mean = np.mean(per_obs)
            obs_mean = np.mean(obs_obs)
            
            per_std = np.std(per_obs)
            obs_std = np.std(obs_obs)
            
            per_dict[i][j].append((obs_mean - per_mean) / np.sqrt(per_std / len(per_obs) + obs_std / len(obs_obs)))
            
        if np.isnan(obs_dict[i][j]).any():
            p_dict1[i][j], p_dict2[i][j] = float('nan'), float('nan')
        else:
            Statistic = [b for b in per_dict[i][j] if b > t] 
            p_dict1[i][j] = len(Statistic) / len(per_dict[i][j])
            p_dict2[i][j] = t

## Plot interaction heatmaps Fig. 6B.
temp = []
for i in p_dict2.values():
    temp += list(i.values())

from matplotlib.colors import LogNorm, Normalize

heatmap_1 = np.array(temp)
heatmap_1 = heatmap_1.reshape(-1,len(Clusters))

mask = heatmap_1
mask = mask < 0.025*16
mask = mask[:-1,:]

df_heatmap = pd.DataFrame(heatmap_1, index=Clusters, columns=Clusters)
df_heatmap = df_heatmap.dropna()
ax = sns.clustermap(df_heatmap, cmap='coolwarm', vmin=-1e3,vmax=1e3,figsize=(8,8),linewidths=1,linecolor="grey")
ax.ax_col_dendrogram.set_visible(False)
ax.set_facetcolor('grey')














