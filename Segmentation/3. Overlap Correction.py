import pandas as pd, tifffile as tf, numpy as np
from read_roi import read_roi_zip
from shapely.geometry import Polygon

data_dir = "C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220322_CODEX_fetalheart_7pcw/"

## Load .csv file with measurements exported from CellProfiler software.
data_path_1 = data_dir + "7pcwCH2FilterObjects.csv"
df_1 = pd.read_csv(data_path_1, sep=',')
# data_path_2 = data_dir + "5pcwCH2_2FilterObjects.csv"
# df_2 = pd.read_csv(data_path_2, sep=',')

# df_CH3_path_1 = data_dir + "5pcwCH3_1FilterObjects.csv"
# df_CH4_path_1 = data_dir + "5pcwCH4_1FilterObjects.csv"
# df_CH3_1 = pd.read_csv(df_CH3_path_1, sep=',')
# df_CH4_1 = pd.read_csv(df_CH4_path_1, sep=',')
# df_CH3_path_2 = data_dir + "5pcwCH3_2FilterObjects.csv"
# df_CH4_path_2 = data_dir + "5pcwCH4_2FilterObjects.csv"
# df_CH3_2 = pd.read_csv(df_CH3_path_2, sep=',')
# df_CH4_2 = pd.read_csv(df_CH4_path_2, sep=',')

df_CH2_cyto_path = data_dir + '7pcwCH2IdentifySecondaryObjects.csv'
df_CH3_cyto_path = data_dir + '7pcwCH3IdentifySecondaryObjects.csv'
df_CH4_cyto_path = data_dir + '7pcwCH4IdentifySecondaryObjects.csv'
df_CH2_cyto_1 = pd.read_csv(df_CH2_cyto_path, sep=',')
df_CH3_cyto_1 = pd.read_csv(df_CH3_cyto_path, sep=',')
df_CH4_cyto_1 = pd.read_csv(df_CH4_cyto_path, sep=',')

# df_CH2_cyto_path = data_dir + '5pcwCH2_2IdentifySecondaryObjects.csv'
# df_CH3_cyto_path = data_dir + '5pcwCH3_2IdentifySecondaryObjects.csv'
# df_CH4_cyto_path = data_dir + '5pcwCH4_2IdentifySecondaryObjects.csv'
# df_CH2_cyto_2 = pd.read_csv(df_CH2_cyto_path, sep=',')
# df_CH3_cyto_2 = pd.read_csv(df_CH3_cyto_path, sep=',')
# df_CH4_cyto_2 = pd.read_csv(df_CH4_cyto_path, sep=',')

## Load nuclei mask.
img_mask_path = data_dir + "Image Processing/Nuclei Objects.tiff"
img_mask = tf.imread(img_mask_path)

## Load ROIs exported from ImageJ.
roi_path = data_dir + "cyc001_reg001_CH1_ROISet.zip"
roi = read_roi_zip(roi_path)

## Retrieve coordinates of nuclei.
roi_keys = list(roi.keys())
roi_x = roi[roi_keys[0]]['x']
roi_y = roi[roi_keys[0]]['y']

## Associate each nuclei in .csv file exported from CellProfiler to ROI exported from ImageJ.
neighbors = dict()

CM_X = np.array(df_1['Location_CenterMassIntensity_X_DAPI'])
CM_Y = np.array(df_1['Location_CenterMassIntensity_Y_DAPI'])
ObjNo = np.array(df_1['ObjectNumber'])

for i in range(len(df_CH2_cyto_1)):
        
    x = (df_1['Location_CenterMassIntensity_X_DAPI'][i] - CM_X)**2
    y = (df_1['Location_CenterMassIntensity_Y_DAPI'][i] - CM_Y)**2
    dist = np.sqrt(x + y)
        
    idx = np.where(dist < 30 + df_1['AreaShape_MaximumRadius'][i])[0]
    
    ObjNo_temp = list()
    for j in idx:
        ObjNo_temp.append(img_mask[round(df_1['AreaShape_Center_Y'][j]),round(df_1['AreaShape_Center_X'][j])])
        
    neighbors[str(img_mask[round(df_1['AreaShape_Center_Y'][i]),round(df_1['AreaShape_Center_X'][i])])] = (np.array(ObjNo_temp), dist[idx])

identities = np.array(list(neighbors.keys()))

## Adjust intensities for each marker.
IntegratedIntensity_Blank1CH2 = list()
IntegratedIntensity_EmptyCyc002CH2 = list()
IntegratedIntensity_CD45 = list()
IntegratedIntensity_DCN = list()
IntegratedIntensity_CK_PAN = list()
IntegratedIntensity_ISL1 = list()
IntegratedIntensity_EmptyCyc007CH2 = list()
IntegratedIntensity_MYH7 = list()
IntegratedIntensity_CD68 = list()
IntegratedIntensity_EmptyCyc010CH2 = list()
IntegratedIntensity_ACTA2 = list()
IntegratedIntensity_JAG1 = list()
IntegratedIntensity_EmptyCyc013CH2 = list()
IntegratedIntensity_Blank2CH2 = list()

IntegratedIntensity_Blank1CH3 = list()
IntegratedIntensity_Notch1 = list()
IntegratedIntensity_MYH6 = list()
IntegratedIntensity_KI67 = list()
IntegratedIntensity_PDPN = list()
IntegratedIntensity_EmptyCyc006CH3 = list()
IntegratedIntensity_TBX20 = list()
IntegratedIntensity_CD31 = list()
IntegratedIntensity_EmptyCyc009CH3 = list()
IntegratedIntensity_CD44 = list()
IntegratedIntensity_EmptyCyc011CH3 = list()
IntegratedIntensity_CD34 = list()
IntegratedIntensity_EmptyCyc013CH3 = list()
IntegratedIntensity_Blank2CH3 = list()

IntegratedIntensity_Blank1CH4 = list()
IntegratedIntensity_EPCAM = list()
IntegratedIntensity_CNP = list()
IntegratedIntensity_GRHL2 = list()
IntegratedIntensity_ACTA1 = list()
IntegratedIntensity_EmptyCyc006CH4 = list()
IntegratedIntensity_FABP5 = list()
IntegratedIntensity_S100B = list()
IntegratedIntensity_HSP90B1 = list()
IntegratedIntensity_NPPA = list()
IntegratedIntensity_WT1 = list()
IntegratedIntensity_EmptyCyc012CH4 = list()
IntegratedIntensity_CTNNB1 = list()
IntegratedIntensity_Blank2CH4 = list()


for i in range(len(identities)):
    
    print(i)
    
    Val_CK_PAN = 0
    Val_CD68 = 0
    Val_ACTA2 = 0
    Val_CD45 = 0
    Val_DCN = 0
    Val_ISL1 = 0
    Val_MYH7 = 0
    Val_JAG1 = 0
    
    Val_CD31 = 0
    Val_TBX20 = 0
    Val_CD34 = 0
    Val_CD44 = 0
    Val_PDPN = 0
    Val_KI67 = 0
    Val_Notch1 = 0
    Val_MYH6 = 0
    
    Val_S100B = 0
    Val_ACTA1 = 0
    Val_GRHL2 = 0
    Val_EPCAM = 0
    Val_CNP = 0
    Val_CTNNB1 = 0
    Val_FABP5 = 0
    Val_NPPA = 0
    Val_HSP90B1 = 0
    Val_WT1 = 0
    
    roi_x = roi[roi_keys[int(identities[i])-1]]['x']
    roi_y = roi[roi_keys[int(identities[i])-1]]['y']
    
    # Define polygon from ROIs
    p = list()
    for j in range(len(roi_x)):
        p.append((roi_x[j], roi_y[j]))
    
    p = Polygon(p)
    
    for j in range(len(neighbors[identities[i]][0])):
        
        if int(identities[i]) == neighbors[identities[i]][0][j]:
            continue
        else:
            
            roi_x = roi[roi_keys[neighbors[identities[i]][0][j]-1]]['x']
            roi_y = roi[roi_keys[neighbors[identities[i]][0][j]-1]]['y']
            
            # Define polygon from ROIs
            q = list()
            for k in range(len(roi_x)):
                q.append((roi_x[k], roi_y[k]))
                    
            q = Polygon(q)
            
            # Find overlap between neighboring polygons
            A = p.intersection(q).area
            
            if A != 0:
            
                point = p.intersection(q).centroid.coords[0]
                identity_overlap = img_mask[round(point[1]),round(point[0])]
                
                if img_mask[round(point[1]),round(point[0])] == int(identities[i]):
                    continue
                else:
                    
                    idx = np.where(identities == str(neighbors[identities[i]][0][j]))[0][0]
                                        
                    #print(df['Intensity_IntegratedIntensity_CK_PAN'][i], 
                    #      df['Intensity_IntegratedIntensity_CK_PAN'][i] - A*df['Intensity_MeanIntensity_CK_PAN'][idx])
                    
                    Val_CK_PAN += A*df_CH2_cyto_1['Intensity_MeanIntensity_CK_PAN'][idx]
                    Val_CD68 += A*df_CH2_cyto_1['Intensity_MeanIntensity_CD68'][idx]
                    Val_ACTA2 += A*df_CH2_cyto_1['Intensity_MeanIntensity_ACTA2'][idx]
                    Val_CD45 += A*df_CH2_cyto_1['Intensity_MeanIntensity_CD45'][idx]
                    Val_DCN += A*df_CH2_cyto_1['Intensity_MeanIntensity_DCN'][idx]
                    Val_ISL1 += A*df_CH2_cyto_1['Intensity_MeanIntensity_ISL1'][idx]
                    Val_MYH7 += A*df_CH2_cyto_1['Intensity_MeanIntensity_MYH7'][idx]
                    Val_JAG1 += A*df_CH2_cyto_1['Intensity_MeanIntensity_JAG1'][idx]
                    
                    Val_CD31 += A*df_CH3_cyto_1['Intensity_MeanIntensity_CD31'][idx]
                    Val_TBX20 += A*df_CH3_cyto_1['Intensity_MeanIntensity_TBX20'][idx]
                    Val_CD34 += A*df_CH3_cyto_1['Intensity_MeanIntensity_CD34'][idx]
                    Val_CD44 += A*df_CH3_cyto_1['Intensity_MeanIntensity_CD44'][idx]
                    Val_PDPN += A*df_CH3_cyto_1['Intensity_MeanIntensity_PDPN'][idx]
                    Val_KI67 += A*df_CH3_cyto_1['Intensity_MeanIntensity_KI67'][idx]
                    Val_MYH6 += A*df_CH3_cyto_1['Intensity_MeanIntensity_MYH6'][idx]
                    Val_Notch1 += A*df_CH3_cyto_1['Intensity_MeanIntensity_Notch1'][idx]
                    
                    Val_S100B += A*df_CH4_cyto_1['Intensity_MeanIntensity_S100B'][idx]
                    Val_ACTA1 += A*df_CH4_cyto_1['Intensity_MeanIntensity_ACTA1'][idx]
                    Val_GRHL2 += A*df_CH4_cyto_1['Intensity_MeanIntensity_GRHL2'][idx]
                    Val_EPCAM += A*df_CH4_cyto_1['Intensity_MeanIntensity_EPCAM'][idx]
                    Val_CNP += A*df_CH4_cyto_1['Intensity_MeanIntensity_CNP'][idx]
                    Val_CTNNB1 += A*df_CH4_cyto_1['Intensity_MeanIntensity_CTNNB1'][idx]
                    Val_FABP5 += A*df_CH4_cyto_1['Intensity_MeanIntensity_FABP5'][idx]
                    Val_NPPA += A*df_CH4_cyto_1['Intensity_MeanIntensity_NPPA'][idx]
                    Val_HSP90B1 += A*df_CH4_cyto_1['Intensity_MeanIntensity_HSP90B1'][idx]
                    Val_WT1 += A*df_CH4_cyto_1['Intensity_MeanIntensity_WT1'][idx]
            else:
                continue
    
    IntegratedIntensity_Blank1CH2.append(df_CH2_cyto_1['Intensity_IntegratedIntensity_blank1'][i])
    IntegratedIntensity_EmptyCyc002CH2.append(df_CH2_cyto_1['Intensity_IntegratedIntensity_EmptyCyc002'][i])
    IntegratedIntensity_CD45.append(max(df_CH2_cyto_1['Intensity_IntegratedIntensity_CD45'][i]-Val_CD45,0))
    IntegratedIntensity_DCN.append(max(df_CH2_cyto_1['Intensity_IntegratedIntensity_DCN'][i]-Val_DCN,0))
    IntegratedIntensity_CK_PAN.append(max(df_CH2_cyto_1['Intensity_IntegratedIntensity_CK_PAN'][i]-Val_CK_PAN,0))
    IntegratedIntensity_ISL1.append(max(df_CH2_cyto_1['Intensity_IntegratedIntensity_ISL1'][i]-Val_ISL1,0))
    IntegratedIntensity_EmptyCyc007CH2.append(df_CH2_cyto_1['Intensity_IntegratedIntensity_EmptyCyc007'][i])
    IntegratedIntensity_MYH7.append(max(df_CH2_cyto_1['Intensity_IntegratedIntensity_MYH7'][i]-Val_MYH7,0))
    IntegratedIntensity_CD68.append(max(df_CH2_cyto_1['Intensity_IntegratedIntensity_CD68'][i]-Val_CD68,0))
    IntegratedIntensity_EmptyCyc010CH2.append(df_CH2_cyto_1['Intensity_IntegratedIntensity_EmptyCyc010'][i])
    IntegratedIntensity_ACTA2.append(max(df_CH2_cyto_1['Intensity_IntegratedIntensity_ACTA2'][i]-Val_ACTA2,0))
    IntegratedIntensity_JAG1.append(max(df_CH2_cyto_1['Intensity_IntegratedIntensity_JAG1'][i]-Val_JAG1,0))
    IntegratedIntensity_EmptyCyc013CH2.append(df_CH2_cyto_1['Intensity_IntegratedIntensity_EmptyCyc013'][i])
    IntegratedIntensity_Blank2CH2.append(df_CH2_cyto_1['Intensity_IntegratedIntensity_blank2'][i])
    
    IntegratedIntensity_Blank1CH3.append(df_CH3_cyto_1['Intensity_IntegratedIntensity_blank1'][i])
    IntegratedIntensity_Notch1.append(max(df_CH3_cyto_1['Intensity_IntegratedIntensity_Notch1'][i]-Val_Notch1,0))
    IntegratedIntensity_MYH6.append(max(df_CH3_cyto_1['Intensity_IntegratedIntensity_MYH6'][i]-Val_MYH6,0))
    IntegratedIntensity_KI67.append(max(df_CH3_cyto_1['Intensity_IntegratedIntensity_KI67'][i]-Val_KI67,0))
    IntegratedIntensity_PDPN.append(max(df_CH3_cyto_1['Intensity_IntegratedIntensity_PDPN'][i]-Val_PDPN,0))
    IntegratedIntensity_EmptyCyc006CH3.append(df_CH3_cyto_1['Intensity_IntegratedIntensity_EmptyCyc006'][i])
    IntegratedIntensity_TBX20.append(max(df_CH3_cyto_1['Intensity_IntegratedIntensity_TBX20'][i]-Val_TBX20,0))
    IntegratedIntensity_CD31.append(max(df_CH3_cyto_1['Intensity_IntegratedIntensity_CD31'][i]-Val_CD31,0))
    IntegratedIntensity_EmptyCyc009CH3.append(df_CH3_cyto_1['Intensity_IntegratedIntensity_EmptyCyc009'][i])
    IntegratedIntensity_CD44.append(max(df_CH3_cyto_1['Intensity_IntegratedIntensity_CD44'][i]-Val_CD44,0))
    IntegratedIntensity_EmptyCyc011CH3.append(df_CH3_cyto_1['Intensity_IntegratedIntensity_EmptyCyc011'][i])
    IntegratedIntensity_CD34.append(max(df_CH3_cyto_1['Intensity_IntegratedIntensity_CD34'][i]-Val_CD34,0))
    IntegratedIntensity_EmptyCyc013CH3.append(df_CH3_cyto_1['Intensity_IntegratedIntensity_EmptyCyc013'][i])
    IntegratedIntensity_Blank2CH3.append(df_CH3_cyto_1['Intensity_IntegratedIntensity_blank2'][i])
    
    IntegratedIntensity_Blank1CH4.append(df_CH4_cyto_1['Intensity_IntegratedIntensity_blank1'][i])
    IntegratedIntensity_EPCAM.append(max(df_CH4_cyto_1['Intensity_IntegratedIntensity_EPCAM'][i]-Val_EPCAM,0))
    IntegratedIntensity_CNP.append(max(df_CH4_cyto_1['Intensity_IntegratedIntensity_CNP'][i]-Val_CNP,0))
    IntegratedIntensity_GRHL2.append(max(df_CH4_cyto_1['Intensity_IntegratedIntensity_GRHL2'][i]-Val_GRHL2,0))
    IntegratedIntensity_ACTA1.append(max(df_CH4_cyto_1['Intensity_IntegratedIntensity_ACTA1'][i]-Val_ACTA1,0))
    IntegratedIntensity_EmptyCyc006CH4.append(df_CH4_cyto_1['Intensity_IntegratedIntensity_EmptyCyc006'][i])
    IntegratedIntensity_FABP5.append(max(df_CH4_cyto_1['Intensity_IntegratedIntensity_FABP5'][i]-Val_FABP5,0))
    IntegratedIntensity_S100B.append(max(df_CH4_cyto_1['Intensity_IntegratedIntensity_S100B'][i]-Val_S100B,0))
    IntegratedIntensity_HSP90B1.append(max(df_CH4_cyto_1['Intensity_IntegratedIntensity_HSP90B1'][i]-Val_HSP90B1,0))
    IntegratedIntensity_NPPA.append(max(df_CH4_cyto_1['Intensity_IntegratedIntensity_NPPA'][i]-Val_NPPA,0))
    IntegratedIntensity_WT1.append(max(df_CH4_cyto_1['Intensity_IntegratedIntensity_WT1'][i]-Val_WT1,0))
    IntegratedIntensity_EmptyCyc012CH4.append(df_CH4_cyto_1['Intensity_IntegratedIntensity_EmptyCyc012'][i])
    IntegratedIntensity_CTNNB1.append(max(df_CH4_cyto_1['Intensity_IntegratedIntensity_CTNNB1'][i]-Val_CTNNB1,0))
    IntegratedIntensity_Blank2CH4.append(df_CH4_cyto_1['Intensity_IntegratedIntensity_blank2'][i])


## Save adjusted intensities.
csv_df = pd.DataFrame({'cell_id:cell_id' : list(df_1['ObjectNumber']), 'Area:Area' : list(df_CH2_cyto_1['AreaShape_Area']), 
                       'x:x' : list(df_CH2_cyto_1['Location_CenterMassIntensity_X_DAPI']), 'y:y' : list(df_CH2_cyto_1['Location_CenterMassIntensity_Y_DAPI']),
                       'MeanRadius:MeanRadius' : list(df_CH2_cyto_1['AreaShape_MeanRadius']), 'MedianRadius' : list(df_CH2_cyto_1['AreaShape_MedianRadius']),
                       'cyc001_CH2:Blank1' : IntegratedIntensity_Blank1CH2, 'cyc001_CH3:Blank1' : IntegratedIntensity_Blank1CH3, 
                       'cyc001_CH4:Blank1' : IntegratedIntensity_Blank1CH4, 'cyc002_CH2:Empty' : IntegratedIntensity_EmptyCyc002CH2, 
                       'cyc002_CH3:Notch1' : IntegratedIntensity_Notch1, 'cyc002_CH4:EPCAM' : IntegratedIntensity_EPCAM, 
                       'cyc003_CH2:CD45' : IntegratedIntensity_CD45, 'cyc003_CH3:MYH6' : IntegratedIntensity_MYH6, 
                       'cyc003_CH4:CNP' : IntegratedIntensity_CNP, 'cyc004_CH2:DCN' : IntegratedIntensity_DCN, 
                       'cyc004_CH3:KI67' : IntegratedIntensity_KI67, 'cyc004_CH4:GRHL2' : IntegratedIntensity_GRHL2, 
                       'cyc005_CH2:CK_PAN' : IntegratedIntensity_CK_PAN, 'cyc005_CH3:PDPN' : IntegratedIntensity_PDPN, 
                       'cyc005_CH4:ACTA1' : IntegratedIntensity_ACTA1, 'cyc006_CH2:ISL1' : IntegratedIntensity_ISL1,
                       'cyc006_CH3:Empty' : IntegratedIntensity_EmptyCyc006CH3, 'cyc006_CH4:Empty' : IntegratedIntensity_EmptyCyc006CH4, 
                       'cyc007_CH2:Empty' : IntegratedIntensity_EmptyCyc007CH2, 'cyc007_CH3:TBX20' : IntegratedIntensity_TBX20, 
                       'cyc007_CH4:FABP5' : IntegratedIntensity_FABP5, 'cyc008_CH2:MYH7' : IntegratedIntensity_MYH7, 
                       'cyc008_CH3:CD31' : IntegratedIntensity_CD31, 'cyc008_CH4:S100B' : IntegratedIntensity_S100B, 
                       'cyc009_CH2:CD68' : IntegratedIntensity_CD68, 'cyc009_CH3:Empty' : IntegratedIntensity_EmptyCyc009CH3, 
                       'cyc009_CH4:HSP90B1' : IntegratedIntensity_HSP90B1, 'cyc010_CH2:Empty' : IntegratedIntensity_EmptyCyc010CH2, 
                       'cyc010_CH3:CD44' : IntegratedIntensity_CD44, 'cyc010_CH4:NPPA' : IntegratedIntensity_NPPA, 
                       'cyc011_CH2:ACTA2' : IntegratedIntensity_ACTA2, 'cyc011_CH3:Empty' : IntegratedIntensity_EmptyCyc011CH3, 
                       'cyc011_CH4:WT1' : IntegratedIntensity_WT1, 'cyc012_CH2:JAG1' : IntegratedIntensity_JAG1, 
                       'cyc012_CH3:CD34' : IntegratedIntensity_CD34, 'cyc012_CH4:Empty' : IntegratedIntensity_EmptyCyc012CH4, 
                       'cyc013_CH2:Empty' : IntegratedIntensity_EmptyCyc013CH2, 'cyc013_CH3:Empty' : IntegratedIntensity_EmptyCyc013CH3, 
                       'cyc013_CH4:CTNNB1' : IntegratedIntensity_CTNNB1, 'cyc014_CH2:Blank2' : IntegratedIntensity_Blank2CH2, 
                       'cyc014_CH3:Blank2' : IntegratedIntensity_Blank2CH3, 'cyc014_CH4:Blank2' : IntegratedIntensity_Blank2CH4})


csv_df.to_csv(data_dir+ 'CODEX_data_matrix.csv', index=False)
