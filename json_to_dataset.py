import json
import pandas as pd
import numpy as np

frames = []
chr_range = np.arange(1, 23).tolist()
chr_range.append('M')
chr_range.append('X')
chr_range.append('Y')
chr_range
path = 'C:\\Users\\io\\Pictures\\res\\'

for i in chr_range:
    with open(path + 'stat2_kidney_snv_updated_chr'+str(i)+'.json') as json_file:
        data = json.load(json_file)
        df1 = pd.DataFrame.from_dict(data, orient='columns')
        df1['target'] = 'kidney'
    frames.append(df1)



for i in chr_range:
    with open(path + 'stat2_kidney_snv_updated_chr'+str(i)+'.json') as json_file:
        data = json.load(json_file)
        df1 = pd.DataFrame.from_dict(data, orient='columns')
        df1['target'] = 'lung'
    frames.append(df1)

result = pd.concat(frames)




### PRINT SIZES #####

# for i in chr_range:
#     if(str(i).isnumeric()):
#         print("kidney chr"+str(i)+": "+str(frames[i-1].shape))
#         print("lung chr"+str(23-i)+": "+str(frames[2*(i-1)].shape)+'\r\n')
# print("kidney chrM:"+str(frames[22].shape))
# print("kidney chrX:"+str(frames[23].shape))
# print("kidney chrY:"+str(frames[24].shape))

# print("lung chrM:"+str(frames[47].shape))
# print("lung chrX:"+str(frames[48].shape))
# print("lung chrY:"+str(frames[49].shape))
        
    
    
    
print("total:" + str(result.shape))
result = result.fillna(0)


########## RENAME PATIENTS NAMES TO 1, 2, 3, ..
# renamed = result.rename(columns={x:y for x,y in zip(result.columns,range(0,len(result.columns)-1))})



##### TO ADD GENE ID AS COLUMN 
# result.insert(0,'Gene ID',result.axes[0].tolist())
# result.reset_index(drop=True, inplace=True)