from collections import defaultdict
import gzip
import pandas as pd
import re
import allel
import csv
import glob
import json

GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame', 'details']
R_SEMICOLON = re.compile(r'\s*;\s*')    
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')
CHROMOSOMES = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']




def dataframe(filename):
    """Open an optionally gzipped GTF file and return a pandas.DataFrame.
    """
    # Each column is a list stored as a value in this dict.
    result = defaultdict(list)

    for i, line in enumerate(lines(filename)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                result[key] = [None] * i

        # Ensure this row has some value for each column.
        for key in result.keys():
            result[key].append(line.get(key, None))

    return pd.DataFrame(result)


def lines(filename):
    """Open an optionally gzipped GTF file and generate a dict for each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open

    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield parse(line)


def parse(line):
    """Parse a single GTF line and return a dict.
    """
    result = {}

    fields = line.rstrip().split('\t')

    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

    for i, info in enumerate(infos, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            result[key] = _get_value(value)

    return result


def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value

def createTSG():
    tsg={}
    for i in CHROMOSOMES:
        tsg[i]=[]
    with open('TSGs.csv', mode='r') as f:
        reader = csv.reader(f, delimiter=';')
        for row in reader:
            try:
                found = re.search('.+Ensembl:(.+?)\|.+', row[3]).group(1)
            except AttributeError:
                found = ''
            tsg['chr'+row[4]].append(found)
    return tsg

def checkPos(df, pos):
    df_ch = df.loc[str(row[0]).split("chr")[1] == new_df['seqname']]
    df_start = df_ch.loc[new_df['start'] <= pos]
    df_end = df_start.loc[df_start['end'] >= pos]
    check=df_end.size
    
    if(check >0):
        print(pos)
        print(set(df_end['gene_id'].tolist()))
#         print(df_end['start'])
#         print(df_end['end'])
        print('----')

        save=df_end['gene_id']
        
        return [check, set(save.tolist())] # to fix
    else:
        return [check, []]

def checkSupp(ch, idG, tsg):
    return (idG in tsg[ch])
    
    
df = dataframe("Homo_sapiens.GRCh38.98.chr.gtf")
#new_df = df[['seqname','feature','start','end','gene_id']].loc[df['feature'] == 'exon']
new_df = df[['seqname','feature','start','end','gene_id']].loc[((df['feature'] == 'gene')| (df['feature'] == 'exon')) & (df['gene_biotype'] == 'protein_coding')]
import os

stat1={}
stat2={}
stat3={}
for i in CHROMOSOMES:
    stat1[i] = {} #it's possible to create all the main keys (chr1, chr2...) before the operation and delete the control if(lista.get(a)==None)
    stat2[i] = {}
    stat3[i] = {}

tsg = createTSG()
new_df.start = pd.to_numeric(new_df.start, errors='coerce')
new_df.end = pd.to_numeric(new_df.end, errors='coerce')
ch = "chr1"
dataframes = {}
j=0
k = 0
for filename in glob.iglob('kidney/*.vcf'): #useful to open all vcf file in the script directory
    allel.vcf_to_csv(f"{filename}", 'example.csv') #convert from vcf to csv 
    with open('example.csv', mode='r') as f:
        reader = csv.reader(f, delimiter=',')
        i = 0
        included_cols = [0, 1, 8]      
        for row in reader:
        
            content = list(row[i] for i in included_cols)
            
            
            if row[8]=='True':
                j = j + 1
                a=row[0]
                if(ch != a):
                    new_line = 1
                    ch = a
                else:
                    new_line = 0
                
                
#                 if(new_line == 1):
#                     print("new chromosome")
                    
                check=checkPos(new_df, int(row[1]))
                
#                 if(check[0] > 0):
#                     print(check[1])
                for key in check[1]:
                    if(key not in stat2[a]):
                        stat2[a][key] = {}
                        stat2[a][key][os.path.basename(filename)] = 1
                        i = i + 1
                    
#                     stat2[a][check[1]][os.path.basename(filename)]=int(check[0])
#                     if(checkSupp(a, check[1], tsg)):
#                         stat3[a][check[1]]=int(check[0])
                    else:
                        if(os.path.basename(filename) not in stat2[a][key]):
                            stat2[a][key][os.path.basename(filename)]=1
                            i = i + 1
                        else:
                            stat2[a][key][os.path.basename(filename)]+=1  
                            i = i + 1
#                     print(pd.DataFrame(stat2))
#                     stat2[a] = {}
                
#         print('end')
        
        
# jsonlist = json.dumps(stat1) #convert "lista" in a json file to be saved
# f = open("stat1_kidney.json","w")
# f.write(jsonlist)
# f.close()

jsonex = json.dumps(stat2) #convert "lista" in a json file to be saved
f = open("stat2_kidney_updated.json","w")
f.write(jsonex)
f.close()

# jsonsupp = json.dumps(stat3) #convert "lista" in a json file to be saved
# f = open("stat3_kidney.json","w")
# f.write(jsonsupp)
# f.close()
# dataframes
for i in CHROMOSOMES:
    df = pd.DataFrame(stat2[i]).fillna(0).transpose()
    df['median'] = df.median(axis=1)
    df.to_json('stat2_kidney_updated_'+str(i)+'.json')
    # display(df)
