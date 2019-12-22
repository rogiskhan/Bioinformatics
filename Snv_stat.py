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
        save=df_end['gene_id']
        return [check, set(save.tolist())]
    else:
        return [check, []]

def checkSupp(ch, idG, tsg):
    return (idG in tsg[ch])

df = dataframe("Homo_sapiens.GRCh38.98.chr.gtf")
new_df = df[['seqname','feature','start','end','gene_id']].loc[((df['feature'] == 'gene')| (df['feature'] == 'exon')) & (df['gene_biotype'] == 'protein_coding')]

import os

statg={}
stat1={}
stat2={}
stat3={}
for i in CHROMOSOMES:
    statg[i] = 0
    stat1[i] = {} 
    stat2[i] = {}
    stat3[i] = {}

tsg = createTSG()
new_df.start = pd.to_numeric(new_df.start, errors='coerce')
new_df.end = pd.to_numeric(new_df.end, errors='coerce')
dataframes = {}

for filename in glob.iglob(os.path.join('kidney', '*.vcf')): #useful to open all vcf file in the script directory
    allel.vcf_to_csv(f"{filename}", 'example.csv') #convert from vcf to csv 
    with open('example.csv', mode='r') as f:
        reader = csv.reader(f, delimiter=',')
        i = 0
        included_cols = [0, 1, 8]
        patient=os.path.basename(filename)
        for row in reader:
        
            content = list(row[i] for i in included_cols)
            
            if row[8]=='True':
                a=row[0]
                    
                check=checkPos(new_df, int(row[1]))
                
                if(check[0] > 0):
                    statg[a]= statg[a]+1        #statistica generica
                    if(patient not in stat1[a]):
                        stat1[a][patient]=1
                    else:
                        stat1[a][patient]+=1
                     
                for key in check[1]:
                    if(key not in stat2[a]):
                        stat2[a][key] = {}
                        stat2[a][key][patient] = 1
                        
                        if(checkSupp(a, key, tsg)):
                            stat3[a][key] = {}
                            stat3[a][key][patient] = 1
                    else:
                        if(patient not in stat2[a][key]):
                            stat2[a][key][patient]=1
                            if(key in stat3[a]):
                                stat3[a][key][patient]=1
                        else:
                            stat2[a][key][patient]+=1
                            if(key in stat3[a]):
                                stat3[a][key][patient]+=1
                                
#print(statg)
#print(stat1)
#print(stat2)
#print(stat3)



jsongen = json.dumps(statg) #convert "lista" in a json file to be saved
f = open("statgeneric_kidney.json","w")
f.write(jsongen)
f.close()

jsonlist = json.dumps(stat1) #convert "lista" in a json file to be saved
f = open("stat1_kidney_snv.json","w")
f.write(jsonlist)
f.close()

jsonex = json.dumps(stat2) #convert "lista" in a json file to be saved
f = open("stat2_kidney_snv.json","w")
f.write(jsonex)
f.close()

jsonsupp = json.dumps(stat3) #convert "lista" in a json file to be saved
f = open("stat3_kidney_snv.json","w")
f.write(jsonsupp)
f.close()

# dataframes


for i in CHROMOSOMES:
    df2 = pd.DataFrame(stat2[i]).fillna(0).transpose()
    df3 = pd.DataFrame(stat3[i]).fillna(0).transpose()
    df2['median'] = df2.median(axis=1)
    df3['median'] = df3.median(axis=1)
    df2.to_json('stat2_kidney_snv_updated_'+str(i)+'.json')
    df3.to_json('stat3_kidney_snv_updated_'+str(i)+'.json')
