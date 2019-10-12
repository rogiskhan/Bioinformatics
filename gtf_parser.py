from collections import defaultdict
import gzip
import pandas as pd
import re


GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')


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

















df = dataframe("........")
df[['feature','start','end']].head(100)



new_df = df[['seqname','feature','start','end']].loc[df['feature'] == 'exon']

new_df.start = pd.to_numeric(new_df.start, errors='coerce')
new_df.end = pd.to_numeric(new_df.end, errors='coerce')
# new_df = new_df.sort_values(by=['start'],ascending=True)
new_df.head(5)


lista={} #it's possible to create all the main keys (chr1, chr2...) before the operation and delete the control if(lista.get(a)==None)
for filename in glob.iglob('D:\\bio\\kidney\\*.vcf'): #useful to open all vcf file in the script directory
    allel.vcf_to_csv(f"{filename}", 'example.csv') #convert from vcf to csv 
    with open('example.csv', mode='r') as f:
        reader = csv.reader(f, delimiter=',')
        included_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8]      
        for row in reader:
            content = list(row[i] for i in included_cols)
            
            if row[8]=='True': 
                a=row[0]              
                if(lista.get(a) == None):
                    lista[a]={}
                    df_ch = new_df.loc[str(row[0]).split("chr")[1] == new_df['seqname']]
                    df_start = df_ch.loc[new_df['start'] <= int(row[1])]
                    df_end = df_start.loc[df_start['end'] >= int(row[1])]
                    print(str(row[0])+ ": "+str(row[1])+': \r\n')
                    print(df_end)
                if(lista[a].get(row[1])==None):
                    lista[a][row[1]]=1
                else:
                    lista[a][row[1]]+=1  #create the dictionary called "lista" with counter of mutations
# json = json.dumps(lista) #convert "lista" in a json file to be saved
# f = open("chr_map.json","w")
# f.write(json)
# f.close()