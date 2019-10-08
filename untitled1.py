# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 10:09:47 2019

@author: monal
"""
import allel
import csv
f1=open('chr1.txt','r')
chr1=f1.read()
f1.close()
f2=open('chr2.txt','r')
chr2=f2.read()
f2.close()
lista={"chrM":chr1, "chr2":chr2}
allel.vcf_to_csv('C:\\Users\\monal\\Documents\\Politecnico\\BIOINFO\\project\\tcga_lung_kidney\\unzipped\\kidney\\0afdd78d-57da-4736-9971-ef1d76d0c3e6.vep.vcf', 'example.csv')
with open('example.csv', mode='r') as f:
    reader = csv.reader(f, delimiter=',')
    included_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8]      
    for row in reader:
        content = list(row[i] for i in included_cols)
        if row[8]=='True': 
            for i in lista.values():
                a=row[0]
                crom=lista.get(f"{a}")
                
            
            
            