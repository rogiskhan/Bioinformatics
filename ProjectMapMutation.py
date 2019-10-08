#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 22:10:43 2019

@author: simonci
"""

import allel
import csv
lista={} #it's possible to create all the main keys (chr1, chr2...) before the operation and delete the control if(lista.get(a)==None)
allel.vcf_to_csv('a.vcf', 'example.csv')
with open('example.csv', mode='r') as f:
    reader = csv.reader(f, delimiter=',')
    included_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8]      
    for row in reader:
        content = list(row[i] for i in included_cols)
        if row[8]=='True': 
            a=row[0]
            if(lista.get(a) == None):
                lista[a]={}
            if(lista[a].get(row[1])==None):
                lista[a][row[1]]=1
            else:
                lista[a][row[1]]+=1
    
