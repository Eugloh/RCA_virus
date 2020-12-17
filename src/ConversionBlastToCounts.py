#!/usr/bin/env python3
import sys  # pour utiliser sys.argv et sys.exit().
import os
import csv
import numpy as np
print('numpy: '+np.version.full_version)
import matplotlib
print('matplotlib: '+matplotlib.__version__)
import matplotlib.pyplot as plt
from os import listdir


# OutputBlastVIRUSES.blast in data folder 
# QueryID	SubjectID	evalue	bitscore	length query	perc id	frames	taxid	kingdom	scientifique name	common name	blast name	title	seq query	startq	stopq
# RCA_Skin5_S17_L00184;clusterid=83;size=797	gi|327195194|gb|JF304769.1|	1.57E-71	279	151	100	1/-1	37957	Viruses	Human papillomavirus type 36	Human papillomavirus type 36	viruses	Human papillomavirus type 36 isolate Muc17.2, complete genome	GCTCCAGACCGTGTGCGAATAGTAGCTCTTGTTCCCAATCTACTAACCCTGACATACCCTGCAGGAGTTGTAGAATATTGGGGACGTCCCAACCGTTGCACATCAAGAAAATCCCTATCTGGTGGTTCCTCAAAAACATCTACATCATGTT	1	151

src_path = os.path.dirname(os.path.realpath(__file__)) # To get the full path to the directory a Python file is contained in 


print(src_path)
listpath = src_path.split("/")
del listpath[-1]
datapath="/".join(["/".join(listpath),"data"])
print(datapath)


files = os.listdir(datapath)


# while 'comptes.csv' in files:
#     del files[files.index('comptes.csv')]


env=["QueryID",	"SubjectID",	"evalue",	"bitscore",	"length query",	"perc id",	"frames",	"taxid",	"kingdom",	"scientifique name"	,"common name"	,"blast name"	,"title",	"seq query",	"startq",	"stopq","origin","size"]

Every={}
listPool=[]
listeverynamevirus=[]
Taxidcorrespond={}
BLAST="/".join([datapath,'OutputBlastVIRUSES.blast'])


with open(BLAST, newline='\n') as csvfile:
    next(csvfile)
    spamreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
    for row in spamreader:
        row.append("_".join(row[0].split("_")[:2])) # RCA_Skin5_S17_L00184;clusterid=83;size=797 => Skin5
        row.append(row[0].split("=")[2]) # 797
        if "_".join(row[0].split("_")[:2]) not in listPool:
            listPool.append("_".join(row[0].split("_")[:2]))
        if row[10] not in listeverynamevirus:
            listeverynamevirus.append(row[10])
            Taxidcorrespond[row[10]]=row[7]
        if row[0] not in Every:
            Every[row[0]]={}
            for el in range(len(env)):
                if env[el] in ["size","evalue","bitscore","length query","perc","startq","stopq"]:
                    Every[row[0]][env[el]]=float(row[el])
                else:
                    Every[row[0]][env[el]]=row[el]
        else:
            if len(row[-5]) > len(Every[row[0]]["seq query"]):
                Every[row[0]]={}
                for el in range(len(env)):
                    if env[el] in ["size","evalue","bitscore","length query","perc","startq","stopq"]:
                        Every[row[0]][env[el]]=float(row[el])
                    else:
                        Every[row[0]][env[el]]=row[el]


# # Every contain keys as : RCA-Vulva-9_S958596;clusterid=58595;size=2

plot1={}
for key,value in Every.items():
    if Every[key]["common name"] not in plot1:
        plot1[Every[key]["common name"]]={}
        for EveryPool in listPool:
            plot1[Every[key]["common name"]][EveryPool]=0
    plot1[Every[key]["common name"]][Every[key]["origin"]]+=Every[key]["size"]
firstraw=["Virus","taxid"]
for id in listPool:
    firstraw.append(id)
rawsplot1=[]
rawsplot1.append(firstraw)
for key in plot1.keys():# les noms de virus
    tmpraw=[key]
    tmpraw.append(Taxidcorrespond[key])
    for value in plot1[key].values():#pour les sous dictionnaires
        tmpraw.append(value)
    rawsplot1.append(tmpraw)

print(rawsplot1)

with open("/".join([datapath,"comptes.csv"]), 'w', newline='' ) as csvfile:
        spamwriter = csv.writer(csvfile, delimiter='\t',quotechar='"')
        for element in rawsplot1:
            spamwriter.writerow(element)

