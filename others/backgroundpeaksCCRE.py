import os
import sys
import pandas as pd
import random
import subprocess
# import matplotlib.pyplot as plt
import pickle as pk
# from giggle import Giggle
# import json
import shutil


bed = pd.read_csv("GRCh38_cCREs_GC.bed", sep="\t", comment="#", header = None)
bed["GC_sum"] = bed.iloc[:,[9,10]].sum(axis=1)
bed = bed.iloc[:,[0,1,2,14,15]]
bed = bed.rename(columns={0: "chr", 1:"start", 2: "end", 14: "length"})

GC_content_record = {}
folder = "backgroundPeaks"
shutil.rmtree(folder)
for i in range(0,999999):
    if i % 1000 == 0:
        subfolder = str(i)
        print(subfolder)
        os.makedirs(folder + "/" + subfolder)
    peaks_number = random.randint(4000, 6000)
    peaks_index = random.sample(range(0, bed.__len__()), peaks_number)
    peaks_index.sort()
    tmpbed = bed.loc[peaks_index,:]
    GC_content = tmpbed.GC_sum.sum()/tmpbed.length.sum()
    GC_content_record[str(i)] = GC_content
    tmpbed = bed.loc[peaks_index,["chr", "start", "end"]]
    output_bed_file = "backgroundPeaks/{subfolder}/{id}.bed".format(id=str(i), subfolder=subfolder)
    tmpbed.to_csv(output_bed_file, sep="\t", header = None, index=0)
    cmd = "sort -k1,1 -k2,2n -k3,3n {bed_file} | bgzip -c > {bed_file}.gz\nrm {bed_file}".format(bed_file = output_bed_file)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    


with open("GC_content_record.pk", "wb") as output_file:
    pk.dump(GC_content_record, output_file)