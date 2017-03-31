#!/usr/bin/env python

import pandas as pd
import numpy as np
import re
from sys import argv

script, infile, filter, mapfile, outfile = argv

fname = re.findall(r"[\w']+", infile)
if fname[-1] == "frq":
    map=pd.read_csv(mapfile, delim_whitespace=True, header=None)
    map=map.rename(columns={0:"CHR", 1:"SNP", 2:"CM", 3:"BP"})
    frq = pd.read_csv(infile, delim_whitespace=True)
    frq["Callrate"]= frq["NCHROBS"]/max(frq["NCHROBS"])
    removed = frq[frq["Callrate"] < 1-float(filter)]
    off = removed[removed["CHR"] != 0]
    off=pd.merge(map, off, on="SNP")
    rm = pd.DataFrame({'Filtering_Step':"SNP_callrate", "CHR":off["CHR_x"], "SNP":off["SNP"], "SNP_Filter_Value": off["Callrate"], "Assay":fname[4], "BP":off["BP"]})
    rm = rm[["CHR", "SNP", "BP","Assay", "Filtering_Step", "SNP_Filter_Value"]]
if fname[-1] == "hwe":
    map=pd.read_csv(mapfile, delim_whitespace=True, header=None)
    map=map.rename(columns={0:"CHR", 1:"SNP", 2:"CM", 3:"BP"})
    hwe = pd.read_csv(infile, delim_whitespace=True)
    removed = hwe[hwe["P"] < float(filter)]
    removed = pd.merge(map, removed, on="SNP")
    rm = pd.DataFrame({'Filtering_Step':"hwe_pvalue", "CHR":removed["CHR_x"], "SNP":removed["SNP"], "SNP_Filter_Value": removed["P"], "Assay":fname[4], "BP":removed["BP"]})
    rm = rm[["CHR", "SNP", "BP","Assay", "Filtering_Step", "SNP_Filter_Value"]]

rm.to_csv(outfile, sep="\t", index=False, header=None)
