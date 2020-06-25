# encoding=utf8
import os
import sys
import re
import time
from random import random
from random import choice
from random import sample
from collections import Counter
import string
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from time import time as t
import math
import scipy.stats as stats
import statistics

where=os.popen("uname -a").read().strip()
if re.match("^Linux",where):
    sys.path.append("/home/chenh19/jb/module")
if re.match("^Darwin",where):
    sys.path.append("/jb/module")
import chctool

d_base_comple = {"A":"T","T":"A","C":"G","G":"C","N":"N","a":"T","t":"A","c":"G","g":"C","n":"N","-":"-"}
d_map = {"stop_retained_variant":"synonymous_variant","stop_lost":"stop_lost","inframe_converted":"inframe","frameshift_converted":"frameshift","inframe_insertion":"inframe","start_lost":"missense_variant","inframe_deletion":"inframe","stop_gained":"stop_gained","frameshift_variant":"frameshift","synonymous_variant":"synonymous_variant","missense_variant":"missense_variant", "inframe":"inframe", "frameshift":"frameshift","nonsynonymous SNV":'missense_variant', 'synonymous SNV':'synonymous_variant','stopgain':'stop_gained','frameshift substitution':'frameshift','nonframeshift substitution':'inframe','stoploss':'stop_lost'}


d_conseq_ratio = {'stop_lost': 0.0006149, 'synonymous_variant': 0.354023, 'frameshift': 0.0129278, 'missense_variant': 0.6102157, 'inframe': 0.0094956, 'stop_gained': 0.012723}

d_codon_ratio = {'r1': 0.3056674, 'r2': 0.2816225, 'r3': 0.3902867, 'frameshift': 0.0129278, 'inframe': 0.0094956}

d_codon_table_complete={'CTT': ['L', 'Leu'], 'ATG': ['M', 'Met'], 'AAG': ['K', 'Lys'], 'AAA': ['K', 'Lys'], 'ATC': ['I', 'Ile'], 'AAC': ['N', 'Asn'], 'ATA': ['I', 'Ile'], 'AGG': ['R', 'Arg'], 'CCT': ['P', 'Pro'], 'ACT': ['T', 'Thr'], 'AGC': ['S', 'Ser'], 'ACA': ['T', 'Thr'], 'AGA': ['R', 'Arg'], 'CAT': ['H', 'His'], 'AAT': ['N', 'Asn'], 'ATT': ['I', 'Ile'], 'CTG': ['L', 'Leu'], 'CTA': ['L', 'Leu'], 'CTC': ['L', 'Leu'], 'CAC': ['H', 'His'], 'ACG': ['T', 'Thr'], 'CAA': ['Q', 'Gln'], 'AGT': ['S', 'Ser'], 'CAG': ['Q', 'Gln'], 'CCG': ['P', 'Pro'], 'CCC': ['P', 'Pro'], 'TAT': ['Y', 'Tyr'], 'GGT': ['G', 'Gly'], 'TGT': ['C', 'Cys'], 'CGA': ['R', 'Arg'], 'CCA': ['P', 'Pro'], 'CGC': ['R', 'Arg'], 'GAT': ['D', 'Asp'], 'CGG': ['R', 'Arg'], 'TTT': ['F', 'Phe'], 'TGC': ['C', 'Cys'], 'GGG': ['G', 'Gly'], 'TAG': ['*', 'STOP'], 'GGA': ['G', 'Gly'], 'TGG': ['W', 'Trp'], 'GGC': ['G', 'Gly'], 'TAC': ['Y', 'Tyr'], 'TTC': ['F', 'Phe'], 'TCG': ['S', 'Ser'], 'TTA': ['L', 'Leu'], 'TTG': ['L', 'Leu'], 'CGT': ['R', 'Arg'], 'GAA': ['E', 'Glu'], 'TCA': ['S', 'Ser'], 'GCA': ['A', 'Ala'], 'GTA': ['V', 'Val'], 'GCC': ['A', 'Ala'], 'GTC': ['V', 'Val'], 'GCG': ['A', 'Ala'], 'GTG': ['V', 'Val'], 'GAG': ['E', 'Glu'], 'GTT': ['V', 'Val'], 'GCT': ['A', 'Ala'], 'ACC': ['T', 'Thr'], 'TGA': ['*', 'STOP'], 'GAC': ['D', 'Asp'], 'TCC': ['S', 'Ser'], 'TAA': ['*', 'STOP'], 'TCT': ['S', 'Ser']}
d_codon_table= {'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'ACG': 'T', 'CAA': 'Q', 'AGT': 'S', 'CAG': 'Q', 'CCG': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CCA': 'P', 'CGC': 'R', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TAG': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TGA': '*', 'GAC': 'D', 'TCC': 'S', 'TAA': '*', 'TCT': 'S'}


# d_mutation_table = {1: {'A': {'ACC': [['C', 'T', 'G'], [0, 0.148, 0.3485, 1]], 'ATG': [['C', 'T', 'G'], [0, 0.098, 0.2308, 1]], 'AAG': [['C', 'T', 'G'], [0, 0.2986, 0.3484, 1]], 'AAA': [['C', 'T', 'G'], [0, 0.236, 0.2831, 1]], 'ATC': [['C', 'T', 'G'], [0, 0.1074, 0.2268, 1]], 'AAC': [['C', 'T', 'G'], [0, 0.248, 0.3801, 1]], 'ATA': [['C', 'T', 'G'], [0, 0.0764, 0.1796, 1]], 'AGG': [['C', 'T', 'G'], [0, 0.2524, 0.3806, 1]], 'AGC': [['C', 'T', 'G'], [0, 0.167, 0.299, 1]], 'ACA': [['C', 'T', 'G'], [0, 0.1144, 0.2511, 1]], 'AGA': [['C', 'T', 'G'], [0, 0.2482, 0.2981, 1]], 'AAT': [['C', 'T', 'G'], [0, 0.2757, 0.4068, 1]], 'ATT': [['C', 'T', 'G'], [0, 0.0876, 0.2035, 1]], 'ACT': [['C', 'T', 'G'], [0, 0.1458, 0.313, 1]], 'ACG': [['C', 'T', 'G'], [0, 0.1074, 0.2954, 1]], 'AGT': [['C', 'T', 'G'], [0, 0.1728, 0.3217, 1]]}, 'C': {'CTT': [['A', 'T', 'G'], [0, 0.1643, 0.7392, 1]], 'CCT': [['A', 'T', 'G'], [0, 0.2006, 0.7714, 1]], 'CGT': [['A', 'T', 'G'], [0, 0.049, 0.9464, 1]], 'CGA': [['A', 'T', 'G'], [0, 0.1368, 0.8771, 1]], 'CCA': [['A', 'T', 'G'], [0, 0.1966, 0.7498, 1]], 'CGC': [['A', 'T', 'G'], [0, 0.0588, 0.953, 1]], 'CGG': [['A', 'T', 'G'], [0, 0.0727, 0.9514, 1]], 'CAT': [['A', 'T', 'G'], [0, 0.1669, 0.8784, 1]], 'CTG': [['A', 'T', 'G'], [0, 0.1181, 0.763, 1]], 'CTA': [['A', 'T', 'G'], [0, 0.1148, 0.6857, 1]], 'CTC': [['A', 'T', 'G'], [0, 0.1238, 0.7521, 1]], 'CAC': [['A', 'T', 'G'], [0, 0.16, 0.8648, 1]], 'CAA': [['A', 'T', 'G'], [0, 0.3244, 0.5841, 1]], 'CAG': [['A', 'T', 'G'], [0, 0.2343, 0.6265, 1]], 'CCG': [['A', 'T', 'G'], [0, 0.1728, 0.7746, 1]], 'CCC': [['A', 'T', 'G'], [0, 0.2219, 0.8228, 1]]}, 'T': {'TAT': [['A', 'C', 'G'], [0, 0.1493, 0.8917, 1]], 'TGA': [['A', 'C', 'G'], [0, 0.2061, 0.8491, 1]], 'TGT': [['A', 'C', 'G'], [0, 0.1416, 0.8064, 1]], 'TCT': [['A', 'C', 'G'], [0, 0.2507, 0.7955, 1]], 'TTT': [['A', 'C', 'G'], [0, 0.1777, 0.7607, 1]], 'TGG': [['A', 'C', 'G'], [0, 0.1416, 0.8072, 1]], 'TGC': [['A', 'C', 'G'], [0, 0.1725, 0.7481, 1]], 'TAG': [['A', 'C', 'G'], [0, 0.0647, 0.817, 1]], 'TAA': [['A', 'C', 'G'], [0, 0.1998, 0.7482, 1]], 'TAC': [['A', 'C', 'G'], [0, 0.149, 0.8787, 1]], 'TTC': [['A', 'C', 'G'], [0, 0.1618, 0.7509, 1]], 'TCG': [['A', 'C', 'G'], [0, 0.1702, 0.6972, 1]], 'TTA': [['A', 'C', 'G'], [0, 0.1777, 0.7385, 1]], 'TTG': [['A', 'C', 'G'], [0, 0.0879, 0.8814, 1]], 'TCC': [['A', 'C', 'G'], [0, 0.2591, 0.7516, 1]], 'TCA': [['A', 'C', 'G'], [0, 0.2581, 0.7369, 1]]}, 'G': {'GCA': [['A', 'C', 'T'], [0, 0.7463, 0.819, 1]], 'GTA': [['A', 'C', 'T'], [0, 0.7802, 0.8804, 1]], 'GCC': [['A', 'C', 'T'], [0, 0.7516, 0.823, 1]], 'GTC': [['A', 'C', 'T'], [0, 0.8332, 0.9108, 1]], 'GCG': [['A', 'C', 'T'], [0, 0.7611, 0.8464, 1]], 'GTG': [['A', 'C', 'T'], [0, 0.7775, 0.8818, 1]], 'GAT': [['A', 'C', 'T'], [0, 0.6689, 0.839, 1]], 'GGT': [['A', 'C', 'T'], [0, 0.7174, 0.8607, 1]], 'GGG': [['A', 'C', 'T'], [0, 0.718, 0.8902, 1]], 'GTT': [['A', 'C', 'T'], [0, 0.7626, 0.8722, 1]], 'GCT': [['A', 'C', 'T'], [0, 0.6674, 0.8011, 1]], 'GGC': [['A', 'C', 'T'], [0, 0.8163, 0.9118, 1]], 'GAG': [['A', 'C', 'T'], [0, 0.7654, 0.9529, 1]], 'GGA': [['A', 'C', 'T'], [0, 0.8343, 0.9452, 1]], 'GAC': [['A', 'C', 'T'], [0, 0.7498, 0.8852, 1]], 'GAA': [['A', 'C', 'T'], [0, 0.7574, 0.9357, 1]]}}, 2: {'A': {'TAT': [['C', 'T', 'G'], [0, 0.0572, 0.1736, 1]], 'AAG': [['C', 'T', 'G'], [0, 0.1796, 0.2826, 1]], 'GAG': [['C', 'T', 'G'], [0, 0.2581, 0.4622, 1]], 'AAA': [['C', 'T', 'G'], [0, 0.1915, 0.2957, 1]], 'AAC': [['C', 'T', 'G'], [0, 0.1276, 0.2237, 1]], 'GAC': [['C', 'T', 'G'], [0, 0.1465, 0.4088, 1]], 'GAT': [['C', 'T', 'G'], [0, 0.1126, 0.4111, 1]], 'CAT': [['C', 'T', 'G'], [0, 0.091, 0.2186, 1]], 'AAT': [['C', 'T', 'G'], [0, 0.0545, 0.1089, 1]], 'TAG': [['C', 'T', 'G'], [0, 0.1654, 0.2371, 1]], 'TAA': [['C', 'T', 'G'], [0, 0.148, 0.2081, 1]], 'TAC': [['C', 'T', 'G'], [0, 0.1237, 0.2824, 1]], 'CAC': [['C', 'T', 'G'], [0, 0.1861, 0.3418, 1]], 'CAA': [['C', 'T', 'G'], [0, 0.227, 0.3303, 1]], 'CAG': [['C', 'T', 'G'], [0, 0.1674, 0.2956, 1]], 'GAA': [['C', 'T', 'G'], [0, 0.2775, 0.4376, 1]]}, 'C': {'ACC': [['A', 'T', 'G'], [0, 0.2605, 0.8219, 1]], 'GCA': [['A', 'T', 'G'], [0, 0.1807, 0.8261, 1]], 'ACA': [['A', 'T', 'G'], [0, 0.1429, 0.8286, 1]], 'ACG': [['A', 'T', 'G'], [0, 0.0381, 0.9615, 1]], 'GCG': [['A', 'T', 'G'], [0, 0.072, 0.9494, 1]], 'CCT': [['A', 'T', 'G'], [0, 0.1277, 0.7444, 1]], 'TCT': [['A', 'T', 'G'], [0, 0.1661, 0.6273, 1]], 'GCC': [['A', 'T', 'G'], [0, 0.1817, 0.8294, 1]], 'GCT': [['A', 'T', 'G'], [0, 0.1193, 0.7661, 1]], 'ACT': [['A', 'T', 'G'], [0, 0.14, 0.7089, 1]], 'TCG': [['A', 'T', 'G'], [0, 0.0223, 0.9421, 1]], 'CCG': [['A', 'T', 'G'], [0, 0.0569, 0.9228, 1]], 'CCA': [['A', 'T', 'G'], [0, 0.1687, 0.8137, 1]], 'TCC': [['A', 'T', 'G'], [0, 0.1828, 0.7682, 1]], 'CCC': [['A', 'T', 'G'], [0, 0.1873, 0.7703, 1]], 'TCA': [['A', 'T', 'G'], [0, 0.0757, 0.8676, 1]]}, 'T': {'CTT': [['A', 'C', 'G'], [0, 0.1633, 0.7468, 1]], 'ATG': [['A', 'C', 'G'], [0, 0.1148, 0.8816, 1]], 'GTC': [['A', 'C', 'G'], [0, 0.2013, 0.8321, 1]], 'ATC': [['A', 'C', 'G'], [0, 0.2677, 0.8779, 1]], 'ATA': [['A', 'C', 'G'], [0, 0.09, 0.9204, 1]], 'GTG': [['A', 'C', 'G'], [0, 0.1414, 0.7884, 1]], 'TTT': [['A', 'C', 'G'], [0, 0.1825, 0.7441, 1]], 'ATT': [['A', 'C', 'G'], [0, 0.0556, 0.9339, 1]], 'CTG': [['A', 'C', 'G'], [0, 0.1694, 0.7877, 1]], 'GTT': [['A', 'C', 'G'], [0, 0.1126, 0.8679, 1]], 'CTA': [['A', 'C', 'G'], [0, 0.1174, 0.8371, 1]], 'CTC': [['A', 'C', 'G'], [0, 0.2223, 0.767, 1]], 'TTC': [['A', 'C', 'G'], [0, 0.2255, 0.7367, 1]], 'TTA': [['A', 'C', 'G'], [0, 0.0858, 0.898, 1]], 'TTG': [['A', 'C', 'G'], [0, 0.0666, 0.6923, 1]], 'GTA': [['A', 'C', 'G'], [0, 0.0969, 0.8768, 1]]}, 'G': {'GGT': [['A', 'C', 'T'], [0, 0.5532, 0.7138, 1]], 'CGG': [['A', 'C', 'T'], [0, 0.9083, 0.953, 1]], 'TGT': [['A', 'C', 'T'], [0, 0.6248, 0.8286, 1]], 'AGG': [['A', 'C', 'T'], [0, 0.6737, 0.8798, 1]], 'CGC': [['A', 'C', 'T'], [0, 0.8809, 0.911, 1]], 'AGC': [['A', 'C', 'T'], [0, 0.6428, 0.8495, 1]], 'AGA': [['A', 'C', 'T'], [0, 0.5239, 0.813, 1]], 'CGA': [['A', 'C', 'T'], [0, 0.9084, 0.9519, 1]], 'TGC': [['A', 'C', 'T'], [0, 0.5826, 0.7533, 1]], 'GGG': [['A', 'C', 'T'], [0, 0.5628, 0.8208, 1]], 'TGA': [['A', 'C', 'T'], [0, 0.7405, 0.8827, 1]], 'GGA': [['A', 'C', 'T'], [0, 0.5991, 0.8055, 1]], 'TGG': [['A', 'C', 'T'], [0, 0.5229, 0.7732, 1]], 'GGC': [['A', 'C', 'T'], [0, 0.5776, 0.7493, 1]], 'AGT': [['A', 'C', 'T'], [0, 0.6317, 0.8583, 1]], 'CGT': [['A', 'C', 'T'], [0, 0.9196, 0.9516, 1]]}}, 3: {'A': {'GCA': [['C', 'T', 'G'], [0, 0.2171, 0.3388, 1]], 'GTA': [['C', 'T', 'G'], [0, 0.0963, 0.2322, 1]], 'ACA': [['C', 'T', 'G'], [0, 0.1674, 0.2891, 1]], 'AAA': [['C', 'T', 'G'], [0, 0.1439, 0.2533, 1]], 'ATA': [['C', 'T', 'G'], [0, 0.1215, 0.2699, 1]], 'CGA': [['C', 'T', 'G'], [0, 0.2332, 0.4252, 1]], 'AGA': [['C', 'T', 'G'], [0, 0.1477, 0.3312, 1]], 'TGA': [['C', 'T', 'G'], [0, 0.2996, 0.4655, 1]], 'GGA': [['C', 'T', 'G'], [0, 0.2008, 0.4145, 1]], 'TAA': [['C', 'T', 'G'], [0, 0.2372, 0.2928, 1]], 'CAA': [['C', 'T', 'G'], [0, 0.0976, 0.1533, 1]], 'TCA': [['C', 'T', 'G'], [0, 0.1965, 0.3437, 1]], 'CCA': [['C', 'T', 'G'], [0, 0.1746, 0.2855, 1]], 'TTA': [['C', 'T', 'G'], [0, 0.1274, 0.2968, 1]], 'GAA': [['C', 'T', 'G'], [0, 0.1902, 0.3239, 1]], 'CTA': [['C', 'T', 'G'], [0, 0.1082, 0.1701, 1]]}, 'C': {'ACC': [['A', 'T', 'G'], [0, 0.1176, 0.8405, 1]], 'GCC': [['A', 'T', 'G'], [0, 0.0837, 0.8704, 1]], 'GTC': [['A', 'T', 'G'], [0, 0.1109, 0.8329, 1]], 'ATC': [['A', 'T', 'G'], [0, 0.1412, 0.8769, 1]], 'AAC': [['A', 'T', 'G'], [0, 0.0914, 0.9253, 1]], 'TGC': [['A', 'T', 'G'], [0, 0.0494, 0.9494, 1]], 'CGC': [['A', 'T', 'G'], [0, 0.1141, 0.8856, 1]], 'CTC': [['A', 'T', 'G'], [0, 0.0784, 0.8018, 1]], 'AGC': [['A', 'T', 'G'], [0, 0.0978, 0.904, 1]], 'TTC': [['A', 'T', 'G'], [0, 0.0862, 0.8934, 1]], 'GGC': [['A', 'T', 'G'], [0, 0.1332, 0.9101, 1]], 'TAC': [['A', 'T', 'G'], [0, 0.0384, 0.9645, 1]], 'CAC': [['A', 'T', 'G'], [0, 0.0804, 0.9027, 1]], 'GAC': [['A', 'T', 'G'], [0, 0.0895, 0.9113, 1]], 'TCC': [['A', 'T', 'G'], [0, 0.0777, 0.896, 1]], 'CCC': [['A', 'T', 'G'], [0, 0.0977, 0.8186, 1]]}, 'T': {'TAT': [['A', 'C', 'G'], [0, 0.0402, 0.9588, 1]], 'CTT': [['A', 'C', 'G'], [0, 0.1289, 0.7217, 1]], 'AGT': [['A', 'C', 'G'], [0, 0.1418, 0.8321, 1]], 'TGT': [['A', 'C', 'G'], [0, 0.0749, 0.8748, 1]], 'GGT': [['A', 'C', 'G'], [0, 0.1802, 0.7833, 1]], 'CCT': [['A', 'C', 'G'], [0, 0.1782, 0.7709, 1]], 'TCT': [['A', 'C', 'G'], [0, 0.1192, 0.7728, 1]], 'GAT': [['A', 'C', 'G'], [0, 0.1672, 0.8634, 1]], 'CAT': [['A', 'C', 'G'], [0, 0.1158, 0.9036, 1]], 'AAT': [['A', 'C', 'G'], [0, 0.128, 0.8995, 1]], 'ATT': [['A', 'C', 'G'], [0, 0.1989, 0.7459, 1]], 'GTT': [['A', 'C', 'G'], [0, 0.159, 0.7021, 1]], 'GCT': [['A', 'C', 'G'], [0, 0.1439, 0.8138, 1]], 'ACT': [['A', 'C', 'G'], [0, 0.1573, 0.7565, 1]], 'TTT': [['A', 'C', 'G'], [0, 0.1009, 0.81, 1]], 'CGT': [['A', 'C', 'G'], [0, 0.1735, 0.8384, 1]]}, 'G': {'ATG': [['A', 'C', 'T'], [0, 0.6922, 0.8217, 1]], 'AAG': [['A', 'C', 'T'], [0, 0.708, 0.866, 1]], 'ACG': [['A', 'C', 'T'], [0, 0.8721, 0.9296, 1]], 'GCG': [['A', 'C', 'T'], [0, 0.8549, 0.9162, 1]], 'GTG': [['A', 'C', 'T'], [0, 0.6988, 0.8461, 1]], 'CTG': [['A', 'C', 'T'], [0, 0.6222, 0.831, 1]], 'TTG': [['A', 'C', 'T'], [0, 0.6155, 0.8219, 1]], 'CGG': [['A', 'C', 'T'], [0, 0.6019, 0.7835, 1]], 'AGG': [['A', 'C', 'T'], [0, 0.6741, 0.8242, 1]], 'GGG': [['A', 'C', 'T'], [0, 0.627, 0.7836, 1]], 'TAG': [['A', 'C', 'T'], [0, 0.7557, 0.926, 1]], 'TGG': [['A', 'C', 'T'], [0, 0.448, 0.7384, 1]], 'GAG': [['A', 'C', 'T'], [0, 0.6982, 0.8851, 1]], 'TCG': [['A', 'C', 'T'], [0, 0.8427, 0.93, 1]], 'CCG': [['A', 'C', 'T'], [0, 0.8636, 0.9222, 1]], 'CAG': [['A', 'C', 'T'], [0, 0.6074, 0.852, 1]]}}, 'inframe': 0.0094956, 'frameshift': 0.0129278}

def get_cadd_score(chr,pos,ref,alt):
    # 1    69100    G    A    2.085694    20.2
    x = os.popen("tabix /home/chenh19/ref/cadd/whole_genome_SNVs.tsv.gz %s:%s-%s|cut -f 4,6"%(chr,pos,pos)).read().split("\n")
    for i in x[:3]:
        tmp = i.split("\t")
        if alt == tmp[0]:
            return float(tmp[1])
    print ('pos= ',pos,x)
    return 'error'

def getseq(chr,pos,strand):
    exec ("x = %s[%s]"%(chr,pos))
    x=x.upper()
    if strand =="+":
        return x
    elif strand == "-":
        return "".join([d_base_comple[i] for i in x[::-1]])
def getsinglebase(chr,pos,strand):
    exec ("x = %s[%s]"%(chr,pos))
    x=x.upper()
    if strand =="+":
        return x
    elif strand == "-":
        return d_base_comple[x]

def getcodon(chr, pos, strand, cds_rg):
    idx = cds_rg.index(pos)

    if strand == "+":
        phase = (idx + 1)%3

        if phase ==1:

            pos_codon1 = pos -1
            pos_codon2 = cds_rg[idx+1] -1
            pos_codon3 = cds_rg[idx+2] -1

            seq_codon1 = getsinglebase(chr, pos_codon1, strand)
            seq_codon2 = getsinglebase(chr, pos_codon2, strand)
            seq_codon3 = getsinglebase(chr, pos_codon3, strand)

            base_ref = seq_codon1
            codon_ref = seq_codon1+seq_codon2+seq_codon3
        elif phase ==2:

            pos_codon1 = cds_rg[idx-1] -1
            pos_codon2 = pos -1
            pos_codon3 = cds_rg[idx+1] -1

            seq_codon1 = getsinglebase(chr, pos_codon1, strand)
            seq_codon2 = getsinglebase(chr, pos_codon2, strand)
            seq_codon3 = getsinglebase(chr, pos_codon3, strand)

            base_ref = seq_codon2
            codon_ref = seq_codon1+seq_codon2+seq_codon3
        elif phase ==0:

            pos_codon1 = cds_rg[idx-2] -1
            pos_codon2 = cds_rg[idx-1] -1
            pos_codon3 = pos -1

            seq_codon1 = getsinglebase(chr, pos_codon1, strand)
            seq_codon2 = getsinglebase(chr, pos_codon2, strand)
            seq_codon3 = getsinglebase(chr, pos_codon3, strand)

            base_ref = seq_codon3
            codon_ref = seq_codon1+seq_codon2+seq_codon3


    if strand == "-":
        phase = (idx + 1)%3
        phase = {1:0,2:2,0:1}[phase]

        if phase ==1:

            pos_codon1 = pos -1
            pos_codon2 = cds_rg[idx-1] -1
            pos_codon3 =  cds_rg[idx-2] -1

            seq_codon1 = getsinglebase(chr, pos_codon1, strand)
            seq_codon2 = getsinglebase(chr, pos_codon2, strand)
            seq_codon3 = getsinglebase(chr, pos_codon3, strand)

            base_ref = seq_codon1
            codon_ref = seq_codon1+seq_codon2+seq_codon3
        elif phase ==2:

            pos_codon1 = cds_rg[idx+1] -1
            pos_codon2 = pos -1
            pos_codon3 = cds_rg[idx-1] -1

            seq_codon1 = getsinglebase(chr, pos_codon1, strand)
            seq_codon2 = getsinglebase(chr, pos_codon2, strand)
            seq_codon3 = getsinglebase(chr, pos_codon3, strand)

            base_ref = seq_codon2
            codon_ref = seq_codon1+seq_codon2+seq_codon3
        elif phase ==0:

            pos_codon1 = cds_rg[idx+2] -1
            pos_codon2 = cds_rg[idx+1] -1
            pos_codon3 = pos -1

            seq_codon1 = getsinglebase(chr, pos_codon1, strand)
            seq_codon2 = getsinglebase(chr, pos_codon2, strand)
            seq_codon3 = getsinglebase(chr, pos_codon3, strand)

            base_ref = seq_codon3
            codon_ref = seq_codon1+seq_codon2+seq_codon3
    phase = {1:1,2:2,0:3}[phase]

    return [base_ref, codon_ref, phase]
def get_prob_inverval(prob_interval, prob):
    # {'A': {1: [['C', 'T', 'G'], [0, 0.1563, 0.2195, 1.0]],}
    # [0, 0.1563, 0.2195, 1.0]
    tmp = [1 for i in prob_interval if i<=prob]
    return len(tmp)-1

def getcdslen(rg):
    x =[i[1]-i[0]+1 for i in rg]
    return sum(x)


def rg_split(rg,codon,strand,returnset = 0):
    # codon should be 1,2,3 the phase of codon, if set as 4, means get all the sites in the range
    rg_all = []
    for tmp in rg:
        s=min(tmp[0],tmp[1])
        e = max(tmp[0],tmp[1])
        rg_small = list(range(s,e+1))
        rg_all.extend(rg_small)

    rg_all = sorted(rg_all)
    codon = int(codon)
    if codon ==4:
        res =  rg_all
    elif strand == "+":
        codon = codon -1
        res = rg_all[codon::3]
    elif strand == "-":
        codon = {1:2,2:1,3:0}[codon]
        res = rg_all[codon::3]
    if returnset :
        return set(res)
    else:
        return res

def getphase(rg, pos,strand):
    pos = int(pos)
    x = [i for i in rg if i[0]<= pos and  i[1]>=pos]
    if len(x)==0:
        return (['error',"not in cds"])

    if len(x)>1:
        return(['error', "multi_match"])

    hit = x[0]
    d = pos - hit[0]

    if strand =="+":
        cdspos = hit[2] + d
    elif strand =="-":
        cdspos = hit[2] - d
    else :
        return (['error',"wrong strand"])

    return [{0:3,1:1,2:2}[cdspos%3],cdspos ]


def get_7mer_context_prob(cdsrg,  strand, cdslen,chr_seq ):
    gene_info = []
    n=0
    ctb = {'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'ACG': 'T', 'CAA': 'Q', 'AGT': 'S', 'CAG': 'Q', 'CCG': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CCA': 'P', 'CGC': 'R', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TAG': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TGA': '*', 'GAC': 'D', 'TCC': 'S', 'TAA': '*', 'TCT': 'S'}
    from chctool import revcompl

    home = os.path.expanduser('~')
    dprob_coding = eval(open('%s/ref/mut_prob/7mer_prob_coding.pdict'%(home)).read())
    dprob_gm = eval(open('%s/ref/mut_prob/7mer_prob_genome.pdict'%(home)).read())
    for i in cdsrg:
        s,e = i[0]-1,i[1]
        for pos in range(s,e):
            n+=1

            context = chr_seq[pos-3:pos+4]

            if strand =='+':
                pos_cds = n
                phase = {0:3,1:1,2:2}[pos_cds%3]
                ref = chr_seq[pos]
                codon = chr_seq[pos-phase+1:pos-phase+4]
                aa1 = ctb[codon]
                x = list('ATCG')
                x.remove(ref)
            elif strand =='-':
                context = revcompl(context)
                pos_cds = cdslen - n +1
                phase = {0:3,1:1,2:2}[pos_cds%3]
                ref = revcompl(chr_seq[pos])
                codon = revcompl(chr_seq[pos+phase-3:pos+phase])
                aa1 = ctb[codon]
                x = list('ATCG')
                x.remove(ref)

            for ialt in x:
                codonalt = (ialt if phase==1 else codon[0]) + (ialt if phase==2 else codon[1]) + (ialt if phase==3 else codon[2])
                aa2 = ctb[codonalt]
                if aa1 == aa2:
                    conseq = 'synonymous_variant'
                elif aa1=='*':
                    conseq = 'stop_lost'
                elif aa2 =='*':
                    conseq = 'stop_gained'
                else:
                    conseq='missense_variant'
                k_coding = context+ialt+str(phase)
                k_gm = context+ialt
                gene_info.append([pos,ref,ialt,codon,phase,context,conseq,dprob_coding[k_coding], dprob_gm[k_gm]])
    return gene_info



def do_permutation_multinominal_multinomial(cycles,  total_count, single_permutation_conseq,sum_permutation_conseq , dm_mut, d_dm_rg, mut_prob_normalized, l_mut_prob):
    '''
    # l_mut_prob :  ([ipos, ref, ialt, codon, prob, conseq ]
    # mut_prob_normalized , got from l_mut_prob,  use the prob / total_count
    # dm_mut,  the observed mut/site count in each domain
    # d_dm_rg, the range of each domain, in set format, used to determin if a site is in the domain
    '''


    for i in range(cycles): # perm 2000 times
        t = np.random.multinomial(total_count,mut_prob_normalized)
        for k1,v1 in (single_permutation_conseq.items()):
            for k2,v2 in v1.items():
                single_permutation_conseq[k1][k2] = 0

        for k1,k2 in zip(t,l_mut_prob):
            if k1>0:
                pos = k2[0]+1
                for idm,dmrg in d_dm_rg.items():
                    if pos in dmrg:
                        conseq = k2[-1]
                        if conseq =='stop_lost':
                            continue
                        single_permutation_conseq[idm][conseq] += k1

        for idm,v1 in single_permutation_conseq.items():
            for iconseq , r in v1.items():
                sum_permutation_conseq[idm][iconseq] .append(r)
    res = {}
    for idm,v1 in sum_permutation_conseq.items():
        res[idm] = {}
        for icons, ct in v1.items():
            ct = np.array(ct)
            res[idm][icons] = 'na'
            if idm not in dm_mut:
                observed = 0
            else:
                observed = dm_mut[idm]['observed'][icons][1]
            perm_mean = ct.mean()
            std = ct.std()
            if std !=0:
                zscore = (perm_mean - observed) / std
            else:
                zscore = 0
            res[idm][icons] = [round(zscore,3),round(perm_mean,3),round(std,3),round(observed,3),]
    return res




























def do_permutation_multinominal_7mer(count_type,cycles,total_snp_count, single_permutation_conseq,sum_permutation_conseq, prob_7mer_normalized, prob_7mer_type , gene_info , dm_mut, d_dm_rg):
    res = {}
    '''
    # gene_info , pos, ref,alt, codon, phase, context, conseq, prob_coding,  prob_gm
    # [[55086970, 'A', 'T', 'ATG', 1, 'GCGATGC', 'missense_variant', 4.04131118096e-05, 0.0013740458015299999],
    # [55086970, 'A', 'C', 'ATG', 1, 'GCGATGC', 'missense_variant', 1.34710372699e-05, 0.000458015267176]]
    # prob_7mer_type,  2 types of prob,  from coding  or from intergenic (gm)
    # prob_7mer_normalized , get from gene_info ,e.g. here use prob_7mer_type=gm, sum_gm = sum([i[-1] for i in info])
    # prob_7mer_normalized =  [i[-1] /sum_gm for i in gene_info],  after norm,  the total of all prob = 1
    # dm_mut,  the observed mut/site count in each domain
    # d_dm_rg, the range of each domain, in set format, used to determin if a site is in the domain
    '''
    select_meth = {'site':1,'mut':0}[count_type]   # choose the right observed count in domain

    for i in range(cycles): # perm 2000 times
        t = np.random.multinomial(total_snp_count,prob_7mer_normalized)
        for k1,v1 in (single_permutation_conseq.items()):
            for k2,v2 in v1.items():
                single_permutation_conseq[k1][k2] = 0
        for k1,k2 in zip(t,gene_info):
            if k1>0:
                pos = k2[0]+1
                for idm,dmrg in d_dm_rg.items():
                    if pos in dmrg:
                        conseq = k2[-3]
                        if conseq =='stop_lost':
                            continue
                        single_permutation_conseq[idm][conseq] += k1

        for idm,v1 in single_permutation_conseq.items():
            for iconseq , r in v1.items():
                sum_permutation_conseq[idm][iconseq] .append(r)

    for idm,v1 in sum_permutation_conseq.items():
        res[idm] = {}
        for icons, ct in v1.items():
            ct = np.array(ct)
            res[idm][icons] = 'na'
            if idm not in dm_mut:
                observed = 0
            else:
                observed = dm_mut[idm]['observed'][icons][select_meth]
            perm_mean = ct.mean()
            std = ct.std()
            if std !=0:
                zscore = (perm_mean - observed) / std
            else:
                zscore = 0
            res[idm][icons] = [round(zscore,3),round(perm_mean,3),round(std,3),round(observed,3),]

    return res


def do_permutation_no_mutatable(gn, permu_cycles, conseq, cdsrg,allrg,  d_cds_phase, d_all_pos_codon_info,   chr, strand,  total_sites, total_mut, dm, d_dm_rg, dm_mut, error,single_permutation_conseq,sum_permutation_conseq,count_type,cadd,d_cadd,d_mutation_table,d_codon_ratio,  d_adjust = 1):
    # total SNP count  - random sample from total CDS range, then the site is randomly choose from 3 directions

    len_total  = len(allrg)

    atcg = list('ATCG')

    # this is because, for some analysis, the total site count is got from the intron adjusted, so the total site may be larger than total cds len
    if total_sites > len_total:
        total_sites = len_total

    set_bad_pos = set()

    if conseq not in ['inframe','frameshift']:
        conseqlist = ["synonymous_variant","missense_variant","stop_gained"]

    for cycle in range(permu_cycles):
        # initiaize single permutation
        for k1,v1 in (single_permutation_conseq.items()):
            for k2,v2 in v1.items():
                single_permutation_conseq[k1][k2] = 0

        leftover = 0
        site_list_raw = sample(allrg, total_sites)
        if count_type  == 'site':
            mut_list = Counter(site_list_raw)
        elif count_type  == 'mut':
            total_mut_raw = total_mut + leftover
            leftover = total_mut_raw - int(total_mut_raw)
            mut_list = [choice(site_list_raw) for _ in range(int(total_mut_raw))]
            mut_list = Counter(mut_list)
        else:
            print ('error from do_perm function, wrong counttype, your input = ',repr(count_type))
            return 'error'

        for idm, dmrg in d_dm_rg.items():
            for snp_perm in mut_list:
                if snp_perm in dmrg:
                    if snp_perm not in d_all_pos_codon_info :
                        print ("not in codon rangegn = %s, pos = %s"%(gn,snp_perm),file=error)
                        if  snp_perm not in set_bad_pos:
                            set_bad_pos.add(snp_perm)
                            print ("bad, fail to get codon, gene = %s\tchr=%s\tpos=%s"%(gn,chr,snp_perm), file=error)
                        continue
                    else :
                        base_ref, codon_ref, phase = d_all_pos_codon_info[snp_perm]

                    # the mutation direction is random
                    possible_allele = [_ for _ in atcg if _ != base_ref]
                    base_mut = choice(possible_allele)

                    codon_mut =  (codon_ref[:phase-1] if phase>1 else "")+ base_mut + (codon_ref[phase:] if phase<3 else  "")

                    aa_ref = d_codon_table[codon_ref]
                    aa_mut = d_codon_table[codon_mut]

                    cadd_score  = 1

                    if aa_ref == aa_mut :
                        single_permutation_conseq[idm]["synonymous_variant"] +=mut_list[snp_perm] * cadd_score

                    elif aa_ref =="*" :
                        pass
                    elif aa_mut == "*":
                        single_permutation_conseq[idm]["stop_gained"] +=mut_list[snp_perm] * cadd_score
                    else :
                        single_permutation_conseq[idm]["missense_variant"] +=mut_list[snp_perm]  * cadd_score

        for idm,v1 in single_permutation_conseq.items():
            for iconseq, r in v1.items():
                sum_permutation_conseq[idm][iconseq].append(r)

    # calc the pvalue
    res = {}
    for i_dm,v1 in sum_permutation_conseq.items():
        res[i_dm] ={}
        for i_conseq, conseq_count in v1.items():
            conseq_count = np.array(conseq_count)
            # s = time()
            res[i_dm][i_conseq] = 'na'

            # for some domain, it's not involved in the exac_count, because there is no  SNP found in the vcf file in this region, for this type, just let the observed  0
            if i_dm not in dm_mut:
                print('i_dm not found', i_dm)
                observed = 0
            else:
                try:
                    # {'observed': {'missense_variant': [23.1, 8.0, 706.27, 239.2],
                    if cadd:
                        observed = dm_mut[i_dm]["observed"][i_conseq][{'bysite':3,'bymut':2,'site':3,'mut':2,'by_site':3,'by_mut':2,}[count_type]]
                    else:
                        observed = dm_mut[i_dm]["observed"][i_conseq][{'bysite':1,'bymut':0,'site':1,'mut':0,'by_site':1,'by_mut':0,}[count_type]]
                except:
                    print ( 'do_permutation fucntion, error get observed count\n\t',dm_mut[i_dm]["observed"][i_conseq], idm, i_conseq)

            groupmean = conseq_count.mean()
            # groupmedian = statistics.median(conseq_count)
            std = conseq_count.std()

            if std !=0:
                #  if expect > observed, zscore is positive
                zscore = (groupmean - observed)/std
            else:
                zscore = 0
            res[i_dm][i_conseq] = [round(zscore,3), round(groupmean,3), round(std,3),round(observed,3)]
    return res



def do_permutation(gn, permu_cycles, conseq, cdsrg,allrg,  d_cds_phase, d_all_pos_codon_info,   chr, strand,  total_sites, total_mut, dm, d_dm_rg, dm_mut, error,single_permutation_conseq,sum_permutation_conseq,count_type,cadd,d_cadd,d_mutation_table,d_codon_ratio,  d_adjust = 1):
    # pass in , total sites, snp counts in gene
    # 1. select `total_sites` <<unique>> sites from the whole cdsrange
    # 2. if conseq == mis/syno/stop_gain,  will expand the snp_site_list, first tally the phase of all selected sites,  the SNP  mutation  count ratio happens in 3 phase are
    #{1:0.326170056,2:0.300468946,3:0.373360998}
    # suppose the site number for codon-1/2/3 = n1/n2/n3
    # find the bigest site number * 1.25 as the base, then duplicate the other 2  ( 1.25 ~  0.3734/0.3, ensure that , after all other 2 codon  site number would be expand)
    # for codon1, n1_new = n2 * (0.3262/0.3004), the number would be added should be n1_new - n1, because the added new sites would be from the original site_list, so the actuall  unique site number wouln't change.

    # if conseq == 'inframe/frameshift', then, just select `total_sites`  sites from  all_rg,  then, inject all the mut to these sites
    # d_codon_ratio = {1:0.326170056,2:0.300468946,3:0.373360998}
    len_total  = len(allrg)

    # this is because, for some analysis, the total site count is got from the intron adjusted, so the total site may be larger than total cds len
    if total_sites > len_total:
        total_sites = len_total


    set_bad_pos = set()
    leftover = 0
    leftover1 = 0
    nerror =0


    if conseq not in ['inframe','frameshift']:
        conseqlist = ["synonymous_variant","missense_variant","stop_gained"]


    for cycle in range(permu_cycles):
        # initiaize single permutation
        for k1,v1 in (single_permutation_conseq.items()):
            for k2,v2 in v1.items():
                single_permutation_conseq[k1][k2] = 0
        # sample
        p1 = [_ for _,v in d_cds_phase.items() if v==1]
        p2 = [_ for _,v in d_cds_phase.items() if v==2]
        p3 = [_ for _,v in d_cds_phase.items() if v==3]

        total_sites_new = total_sites + leftover1
        ct1 = int(round(total_sites_new*d_codon_ratio[1],0))
        ct2 = int(round(total_sites_new*d_codon_ratio[2],0))
        ct3 = int(round(total_sites_new*d_codon_ratio[3],0))
        leftover1 = total_sites - ct1 -ct2 -ct3

        site_list_raw = []
        site_list_raw +=sample(p1,ct1)
        site_list_raw +=sample(p2,ct2)
        site_list_raw +=sample(p3,ct3)

        # this is the previous version to mimic the mutation ratio on each phase,
        # this is not accurate, because this would change the actuall total mutation count, and  this number is all that matter to calculate the z score
        # so we need to sample from differnt phase, rather than from total cds region
        ##################
        # if conseq not in ['inframe','frameshift']:
        #     d_site_phase_count={1:0,2:0,3:0}
        #     site_list_phase = {1:[],2:[],3:[]}
        #     for i in site_list_raw:
        #         phase = d_cds_phase[i]
        #         d_site_phase_count[phase] +=1
        #         site_list_phase[phase].append(i)

        #     site_phase_count = sorted(d_site_phase_count.items(),key=lambda x:x[1])

        #     # choose the starting site num
        #     # starter / ratio_higest_count  * lowest_ratio  > second highest site num
        #     starter = site_phase_count[1][1] * d_codon_ratio[site_phase_count[-1][0]]  / min(d_codon_ratio.values())
        #     starter  = max(int(starter) +1, site_phase_count[-1][1])

        #     for i in range(1,4):
        #         tmp = d_codon_ratio[i]/d_codon_ratio[site_phase_count[-1][0]]
        #         new_site_len = int(tmp* starter)
        #         gap = new_site_len - d_site_phase_count[i]

        #         # this step, will modify the situation, when  phase_range_cds[i] is empty, the coice will error
        #         if d_site_phase_count[i] ==0:
        #             choice_new = rg_split(cdsrg, i,strand)
        #         else:
        #             choice_new = site_list_phase[i]

        #         expand = [choice(choice_new) for _ in range(gap)]
        #         site_list_raw.extend(expand)

        if count_type  == 'site':
            mut_list = Counter(site_list_raw)
        elif count_type  == 'mut':
            total_mut_raw = total_mut + leftover
            leftover = total_mut_raw - int(total_mut_raw)
            mut_list = [choice(site_list_raw) for _ in range(int(total_mut_raw))]
            mut_list = Counter(mut_list)
        else:
            print ('error from do_perm function, wrong counttype, your input = ',repr(count_type))
            return 'error'

        if conseq == 'frameshift':
            for idm, v in dm.items():
                if conseq not in dm_mut[idm]['observed']:
                    continue
                for snp_perm in mut_list:
                    dmrg = v[0]
                    if test_protein_trunc(dmrg, strand, snp_perm):
                        single_permutation_conseq[idm]['frameshift'] += mut_list[snp_perm]

        elif conseq == 'inframe' :
            for idm, dmrg in d_dm_rg.items():
                if conseq not in dm_mut[idm]['observed']:
                    continue
                for snp_perm in mut_list:
                    if snp_perm in dmrg:
                        single_permutation_conseq[idm]['inframe'] += mut_list[snp_perm]

        else:
            for idm, dmrg in d_dm_rg.items():
                for snp_perm in mut_list:
                    if snp_perm in dmrg:
                        if snp_perm not in d_all_pos_codon_info :
                            print ("not in codon rangegn = %s, pos = %s"%(gn,snp_perm),file=error)
                            if  snp_perm not in set_bad_pos:
                                set_bad_pos.add(snp_perm)
                                print ("bad, fail to get codon, gene = %s\tchr=%s\tpos=%s"%(gn,chr,snp_perm), file=error)
                            continue
                        else :
                            base_ref, codon_ref, phase = d_all_pos_codon_info[snp_perm]

                        # generate a random number in [0,1)
                        direction = random()
                        mutation_table = d_mutation_table[phase][base_ref][codon_ref]
                        mutation_table_prob = mutation_table[1]    # [0, 0.1563, 0.2195, 1.0]
                        mutation_table_seq = mutation_table[0]  # ['C', 'T', 'G']

                        idx_tmp = get_prob_inverval(mutation_table_prob, direction)

                        base_mut = mutation_table_seq[idx_tmp]
                        codon_mut =  (codon_ref[:phase-1] if phase>1 else "")+ base_mut + (codon_ref[phase:] if phase<3 else  "")

                        aa_ref = d_codon_table[codon_ref]
                        aa_mut = d_codon_table[codon_mut]

                    # this step may be slow, may need to modify into batch tabix mode -R
                        if  cadd:
                            basemut1 = base_mut if strand =='+' else d_base_comple[base_mut]
                            try:
                                cadd_score = d_cadd[snp_perm][basemut1]
                            except:
                                print ('error, site not found in cadd dict, site = %s, strand = %s, ref, alt = %s,%s '%(snp_perm,strand,base_ref,base_mut),file=error  )
                                continue
                        else:
                            cadd_score = 1
                        if aa_ref == aa_mut :
                            single_permutation_conseq[idm]["synonymous_variant"] +=mut_list[snp_perm] * cadd_score

                        elif aa_ref =="*" :
                            pass
                        elif aa_mut == "*":
                            single_permutation_conseq[idm]["stop_gained"] +=mut_list[snp_perm] * cadd_score
                        else :
                            single_permutation_conseq[idm]["missense_variant"] +=mut_list[snp_perm]  * cadd_score

        for idm,v1 in single_permutation_conseq.items():
            for iconseq, r in v1.items():
                sum_permutation_conseq[idm][iconseq].append(r)

    # calc the pvalue
    res = {}

    for i_dm,v1 in sum_permutation_conseq.items():
        res[i_dm] ={}
        for i_conseq, conseq_count in v1.items():
            conseq_count = np.array(conseq_count)
            # s = time()
            res[i_dm][i_conseq] = 'na'

            # for some domain, it's not involved in the exac_count, because there is no  SNP found in the vcf file in this region, for this type, just let the observed  0
            if i_dm not in dm_mut:
                observed = 0
            else:
                try:
                    # {'observed': {'missense_variant': [23.1, 8.0, 706.27, 239.2],
                    if cadd:
                        observed = dm_mut[i_dm]["observed"][i_conseq][{'bysite':3,'bymut':2,'site':3,'mut':2,'by_site':3,'by_mut':2,}[count_type]]
                    else:
                        observed = dm_mut[i_dm]["observed"][i_conseq][{'bysite':1,'bymut':0,'site':1,'mut':0,'by_site':1,'by_mut':0,}[count_type]]
                except:
                    print ( 'do_permutation fucntion, error get observed count\n\t',dm_mut[i_dm]["observed"][i_conseq], idm, i_conseq)

            groupmean = conseq_count.mean()
            # groupmedian = statistics.median(conseq_count)
            std = conseq_count.std()

            if std !=0:
                #  if expect > observed, zscore is positive
                zscore = (groupmean - observed)/std
            else:
                zscore = 0
            res[i_dm][i_conseq] = [round(zscore,3), round(groupmean,3), round(std,3),round(observed,3)]
    return res
            # if observed > simulation mean, pvalue = negative, means this site is more trend to mut
            # if observed < simulation mean, pvalue = positive, means this site is more conserved`
            # if groupmean - observed <0:
            #     # z-score is negative
            #     if     stats.norm.cdf(zscore)==0:
            #         # if z score <-37, then  pvalue would be 0 (z-score=37, pvalue = 5e-300)
            #         pvalue =300
            #     else:
            #         pvalue = math.log(stats.norm.cdf(zscore),10)
            # else:
            #     # z-score >=0
            #     if     stats.norm.cdf(zscore)==0:
            #         # if z score <37, then  pvalue would be 0 (z-score=37, pvalue = 5e-300)
            #         pvalue =-300
            #     else:
            #         pvalue = -math.log(stats.norm.cdf( -zscore),10)



def test_protein_trunc(rg,strand,pos):
    if strand == "+":
        point = max(rg[-1][1], rg[0][1])
        rs = (0,1)[pos <point]
    elif strand == "-":
        point = min(rg[0][0],rg[-1][0])
        rs = (0,1)[pos >point]
    return rs

def get_cds_len(rg):
    x =[i[1]-i[0]+1 for i in rg]
    return sum(x)


def get_all_codon_in_cds(chr_seq, cds_rg, strand):
    # chr is the real total chrom sequence

    d_all_codon={}   # store all the codon result
    allrg = rg_split(cds_rg, 4,strand)

    # print (allrg[:10])

    cdslen = get_cds_len(cds_rg)
    chr = chr_seq

    seqall = ""
    flag = "good"

    stop_codon ="TAG,TGA,TAA".split(",")
    if strand not in ["-","+"]:
        return ["error","strand not in + -"]

    # step1, get expected cds len, and first 4bp, last 4nt

    seq1 = [chr[i-1] for i in allrg[0:3]]   # expected ATG
    pos1a = allrg[0]-1
    seq1a = chr[pos1a -1]   # pos -1
    pos1b = allrg[3]
    seq1b = chr[pos1b -1]   # pos 4

    seq2 = [chr[i-1] for i in  allrg[-3:]]  # expected stop codon
    pos2a = allrg[-1] +1
    seq2a = chr[pos2a -1]  # end pos +1
    pos2b = allrg[-4]
    seq2b = chr[pos2b -1]  # 4th to the last

    pos1 = allrg[0]
    pos2 = allrg[-1]


    if strand =='-':
        rev_seq1 = [d_base_comple[i] for i in seq2][::-1]
        rev_seq2 = [d_base_comple[i] for i in seq1][::-1]

        rev_seq1a = d_base_comple[seq2a]
        rev_seq1b = d_base_comple[seq2b]

        rev_seq2a = d_base_comple[seq1a]
        rev_seq2b = d_base_comple[seq1b]

        revpos1a = pos2a
        revpos1b = pos2b
        revpos2a = pos1a
        revpos2b = pos1b

        seq1 = rev_seq1
        seq2 = rev_seq2
        seq1a = rev_seq1a
        seq1b = rev_seq1b
        seq2a = rev_seq2a
        seq2b = rev_seq2b

        pos1a = revpos1a
        pos1b = revpos1b
        pos2a = revpos2a
        pos2b = revpos2b

        pos1 = allrg[-1]
        pos2 = allrg[0]

    # in case , some ts cds_range, won't include the stopcodon in the seq, check the trailing 3bp
    if strand =='+':
        check_stop = ''.join([chr[i-1] for i in [allrg[-1]+1, allrg[-1]+2, allrg[-1]+3]])
    if strand =='-':
        check_stop = [chr[i-1] for i in [allrg[0]-1, allrg[0]-2, allrg[0]-3]]
        check_stop = [d_base_comple[i] for i in check_stop]
        check_stop = ''.join(check_stop)


    shift_start, shift_end = 100,100


    if ''.join(seq1) == 'ATG':
        shift_start = 0
    if seq1a+seq1[0]+seq1[1] =='ATG':  # left shift, if rescue,  will add pos1a, total_len +1
        shift_start = 1
    elif seq1[1]+seq1[2]+seq1b =='ATG':  # right shift, if rescue , will remove seq1[0], total_len -1
        shift_start = -1


    if ''.join(seq2) in stop_codon:
        shift_end = 0
    if seq2b+seq2[0]+seq2[1] in  stop_codon:  # left shift, if rescue,  will remove seq2[2], total_len -1
        shift_end = -1
    elif seq2[1]+seq2[2]+seq2a in  stop_codon:  # right shift, if rescue , will add pos2a, total_len +1
        shift_end = 1
    elif check_stop in  stop_codon:
        shift_end = 0

#     rg_for_calc = allrg[:]
    if shift_start + shift_end ==200:
        return ['error', 'both start and stop codon unfound']
    elif shift_start ==100:
        return ['error','start codon not found']
    elif shift_end == 100:
        return ['error','stop codon not found']
    elif shift_start == 0 and shift_end ==0:
        if cdslen%3 !=0:
            return ['error','start and stop codon all found, but, cdslen(%s) %% 3 !=0  (%s)'%(cdslen, cdslen%3)]
    elif  shift_start == 0 and shift_end ==1:   # add one base
        if (cdslen +1)%3 !=0:
            return ['error', 'extend 1bp find stop codon, but new_cdslen%s %%3 !=0  (%s)'%(cdslen +1,(cdslen +1)%3 )]
        else:
            allrg.append(pos2a)

    elif shift_start == 0 and shift_end == -1:   # del one base
        if (cdslen -1)%3 !=0:
            return ['error', 'del 1bp find stop codon, but new_cdslen%s %%3 !=0  (%s)'%(cdslen -1,(cdslen -1)%3 )]
        else:
            allrg.remove(pos2)

    elif  shift_start == 1 and shift_end ==0:   # add one base
        if (cdslen +1)%3 !=0:
            return ['error', 'extend 1bp find start codon, but new_cdslen%s %%3 !=0  (%s)'%(cdslen +1,(cdslen +1)%3 )]
        else:
            allrg.append(pos1a)

    elif shift_start == 1 and shift_end ==1:   # del one base
        if (cdslen +2 )%3 !=0:
            return ['error', 'expand 1bp find start codon, expand 1bp find stop codon, but new_cdslen%s %%3 !=0  (%s)'%(cdslen +2,(cdslen +2)%3 )]
        else:
            allrg.append(pos2a)
            allrg.append(pos1a)

    elif  shift_start == 1 and shift_end ==-1:   # add one base
        if (cdslen )%3 !=0:
            return ['error', 'globally left shifted, but new_cdslen%s %%3 !=0  (%s)'%(cdslen +1,(cdslen +1)%3 )]
        else:
            rgnew = [[i[0]-1,i[1]-1 ] for i in cds_rg]
            allrg =  rg_split(rgnew, 4, strand)

    elif  shift_start == -1 and shift_end ==0:   # add one base
        if (cdslen -1)%3 !=0:
            return ['error', 'del 1bp find start codon, but new_cdslen%s %%3 !=0  (%s)'%(cdslen -1,(cdslen -1)%3 )]
        else:
            allrg.remove(pos1)

    elif shift_start == -1 and shift_end ==1:   # del one base
        if (cdslen )%3 !=0:
            return ['error', 'globally right shifted, but new_cdslen%s %%3 !=0  (%s)'%(cdslen,(cdslen)%3 )]
        else:
            rgnew = [[i[0]+1,i[1]+1 ] for i in cds_rg]
            allrg =  rg_split(rgnew, 4, strand)

    elif  shift_start == -1 and shift_end ==-1:   # add one base
        if (cdslen -2)%3 !=0:
            return ['error', 'del 1bp find start codon, del 1bp find stop codon, but new_cdslen%s %%3 !=0  (%s)'%(cdslen -2,(cdslen -2)%3 )]
        else:
            allrg.remove(pos1)
            allrg.remove(pos2)


    if strand =="+":

        allrg = sorted(allrg)
        allseq = [chr[i-1] for i in allrg]
        n_aa = int(len(allrg)/3)
        allcodon = [allseq[3*i]+allseq[3*i+1]+allseq[3*i+2] for i in range(n_aa)]
        phase = [1,2,3]* n_aa

        for n,pos in enumerate(allrg):
            d_all_codon[pos]=[allseq[n],allcodon[int(n/3)], phase[n]]

    elif strand =="-":

        allrg = sorted(allrg,reverse=True)

        allseq = [d_base_comple[chr[i-1]] for i in allrg]
        n_aa = int(len(allrg)/3)
        allcodon = [allseq[3*i]+allseq[3*i+1]+allseq[3*i+2] for i in range(n_aa)]
        phase = [1,2,3]* n_aa

        for n,pos in enumerate(allrg):
            d_all_codon[pos]=[allseq[n],allcodon[int(n/3)], phase[n]]

    return ['good',d_all_codon]
