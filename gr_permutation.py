# previous version, the total permuation num = total mutations  (one site could have multiple times mutation)
# in this version , the perm num = site count


import re
import os
import sys
import gzip
import importlib

where = os.popen("uname -a").read().strip()
if re.match("^Linux", where):
    sys.path.append("/home/chenh19/jb/module")
if re.match("^Darwin", where):
    sys.path.append("/jb/module")

from random import choice
from random import random
import scipy.stats as stats
import statistics
import math
from time import time

import permutation_functions
import chctool
importlib.reload(permutation_functions)
importlib.reload(chctool)

from chc_tool import f_make_dirs
from chc_tool import merge_range


from permutation_functions import *
import argparse
import pprint
p = pprint.pprint
chrmap = chctool.chrmap()


# total mut count in exac_non_tcga = 19634231.6
# total missense in exac_non_tcga = 11754032.9
# total synonymous in exac_non_tcga = 7233291.6

# genome overall site count
# frameshift              82549
# inframe                 33794
# missense_variant      2793399
# stop_gained             77214
# stop_lost                3028
# synonymous_variant    1401822

# frameshift             1.88%
# inframe                0.77%
# missense_variant       63.6%
# stop_gained            1.76%
# stop_lost              0.07%
# synonymous_variant    31.92%

# arguments = 4.  1=genelist, 2=project_name, 3= chr using (optional, def = all chrs, 4=count_type, site/mut, default = by site   --perm_cycles, --n,-n, -perm, the permutation times default = 1000
ps = argparse.ArgumentParser()
ps.add_argument('gene_list', help='the genelist used for permutation')
ps.add_argument('project_name', help='the output project name, will make a dir under ~/gr/')
ps.add_argument(
    'chrin', help='the chrom sequence need to be read in, default = readin all the chromosomes, could be multiple',
    nargs='*', choices=list(chrmap.keys()) + ['all'])
ps.add_argument('--counttype', '-count', '-counttype', dest='counttype', help='count type, by site or by mut',
                default='site', choices=['site', 'mut', 'by_site', 'by_mut', 'bysite', 'bymut'])
ps.add_argument('--site', '-s', '-site', '--s', dest='site', help='use the count site method', action='store_true')
ps.add_argument('--mut', '-m', '-mut', '--m', dest='mut', help='use the count mut method', action='store_true')
ps.add_argument('--permu_cycles', '--n', '-n', '--perm', '-perm', dest='permu_cycles',
                help='permutation times, default = 2500', default=3000)
ps.add_argument('--cadd', '-cadd', dest='cadd',
                help='run permutation with cadd weight added, default = False', action='store_true')
ps.add_argument('--ajust_intron', '-intron', dest='adjust_intron',
                help='use the site count/ mut count after ajusted by intron = False', action='store_true')
ps.add_argument(
    '--dict', '--pdict', '--countpdict', '-countpdict', '-pdict', '-dict', dest='pdict',
    help='specify the variant count pdict, default = /home/chenh19/gr/snpcount/mgrb.all.lite.maf_0.001.count.pdict',
    default='/home/chenh19/gr/snpcount/mgrb.all.lite.maf_0.001.count.pdict')
ps.add_argument('--mut_table', '-mut_table', '-muttab', '-muttable', dest='mut_table',
                help='specify the mutation table, default = coresponding table with pdict in gr/mut_table', default='na')
ps.add_argument(
    '--codon_ratio', '-codon_ratio', '-codonratio', '-codon', dest='codon_ratio',
    help='specify the  codon 3 phase ratio,, default = coresponding table with pdict in gr/mut_table ', default='na')

ps.add_argument('--name_appendix', '-suf', '-suffix', dest='suf',
                help='the  suffix for the final outpu tile', default='')

args = ps.parse_args()

permu_cycles = int(args.permu_cycles)
gene_list = args.gene_list
project_name = args.project_name
counttype = args.counttype
counttype = 'site' if counttype.find('site') > -1 else 'mut'
if args.site:
    counttype = 'site'
elif args.mut:
    counttype = 'mut'

chrin = args.chrin
if chrin == ['all']:
    chrin = list(map(str, range(1, 25)))
    print('you choose to load all the chrom sequence')
    p(chrin)
else:
    chrin = [chrmap[i] for i in chrin]

defaut_mut_table = "/home/chenh19/gr/mut_table/gnomad.genome.maf_0.001.mutation_table.txt"
# defaut_codon_ratio = "/home/chenh19/gr/mut_table/gnomad.gm.control.maf_0.001.codon_ratio.txt"
defaut_codon_ratio = "/home/chenh19/gr/mut_table/gnomad_gm.maf_0.01.codon_ratio_bysite.txt"

pdict_prefix = args.pdict.split(".count")[0].rsplit("/", 1)[-1]

if args.codon_ratio == 'na':
    if os.path.exists("/home/chenh19/gr/mut_table/%s.codon_ratio.txt" % (pdict_prefix)):
        args.codon_ratio = "/home/chenh19/gr/mut_table/%s.codon_ratio.txt" % (pdict_prefix)
    else:
        args.codon_ratio = defaut_codon_ratio


if args.mut_table == 'na':
    if os.path.exists("/home/chenh19/gr/mut_table/%s.mutation_table.txt" % (pdict_prefix)):
        args.mut_table = "/home/chenh19/gr/mut_table/%s.mutation_table.txt" % (pdict_prefix)
    else:
        args.mut_table = defaut_mut_table

d_mutation_table = eval(open(args.mut_table).read())
d_codon_ratio = eval(open(args.codon_ratio).read())

adjust_intron = args.adjust_intron
cadd = args.cadd

af_max_threshold = 100

fn_gn_info = "/home/chenh19/ref/basic/gene_full_info.pdict"
# fn_exac_count = "/home/chenh19/gr/snpcount/mgrb.count.pdict"
fn_exac_count = args.pdict

gene_list_in = set()

with open(gene_list, 'r') as fp:
    for i in fp:
        i = i.strip()
        gene_list_in.add(i)

suf = '' if args.suf == '' else '_%s' % (args.suf)
gene_list_in = [_f for _f in gene_list_in if _f]
short_name = gene_list.rsplit("/", 1)[-1].rsplit(".")[0] + suf
pw = "/home/chenh19/gr/perm/%s" % (project_name)


f_make_dirs(pw, ["log", "pbs", "pdict", "candicate_gene_list"])


error = open("%s/log/%s_error_during_gr_permutation.txt" % (pw, short_name), 'w')
flog = open("%s/log/%s_run_status.txt" % (pw, short_name), 'w')
fl_result_pdict = "%s/pdict/%s.result.pdict" % (pw, short_name)

print('mutation table = %s' % (args.mut_table))
print("codon ratio table = %s" % (args.codon_ratio))
print('mutation table = %s' % (args.mut_table), file=flog)
print("codon ratio table = %s" % (args.codon_ratio), file=flog)
print(args, file=flog)

# fh_exac_pdict = open(fl_exac_pdict,'w')
fh_result_pdict = open(fl_result_pdict, 'w')

conseq_list = ['stop_lost', 'synonymous_variant', 'frameshift', 'missense_variant', 'inframe', 'stop_gained']

# read in the gene full info and exac snp count
exec("d1 = %s" % (open(fn_gn_info).read()))

# load the exac snp count
d_exac_snp_count = eval(open(fn_exac_count).read())
# exec("d_exac_snp_count=%s"%(open(fn_exac_count).read()) )

print("start to load sequence", file=flog)
d_seq = {}
fa = """/home/chenh19/ref/hg19/oneline/chr3.txt
/home/chenh19/ref/hg19/oneline/chr10.txt
/home/chenh19/ref/hg19/oneline/chrY.txt
/home/chenh19/ref/hg19/oneline/chr13.txt
/home/chenh19/ref/hg19/oneline/chr4.txt
/home/chenh19/ref/hg19/oneline/chr14.txt
/home/chenh19/ref/hg19/oneline/chr18.txt
/home/chenh19/ref/hg19/oneline/chr5.txt
/home/chenh19/ref/hg19/oneline/chr21.txt
/home/chenh19/ref/hg19/oneline/chr6.txt
/home/chenh19/ref/hg19/oneline/chr12.txt
/home/chenh19/ref/hg19/oneline/chr16.txt
/home/chenh19/ref/hg19/oneline/chr19.txt
/home/chenh19/ref/hg19/oneline/chr7.txt
/home/chenh19/ref/hg19/oneline/chr22.txt
/home/chenh19/ref/hg19/oneline/chr1.txt
/home/chenh19/ref/hg19/oneline/chr8.txt
/home/chenh19/ref/hg19/oneline/chr20.txt
/home/chenh19/ref/hg19/oneline/chr15.txt
/home/chenh19/ref/hg19/oneline/chr11.txt
/home/chenh19/ref/hg19/oneline/chr2.txt
/home/chenh19/ref/hg19/oneline/chr9.txt
/home/chenh19/ref/hg19/oneline/chrX.txt
/home/chenh19/ref/hg19/oneline/chr17.txt
/home/chenh19/ref/hg19/oneline/chr23.txt
/home/chenh19/ref/hg19/oneline/chr24.txt
/home/chenh19/ref/hg19/oneline/chr25.txt""".split("\n")

for i in fa:
    chr = i.rsplit("/", 1)[-1].split(".")[0]

    for ichr in chrin:
        if chr == ichr:
            with open(i, 'r') as fp:
                print("\tloading %s" % (chr), file=flog)
                exec("%s = fp.read()" % (chr))


print(" sequence loaded", file=flog)
# print(" sequence loaded")
flog.flush()

#  build the exac_lite, which include only the gene in genelist
exac_lite = {}

exac_perm_result = {}

for i in gene_list_in:
    if i in d_exac_snp_count:
        exac_lite[i] = d_exac_snp_count[i]

print("now doing %s" % (gene_list), file=flog)
print("total genes input =  %s" % (len(gene_list_in)), file=flog)
print('there are **%s** genes found in vcf file SNP count' % (len(exac_lite)), file=flog)
not_in_exac = set(gene_list_in) - set(exac_lite)

if len(not_in_exac) > 0:
    print('\n'.join(['%s  ***  not included in SNP count pdict' % (i) for i in not_in_exac]), file=flog)
    print('\n'.join(['%s  ***  not included in SNP count pdict' % (i) for i in not_in_exac]), file=error)

flog.flush()
# {'total_snp_in_cds': 2116.0,
#  'total_sites_in_cds': 80.0,
#  'expected_sites': 79.81,
#  'expected_mut': 1999.84,
#  'cds_conseq_observed_mut': {'synonymous_variant': 743.0,
#   'missense_variant': 1337.0,
#   'stop_gained': 36.0},
#  'cds_conseq_observed_site': {'synonymous_variant': 24.0,
#   'missense_variant': 54.0,
#   'stop_gained': 2.0},
#  'dm': {'PF00207': {'observed': {'synonymous_variant': [36.0, 2.0],
#     'missense_variant': [89.0, 4.0],
#     'stop_gained': [0.0, 0.0]},
#    'pvalue': {'synonymous_variant': 'na',
#     'missense_variant': 'na',
#     'stop_gained': 'na'}},


for gn in sorted(exac_lite.keys()):
    exac_info = exac_lite[gn]
    print("now running %s\t" % (gn), os.popen("date").read(), file=flog)
    # print ("now running %s\t"%(gn),os.popen("date").read())
    if gn not in d1:
        print("gene = %s, not found in gene_full_info" % (gn))
        continue
    # 0=cdsrg, 1=dm, 2='-', '19', 1488, 'P04217', 'NM_130786']
    chr = d1[gn][3]
    chr_seq = "chr%s" % (chr)
    strand = d1[gn][2]
    dm = d1[gn][1]
    cdsrg = d1[gn][0]
    allrg = rg_split(cdsrg, 4, strand)

    if chr_seq not in globals():
        print("gene=%s in wrong chr, now loading %s" % (gn, chr_seq), file=flog)
        tmpfl = "/home/chenh19/ref/hg19/oneline/%s.txt" % (chr_seq)
        with open(tmpfl, 'r') as fp:
            exec("%s = fp.read()" % (chr_seq))
    conseq = 'all'

    # cdslen = getcdslen(cdsrg)
    chr_seq_usable = eval(chr_seq)

    d_cds_phase = {i: getphase(cdsrg, i, strand)[0] for i in allrg}
    # getphase passin = gm position,
    # return = (phase, cds_position), phase is in 1,2,3


    d_all_pos_codon_info = get_all_codon_in_cds(chr_seq_usable, cdsrg, strand)

    if d_all_pos_codon_info[0] == 'error':
        print("gn %s can't get the codon\terror=%s" % (gn, d_all_pos_codon_info[1]), file=error)
        continue
    else:
        d_all_pos_codon_info = d_all_pos_codon_info[1]

    # print (d_all_pos_codon_info[140449117])

    # 'cds_count': [219.6, 71.0, 276.88, 84.29],
    if adjust_intron:
        total_sites = int(exac_info['cds_count'][3])
        total_mut = exac_info['cds_count'][2]
    else:
        total_sites = int(exac_info['cds_count'][1])
        total_mut = exac_info['cds_count'][0]

    d_dm_rg = {idm: rg_split(dm[idm][0], 4, strand, 1) for idm in dm}

    dm_mut = exac_info['dm']

    sum_permutation_conseq = {}
    single_permutation_conseq = {}

    # generate the cadd dict
    #  {'PF00757': [[[55214427, 55214433, 553, 559],
#            [55218987, 55219055, 560, 628],
#            [55220239, 55220357, 629, 747],
#            [55221704, 55221845, 748, 889],
#            [55223523, 55223639, 890, 1006],
    #    [55224226, 55224233, 1007, 1014]],
    d_cadd = {}
    if cadd:
        for idm in dm:
            tmprg = dm[idm][0]
            for k1 in tmprg:
                s, e = k1[0], k1[1]
                chr_cadd = chr if chr not in ['23', '24'] else {'23': 'X', '24': 'Y'}[chr]
                res_cadd = os.popen("tabix /home/chenh19/ref/cadd/whole_genome_SNVs.tsv.gz %s:%s-%s" %
                                    (chr_cadd, s, e)).read().strip().split("\n")
                # 7	140477796	T	A	0.638154	10.59
                for i1 in res_cadd:
                    a1 = i1.split("\t")
                    try:
                        ipos1 = int(a1[1])
                    except:
                        print(f's={s},e={e}, res = {i1}')

                    if ipos1 not in d_cadd:
                        d_cadd[ipos1] = {}

                    d_cadd[ipos1][a1[3]] = float(a1[5])
    # print (len(d_cadd))

    for idm in dm:
        single_permutation_conseq[idm] = {"missense_variant": 0, "synonymous_variant": 0, "stop_gained": 0}
        sum_permutation_conseq[idm] = {"missense_variant": [], "synonymous_variant": [], "stop_gained": []}

    # do the permutation
    res = do_permutation(gn, permu_cycles, conseq, cdsrg, allrg, d_cds_phase, d_all_pos_codon_info, chr, strand,
                         total_sites, total_mut, dm, d_dm_rg, dm_mut, error, single_permutation_conseq,
                         sum_permutation_conseq, counttype, cadd, d_cadd, d_mutation_table, d_codon_ratio)

    if res == 'error':
        print('\terror during permutation\t', gn)
        continue
    exac_perm_result[gn] = res

error.close()
print(exac_perm_result, file=fh_result_pdict)
fh_result_pdict.close()
flog.close()

# exac count info
# {'total_snp_in_cds': 2089.48,
#  'total_sites_in_cds': 559,
#  'cds_conseq_observed_mut': {'synonymous_variant': 1039.91,
#   'missense_variant': 1045.59,
#   'stop_gained': 2.09,
#   'inframe': 0.94,
#   'frameshift': 0.94},
#  'cds_conseq_observed_site': {'synonymous_variant': 220,
#   'missense_variant': 335,
#   'stop_gained': 2,
#   'inframe': 1,
#   'frameshift': 1},
#  'dm': {'PF00757': {'observed': {'synonymous_variant': [177.15, 32],
#     'missense_variant': [121.67, 42],
#     'stop_gained': [0.94, 1],
#     'inframe': 0,
#     'frameshift': 0},
#    'pvalue': {'synonymous_variant': 'na',
#     'missense_variant': 'na',
#     'stop_gained': 'na',
#     'inframe': 'na',
#     'frameshift': 'na'}},
#   'PF01030': {'observed': {'synonymous_variant': [102.7, 34],
#     'missense_variant': [70.65, 46],
#     'stop_gained': 0,
#     'inframe': 0,
#     'frameshift': 0},
#    'pvalue': {'synonymous_variant': 'na',
#     'missense_variant': 'na',
#     'stop_gained': 'na',
#     'inframe': 'na',
#     'frameshift': 'na'}},
#   'PF07714': {'observed': {'synonymous_variant': [360.17, 46],
#     'missense_variant': [274.16, 64],
#     'stop_gained': 0,
#     'inframe': 0,
#     'frameshift': 0},
#    'pvalue': {'synonymous_variant': 'na',
#     'missense_variant': 'na',
#     'stop_gained': 'na',
#     'inframe': 'na',
#     'frameshift': 'na'}},
#   'PF14843': {'observed': {'synonymous_variant': [131.04, 30],
#     'missense_variant': [138.49, 45],
#     'stop_gained': 0,
#     'inframe': 0,
#     'frameshift': 0},
#    'pvalue': {'synonymous_variant': 'na',
#     'missense_variant': 'na',
#     'stop_gained': 'na',
#     'inframe': 'na',
#     'frameshift': 'na'}},
#   'outside_domain': {'observed': {'synonymous_variant': [268.86, 78],
#     'missense_variant': [440.63, 138],
#     'stop_gained': [1.15, 1],
#     'inframe': [0.94, 1],
#     'frameshift': [0.94, 1]},
#    'pvalue': {'synonymous_variant': 'na',
#     'missense_variant': 'na',
#     'stop_gained': 'na',
#     'inframe': 'na',
#     'frameshift': 'na'}}}}


# exac info
# {'total_snp_in_cds': '2089.4778429379794', '
# total_sites_in_cds': '558',
# 'cds_conseq_observed': {'missense_variant': '1045.593556434887', 'stop_gained': '2.0881009764169423', 'synonymous_variant': '1039.9123961924395', 'frameshift': '0.942258404944972', 'inframe': '0.9415309292910272'},

# 'dm': {
# 'PF00757':
# 	{'observed': {'missense_variant': '121.66604408205137', 'stop_gained': '0.9415486592347092', 'synonymous_variant': '177.15093849986323'},
# 	'pvalue': {'missense_variant': 'na', 'stop_gained': 'na', 'synonymous_variant': 'na'}},
# 'PF01030': {'observed': {'missense_variant': '70.64726611150668', 'synonymous_variant': '102.69789329530236'},
# 'pvalue': {'missense_variant': 'na', 'synonymous_variant': 'na'}},

# 'PF07714': {'observed': {'missense_variant': '274.16135027650887', 'synonymous_variant': '360.1661615496443'},
# 'pvalue': {'missense_variant': 'na', 'synonymous_variant': 'na'}},

# 'PF14843': {'observed': {'missense_variant': '138.49112195620094', 'synonymous_variant': '131.03921325613084'},
# 'pvalue': {'missense_variant': 'na', 'synonymous_variant': 'na'}},

# 'outside_domain': {'observed': {'missense_variant': '440.62777400861904', 'stop_gained': '1.1465523171822332', 'synonymous_variant': '268.85818959149873', 'frameshift': '0.942258404944972', 'inframe': '0.9415309292910272'},
# q'pvalue': {'missense_variant': 'na', 'stop_gained': 'na', 'synonymous_variant': 'na', 'frameshift': 'na', 'inframe': 'na'}}}}


# gene_info: (egfr)
# [[[55086971, 55087058, 1, 88],
#   [55209979, 55210130, 89, 240],
#   [55210998, 55211181, 241, 424],
#   [55214299, 55214433, 425, 559],
#   [55218987, 55219055, 560, 628],
#   [55220239, 55220357, 629, 747],
#   [55221704, 55221845, 748, 889],
#   [55223523, 55223639, 890, 1006],
#   [55224226, 55224352, 1007, 1133],
#   [55224452, 55224525, 1134, 1207],
#   [55225356, 55225446, 1208, 1298],
#   [55227832, 55228031, 1299, 1498],
#   [55229192, 55229324, 1499, 1631],
#   [55231426, 55231516, 1632, 1722],
#   [55232973, 55233130, 1723, 1880],
#   [55238868, 55238906, 1881, 1919],
#   [55240676, 55240817, 1920, 2061],
#   [55241614, 55241736, 2062, 2184],
#   [55242415, 55242513, 2185, 2283],
#   [55248986, 55249171, 2284, 2469],
#   [55259412, 55259567, 2470, 2625],
#   [55260459, 55260534, 2626, 2701],
#   [55266410, 55266556, 2702, 2848],
#   [55268009, 55268106, 2849, 2946],
#   [55268881, 55269048, 2947, 3114],
#   [55269428, 55269475, 3115, 3162],
#   [55270210, 55270318, 3163, 3271],
#   [55272949, 55273310, 3272, 3633]],
#  {'PF00757': [[[55214427, 55214433, 553, 559],
#     [55218987, 55219055, 560, 628],
#     [55220239, 55220357, 629, 747],
#     [55221704, 55221845, 748, 889],
#     [55223523, 55223639, 890, 1006],
#     [55224226, 55224233, 1007, 1014]],
#    'Furin-like'],
#   'PF01030': [[[55210059, 55210130],
#     [55210998, 55211181],
#     [55214299, 55214375],
#     [55224300, 55224352],
#     [55224452, 55224525],
#     [55225356, 55225446],
#     [55227832, 55227973]],
#    'Recep_L_domain'],
#   'PF07714': [[[55241692, 55241736, 2140, 2184],
#     [55242415, 55242513, 2185, 2283],
#     [55248986, 55249171, 2284, 2469],
#     [55259412, 55259567, 2470, 2625],
#     [55260459, 55260534, 2626, 2701],
#     [55266410, 55266556, 2702, 2848],
#     [55268009, 55268055, 2849, 2895]],
#    'Pkinase_Tyr'],
#   'PF14843': [[[55229206, 55229324, 1513, 1631],
#     [55231426, 55231516, 1632, 1722],
#     [55232973, 55233130, 1723, 1880],
#     [55238868, 55238895, 1881, 1908]],
#    'GF_recep_IV']},
#  '+',
#  '7',
#  3633,
#  'P00533',
#  '',
#  'NM_005228']
