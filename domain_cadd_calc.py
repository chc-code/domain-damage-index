# extract all the bases in the domain
# for each base, there are 3 values.  we just average all the values
# e.g  chr1:101-200,   200bp, the result would be 600 lines, and we just get the raw_cadd score , and sum it, and then divide it by 600, to be the score for this domain
# didn't consider the maf ...


# !!!!!! the output result, total CADD score count is 3 times of the real domain length.  eg, count = 999, the real domain len = 333. but the average stays the same
import os
import sys
import gzip
# import pandas as pd


import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('gnlist', help="""only run the given genelist, otherwize, would run all the genes. \nthis is very useful when you need to run this script parallelly""", nargs='*', default=None)
ps.add_argument('-listgene', '-genelist', help="""export all genename to file genelist.txt, used for parallel""", action='store_true')
ps.add_argument('-suf', '-suffix', dest='suffix', help="""add the suffix to the output file, used when run in parallel, avoid the file write confilict, usaully set as the thread id""", default=None)
args = ps.parse_args()

if args.listgene:
    info = eval(open("/home/chenh19/ref/basic/gene_full_info.pdict").read())
    with open('/home/chenh19/ref/cadd/genelist.txt', 'w') as out:
        print('\n'.join(info.keys()), file=out)
    print('export genelist done')
    sys.exit()

column_to_sum = 6
fn = "/home/chenh19/ref/cadd/whole_genome_SNVs.tsv.gz"
fnout = "/home/chenh19/ref/cadd/all_domain_casdd_score.tsv"
try:
    fnout += args.suffix
except:
    pass


# Chr	Pos	Ref	Alt	RawScore	PHRED
# 1	10001	T	A	0.118631	4.575

if not os.path.exists(fnout):
    out = open(fnout, 'w')
    print('gn,dm, dmname,cadd_sum,cadd_query_count'.replace(",", '\t'), file=out)
else:
    out = open(fnout, 'a')

gnlist = args.gnlist or list('abcd')
# print(gnlist, file=out)
# print('end', file=out)
# sys.exit()

info = eval(open("/home/chenh19/ref/basic/gene_full_info.pdict").read())
for gn in gnlist:
    try:
        v1 = info[gn]
    except:
        print(f'gn not found {gn}')
        continue
    dm = v1[1]
    chr_ = v1[3].replace("chr", '')
    for v2 in dm.values():
        rg = v2[0]
        dmname = v2[3]
        pid = v2[1]
        if pid == 'outdm':
            continue
        sum1, n1 = 0, 0
        for irg in rg:
            # export the sum of cadd and count of lines
            # if the domain len = n, the cadd query result = 3n
            # each irg, would transfer the sum1, n1 to the awk, so  the total sum / number in all irg would be added
            tmp = os.popen(
                "tabix %s %s:%s-%s|awk 'BEGIN{sum=0;sum1=%s;n1=%s}{sum+=$%s}END{print sum+sum1,n1+NR}' 2>/dev/null" %
                (fn, chr_, irg[0], irg[1],
                 sum1, n1, column_to_sum)).read().strip().split(" ")
            if len(tmp) != 2:
                print('error get domain small range value:\t%s\t%s\t%s' % (gn, pid, irg))
                continue
            sum1, n1 = tmp[0], tmp[1]

        print("%s\t%s\t%s\t%s\t%s" % (gn, pid, dmname, sum1, n1), file=out, flush=True)

out.close()
