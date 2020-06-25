# all the tracks comes from 

# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/UP000005640_9606_beds//


# file lookes like
# 1    chr1
# 2    6825162
# 3    7640453
# 4    Q9Y6Y1
# 5    0
# 6    +
# 7    6825162
# 8    7640453
# 9    0,153,51
# 10    1
# 11    815291
# 12    0
# 13    .
# 14    K63-C188; CG-1

# the bed file is 0-index, half-closed

import os,sys,re

pwref = "/work/ref"
fn=sys.argv[1]

sys.path.append("/jb/module")
from chctool import map_aa_to_gm


ts  = eval(open("%s/pdict/all_transcript_cds_range.pdict"%(pwref)).read())
# {'NR_024540': ['WASH7P', '1', '-', [[14362, 14829, 1769, 1302], [14970, 15038, 1301, 1233], [15796, 15947, 1232, 1081], [16607, 16765, 1080, 922], [16858, 17055, 921, 724], [17233, 17368, 723, 588], [17606, 17742, 587, 451], [17915, 18061, 450, 304], [18268, 18366, 303, 205], [24738, 24891, 204, 51], [29321, 29370, 50, 1]], 
# '585\tNR_024540\tchr1\t-\t14361\t29370\t29370\t29370\t11\t14361,14969,15795,16606,16857,17232,17605,17914,18267,24737,29320,\t14829,15038,15947,16765,17055,17368,17742,18061,18366,24891,29370,\t0\tWASH7P\tunk\tunk\t-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,']

uid2ts  = eval(open("%s/pdict/canonical_uniprot_id_info_include_transcript_isoform.pdict"%(pwref)).read())
# 'P31946': {'ts': ['NM_003404', 'NM_139323', 'XM_017028039', 'ENST00000353703', 'ENST00000372839'], 'primary_uid': 'P31946', 'isoform': 'P31946-1', 'all_uid': ['P31946', 'A8K9K2', 'E1P616'], 'uacc': '1433B_HUMAN'},

print ('dict loaded')


fn = fn.split(",")
for ifn in fn:
    ifn = ifn.strip()
    lb = ifn.rsplit("/",1)[-1].rsplit(".",1)[0]

    out = open("%s_hg19_pos_half_close.bed"%(lb),'w')
    err = open("err_%s_get_hg19_pos.txt"%(lb),'w')
    uidnotfound = set()
    uidnots = set()
    wrongts  = set()

    with open(ifn) as fp:
        print ("running %s"%(ifn))
        for i in fp:
            i=i.strip()
            a = i.split("\t")
            chr,uid,strand,pos = a[0],a[3],a[5],a[13]
            pos = pos.split(";")[0].strip()
            pos = pos.split(",")
            posnew = []
            # D111-E122,D150-E161,D246-E257,D291-E302,D327-E338; calcium-binding region_111-122_150-161_246-257_291-302_327-338
            for ipos in pos: 
                ipos=ipos.strip()
                ipos = ipos.split("-")
                ipos = [re.sub("[A-Za-z*]+",'',_) for _ in ipos]

                try:
                    ipos = [int(_) for _ in ipos]
                except:
                    print ("invalid aa pos,",uid,ipos, a[13],file=err)    
                    continue

                if len(ipos ) ==1:
                    posnew.append(ipos[0])
                elif len(ipos) >2:
                    print ('bad, wrong aa range',uid,a[13],file=err)
                    continue
                else:
                    posnew.append(ipos)    

            if uid not in uid2ts:
                uidnotfound.add(uid)
            else:
                res1 = uid2ts[uid]
                count = 0
                for its in res1['ts']:
                    if its in ts:
                        count = 1
                        res = ts[its]
                        gn,chr1,cdsrg,strand1 = res[0],res[1],res[3],res[2]
                        break
                if count == 0:
                    uidnots.add(uid)
                    continue

                for ipos in posnew:    
                    gm_pos = map_aa_to_gm(cdsrg, strand1, ipos , 'aa')
                    if gm_pos[0] =='error':
                        wrongts.add("%s\t%s"%(uid,its))
                        continue

                    for irg in gm_pos:
                        s = min(irg[0],irg[1])-1   # the bed file is 0-index, half-closed    
                        e = max((irg[0],irg[1]))
                        print ('%s\t%s\t%s\t%s\t%s\t%s\t%s'%(chr1,s,e,strand,gn,ipos,gm_pos),file=out)    

    out.close()                
    if len(uidnotfound)>0:
        print ("uid not found",file=err)
        print("\n".join(uidnotfound),file=err)
        print ("*********\n",file=err)
    if len(uidnots)>0:
        print ("transcript not found",file=err)
        print("\n".join(uidnots),file=err)
        print ("*********\n",file=err)
    if len(wrongts)>0:
        print ("wrong,transcript requested > max ts len",file=err)
        print("\n".join(wrongts),file=err)
        print ("*********\n",file=err)        

















