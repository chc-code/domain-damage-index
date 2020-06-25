import os,sys,re
import pandas as pd
home = os.path.expanduser('~')
sys.path.append("%s/jb/module"%(home))

from chctool import newfn
fnraw = sys.argv[1]
pw,fn,ext = newfn(fnraw)
fn = fn+'.'+ext if ext else fn



m = re.match('(.+).in_cds.+anno.refined.cadd.txt.gz',fn)
if m:
    lb = m.groups()[0]
else:
    print ('bad file name pattern')
    sys.exit()




anno = pd.read_csv(fnraw,sep='\t')
anno['s'] = anno['s'] -1
anno['chr,s,e,ref,alt,af,conseq,gn,cadd'.split(",")].to_csv("{0}{1}.incds.anno.cadd.bed.gz".format(pw,lb),index=False,header=None,na_rep='NA',sep='\t')

os.system("intersectBed -a {0}{1}.incds.anno.cadd.bed.gz -b {2}/ref/basic/gene_full_info_domain_range.bed -wb|cut -f 1-7,9,13-15|  gzip >{0}{1}.in_domain.anno.cadd.bed.gz".format(pw,lb,home))




