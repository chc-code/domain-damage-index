#%% 
import pandas as pd
import os,sys,gzip,re
home=os.path.expanduser('~')
sys.path.append('%s/jb/module'%(home))
# import imp
# import chctool
# imp.reload(chctool)

from chctool import newfn
from chctool import addfn
from chctool import intersect
from chctool import intersectsimple

import argparse
ps=argparse.ArgumentParser('prepare positive set')
ps.add_argument('fn', help = 'filename of the postive set, should be like clinvar.in_domain.anno.cadd.bed.gz(by intersect), or .in_domain.anno.txt.gz (by pandas)')
ps.add_argument('-gnomad','-gmaf','-gm',dest='gm', help = 'add the gnomad af or not',action='store_true')
args = ps.parse_args()

fn = args.fn
revel = '/work/ref/revel/revel_all_chromosomes.csv.gz'

subrvis = '/Users/zeiss/ref/subrvis/subrvis_domain_value.lite.bed'
ccr = '/work/ref/ccr/ccr.all.lite.bed.gz'
mtr = '/work/ref/mtr/mtr.all.lite.bed.gz'  # just read to pandas
gnomad = '/work/ref/gnomad/gnomad.exome.txt.gz'

# fn = '/work/gr/roc/1kg.in_domain.anno.cadd.bed.gz'
tmp = newfn(fn)
try:
    os.mkdir("/Users/files/nogit/generisk/roc/roc_dataset/")
except:
    pass    

fn_prefix = "/Users/files/nogit/generisk/roc/roc_dataset/"+re.sub("\.(bed|txt)",'',tmp[1])


# dedup

print ('step1 , dedup')

fn_dedup = fn_prefix +".missense.dedup.pk"

if os.path.exists(fn_dedup):
#       index  0          1          2  3  4         5                 6      7               dm
# 0   0  1   19203991   19203992  G  A  0.000096  missense_variant  25.60  ALDH4A1:PF00171
    d1 = pd.read_pickle(fn_dedup)
else:
    if fn.find(".anno.cadd.bed.gz")>0:
        d = pd.read_csv(fn,sep='\t',header=None)
        # 1       949421  949422  G       A       Benign  missense_variant        9.188   ISG15   PF00240 ubiquitin
        d = d.loc[d[6] =='missense_variant']
        s = d.groupby(list(range(5))).apply(lambda x:x[8].str.cat(x[9],sep=':').str.cat(sep="*")).rename('dm')

        d1 = d.drop([8,9,10],axis=1).drop_duplicates(subset=range(5)).join(s,on=[0,1,2,3,4]).reset_index()
        d1.to_pickle(fn_dedup)
    elif fn.find(".anno.txt.gz")> 0:
        # with header!
        # 1  949422  949422   G   A                  Benign  ISG15  PF00240  ubiquitin  missense_variant   9.188    
        d = pd.read_csv(fn,sep='\t')
        d = d.loc[d['conseq']== 'missense_variant' ] 
        d.columns=[0,1,2,3,4,5,8,9,10,6,7]
        s = d.groupby(list(range(5))).apply(lambda x:x[8].str.cat(x[9],sep=':').str.cat(sep="*")).rename('dm')
        d1 = d.drop([8,9,10],axis=1).drop_duplicates(subset=range(5)).join(s,on=[0,1,2,3,4]).reset_index()
        d1.to_pickle(fn_dedup)

# revel
print ('step2 , revel')
fn_pre_revel = fn_prefix + ".for_revel.gz"
fn_revel = fn_prefix + ".revel.gz"
if os.path.exists(fn_revel):
    d_revel = pd.read_csv(fn_revel, sep='\t',header=None)
else:
    d1[[0,2,3,4,'index']].to_csv(fn_pre_revel,sep='\t',header=None,index=False)

    d_revel = intersect(fn_pre_revel, revel ,colfn = [2,3,4],coldb=[2,3,6] )
    # d_revel.columns=[0,'s',1,2,3,4,5]
    d_revel = d_revel.sort_values([6]).drop_duplicates(subset=[0,1,2,3],keep='last').sort_index()
    # d_revel.to_pickle(fn_revel)
    d_revel.to_csv(fn_revel,sep='\t',  index=False, header=None, na_rep='NA')
# MTR
print ('step3 , MTR')
fn_mtr = fn_prefix + ".revel.mtr.bed.gz"
# fn_mtr_pk = fn_prefix + ".revel.mtr.pk"
if os.path.exists(fn_mtr):
    dm_mtr =  pd.read_csv(fn_mtr, sep='\t',header=None)
else:
    dm_mtr = intersectsimple(fn_revel,mtr,colfn=[0,2,3,4,5,6],coldb=[0,1,2])
    # remove duplicate
    dm_mtr = dm_mtr.sort_values(['2_y']).drop_duplicates(subset=[0,2,3,4],keep='last').sort_index()

    # dm_mtr.to_pickle(fn_mtr_pk)
    dm_mtr.to_csv(fn_mtr,sep='\t',  index=False, header=None, na_rep='NA')

# print (dm_mtr.head())


# subrvis
print ('step4 , subRVIS')
fn_subrvis = fn_prefix + ".revel.mtr.subrvis.bed.gz"
if not os.path.exists(fn_subrvis):
    dm_subrvis = pd.read_csv(os.popen("intersectBed -a %s -b %s -wao|cut -f 1-8,12 "%(fn_mtr,subrvis)),sep='\t',header=None,na_values=['.','NA','na']).sort_values(8).drop_duplicates(subset=5,keep='last').sort_index()

    dm_subrvis.to_csv(fn_subrvis,sep='\t',header=None,index=False,na_rep='NA')
else:
    dm_subrvis = pd.read_csv(fn_subrvis,sep='\t',header=None)

# print (dm_subrvis.head())

# now got 3 values, revel, mtr, subrvis
#%%
# ccr
print ('step5 , CCR')
fn_ccr = fn_prefix + ".revel.mtr.subrvis.ccr.bed.gz"

if not os.path.exists(fn_ccr):
    dm_ccr = pd.read_csv(os.popen("intersectBed -a %s -b %s -wao|cut -f 1-9,13"%(fn_subrvis, ccr)),sep='\t',header=None,na_values=['.','NA','na']).sort_values(9).drop_duplicates(subset=5,keep='last').sort_index()

    dm_ccr.to_csv(fn_ccr,sep='\t',header=None,index=False,na_rep='NA')
else:
    dm_ccr= pd.read_csv(fn_ccr,sep='\t',  header=None)

# print (dm_ccr.head())


# add gnomad AF
# some dataset may not need this step

if args.gm :
    print ('step6 , add gnomad')
    fn_gm = fn_prefix + ".revel.mtr.subrvis.ccr.gnomad_af.bed.pk"

    if os.path.exists(fn_gm):
        dm_gnomad = pd.read_pickle(fn_gm)
    else:
        dm_gnomad = intersect(fn_ccr, gnomad ,colfn = [3,4,5,6,7,8,9],coldb=[2,3,4],inputbed = True )
        # print (dm_gnomad.head())

        # dm_gnomad.columns = [0,'s',2,3,4,5,6,7,8,9,10]
        dm_gnomad = dm_gnomad.sort_values(['4_y']).drop_duplicates(subset=[0,2,3,4],keep='last').sort_index()
        dm_gnomad.to_pickle(fn_gm)
    dm_pre_final = dm_gnomad
else:
    dm_pre_final = dm_ccr


#  release the domain/gn
print ('step6 , final')
if args.gm :
    fn_final = fn_prefix + ".revel.mtr.subrvis.ccr.gnomad.final.bed.pk"
    fn_final_bed = fn_prefix + ".revel.mtr.subrvis.ccr.gnomad.final.bed.gz"
else:
    fn_final = fn_prefix + ".revel.mtr.subrvis.ccr.final.bed.pk"
    fn_final_bed = fn_prefix + ".revel.mtr.subrvis.ccr.final.bed.gz"

if not os.path.exists(fn_final):
    d2 = pd.merge(d1,dm_pre_final.loc[:,dm_pre_final.columns[5:]],left_on=['index'],right_on=[5]).drop(['index','5_y'],axis=1)

    # print (d2.head())
    s = d2['dm'].str.split("*",expand=True).stack().reset_index(level=1,drop=True).str.split(":",expand=True)
    s.columns=['gn','dm']
    d2 = d2.drop('dm',axis=1)
    d2 = d2.join(s)

    col_end = ['gnomad_maf' ,'gn','dm'] if args.gm else ['gn','dm']
    d2.columns = 'chr,s,e,ref,alt,sig,conseq,cadd,revel,mtr,subrvis,ccr'.split(",")+col_end
    d2.to_pickle(fn_final)
    d2 . to_csv(fn_final_bed,index=False, sep='\t')
