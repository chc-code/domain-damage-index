# previously, the total mut count of the gene used for permutation = actual total mut count.
# but for those which the ratio of syno/missense largely devirate from the average ratio,
# both the missense and syno zscore would be deviate from 0
# here we assume that syno mutation is under the average selection pressure, 
# so we used the  n_syno_mut / average_syno_ratio  as the total mutation count
# the averge syno_mut ratio = 1/3

import sys,os,re,gzip
home=os.popen("echo ~").read().strip()
sys.path.append('%s/jb/module'%(home))

import pandas as pd;import numpy as np;
import matplotlib.pyplot as plt
import sklearn.metrics as sk
import seaborn as sns

# the exac is the one permuated with the new  total mut count


pw = "%s/generisk"%(home)



# load the compare

fexac_zscore = "%s/exac/exac_zscore.lite_gn_checked.txt"%(pw)
fgdi = "%s/gdi/gdi_score_pnas_gn_checked.txt"%(pw)
fperm_z = "%s/exac/exac_201907_adjusted_gene_just_site_count_perm_result.tsv"%(pw)
ct = eval(open('%s/exac/exac_count_final.pdict'%(pw)).read())
omim = pd.read_csv("/Users/zeiss/gr/omim/omim_dedup.tsv",sep='\t')


exacz = pd.read_csv(fexac_zscore,sep='\t')
exacz=exacz[['gn','exac_z']]
exacz['exac_zrank'] = round(exacz['exac_z'].rank()/len(exacz.index)*100,2)

gdi = pd.read_csv(fgdi,sep="\t")
gdi = gdi[['gn','gdi','gdi_phred_raw']]
gdi['gdi_rank'] = 100- round(gdi['gdi'].rank()/len(gdi.index)*100,2)

z = pd.read_csv(fperm_z,sep='\t')
z = z.loc[z.conseq=='missense_variant']
z['zrank'] = round(z['zscore'].rank()/len(z.index)*100,2)



# z.describe()
# zscore    perm_mean    perm_std    observ
# count    45043.000000    45043.000000    45043.000000    45043.000000
# mean    -0.442758    283.086033    27.482658    260.533714
# std    6.828854    424.159458    16.187950    431.084395
# min    -443.330000    0.000000    0.000000    0.000000
# 25%    -2.360000    74.315000    16.510000    49.890000
# 50%    0.910000    160.110000    24.510000    138.780000
# 75%    2.950000    332.030000    35.090000    313.075000
# max    33.900000    14648.650000    213.960000    18941.130000

# print(z.shape)
# print(gdi.shape)
# print(exacz.shape)

# (45043, 8)
# (18492, 4)
# (17874, 3)

z1 = z.loc[z['dm'] !='outside_domain']  #27849
z1['zrank_dm'] = round(z1['zscore'].rank()/len(z1.index) *100,2)

dm1 = pd.merge(z1,exacz,on='gn')
ddi = pd.merge(dm1,gdi,on='gn')


# negative set 
neg = pd.read_csv("/work/ref/exac/nontcga_filtered.halfclose.annotated.refined.map_to_domain.bed",sep='\t',header=None)
neg= neg.loc[neg[9] =='missense_variant']
neg = neg[[13,14]].rename(columns={13:'gn',14:'dm'})
neg_ddi = pd.merge(neg, ddi, on=['gn','dm'])


# positive set clinvar

clinvar = pd.read_csv("/Users/zeiss/gr/clinvar/clinvar_map_to_domain.tsv",sep='\t',header=None)
clinvar = clinvar[[5,11,12]].rename(columns={5:'sig',11:'gn',12:'dm'})
snp = clinvar
snp_ddi = pd.merge(snp, ddi,on=['gn','dm'])
snp_ddi_omim = pd.merge(snp_ddi, omim, on = 'gn')


# find best way to combine the gene/ domain zscore


pos1 = snp_ddi_omim.loc[(snp_ddi_omim['sig']=='Pathogenic') & (snp_ddi_omim['inher'] !='AR'),['gn','dm','gdi_rank','zrank','exac_zrank']]
neg1 = neg_ddi[['gn','dm','gdi_rank','zrank','exac_zrank']]


pos1= pos1.loc[pos1['dm']!= 'outside_domain']
neg1= neg1.loc[neg1['dm']!= 'outside_domain']


pos1['flag'] = 1
neg1['flag'] = 0

plt.plot([0,1],[0,1],color='navy',linestyle='--')
plt.xlabel('FPR')
plt.ylabel('TPR')

plt.title("ROC Curve clinvar-AD/X-link_gene + exac_non_tcga")
roc1 = pd.concat([pos1,neg1])

roc1['new'] = pd.np.nan


# base
fpr,tpr,thres = sk.roc_curve(roc1['flag'],roc1['exac_zrank'])
auc2 = sk.auc(fpr,tpr)
plt.plot(fpr,tpr,color ='navy',label = "exac_z AUC = %.4f"% auc2 )


fpr,tpr,thres = sk.roc_curve(roc1['flag'],roc1['zrank'])
auc2 = sk.auc(fpr,tpr)
plt.plot(fpr,tpr,color ='red',label = "zrank AUC = %.3f"% auc2 )


# begin combin

left,mid,mid2,right = 10,30,70,90
ratio1 = (0.2,0.8)  # left - mid, zrank_exac percent, zrank percent
ratio3 = (0.5,0.5)  # mid - mid1 zrank_exac percent, zrank percent
ratio2 = (0.8,0.2)  # mid1 - right, zrank_exac percent, zrank percent
#rs = 0.


# # 0-left  = zrank
# # left - mid  = ratio1
# # mid - mid2 = min(exac_z, zrank)
# # mid2 - right  = ratio2
# # right - 100 = exac_z


roc1.loc[(roc1['exac_zrank']<left),'new'] = roc1.loc[(roc1['exac_zrank']<left),'zrank']
roc1.loc[(roc1['exac_zrank']>right),'new'] = roc1.loc[(roc1['exac_zrank']>right),'exac_zrank']



roc1.loc[(roc1['exac_zrank']>=left) & (roc1['exac_zrank']<=mid)  ,'new'] = roc1.loc[(roc1['exac_zrank']>=left) & (roc1['exac_zrank']<=mid),'exac_zrank']*ratio1[0] +roc1.loc[(roc1['exac_zrank']>=left) & (roc1['exac_zrank']<=mid),'zrank']*ratio1[1]
# roc1.loc[(roc1['exac_zrank']>mid) & (roc1['exac_zrank']<=mid2)  ,'new'] = \
# roc1.loc[(roc1['exac_zrank']>mid) & (roc1['exac_zrank']<=mid2),['exac_zrank','zrank_adj']].max(axis=1)

roc1.loc[(roc1['exac_zrank']>mid) & (roc1['exac_zrank']<=mid2)  ,'new'] = \
    roc1.loc[(roc1['exac_zrank']>mid) & (roc1['exac_zrank']<=mid2),'exac_zrank'] * ratio3[0]  + \
    roc1.loc[(roc1['exac_zrank']>mid) & (roc1['exac_zrank']<=mid2),'zrank'] * ratio3[1]


roc1.loc[(roc1['exac_zrank']>mid2) & (roc1['exac_zrank']<=right)  ,'new'] = roc1.loc[(roc1['exac_zrank']>mid2) & (roc1['exac_zrank']<=right),'exac_zrank']*ratio2[0] + roc1.loc[(roc1['exac_zrank']>mid2) & (roc1['exac_zrank']<=right),'zrank']  *ratio2[1] 

roc2 = roc1.loc[pd.notna(roc1['new'])]

print(roc1.shape,roc2.shape)


fpr,tpr,thres = sk.roc_curve(roc2['flag'],roc2['new'])
auc2 = sk.auc(fpr,tpr)
plt.plot(fpr,tpr,color ='pink',label = "new%s_%s AUC = %.4f"% (left,right,auc2) )


print (round(auc2,4))

plt.legend(loc="lower right") 
