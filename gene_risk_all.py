# import package
import sys
import os
import pandas as pd
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib

home = os.path.expanduser('~')
sys.path.append('%s/jb/module' % (home))

from chctools import chrmap
chrmap = chrmap()
pw_final = f'{home}/gr/final'

# define

pathlist = set(['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic'])
pathlist1 = set(['Pathogenic'])
blist = ['Likely_benign', 'Benign', 'Benign/Likely_benign']


# load data
ginfo = eval(open('%s/ref/basic/gene_full_info.pdict' % (home)).read())
d_count = eval(open('%s/gr/snpcount/gnomad.genome.maf_0.001.count.adj_intron.pdict' % (home)).read())
pid = pd.read_csv('/Users/zeiss/ref/basic/pfamDesc.txt.gz', sep='\t',   # domain full name
                  header=None, names=['dm', 'dmname', 'dm_fullname'])
# zscore
perm_file = '/Users/zeiss/gr/perm_result/gnomad_gm_site_intron_nocadd_perm_result.tsv'
z = pd.read_pickle('/ref/basic/gene_risk_zscore_gnomad_gm_site_intron_no_cadd.pk')

# gene level index, out idx , omim
out_idx = pd.read_pickle(f'{home}/gr/final/out_idx.pk')

# clinvar
# (3982, 22)
clinvar_count = pd.read_csv("/Users/zeiss/work/gr/final/clinvar_enrichment_for_heatscatter.csv")

# 74.3% missense pathogenic variants are in domain  (25793/34726)
cling = pd.read_pickle(f'{home}/ref/clinvar/clinvar_cds.pk')  # (274964, 7)
indm = pd.read_pickle(f'{home}/ref/clinvar/clinvar_domain.pk')  # (143972, 8)
cling_path = cling.loc[cling.sig.isin(pathlist)]  # (34726, 7)
indm_path = indm.loc[indm.sig.isin(pathlist)]  # (25793, 8)
#


# how many genes in a certain domain family
zraw = z
dm_type = z.groupby('dm').apply(lambda x: len(x.index)).reset_index()
x = len(dm_type.loc[dm_type[0] == 1])
y = x / len(dm_type)
print('domain in only 1 gene', x, y)   # 0.4924

x = len(dm_type.loc[dm_type[0] == 2])
y = x / len(dm_type)
print('domain in only 2 gene', x, y)  # 0.1798

x = len(dm_type.loc[dm_type[0] > 2])
y = x / len(dm_type)
print('domain in more than gene', x, y)  # 0.3277

# domain in only 1 gene 3129 0.49244570349386213
# domain in only 2 gene 1143 0.17988668555240794
# domain in more than gene 2082 0.3276676109537299


# OR >1   gnomad vs pathogenic enrichment in domain
n_or_gt_1 = len(clinvar_count.loc[clinvar_count['or'] > 1])
print(f'OR >1   gnomad vs pathogenic enrichment in domain {n_or_gt_1}, \nratio = {n_or_gt_1/len(clinvar_count)}')
# OR >1   gnomad vs pathogenic enrichment in domain 2025,
# ratio = 0.7596685082872928

# show the OR distribution
max_or = max(clinvar_count['or'])
# bins = pd.IntervalIndex.from_tuples([(0,0.1),(0.1,0.5),(0.5,1),(1,2),(2,10),(10,100),(100,max_or)],  closed='left')
clinvar_count['or_bin'] = pd.cut(clinvar_count['or'],[0, 0.1, 0.5, 1, 2, 10, 100, max_or+1], labels='0-0.1,0.1-0.5,0.5-1,1-2,2-10,10-100,>100'.spli(","), right=False)

tmp = clinvar_count.groupby('or_bin').apply(lambda x:len(x.index)).reset_index()

plt.figure(figsize=(5,5),dpi=200)
plt.bar(tmp.or_bin, tmp[0], capsize=2)
plt.tick_params(axis='y', labelsize=14)
plt.tick_params(axis='x', labelsize=13, labelrotation=30 )

fontdict = {'fontweight':'bold', 'size':14}
plt.ylabel("Domain Count", fontdict=fontdict)
plt.xlabel("Odds Ratio bin", fontdict=fontdict)
plt.savefig('/Users/files/work/generisk/final/domain_odds_ratio_distribution.png')



# load cosmic patho count
cos_path_ct = pd.read_csv("/Users/zeiss/work/gr/final/cosmic_enrichment_for_heatscatter.csv")



# how many genes, the pathogenic is exclusively in domains
total_gene_count = len(clinvar_count.gn.unique())  # 2868
r100 = clinvar_count.loc[clinvar_count.ratio_patho > 0.99]   # 1237
r80 = clinvar_count.loc[clinvar_count.ratio_patho > 0.8]   # 1463

print(f'{len(r100)} genes , ratio= {len(r100) / total_gene_count}, pathogenic variants are exclusively in 1 domain')
print(f'{len(r80)} genes , ratio= {len(r80) / total_gene_count}, > 80 pathogenic variants are in 1 domain')

# 1237 genes , ratio= 0.4313110181311018, pathogenic variants are exclusively in 1 domain
# 1463 genes , ratio= 0.5101115760111576, > 80 pathogenic variants are in 1 domain


# plot - mos constrained domains across genes

# load the most_least_constrained_domains
d_constrained = pd.read_pickle('/Users/zeiss/work/gr/final/constrained_domains_across_gene_pvalue.pk')

set_most_constrained = set(d_constrained.loc[d_constrained['utest-less'] < 2.5e-5, 'dm'])  # 25
set_least_constrained = set(d_constrained.loc[d_constrained['utest-more'] < 2.5e-5, 'dm'])  # 16

# the actual zrank value
most_constrained = z.loc[z.dm.isin(set_most_constrained)].sort_values('zscore_median')[['dmname', 'zscore', 'zrank']]
least_constrained = z.loc[z.dm.isin(set_least_constrained)].sort_values('zscore_median')[['dmname', 'zscore', 'zrank']]

# plot most constrained domains
plt.figure(figsize=(8, 4), dpi=200)
# print(most_constrained.dmname.value_counts())
sns.boxplot(x='dmname', y='zscore', data=most_constrained, color='silver')
plt.tick_params(axis='y', labelsize=14)
plt.tick_params(axis='x', labelsize=10, labelrotation=45)
fontdict = {'fontweight': 'bold', 'size': 14}
plt.ylabel("DDI", fontdict=fontdict)
plt.xlabel("Domain Name", fontdict=fontdict)
plt.savefig(f'{pw_final}/box_plot_most_constrained_domains_across_genes.png')


# plot least constrained domains
plt.figure(figsize=(8, 2.5), dpi=200)
# print(most_constrained.dmname.value_counts())
sns.boxplot(x='dmname', y='zrank', data=least_constrained, color='silver')
plt.tick_params(axis='y', labelsize=14)
plt.tick_params(axis='x', labelsize=10, labelrotation=45)
fontdict = {'fontweight': 'bold', 'size': 12}
plt.ylabel("DDI percentile", fontdict=fontdict)
plt.xlabel(" ", fontdict=fontdict)
plt.savefig(f'{pw_final}/box_plot_least_constrained_domains_across_genes.png')

# most least constrained domains across gene





# gnomad variants distribution in CACNA1A
gene_full_info = eval(open('/Users/files/ref_basic/gene_full_info.pdict').read())
gnomad = pd.read_csv('/Users/files/ref/gnomad/gnomad.genome.in_cds.maf_0.001.anno.refined.cadd.txt.gz', sep='\t')
domain_aa_range = eval(open('/Users/files/ref_basic/domain_aa_range.pdict').read())

from matplotlib import rcParams, rc


gn = 'CACNA1A'
gene_info = gene_full_info[gn]
ca = gnomad.loc[(gnomad.conseq=='missense_variant') & (gnomad.gn=='CACNA1A')]
ca['aa_pos'] = ca['s'].apply(lambda x: genomic_to_aa(x, gene_info))

uid = gene_info[5]
try:
    domain_aa_info = domain_aa_range[uid]
    domain_aa_info = sorted(domain_aa_info, key=lambda x:x[2])
except:
    print(f'this gene\' uid is not found in pfam  {gn}, {uid}')


total_aa = int(gene_info[4] / 3)
bin_range = list(range(0, int(total_aa) + 1, 10))
if bin_range[-1] < total_aa:
    bin_range.append(total_aa)
ca['aa_bin'] = pd.cut(ca.aa_pos, bin_range, labels=bin_range[1:])

ca.to_csv(f'{gn}.gnomad.missense.csv', index=False)

tmp = ca.groupby('aa_bin').apply(lambda x:len(x.index)).reset_index()

rc('font', weight='bold')
fig = plt.figure(figsize=(15,3),dpi=200)
rcParams.update({'figure.autolayout': True})
ax = fig.add_subplot(111)
plt.xlim(-10, total_aa+1)
plt.ylim(0, 8)
plt.bar(tmp.aa_bin, tmp[0], capsize=2, color='black')
plt.tick_params(axis='y', labelsize=15, )
plt.tick_params(axis='x', labelsize=15)
plt.xticks(range(0, total_aa+1, 200))

fontdict = {'fontweight':'bold', 'size':18}
fontdict1 = {'fontweight':'bold', 'size':14}
fontdict2 = {'fontweight':'bold', 'size':12}
plt.ylabel("gnomAD Variant Count", fontdict=fontdict1)
plt.xlabel("CDS position", fontdict=fontdict)

# add patches
s_lengend = set()
cmap = matplotlib.cm.get_cmap('Spectral')
all_dm = []
[all_dm.append(_[1]) for _ in domain_aa_info if _[1] not in all_dm]
green = [48, 207, 0]
red = [255, 83, 83]
blue = [90, 92, 255]

d_color = {idm: [_/255 for _ in icolor] for idm, icolor in zip(all_dm, [green, red, blue])}


for idm in domain_aa_info:
    dm_id, dm_name, s, e = idm
    if dm_name not in s_lengend:
        s_lengend.add(dm_name)
        rect = matplotlib.patches.Rectangle((s, 0), e-s, 6, fill=True, alpha=0.25, color=d_color[dm_name], label=dm_name)
    else:
        rect = matplotlib.patches.Rectangle((s, 0), e-s, 6, fill=True, alpha=0.25, color=d_color[dm_name],)
    ax.add_patch(rect)
plt.legend(loc='upper right', prop={'size' :12})
plt.savefig(f'/Users/files/work/generisk/final/{gn}.gnomad.misssense.png',bbox_inches = "tight")






def most_least_constrained_domains():
    import statsmodels.stats.multitest
    os.chdir(f'{home}/gr/final')
    tmp = z.loc[z.dm_count > 20, 'dm'].unique()  # 202

    # olfactory receptor = PF13853,  count = 373 genes have this domain
    # utest pvalue, alternative = zscore less than average, pvalue = 1.4e-35

    d = pd.DataFrame({'dm': [], 'utest-less': [], 'utset-more': [], 'utest-2side': []})
    t = z.zscore  # universal zscore
    for i in tmp:
        t1 = z.loc[z.dm == i, 'zscore']
        p1 = st.mannwhitneyu(t, t1, alternative='less').pvalue
        p2 = st.mannwhitneyu(t, t1, alternative='greater').pvalue
        p3 = st.mannwhitneyu(t, t1, alternative='two-sided').pvalue

        d = d.append(pd.Series([i, p1, p2, p3], index=['dm', 'utest-less',
                                                       'utest-more', 'utest-2side']), ignore_index=True)

    d1 = d.merge(pid, on='dm')  # 200,  2 are lost
    tmp = z[['dm', 'dm_count', 'zscore_std', 'zscore_mean']].drop_duplicates()
    d1 = d1.merge(tmp, on='dm')   # 200
    d1 = d1.sort_values('utest-2side')
    d1['FDR_less'] = statsmodels.stats.multitest.multipletests(d1['utest-less'], method = 'fdr_bh')[1]
    d1['FDR_more'] = statsmodels.stats.multitest.multipletests(d1['utest-more'], method = 'fdr_bh')[1]

    d1.to_excel('constrained_domains_across_gene_pvalue.xlsx', index=False)
    d1.to_pickle('constrained_domains_across_gene_pvalue.pk')


def highest_lowest_ddi():
    # highest and lowest domains
    os.chdir('/Users/zeiss/gr/final')

    z_highest = z.loc[z.zrank > 99].sort_values('zrank', ascending=False)
    z_lowest = z.loc[z.zrank < 1].sort_values('zrank', ascending=True)

    z_highest.to_excel("~/gr/final/least_damaged_ddi_gt_0.99.xlsx", index=False, na_rep='NA')
    z_lowest.to_excel("~/gr/final/most_damaged_ddi_lt_0.01.xlsx", index=False, na_rep='NA')


def get_cosg():
    cosg = pd.read_csv('/Users/files/ref/cosmic/cosmic.in_cds.anno.refined.cadd.txt.gz', sep='\t',
                       header=None, names='chr,s,e,ref,alt,sig,_1,gn,_2,cadd'.split(','))
    cosg1 = cosg.loc[cosg._1 == 'missense_variant'].drop('_1,_2,cadd'.split(','), axis=1).drop_duplicates()  # (84598, 7
    cosdm = pd.read_csv('/Users/files/ref/cosmic/cosmic.in_domain.anno.cadd.bed.gz', sep='\t',
                        header=None, names='chr,s,e,ref,alt,sig,_1,cadd,gn,dm,dmname'.split(','))
    cosdm1 = cosdm.loc[cosdm._1 == 'missense_variant'].drop(
        '_1,cadd,dmname'.split(','), axis=1).drop_duplicates()  # (41535, 8)
    # pathogenic
    cosg_path = cosg1.loc[cosg1.sig == 'PATHOGENIC']  # (66297, 7)
    cosdm_path = cosdm1.loc[cosdm1.sig == 'PATHOGENIC']  # (35277, 8)
    # length
    len_dm = pd.read_pickle('/ref/basic/length_dm_cds.pk')  # 46326, include len_cds, len_dm
    len_cds = pd.read_csv("/ref/basic/length_gene_cds.tsv", sep='\t', names=['gn', 'len_cds'])  # 19622, only len_cds

    # gene / domain list
    gnlist = cosg_path[['gn']].drop_duplicates()  # 6343
    gnlist = gnlist.merge(len_cds, on='gn')
    dmlist = cosdm_path[['gn', 'dm']].drop_duplicates()  # 5867
    dmlist = dmlist.merge(len_dm, on=['gn', 'dm'])

    total_leng = gnlist.len_cds.sum()
    total_lendm = dmlist.len_dm.sum()
    print(f"""total gene num - genelist = {len(gnlist.index)}
    total dm num = {len(dmlist.index)}
    total gene len = {total_leng}
    total domain len = {total_lendm}
    len ratio ={total_lendm/total_leng} """)
    # total gene num - genelist = 6342
    # total dm num = 5831
    # total gene len = 16364160
    # total domain len = 4498491
    # len ratio =0.2748989865657632

    n_g = len(cosg_path.index)
    n_dm = len(cosdm_path.index)
    print(f"""pathogenic, missense, SNP count
    in gene = {n_g}
    in domain = {n_dm}
    ratio ={n_dm/n_g} """)

    # pathogenic, missense, SNP count
    # in gene = 66297
    # in domain = 35277
    # ratio =0.532105525136884

    #  count pathogenic varints in domain and gene
    ct_dm = cosdm_path.groupby(['gn', 'dm']).apply(lambda x: len(x.index)).reset_index()  # 5867
    ct_dm.columns = ['gn', 'dm', 'path_count_dm']
    ct_g = cosg_path.groupby(['gn']).apply(lambda x: len(x.index)).reset_index()  # 6343
    ct_g.columns = ['gn', 'path_count_gene']

    ct_dm2 = ct_dm.merge(ct_g, on='gn')  # 5837
    ct_dm3 = ct_dm2.merge(len_dm, on=['gn', 'dm'])  # 5801
    ct_dm3['ratio_len'] = round(ct_dm3['len_dm'] / ct_dm3['len_cds'], 6)
    ct_dm3['ratio_patho'] = round(ct_dm3['path_count_dm'] / ct_dm3['path_count_gene'], 6)
    ct_dm3['path_count_out'] = ct_dm3.path_count_gene - ct_dm3.path_count_dm
    ct_dm3['len_out'] = ct_dm3.len_cds - ct_dm3.len_dm

    ct_gmad = eval(open('/Users/files/work/generisk/snpcount/gnomad.genome.maf_0.001.count.adj_intron.pdict').read())

    gmad_dm = []
    gmad_cds = []
    for i in ct_dm3.itertuples(False):
        gn = i.gn
        dm = i.dm
        if gn not in ct_gmad:
            gmad_dm.append(np.nan)
            gmad_cds.append(np.nan)

        elif dm not in ct_gmad[gn]['dm']:
            gmad_dm.append(np.nan)
            gmad_cds.append(ct_gmad[gn]['cds_conseq']['missense_variant'][1])
        else:
            gmad_dm.append(ct_gmad[gn]['dm'][dm]['observed']['missense_variant'][1])
            gmad_cds.append(ct_gmad[gn]['cds_conseq']['missense_variant'][1])

    ct_dm3['gmad_cds'] = gmad_cds
    ct_dm3['gmad_dm'] = gmad_dm
    ct_dm3['ratio_gmad'] = ct_dm3['gmad_dm'] / ct_dm3['gmad_cds'] # 5801

    # calculate the domain odds ratio
    # Haldane-Anscombe correction
    # https://www.researchgate.net/post/How_to_calculate_OR_odd_ratio_if_one_of_groups_is_0_in_a_case-control_study
    ct_dm3['or'] = ((ct_dm3['path_count_dm'] + 0.5) / (ct_dm3['path_count_out'] + 0.5)) / \
        ((ct_dm3['gmad_dm'] + 0.5) / (ct_dm3['gmad_cds'] - ct_dm3['gmad_dm'] + 0.5))   # 5801

    # merge with zscore
    # 5745,   56 domains are lost
    ct_dm4 = ct_dm3.merge(z.drop(['len_dm', 'len_cds'], axis=1), on=['gn', 'dm'])
    ct_dm4.to_csv("/Users/zeiss/work/gr/final/cosmic_enrichment_for_heatscatter.csv", index=False)

    return ct_dm4


def get_clinvar_count():
    cling = pd.read_pickle(f'{home}/ref/clinvar/clinvar_cds.pk')  # (274964, 7)
    indm = pd.read_pickle(f'{home}/ref/clinvar/clinvar_domain.pk')  # (143972, 8)

    # pathogenic
    cling_path = cling.loc[cling.sig.isin(pathlist)]  # (34726, 7)
    indm_path = indm.loc[indm.sig.isin(pathlist)]  # (25793, 8)

    # length
    len_dm = pd.read_pickle('/ref/basic/length_dm_cds.pk')  # 46326, include len_cds, len_dm
    len_cds = pd.read_csv("/ref/basic/length_gene_cds.tsv", sep='\t', names=['gn', 'len_cds'])  # 19622, only len_cds

    # gene / domain list
    gnlist = cling_path[['gn']].drop_duplicates()  # 3533
    gnlist = gnlist.merge(len_cds, on='gn')  # 3533
    dmlist = indm_path[['gn', 'dm']].drop_duplicates()  # 4028
    dmlist = dmlist.merge(len_dm, on=['gn', 'dm'])  # 4028

    # calculate the domain length ratio
    total_leng = gnlist.len_cds.sum()
    total_lendm = dmlist.len_dm.sum()
    print(f"""total gene num - genelist = {len(gnlist.index)}
    total dm num = {len(dmlist.index)}
    total gene len = {total_leng}
    total domain len = {total_lendm}
    len ratio ={total_lendm/total_leng} """)

    # total gene num - genelist = 3533
    # total dm num = 4028
    # total gene len = 8440629
    # total domain len = 2995224
    # len ratio =0.35485791402512773

    n_g = len(cling_path.index)
    n_dm = len(indm_path.index)
    print(f"""pathogenic, missense, SNP count
    in gene = {n_g}
    in domain = {n_dm}
    ratio ={n_dm/n_g} """)

    # pathogenic, missense, SNP count
    # in gene = 34726
    # in domain = 25793
    # ratio =0.7427575879744284

    #  count pathogenic varints in domain and gene
    ct_dm = indm_path.groupby(['gn', 'dm']).apply(lambda x: len(x.index)).reset_index()  # 4028
    ct_dm.columns = ['gn', 'dm', 'path_count_dm']
    ct_g = cling_path.groupby(['gn']).apply(lambda x: len(x.index)).reset_index()  # 3533
    ct_g.columns = ['gn', 'path_count_gene']

    ct_dm2 = ct_dm.merge(ct_g, on='gn')  # 4027
    ct_dm3 = ct_dm2.merge(len_dm, on=['gn', 'dm'])  # 4027
    ct_dm3['ratio_len'] = round(ct_dm3['len_dm'] / ct_dm3['len_cds'], 6)
    ct_dm3['ratio_patho'] = round(ct_dm3['path_count_dm'] / ct_dm3['path_count_gene'], 6)
    ct_dm3['path_count_out'] = ct_dm3.path_count_gene - ct_dm3.path_count_dm
    ct_dm3['len_out'] = ct_dm3.len_cds - ct_dm3.len_dm

    ct_gmad = eval(open('/Users/files/work/generisk/snpcount/gnomad.genome.maf_0.001.count.adj_intron.pdict').read())

    gmad_dm = []
    gmad_cds = []
    for i in ct_dm3.itertuples(False):
        gn = i.gn
        dm = i.dm
        if gn not in ct_gmad:
            gmad_dm.append(np.nan)
            gmad_cds.append(np.nan)

        elif dm not in ct_gmad[gn]['dm']:
            gmad_dm.append(np.nan)
            gmad_cds.append(ct_gmad[gn]['cds_conseq']['missense_variant'][1])
        else:
            gmad_dm.append(ct_gmad[gn]['dm'][dm]['observed']['missense_variant'][1])
            gmad_cds.append(ct_gmad[gn]['cds_conseq']['missense_variant'][1])

    ct_dm3['gmad_cds'] = gmad_cds
    ct_dm3['gmad_dm'] = gmad_dm
    ct_dm3['ratio_gmad'] = ct_dm3['gmad_dm'] / ct_dm3['gmad_cds']

    # calculate the domain odds ratio
    # Haldane-Anscombe correction
    # https://www.researchgate.net/post/How_to_calculate_OR_odd_ratio_if_one_of_groups_is_0_in_a_case-control_study
    ct_dm3['or'] = ((ct_dm3['path_count_dm'] + 0.5) / (ct_dm3['path_count_out'] + 0.5)) / \
        ((ct_dm3['gmad_dm'] + 0.5) / (ct_dm3['gmad_cds'] - ct_dm3['gmad_dm'] + 0.5))   # 4027

    # merge with zscore
    # 3982,   45 domains are lost
    ct_dm4 = ct_dm3.merge(z.drop(['len_dm', 'len_cds'], axis=1), on=['gn', 'dm'])
    ct_dm4.to_csv("/Users/zeiss/work/gr/final/clinvar_enrichment_for_heatscatter.csv", index=False)

    return ct_dm4


# plot lollipop

def low_gdi_high_ddi():
    zm = z.merge(out_idx, on=['gn', 'conseq'])
    zm = zm.drop('gdi,gdi_phred_raw,rvis'.split(","), axis=1)

    # select low gdi, high ddi
    z_low_gdi = zm.loc[(zm.rvis_rk < 20) | (zm.exac_zrank < 20) | (zm.gdi_rank < 20)]
    z_low_gdi_zrank90 = z_low_gdi.loc[z_low_gdi.zrank > 90]  # (265, 15)

    z_low_gdi_zrank90_clinvar = z_low_gdi_zrank90.merge(ct_cds, on='gn')  # intersect with clinvar count, 79 genes left

    # merge with clinvar count in dm, then fillna
    z_low_gdi_zrank90_clinvar = z_low_gdi_zrank90_clinvar.merge(ct_dm, on=['gn', 'dm'], how='left')
    z_low_gdi_zrank90_clinvar['count_dm'] = z_low_gdi_zrank90_clinvar['count_dm'].fillna(0)

    z_low_gdi_zrank90_clinvar.sort_values(['count_cds', 'count_dm'], ascending=False, inplace=True)
    z_low_gdi_zrank90_clinvar = z_low_gdi_zrank90_clinvar.merge(dmname, on='dm')

    z_low_gdi_zrank90_clinvar.to_csv(
        "/Users/zeiss/generisk/final/low_gdi_high_ddi.clinvar_count.tsv", sep='\t', index=False)


def plot_lollipop(clinvar_path_count):
    """
    clinvar_path_count is a df, the return value of data_for_lollipop
    """
    os.chdir("/Users/Files/nogit/generisk/lollipop/figure")
    app = '/Users/files/work/bin/lollipops'
    for i in x.itertuples(False, None):
        # uid, variants, gn
        os.system("%s -domain-labels=fit  -U %s -o %s.svg %s" % (app, i[0], i[2], i[1]))
        #   os.system("rsvg-convert   %s.svg -o %s_rsvg.png"%(i[2],i[2]))
        os.system("convert -density 500 -units pixelsperinch %s.svg %s.png" % (i[2], i[2]))


def get_out_idx():
    """
    get the gene level score, RVIS, GDI, exacz
    merge method = outer
    """
    exacz = pd.read_csv(f'{home}/ref/exac/exac_zscore_mimssense+stopgain_gn_checked.txt', sep='\t')
    exacz = exacz[['gn', 'conseq', 'exac_z', 'exac_zrank']]

    gdi = pd.read_csv(f'{home}/work/generisk/gdi/gdi_score_pnas_gn_checked.txt', sep="\t")
    gdi = gdi[['gn', 'gdi', 'gdi_phred_raw']]
    gdi['gdi_rank'] = 100 - round(gdi['gdi'].rank() / len(gdi.index) * 100, 2)

    rvis = pd.read_csv(f"{home}/ref/rvis/rvis_lite.txt", sep='\t')

    out_idx = pd.merge(exacz, gdi, on='gn', how='outer')
    out_idx = pd.merge(out_idx, rvis, on='gn', how='outer')

    # merge with omim
    omim = pd.read_csv(f"{home}/ref/omim/omim_dedup.tsv", sep='\t', usecols='gn,inher'.split(","))
    out_idx = pd.merge(out_idx, omim, on='gn', how='left')
    out_idx['inher'] = out_idx['inher'].fillna('na')

    # 18090
    out_idx = out_idx.loc[out_idx.conseq == 'missense_variant'].drop('conseq', axis=1)

    out_idx.to_pickle(f'{home}/gr/final/out_idx.pk')

    return out_idx


def data_for_lollipop():
    os.chdir("/Users/files/ref/clinvar")

    c = pd.read_csv(
        "clinvar.in_cds.anno.exonic_variant_function.gz", sep='\t', header=None,
        names='_1,conseq,gn1,chr,s,e,ref,alt,sig,gn_bedtools'.split(",")).drop_duplicates(
        ['chr', 's', 'ref', 'alt', 'gn1'])
    c = c.replace({
        'conseq':
        {'nonsynonymous SNV': 'missense_variant', 'synonymous SNV': 'synonymous_variant',
         'frameshift substitution': 'frameshift', 'nonframeshift substitution': 'inframe',
         'stopgain': 'stop_gained', 'stoploss': 'stop_lost'}})
    c['gn,ts,aachange'.split(",")] = c.gn1.str.split("[:,]", expand=True)[[0, 1, 4]]
    c.aachange = c.aachange.str.replace("p.", '')

    geneinfo_basic = pd.read_csv("/Users/files/ref/basic/gene_full_info.tsv", sep='\t',
                                 header=None, names='gn,strand,chr,len,uid,ts'.split(','), usecols='gn,uid'.split(","))
    c1 = c.merge(geneinfo_basic, on='gn')
    c2 = c1['gn,sig,chr,s,e,ref,alt,conseq,aachange,ts,uid'.split(",")]

    sig2 = ['Pathogenic', 'Likely_pathogenic']
    sig3 = ['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic', 'drug_response']

    c3 = c2.loc[(c2.sig.isin(sig2)) & (c2.conseq == 'missense_variant')]

    # prepare the command for lollipops

    x = c3.groupby('uid')['aachange'].apply(lambda x: x.str.cat(sep=' ')).reset_index()
    clinvar_path_count = x.merge(c3[['gn', 'uid']].drop_duplicates(), on=['uid'])
    clinvar_path_count.to_csv(
        "/Users/Files/nogit/generisk/lollipop/clinvar.pathogenic.sig2.missense.all.tsv", sep='\t', index=False)

    return clinvar_path_count





def genomic_to_cds(pos_gm, gene_info):
    """
    gene_info = the value of gene_full_info[gn]
    """
    strand = gene_info[2]
    rg = gene_info[0]
    head_tail = rg[0][:2] + rg[-1][:2]
    s, e = min(head_tail), max(head_tail)
    irg = [_ for _ in rg if _[0] <= pos_gm <= _[1]]

    if len(irg) < 0:
        return ['err', f'{pos_gm} not in cds range head={s}, tail={e}']
    if len(irg) > 1:
        return ['err', f'multiple match found: {pos_gm}  \n{irg}']

    irg = irg[0]
    len1 = irg[1] - pos_gm
    if len1 < 0:
        return ['err', f'wrong range order: strand={strand}, gm_pos={pos_gm}, irg={irg}']

    if strand == '+':
        if irg[2] > irg[3]:
            return ['err', f'wrong range cds order: strand = plus, irg={irg}']
        return irg[2] + len1

    elif strand == '-':
        if irg[2] < irg[3]:
            return ['err', f'wrong range cds order: strand = minus, irg={irg}']

        return irg[3] + len1

    else:
        return ['err', f'wrong strand {strand}']


def genomic_to_aa(pos_gm, gene_info):
    res =  genomic_to_cds(pos_gm, gene_info)
    try:
        return res//3 + 1
    except:
        return res



def get_clinvar():

    cling = pd.read_csv("/Users/files/ref/clinvar/clinvar20200415/clinvar.in_cds.anno.refined.nocadd.txt.gz", sep='\t')
    cling = cling[cling['conseq'] == 'missense_variant']  # (275400, 9)
    cling = cling.drop(['conseq', 'gn_bedtools'], axis=1)
    cling.rename(columns={'af': 'sig'}, inplace=True)
    cling = cling.drop_duplicates(subset=['chr', 's', 'ref', 'alt'])  # (274964, 9)

    indm = pd.read_csv("/Users/files/ref/clinvar/clinvar20200415/clinvar.in_domain.anno.txt.gz", sep='\t')
    indmraw = indm[indm['conseq'] == 'missense_variant']  # (270150, 11)
    indm = indmraw.loc[indmraw.dm != 'outdm']   # (145175, 11)
    indm = indm.drop_duplicates(subset=['chr', 's', 'ref', 'alt'])  # (143972, 11)
    indm = indm.drop(['conseq', 'dmname', 'cadd'], axis=1)
    indm.rename(columns={'af': 'sig'}, inplace=True)

    cling.to_pickle('/ref/clinvar/clinvar_cds.pk')
    indm.to_pickle('/ref/clinvar/clinvar_domain.pk')


def get_zscore():

    perm_file = '/Users/zeiss/gr/perm_result/gnomad_gm_site_intron_nocadd_perm_result.tsv'
    z = pd.read_csv(perm_file, sep='\t')
    z = z.loc[z.conseq == 'missense_variant']
    z.drop('conseq,perm_mean,perm_std,observ'.split(','), axis=1, inplace=True)  # (28507, 3)
    z = z.merge(pid, on='dm', how='left')  # (28507, 5)
    z['zrank'] = round(z['zscore'].rank() / len(z.index) * 100, 2)
    z['dm_count'] = z.groupby('dm')['gn'].transform('count')
    z['zscore_median'] = z.groupby('dm')['zscore'].transform('median')
    z['zscore_std'] = z.groupby('dm')['zscore'].transform('std')
    z['zscore_mean'] = z.groupby('dm')['zscore'].transform('mean')
    z = z.merge(len_dm, on=['gn', 'dm'])  # (28384, 12)
    z.to_pickle('/ref/basic/gene_risk_zscore_gnomad_gm_site_intron_no_cadd.pk')


def gene_element_len():
    d_dm_len = {}
    dm_bad = {}  # 0, no bad found

    # get cds len
    hgnc = pd.read_csv("/Users/files/ref_basic/hgnc.coding+pseudo.txt", sep='\t')
    uniprot = eval(open('/Users/files/ref_basic/uniprot.pdict').read())
    all_to_uidfull = eval(open('/ref/basic/all_2_uidfullname.pdict').read())

    bad_uid = []
    cds_len = {}
    dup = []
    no_uid = []
    for gn, uid, mol_type in zip(hgnc.symbol, hgnc.uniprot, hgnc.locus_type):
        if uid is np.nan:
            if mol_type == 'gene with protein product':
                no_uid.append(gn)
            continue
        try:
            uid_full = all_to_uidfull['uid'][uid][0]
            aalen = uniprot[uid_full]['aa']
            len1 = aalen * 3 + 3
            if gn not in cds_len:
                cds_len[gn] = len1
            elif len1 == cds_len[gn]:
                continue
            elif len1 < cds_len[gn]:
                dup.append([gn, len1, cds_len[gn]])
            else:
                dup.append([gn, len1, cds_len[gn]])
                cds_len[gn] = len1
        except:
            bad_uid.append([gn, uid, all_to_uidfull['uid'].get(uid)])
    # cds_len = 19622
    with open('/ref/basic/length_gene_cds.tsv', 'w') as out:
        for k, v in cds_len.items():
            print(f'{k}\t{v}', file=out)

    # get_dm len
    for ign, v in ginfo.items():
        d_dm_len[ign] = {}

        for v2 in v[1].values():
            irg = v2[0]
            idm = v2[1]
            i_len = 0
            for _ in irg:
                i_len += abs(_[3] - _[2]) + 1
            if i_len < 1:
                dm_bad[ign][idm] = v2
            else:
                try:
                    d_dm_len[ign][idm] += i_len
                except:
                    d_dm_len[ign][idm] = i_len

    with open('/ref/basic/length_domain.tsv', 'w') as out:
        for ign, v in d_dm_len.items():
            for idm, v2 in v.items():
                print(f'{ign}\t{idm}\t{v2}', file=out)

    len_cds = pd.read_csv("/ref/basic/length_gene_cds.tsv", sep='\t', names=['gn', 'len_cds'])  # 19237
    len_dm = pd.read_csv("/ref/basic/length_domain.tsv", sep='\t', names=['gn', 'dm', 'len_dm'])  # 46333
    d2 = len_dm.merge(len_cds, on='gn')  # 46326
    d2.to_pickle('/ref/basic/length_dm_cds.pk')


def get_clinvar_gene_list():
    gnlist = cling_path[['gn']].drop_duplicates()  # 3533
    gnlist = gnlist.merge(len_dm[['gn', 'len_cds']].drop_duplicates(), on='gn')  # 3430
