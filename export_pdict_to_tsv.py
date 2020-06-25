import os,sys
home=os.path.expanduser('~')
sys.path.append('%s/jb/module'%(home))

import chctool
newlist = chctool.newlist

project_name = sys.argv[1]
project_name = project_name.rsplit("/",1)[-1]
pw = "/home/chenh19/gr/perm/%s"%(project_name)
pwsum = "/home/chenh19/gr/perm/%s/pdict"%(project_name)

fl_out = "%s/%s_perm_result.tsv"%(pw,project_name)

d_all ={}

x = os.popen("ls %s/*.pdict"%(pwsum)).read().split("\n")
x = newlist(x)

for i in x:
    try:
        dtmp = eval(open(i).read())
        d_all.update(dtmp)
    except:
        print    ('invalid pdict  ',i)

with open("%s/%s_perm_result.pdict"%(pw,project_name),'w') as out:
    print (d_all,file=out)


# dtmp = eval(open("%s/gene_list_fail_to_get_codon.result.pdict"%(pwsum)).read())
# print(len(dtmp))

# d_all.update(dtmp)
print ('project = ', project_name)
print (f'total genes in pdict = {len(d_all)}')
out = open(fl_out,'w')
print ('gn,dm,conseq,zscore,perm_mean,perm_std,observ'.replace(",",'\t'),file=out)


# {'ZMYM5': {'PF06467': {'missense_variant': [-0.047, 248.427, 81.302, 252.26], 'synonymous_variant': [1.559, 95.807, 48.818, 19.69], 'stop_gained': [0.194, 9.345, 22.139, 5.04]}}, 'ZDHHC20': {'PF01529': {'missense_variant': [0.841, 233.221, 79.66, 166.26], 'synonymous_variant': [-1.875, 45.731, 25.415, 93.38], 'stop_gained': [0.409, 10.699, 26.153, 0.0]}}, 'ZAR1L': {'PF13695': {'missense_variant': [0.729, 179.531, 89.952, 114.0], 'synonymous_variant': [0.682, 43.353, 34.178, 20.06], 'stop_gained': [0.388, 13.541, 34.929, 0.0]}}, 'ZMYM2': {'PF06467': {'missense_variant': [1.654, 482.479, 118.635, 286.23], 'synonymous_variant': [-0.803, 107.273, 41.613, 140.69], 'stop_gained': [0.55, 21.772, 33.057, 3.58]}, 'PF12012': {'missense_variant': [1.162, 258.884, 90.324, 153.89], 'synonymous_variant': [0.11, 57.164, 29.636, 53.89], 'stop_gained': [0.342, 12.663, 26.249, 3.68]}}, 'ZC3H13': {'PF00642': {'missense_variant': [0.612, 53.885, 44.185, 26.83], 'synonymous_variant': [0.099, 13.948, 16.409, 12.32], 'stop_gained': [0.239, 2.997, 12.543, 0.0]}}, 'WDFY2': {'PF01363': {'missense_variant': [1.242, 216.829, 84.562, 111.8], 'synonymous_variant': [1.034, 39.421, 25.238, 13.33], 'stop_gained': [0.464, 12.57, 27.066, 0.0]}, 'PF00400': {'missense_variant': [-0.797, 309.818, 94.149, 384.86], 'synonymous_variant': [1.099, 62.398, 30.082, 29.34], 'stop_gained': [0.277, 19.313, 32.982, 10.17]}}, 'ZIC5': {'PF18366': {'missense_variant': [1.142, 102.037, 68.821, 23.43], 'synonymous_variant': [0.427, 37.56, 37.098, 21.73], 'stop_gained': [0.101, 0.564, 5.599, 0.0]}, 'PF00096': {'missense_variant': [1.663, 185.696, 98.597, 21.69], 'synonymous_variant': [0.447, 46.902, 39.384, 29.28], 'stop_gained': [0.241, 3.848, 15.984, 0.0]}}, 'ZIC2': {'PF18366': {'missense_variant': [1.547, 136.706, 79.291, 14.01], 'synonymous_variant': [-0.139, 32.015, 28.727, 36.01], 'stop_gained': [0.29, 5.876, 20.294, 0.0]}, 'PF00096': {'missense_variant': [2.035, 221.753, 100.539, 17.15], 'synonymous_variant': [0.694, 57.347, 41.041, 28.88], 'stop_gained': [0.188, 2.653, 14.148, 0.0]}}, 'XPO4': {'PF08767': {'missense_variant': [-0.94, 221.133, 78.427, 294.88], 'synonymous_variant': [0.444, 45.498, 24.635, 34.56], 'stop_gained': [0.379, 8.098, 21.366, 0.0]}}}


for gn,v1 in d_all.items():
    for idm,v2 in v1.items():
        for cons,v3 in v2.items():
            res = [gn,idm,cons] + v3
            res = list(map(str,res))
            print ('\t'.join(res),file = out)
            # print('\t'.join([gn,idm,cons]+[str(round(_,2)) for _ in v3]),file=out)
out.close()

os.system("ln -s -f {0}/{1}_perm_result.tsv /home/chenh19/gr/perm/pdict".format(pw,project_name))
os.system("ln -s -f {0}/{1}_perm_result.pdict /home/chenh19/gr/perm/pdict".format(pw,project_name))
