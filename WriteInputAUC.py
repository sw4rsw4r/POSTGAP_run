import cPickle
import glob
import sys

#fLst = glob.glob('/data/hwangsw/gg4_hwangse/master/postgap/'+sys.argv[1]+"/res[0-9]*.pkl")
#fLst = glob.glob('/home/msung/hwangsw/Projects/postgap/simulation_POSTGAP/'+sys.argv[1]+"/res[0-9]*.pkl")
# fLst = glob.glob('/Users/msung/projects/p_postgap/postgap/simulation_POSTGAP/'+sys.argv[1]+"/res[0-9]*.pkl")
fLst = glob.glob('/Users/msung/projects/p_postgap/postgap/Simple_POSTGAP/'+sys.argv[1]+"/res[0-9]*.pkl")
#fLst = [f for f in fLst if not 'test' in f ]
fLst = [f for f in fLst if not 'TEST' in f ]
fLst = [f for f in fLst if not 'clusters' in f ]

for idx, f_path in enumerate(fLst):
    with open(f_path,'r') as f:
        res_out_lst = cPickle.load(f)
        for res_out in res_out_lst:
            for tissue, out in res_out.collocation_posterior.items():
                gwas_PIP, CLPP = out
                print CLPP
