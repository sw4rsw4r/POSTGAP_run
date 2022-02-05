#import cPickle
import cPickle as pickle
import glob
import sys

fLst = glob.glob(sys.argv[1])
fLst = [f for f in fLst if not 'test' in f ]
fLst = [f for f in fLst if not 'clusters' in f ]

for idx, f_path in enumerate(fLst):
    with open(f_path,'r') as f:
        res_out_lst = pickle.load(f)
        for res_out in res_out_lst:
            for tissue, out in res_out.collocation_posterior.items():
                gwas_PIP, CLPP = out
                #print gwas_PIP.tolist()
#                PIP = gwas_PIP.tolist()
                for pip in gwas_PIP:
                    print(float(pip))
