





def finemap_v1(z_scores, beta_scores, cov_matrix, n, labels, sample_label, lambdas, mafs, annotations, kstart=1, kmax=1, corr_thresh=0.9, max_iter=300, output="configuration", prior="independence_robust", v_scale=0.0025, g="BRIC", eigen_thresh=0.1, verbose=False, isGWAS=False):

    # shotgun stochastic search
    if kstart < kmax:
        n_all_pos_config = sum(scipy.special.comb(
            len(z_scores), r, exact=True) for r in range(kstart, kmax+1))


 


