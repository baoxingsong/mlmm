"""
Created: 25/12/2017  00:40:39
Author: Baoxing Song
Email: songbaoxing168@163.com
This source code is partially adopted from https://github.com/bvilhjal/mixmogam/blob/master/linear_models.py
And the R source code of emma and the intel version c code of emmax is referred
The GEMMA paper is referred for the EM and NR algorithm
And the SVS gives a lot of details could not be found any where else
http://doc.goldenhelix.com/SVS/latest/svsmanual/mixedModelMethods/overview.html#overviewofmixedlinearmodels
"""

import sys
this_path = "/Users/song/Dropbox/PycharmProjects/myEmmax"
sys.path = [this_path] + sys.path

import numpy as np
import scipy as sp
from numpy import linalg
from scipy import stats
from scipy import optimize
import warnings
import kinship
import re
import inspect, os
import pickle
from dataParsers import parse_plink_tped_file
from dataParsers import parse_plink_fam_phenotype_file

def qr_decomp(X):
    """
    QR decomposition.  Adaptations to changes between scipy versions.
    """
    ver_list = tuple(map(int, (sp.__version__).split('.')[:2]))
    if ver_list >= (0, 9):
        return linalg.qr(X, mode='reduced')  # Do the QR-decomposition for the Gram-Schmidt process.
    else:
        return linalg.qr(X, econ=True)  # Do the QR-decomposition for the Gram-Schmidt process.

class LinearMixedModel:
    """
    A class for linear mixed models
    """
    def __init__(self, Y=None, k=None, dtype='double'):
        """
        The fixed effects should be a list of fixed effect lists (SNPs)
        """
        self.n = len(Y)
        self.y_var = np.var(Y, ddof=1, dtype=dtype)
        self.Y = np.matrix(Y, dtype=dtype)
        self.Y.shape = (self.n, 1)
        self.X = np.matrix(np.ones((self.n, 1), dtype=dtype))  # The intercept
        self.p = 1
        self.beta_est = None
        self.K = k
        # A list of random effect type, and the cov matrix.
        self.random_effects = [('normal', np.matrix(np.identity(self.n)))]  # The first random effect is the IID error.

    def add_factor(self, x, lin_depend_thres=1e-4):
        """
        Adds an explanatory variable to the X matrix.
        """
        # Checking whether this new cofactor in linearly independent.
        new_x = np.array(x)
        new_x.shape = len(x)
        (beta, rss, rank, sigma) = linalg.lstsq(self.X, new_x)
        if float(rss) < lin_depend_thres:
            warnings.warn('A factor was found to be linearly dependent on the factors already in the X matrix. '
                          'Hence skipping it!')
            return False
        new_x.shape = (self.n, 1)
        self.X = np.hstack([self.X, new_x])
        self.p += 1
        return True

    def _get_eigen_L_(self, K=None, dtype='double'): # in the source code of emma, there is a matrix Z,
        # I think we should never use Z here
        if K is None:
            K = self.K
        evals, evecs = linalg.eigh(K)
        evals = np.array(evals, dtype=dtype)
        return {'values': evals, 'vectors': np.mat(evecs, dtype=dtype).T}

    def _get_eigen_R_(self, X=None, K=None, hat_matrix=None, dtype='double'):
        if X is None:
            X = self.X
        q = X.shape[1]
        if not hat_matrix:
            X_squared_inverse = linalg.pinv(X.T * X)  # (X.T*X)^{-1}
            # linalg.pinv: Calculate the generalized inverse of a matrix using
            # its singular-value decomposition (SVD) and including all large singular values.
            hat_matrix = X * X_squared_inverse * X.T
        if K is None:
            K = self.K
        S = np.mat(np.identity(self.n)) - hat_matrix  # S=I-X(X'X)^{-1}X'
        M = np.mat(S * (K + np.matrix(np.identity(self.n))) * S, dtype='double')
        evals, evecs = linalg.eigh(M)  # eigen of S(K+I)S
        eig_values = np.array(evals[q:], dtype=dtype) - 1  # Don't know what minute one, the emma R code did this
        return {'values': eig_values, 'vectors': (np.mat(evecs, dtype=dtype).T[q:])} # here it seems different with R,
        # it is due to the eigen values result from R is decreasing while from numpy it is increasing
        # if the X is None, this result should be same with _get_eigen_L

    def _ll_(self, delta, eig_vals, eig_vals_L, sq_etas):
        n = self.n
        c_1 = 0.5 * n * (sp.log(n / (2.0 * sp.pi)) - 1)
        v1 = eig_vals + delta
        v2 = eig_vals_L + delta
        res = c_1 - 0.5 * (n * sp.log(sp.sum(sq_etas.flatten() / v1)) + sp.sum(sp.log(v2)))
        return res # log-likelihoods (eq. 6 from paper), this is the emma.delta.ML.LL.wo.Z function in R emma

    def _dll_(self, delta, eig_vals, eig_vals_L, sq_etas):
        v1 = eig_vals + delta
        v2 = sq_etas.flatten() / v1
        v3 = eig_vals_L + delta
        res = (self.n * sp.sum(v2 / v1) / sp.sum(v2) - sp.sum(1.0 / v3))
        return res # diffrentiated log-likelihoods (eq. 8 from paper) the emma.delta.ML.dLL.wo.Z function in R emma

    def _rell_(self, delta, eig_vals, sq_etas):
        num_eig_vals = len(eig_vals)
        c_1 = 0.5 * num_eig_vals * (sp.log(num_eig_vals / (2.0 * sp.pi)) - 1)
        v = eig_vals + delta
        res = c_1 - 0.5 * (num_eig_vals * sp.log(sp.sum(sq_etas.flatten() / v)) + sp.sum(sp.log(v)))
        return res  # log-likelihoods (eq. 7 from paper)

    def _redll_(self, delta, eig_vals, sq_etas):
        num_eig_vals = len(eig_vals)
        v1 = eig_vals + delta
        v2 = sq_etas.flatten() / v1
        res = (num_eig_vals * sp.sum(v2 / v1) / sp.sum(v2) - sp.sum(1.0 / v1))
        return res  # diffrentiated log-likelihoods (*2) (eq. 9 from paper)

    def get_estimates_another_matrix(self, k2, xs=None, ngrids=[5, 5, 5, 5, 5], llim=-10, ulim=10, method='REML'):
        """
        Handles two K matrices, and one I matrix.
        Methods available are 'REML', and 'ML'
        """
        print("this function is not well tested")
        return 0
        if xs != None:
            X = sp.hstack([self.X, xs])
        else:
            X = self.X

        for it_i in range(len(ngrids)):
            delta = float(ulim - llim) / ngrids[it_i]
            # narrow this range, to get a roughly optimized value. So the final value is not the exact best value
            log_k_ratio = llim
            lls = []
            res_list = []
            for i in range(ngrids[it_i] + 1):
                k_ratio = sp.exp(log_k_ratio)
                a = k_ratio / (k_ratio + 1.0)
                K = a * self.K + (1 - a) * k2
                eig_L = self._get_eigen_L_(K=K)
                # Now perform EMMA
                res_dict = self.get_estimates(eig_L, K=K, xs=xs, ngrids=10, method=method, llim=-10, ulim=10)
                res_list.append(res_dict)
                lls.append(res_dict['max_ll'])
                log_k_ratio += delta
            max_ll_i = sp.argmax(lls)
            ulim = llim + delta * (max_ll_i + 1) # This range is updated
            llim = llim + delta * (max_ll_i - 1)

        opt_k_ratio = sp.exp(log_k_ratio)
        a = opt_k_ratio / (opt_k_ratio + 1)
        opt_k = kinship.scale_k(a * self.K + (1 - a) * k2)
        res_dict = self.get_estimates(eig_L, K=opt_k, xs=xs)
        res_dict['opt_k'] = opt_k
        res_dict['opt_k_ratio'] = opt_k_ratio
        res_dict['perc_var1'] = a * res_dict['pseudo_heritability']
        res_dict['perc_var2'] = (1 - a) * res_dict['pseudo_heritability']
        # this function is computationally density, should not use emma method to get the p_val
        # try to get the emmax method
        return res_dict

    def get_estimates(self, eig_L, K=None, xs=None, ngrids=50, llim=-10, ulim=10, esp=1e-6,
                      return_pvalue=False, return_f_stat=False, method='REML',
                      dtype='double', eig_R=None):
        """
        Get ML/REML estimates for the effect sizes, as well as the random effect contributions.
        Using the EMMA algorithm (Kang et al., Genetics, 2008).
        Methods available are 'REML', and 'ML'
        """

        if xs is None:
            X = self.X
        else:
            X = np.hstack([self.X, xs])

        if not (eig_R and xs != None):
            eig_R = self._get_eigen_R_(X=X, K=K)

        q = X.shape[1]  # number of fixed effects
        n = self.n  # number of individuls
        p = n - q
        m = ngrids + 1
        print ("eigen_L_V ", eig_L['values'])
        print ("eigen_R_V ", eig_R['values'])

        print ("eigen_L ", eig_L['vectors'])
        print ("eigen_R ", eig_R['vectors'])
        print ("self.Y ", self.Y)
        etas = np.array(eig_R['vectors'] * self.Y, dtype=dtype)
        print ("etas ", etas)
        sq_etas = etas * etas
        print ("sq_etas ", sq_etas)
        log_deltas = (np.arange(m, dtype=dtype) / ngrids) * (ulim - llim) + llim  # a list of deltas to search
        deltas = np.exp(log_deltas)
        eig_vals = np.array(eig_R['values'], dtype=dtype)
        lambdas = np.reshape(np.repeat(eig_vals, m), (p, m)) + np.reshape(np.repeat(deltas, p), (m, p)).T
        s1 = np.sum(sq_etas / lambdas, axis=0)
        print ("line204 ", s1)
        if method == 'REML':
            s2 = np.sum(np.log(lambdas), axis=0)
            print ("line284 ", s2)
            lls = 0.5 * (p * (sp.log((p) / (2.0 * sp.pi)) - 1 - sp.log(s1)) - s2)# log likelihoods (eq. 7 from paper)
            print ("line288 ", lls)
            s3 = np.sum(sq_etas / (lambdas * lambdas), axis=0)
            print ("s3", s3)
            s4 = np.sum(1 / lambdas, axis=0)
            dlls = 0.5 * (p * s3 / s1 - s4) # this is the derivation of log likelihood (eq. 9 from paper)
            print("s4", s4)
            print ("dlls", dlls)
        elif method == 'ML': # this part is the function emma.MLE in R emma
            eig_vals_L = sp.array(eig_L['values'], dtype=dtype)
            xis = sp.reshape(sp.repeat(eig_vals_L, m), (n, m)) + \
                  sp.reshape(sp.repeat(deltas, n), (m, n)).T
            s2 = np.sum(np.log(xis), axis=0)
            lls = 0.5 * (n * (np.log((n) / (2.0 * np.pi)) - 1 - np.log(s1)) - s2)# log likelihoods (eq. 6 from paper)
            s3 = sp.sum(sq_etas / (lambdas * lambdas), axis=0)
            s4 = sp.sum(1 / xis, axis=0)
            dlls = 0.5 * (n * s3 / s1 - s4)            # this is the derivation of log likelihood (eq. 8 from paper)

        max_ll_i = sp.argmax(lls)
        max_ll = lls[max_ll_i]

        last_dll = dlls[0]
        last_ll = lls[0]

        print("max_ll_i ",  max_ll_i)
        print("max_ll " , max_ll )
        print("last_dll " , last_dll )
        print("last_ll " , last_ll )

        zero_intervals = []
        for i in range(1, len(dlls)):
            if dlls[i] < 0 and last_dll > 0:
                zero_intervals.append(((lls[i] + last_ll) * 0.5, i))# There is likelihoods value peak in
                # thie interval, go to this interval search the maximum likelihood
            last_ll = lls[i]
            last_dll = dlls[i]

        if len(zero_intervals) > 0:
            opt_ll, opt_i = max(zero_intervals)# what does this max mean?
            opt_delta = 0.5 * (deltas[opt_i - 1] + deltas[opt_i])
            # Newton-Raphson
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    if method == 'REML':
                        new_opt_delta = optimize.newton(self._redll_, opt_delta, args=(eig_vals, sq_etas), tol=esp,
                                                        maxiter=100)
                    elif method == 'ML':
                        new_opt_delta = optimize.newton(self._dll_, opt_delta, args=(eig_vals, eig_vals_L, sq_etas),
                                                        tol=esp, maxiter=100)
            except Exception:
                new_opt_delta = opt_delta
            # Validating the delta
            if opt_i > 1 and deltas[opt_i - 1] - esp < new_opt_delta < deltas[opt_i] + esp:
                opt_delta = new_opt_delta
                opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
            # Cheking lower boundary
            elif opt_i == 1 and 0.0 < new_opt_delta < deltas[opt_i] + esp:
                opt_delta = new_opt_delta
                opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
            # Cheking upper boundary
            elif opt_i == len(deltas) - 1 and new_opt_delta > deltas[opt_i - 1] - esp \
                    and not np.isinf(new_opt_delta):
                opt_delta = new_opt_delta
                opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)

            if method == 'REML':
                opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
            elif method == 'ML':
                opt_ll = self._ll_(opt_delta, eig_vals, eig_vals_L, sq_etas)

            if opt_ll < max_ll:
                opt_delta = deltas[max_ll_i]
        else:
            opt_delta = deltas[max_ll_i]
            opt_ll = max_ll

        #        likelyhood = self._ll_(opt_delta, eig_vals, eig_L['values'], sq_etas) # this is used for LRT
        print("line 418 opt_delta ", opt_delta)
        print("line 419 opt_ll ", opt_ll)
        print("sq_etas ", sq_etas)
        print("eig_vals ", eig_vals)
        print("opt_delta ", opt_delta)
        R = sq_etas / (eig_vals + opt_delta)
        print("R ", R)
        opt_vg = np.sum(R) / p  # vg, p = n-q q is the number of fixed effects.
        print("opt_vg ", opt_vg)
        # This is the REML estimation. the ML estimation is np.sum(R) / n
        opt_ve = opt_vg * opt_delta  # ves
        print("opt_ve ", opt_ve)
        # the R emma.MLE and emma.REMLE code stopped here

        # the solve of mixed model equation is mentioned in
        # http://doc.goldenhelix.com/SVS/latest/svsmanual/mixedModelMethods/overview.html#overviewofmixedlinearmodels
        H_sqrt_inv = np.mat(np.diag(1.0 / np.sqrt(eig_L['values'] + opt_delta)), dtype=dtype) * eig_L['vectors'] # this
        print("H_shape ", np.shape(H_sqrt_inv))
        print("H_sqrt_inv ", H_sqrt_inv)
        # is the U value from R emma code
        # V = opt_vg * K + opt_ve * sp.eye(len(K))
        # H_sqrt = cholesky(V).T
        # H_sqrt_inv = H_sqrt.I
        X_t = H_sqrt_inv * X
        Y_t = H_sqrt_inv * self.Y
        (beta_est, mahalanobis_rss, rank, sigma) = linalg.lstsq(X_t, Y_t, rcond=-1.0)
        x_beta = X * beta_est
        print("x_beta ", x_beta)
        print("Y_t ", Y_t)
        print("Y_t - x_beta ", Y_t - x_beta)
        print("Y_t - x_beta test ", np.sum((Y_t - X_t*beta_est).T*(Y_t - X_t*beta_est)))
        residuals = self.Y - x_beta
        print("residuals ", residuals)
        rss = residuals.T * residuals
        print("rss ", rss)
        print("mahalanobis_rss ", mahalanobis_rss)
        # x_beta_var = sp.var(x_beta, ddof=1)
        # var_perc = x_beta_var / self.y_var

        # get the likelyhood value for LRT
        opt_ll = self._ll_(opt_delta, eig_vals, eig_L['values'], sq_etas) # recalculate likelyhood
        print("opt_ll", opt_ll)
        opt_dll = self._dll_(opt_delta, eig_vals, eig_L['values'], sq_etas) # recalculate likelyhood
        print("opt_dll", opt_dll)
        opt_rell = self._rell_(opt_delta, eig_vals, sq_etas)
        print("opt_rell", opt_rell)
        opt_redll = self._redll_(opt_delta, eig_vals, sq_etas)
        print("opt_redll", opt_redll)


        res_dict = {'max_ll': opt_ll, 'delta': opt_delta, 'beta': beta_est, 've': opt_ve, 'vg': opt_vg,
                    'rss': rss, 'mahalanobis_rss': mahalanobis_rss, 'H_sqrt_inv': H_sqrt_inv,
                    'pseudo_heritability': 1.0 / (1 + opt_delta)}

        if (xs is not None) and return_f_stat:
            #            rss_ratio = h0_rss / rss_list
            #            var_perc = 1 - 1 / rss_ratio
            #            f_stats = (rss_ratio - 1) * n_p / float(q)

            h0_X = H_sqrt_inv * self.X
            (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y_t, rcond=-1)
            print("h0_rss ", h0_rss)
            f_stat = (h0_rss / mahalanobis_rss - 1) * p / xs.shape[1]
            res_dict['var_perc'] = 1.0 - mahalanobis_rss / h0_rss
            res_dict['f_stat'] = float(f_stat)
            print("f_stat ", f_stat)
        if return_pvalue:
            p_val = stats.f.sf(f_stat, (xs.shape[1]), p)
            res_dict['p_val'] = float(p_val)
            print("p_val ", p_val)
        return res_dict  # , lls, dlls, sp.log(deltas)

    def emma_test(self, genotypes, f_test=True, lrt_test=True, ngrids=50, llim=-10, ulim=10, esp=1e-6, eig_L=None,
                  maf=0.05, missing_rate=0.2, outPutFile="emmax_f_test.ps"):
        """
        Single SNP analysis, i.e. EMMA, but faster than R EMMA.
        """
        f = open(outPutFile, 'w')
        K = self.K
        if eig_L is None:
            eig_L = self._get_eigen_L_(K)
        num_snps = len(genotypes.snps[0])
        # f_stats = np.full(num_snps, None)
        # vgs = np.full(num_snps, None)
        # ves = np.full(num_snps, None)
        # max_lls = np.full(num_snps, None)
        # var_perc = np.full(num_snps, None)
        # chr_list = np.empty(num_snps, dtype='U20')
        # position_list = np.empty(num_snps, dtype='u4')
        # ids_list = np.empty(num_snps, dtype='U20')
        # rss_list = np.full(num_snps, None)
        # betas = []
        # f_test_p_vals = np.full(num_snps, None, dtype='double')
        # lrt_test_p_vals = np.full(num_snps, None, dtype='double')

        # Run null model....
        if lrt_test:
            null_log_likeHoohd = self.get_estimates(eig_L=eig_L, xs=None, ngrids=ngrids, llim=llim, ulim=ulim,
                                                    esp=esp, return_pvalue=False, return_f_stat=False)['max_ll']
        # if the count of minor allele smaller than this number, do not care about it
        ma_count = maf * self.n
        la_count = self.n - ma_count
        missing_count = missing_rate * self.n

        the_snp_list = genotypes.snps.T.tolist() # this data structure should be faster than numpy array
        for i in np.arange(0, 1):
            if not re.findall("b$", genotypes.variants['id'][i][
                0]):  # take them as category variables, this is the specific format from Irisas
                snp = the_snp_list[i]
                a = np.array(np.unique(snp, return_counts=True)).T
                b = a[np.argsort(a[:, 1], kind='mergesort')]
                number = len(a)
                if (number > 1) and (a[0][0] != -1 or a[0][1] < missing_count):

                    this_variable = np.zeros(0)
                    for j in np.arange(number):
                        # if the frequency is too small, do not test it. If the frequency is too large, do not test it
                        if (b[j][1] >= ma_count) and (b[j][1] <= la_count and b[j][0] != -1): # filtering missing values
                            g = np.where(snp == b[j][0], 1, 0)
                            c = np.where(snp == b[j][0])
                            m = np.where(snp == -1)
                            mean = len(c[0]) / (len(g) - len(m[0]))
                            this_variable = np.append(this_variable, np.array(np.where((snp != -1), g, mean)))
                    this_variable = this_variable.reshape(int(len(this_variable) / self.n), int(self.n))
                    if this_variable.shape[0]==1 or (this_variable.shape[0]==2 and number==2):
                        if this_variable.shape[0] == 2:
                            this_variable = this_variable[0]
                            this_variable = this_variable.reshape(1, int(self.n))

                        if f_test:
                            res = self.get_estimates(eig_L=eig_L, xs=np.matrix(this_variable).T, ngrids=ngrids,
                                                     llim=llim,
                                                     ulim=ulim, esp=esp, return_pvalue=True, return_f_stat=True)
                        elif lrt_test:
                            res = self.get_estimates(eig_L=eig_L, xs=np.matrix(this_variable).T, ngrids=ngrids,
                                                     llim=llim, ulim=ulim,
                                                     esp=esp, return_pvalue=False, return_f_stat=True)
                        if f_test or lrt_test:
                            # f_stats[i] = res['f_stat']
                            # vgs[i] = res['vg']
                            # ves[i] = res['ve']
                            # max_lls[i] = res['max_ll']
                            # var_perc[i] = res['var_perc']
                            # betas.append(map(float, list(res['beta'])))
                            # rss_list[i] = res['rss']
                            # log_likeHood = res['max_ll']
                            # f_stats = res['f_stat']
                            # vgs = res['vg']
                            # ves = res['ve']
                            # max_lls = res['max_ll']
                            # var_perc = res['var_perc']
                            # betas.append(map(float, list(res['beta'])))
                            # rss_list[i] = res['rss']
                            log_likeHood = res['max_ll']
                        chr = genotypes.variants[i]['chromosome'][0]
                        pos = genotypes.variants[i]['position'][0]
                        id = genotypes.variants[i]['id'][0]
                        f.write("%s\t%d\t%s" % (chr, pos, id))
                        if lrt_test:
                            stat = 2 * (log_likeHood-null_log_likeHoohd)
                            # lrt_test_p_vals[i] = sp.stats.distributions.chi2.sf(stat, this_variable.shape[0])
                            lrt_test_p_vals = sp.stats.distributions.chi2.sf(stat, this_variable.shape[0])
                            f.write("\t%10.20f" % (lrt_test_p_vals))
                        if f_test:
                            # f_test_p_vals[i] = res['p_val']
                            f_test_p_vals = res['p_val']
                            f.write("\t%10.20f" % (f_test_p_vals))
                        f.write("\n")
                    elif this_variable.shape[0] >= 2:  # test one by one
                        for this_variable_index in range(this_variable.shape[0]):
                            this_variable_2 = this_variable[this_variable_index,:]
                            this_variable_2 = this_variable_2.reshape(1, int(self.n))
                            if f_test:
                                res = self.get_estimates(eig_L=eig_L, xs=np.matrix(this_variable_2).T, ngrids=ngrids,
                                                         llim=llim,
                                                         ulim=ulim, esp=esp, return_pvalue=True, return_f_stat=True)
                            elif lrt_test:
                                res = self.get_estimates(eig_L=eig_L, xs=np.matrix(this_variable_2).T, ngrids=ngrids,
                                                         llim=llim, ulim=ulim,
                                                         esp=esp, return_pvalue=False, return_f_stat=True)
                            if f_test or lrt_test:
                                log_likeHood = res['max_ll']
                            chr = genotypes.variants[i]['chromosome'][0]
                            pos = genotypes.variants[i]['position'][0]
                            id = genotypes.variants[i]['id'][0]
                            f.write("%s\t%d\t%s" % (chr, pos, id))
                            if lrt_test:
                                stat = 2 * (log_likeHood - null_log_likeHoohd)
                                # lrt_test_p_vals[i] = sp.stats.distributions.chi2.sf(stat, this_variable.shape[0])
                                lrt_test_p_vals = sp.stats.distributions.chi2.sf(stat, this_variable.shape[0])
                                f.write("\t%10.20f" % (lrt_test_p_vals))
                            if f_test:
                                # f_test_p_vals[i] = res['p_val']
                                f_test_p_vals = res['p_val']
                                f.write("\t%10.20f" % (f_test_p_vals))
                            f.write("\n")
            else:
                snp = the_snp_list[i]
                a = np.array(np.unique(snp, return_counts=True)).T
                number = len(a)
                if number > 1:
                    this_variable = np.array(snp)
                    this_variable = this_variable.reshape(1, int(self.n))
                    if f_test:
                        res = self.get_estimates(eig_L=eig_L, xs=np.matrix(this_variable_2).T, ngrids=ngrids,
                                                 llim=llim,
                                                 ulim=ulim, esp=esp, return_pvalue=True, return_f_stat=True)
                    elif lrt_test:
                        res = self.get_estimates(eig_L=eig_L, xs=np.matrix(this_variable_2).T, ngrids=ngrids,
                                                 llim=llim, ulim=ulim,
                                                 esp=esp, return_pvalue=False, return_f_stat=True)
                    if f_test or lrt_test:
                        log_likeHood = res['max_ll']
                    chr = genotypes.variants[i]['chromosome'][0]
                    pos = genotypes.variants[i]['position'][0]
                    id = genotypes.variants[i]['id'][0]
                    f.write("%s\t%d\t%s" % (chr, pos, id))
                    if lrt_test:
                        stat = 2 * (log_likeHood - null_log_likeHoohd)
                        # lrt_test_p_vals[i] = sp.stats.distributions.chi2.sf(stat, this_variable.shape[0])
                        lrt_test_p_vals = sp.stats.distributions.chi2.sf(stat, this_variable.shape[0])
                        f.write("\t%10.20f" % (lrt_test_p_vals))
                    if f_test:
                        # f_test_p_vals[i] = res['p_val']
                        f_test_p_vals = res['p_val']
                        f.write("\t%10.20f" % (f_test_p_vals))
                    f.write("\n")
            #     else:
            #         betas.append(None)
            # else:
            #     betas.append( None )

            # chr_list[i] = genotypes.variants[i]['chromosome'][0]
            # position_list[i] = genotypes.variants[i]['position'][0]
            # ids_list[i] = genotypes.variants[i]['id'][0]
        f.close()
        # this_result = {'f_stats': f_stats, 'vgs': vgs, 'ves': ves, 'var_perc': var_perc,
        #         'max_lls': max_lls, 'betas': betas, 'rss': rss_list, 'id': ids_list,
        #                'chr': chr_list, 'position': position_list}
        # if f_test:
        #     this_result['f_p_value'] = f_test_p_vals
        # if lrt_test:
        #     this_result['lrt_p_value'] = lrt_test_p_vals
        return 0

    # the emmax part
    def emmax_f_test(self, genotypes, method='REML',
                     eig_L=None, eig_R=None, maf=0.05, outPutFile="emmax_f_test.ps"):
        """
        EMMAX implementation (in python)
        Single SNPs
        """
        if not eig_L:
            eig_L = self._get_eigen_L_()
        if not eig_R:
            eig_R = self._get_eigen_R_(X=self.X)

        res = self.get_estimates(eig_L, method=method, eig_R=eig_R)  # Get the variance estimates...
        # similar with the 463 line of emmax code
        r = self._emmax_f_test_(genotypes, res['H_sqrt_inv'],  maf=maf, outPutFile=outPutFile)
        # r['pseudo_heritability'] = res['pseudo_heritability']
        # r['ve'] = res['ve']
        # r['vg'] = res['vg']
        # r['max_ll'] = res['max_ll']
        return r

    def _emmax_f_test_(self, genotypes, H_sqrt_inv, maf=0.05, missing_rate=0.2, outPutFile="emmax_f_test.ps"):
        """
        EMMAX implementation (in python)
        Single SNPs

        Methods:
            normal - Normal regression
            qr - Uses QR decomposition to speed up regression with many co-factors.
        """
        f = open(outPutFile, 'w')

        dtype = 'double'

        # if the count of minor allele smaller than this number, skip it
        ma_count = maf * self.n
        la_count = self.n - ma_count

        missing_count = missing_rate * self.n

        n = self.n
        h0_X = np.mat(H_sqrt_inv * self.X, dtype=dtype)
        Y = H_sqrt_inv * self.Y  # The transformed outputs.
        (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y, rcond=None)
        Y = sp.mat(Y - h0_X * h0_betas, dtype=dtype) # remove the effect of intercept and covariates to
        # accelerate computing
        # this is the MB^(-1)y in the SVS document
        Q, R = qr_decomp(h0_X)  # Do the QR-decomposition for the Gram-Schmidt process.
        Q = sp.mat(Q)
        Q2 = Q * Q.T
        M = sp.mat(H_sqrt_inv.T * (sp.eye(n) - Q2), dtype=dtype)

        the_snp_list = genotypes.snps.T.tolist()  # this data structure should be faster than numpy array
        num_snps = len(genotypes.snps[0])

        # chr_list = np.empty(num_snps, dtype='U20')
        # position_list = np.empty(num_snps, dtype='u4')
        # ids_list = np.empty(num_snps, dtype='U20')
        # f_stats = np.full(num_snps, None)
        # rss_list = np.full(num_snps, None)
        # p_vals = np.full(num_snps, None, dtype='double')

        for i in np.arange(0, num_snps):
            if not re.findall("b$", genotypes.variants['id'][i][0]): # take them as category variables, this is the
                # specific format from Irisas
                snp = the_snp_list[i]
                a = np.array(np.unique(snp, return_counts=True)).T
                b = a[np.argsort(a[:, 1], kind='mergesort')]
                number = len(a)
                if (number > 1) and (a[0][0] != -1 or a[0][1]<missing_count):
                    this_variable = np.zeros(0)
                    for j in np.arange(number):
                        # if the frequency is too small, do not test it. If the frequency is too large, do not test it
                        if (b[j][1] >= ma_count) and (b[j][1] <= la_count and b[j][0] != -1): # filter missing value
                            g = np.where(snp == b[j][0], 1, 0)
                            c = np.where(snp == b[j][0])
                            m = np.where(snp == -1)
                            mean = len(c[0]) / (len(g) - len(m[0]))
                            this_variable = np.append( this_variable, np.array(np.where((snp != 0), g, mean)) )
                    this_variable = this_variable.reshape(int(len(this_variable) / self.n), int(self.n))
                    if this_variable.shape[0]==1 or (this_variable.shape[0]==2 and number==2):
                        if this_variable.shape[0] == 2:
                            this_variable = this_variable[0,:]
                            this_variable = this_variable.reshape(1, int(self.n))
                        q = 1 # line 359/367 in the emmax intel version c code
                        p = len(self.X.T) + q
                        n_p = n - p
                        xs = np.matrix(this_variable)
                        xs = sp.mat(xs * M, dtype=dtype)
                        (betas, rss, rank, sigma) = linalg.lstsq(xs.T, Y, rcond=-1)
                        if rss:
                            #
                            # rss_list[i] = rss[0]
                            # rss_ratio = h0_rss / rss_list[i]
                            # f_stats[i] =(rss_ratio - 1) * n_p / float(q)
                            # p_vals[i] = stats.f.sf(f_stats[i], q, n_p)
                            rss_list = rss[0]
                            rss_ratio = h0_rss / rss_list
                            f_stats = (rss_ratio - 1) * n_p / float(q)
                            pvalue = stats.f.sf(f_stats, q, n_p)
                            chr = genotypes.variants[i]['chromosome'][0]
                            pos = genotypes.variants[i]['position'][0]
                            id = genotypes.variants[i]['id'][0]
                            f.write("%s\t%d\t%s\t%10.20f\n" % (chr, pos, id, pvalue))
                    elif this_variable.shape[0] >= 2: # test one by one
                        for this_variable_index in range(this_variable.shape[0]):
                            this_variable_2 = this_variable[this_variable_index,:]
                            this_variable_2 = this_variable_2.reshape(1, int(self.n))
                            q = 1 # line 359/367 in the emmax intel version c code
                            p = len(self.X.T) + q
                            n_p = n - p
                            xs = np.matrix(this_variable_2)
                            xs = sp.mat(xs * M, dtype=dtype)
                            (betas, rss, rank, sigma) = linalg.lstsq(xs.T, Y, rcond=-1)
                            if rss:
                                rss_list = rss[0]
                                rss_ratio = h0_rss / rss_list
                                f_stats = (rss_ratio - 1) * n_p / float(q)
                                pvalue = stats.f.sf(f_stats, q, n_p)
                                chr = genotypes.variants[i]['chromosome'][0]
                                pos = genotypes.variants[i]['position'][0]
                                id = genotypes.variants[i]['id'][0]
                                f.write("%s\t%d\t%s\t%10.20f\n" % (chr, pos, id, pvalue))
            else: # the boundary regions are token as quantitive variables
                snp = the_snp_list[i]
                a = np.array(np.unique(snp, return_counts=True)).T
                number = len(a)
                if number > 1:
                    this_variable = np.array(snp)
                    this_variable = this_variable.reshape(1, int(self.n))
                    q = 1  # line 359/367 in the emmax intel version c code
                    p = len(self.X.T) + q
                    n_p = n - p
                    xs = np.matrix(this_variable)
                    xs = sp.mat(xs * M, dtype=dtype)
                    (betas, rss, rank, sigma) = linalg.lstsq(xs.T, Y, rcond=-1)
                    if rss:
                        rss_list = rss[0]
                        rss_ratio = h0_rss / rss_list
                        f_stats = (rss_ratio - 1) * n_p / float(q)
                        pvalue = stats.f.sf(f_stats, q, n_p)
                        chr = genotypes.variants[i]['chromosome'][0]
                        pos = genotypes.variants[i]['position'][0]
                        id = genotypes.variants[i]['id'][0]
                        f.write("%s\t%d\t%s\t%10.20f\n" % (chr, pos, id, pvalue))
            #
        #     chr_list[i] = genotypes.variants[i]['chromosome'][0]
        #     position_list[i] = genotypes.variants[i]['position'][0]
        #     ids_list[i] = genotypes.variants[i]['id'][0]
        # res_d = {'ps': p_vals, 'f_stats': f_stats, 'rss': rss_list, 'id': ids_list,
        #          'chr': chr_list, 'position': position_list}
        # return res_d
        f.close()
        return 0


def emma(genotypes, phenotypes, k, cofactors=None, f_test=True, lrt_test=True, maf=0.05, outPutFile="emmax_f_test.ps"):
    """
    Run EMMA
    """
    lmm = LinearMixedModel(phenotypes, k=k)
    if cofactors:
        for cofactor in cofactors:
            lmm.add_factor(cofactor)
    res = lmm.emma_test(genotypes, f_test=f_test, lrt_test=lrt_test, maf=maf, outPutFile=outPutFile)
    return res

ped = parse_plink_tped_file("/Users/song/Dropbox/mlmm_cpp/src/tests/testData/orf")
k1 = np.loadtxt("/Users/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf") # this is the kinship matrix
k1 = np.mat(k1)
#k = scale_k(k)
phe = parse_plink_fam_phenotype_file("/Users/song/Dropbox/mlmm_cpp/src/tests/testData/phenotype.tfam")

n_individual = len(phe.individuals)
individual_ids = []
for i in range(n_individual):
    individual_ids.append( phe.individuals[i]['individual_id'][0] )
ped.extractAccessions(individual_ids)

phenotype = phe.get_phenotypes_with_individuals_id(ped.individuals['individual_id'])
outputfile="test.ps"
emmax_result = emma(ped, phenotype, k1, maf=0.1, outPutFile=outputfile)

