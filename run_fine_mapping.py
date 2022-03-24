import pandas as pd
import numpy as np
import scipy.linalg
import time
import gzip
import logging
import os
import sys
from susieinf import susie,cred
from finemapinf import finemap
import argparse
import pkg_resources

def read_sumstats_from_file(z_file, args):
    logging.info('Reading summary statistics from file %s'%(z_file))
    z_df = pd.read_csv(z_file, delimiter='\s+')
    if args.z_col_name:
        z = z_df[args.z_col_name].values.flatten()
    elif args.beta_col_name and args.se_col_name:
        z = z_df[args.beta_col_name]/z_df[args.se_col_name]
        z = z.values.flatten()
    else:
        raise ValueError('Either column name for z score or column names for both marginal effect size and standard error must be provided')
    logging.info('%s SNPs in summary statistics file'%(len(z)))
    return z,z_df

def read_large_file(file_name):
    file_ext = os.path.splitext(file_name)[1]
    if file_ext == '.npy':
        return np.load(file_name)
    elif file_ext == '.npz':
        npz_file = np.load(file_name)
        file_keys = [k for k in npz_file.keys()]
        if len(file_keys) == 1:
            return npz_file[file_keys[0]]
        else:
            raise ValueError('%s contains multiple keys'%(file_name))
    elif file_ext in ['.gz', '.bgz']:
        with gzip.open(file_name, 'rb') as f:
            return np.loadtxt(f)
    else:
        raise ValueError('File extension %s of file %s currently not supported'%(file_ext, file_name))

def read_V_Dsq_from_file(V_file, Dsq_file):
    logging.info('Reading V and Dsq from file %s and %s'%(V_file,Dsq_file))
    t0 = time.time()
    V = read_large_file(V_file)
    Dsq = read_large_file(Dsq_file)
    logging.info('Reading V and Dsq files took  %0.2f seconds'%(time.time() - t0))
    # check data shape
    V_shape = V.shape
    Dsq_shape = Dsq.shape
    if V_shape[0]!=V_shape[1] or V_shape[0]!=Dsq_shape[0]:
        raise ValueError('Incorrect data shape for V or Dsq')
    return V,Dsq

def process_output(method_name, output_dict, df, output_prefix):
    if method_name == 'susieinf':
        df['prob'] = 1 - (1 - output_dict['PIP']).prod(axis=1)
        L = output_dict['PIP'].shape[1]
        alpha_cols = ['alpha{}'.format(i) for i in range(1,L+1)]
        mu_cols = ['mu{}'.format(i) for i in range(1,L+1)]
        df[alpha_cols + mu_cols + ['omega{}'.format(i) for i in range(1,L+1)]] = pd.DataFrame(np.concatenate((output_dict['PIP'], output_dict['mu'], output_dict['omega']), axis=1))
        df['alpha'] = output_dict['alpha']
        df['post_mean'] = np.sum(df[mu_cols].values * df[alpha_cols].values, axis=1) + df['alpha']
        df['tausq'] = output_dict['tausq']
        df['sigmasq'] = output_dict['sigmasq']
        df['cs'] = -1
        if len(output_dict['cred'])>0:
            df = df.reset_index(drop=True)
            for i,x in enumerate(output_dict['cred']): df.loc[x,'cs'] = i+1
        out_file = output_prefix+'.susieinf.gz'
        logging.info('Saving output to %s'%(out_file))
        df.to_csv(out_file, sep='\t', index=False, compression='gzip')
    elif method_name == 'finemapinf':
        df['prob'] = output_dict['PIP']
        df['post_mean_cond'] = output_dict['beta']
        df['post_sd_cond'] = output_dict['se']
        df['alpha'] = output_dict['alpha']
        df['post_mean'] = df['post_mean_cond'] * df['prob'] + df['alpha']
        df['tausq'] = output_dict['tausq']
        df['sigmasq'] = output_dict['sigmasq']
        out_file = output_prefix + '.finemapinf.gz'
        logging.info('Saving output to %s'%(out_file))
        df.to_csv(out_file, sep='\t', index=False, compression='gzip')
        config_df = pd.DataFrame(output_dict['models'], columns=['prob','config'])
        out_file = output_prefix + '.finemapinf.config.gz'
        logging.info('Saving FINEMAP-inf configurations to %s'%(out_file))
        config_df.to_csv(out_file, sep='\t', index=False, compression='gzip')

def susieinf_splash_screen():
    logging.info('*********************************************************************')
    logging.info('* SuSiE-inf')
    logging.info('* Version {}'.format(pkg_resources.require("susieinf")[0].version))
    logging.info('* (C) Zhou Fan')
    logging.info('*********************************************************************')

def finemapinf_splash_screen():
    logging.info('*********************************************************************')
    logging.info('* FINEMAP-inf')
    logging.info('* Version {}'.format(pkg_resources.require("finemapinf")[0].version))
    logging.info('* (C) Zhou Fan')
    logging.info('*********************************************************************')



if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # input file params
    parser.add_argument('--sumstats', type=str, required=True, help='Name of summary statistics file')
    parser.add_argument('--z-col-name', type=str, help='Name of z score column name')
    parser.add_argument('--beta-col-name', type=str, help='Marginal effect size column name, if z score column name is not provided. --se-col-name must also be specified')
    parser.add_argument('--se-col-name', type=str, help='Marginal standard error column name, if z score column name is not provided. --beta-col-name must also be specified')
    parser.add_argument('--ld-file', type=str, help='Name of LD matrix file')
    parser.add_argument('--V-file', type=str, help='Name of file containing eigenvectors of XtX')
    parser.add_argument('--Dsq-file', type=str, help='Name of file containing eigenvalues of XtX')

    # input algorithm params
    parser.add_argument('--meansq', type=float, default=1.0, help='Average squared magnitude of y (equal to yT.dot(y)/n)')
    parser.add_argument('--n', type=int, required=True, help='Sample size')
    parser.add_argument('--num-sparse-effects', type=int, default=10, help='Maximum number of sparse large effects in the region')
    parser.add_argument('--method', type=str, default='susieinf', help='Comma delimited fine-mapping methods to run, e.g. susieinf,finemapinf')
    parser.add_argument('--use-susieinf-tausq', action='store_true', help='Use SuSiE-inf generated tau-squared for FINEMAP-inf')
    parser.add_argument('--empirical-Bayes-method', type=str, default='moments', help='One of {"moments","MLE"}. Sigma-squared and tau-squared will be estimated using either method-of-moments or MLE')
    parser.add_argument('--prior', type=str, help='File name for user specified causal priors. Values must sum to 1. Algorithm default uses uniform prior.')
    parser.add_argument('--num-epochs', type=int, default=5, help='Number of epochs for FINEMAP-inf. Tau-squred and sigma-squared estimates can potentially gain accuracy with more epochs.')
    parser.add_argument('--coverage', type=float, default=0.95, help='Credible set coverage')
    parser.add_argument('--purity', type=float, default=0.5, help='Credible set purity threshold')

    # output params
    parser.add_argument('--eigen-decomp-prefix', type=str, help='Save V and Dsq from decomposing of XtX to .npz files with specified file prefix')
    parser.add_argument('--save-npz', action='store_true', help='Save output dictionary to .npz file')
    parser.add_argument('--save-tsv', action='store_true', help='Save output as formatted .tsv file')
    parser.add_argument('--output-prefix', type=str, required=True, help='File prefix for output files')
    args = parser.parse_args()

    # configure logger
    logging.basicConfig(filename=args.output_prefix+'.log', filemode='w', level=logging.DEBUG,
            format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    # read in summary statistics
    z,z_df = read_sumstats_from_file(args.sumstats, args)

    # read in V and Dsq or reading in LD and performing eigen decomposition
    if args.ld_file is None:
        if args.V_file is None or args.Dsq_file is None:
            raise ValueError('--V-file and --Dsq-file must be provided if --ld-file is not provided')
        else:
            V,Dsq = read_V_Dsq_from_file(args.V_file, args.Dsq_file)
    if args.ld_file:
        if args.V_file is not None or args.Dsq_file is not None:
            logging.warning('Using --V-file and --Dsq-file instead of --ld-file. Please do not specify --ld-file if V and Dsq are available')
            V,Dsq = read_V_Dsq_from_file(args.V_file, args.Dsq_file)
        else:
            logging.info('Reading in LD matrix from file %s'%(args.ld_file))
            t0 = time.time()
            LD = read_large_file(args.ld_file)
            logging.info('Reading in LD matrix took %0.2f seconds'%(time.time() - t0))
            if len(LD)!=len(z): raise ValueError('Size of LD matrix does not match summary statistics')
            logging.info('Performing eigen decomposition')
            t0 = time.time()
            eigenvals,V = scipy.linalg.eigh(LD, driver='evd')
            logging.info('Eigen decomposition took %0.2f seconds'%(time.time() - t0))
            Dsq = args.n * eigenvals
            if args.eigen_decomp_prefix is not None:
                np.savez_compressed(args.eigen_decomp_prefix+'.Dsq.npz', Dsq)
                np.savez_compressed(args.eigen_decomp_prefix+'.V.npz', V)
    if len(V)!=len(z): raise ValueError('Size of V and Dsq do not match summary statistics')

    # read in priors
    if args.prior:
        prior_df = pd.read_csv(args.prior, delimiter='\s+')
        pi0 = prior_df.values.flatten()
        if len(pi0)!=len(z): raise ValueError('Length of prior file does not match length of summary statistics')
    else:
        pi0=None

    # run fine-mapping algorithms
    methods = args.method.split(',')
    logging.info('Running %s'%(args.method))
    if 'susieinf' in methods:
        t0_susieinf = time.time()
        susieinf_splash_screen()
        susie_output = susie(z, args.meansq, args.n, args.num_sparse_effects, LD=None, V=V, Dsq=Dsq,
                est_ssq=True,ssq=None,ssq_range=(0,1),pi0=pi0, method=args.empirical_Bayes_method,
                sigmasq_range=None,tausq_range=None,PIP=None,mu=None,maxiter=100,PIP_tol=1e-3,verbose=True)
        susie_output['cred'] = cred(susie_output['PIP'], coverage=args.coverage, purity=args.purity, LD=None,V=V, Dsq=Dsq, n=args.n)
        logging.info('Running SuSiE-inf took %0.2f seconds'%(time.time() - t0_susieinf))
        if args.save_npz:
            out_file = args.output_prefix+'.susieinf.npz'
            logging.info('Saving output dictionary to %s'%(out_file))
            np.savez_compressed(out_file, **susie_output)
        if args.save_tsv:
            process_output('susieinf', susie_output, z_df.copy(), args.output_prefix)
    if 'finemapinf' in methods:
        t0_finemapinf = time.time()
        if args.use_susieinf_tausq:
            logging.info('Using tau-squared generated by SuSiE-inf: %.2E'%(susie_output['tausq']))
            finemapinf_splash_screen()
            finemap_output = finemap(z, args.meansq, args.n, args.num_sparse_effects, LD=None, V=V, Dsq=Dsq,
                sigmasq=1, tausq=susie_output['tausq'], pi0=pi0, S=[0.05**2], sched_sss=[10000], n_conv_sss=100,
                nmodels=5000, prob_tol_sss=1e-3, seed=21,models_nbins=None,pairs_nbins=None,verbose=1)
        else:
            finemapinf_splash_screen()
            finemap_output = finemap(z, args.meansq, args.n, args.num_sparse_effects, LD=None, V=V, Dsq=Dsq,
                    sigmasq=1, tausq=0, pi0=pi0, S=[0.05**2], sched_sss=[100]*(args.num_epochs-1)+[10000], n_conv_sss=100,
                    nmodels=5000, prob_tol_sss=1e-3, seed=21,models_nbins=None,pairs_nbins=None,verbose=1)
        logging.info('Running FINEMAP-inf took %0.2f seconds'%(time.time() - t0_finemapinf))
        if args.save_npz:
            out_file = args.output_prefix+'.finemapinf.npz'
            logging.info('Saving output dictionary to %s'%(out_file))
            finemap_output['models'] = np.asanyarray(finemap_output['models'], dtype=object)
            np.savez_compressed(out_file, **finemap_output)
        if args.save_tsv:
            process_output('finemapinf', finemap_output, z_df.copy(), args.output_prefix)

