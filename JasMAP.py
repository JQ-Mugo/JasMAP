#!/usr/bin/python

import numpy as np
import pandas as pd
import os, sys, fileinput,time,datetime
import logging
import warnings
logger = logging.getLogger()
logger.setLevel(logging.INFO)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logger.addHandler(console)
warnings.filterwarnings("ignore")
import scipy.stats as sts
from scipy import linalg,stats
from scipy import optimize

class jasMapUtils:
    '''     Useful functions    '''
    def __init__(self):
        return None

    def logical(self, Value):
        '''Convert a string to logical variable 
        Value: a string, NO or YES
        Return a logical variable'''

        if Value == "YES":
                return True
        else:
                return False

    def safe_open(self,filename):
        '''Opens a file in a safe manner or raise an exception if the file doesn't exist filename: a string, name of a file
        Return: file object'''

        if os.path.exists(filename):
                return open(filename)
        else:
                self.printHead()
                sys.stderr.write("ERROR ==> No such file or directory\"%s\"\n"%format(filename))
                return False

    def check_files(self,params):
        '''Check existed of files in a dict or a list of  params: a dictionary, of option:value '''
        if len(params) == 1:
                if os.path.exists(params[0]) == False:
                        sys.stderr.write("ERROR ==> No such file or directory \"%s\"\n\n" %params[0])
                        sys.exit(1)
        else:
                if params[0] in ['covariates_file'] and os.path.exists(params[1]) == True:
                        self.covariates = True
                else:
                        self.covariates = False

    def terminate(self):
        '''Terminate the process'''
        log_rem = os.getcwd()
        try:
                pkl_files = [fil for fil in os.listdir(self.outfolder) if fil.endswith('pkl')]
                for fil in pkl_files:
                        pass
                os.system("cp"+" "+log_rem+'/'+self.logFile+" "+self.outfolder)
        except:
                filenames = os.listdir(os.curdir)
                for filename in filenames:
                        if os.path.isfile(filename) and filename.endswith('.log'):
                                os.system("rm"+" "+ log_rem+"/*.log")
        finally:
                log = log_rem+'/'+self.logFile
                logger.info('Log file generated in '+log+'\nHave a nice day!\n')

                filenames = os.listdir(os.curdir)
                for filename in filenames:
                        if os.path.isfile(filename) and filename.endswith('.log'):
                                pass
                sys.exit(1)


class jasMapInit(jasMapUtils):

    try:
        logger.removeHandler(logger.handlers[1])
    except:
        pass
    log_rem = os.getcwd()
    logger.setLevel(logging.INFO)
    filenames = os.listdir(os.curdir)
    for filename in filenames:
        if os.path.isfile(filename) and filename.endswith('.log'):
                os.system("rm"+" "+log_rem+"/*.log")

    logFile = 'jasmap-'+str(time.time()).split('.')[0]+'.log'
    fh = logging.FileHandler(logFile, mode='w')
    logger.addHandler(fh)


    logger.info("\n***********************************************************************************************************************")
    logger.info("               JASMAP: A Joint Ancestry and SNP Mapping Method for Multi-way Admixed Populations                        ")
    logger.info("                                        Computational Biology Group                                                    ")
    logger.info("                               Intergrative Biomedical Sciences Department (IBMS)                                     ")
    logger.info("                                   2023, University of Cape Town, South Africa                                       ")
    logger.info("                                            Verson 1.0 Beta                                                         ")
    logger.info(  "***********************************************************************************************************************\n")

    def __init__(self, argv, logFile=logFile):
        '''Initializing JasMAP by reading the parameter file'''
        self.argv = [argv]
        self.logFile = logFile
        popnames = []

        if len(self.argv) == 0 or self.argv == ['']:
                logger.info('Command line usage: %s <parameter file>  ' % sys.argv[0])
                logger.info('eg: python JasMAP.py parameter_file.txt\n')
                sys.exit(1)
        elif len(self.argv) == 1:
                try:
                        # Check if parameter file is available
                        self.paramFile = self.argv[0]
                        if os.path.exists(os.getcwd()+'/'+self.paramFile):
                                inFile = open(os.getcwd()+'/'+self.paramFile)
                        elif os.path.exists(self.paramFile):
                                inFile = open(self.paramFile)
                        else:
                                logger.info('\nERROR ==> Failed to process the input; No parameter file!\n')
                                self.terminate()
                                sys.exit()

                        # Read in the parameter file and set logical parameters true/false
                        self.params_dict = {}
                        keyList = []; logicalList = []
                        for line in inFile:
                                if not line.startswith("#"):
                                        data = line.split()
                                        if len(data) != 0:
                                                data1,data2 = data[0].split(":")
                                                if data1 in logicalList:
                                                        data2 = self.logical(data2)
                                                self.params_dict[data1] = data2
                                                keyList.append(data1)

                        # check all parameters are listed and valid
                        if len(self.params_dict) != 9 :
                                sys.stderr.write('ERROR ==> Missing parameters!!! Failed to process the input, check the parameters!\n')
                                sys.exit()

                        expected_params_list = ['infolder','lanc_file','geno_file','pheno_file','covariates_file','no_of_ancs','AM_eff_no_tests','no_of_snps','outfolder']
                        unknown_params_list = []
                        for i in range(len(keyList)):
                                if keyList[i] not in expected_params_list:
                                        unknown_params_list.append(keyList[i])

                        if len(unknown_params_list) != 0:
                                sys.stderr.write('\nInvalid option: '+','.join(unknown_params_list))
                                sys.stderr.write('\nERROR ==> Failed to process the input, check the parameters!\n')
                                self.terminate()
                                sys.exit()


                        #Check/create output folder
                        try:
                                self.outfolder = self.params_dict['outfolder']+"/"
                                path = os.getcwd()
                                if os.path.exists(self.outfolder):
                                        pass
                                else:
                                        os.makedirs(self.outfolder)
                        except IndexError:
                                sys.stderr.write('ERROR ==> Cannot create directory. Please create a directory OUT in your working directory\n')
                                sys.exit(1)


                        #assign complete input-file-paths to the input files and check existence
                        self.infolder = self.params_dict['infolder']
                        self.check_files([self.infolder])
                        for file in ['lanc_file','geno_file','pheno_file']: #compulsory files
                                self.check_files([self.infolder+"/"+self.params_dict[file]])
                                self.params_dict[file] = self.infolder+"/"+self.params_dict[file]

                        #Additional files (covariates and other factors)
                        for file2 in ['covariates_file']:
                                self.check_files(['covariates_file',self.infolder+"/"+self.params_dict[file2]])
                                self.params_dict[file2] = self.infolder+"/"+self.params_dict[file2]

                        #create global variable options
                        self.lanc_file = self.params_dict['lanc_file']
                        self.geno_file = self.params_dict['geno_file']
                        self.pheno_file = self.params_dict['pheno_file']
                        self.covar_file = self.params_dict['covariates_file']
                        self.no_of_ancs = int(self.params_dict['no_of_ancs'])
                        self.no_of_snps = int(self.params_dict['no_of_snps'])
                        self.AM_eff_no_tests = float(self.params_dict['AM_eff_no_tests'])
                        self.lanc_labels_list = [i+1 for i in range(self.no_of_ancs)]

                        #create a final printing dictionary
                        self.Params = {}
                        i=0
                        for param in keyList:
                                self.Params[i] = [param,self.params_dict[param]]
                                i+=1

                except (IndexError, TypeError):
                        sys.stderr.write('ERROR ==> Failed to process the input, check the parameters!\n\n')
                        self.terminate()
                        sys.exit(1)
        else:
                logger.info('Command line usage: %s <parameter file>  ' % sys.argv[0])
                logger.info('eg. python JasMAP.py parameters_file_jasmap.txt\n')
                self.terminate()
                sys.exit(1)


class LMM:
    """
         The LMM class is from pylmm by Nicholas A. Furlotte, used and modified under the the terms of the GNU Affero General Public License 
         # pylmm is a python-based linear mixed-model solver with
         # applications to GWAS Copyright (C) 2015 Nicholas A. Furlotte
         # (nick.furlotte@gmail.com)

         #    This program is free software: you can redistribute it and/or modify
         #    it under the terms of the GNU Affero General Public License as published by
         #    the Free Software Foundation, either version 3 of the License, or
         #    (at your option) any later version.

         #    This program is distributed in the hope that it will be useful,
         #    but WITHOUT ANY WARRANTY; without even the implied warranty of
         #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
         #    GNU General Public License for more details.

         #    You should have received a copy of the GNU General Public License
         #    along with this program.  If not, see <http://www.gnu.org/licenses/>.
      	 
         This is a simple version of EMMA/fastLMM.
         The main purpose of this module is to take a phenotype vector (Y), a set of covariates (X) and a kinship matrix (K)
         and to optimize this model by finding the maximum-likelihood estimates for the model parameters.
         There are three model parameters: heritability (h), covariate coefficients (beta) and the total
         phenotypic variance (sigma).
         Heritability as defined here is the proportion of the total variance (sigma) that is attributed to
         the kinship matrix.

         For simplicity, we assume that everything being input is a numpy array.
         If this is not the case, the module may throw an error as conversion from list to numpy array
         is not done consistently.

    """
   
    def __init__(self,Y,K,Kva=[],Kve=[],X0=None,verbose=False):
      
        """
            The constructor takes a phenotype vector or array of size n.
            It takes a kinship matrix of size n x n.  Kva and Kve can be computed as Kva,Kve = linalg.eigh(K) and cached.
            If they are not provided, the constructor will calculate them.
            X0 is an optional covariate matrix of size n x q, where there are q covariates.
            When this parameter is not provided, the constructor will set X0 to an n x 1 matrix of all ones to represent a mean effect.
        """

        if str(X0) == str(None):
                X0 = np.ones(len(Y)).reshape(len(Y),1)
        self.verbose = verbose
      
        x = True ^ np.isnan(Y)
        x = x.reshape(-1,)
        if not x.sum() == len(Y):
                if self.verbose: sys.stderr.write("Removing %d missing values from Y\n" % ((True ^ x).sum()))
                Y = Y[x]
                K = K[x,:][:,x]
                X0 = X0[x,:]
                Kva = []
                Kve = []
        self.nonmissing = x

        if len(Kva) == 0 or len(Kve) == 0:
                if self.verbose: sys.stderr.write("Obtaining eigendecomposition for %dx%d matrix\n" % (K.shape[0],K.shape[1]) )
                begin = time.time()
                if np.isnan(K).any():K = np.nan_to_num(K)
                Kva,Kve = linalg.eigh(K)
                end = time.time()
                if self.verbose: sys.stderr.write("Total time: %0.3f\n" % (end - begin))

        self.K = K
        self.Kva = Kva
        self.Kve = Kve
        self.N = self.K.shape[0]
        self.Y = Y.reshape((self.N,1))
        self.X0 = X0

        if sum(self.Kva < 1e-6):
                if self.verbose: sys.stderr.write("Cleaning %d eigen values\n" % (sum(self.Kva < 0)))
                self.Kva[self.Kva < 1e-6] = 1e-6
        self.transform()

    def transform(self):

        """
            Computes a transformation on the phenotype vector and the covariate matrix.
            The transformation is obtained by left multiplying each parameter by the transpose of the
            eigenvector matrix of K (the kinship).
        """

        self.Yt = np.matmul(self.Kve.T, self.Y)
        self.X0t = np.matmul(self.Kve.T, self.X0)
        self.X0t_stack = np.hstack([self.X0t, np.ones((self.N,1))])
        self.q = self.X0t.shape[1]

    def getMLSoln(self,h,X,beta):
        """
            Obtains the maximum-likelihood estimates for the covariate coefficients (beta),
            the total variance of the trait (sigma) and also passes intermediates that can
            be utilized in other functions. The input parameter h is a value between 0 and 1 and represents
            the heritability or the proportion of the total variance attributed to genetics.  The X is the
            covariate matrix.
        """

        S = 1.0/(h*self.Kva + (1.0 - h))
        Xt = X.T*S
        XX = np.matmul(Xt,X)
        if np.isnan(XX).any(): XX = np.nan_to_num(XX)
        try:
                XX_i = linalg.inv(XX)
        except:
	        return "pass"
        if type(beta) == type(None):
                beta = np.matmul(np.matmul(XX_i,Xt),self.Yt)
        else:
                Xf = self.X0t
                Xft = Xf.T*S
                XXf = np.matmul(Xft,Xf)
                XXf_i = linalg.inv(XXf)
                betaf = np.matmul(np.matmul(XXf_i,Xft),self.Yt)
                beta = np.append(betaf,beta,axis=0)
        Yt = self.Yt - np.matmul(X,beta)
        Q = np.dot(Yt.T*S,Yt)
        sigma = Q * 1.0 / (float(self.N) - float(X.shape[1]))
        return beta,sigma,Q,XX_i,XX

    def LL_brent(self,h,X=None,REML=False):
        if h < 0: return 1e6
        return -self.LL(h,X,stack=False,REML=REML)[0]


    def LL(self,h,X=None,beta=None,stack=True,REML=False):
        """
            Computes the log-likelihood for a given heritability (h).  If X==None, then the
            default X0t will be used.  If X is set and stack=True, then X0t will be matrix concatenated with
            the input X.  If stack is false, then X is used in place of X0t in the LL calculation.
            REML is computed by adding additional terms to the standard LL and can be computed by setting REML=True.
        """

        if str(X) == str(None): X = self.X0t
        elif stack:
                self.X0t_stack[:,(self.q)] = np.matmul(self.Kve.T,X)[:,0]
                X = self.X0t_stack
        n = float(self.N)
        q = float(X.shape[1])
        MLSoln = self.getMLSoln(h,X,beta)
        if MLSoln == "pass":
                return "pass"
        else:
                beta,sigma,Q,XX_i,XX = MLSoln

        LL = n*np.log(2*np.pi) + np.log(h*self.Kva + (1.0-h)).sum() + n + n*np.log(1.0/n * Q)
        LL = -0.5 * LL

        if REML:
                LL_REML_part = q*np.log(2.0*np.pi*sigma) + np.log(linalg.det(np.matmul(X.T,X))) - np.log(linalg.det(XX))
                LL = LL + 0.5*LL_REML_part

        LL = LL.sum()

        return LL,beta,sigma,XX_i

    def getMax(self,H, X=None,REML=False):

        """
              Helper functions for .fit(...).
              This function takes a set of LLs computed over a grid and finds possible regions
              containing a maximum.  Within these regions, a Brent search is performed to find the
              optimum.
        """
        n = len(self.LLs)
        HOpt = []
        for i in range(1,n-2):
                if self.LLs[i-1] < self.LLs[i] and self.LLs[i] > self.LLs[i+1]:
                        HOpt.append(optimize.brent(self.LL_brent,args=(X,REML),brack=(H[i-1],H[i+1])))
                        if np.isnan(HOpt[-1]): HOpt[-1] = H[i-1]

        if len(HOpt) > 1:
                if self.verbose: sys.stderr.write("NOTE: Found multiple optima.  Returning first...\n")
                return HOpt[0]
        elif len(HOpt) == 1: return HOpt[0]
        elif self.LLs[0] > self.LLs[n-1]: return H[0]
        else: return H[n-1]

    def fit(self,X=None,ngrids=100,REML=True):

        """
             Finds the maximum-likelihood solution for the heritability (h) given the current parameters.
             X can be passed and will transformed and concatenated to X0t.  Otherwise, X0t is used as
             the covariate matrix.

             This function calculates the LLs over a grid and then uses .getMax(...) to find the optimum.
             Given this optimum, the function computes the LL and associated ML solutions.

        """

        if str(X) == str(None): X = self.X0t
        else:
                self.X0t_stack[:,(self.q)] = np.matmul(self.Kve.T,X)[:,0]
                X = self.X0t_stack
        H = np.array(range(ngrids)) / float(ngrids)
        L = np.array([self.LL(h,X,stack=False,REML=REML)[0] for h in H])
        self.LLs = L

        hmax = self.getMax(H,X,REML)

        L,beta,sigma,betaSTDERR = self.LL(hmax,X,stack=False,REML=REML)

        m = np.shape(self.K)[0]

        self.H = H
        self.optH = hmax.sum()
        self.optLL = L
        self.optBeta = beta

        self.optSigma = sigma.sum()
        return hmax,beta,sigma,L

    def association(self,X, h = None, beta=None, stack=True,REML=True, returnBeta=False):

        """
            Calculates association statitics for the SNPs encoded in the vector X of size n.
            If h == None, the optimal h stored in optH is used.

        """
        if stack:
                self.X0t_stack[:,(self.q)] = np.matmul(self.Kve.T,X)[:,0]
                X = self.X0t_stack
        if h == None: h = self.optH

        L,beta,sigma,betaVAR = self.LL(h,X,beta,stack=False,REML=REML)
        q  = len(beta)

        if beta == "a": return "nan","nan","nan","nan"
        
        stderr,ts,ps = self.tstat(beta[q-1],betaVAR[q-1,q-1],sigma,q)

        if returnBeta: return beta[q-1].sum(),stderr,ts,ps

        return beta[q-1].sum(),stderr,ts,ps
      
    def tstat(self,beta,var,sigma,q,log=False):

        """
            Calculates a t-statistic and associated p-value given the estimate of beta and its standard error.
            This is actually an F-test, but when only one hypothesis is being performed, it reduces to a t-test.
        """
        stderr = np.sqrt(var * sigma)
        ts = beta / stderr
        pvalue = 2.0*(stats.norm.sf(np.abs(ts)))

        if log:
                ps = 2.0 + (stats.t.logsf(np.abs(ts), self.N-q))
        else:
                ps = 2.0*(stats.t.sf(np.abs(ts), self.N-q))
        if not len(ts) == 1 or not len(ps) == 1:
                raise Exception("Something bad happened :(")
        return stderr,ts.sum(),ps.sum()

    def matrixMult(self,A,B):
        "A function to multiply matrices"
        #If there is no fblas then we will revert to np.dot()
        try:
                linalg.fblas
        except AttributeError:
                return np.dot(A,B)

        # If the matrices are in Fortran order then the computations will be faster
        # when using dgemm.  Otherwise, the function will copy the matrix and that takes time.
        if not A.flags['F_CONTIGUOUS']:
                AA = A.T
                transA = True
        else:
                AA = A
                transA = False

        if not B.flags['F_CONTIGUOUS']:
                BB = B.T
                transB = True
        else:
                BB = B
                transB = False
        return linalg.fblas.dgemm(alpha=1.,a=AA,b=BB,trans_a=transA,trans_b=transB)

    def calculateKinship(self,K,center=False):
        """
        W is an n x m matrix encoding SNP minor alleles.

        This function takes a matrix oF SNPs and obtaining the Kinship matrix 
        """
        starting_time = datetime.datetime.today()
        W = np.array([i for i in K])
        n = W.shape[0]
        m = W.shape[1]
        keep = []
        for i in range(m):
                mn = W[True ^ np.isnan(W[:,i]),i].mean()
                W[np.isnan(W[:,i]),i] = mn
                vr = W[:,i].var()
                if vr == 0: continue

                keep.append(i)
                W[:,i] = (W[:,i] - mn) / np.sqrt(vr)

        W = W[:,keep]
        K = self.matrixMult(W,W.T) * 1.0/float(m)
        if center:
                P = np.diag(np.repeat(1,n)) - 1/float(n) * np.ones((n,n))
                S = np.trace(self.matrixMult(self.matrixMult(P,K),P))
                K_n = (n - 1)*K / S
                finishing_time = datetime.datetime.today()
                logger.info('Obtaining the kinship matrix took :%s' % str(finishing_time-starting_time))
                return K_n
        return K

class  Joint_Association(jasMapInit,LMM):

    '''This module reads the input files and implements the joint association.'''

    def read_lanc(self):
        '''Read the local ancestry inference input and split into different ancestral populations matrices'''
        starting_time = datetime.datetime.today()
        ###self.split_anc_per_SNP_file = open(self.outfolder+"ancestry_per_SNP.txt","wt")
        self.lanc_dict = {} ; self.genotype_strata = {} ; self.phenotype_strata = {}


        for line in fileinput.input(self.lanc_file):
                data = line.rstrip().split(" ")
                if len(data[0])>1: data = list(line.rstrip())
                if all(int(i) <= self.no_of_ancs and int(i)>0 for i in data) == False:
                        sys.stderr.write('ERROR ==> Unknown anc label in line %i!\n\n'%(fileinput.lineno()))
                        self.terminate()
                else:
                        self.lanc_dict[fileinput.lineno()] = {}
                        for idx in range(self.no_of_ancs):
                                self.lanc_dict[fileinput.lineno()][idx+1] = []
                        ind = 0
                        for hap in range(0,len(data),2):
                                if data[hap] == data[hap+1]:
                                        for anc in self.lanc_labels_list:
                                                if int(data[hap]) == anc:
                                                        self.lanc_dict[fileinput.lineno()][anc].append(2)
                                                else:
                                                        self.lanc_dict[fileinput.lineno()][anc].append(0)
                                else:
                                        for anc2 in self.lanc_labels_list:
                                                if int(data[hap]) == anc2 or int(data[hap+1]) == anc2 :
                                                        self.lanc_dict[fileinput.lineno()][anc2].append(1)
                                                else:
                                                        self.lanc_dict[fileinput.lineno()][anc2].append(0)
        finishing_time = datetime.datetime.today()
        logger.info("Reading/converting ancestry file took %s time to finish.\n\n\n"%str(finishing_time-starting_time))

    def get_ancestry_AM_instances(self,X0s):
        '''Splitting the self.lanc_dict per anc, Get AM LMM instances, ancestry relatedness matrices and ancestry 
           heritability'''
        self.lanc_dict_per_anc = {}
        for anc in range(1,self.no_of_ancs+1):
                self.lanc_dict_per_anc[anc] = []
        for snp in self.lanc_dict:
                for anc_idx in self.lanc_dict[snp]:
                        self.lanc_dict_per_anc[anc_idx].append(self.lanc_dict[snp][anc_idx]) 

        self.lai_h_dict = {} ; self.lai_L_instances = {}
        for anc2 in range(1,self.no_of_ancs+1):
                lai_dict_per_anc_df = pd.DataFrame(self.lanc_dict_per_anc[anc2])
                lai_dict_per_anc_edited_df = lai_dict_per_anc_df.drop_duplicates()
                lai_dict_per_anc = [list(lai_dict_per_anc_edited_df.iloc[i]) for i in list(range(len(lai_dict_per_anc_edited_df.index.values)))]
                lanc_dict_per_anc_T = np.array(lai_dict_per_anc).T
                lai_kinship = self.calculateKinship(K=lanc_dict_per_anc_T)
                lai_L = LMM(self.phenotypes,lai_kinship,X0s)
                lai_hmax,lai_beta,lai_sigma,lai_LL_null = lai_L.fit()
                self.lai_h_dict[anc2] = lai_hmax ; self.lai_L_instances[anc2] = lai_L
                ##logger.info('Pseudo-heritability for ancestry %s = %s.\n\n\n'%(str(anc2),str(lai_hmax)))

    def read_covar_file(self):
        '''Reading in the covariates file PCs + others'''
        covar_list = []
        if self.covariates:
                for line in fileinput.input(self.covar_file):
                        if fileinput.lineno() > 0:
                                data = line.rstrip().split()[2:]
                                covar_list.append(np.array([float(i) for i in data]))
                self.covar_list = np.array(covar_list)
        else:
                self.covar_list = np.ones(self.inds_size)
        logger.info("Reading PCA file done.")

    def read_geno_file(self):
        '''Reading the genotype information provided {0,1,2}'''
        #print("Reading genotypes...")
        self.genotypes = {}
        for line in fileinput.input(self.geno_file):
                #print("Reading geno file...",fileinput.lineno())
                data = line.rstrip().split()
                if len(data)==1:
                        data = list(data[0])
                geno = np.array([int(i) for i in data])
                geno_mean = np.mean(geno) ; geno_var = np.var(geno)
                self.genotypes[fileinput.lineno()] = (geno - geno_mean)/geno_var
        self.genotypes_T = np.array([list(i) for i in self.genotypes.values()]).T
        logger.info("Reading genotype file done.")

    def read_pheno_file(self):
        '''Read in the phenotype file format {2,1} in one column'''
        phenotypes = [] ; self.case_idx = [] ; self.control_idx = []
        for line in fileinput.input(self.pheno_file):
                data = line.rstrip()
                phenotypes.append(int(data))
                if data == "1":
                        self.control_idx.append(fileinput.lineno()-1)
                else:
                        self.case_idx.append(fileinput.lineno()-1)
        self.no_of_cases = phenotypes.count(2) ; self.no_of_controls = phenotypes.count(1)
        self.phenotypes = np.array(phenotypes)
        logger.info("Reading phenotype file done.")


    def get_ll(self,L,beta,hmax,true_geno):
        '''This function returns the likelihood of the data given the parameter beta'''

        X_t = list(true_geno)
        Xt = np.array([X_t]).T
        q = int(Xt.shape[1])
        beta = np.array([beta])
        out = L.LL(h=hmax,X=Xt,beta=beta,stack=True,REML=True)
        if out == "pass":
                return "pass"
        else:
                ll = out[0] ; stderr = np.sqrt(out[-2]*out[-1][q,q])
                return ll,stderr

    def run_joint_analysis(self):
        '''The main function where estimation of joint beta is done and the corresponding information written on the output file'''

        starting_time = datetime.datetime.today()

        #Output file specified, and intiliazed for the summary stats
        joint_fout = open(self.outfolder+"Jasmap.Joint.assoc","wt")
        joint_fout.writelines("SNP"+"\t"+"anc"+"\t"+"AM_beta"+"\t"+"AM_Pval"+"\t"+"SNP_beta"+"\t"+"SNP_Pval"+"\t"+"Joint_PPA"+"\n")

        anc_prior_fout_dict={}
        for idx in range(1,self.no_of_ancs+1):
                fout = open(self.outfolder+"Jasmap.ANC"+str(idx)+".AM.assoc","wt")
                anc_prior_fout_dict[idx] = fout
                anc_prior_fout_dict[idx].writelines("SNP"+"\t"+"beta"+"\t"+"Pval"+"\n")

        #Read other inputs, genotypes, phenotype and covariates
        self.read_geno_file()
        self.read_pheno_file() ; self.inds_size = len(self.phenotypes)
        self.read_covar_file() ; X0 = self.covar_list.T
        X0s = self.covar_list #np.ones((self.inds_size,1))
        #np.append(X0s,self.covar_list)

        #Read local ancestry file and separate to different ancestries and ancestry matrices
        self.read_lanc()
        self.get_ancestry_AM_instances(X0s)

        #Obtain the kinship matrix
        kS_time = datetime.datetime.today()
        self.kinship = self.calculateKinship(K=self.genotypes_T)
        kF_time = datetime.datetime.today()
        logger.info('Obtaining kinship matrix took :%s' % str(kF_time-kS_time))

        #Call an instance of the LMM class.
        logger.info("Calling instance of LMM ...")
        #X0s = np.ones((self.inds_size,1))
        #np.append(X0s,self.covar_list)
        L = LMM(self.phenotypes,self.kinship,X0s)

        #Obtain the pseudo heritability parameter (rho:p) from the function fit() in lmm class
        hmax,beta_null,sigma,LL_null = L.fit()
        ##logger.info('Pseudo-heritability for genotype = %s.\n\n\n'%str(hmax))

        ###########################################################################################
        #                    MAIN ANALYSIS
        ###########################################################################################

        # The analysis runs per SNPs by reading through the LAI dictionary with the LAI inference per SNP per ancestry
        # structure: self.lanc_dict = {SNP_1:{ANC_1:[LAI_1],ANC_2:[LAI_2],...,ANC_K:[LAI_K]},...,SNP_K:{...}}

        logger.info("Main analysis running ...")

        for snp in self.lanc_dict:
                # The true genotype for that SNP
                true_geno = self.genotypes[snp]

                #Initial dicts and lists for storage
                beta_AM_dict = {} ; stderr_dict = {} ; pval_dict = {}
                beta_AM_list = [] ; stderr_list = [] ; pval_list = [] ; norm_anc_list = []

                ############## PART I: ADMIXTURE MAPPING  #####################################

                # We loop through the different ancestries
                for anc in range(1,self.no_of_ancs+1):

                        '''First: For each SNP j perform admixture mapping case\control in the prior and obtain corresponding p-values'''
                        # Obtain the phenotype as 1's and 0's, we read them in as 2's and 1's
                        # Will edit this later to require 1's and 0's input
                        pheno = np.array(self.phenotypes)-1

                        # Normalize ancestry k, for SNP j using local ancestry inferences averages (lai_ave) and its variance
                        #norm_anc_list = [] ; beta_prob_list = []
                        lai_ave = np.sum(self.lanc_dict[snp][anc])/float(2*len(self.lanc_dict[snp][anc]))
                        lai_var = 2*lai_ave*(1-lai_ave)


                        # If the lai-average for that SNP <=0.01, or var is 0 then we skip performing association, to avoid spurious
                        # associations
                        if lai_ave < 0.01 or lai_var==0.0:
                                beta_AM = np.nan ;  beta_pval = np.nan ; stderr = np.nan
                        else:
                                norm_anc = (self.lanc_dict[snp][anc] - lai_ave)/lai_var

                                # Call the logistic prior function, with phenotype, normalized ancestry and
                                # the initially calculated intecept as input.
                                At = norm_anc.T ; At.shape = (len(At),1); lai_hmax = self.lai_h_dict[anc] ; lai_L = self.lai_L_instances[anc]
                                Abeta,Astderr,Ats,Aps = lai_L.association(X=At, h=lai_hmax, beta=None, stack=True,REML=True, returnBeta=False)
                                if Aps == 'nan':
                                        beta_AM = np.nan ;  beta_pval = np.nan ; stderr = np.nan
                                else:
                                        beta_AM = Abeta ;  beta_pval = Aps ; stderr = Astderr

                        # Store each ancestry stats in list; the obtained betas, stderr and p_value
                        beta_AM_list.append(beta_AM) ; stderr_list.append(stderr) ; pval_list.append(beta_pval)

                        # Output the Admixture mapping results for each SNP for each ancestry

                        anc_prior_fout_dict[anc].writelines(str(snp)+"\t"+str(beta_AM)+"\t"+str(beta_pval)+"\n")


                ############# PART II: Joint Admixture and SNP association  #####################

                # Obtain the ancestry from the logistic Admixture mapping with lowest p-value, i.e most significant.
                nan_pval_idx = [i for i,p in enumerate(pval_list) if np.isnan(p)]
                pval_list_updated = [j for j in pval_list if not np.isnan(j)]
                beta_AM_list_updated = [b for b in beta_AM_list if beta_AM_list.index(b) not in nan_pval_idx]

                #If none of the ancestries for that SNP produces a p-value. Then we continue with no ancestry effect p_value ~ 1
                # and record 'nan' as the most informative ancestry and proceed. Otherwise record the choosen ancestry and assiciated
                # summary stats.
                if len(pval_list_updated) == 0:
                        min_pval = 0.9999
                        joint_fout.writelines(str(snp)+"\t"+"nan"+"\t"+"nan"+"\t"+"nan"+"\t")
                else:
                     min_pval = np.min(pval_list_updated)
                     min_pval_idx = pval_list.index(min_pval)
                     am_beta = beta_AM_list[min_pval_idx]
                     joint_fout.writelines(str(snp)+"\t"+str(min_pval_idx+1)+"\t"+str(am_beta)+"\t"+str(min_pval)+"\t")

                ### 2.1 Convert the chosen p-value to chisq stats, then to pdf under the null and under the alternative hypothesis.
                am_burden = self.AM_eff_no_tests  #Input by user
                am_lambda = (sts.norm.ppf(1-0.05/float(am_burden)/2) + sts.norm.ppf(0.8))**2
                am_prior = 1/float(am_burden)
                am_chisq = sts.chi2.ppf(1-np.array(min_pval),1) #convert to chisq
                if np.isinf(am_chisq): am_chisq=68.7632522116684
                # Convert chisq to density functions
                am_pdf_lambda = sts.ncx2.pdf(am_chisq,1,nc=am_lambda) ; am_pdf_0 = sts.chi2.pdf(am_chisq,1)

                #Obtain admixture mapping posterior probability of association (AM_POST)
                am_post = (am_pdf_lambda*am_prior)/float((am_pdf_lambda*am_prior)+(am_pdf_0*(1-am_prior)))

                #### 2.2 Perform the SNP association using LMM
                X_t = list(true_geno)
                Xt = np.array([X_t]).T
                beta,stderr,ts,ps = L.association(X=Xt, h=hmax, beta=None, stack=True,REML=True, returnBeta=False)

                if ps == "nan":
                        joint_post = "nan"
                else:
                        # Convert obtained p-value to chisq then to association density functions
                        if ps == 1.0: ps = 0.99999999
                        assoc_lambda = (sts.norm.ppf(1-0.05/float(self.no_of_snps)/2.0) + sts.norm.ppf(0.8))**2
                        assoc_chisq = sts.chi2.ppf(1-np.array(ps),1)
                        if np.isinf(assoc_chisq): assoc_chisq=68.7632522116684
                        joint_pdf_lambda = sts.ncx2.pdf(assoc_chisq,1,nc=assoc_lambda) ; joint_pdf_0 = sts.chi2.pdf(assoc_chisq,1)

                        # Obtain Joint posterior probability of assocation with AM_POST as prior
                        joint_post = (joint_pdf_lambda*am_post)/float((joint_pdf_lambda*am_post)+(joint_pdf_0*(1-am_post)))
                        if joint_post == 1.0: joint_post = 0.99999999

                #Record the SNP-only pvalue and the Joint PPA.
                joint_fout.writelines(str(beta)+"\t"+str(ps)+"\t"+str(joint_post)+"\n")

        #Close the files opened for writting
        joint_fout.close()
        for idx2 in anc_prior_fout_dict:
                anc_prior_fout_dict[idx2].close()

class JasMAP(Joint_Association):
    '''Performs the overrall steps of JasMAP as described in the JasMAP method'''

    def jasmap(self):
        ''' Running JasMAP '''

        overall_starting_time = datetime.datetime.today()
        logger.info('Starting at time:%s' % str(overall_starting_time))
        self.res = run
        logger.info("Loading parameters from %s ..."%os.path.abspath(argv))
        logger.info("Options in effect:")

        for param in self.res.Params:
                logger.info('             '+str(self.res.Params[param][0])+': '+str(self.res.Params[param][1]))

        self.run_joint_analysis()

        overall_finishing_time = datetime.datetime.today()
        logger.info("\n\nFinish at time:%s"%str(overall_finishing_time))
        logger.info("Analysis took %s time to finish.\n\n\n"%str(overall_finishing_time-overall_starting_time))
        self.res.terminate()


if __name__ == '__main__':

    try:
        global gpv1
        argv = sys.argv[1]
    except IndexError:
        argv = ''
    finally:
        run = JasMAP(argv)
        run.jasmap()






