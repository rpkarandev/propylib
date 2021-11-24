
#### Junk for now
#### todo: organize and clean



import numpy as np
import pandas as pd
import propylib as plb
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
import os
import sys
from scipy import stats 
import warnings


path = '/lib/forcefields/'

## using adjusted boxplot to mark outliers
#from statsmodels.stats import stattools 
from robustats import medcouple


def getFromDict1(dic, x):
    if x in dic.keys():
        return dic[x]
    else:
        return np.nan
    
def getFromDict2(dic, x):
    if x in dic.keys():
        return dic[ x ]
    else:
        return ''
    
def get_mZscore(x0, dist='n', returnCutoff=False):
    x = x0
    #x = (x0 - np.nanmin(x0) )/ ( np.nanmax(x0) - np.nanmin(x0) ) * 100 ## does not affect
    med = np.nanmedian( x )
    ma = stats.median_abs_deviation( x, nan_policy='omit' )
    ma_med = np.nanmedian( ma )
    b = 0.67449 # need 0.75 quantile from underlying distribution
    ## estimating underlying distribution is not an easy task
    #b = np.quantile(x1/100, 0.75)
    
    z = ( b * (x - med) ) / ma_med
    print( 'm:{0:4.1f}, m(mad):{1:4.3f}, range: {2:4.3f}-{3:4.3f}, outlierfrac: {4:4.3f},{5:4.3f}'.format( med, ma_med, np.nanmin(z), np.nanmax(z), np.nanmean(z < -3.5),  np.nanmean(z > 3.5) ) )
    if not returnCutoff:
        return z
    else:
        return z, dict(zip( ['b', 'med', 'ma_med'], [b, med, ma_med] ))

def transform_and_get_mZscore(x0, method='yeojohnson', ax=None, title='Probplot '):
    #x = x0
    ## transforming into Normal distribution
    if method == 'yeojohnson':
        x, lmbda = stats.yeojohnson(x0)
    else:
        x0 = x0 - np.nanmin(x0) + 0.1 ## must be of positive numbers
        x, lmbda = stats.boxcox(x0)
        
    if ax:
        prob = stats.probplot( x, dist=stats.norm, plot=ax)
        ax.set_title(title + method)
        
    med = np.nanmedian( x )
    ma = stats.median_abs_deviation( x, nan_policy='omit' )
    ma_med = np.nanmedian( ma )
    b = 0.67449 # need 0.75 quantile from underlying distribution
    ## estimating underlying distribution is not an easy task
    #b = np.quantile(x1/100, 0.75)
    
    z = ( b * (x - med) ) / ma_med
    print( 'm:{0:4.1f}, m(mad):{1:4.3f}, range: {2:4.3f}-{3:4.3f}, outlierfrac: {4:4.3f},{5:4.3f}'.format( med, ma_med, np.nanmin(z), np.nanmax(z), np.nanmean(z < -3.5),  np.nanmean(z > 3.5) ) )
    return z, lmbda


def getOutliers_adjBox(x, d=1.5, returnCutoff=False):
    #if len(x) < 10000: ## compare and matches ; ;; ; they are close
    #    mc0 = stattools.medcouple( x )
    #    print( 'mc0:{}'.format(mc0) )
    mc = medcouple( x[~x.isnull()].values )
    Q1 = np.nanquantile(x, 0.25)
    Q3 = np.nanquantile(x, 0.75)
    Q2 = np.nanquantile(x, 0.5)
    IQR = Q3 - Q1
    Ql = Q1 - d * IQR * np.exp(-3 * mc)
    Qu = Q3 + d * IQR * np.exp(4 * mc)
    
    y = np.zeros( x.shape )
    y[ x < Ql ] = -1
    y[ x > Qu ] = 1

    Qlf = np.nanmean(y == -1 )
    Quf = np.nanmean(y == 1 )
    print( 'mc:{}, Q1:{}, Q2:{}, Q3:{}, IQR:{}, Ql:{}, Qu:{}, Qlf:{}, Quf:{}'.format(mc, Q1, Q2, Q3, IQR, Ql, Qu, Qlf, Quf ) )
    if not returnCutoff:
        return pd.Series( y, index=x.index)
    else:
        return pd.Series( y, index=x.index), dict(zip( ['mc', 'Q1', 'Q2', 'Q3', 'IQR', 'Ql', 'Qu', 'Qlf', 'Quf'], [mc, Q1, Q2, Q3, IQR, Ql, Qu, Qlf, Quf] ))
    
def getStats(x):
    d = {}
    d['count'] = len(x)
    d['mean'] = np.nanmean(x)
    d['median'] = np.nanmedian(x)
    d['median'] = np.nanmedian(x)
    d['max'] = np.nanmax(x)
    d['min'] = np.nanmin(x)
    d['std'] = np.nanstd(x)
    d['se'] = np.nanstd(x)/np.sqrt( len(x) )
    d['per025'] = np.nanpercentile(x, 2.5)
    d['per10'] = np.nanpercentile(x, 10)
    d['per90'] = np.nanpercentile(x, 90)
    d['per95'] = np.nanpercentile(x, 95)
    d['per975'] = np.nanpercentile(x, 97.5)
    return pd.Series(d)

def pvalue2code(x):
    if x > 0.05:
        return 0
    elif x > 0.01:
        return 1
    elif x > 0.001:
        return 2
    elif x > 0.00001:
        return 3
    elif x <= 0.00001:
        return 4
    else:
        -1
pvalue2code = np.vectorize( pvalue2code )

# Create models from data
def best_fit_distribution(data, bins=200, ax=None):
    """Model data by finding best fit distribution to data"""
    # Get histogram of original data
    y, x = np.histogram(data, bins=bins, density=True)
    x = (x + np.roll(x, -1))[:-1] / 2.0

    percentile_bins = np.linspace(0,100,11)
    percentile_cutoffs = np.percentile(data, percentile_bins)
    observed_frequency, bins = (np.histogram(data, bins=percentile_cutoffs))
    cum_observed_frequency = np.cumsum(observed_frequency)
    size = len(data)
    
    # Distributions to check
    DISTRIBUTIONS = [        
        stats.alpha,stats.anglit,stats.arcsine,stats.beta,stats.betaprime,stats.bradford,stats.burr,stats.cauchy,stats.chi,stats.chi2,stats.cosine,
        stats.dgamma,stats.dweibull,stats.erlang,stats.expon,stats.exponnorm,stats.exponweib,stats.exponpow,stats.f,stats.fatiguelife,stats.fisk,
        stats.foldcauchy,stats.foldnorm,stats.genlogistic,stats.genpareto,stats.gennorm,stats.genexpon,
        stats.genextreme,stats.gausshyper,stats.gamma,stats.gengamma,stats.genhalflogistic,stats.gilbrat,stats.gompertz,stats.gumbel_r,
        stats.gumbel_l,stats.halfcauchy,stats.halflogistic,stats.halfnorm,stats.halfgennorm,stats.hypsecant,stats.invgamma,stats.invgauss,
        stats.invweibull,stats.johnsonsb,stats.johnsonsu,stats.ksone,stats.kstwobign,stats.laplace,stats.levy,stats.levy_l,stats.levy_stable,
        stats.logistic,stats.loggamma,stats.loglaplace,stats.lognorm,stats.lomax,stats.maxwell,stats.mielke,stats.nakagami,stats.ncx2,stats.ncf,
        stats.nct,stats.norm,stats.pareto,stats.pearson3,stats.powerlaw,stats.powerlognorm,stats.powernorm,stats.rdist,stats.reciprocal,
        stats.rayleigh,stats.rice,stats.recipinvgauss,stats.semicircular,stats.t,stats.triang,stats.truncexpon,stats.truncnorm,stats.tukeylambda,
        stats.uniform,stats.vonmises,stats.vonmises_line,stats.wald,stats.weibull_min,stats.weibull_max,stats.wrapcauchy
    ]
    DISTRIBUTIONS = ['gamma', 'alpha', 'beta', 'rayleigh', 'norm', 'pareto', 'invgauss'
                    'rayleigh', 'weibull_min','weibull_max', 'rayleigh','powerlaw', 'dweibull', 
                     'pearson3','chi', 'chi2','f', 't', 'exponweib', 'laplace', 'invweibull',
                    'powerlognorm', 'powerlaw','powernorm', 'laplace','invgamma', 'logistic',
                     'genlogistic','gennorm', 'genexpon', 'exponpow','exponnorm', 'expon', 'maxwell',
                    'invweibull', 'loggamma', 'lognorm', 'rdist', 'vonmises', 'betaprime', 'erlang']
    # Best holders
    best_distribution = stats.norm
    best_params = (0.0, 1.0)
    best_sse = np.inf

    results_stats = []
    # Estimate distribution parameters from data
    for dist_name in DISTRIBUTIONS:
        #print(dist_name)
        # Try to fit the distribution
        try:
            # Ignore warnings from data that can't be fit
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')

                # fit dist to data
                distribution = getattr(stats, dist_name)
                params = distribution.fit(data)

                # Separate parts of parameters
                arg = params[:-2]
                loc = params[-2]
                scale = params[-1]

                # Calculate fitted PDF and error with fit in distribution
                pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)
                sse = np.sum(np.power(y - pdf, 2.0))

                # if axis pass in add to plot
                try:
                    if ax:
                        pd.Series(pdf, x).plot(ax=ax)
                    end
                except Exception:
                    pass

                # Get expected counts in percentile bins
                # cdf of fitted sistrinution across bins
                cdf_fitted = distribution.cdf(percentile_cutoffs, *params)
                expected_frequency = []
                for bin in range(len(percentile_bins)-1):
                    expected_cdf_area = cdf_fitted[bin+1] - cdf_fitted[bin]
                    expected_frequency.append(expected_cdf_area)

                # Chi-square Statistics
                expected_frequency = np.array(expected_frequency) * size
                cum_expected_frequency = np.cumsum(expected_frequency)
                ss = sum (((cum_expected_frequency - cum_observed_frequency) ** 2) / cum_observed_frequency)

                ## KS test
                #print(distribution, params)
                D, p = stats.kstest(data, dist_name, args=params)
                results_stats.append( [dist_name,  params, sse, ss, D, p] )

                # identify if this distribution is better
                if best_sse > sse > 0:
                    best_distribution = distribution
                    best_params = params
                    best_sse = sse

        except Exception:
            pass
    dfResult = pd.DataFrame( results_stats )
    dfResult.columns = ['dist', 'params', 'sse', 'chi2', 'KS_D', 'KS_pvalue']
    return (dfResult.sort_values(['chi2', 'KS_D']), best_distribution.name, best_params)

def make_pdf(dist, params, size=10000):
    """Generate distributions's Probability Distribution Function """

    # Separate parts of parameters
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Get sane start and end points of distribution
    start = distats.ppf(0.01, *arg, loc=loc, scale=scale) if arg else distats.ppf(0.01, loc=loc, scale=scale)
    end = distats.ppf(0.99, *arg, loc=loc, scale=scale) if arg else distats.ppf(0.99, loc=loc, scale=scale)

    # Build PDF and turn into pandas Series
    x = np.linspace(start, end, size)
    y = distats.pdf(x, loc=loc, scale=scale, *arg)
    pdf = pd.Series(y, x)
    
    return pdf

######################################################################################################

###############################
atom_nomen = pd.read_csv(path + 'atom_nom_mod.tbl', sep='\t', comment='#', header=None)
atom_nomen.columns = 'aa1		BMRB	SC	PDB	UCSF	MSI	XPLOR	SYBYL*	MIDAS*	DIANA'.split('\t')
atom_nomen['aa3'] = atom_nomen.aa1.apply(lambda x: aaA2AAA[x])
###############################
amass = pd.read_csv(path + 'amber/parm14ipq_mass.dat', keep_default_na=False, engine='python', sep='\s{2,}', header=None)
amass.columns = ['atom', 'amass', 'atpol', 'comment']
amass.set_index('atom', drop=False, inplace=True)
'''
KNDSYM     The unique atom symbol used in the system.

 AMASS      Atomic mass of the center having the symbol "KNDSYM".

 ATPOL      The atomic polarizability for each atom (in A**3)
            This is the type of polarizability used in sander
            and gibbs. No parameters are supplied for this since
            the feature is still in development (Amber 4.1).
'''
###############################
aradius = pd.read_csv(path + 'amber/parm14ipq_radius.dat', keep_default_na=False, engine='python', sep='\s{2,}', header=None)
aradius.columns = ['atom', 'vradius', 'potwell', 'comment']
aradius.set_index('atom', drop=False, inplace=True)
'''
'RE'       van der Waals radius and the potential well depth
                     parameters are read.
'''
###############################
libfiles = glob(path + 'amber/dat/leap/lib/amino*14ipq.lib')
lines0 = []
for libfile in libfiles:
    with open(libfile) as Famino:
        print('Reading amber library file {}'.format(libfile))
        lines0 += Famino.readlines()

if len(lines0) > 0:
    lines = list(lines0)
    aalist = [] 
    i = 0
    aa30 = ''
    print( 'no of lines:', len(lines) )
    line = re.split(' +', lines.pop(0).strip())
    tabledict = {}
    arraydicts = {}
    valuedicts = {}
    while lines and line[0].startswith('!'):
        title = line[0].split('.')
        if title[0] == '!entry':
            #print(title)
            datatitle = title[3]
            aa3 = title[1]
            if aa30 != aa3:
                if aa30 != '':
                    arraydicts[aa30] = arraydict
                    valuedicts[aa30] = valuedict
                aa30 = aa3
                arraydict = {}
                valuedict = {}
               
            if line[1] == 'table':
                cols = [ line[ci] for ci in range(3, len(line), 2 ) ]
                temparray = []
                line = re.split(' +', lines.pop(0).strip())
                i = i + 1
                while lines and not line[0].startswith('!'):
                    temparray.append( line ) 
                    line = re.split(' +', lines.pop(0).strip())
                    i = i + 1                
                tempdf = pd.DataFrame(temparray, columns=cols)
                tempdf.reset_index(drop=False, inplace=True)
                tempdf['aa3'] = aa3
                if datatitle in tabledict.keys():
                    tabledict[datatitle].append(tempdf)
                else:
                    tabledict[datatitle] = [tempdf]
            elif line[1] == 'array':
                temparray = []
                line = re.split(' +', lines.pop(0).strip())
                i = i + 1
                while lines and not line[0].startswith('!'):
                    temparray.append( line[0] ) 
                    line = re.split(' +', lines.pop(0).strip())
                    i = i + 1                
                arraydict[ datatitle ] = temparray
            elif line[1] == 'single':
                line = re.split(' +', lines.pop(0).strip())
                i = i + 1
                while lines and not line[0].startswith('!'):
                    valuedict[ datatitle ] = line[0]   
                    line = re.split(' +', lines.pop(0).strip())
                    i = i + 1                
            else:
                line = re.split(' +', lines.pop(0).strip())
                i = i + 1
        elif title[0] == '!!index':
            line = re.split(' +', lines.pop(0).strip())
            i = i + 1
            while lines and not line[0].startswith('!'):
                aalist.append(line)
                line = re.split(' +', lines.pop(0).strip())
                i = i + 1
        else:
            line = re.split(' +', lines.pop(0).strip())
            i = i + 1
    singlevalues = pd.DataFrame.from_dict(valuedicts).T
    arrayvalues = pd.DataFrame.from_dict(arraydicts).T
    tables = {x:pd.concat( tabledict[x] ) for x in tabledict.keys() }
    atoms = tables['atoms']
    atoms['name'] = atoms['name'].str.replace('"', '')
    atoms['type'] = atoms['type'].str.replace('"', '')
    atoms['mass'] = atoms['type'].apply(lambda x: amass['amass'][x])
    atoms['atpol'] = atoms['type'].apply(lambda x: amass['atpol'][x])
    atoms['vradius'] = atoms['type'].apply(lambda x: aradius['vradius'][x])
    atoms['potwell'] = atoms['type'].apply(lambda x: aradius['potwell'][x])
atoms['chg'] = atoms.chg.astype('float')   
atoms['aa30'] = atoms['aa3'].apply(lambda x: x[-3:])
atoms['aa1'] = atoms.aa30.apply(lambda x: aaAAA2A[x])
 
###########################
tables['positions']['x'] = tables['positions']['x'].astype(float)
tables['positions']['y'] = tables['positions']['y'].astype(float)
tables['positions']['z'] = tables['positions']['z'].astype(float)

atomsi = atoms.set_index(['aa3', 'index'], drop=False)
atom_positions = tables['positions'].set_index(['aa3', 'index'], drop=False)
for col in ['x', 'y', 'z']:
    atoms[col] =  atoms.apply(lambda x: atom_positions.loc[ (x['aa3'], x['index']), col], axis=1)

###########################

###########################
atom_nomen['amber'] = atom_nomen.iloc[:,2]
#atom_nomen.loc[atom_nomen.aa1 =='X' ,'amber'] = []
x = atoms[atoms.aa30 == 'TYR']
(atoms.name.isin(atom_nomen.iloc[:,2])).mean()
x = atoms[~atoms.name.isin(atom_nomen.iloc[:,2]) ]
{x: (atoms[atoms.aa1==x].name.isin(atom_nomen[atom_nomen.aa1.isin([x,'X'])].iloc[:,2])).mean() for x in aa1code}
xa = 'N'
x = atoms[~atoms.name.isin(atom_nomen[atom_nomen.aa1.isin([xa,'X'])].iloc[:,2]) ]
x[x.aa1=='N']

###########################
atoms.to_csv('./atom_details_amber.csv', index=False)





#######################################################################################################
'Sequence Mutation - GAP - others'
'Single mutation of sequence '
#Functions: mutate, tangorun  ..  creates single fasta file with all mutants
def mutate(infile, start = 1, end = 0, mutAA = 'A'):
	#Input Parameters
	#start = 1 # Default 1
	#end = 0   # 0 not defined: auto set to max (total length)
	#mutAA = 'A'

	#Verifying the argument
	IFILE = infile
	temp = IFILE.split('/')
	temp = temp[len(temp)-1].split('.')
	if temp[1].lower() != 'fasta' :
		print("try executing the script again, with just fasta file !!")
		exit()

	print('Starting script ... provided "' + temp[0] + "." + temp[1] + '" as input fasta file')
	pro = temp[0] + "_"

	#Reading Input file 
	F = open(IFILE, 'r')
	print("Reading file . . .. ")
	content = F.read().splitlines()

	i = 1
	seq = ""
	title =  ""
	#Extracting Sequence and Protein Title from input Fasta sequence
	for line in content:			
		if i == 1 and line[0] == '>':
			title = line
		elif i != 1 and line[0] == '>':
			print("Only one fasta sequence is expected. Plz check the file content for multiple sequence or blank lines !!")
			break
		else:
			seq = seq + line
			print(i-1, "  ", line)
		i = i + 1
		F.close()
	print("Input sequence file: ", IFILE, " read successfully.\nThe sequence titled as ", title)
	print("The sequence (", len(seq), "AA ) is ", seq)

	i = 1	
	#Creating folder and Defining Prefix for output file names
	temp = pro + 'mutants_' + mutAA + ".fasta" 
	while os.path.exists(temp):
		temp = pro + 'mutants_'+  mutAA + str(i) + ".fasta" 
		i = i + 1
	OFILE = temp

	mutno = 0
	#Creating Mutant Sequences
	i = start-1;
	if end == 0:
		end = len(seq)


	tmp = "X"
	F = open(OFILE, 'w')
	F.write(">mutant_" + str(mutno) + "_res" + str(0) + tmp + "\n")
	F.write(seq+ "\n")	
	
	while i < end:
#		print i;	
		#print "started"
		temp = list(seq)
		tmp = 'X'
		if (temp[i] != mutAA) & (i > -1):
		#Mutating the original sequence	using Alanine at position i			
			tmp = temp[i]
			temp[i] = mutAA
			mutseq = "".join(temp)
			mutno = mutno + 1
		#Writing mutated sequence to file
			#temp = title + " Mutant_" + str(i+1) + "\n"
			F.write(">mutant_" + str(mutno) + "_res" + str(i+1) + tmp + "\n")
			F.write(mutseq+ "\n")
		i = i + 1 
	F.close()
	print("job completed successfully !!", mutno,"Mutant sequences have been created in file: ", OFILE, 'by replacing residues with amino acid "'+mutAA+'"'			)
	return OFILE	

def tangorun(infile, ct="N",nt="N",ph=7.4,te=303,io=0.1):
	#Input Parameters
	curdir = os.getcwd();
	epath = 'C:/Users/rpkaran/Documents/GoogleDrive/Project/mylib/mutation_executables/Tango.exe'
	#ct="N"		# Protection at the C-terminus: can be N for no or Y for amidated
	#nt="N"  	# Protection at the N-terminus: can be N for no, A for acetylated or S for succinilated
	#ph="7.4" 	# pH
	#te="303" 	# Temperature in Kelvin
	#io="1.05" 	# Ionic strength M
	seq = ""
	param = ct + " " + nt + " " + str(ph)  + " " + str(te)  + " " + str(io)
#	param = "N N 7 298 0.1"
	app = "tango"
	titlelist = list()

	#Verifying the argument 
	
	IFILE = infile
	temp = IFILE.split('/')
	temp = temp[len(temp)-1].split('.')
	pro = temp[0]
	if temp[1].lower() != 'fasta' :
		print("try executing the script again, with just fasta file !!")
		exit()

	#Reading Input File and Creating output file (Tango input file)
	print('Starting script ... provided "' + temp[0] + "." + temp[1] + '" as input fasta file')
	OFILE = "mutants_inputs" + ".txt"

	IN = open(IFILE, 'r')
	print("Reading file . . .. ")
	content = IN.read().splitlines()


	print("Creating output directory ...")
	#Creating folder and Defining Prefix for output file names
	i = 1
	temp = app + "_" + pro
	while os.path.exists(temp):
		temp = app  + str(i) + pro + "_" 
		i = i + 1
	os.mkdir(temp)
	os.chdir(temp)
	print("Created output folder: ", temp)
	OUT = open(OFILE, 'w')

	i = 1
	seq = ""
	seqno = 0
	title =  ""
	#Extracting Sequence and Protein Title from input Fasta sequence
	for line in content:	
		if i == 1 and line[0] == '>':
			title = line[1:len(line)]
		elif i != 1 and line[0] == '>':
			OUT.write(title + " " + param + " ")
			OUT.write(seq + "\n")
			seqno = seqno + 1
			seq = ""
			titlelist.append(title)
			#print "Processed sequence #",seqno
			title = line[1:len(line)]		
		else:
			seq = seq + line		
		i = i + 1
		IN.close()
	OUT.write(title + " " + param + " ")
	OUT.write(seq + "\n")
	seqno = seqno + 1
	titlelist.append(title)
	#print "Processed sequence #",seqno,"\n"
	OUT.close()

	print("Input sequence file: ", IFILE, " read successfully.",seqno,"sequences were found !!")
	print("Input file for",app,"created successfully: ", OFILE, "In folder ", temp)

	cmd = epath + " -inputfile=" + OFILE + " > tango.log"
	print("executing ",cmd)
	status = subprocess.call(cmd, shell=True)
	os.chdir(curdir)
	if status == 0: 
		print(app,"successfully executed !!")
		#os.chdir(curdir)  #os.chdir('../')
		#print curdir 	print os.getcwd()
	else:
		print("Error in execution of", app)
	return 0	

def GAPrun(infile):
	#Input Parameters
	#epath = '/home/user1/Documents/Project/mylib/mutation_executables/tango_x86_64_release'
	#ct="N"		# Protection at the C-terminus: can be N for no or Y for amidated
	#nt="N"  	# Protection at the N-terminus: can be N for no, A for acetylated or S for succinilated
	#ph="7.4" 	# pH
	#te="303" 	# Temperature in Kelvin
	#io="1.05" 	# Ionic strength M
	seq = ""
#	param = ct + " " + nt + " " + str(ph)  + " " + str(te)  + " " + str(io)
	param = ""
#	param = "N N 7 298 0.1"
	app = "GAP"
	titlelist = list()

	#Verifying the argument 
	
	IFILE = infile
	temp = IFILE.split('/')
	temp = temp[len(temp)-1].split('.')
	pro = temp[0]
	if temp[1].lower() != 'fasta' :
		print("try executing the script again, with just fasta file !!")
		exit()

	#Reading Input File ##and Creating output file (Tango input file)
	print('Starting script ... provided "' + temp[0] + "." + temp[1] + '" as input fasta file')
	#OFILE = "mutants_inputs" + ".txt"

	IN = open(IFILE, 'r')
	print("Reading file . . .. ")
	content = IN.read().splitlines()


	print("Creating output directory ...")
	#Creating folder and Defining Prefix for output file names
	i = 1
	temp = app + "_" + pro
	while os.path.exists(temp):
		temp = app  + str(i) + pro + "_" 
		i = i + 1
	os.mkdir(temp)
	os.chdir(temp)
	print("Created output folder: ", temp)
	#OUT = open(OFILE, 'w')

	i = 1
	seq = ""
	seqno = 0
	title =  ""
	#Extracting Sequence and Protein Title from input Fasta sequence
	for line in content:	
		if i == 1 and line[0] == '>':
			title = line[1:len(line)]
		elif i != 1 and line[0] == '>':
			#OUT.write(title + " " + param + " ")
			#OUT.write(seq + "\n")
			seqno = seqno + 1
			seq = ""
			titlelist.append(title)
			#print "Processed sequence #",seqno
			title = line[1:len(line)]		
		else:
			seq = seq + line		
		i = i + 1
		IN.close()
	#OUT.write(title + " " + param + " ")
	#OUT.write(seq + "\n")
	seqno = seqno + 1
	titlelist.append(title)
	#print "Processed sequence #",seqno,"\n"
	#OUT.close()

	print("Input sequence file: ", IFILE, " read successfully.",seqno,"sequences were found !!")
	#print "Input file for",app,"created successfully: ", OFILE, "In folder ", temp

	#cmd = epath + " -inputfile=" + OFILE + " > GAP.log"
#	print "executing ",cmd
	#status = subprocess.call(cmd, shell=True)
	status = accessGAP(infile, titlelist)
	if status == 0:
		print(app,"successfully executed !!")
		os.chdir('../')
		return titlelist
	else:
		print("Error in execution of", app)
		return -1
		
		
def accessGAP(infile, titlelist):
	outfile1 = "GAPout.html"
	outfile2 = "GAPout.txt"
	#os.mkdir("test")
	F = open("../"+infile, 'r')
	Seq = F.read()
	F.close()

	#Seq = '>sp|P37840|SYUA_HUMAN Alpha-synuclein OS=Homo sapiens GN=SNCA PE=1 SV=1' + '\n' + 'MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'
	params = urllib.urlencode({'sequence': Seq, 'predict': 'Predict'})
	#, '@type': 'submit'
	headers = {"Content-type": "application/x-www-form-urlencoded", "Accept": "text/html","Accept-Encoding":"gzip, deflate", "Connection":"keep-alive","Host":"www.iitm.ac.in","Origin":"http://www.iitm.ac.in","Referer":"http://www.iitm.ac.in/bioinfo/GAP/"}
	conn = httplib.HTTPConnection("10.93.249.92:80")
	#"Accept":"text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8"

	#print(Seq)
	#print(params)
	#print(headers)

	conn.request("POST", 'http://www.iitm.ac.in/bioinfo/cgi-bin/GAP/cal_fasta.pl', params, headers)
	response = conn.getresponse()
	print(response.status, response.reason)
	#if response.headers().get('Content-Encoding') == 'gzip':
	buf = StringIO( response.read())
	f = gzip.GzipFile(fileobj=buf)
	conn.close()

	html = f.read()
	#html = html.decode('utf8')

	#html = open("foobar.html").read()
	#data = html2text.html2text(html)
	#print(data)
	#soup = BeautifulSoup(html)

	F = open(outfile1, 'w')
	F.write(html)
	F.close()

#	F = open(outfile2, 'w')
#	F.write(data)
#	F.close()
	
	html2dat(outfile1, titlelist)
	return 0
	
def html2dat(infile, titlelist):
	
	#infile = "out.html"
	outfile1 = "outtable1.txt"
	outfile2 = "outtable2.txt"
	outfile3 = "outtable3.txt"

	PATnumSeq = re.compile('(?<=Number of sequences: )\d*',re.IGNORECASE)
	PATSeq = re.compile('(?<=Input sequence</u>: )\w*',re.IGNORECASE)
	PATout = re.compile('(?<=<br>The given peptide is Amyloid \()\d*.\d*',re.IGNORECASE)

	PATtable = re.compile('<table.*?</table>',re.IGNORECASE)
	PATtr = re.compile('<tr>.*?</tr>',re.IGNORECASE)
	PATtd = re.compile('<td.*?</td>',re.IGNORECASE)
	PATtag = re.compile('<.*?>',re.IGNORECASE)

	theader = ("No","Peptide","Amyloid probability","Non-amyloid probability")
	gap = 6 - 1
	txt = ""
	ttxt = "\t" + " No " + "\t" + "Peptide" + "\t" + "Amyloid probability" + "\t" + "Non-amyloid probability" + '\t'
	t = list()
	data = list()

	F = open(infile, 'r')
	Seq = F.read()
	F.close()

	NumSeqs = int(PATnumSeq.findall(Seq)[0])
	Seqs = PATSeq.findall(Seq)
	outs = PATout.findall(Seq)

	print(NumSeqs)
	#print(Seqs)
	#print(outs)

	for i in range(0,NumSeqs):
		#print i
		txt = txt + titlelist[i] + '\t' + str(i) + '\t' + outs[i] + '\t' + Seqs[i] + '\n'
		#print txt
	
	with open(outfile1, "w") as output:
		output.write(txt)

	txt = ""
	seltxt = ""
	
	#tables = pat1.match(Seq)
	tables = PATtable.findall(Seq)
	if tables:
	#	print tables
	#	print tables.groups()
		i = -1
		for table in tables:		
			i = i + 1
			mutation = titlelist[i]
			#mutant_0_res0X
			sel = int(re.findall('(?<=res)\d*',mutation)[0])
	#		print i,table,'\n'
			trows = PATtr.findall(table)
			if trows:
				j = -1
				for trow in trows:
					j = j + 1
	#				print i,j,PATtag.sub("\t",trow) ,i
					temp = PATtag.sub("\t",trow) 
					temp = re.sub("\t\t","\t",temp)
					if (temp != ttxt):
						txt = txt + mutation + temp + '\n'
						if (j in range(sel - gap,sel+1)) or (sel == 0):
							seltxt = seltxt + mutation + temp + '\n'				
			
	if (NumSeqs != i+1):
		print("Warn: Out of ",NumSeq,", ",i+1," outputs were retrieved !!!. #########")

	#print(data)
	#print(type(txt))		
	#F.write(str(data))	
	#with open(outfile2, "w") as output:
	#	writer = csv.writer(output, lineterminator='\n')
	#	writer.writerows(data)			
	with open(outfile2, "w") as output:
		output.write(txt)
	with open(outfile3, "w") as output:
		output.write(seltxt)
	return(0)	
	
def PASTA2run(infile):
	perlf="PastaPairs.pl"
	pot="pot_pasta.dat"
	abspath = os.path.dirname(os.path.abspath(__file__)) + "/PASTA2II/"
	path= abspath 
	top=10
	flag= 2
	typ = "self"
	cutoff = -5.5
	base = ""
	#Input Parameters
	epath = '/home/user1/Documents/Project/mylib/pasta_exe/'
	#ct="N"		# Protection at the C-terminus: can be N for no or Y for amidated
	#nt="N"  	# Protection at the N-terminus: can be N for no, A for acetylated or S for succinilated
	#ph="7.4" 	# pH
	#te="303" 	# Temperature in Kelvin
	#io="1.05" 	# Ionic strength M
	app = "PASTA2"
	titlelist = list()

	#Verifying the argument 
	
	IFILE = infile
	temp = IFILE.split('/')
	temp = temp[len(temp)-1].split('.')
	pro = temp[0]
	if temp[1].lower() != 'fasta' :
		print("try executing the script again, with just fasta file !!")
		exit()

	#Reading Input File and Creating output file (Tango input file)
	print('Starting script ... provided "' + temp[0] + "." + temp[1] + '" as input fasta file')
	OFILE = "mutants_inputs" + ".txt"

	IN = open(IFILE, 'r')
	print("Reading file . . .. ")
	content = IN.read().splitlines()


	print("Creating output directory ...")
	#Creating folder and Defining Prefix for output file names
	i = 1
	temp = app + "_" + pro
	while os.path.exists(temp):
		temp = app  + str(i) + pro + "_" 
		i = i + 1
	os.mkdir(temp)
	os.chdir(temp)
	path = path + temp + '/'
	print("Created output folder: ", temp)
	

	i = 1
	seq = ""
	seqno = 0
	title =  ""
	#Extracting Sequence and Protein Title from input Fasta sequence
	status = 0
	for line in content:	
		if i == 1 and line[0] == '>':
			title = line[1:len(line)]
		elif i != 1 and line[0] == '>':			
			base = title
			os.mkdir(base)
			os.chdir(base)
			OFILE = "mutants_inputs" + str(seqno) + ".fasta"
			OUT = open(OFILE, 'w')			
			OUT.write(">" + title + "\n" + seq )
			OUT.close()			
			cmd = "(cd " + epath + " && perl " + perlf + " " + pot + " " + path + base + '/ ' + str(top) + " " + str(flag) + " " + str(typ) + " " + str(cutoff) + " " + base + ")" + ' > PASTA2.log'
			#print cmd
			print("Executing ", app,": mutant seqno -> ", seqno)
			status = subprocess.call(cmd, shell=True)
			if status != 0:				
				print("Error in execution of", app, ": sequenceno",seqno)
				break
			os.chdir(path)
			seq = ""
			seqno = seqno + 1
			titlelist.append(title)
			#print "Processed sequence #",seqno
			title = line[1:len(line)]		
		else:
			seq = seq + line		
			i = i + 1
		IN.close()
	if status == 0:		
		base = title
		os.mkdir(base)
		os.chdir(base)
		OFILE = "mutants_inputs" + str(seqno) + ".fasta"
		OUT = open(OFILE, 'w')
		OUT.write(">" + title + "\n" + seq )
		OUT.close()
		cmd = "(cd " + epath + " && perl " + perlf + " " + pot + " " + path + base + '/ ' + str(top) + " " + str(flag) + " " + str(typ) + " " + str(cutoff) + " " + base + ")" + ' > PASTA2.log'
		status = subprocess.call(cmd, shell=True)
		if status != 0:				
			print("Error in execution of", app, ": sequenceno",seqno)
		seqno = seqno + 1
		titlelist.append(title)	
	
	print("Input sequence file: ", IFILE, " read successfully.",seqno,"sequences were found !!")
	print("Input file for",app,"created successfully in folder ", temp)
	os.chdir(abspath)
	return status	


def waltzrun(infile):
	#Input Parameters
	#epath = '/home/user1/Documents/Project/mylib/mutation_executables/tango_x86_64_release'
	#ct="N"		# Protection at the C-terminus: can be N for no or Y for amidated
	#nt="N"  	# Protection at the N-terminus: can be N for no, A for acetylated or S for succinilated
	#ph="7.4" 	# pH
	#te="303" 	# Temperature in Kelvin
	#io="1.05" 	# Ionic strength M
	seq = ""
#	param = ct + " " + nt + " " + str(ph)  + " " + str(te)  + " " + str(io)
	param = ""
#	param = "N N 7 298 0.1"
	app = "waltz"
	titlelist = list()

	#Verifying the argument 
	
	IFILE = infile
	temp = IFILE.split('/')
	temp = temp[len(temp)-1].split('.')
	pro = temp[0]
	if temp[1].lower() != 'fasta' :
		print("try executing the script again, with just fasta file !!")
		exit()

	#Reading Input File ##and Creating output file (Tango input file)
	print('Starting script ... provided "' + temp[0] + "." + temp[1] + '" as input fasta file')
	#OFILE = "mutants_inputs" + ".txt"

	IN = open(IFILE, 'r')
	print("Reading file . . .. ")
	content = IN.read().splitlines()


	print("Creating output directory ...")
	#Creating folder and Defining Prefix for output file names
	i = 1
	temp = app + "_" + pro
	while os.path.exists(temp):
		temp = app  + str(i) + pro + "_" 
		i = i + 1
	os.mkdir(temp)
	os.chdir(temp)
	print("Created output folder: ", temp)
	#OUT = open(OFILE, 'w')

	i = 1
	seq = ""
	seqno = 0
	title =  ""
	#Extracting Sequence and Protein Title from input Fasta sequence
	for line in content:	
		if i == 1 and line[0] == '>':
			title = line[1:len(line)]
		elif i != 1 and line[0] == '>':
			#OUT.write(title + " " + param + " ")
			#OUT.write(seq + "\n")
			seqno = seqno + 1
			seq = ""
			titlelist.append(title)
			#print "Processed sequence #",seqno
			title = line[1:len(line)]		
		else:
			seq = seq + line		
		i = i + 1
		IN.close()
	#OUT.write(title + " " + param + " ")
	#OUT.write(seq + "\n")
	seqno = seqno + 1
	titlelist.append(title)
	#print "Processed sequence #",seqno,"\n"
	#OUT.close()

	print("Input sequence file: ", IFILE, " read successfully.",seqno,"sequences were found !!")
	#print "Input file for",app,"created successfully: ", OFILE, "In folder ", temp

	#cmd = epath + " -inputfile=" + OFILE + " > GAP.log"
#	print "executing ",cmd
	#status = subprocess.call(cmd, shell=True)
	status = accesswaltz(infile, titlelist)
	if status == 0:
		print(app,"successfully executed !!")
		os.chdir('../')
		return titlelist
	else:
		print("Error in execution of", app)
		return -1
		
		
def accesswaltz(infile, titlelist):
	outfile1 = "waltzout.html"
	outfile2 = "waltzout.txt"
	#os.mkdir("test")
	F = open("../"+infile, 'r')
	Seq = F.read()
	F.close()

	#Seq = '>sp|P37840|SYUA_HUMAN Alpha-synuclein OS=Homo sapiens GN=SNCA PE=1 SV=1' + '\n' + 'MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'
	params = urllib.urlencode({'sequence': Seq, 'threshold': '92', 'custom_value' : 0, 'ph': 7.0, 'output': 'text_short', 'Submit': 'Submit sequences'})
	#, '@type': 'submit'
	headers = {"Content-type": "application/x-www-form-urlencoded", "Accept": "text/html","Accept-Encoding":"gzip, deflate", "Connection":"keep-alive","Host":"waltz.switchlab.org","Origin":"http://waltz.switchlab.org","Referer":"http://waltz.switchlab.org/results.cgi"}
	conn = httplib.HTTPConnection("134.58.50.6:80")
	#"Accept":"text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8"

	#print(Seq)
	#print(params)
	#print(headers)

	conn.request("POST", 'http://waltz.switchlab.org/results.cgi', params, headers)
	response = conn.getresponse()
	print(response.status, response.reason)
	#if response.headers().get('Content-Encoding') == 'gzip':
	buf = StringIO( response.read())
	f = gzip.GzipFile(fileobj=buf)
	conn.close()

	html = f.read()
	#html = html.decode('utf8')
	PATURL = re.compile('<a href="\.\.(.*?)">here</a>',re.IGNORECASE)
	URL = PATURL.findall(html)[0]		

	URL = 'http://waltz.switchlab.org/output' + URL
	print(URL)
	conn = httplib.HTTPConnection("134.58.50.6:80")
	conn.request("GET", URL)
	response = conn.getresponse()	
	print(response.status, response.reason)
#	if response.headers().get('Content-Encoding') == 'gzip':
#	buf = StringIO( response.read())
#	f = gzip.GzipFile(fileobj=buf)
#	html = f.read()
#	else:
	html = response.read()	
	conn.close()	

	
	F = open(outfile1, 'w')
	F.write(html)
	F.close()
	
	
	html2datwaltz(outfile1, titlelist)
	return 0
	
def html2datwaltz(infile, titlelist):	
	#infile = "out.html"
	outfile1 = "outtable1.txt"
	outfile2 = "outtable2.txt"
	outfile3 = "outtable3.txt"
	
	PATdiv = re.compile('<div id="bodyText">.*?</div>', re.DOTALL)
	PATtable = re.compile('<table.*?</table>',re.DOTALL)
	PATtr = re.compile('<tr>.*?</tr>', re.DOTALL)
	PATtd = re.compile('<td.*?</td>',re.IGNORECASE)
	PATtag = re.compile('<.*?>',re.IGNORECASE)
		
	txt = ""
	F = open(infile, 'r')
	Seq = F.read()
	F.close()
	theader = ("Job","Regions","Description")
	#print(Seq) 
	#print("######## seq ")
	div = PATdiv.findall(Seq)[0]
	#print div
	if div:
		table = PATtable.findall(div)[0]
	#print table				
	if table:
		trows = PATtr.findall(table)
	if trows:
		j = -1
		for trow in trows:
			j = j + 1			
			tcells = PATtd.findall(trow)
			if tcells:
				i = 0
				for tcell in tcells:				
					tcell = PATtag.sub("",tcell)					
					txt = txt + tcell + "\t"
					
				txt = txt + "\n"
							
	print("######", j ," outputs were retrieved !!!. #########")

	#print(data)
	#print(type(txt))		
	#F.write(str(data))	
	#with open(outfile2, "w") as output:
	#	writer = csv.writer(output, lineterminator='\n')
	#	writer.writerows(data)			
	with open(outfile1, "w") as output:
		output.write(txt)
	return(0)

####################################################################
def read_amylpred(filename):
    with open(filename, 'r') as ifile:
        tdict = {}
        line = ifile.readline()
        while line:
            if line[:16] == 'Sequence length:' :
                l = int(line.replace(' ', '').split(':')[1])
            if line[0:4] == 'HITS':
                line = ifile.readline()
                line = ifile.readline()
                temp = line.strip().split(':')
                while len(temp) == 2:                
                    #print(temp)
                    tdict[temp[0]] = temp[1]
                    line = ifile.readline()   
                    temp = line.strip().strip('>--> ').replace(' ', '').split(':')
            line = ifile.readline()
    tdict['seqlen'] = l
    tdict['filename'] = filename
    return tdict

def read_amylpred2(ifilename):
    pattern = re.compile("(\w+\s*\w*): ([0-9-,:# ]+)")
    flag = False
    dictlist = []
    dict1 = {}
    with open(ifilename, 'r') as ifile:
        line = ifile.readline()
        while line :
            match = pattern.findall(line)
            if len(match) > 0:
                #print( match )
                match = match[0]
                if match[0] == "Job ID":
                    if flag:
                        dictlist.append(dict1)
                        dict1 = {}
                    else:
                        flag = True
               # elif match[0] == "Sequence length":
                #    print(match)
                dict1[match[0]] = match[1]                
            line = ifile.readline()
        dictlist.append(dict1)
    return dictlist
