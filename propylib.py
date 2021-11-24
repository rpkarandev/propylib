#!/user/bin/python
## load core library
import glob
import math 
import numpy as np
import pandas as pd

### todo
## why not replace math with np: todo
## organize and split the file into smaller ones
## remove unwanted package

## load auxillary
#import csv
import gzip
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#from mpl_toolkits.mplot3d import axes3d
import os
import re
#import scipy
#import scipy.stats
#from sklearn.model_selection import KFold
#from sklearn.cluster import KMeans
import sys
import subprocess
#import httplib, urllib
#from StringIO import StringIO

#####################################################################################
### identify continuous stretches from a list of numbers 
### join stretches if they closer than "win" length
### input is a list
def find_uninterruptedseq(resnums0, win=0):
    resnums = resnums0.copy()
    #resnums = list(set(resnums))
    nlists = []
    while len(resnums) > 0: 
        s1 = resnums.pop(0)
        if (len(resnums) > 0) and resnums[0] == s1 + 1:           
            nlist = [s1, resnums[0]]
            s1 = resnums.pop(0)
            while (len(resnums) > 0) and (resnums[0] == s1 + 1):
                s1 = resnums.pop(0)
                nlist.append(s1)         
            if len(nlist) > 0:
                nlists.append(nlist)
        else:
            nlists.append([s1])
    if win >  0 :
        nlists1 = continuous(nlists, win)
        return nlists1
    else:
        return nlists 

### join segments/stretches
### seglist is a list of list
def continuous(seglist, win):
    if len(seglist) > 1 :
        seglistn = []
        for i in range(1,len(seglist)):
            if min(seglist[i]) - max(seglist[i-1]) < win:
                seglist[i] = seglist[i-1] + np.arange(max(seglist[i-1])+1, min(seglist[i])).tolist() + seglist[i]
            else:
                seglistn.append(seglist[i-1])
        seglistn.append(seglist[i])
        return seglistn
    else:
        return seglist
 
### covert boolean list to index list
def bool2ind(predicted, b=True):
    return np.where( np.array(predicted) == b)[0].tolist()    

### returns index of A which are present in B
def ismember(A, B):
    AinB = [x for x in A if x in B]
    AinBi = [i for i, x in enumerate(A) if x in AinB]
	return AinBi


############################################################################################
#### 'common amino acid residue values '
AAlist = list("GAVLIMFWPSTCYNQDEKRH")
hphoAAlist = list("GAVLIMFWP") #[ "G", "A", "V", "L", "I", "M", "F", "W", "P" ]
hphiAAlist = list("STCYNQDEKRH")

aa1code = ['A', 'D', 'C', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
aa3code = ['ALA', 'ASP', 'CYS', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
#aa1code = ['A', 'D', 'C', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
aapair = np.array([ [y+x for x in aa1code] for y in aa1code]).flatten()

aaA2AAA = {'A' :'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER', 'T':'THR', 'V': 'VAL', 'W' : 'TRP', 'Y' : 'TYR', 'X': 'UNK' } #, 'B': 'MSE', 'J' : 'PCA', 'O': 'SL5', 'J': 'SWG'
aaAAA2A = {'ALA':'A', 'CYS':'C', 'ASP' : 'D',  'GLU' : 'E' , 'PHE' : 'F', 'GLY' : 'G', 'HIS' : 'H', 'ILE' : 'I', 'LYS': 'K' , 'LEU':'L', 'MET':'M' , 'ASN':'N', 'PRO':'P' , 'GLN':'Q' , 'ARG':'R' , 'SER':'S' , 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y', 'UNK' :'X'} #, 'MSE':'B', 'PCA' : 'B', 'SL5' : 'B', 'SWG' : 'B' , 'TPO' : 'B', 'MIR' : 'B', 'PTR' : 'B', 'PIA' : 'B', 'CRF' : 'B', 'CZZ' : 'B'  

aa2chph74 = {'A' :0, 'C':0, 'D':-1, 'E':-1, 'F':0, 'G':0, 'H':0.5, 'I':0, 'K':1, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':1, 'S':0, 'T':0, 'V':0, 'W' : 0, 'Y' :0, 'X':0} 
aa2chph26 = {'A' :0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':1, 'I':0, 'K':1, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':1, 'S':0, 'T':0, 'V':0, 'W' : 0, 'Y' :0, 'X':0} 

## simplify DSSP/STRIDE secstr code
dssp2ss1 = {'H':'H', 'G':'H', 'I':'H', 'E':'E', 'B':'E', 'S':'T', 'T':'T', '':'C'}
dssp2ss2 = {'H':'H', 'G':'H', 'I':'H', 'E':'E', 'B':'E', 'S':'C', 'T':'C', '':'C'}
stride2ss = {'H':'H', 'E':'E', 'T':'C', 'C':'C', 'G':'H', 'I':'H', 'b':'E', 'B':'E'} 

aatASA = [110.2, 144.1, 140.4, 174.7, 200.7, 78.7, 181.9, 185, 205.7, 183.1, 200.1, 146.4, 141.9, 178.6, 229, 117.2, 138.7, 153.7, 240.5, 213.7];
aaShphobic = [0.81, 0.88, 0.72, 0.46, 0.95, 0.77, 0.54, 1, 0.26, 0.92, 0.81, 0.45, 0.68, 0.43, 0, 0.6, 0.63, 0.92, 0.85, 0.71] 
aaShphobic_c = [x-0.77 for x in aaShphobic]
aa3toaa1 = { a:[b,c,d] for a,b,c,d in zip(aa3code,aa1code,aatASA,aaShphobic)}
aaA2tASA = dict(zip(aa1code, aatASA)) 
aaA2hpo = dict(zip(aa1code, aaShphobic))       
aaA2hpo_c = dict(zip(aa1code, aaShphobic_c))      

def aminoA2tASA(aa):
	try:
        	return aaA2tASA[aa]
	except:
        	return 0
        	
def aminoA2hpo(aa):
	try:
        	return aaA2hpo[aa]
	except:
        	return 0

def aminoA2hpo_c(aa):
	try:
        	return aaA2hpo_c[aa]
	except:
        	return 0

def aminoAAA2A(res):
	try:
		return aaAAA2A[res]
	except:
		return 'B'

def aminoA2AAA(res):
	try:
		return aaA2AAA[res]
	except:
		return 'NSA'

def aminoA2chph74(res):
	try:
		return aa2chph74[res]
	except:
		return 0

def aminoA2chph26(res):
	try:
		return aa2chph26[res]
	except:
		return 0
				
aminoAAA2A = np.vectorize(aminoAAA2A)
aminoA2AAA = np.vectorize(aminoA2AAA)
aminoA2hpo = np.vectorize(aminoA2hpo)
aminoA2hpo_c = np.vectorize(aminoA2hpo_c)
aminoA2tASA = np.vectorize(aminoA2tASA)
aminoA2chph74 = np.vectorize(aminoA2chph74)
aminoA2chph26 = np.vectorize(aminoA2chph26)

#### checks a filename exists and accessible
def file_len(fname):
    i = 0
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


###########################################################################################
### Sequence-related functions
###########################################################################################
def readfasta(filename):    
    header = []; seqlist = []; flag = False; seq = [];
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>'):                
                header.append(line[1:-1].strip() )
                if flag :
                    seqlist.append(''.join(seq))
                    seq = [] 
                else:
                    flag = True
                    seq = []  
            else:                    
                seq.append(line[:-1].replace(' ', ''))
    seqlist.append(''.join(seq))
    seq = pd.DataFrame(header); seq['seq'] = seqlist
    seq.columns = ['header', 'seq']
    return seq

def writefasta(filename, header, seqlist, mode="w") :
	if len(header) == 0 :
		header = ['seq' + str(i) for i in range(len(seqlist))]
	elif len(header) == 1 :
		header = [header + str(i) for i in range(len(seqlist))]	
	file = open(filename, mode)
	for i in range(len(seqlist)):
		strr = '>' + header[i] + '\n' + seqlist[i] + '\n'
		file.write(strr)
	file.close()
	return True

###########################################################################################
# 'Old conservation score program for MSA sequences'
def conservation_score(list_seq):
    f_matrix = [[0 for c in range(20)] for r in range(len(list_seq[0]))]
    r_list = ['G','A','P','V','L','I','M','F','Y','W','S','T','C','N','Q','K','H','R','D','E']
    len_seq = len(list_seq[0])
    seq_count = len(list_seq)
    for j in range(0,len_seq):
        for k in range(0,20):
            for i in range(0,seq_count):
                if list_seq[i][j] == r_list[k]:
                    f_matrix[j][k] = f_matrix[j][k] + 1
    for j in range(0,len_seq):
        for k in range(0,20):
            f_matrix[j][k] = (f_matrix[j][k]/seq_count)
    score = [0]*(len_seq)
    for j in range(0,len_seq):
        for k in range(0,20):
            if f_matrix[j][k] != 0:
                score[j] = score[j] + (f_matrix[j][k]*(math.log(f_matrix[j][k])))
    return score    

def binary2indseg(binseq, cl):
    seglist = []; i = 0
    while i < len(binseq) :
        if binseq[i] == cl:
            seg = []
            while i < len(binseq) and binseq[i] == cl :
                seg.append(i)
                i = i + 1
            seglist.append(seg)
        i = i + 1
    return seglist
def seglist2boolind(seglist, length, sep1=',', sep2='-'): 
    ind = []
    if len(seglist) == 0 :
        segboollist = np.zeros(length, bool)
    else:
        for seg in seglist.split(sep1):
            seg = seg.split(sep2)            
            try:
                st = int(seg[0])
                ed = int(seg[1])
                if st <= ed :
                    ind.extend( [x for x in range( st-1, ed ) ] )
                else:
                    ind.extend( [x for x in range( ed-1, st ) ] )
            except:
                print('seglist2boolind: Bad segment format !!  ',seg)
                return -1
        segboollist = np.zeros(length, bool)
        segboollist[ind] = True
    return segboollist

def ind2bool(ind, length):
    bind = [0 for i in range(length) ]
    for i in ind:
        bind[i] = 1
    return bind

def boolind2seglist(segboollist, sep1=',', sep2='-'): 
    i = 0 ; Flag = 0; seglist = []
    if len(segboollist) <= 0:
        seglist = ''
    else:
        while i < len(segboollist) :
            if segboollist[i] == 1 and Flag == 0:
                Flag = 1; st = i+1
            elif segboollist[i] == 0 and Flag == 1:
                Flag = 0; ed = i;
                seglist.append( str(st)+sep2+str(ed))
            i += 1 
        if Flag == 1:
            ed = i;
            seglist.append( str(st)+sep2+str(ed))
        seglist = ','.join(seglist)
    return seglist    


#posvec = [ [0,1,2,3,4,0,1,2,3], [1,2,3,4,5,2,3,4,5] ]
#posvec = [ [0,0], [1,2]]
def pairpreference(seqlist, mode):
    aa1code = ['A', 'D', 'C', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'] 
    aapair = np.array([ [y+x for x in aa1code] for y in aa1code]).flatten()
    pairpref = pd.DataFrame( index = aapair)
    #pairpref['aa1'] = pairpref.index.str.slice(0,1)
    #pairpref['aa2'] = pairpref.index.str.slice(1,2)
    if mode == 0 : 
        pairpref['f12'] = 0.0
        pairpref['f13'] = 0.0        
        for seq in seqlist:
            for i in range(len(seq)-2):
                try:                     
                    pairpref.ix[ seq[i]+seq[i+1],'f12'] += 1
                    pairpref.ix[ seq[i]+seq[i+2],'f13'] += 1         
                except :#KeyError as exc:
                    #print("Caught KeyError: {}".format(exc))
                    print(seq,i)
            pairpref.ix[ seq[i+1]+seq[i+2],'f12'] += 1
    elif mode == 1:
        posvec = [ [0,1,2,3,4,0,1,2,3], [1,2,3,4,5,2,3,4,5] ]
        clist = []
        for i in range(len(posvec[0])):
            clist.append('f'+str(posvec[0][i]+1)+str(posvec[1][i]+1) )                
            pairpref[clist[i]] = 0.0;       
        for seq in seqlist:
            for i in range(len(posvec[0])):                
                pairpref.ix[ seq[posvec[0][i]]+seq[posvec[1][i]],clist[i]] += 1                               
    pairpref /= pairpref.sum()    
    return pairpref

def calscore_pairpref(seqlist, pairpref, mode):
    if mode == 0 : 
        score = np.zeros((len(seqlist), 2))
        k = -1
        for seq in seqlist:
            k = k + 1; l = len(seq)
            score[k, 0] = np.sum(pairpref.loc[[seq[x]+seq[x+1] for x in np.arange(l-1)],'f12']) 
            score[k, 1] = np.sum(pairpref.loc[[seq[x]+seq[x+2] for x in np.arange(0,l-2)], 'f13'])
    elif mode == 1:
        posvec = [ [0,1,2,3,4,0,1,2,3], [1,2,3,4,5,2,3,4,5] ]
        score = np.zeros((len(seqlist), 9))
        clist = []
        for i in range(len(posvec[0])):
            clist.append('f'+str(posvec[0][i]+1)+str(posvec[1][i]+1) )                
        k = -1   
        for seq in seqlist:
            k = k + 1
            for i in range(len(posvec[0])):
                #seq[posvec[0][i]]+seq[posvec[1][i]]                #clist[i]      
                pairpref.loc[ seq[posvec[0][i]]+seq[posvec[1][i]],clist[i]] 
                score[k, i] = pairpref.loc[ seq[posvec[0][i]]+seq[posvec[1][i]],clist[i]]  
    elif mode == 2:
        fl = pairpref.shape[1]
        score = np.zeros((len(seqlist), fl))
        k = -1
        for seq in seqlist:
            k = k + 1; l = len(seq)
            for f in range(fl):
                score[k, f] = np.sum(pairpref.ix[[seq[x]+seq[x+1] for x in np.arange(l-1)], f]) 
                score[k, f] = np.sum(pairpref.ix[[seq[x]+seq[x+2] for x in np.arange(0,l-2)], f])                      
    return score

def calcompos(seqlist):
    mat = np.zeros((len(seqlist), 20))
    k = -1 
    for seq in seqlist :
        k = k +1
        mat[k, :] = [seq.count(x) for x in aa1code]
    mat = mat/np.array([len(x) for x in seqlist]).reshape((len(seqlist),1))
    return mat


###########################################################################################
### structure-related functions
###########################################################################################

###########################################################################3333############
### 'Protein Structure related Functions'
def assign_mass(s):
    if (s=="C"):
        return 12.0107
    elif (s=="H"):
        return 1.00794
    elif (s=="N"):
        return 14.0067
    elif (s=="O"):
        return 15.9994
    elif (s=="P"):
        return 30.973762
    elif (s=="S"):
        return 32.065
    else:
        print("unidentified atom: error mass not assigned", s)
        return 0
    
def distmatrix_euclid(coords):
    l = coords.shape[0]
    mat = np.zeros((l,l))
    for i in range(l):
        mat[i,:] = np.sum(np.power( coords - coords[i,:], 2), axis=1)    
    return np.power( mat, 0.5)

def dist_euclid(coords1, coords2):
    l1 = coords1.shape[0]
    l2 = coords2.shape[0]
    if l1  > 10000 or l2 > 10000:
        print('calulating distmat using np.float32')
        mat = np.zeros( (l1, l2), dtype=np.float32 )
    else:
        mat = np.zeros((l1,l2) )
    for i in range(l2):
        mat[:,i] = np.sum(np.power( coords1 - coords2[i,:], 2), axis=1)    
    return np.power( mat, 0.5)

def circular_stat(angs, degree=True):
    if degree:
        angs = [np.deg2rad(x) for x in angs]        
    angs_cos = np.sum( [np.cos(a) for a in angs] )
    angs_sin = np.sum( [np.sin(a) for a in angs] )
    n = len(angs)
    angs_cos_n = angs_cos/n
    angs_sin_n = angs_sin/n
    #R = np.sqrt( angs_cos * angs_cos + angs_sin * angs_sin  )
    R_n = np.sqrt( angs_cos_n * angs_cos_n + angs_sin_n*  angs_sin_n  )
    T = np.arctan2( angs_sin_n, angs_cos_n ) 
    sd = np.sqrt( -2 * np.log(R_n) )
    if degree:
        T = np.rad2deg( T )
        sd = np.rad2deg( sd )
    #np.arctan( angs_sin_n/angs_cos_n )
    
    return T, 1-R_n, sd # circular Mean, circular Variance, stddev

def get_entropy_1dang( angs0, nbins=36 ):
    #angs0 = np.arange(-360 , 720, 60 )
    angs = []
    for a in angs0:
        if a < 0:
            angs.append( a + 360 )
        else:
            while a >= 360:
                a -= 360
            angs.append( a )            
    binnos = np.digitize( angs, np.arange(0, 360.5, 360/nbins) )
    binx, freqy = np.unique( binnos, return_counts=True)
    proby = freqy/freqy.sum()
    t = np.log(nbins)
    H = -np.sum( [x * np.log(x) for x in proby] )
    return H/t
    

### 'Read Stride output and write in TSV format'
def rwstride(ifile, ofile):	#script to reformat stride output 
	with open(ifile, 'r') as infile :
		with open(ofile, 'wb') as outfile :
			outfile.write("res\tch\tresno\tsno\tss\tsecstruct\tphi\tpsi\tarea\tpdb\n")
			line = infile.readline()
			while line :
				outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (line[5:8].strip(), line[8:10].strip(), line[10:15].strip(), line[15:20].strip(), line[20:25].strip(), line[25:39].strip(), line[39:49].strip(), line[49:59].strip(), line[59:69].strip(), line[69:79].strip() ))
				line = infile.readline()
	return 1

def readdssp(fname, mode=0, pdbid=''):
    if not os.path.isfile(fname):
        raise IOError('{0} is not a valid file path'.format(fname))   
    else:    
        dssp = open(fname); 
        nl = file_len(fname);  i = 0
        for line in dssp:
            i += i
            if line.startswith('  #  RESIDUE'):
                break
        #### Initialize (efficient for large data)            
        n_atoms =  nl - i   
        data = pd.DataFrame(); k = - 1 
        NUMBER0 = np.ones(n_atoms, int) * -100                  
        NUMBER = np.zeros(n_atoms, int)
        CHAIN = np.zeros(n_atoms, dtype=np.dtype(str))
        SS = np.zeros(n_atoms, dtype=np.dtype(str))
        RES = np.zeros(n_atoms, dtype=np.dtype(str))
        SHEETLABEL = np.zeros(n_atoms, np.array(['a']).dtype.char + '1')
        ACC = np.zeros(n_atoms, float)
        KAPPA = np.zeros(n_atoms, float)
        ALPHA = np.zeros(n_atoms, float)
        PHI = np.zeros(n_atoms, float)
        PSI = np.zeros(n_atoms, float)
        if mode == 1:
            BP1 = np.zeros(n_atoms, int)
            BP2 = np.zeros(n_atoms, int)
            NH_O_1 = np.zeros(n_atoms, int)
            NH_O_1_nrg = np.zeros(n_atoms, float)
            O_HN_1 = np.zeros(n_atoms, int)
            O_HN_1_nrg = np.zeros(n_atoms, float)
            NH_O_2 = np.zeros(n_atoms, int)
            NH_O_2_nrg = np.zeros(n_atoms, float)
            O_HN_2 = np.zeros(n_atoms, int)
            O_HN_2_nrg = np.zeros(n_atoms, float)
            TCO = np.zeros(n_atoms, float)
    
        #ag.setSecstrs(np.zeros(n_atoms, dtype=ATOMIC_FIELDS['secondary'].dtype))    
        for line in dssp:
            if line[13] == '!':
                continue
            resid = line[11]+line[5:10]+line[10].strip()
            #print(resid)
            if resid is None:
                continue
            k = k + 1
            NUMBER0[k] = k+1
            CHAIN[k] = line[10:12].strip(); 
            SS[k] = line[16].strip()
            RES[k] = line[13]
            NUMBER[k] = int(line[5:10])
            SHEETLABEL[k] = line[33].strip()
            ACC[k] = int(line[35:38])
            KAPPA[k] = float(line[91:97])
            ALPHA[k] = float(line[97:103])
            PHI[k] = float(line[103:109])
            PSI[k] = float(line[109:115])
    
            if mode == 1:
                BP1[k] = int(line[25:29])
                BP2[k] = int(line[29:33])
                NH_O_1[k] = int(line[38:45])
                NH_O_1_nrg[k] = float(line[46:50])
                O_HN_1[k] = int(line[50:56])
                O_HN_1_nrg[k] = float(line[57:61])
                NH_O_2[k] = int(line[61:67])
                NH_O_2_nrg[k] = float(line[68:72])
                O_HN_2[k] = int(line[72:78])
                O_HN_2_nrg[k] = float(line[79:83])
                TCO[k] = float(line[85:91])
        data['sno'] = NUMBER0        
        data['ch'] =  CHAIN
        data['resnum'] =  NUMBER
        data['res'] =  RES    
        data['sheet_label'] = SHEETLABEL
        data['ss'] =  SS
        data['acc'] = ACC
        data['kappa'] = KAPPA
        data['alpha'] = ALPHA
        data['phi'] = PHI
        data['psi'] = PSI
        if mode == 1:
            data['bp1'] = BP1
            data['bp2'] = BP2
            data['NH_O_1_index'] = NH_O_1
            data['NH_O_1_energy'] = NH_O_1_nrg
            data['O_NH_1_index'] = O_HN_1
            data['O_NH_1_energy'] = O_HN_1_nrg
            data['NH_O_2_index'] = NH_O_2
            data['NH_O_2_energy'] = NH_O_2_nrg
            data['O_NH_2_index'] = O_HN_2
            data['O_NH_2_energy'] = O_HN_2_nrg
            data['tco'] = TCO
        data['ss1'] = data['ss'].apply(lambda x: dssp2ss1[x])
        data['ss2'] = data['ss'].apply(lambda x: dssp2ss2[x])
        data = data[data.sno != -100]  
        temp = data['res'].apply(aminoA2tASA)
        data['racc'] = data['acc']/temp
        if pdbid != '' :
        	data['pdbid'] = pdbid         
    return data    
'''
to compile: python -m py_compile propylib.py
'''

###############################################################################################################
'Coordinate system related'

def dihedral2SpatialCoords(dihangles):	
#	dihanglesCoords = np.zeros(np.shape(dihangles), np.float)
	dihanglesCoords = [ np.array([ polar2cartesian(1, theta) for theta in row ]).flatten() for row in dihangles ] 
	return np.array(dihanglesCoords)

def polar2cartesian( r, theta):
	x = r * math.sin(theta) 
	y = r * math.cos(theta) 
	return [x,y]

def spherical2cartesian( r, phi, theta):
	x = r * math.sin(phi) * math.cos(theta)
	y = r * math.sin(phi) * math.sin(theta)
	z = math.cos(phi)
	return [x,y,z]
	
def threespherical2cartesian( r, psi, theta, phi):
	x1 = r * math.cos(psi)
	x2 - r * math.sin(psi) * math.cos(theta)
	x3 = r * math.sin(psi) * math.sin(theta) * math.cos(phi)
	x4 = r * math.sin(psi) * math.sin(theta) * math.sin(phi)
	return [x1, x2, x3, x4]

'RMSD in 2 dimension between two lists of coordinates'
def rmsd2(coords2list1, coords2list2): # same length
	n = len(coords2list1)
	#dist = np.empty((n,1), dtype=np.float) 
	rmsd = np.power(np.mean(np.power((coords2list1[:,0] - coords2list2[:,0]),2) + np.power(coords2list1[:,1] - coords2list2[:,1],2) ), 0.5)
	return rmsd
'RMSD in 6 dimension between two lists of coordinates	'
def rmsd6(coords6list1, coords6list2): # same length
	n = len(coords6list1) ; rmsd = 0
	#dist = np.empty((n,1), dtype=np.float) 
	rmsd =  np.power(coords6list1[:,0] - coords6list2[:,0],2)
	for i in range(1, np.shape(coords6list1)[1] ) :
		rmsd =  rmsd + np.power(coords6list1[:,1] - coords6list2[:,1],2)
	rmsd = np.power(np.mean(rmsd ), 0.5)
	return rmsd
'Align - min RMSD position in 2 dimension  - useles '
def aligncoords2(coords2list1, coords2list2):	
	n1 = len(coords2list1)
	n2 = len(coords2list2)
	dist = []
	if n2 > n1 :
		n = n1; ch = 0+1
		for i in range(n2-n1+1) :			
			dist.append(rmsd2(coords2list1, coords2list2[i:n+i,:]))
		rmsd = min(dist)
		ii = [x+1 for x in range(len(dist)) if dist[x] == rmsd]
	elif n1 > n2 :
		n = n2; ch = 0+2
		for i in range(n1-n2+1) :			
			dist.append(rmsd2(coords2list2, coords2list1[i:n+i,:]))
		rmsd = min(dist)
		ii = [x+1 for x in range(len(dist)) if dist[x] == rmsd]
	elif n1 == n2 :
		rmsd = rmsd2(coords2list1, coords2list2)
		ii = 0+1; n = n1; ch = 0+1
	return [rmsd, ch, ii, n ]
'Align - min RMSD position in 6 dimension'
def aligncoords6(coords6list1, coords6list2):	
	n1 = len(coords6list1)
	n2 = len(coords6list2)
	dist = []
	if n2 > n1 :
		n = n1; ch = 0+1
		for i in range(n2-n1+1) :			
			dist.append(rmsd6(coords6list1, coords6list2[i:n+i,:]))
		rmsd = min(dist)
		ii = [x+1 for x in range(len(dist)) if dist[x] == rmsd]
	elif n1 > n2 :
		n = n2; ch = 0+2
		for i in range(n1-n2+1) :			
			dist.append(rmsd6(coords6list2, coords6list1[i:n+i,:]))
		rmsd = min(dist)
		ii = [x+1 for x in range(len(dist)) if dist[x] == rmsd]
	elif n1 == n2 :
		rmsd = rmsd6(coords6list1, coords6list2)
		ii = 0+1; n = n1; ch = 0+1
	return [rmsd, ch, ii, n ]
		
###########################################################################################
'Calculates dihedral angle based on given four points'    
def calcdihedral(coords):
	#coords = [[-7.859, -12.441, 4.87], [-7.698, -12.22, 3.413] , [-6.758, -11.046, 3.157], [-5.621, -11.068, 3.798]]
	#coords = [[-5.621, -11.068, 3.798],[-4.654, -9.957, 3.601] ,[-4.323 ,-9.78 ,2.126], [-3.827, -8.623, 1.782]] 
	#coords = [[-3.827, -8.623, 1.782], [-3.48, -8.37, 0.361], [-2.491, -9.407, -0.147], [-1.874, -10.114, 0.627]] 
	#coords = np.array(coords)
	BA = coords[0,:] - coords[1,:]; CD = coords[3,:] - coords[2,:]
	BC = coords[2,:] - coords[1,:]; bc = BC/np.linalg.norm(BC) 
	ba = BA - np.dot(BA, bc) * bc; 
	cd = CD - np.dot(CD, bc) * bc
	cosang = np.dot(ba,cd)/(np.linalg.norm(ba)*np.linalg.norm(cd))
	N3 = np.cross(ba,cd); 	n3 = N3/np.linalg.norm(N3)
#	n1 = np.cross(bc,BA)
#	n2 = np.cross(bc, CD)
#	cosang = np.dot(n1,n2)/(np.linalg.norm(n1)*np.linalg.norm(n2))
#	N3 = np.cross(n1,n2)
	n3 = N3/np.linalg.norm(N3)
	if np.dot(bc,n3) < 0 :
		ANG =  -np.arccos(cosang)
	else:
		ANG =  np.arccos(cosang)	 
	ANG = math.degrees(ANG)
#	print coords
#	print ANG
	return ANG

###########################################################################################
### PDB object
###########################################################################################

class coord():
	'Special class for coordinates and related function'
	def __init__(self, d=3):
		self.n = 0
		self.dim = d
		self.data = np.empty((self.n, self.dim), dtype=np.float) #### peformance tuning need to be done 
		
	def append(self, newcoord):
		newcoord = np.reshape(newcoord, (1, self.dim))
		self.data = np.append(self.data, newcoord, axis=0)
		self.n +=1 
   
	def extend(self, newcoords):
		newcoords = np.array(newcoords)
		self.data = np.append(self.data, newcoords, axis=0)
		self.n += newcoords.shape[0]
        
	def __getitem__(self, index):
		index1, index2 = index
		if isinstance(index1[0], bool) :
			index1 = [x for x in range(len(index1)) if index1[x] == True ]
		if isinstance(index2[0], bool) :
			index2 = [x for x in range(len(index2)) if index2[x] == True ]						
		return self.data[np.ix_(index1, index2)]
			
	def transform(self, tcoord):
		self.data = self.data + tcoord
        
	def __len__(self):
		return self.n
	# include rotation

class atomgrp_o:
	'Structure for storing pdb data of atoms/hetatoms'
	def __init__(self):
		self.atomno = []		
		self.atomname =[]
		self.atomtype = []
		self.occupancy = []
		self.bfac = []

		self.resno = [] 
		self.resname = []
		self.chain = []
		self.coords = coord()
		
	def add(self, record):
		self.atomno.append(record[0])
		self.atomname.append(record[1])
		self.resname.append(record[2])
		self.chain.append(record[3])
		self.resno.append(record[4])
		self.coords.append( record[5] )
		self.occupancy.append(record[6])
		self.bfac.append(record[7])			
		self.atomtype.append(record[8])
class atomgrp:
	'Structure for storing pdb data of atoms/hetatoms'
	def __init__(self):
		self.atomno = np.empty((0), dtype=int)		
		self.atomname =np.empty((0), dtype='<U3')
		self.atomtype = np.empty((0), dtype='<U2')
		self.occupancy = np.empty((0), dtype=float)
		self.bfac = np.empty((0), dtype=float)

		self.resno = np.empty((0), dtype=int) 
		self.resname = np.empty((0), dtype='<U3')
		self.chain = np.empty((0), dtype='<U2')
		self.coords = coord()
		
	def add(self, record):
		self.atomno = np.append(self.atomno, record[0])
		self.atomname = np.append(self.atomname, record[1])
		self.resname = np.append(self.resname, record[2])
		self.chain = np.append(self.chain, record[3])
		self.resno = np.append(self.resno, record[4])
		self.coords.append( record[5] )
		self.occupancy = np.append(self.occupancy, record[6])
		self.bfac = np.append(self.bfac, record[7])			
		self.atomtype = np.append(self.atomtype, record[8])
		
        
###########################################################################3333############
#pdbfile = '/home/prabakar/Documents/mylib/compile_propylib/4zyo.pdb'
#pdb = PDB(pdbfile, 0, '')
#pdbfile = pdbfile.split('/')[-1].split('.')[0]
#coords = [pdb.atom.coords.data[x,:] for x in range(len(pdb.atom.atomno)) if pdb.atom.chain[x] == ch and pdb.atom.atomtype[x] != 'H' ]
#Data = pd.DataFrame(columns =['pdb', 'ch', 'res1', 'res2', 'dist'])
#for i in coords:
#    data = pd.DataFrame(columns =['pdb', 'ch', 'res1', 'res2', 'dist'])
#    data['dist'] = np.sqrt(np.sum(np.power((pdb.atom.coords.data - pdb.atom.coords.data[i, :]),2),axis = 1))
#    data['res2'] = pdb.atom.resno
#    data['res1'] = i
#    data['pdb'] = pdbfile  
#    data['ch'] = ch
#    Data = pd.concat((Data, data.loc[data.dist <= 8]), axis=0)
#Data.to_csv('')        
    
#status = pdb.listresstretchesCA('A', 4, 0)
#print(status)
#dihangles = pdb.dihedralangles('A')
#coords = dihedral2SpatialCoords(dihangles[0])
#print(coords)
#for row in dihangles[0]:
#	print( row )	

class PDB:
	'Class for storing PDB data'
	def __init__(self, filename, modelno=0, chainid=''):
		self.fname = ((filename.split('/'))[-1].split('.'))[0]
		#print( 'chainid', chainid )
		if len(chainid) == 1 :
			chainid2 = ''.join([' ',chainid])
		self.missingres_seqno = [] 
		self.chains={}
		self.resnos = []
		self.resnames = []
		self.resnames1 = []

		self.atom = atomgrp()
		self.hetatom = atomgrp()

		self.indCalpha = []
		self.indBB = []
		self.indnH = []
		self.res_many_occurence = []
		modelflag = False
		modelcounter = -1
		with open(filename, 'r') as pdbfile :
			line = pdbfile.readline()
			if chainid == '' :			
				while line :
					#print(line) 
	## Read Missing residues			
					if line[:27] == 'REMARK 465   M RES C SSSEQI':
						line = pdbfile.readline()
						while line[:10] == 'REMARK 465' :
							temp = line[21:26].strip()						
							self.missingres_seqno.append(int(temp))							
							line = pdbfile.readline()## Read seqres				
	## Read SeqRes				
					while line[:6] == 'SEQRES':	
						if line[10:12].strip() in self.chains.keys() : 
							self.chains[line[10:12].strip()]  = self.chains[line[10:12].strip()] + line[19:80].split()
						else :
							self.chains[line[10:12].strip()] = line[19:80].split()
							#print(line[10:12].strip())		
							#print(line		)
						line = pdbfile.readline()
					if line[:5] == 'MODEL': 
						print(modelcounter)
						modelcounter += 1
						if modelcounter == modelno:
						    modelflag = True 
					if modelflag and line[:6] == 'ENDMDL' :
						break 
                       
					if modelflag or modelcounter == -1: 	
						#print(line) 
## Read atom coordinates					
						if line[:6] == 'ATOM  ':
							t = [int(line[6:11].strip()), line[12:16].strip(), line[17:20].strip(), line[21:22].strip(), int(line[22:26].strip()), [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())] ]
							for ti in [ line[54:60].strip(), line[60:66].strip() ]:
								if ti != '':
									t.append( float( ti ) )
								else:
									t.append( 0 )
							t.append( line[76:78].strip() )
							self.atom.add( t ) 
									
		## Read hetatom coordinates
						if line[:6] == 'HEATOM':
							t = [int(line[6:11].strip()), line[12:16].strip(), line[17:20].strip(), line[21:22].strip(), int(line[22:26].strip()), [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())] ]
							for ti in [ line[54:60].strip(), line[60:66].strip() ]:
								if ti != '':
									t.append( float( ti ) )
								else:
									t.append( 0 )
							t.append( line[76:78].strip() )							
							self.hetatom.add( t )						
	
					line = pdbfile.readline()
			else:
				while line :			
					if line[:27] == 'REMARK 465   M RES C SSSEQI':
						line = pdbfile.readline()
						while line[:10] == 'REMARK 465' and line[18:20] == chainid2:
							temp = line[21:26].strip()						
							self.missingres_seqno.append(int(temp))							
							line = pdbfile.readline()## Read seqres				
	
					while line[:6] == 'SEQRES':	
						if line[10:12].strip() in self.chains.keys() : 
#							self.chains[line[10:12].strip()].append(line[19:80].split()) 
							self.chains[line[10:12].strip()]  = self.chains[line[10:12].strip()] + line[19:80].split()
						elif line[10:12] == chainid2:
							self.chains[line[10:12].strip()] = line[19:80].split()				
						line = pdbfile.readline()
                        
					if line[:5] == 'MODEL':
						modelflag = True 
					if modelflag and line[:6] == 'ENDMDL' :
						break                     
					if line[:6] == 'ATOM  ' and line[20:22] == chainid2 :	
						t = [int(line[6:11].strip()), line[12:16].strip(), line[17:20].strip(), line[21:22].strip(), int(line[22:26].strip()), [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())] ]
						for ti in [ line[54:60].strip(), line[60:66].strip() ]:
							if ti != '':
								t.append( float( ti ) )
							else:
								t.append( 0 )
						t.append( line[76:78].strip() )
						self.atom.add( t ) 

					if line[:6] == 'HEATOM':
						t = [int(line[6:11].strip()), line[12:16].strip(), line[17:20].strip(), line[21:22].strip(), int(line[22:26].strip()), [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())] ]
						for ti in [ line[54:60].strip(), line[60:66].strip() ]:
							if ti != '':
								t.append( float( ti ) )
							else:
								t.append( 0 )
						t.append( line[76:78].strip() )
						self.hetatom.add( t ) 
					line = pdbfile.readline()
		
		self.chains = np.unique( self.atom.chain ).tolist()
		#self.atom.coords.data = self.atom.coords.data[1:,:]
		temp0 = [[ x for x in range(len(self.atom.atomname)) if (self.atom.atomname[x] == 'CA') and (self.atom.chain[x] == y)] for y in self.chains]	
		## removing multiple occuring residues - choose first one
		for chainx in temp0:			
			temp1 = [ self.atom.resno[y] for y in chainx ]
			temp3 = {}
			for rsi in range(len(temp1)) :
				resno = temp1[rsi]		
				if resno in temp3.keys():
					self.res_many_occurence.append(resno)
				else:
					temp3[resno] = chainx[rsi]
			self.indCalpha.append(list(temp3.values()))
		
		self.resnos = [ [ self.atom.resno[y] for y in x ] for x in self.indCalpha]
		self.resnames = [ [ self.atom.resname[y] for y in x ] for x in self.indCalpha]
		self.resnames1 = [ [ aminoAAA2A(y).tolist() for y in x ] for x in self.resnames]
		self.seqs = [''.join(x) for x in self.resnames1]
		self.indBB = [[x for x in range(len(self.atom.atomname)) if (self.atom.atomname[x] == 'N') or (self.atom.atomname[x] == 'CA') or (self.atom.atomname[x] == 'C')] for y in self.chains]	
		self.indnH = [[x for x in range(len(self.atom.atomname)) if (self.atom.atomname[x] != 'H') or (self.atom.atomname[x] != 'D') ] for y in self.chains]	
#		print(self.chains)
#		print(self.resnos)
#		print(self.resnames)		
	'writes PDB seqres in fasta format'
	'Appends to the file given'
	def fasta(self, fname, ch='', mode='w'):
		with open(fname, mode) as fid :
			if ch == '' :				
				chs = list(self.chains.keys())
			else:
				chs = [ch]
			for ch in chs:
				try:
					seq = ''.join([aminoAAA2A(y) for y in self.chains[ch]])
					fid.write(">%s_%s_%d\n%s\n" % (self.fname,ch,len(seq),seq) )
				except:
					print('Couldn\'t fetch the sequence of the chain: ' + ch)				
		return len(chs)

	'Calculates all the dihedral angles in the PDB and returns n x 3 np.array for each chain in list'
	' tested with dssp results . ... perfect :-)'
	def dihedralangles( self, ch=''):
		dihanglestotal = []; chs0 = [x for x in self.chains.keys()]
		if ch == '' :				
			chs = chs0
		else:
			chs = [ch]
		for ch in chs:
			chno = ([x for x in range(len(chs0)) if chs0[x] == ch ])[0]
			chno = ([x for x in range(len(chs0)) if chs0[x] == ch ])[0]
			bbresno	= [self.atom.resno[x] for x in self.indBB[chno] ]
			bbind = self.indBB[chno]
			bbcoords = self.atom.coords.data[bbind, :]
			bbatomname = [self.atom.atomname[x] for x in bbind]
			allresnos = np.arange( min(self.resnos[chno]), max(self.resnos[chno])+1 )
			dihangles = np.ones(( len(allresnos) ,3), dtype=np.float) * 360.0
			#for i in allresnos:
			#	bbcoords
			for i in range(len(bbind)-3):
				ri = [x for x in allresnos if x == bbresno[i+1] ]
				patternresno = [ bbresno[i], bbresno[i+1], bbresno[i+2], bbresno[i+3] ]
				patternresno = [x-min(patternresno) for x in patternresno]		
				resnoorder =  sum([x == 1 for x in patternresno])
				coords = bbcoords[i:i+4,:]
				if (resnoorder == 1) and (bbatomname[i] == 'N') :
					dihangles[ri, 1] = calcdihedral( coords )		
				elif (resnoorder == 2) and (bbatomname[i] == 'CA') :
					dihangles[ri, 2] = calcdihedral( coords )
				elif (resnoorder == 3) and (bbatomname[i] == 'C') :
					dihangles[ri, 0] = calcdihedral( coords )
				#print(resnoorder, patternresno, coords)
			dihanglestotal.append(dihangles)			
		return dihanglestotal
			# to be continued ...
			
	'Calculates distance from a residues to list of residue in a protein'
	def resdistCA(self, ch1, ch2, res1, reslist2):
		if type(reslist2) == int :		
			reslist2 = [reslist2]
		chno1 = [ x for x,y in enumerate(self.chains.keys()) if y == ch1 ]
		chno2 = [ x for x,y in enumerate(self.chains.keys()) if y == ch2 ]		
		if len(chno1) == 0 or len(chno2) == 0 :
			print("resdistCA -> Error: chains not found in pdb")
			return [], []
		else:
			chno1 = chno1[0]
			chno2 = chno2[0]			
		ind1 = [ x  for x in self.indCalpha[chno1] if self.atom.resno[x] == res1 ]		
		reslist2 = [x for x in reslist2 if  x in self.resnos[chno2] ]
		ind2 = [ x  for x in self.indCalpha[chno2] if self.atom.resno[x] in reslist2 ]
		if len(ind1) == 0 or len(ind2) == 0 :
			print("resdistCA -> Error: residues not found in pdb")
			return [], []
		else:
			ind1 = ind1[0]
		coord1 = self.atom.coords.data[ind1,:]		
		coords2 = np.array([self.atom.coords.data[x,:] for x in ind2])		
		resdist = np.sqrt(np.sum(np.power((coords2 - coord1),2),axis = 1))
#		print(resdist, reslist2)
		return resdist, reslist2

	'Gives list( of indices) of all possible residue stretches of given length'
	'c == 1, give the centre residues, t == 1, give coords instead of indicies, t == 2 gives resno'
	def listresstretchesCA(self, ch='A', l=6, c=0, t=0):
		ll = l - 1 ;
		chno = [ x for x,y in enumerate(self.chains.keys()) if y == ch ][0]
		cc = l/2 - 1
		if c == 0 :		
			indresstretchlist = [self.indCalpha[chno][ii:ii+l] for ii in range(len(self.indCalpha[chno])-ll) ]
		else :
			indresstretchlist = [self.indCalpha[chno][ii+cc] for ii in range(len(self.indCalpha[chno])-ll) ]			

		if t == 0 :
			return indresstretchlist 
		elif t == 1:
			return np.array([self.atom.coords.data[x,:] for x in indresstretchlist ])
		elif t == 2:
			return [self.atom.resno[x] for x in indresstretchlist ]
		else :
			return -1
			
	'calculates contact map between two chains'
	def contactmapCA(self, filename, chain1, chain2):
		if ((not chain1 in self.chains.keys()) and  (chain1 != -1)) or ((not chain2 in self.chains.keys()) and  (chain2 != -1)):			
			print('Chain(s) were not found in the PDB') ;
			return -1
		if chain1 == chain2 :
			if chain1 == -1 :
				coords1 = self.atom.coords.data 
				coords2 = coords1 				
			else:
				index1 =  [x==chain1 for x in self.atom.chain]
				coords1 = self.atom.coords[(index1, [0, 1, 2])]
				index2 =  [x==chain2 for x in self.atom.chain]
				coords2 = self.atom.coords[(index2, [0, 1, 2])]
		else:
			if chain1 ==coord -1:
				coords1 = self.atom.coords.data
			else:
				index1 =  [x==chain1 for x in self.atom.chain]
				coords1 = self.atom.coords[(index1, [0, 1, 2])]
			if chain2 == -1:
				coords2 = self.atom.coords.data	
			else:
				index2 =  [x==chain2 for x in self.atom.chain]
				coords2 = self.atom.coords[(index2, [0, 1, 2])]		
		
		size1 = np.shape(coords1)[0]
		size2 = np.shape(coords2)[0]		
		cont_dist = np.zeros((size1,size2))
		cont_bin = np.zeros((size1,size2))
		i = 0
		for c1 in coords1:
			cont_dist[i,:] = np.sqrt(np.sum(np.power((coords2 - c1),2),axis =1))
			i = i+1
		
		res_number1 =  [self.atom.resno[x] for x in range(len(index1)) if index1[x] == True ] if 'index1' in locals() else self.atom.resno
		res_number2 =  [self.atom.resno[x] for x in range(len(index2)) if index2[x] == True ] if 'index2' in locals() else self.atom.resno
		
		contactfile = open(filename,'w')
		contactfile.write('\t'+'\t'.join([str(x) for x in res_number2])+'\n')
		for i in range(size1):
			contactfile.write(str(res_number1[i])+'\t')
			for j in range(size2):
				if cont_dist[i,j] > 6:
					cont_bin[i,j] = 0
				else:
					cont_bin[i,j] = 1
				contactfile.write (str(int(cont_bin[i,j]))+'\t')				
			contactfile.write ('\n')
		
		print(cont_bin)
#		plt.imshow(cont_bin[size1:1:-1,:], cmap=cm.hot)
		plt.imshow(cont_bin)
		plt.colorbar()
		plt.show()		
#		plt.set_xticklabels(res_number2)		
		return cont_dist
	'calculates contact map between two chains'
	def contactmapCA1(self, chain, distcutoff):
		if ((not chain in self.chains.keys()) and  (chain != -1)) :			
			print('Chain(s) were not found in the PDB') ;
			return -1
		index = [(x==chain) and (y!='H') for x,y in zip(self.atom.chain,self.atom.atomtype)]
		coords1 = self.atom.coords[(index, [0, 1, 2])]
		coords2 = coords1        
		resnum = self.atom.resno[index]
		reslist = np.unique(resnum).tolist()
		size1 = np.shape(coords1)[0]
		size2 = np.shape(coords2)[0]		
		cont_dist = np.zeros((size1,size2))
		cont_bin = np.zeros((len(reslist),len(reslist)))
		i = 0
		for c1 in coords1:
			cont_dist[i,:] = np.sqrt(np.sum(np.power((coords2 - c1),2),axis =1))
			i = i+1
		
#		contactfile = open(filename,'w')
#		contactfile.write('\t'+'\t'.join([str(x) for x in res_number2])+'\n')
		for i in range(len(reslist)):
#			contactfile.write(str(res_number1[i])+'\t')
			ind1 = np.where(resnum==reslist[i])
			for j in range(len(reslist)):
			    ind2 = np.where(resnum==reslist[j])
			    cont_bin[i,j] = np.sum(cont_dist[np.meshgrid( ind1, ind2)] <= distcutoff)
			#	contactfile.write (str(int(cont_bin[i,j]))+'\t')				
#			contactfile.write ('\n')
	
		return cont_dist,cont_bin
	'Calculates network based on distance between residues in a chain'
#	def networkCA(self, ch1, d):
#		if type(reslist2) == int :		
#			reslist2 = [reslist2]
#		chno1 = [ x for x,y in enumerate(self.chains.keys()) if y == ch1 ]
#		chno2 = [ x for x,y in enumerate(self.chains.keys()) if y == ch2 ]		
#		if len(chno1) == 0 or len(chno2) == 0 :
#			print("resdistCA -> Error: chains not found in pdb")
#			return [], []
#		else:
#			chno1 = chno1[0]
#			chno2 = chno2[0]			
#		ind1 = [ x  for x in self.indCalpha[chno1] if self.atom.resno[x] == res1 ]		
#		reslist2 = [x for x in reslist2 if  x in self.resnos[chno2] ]
#		ind2 = [ x  for x in self.indCalpha[chno2] if self.atom.resno[x] in reslist2 ]
#		if len(ind1) == 0 or len(ind2) == 0 :
#			print("resdistCA -> Error: residues not found in pdb")
#			return [], []
#		else:
#			ind1 = ind1[0]
#		coord1 = self.atom.coords.data[ind1,:]		
#		coords2 = np.array([self.atom.coords.data[x,:] for x in ind2])		
#		resdist = np.sqrt(np.sum(np.power((coords2 - coord1),2),axis = 1))
##		print(resdist, reslist2)
#		return resdist, reslist2

class PDB2:
	'Class for storing PDB data'
	def __init__(self, filename, modelno=0, chainid=''):
		self.fname = ((filename.split('/'))[-1].split('.'))[0]
		print(chainid)
		if len(chainid) == 1 :
			chainid2 = ''.join([' ',chainid])
		self.missingres_seqno = [] 
		self.chains={}
		self.resnos = []
		self.resnames = []
		self.resnames1 = []

		self.atom = atomgrp()
		self.hetatom = atomgrp()

		self.indCalpha = []
		self.indBB = []
		self.indnH = []
		self.res_many_occurence = []
		modelflag = False
		with open(filename, 'r') as pdbfile :
			line = pdbfile.readline()
			if chainid == '' :
				while line :
	## Read Missing residues			
					if line[:27] == 'REMARK 465   M RES C SSSEQI':
						line = pdbfile.readline()
						while line[:10] == 'REMARK 465' :
							temp = line[21:26].strip()						
							self.missingres_seqno.append(int(temp))							
							line = pdbfile.readline()## Read seqres				
	## Read SeqRes				
					while line[:6] == 'SEQRES':	
						if line[10:12].strip() in self.chains.keys() : 
							self.chains[line[10:12].strip()]  = self.chains[line[10:12].strip()] + line[19:80].split()
						else :
							self.chains[line[10:12].strip()] = line[19:80].split()
							#print(line[10:12].strip())		
							#print(line		)
						line = pdbfile.readline()
					if line[:5] == 'MODEL':
						modelflag = True 
					if modelflag and line[:6] == 'ENDMDL' :
						break 					 
	## Read atom coordinates					
					if line[:6] == 'ATOM  ':
						self.atom.add([int(line[6:11].strip()), line[12:16].strip(), line[16:20].strip(), line[20:22].strip(), int(line[22:26].strip()), [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())], float(line[54:62].strip()), float(line[62:68].strip()), line[76:78].strip() ]  ) 
								
	## Read hetatom coordinates
				#	if line[:6] == 'HEATOM':
				#		self.hetatom.add([int(line[6:11].strip()), line[12:16].strip(), line[17:20].strip(), line[21:22].strip(), int(line[22:26].strip()), [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())], float(line[54:60].strip()), float(line[60:66].strip()), line[76:78].strip() ] )						

					line = pdbfile.readline()
			else:
				while line :			
					if line[:27] == 'REMARK 465   M RES C SSSEQI':
						line = pdbfile.readline()
						while line[:10] == 'REMARK 465' and line[18:20] == chainid2:
							temp = line[21:26].strip()						
							self.missingres_seqno.append(int(temp))							
							line = pdbfile.readline()## Read seqres				
	
					while line[:6] == 'SEQRES':	
						if line[10:12].strip() in self.chains.keys() : 
#							self.chains[line[10:12].strip()].append(line[19:80].split()) 
							self.chains[line[10:12].strip()]  = self.chains[line[10:12].strip()] + line[19:80].split()
						elif line[10:12] == chainid2:
							self.chains[line[10:12].strip()] = line[19:80].split()				
						line = pdbfile.readline()
                        
					if line[:5] == 'MODEL':
						modelflag = True 
					if modelflag and line[:6] == 'ENDMDL' :
						break                     
					if line[:6] == 'ATOM  ' and line[20:22] == chainid2 :						
						self.atom.add([int(line[6:11].strip()), line[12:16].strip(), line[17:20].strip(), line[20:22].strip(), int(line[22:26].strip()), [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())], float(line[54:60].strip()), float(line[60:66].strip()), line[76:78].strip() ]  ) 

					if line[:6] == 'HEATOM':
						self.hetatom.add([int(line[6:11].strip()), line[12:16].strip(), line[17:20].strip(), line[20:22].strip(), int(line[22:26].strip()), [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())], float(line[54:60].strip()), float(line[60:66].strip()), line[76:78].strip() ] )						
					line = pdbfile.readline()
				
		#self.atom.coords.data = self.atom.coords.data[1:,:]
		temp0 = [[ x for x in range(len(self.atom.atomname)) if (self.atom.atomname[x] == 'CA') and (self.atom.chain[x] == y)] for y in self.chains]	
		## removing multiple occuring residues - choose first one
		for chainx in temp0:			
			temp1 = [ self.atom.resno[y] for y in chainx ]
			temp3 = {}
			for rsi in range(len(temp1)) :
				resno = temp1[rsi]		
				if resno in temp3.keys():
					self.res_many_occurence.append(resno)
				else:
					temp3[resno] = chainx[rsi]
			self.indCalpha.append(list(temp3.values()))
		
		self.resnos = [ [ self.atom.resno[y] for y in x ] for x in self.indCalpha]
		self.resnames = [ [ self.atom.resname[y] for y in x ] for x in self.indCalpha]
		self.resnames1 = [ [ aminoAAA2A(y) for y in x ] for x in self.resnames]
		self.indBB = [[x for x in range(len(self.atom.atomname)) if (self.atom.atomname[x] == 'N') or (self.atom.atomname[x] == 'CA') or (self.atom.atomname[x] == 'C')] for y in self.chains]	
		self.indnH = [[x for x in range(len(self.atom.atomname)) if (self.atom.atomname[x] != 'H') or (self.atom.atomname[x] != 'D') ] for y in self.chains]	
#		print(self.chains)
#		print(self.resnos)
#		print(self.resnames)		
	'writes PDB seqres in fasta format'
	'Appends to the file given'
	def fasta(self, fname, ch='', mode='w'):
		with open(fname, mode) as fid :
			if ch == '' :				
				chs = list(self.chains.keys())
			else:
				chs = [ch]
			for ch in chs:
				try:
					seq = ''.join([aminoAAA2A(y) for y in self.chains[ch]])
					fid.write(">%s_%s_%d\n%s\n" % (self.fname,ch,len(seq),seq) )
				except:
					print('Couldn\'t fetch the sequence of the chain: ' + ch)				
		return len(chs)

	'Calculates all the dihedral angles in the PDB and returns n x 3 np.array for each chain in list'
	' tested with dssp results . ... perfect :-)'
	def dihedralangles( self, ch=''):
		dihanglestotal = []; chs0 = [x for x in self.chains.keys()]
		if ch == '' :				
			chs = chs0
		else:
			chs = [ch]
		for ch in chs:
			chno = ([x for x in range(len(chs0)) if chs0[x] == ch ])[0]
			chno = ([x for x in range(len(chs0)) if chs0[x] == ch ])[0]
			bbresno	= [self.atom.resno[x] for x in self.indBB[chno] ]
			bbind = self.indBB[chno]
			bbcoords = self.atom.coords.data[bbind, :]
			bbatomname = [self.atom.atomname[x] for x in bbind]
			allresnos = np.arange( min(self.resnos[chno]), max(self.resnos[chno])+1 )
			dihangles = np.ones(( len(allresnos) ,3), dtype=np.float) * 360.0
			#for i in allresnos:
			#	bbcoords
			for i in range(len(bbind)-3):
				ri = [x for x in allresnos if x == bbresno[i+1] ]
				patternresno = [ bbresno[i], bbresno[i+1], bbresno[i+2], bbresno[i+3] ]
				patternresno = [x-min(patternresno) for x in patternresno]		
				resnoorder =  sum([x == 1 for x in patternresno])
				coords = bbcoords[i:i+4,:]
				if (resnoorder == 1) and (bbatomname[i] == 'N') :
					dihangles[ri, 1] = calcdihedral( coords )		
				elif (resnoorder == 2) and (bbatomname[i] == 'CA') :
					dihangles[ri, 2] = calcdihedral( coords )
				elif (resnoorder == 3) and (bbatomname[i] == 'C') :
					dihangles[ri, 0] = calcdihedral( coords )
				#print(resnoorder, patternresno, coords)
			dihanglestotal.append(dihangles)			
		return dihanglestotal
			# to be continued ...
			
	'Calculates distance from a residues to list of residue in a protein'
	def resdistCA(self, ch1, ch2, res1, reslist2):
		if type(reslist2) == int :		
			reslist2 = [reslist2]
		chno1 = [ x for x,y in enumerate(self.chains.keys()) if y == ch1 ]
		chno2 = [ x for x,y in enumerate(self.chains.keys()) if y == ch2 ]		
		if len(chno1) == 0 or len(chno2) == 0 :
			print("resdistCA -> Error: chains not found in pdb")
			return [], []
		else:
			chno1 = chno1[0]
			chno2 = chno2[0]			
		ind1 = [ x  for x in self.indCalpha[chno1] if self.atom.resno[x] == res1 ]		
		reslist2 = [x for x in reslist2 if  x in self.resnos[chno2] ]
		ind2 = [ x  for x in self.indCalpha[chno2] if self.atom.resno[x] in reslist2 ]
		if len(ind1) == 0 or len(ind2) == 0 :
			print("resdistCA -> Error: residues not found in pdb")
			return [], []
		else:
			ind1 = ind1[0]
		coord1 = self.atom.coords.data[ind1,:]		
		coords2 = np.array([self.atom.coords.data[x,:] for x in ind2])		
		resdist = np.sqrt(np.sum(np.power((coords2 - coord1),2),axis = 1))
#		print(resdist, reslist2)
		return resdist, reslist2

	'Gives list( of indices) of all possible residue stretches of given length'
	'c == 1, give the centre residues, t == 1, give coords instead of indicies, t == 2 gives resno'
	def listresstretchesCA(self, ch='A', l=6, c=0, t=0):
		ll = l - 1 ;
		chno = [ x for x,y in enumerate(self.chains.keys()) if y == ch ][0]
		cc = l/2 - 1
		if c == 0 :		
			indresstretchlist = [self.indCalpha[chno][ii:ii+l] for ii in range(len(self.indCalpha[chno])-ll) ]
		else :
			indresstretchlist = [self.indCalpha[chno][ii+cc] for ii in range(len(self.indCalpha[chno])-ll) ]			

		if t == 0 :
			return indresstretchlist 
		elif t == 1:
			return np.array([self.atom.coords.data[x,:] for x in indresstretchlist ])
		elif t == 2:
			return [self.atom.resno[x] for x in indresstretchlist ]
		else :
			return -1
			
	'calculates contact map between two chains'
	def contactmapCA(self, filename, chain1, chain2):
		if ((not chain1 in self.chains.keys()) and  (chain1 != -1)) or ((not chain2 in self.chains.keys()) and  (chain2 != -1)):			
			print('Chain(s) were not found in the PDB') ;
			return -1
		if chain1 == chain2 :
			if chain1 == -1 :
				coords1 = self.atom.coords.data 
				coords2 = coords1 				
			else:
				index1 =  [x==chain1 for x in self.atom.chain]
				coords1 = self.atom.coords[(index1, [0, 1, 2])]
				index2 =  [x==chain2 for x in self.atom.chain]
				coords2 = self.atom.coords[(index2, [0, 1, 2])]
		else:
			if chain1 ==coord -1:
				coords1 = self.atom.coords.data
			else:
				index1 =  [x==chain1 for x in self.atom.chain]
				coords1 = self.atom.coords[(index1, [0, 1, 2])]
			if chain2 == -1:
				coords2 = self.atom.coords.data	
			else:
				index2 =  [x==chain2 for x in self.atom.chain]
				coords2 = self.atom.coords[(index2, [0, 1, 2])]		
		
		size1 = np.shape(coords1)[0]
		size2 = np.shape(coords2)[0]		
		cont_dist = np.zeros((size1,size2))
		cont_bin = np.zeros((size1,size2))
		i = 0
		for c1 in coords1:
			cont_dist[i,:] = np.sqrt(np.sum(np.power((coords2 - c1),2),axis =1))
			i = i+1
		
		res_number1 =  [self.atom.resno[x] for x in range(len(index1)) if index1[x] == True ] if 'index1' in locals() else self.atom.resno
		res_number2 =  [self.atom.resno[x] for x in range(len(index2)) if index2[x] == True ] if 'index2' in locals() else self.atom.resno
		
		contactfile = open(filename,'w')
		contactfile.write('\t'+'\t'.join([str(x) for x in res_number2])+'\n')
		for i in range(size1):
			contactfile.write(str(res_number1[i])+'\t')
			for j in range(size2):
				if cont_dist[i,j] > 6:
					cont_bin[i,j] = 0
				else:
					cont_bin[i,j] = 1
				contactfile.write (str(int(cont_bin[i,j]))+'\t')				
			contactfile.write ('\n')
		
		print(cont_bin)
#		plt.imshow(cont_bin[size1:1:-1,:], cmap=cm.hot)
		plt.imshow(cont_bin)
		plt.colorbar()
		plt.show()		
#		plt.set_xticklabels(res_number2)		
		return cont_dist
	'calculates contact map between two chains'
	def contactmapCA1(self, chain, distcutoff):
		if ((not chain in self.chains.keys()) and  (chain != -1)) :			
			print('Chain(s) were not found in the PDB') ;
			return -1
		index = [(x==chain) and (y!='H') for x,y in zip(self.atom.chain,self.atom.atomtype)]
		coords1 = self.atom.coords[(index, [0, 1, 2])]
		coords2 = coords1        
		resnum = self.atom.resno[index]
		reslist = np.unique(resnum).tolist()
		size1 = np.shape(coords1)[0]
		size2 = np.shape(coords2)[0]		
		cont_dist = np.zeros((size1,size2))
		cont_bin = np.zeros((len(reslist),len(reslist)))
		i = 0
		for c1 in coords1:
			cont_dist[i,:] = np.sqrt(np.sum(np.power((coords2 - c1),2),axis =1))
			i = i+1
		
#		contactfile = open(filename,'w')
#		contactfile.write('\t'+'\t'.join([str(x) for x in res_number2])+'\n')
		for i in range(len(reslist)):
#			contactfile.write(str(res_number1[i])+'\t')
			ind1 = np.where(resnum==reslist[i])
			for j in range(len(reslist)):
			    ind2 = np.where(resnum==reslist[j])
			    cont_bin[i,j] = np.sum(cont_dist[np.meshgrid( ind1, ind2)] <= distcutoff)
			#	contactfile.write (str(int(cont_bin[i,j]))+'\t')				
#			contactfile.write ('\n')
	
		return cont_dist,cont_bin
    
#########################################################################################################
### Machine learning related - use sklearn instead
#########################################################################################################

def getScores(estimator, X, y):
    yp = estimator.predict(X)      
    tp = sum( yp[y==1] == 1 )
    tn = sum( yp[y==0] == 0 )
    Npos =  sum(y==1) ; Nneg = sum(y==0);
    nan = float('nan')   
    TP = tp; TN = tn; FN = Npos-TP ; FP = Nneg-TN  
    mcc = ( TP * TN - FP * FN )/( (TP+FP) * (TP+FN) * (TN+FP) * (TN+FN) )**0.5          
    if Npos == 0 and Nneg == 0:
        return [ tp, tn, Npos, Nneg, nan, nan,nan,  mcc]
    if Npos == 0:
        return [ tp, tn, Npos, Nneg, nan, tn/Nneg * 100,(tp+tn)/(Npos+Nneg)*100,  mcc]
    elif Nneg == 0:
        return [ tp, tn, Npos, Nneg, tp/Npos * 100, nan,(tp+tn)/(Npos+Nneg)*100,  mcc]
    else:
        return [ tp, tn, Npos, Nneg, tp/Npos * 100, tn/Nneg * 100,(tp+tn)/(Npos+Nneg)*100,  mcc]
def getScores_direct(y, yp):
    y = np.array(y)
    yp = np.array(yp)
    tp = sum( yp[y==1] == 1 )
    tn = sum( yp[y==0] == 0 )
    Npos =  sum(y==1) ; Nneg = sum(y==0);
    nan = float('nan')   
    TP = tp; TN = tn; FN = Npos-TP ; FP = Nneg-TN  
    mcc = ( TP * TN - FP * FN )/( (TP+FP) * (TP+FN) * (TN+FP) * (TN+FN) )**0.5          
    tpr = tp/Npos * 100
    tnr = tn/Nneg * 100
    if Npos == 0 and Nneg == 0:
        return [ tp, tn, Npos, Nneg, nan, nan,nan,  mcc]
    if Npos == 0:
        return [ tp, tn, Npos, Nneg, nan, tnr, tnr/2, (tp+tn)/(Npos+Nneg)*100,  mcc]
    elif Nneg == 0:
        return [ tp, tn, Npos, Nneg, tpr, nan, tpr/2, (tp+tn)/(Npos+Nneg)*100,  mcc]
    else:
        return [ tp, tn, Npos, Nneg, tpr, tnr, (tpr + tnr)/2, (tp+tn)/(Npos+Nneg)*100,  mcc]
def getScores_singlevalue(estimator, X, y):
    yp = estimator.predict(X)      
    tp = sum( yp[y==1] == 1 )
    tn = sum( yp[y==0] == 0 )
    Npos =  sum(y==1) ; Nneg = sum(y==0);
    nan = float('nan')   
    TP = tp; TN = tn; FN = Npos-TP ; FP = Nneg-TN  
    mcc = ( TP * TN - FP * FN )/( (TP+FP) * (TP+FN) * (TN+FP) * (TN+FN) )**0.5          
    return mcc  
def cross_val(classifier, X, y, k):
    kf = KFold(n_splits=k,shuffle=True); score = []
    Xpos = X[y.ravel()==1, :]; ypos = y[y==1];
    Xneg = X[y.ravel()==0, :]; yneg = y[y==0];        
    for (indtrp, indtep), (indtrn, indten) in zip(kf.split(Xpos), kf.split(Xneg)):
       # print(len(indtrp), len(indtep), len(indtrn), len(indten))
        Xtr = np.concatenate( (Xpos[indtrp, :],Xneg[indtrn, :] ), axis=0 )
        ytr = np.concatenate( (ypos[indtrp], yneg[indtrn] ), axis=0 )
        Xte = np.concatenate( (Xpos[indtep, :],Xneg[indten, :] ), axis=0 )
        yte = np.concatenate( (ypos[indtep], yneg[indten] ), axis=0 )
        estimator = classifier.fit(Xtr, ytr.ravel())
        score.append( getScores(estimator, Xte, yte.ravel()) )
    return np.array(score)

def jackknife(classifier, X, y):
    k = len(y)-1
    kf = KFold(n_splits=k,shuffle=True); score = []            
    for indtr, indte in kf.split(X):        
        estimator = classifier.fit(X[indtr, :], y[indtr].ravel())
        score.append( getScores(estimator, X[indte, :], y[indte].ravel() ) )
    return np.array(score)

def pandas2arff(df,filename,wekaname = "pandasdata",cleanstringdata=True,cleannan=True):
    """
    converts the pandas dataframe to a weka compatible file
    df: dataframe in pandas format
    filename: the filename you want the weka compatible file to be in
    wekaname: the name you want to give to the weka dataset (this will be visible to you when you open it in Weka)
    cleanstringdata: clean up data which may have spaces and replace with "_", special characters etc which seem to annoy Weka. 
                     To suppress this, set this to False
    cleannan: replaces all nan values with "?" which is Weka's standard for missing values. 
              To suppress this, set this to False
    """
    import re
    
    def cleanstring(s):
        if s!="?":
            return re.sub('[^A-Za-z0-9]+', "_", str(s))
        else:
            return "?"
            
    dfcopy = df #all cleaning operations get done on this copy

    
    if cleannan!=False:
        dfcopy = dfcopy.fillna(-999999999) #this is so that we can swap this out for "?"
        #this makes sure that certain numerical columns with missing values don't get stuck with "object" type
 
    f = open(filename,"w")
    arffList = []
    arffList.append("@relation " + wekaname + "\n")
    #look at each column's dtype. If it's an "object", make it "nominal" under Weka for now (can be changed in source for dates.. etc)
    for i in range(df.shape[1]):
        if dfcopy.dtypes[i]=='O' or (df.columns[i] in ["Class","CLASS","class"]):
            if cleannan!=False:
                dfcopy.iloc[:,i] = dfcopy.iloc[:,i].replace(to_replace=-999999999, value="?")
            if cleanstringdata!=False:
                dfcopy.iloc[:,i] = dfcopy.iloc[:,i].apply(cleanstring)
            _uniqueNominalVals = [str(_i) for _i in np.unique(dfcopy.iloc[:,i])]
            _uniqueNominalVals = ",".join(_uniqueNominalVals)
            _uniqueNominalVals = _uniqueNominalVals.replace("[","")
            _uniqueNominalVals = _uniqueNominalVals.replace("]","")
            _uniqueValuesString = "{" + _uniqueNominalVals +"}" 
            arffList.append("@attribute " + df.columns[i] + _uniqueValuesString + "\n")
        else:
            arffList.append("@attribute " + df.columns[i] + " real\n") 
            #even if it is an integer, let's just deal with it as a real number for now
    arffList.append("@data\n")           
    for i in range(dfcopy.shape[0]):#instances
        _instanceString = ""
        for j in range(df.shape[1]):#features
                if dfcopy.dtypes[j]=='O':
                    _instanceString+="\"" + str(dfcopy.iloc[i,j]) + "\""
                else:
                    _instanceString+=str(dfcopy.iloc[i,j])
                if j!=dfcopy.shape[1]-1:#if it's not the last feature, add a comma
                    _instanceString+=","
        _instanceString+="\n"
        if cleannan!=False:
            _instanceString = _instanceString.replace("-999999999.0","?") #for numeric missing values
            _instanceString = _instanceString.replace("\"?\"","?") #for categorical missing values
        arffList.append(_instanceString)
    f.writelines(arffList)
    f.close()
    del dfcopy
    return True

#java -cp weka.classifiers.meta.RotationForest -T focuslib_features.arff -l ../../model_RotationForest_dEnCcv10_02Jun17.model -p 0 > weka_focuslib_modelRF02Jun17.out

## Read weka output
def readwekaout(filename, mode):
    if mode == 0:
        wekadata = pd.read_csv(filename, skiprows=4, sep='\t')
        #wekadata.iloc[:,0].str.extract( '[\s]*([\d]+)[\s]*([\d]):([\d])[\s]*([\d]):([\d])[\s]*([+]*)[\s]*([-.\d]*)[\s]*' )  
        wekadata = wekadata.iloc[:,0].str.extract( '([\d]+)[\s]*([\d]):([A-Za-z?]*)[\s]*([\d]):([A-Za-z?]*)[\s]*([+]*)[\s]*([-.\d]*)[\s]*' )  
        wekadata.columns = ['sno', 'actualno', 'actualcl', 'predno', 'predcl', 'iserror', 'prob']
        wekadata['prob'] = wekadata.prob.astype('float32')
    return wekadata

def forward_selection(classifier, X, y, k):
    k = k + 1
    remaining = list(np.arange( X.shape[1] ))     
    selected = []
    current_score, best_new_score = 0.0, 0.0
    while len(remaining) > 0 and len(selected) < k and current_score == best_new_score:        
        scores1 = [] ; candidates1 = []; scoresx1 = [] ;
        for candidate in remaining :
            selected1 = list(selected)
            selected1.append(candidate)
            #scorex = np.nanmean(jackknife(classifier, X[:,selected1], y), axis=0)[4:8]
            scorex = np.nanmean(cross_val(classifier, X[:,selected1], y, 5), axis=0)[4:8]
            score = scorex[3]
            #score = -np.sum(np.power( 1.0 - scorex[0:2], 2))
            #score = (100.0-scorex[0])**2 + (100.0-scorex[1])**2                        
            scores1.append(score); candidates1.append(candidate); scoresx1.append(scorex);
        best_new_score = max(scores1)
        best_candidate = [candidates1[x] for x in range(len(candidates1)) if scores1[x] == best_new_score][0]        
        best_scorex = [scoresx1[x] for x in range(len(candidates1)) if scores1[x] == best_new_score][0]        
        if current_score < best_new_score:
            remaining.remove(best_candidate)
            selected.append(best_candidate)
            current_score = best_new_score
        print(current_score, best_new_score, best_scorex, len(selected), selected)    
    return selected



###########################################################################################
#### Experimental
###########################################################################3333############ 

###########################################################################3333############
def findoptk(coords, limit):
    mstd = []
    for k in range(1,limit):
        kmeans = KMeans(n_clusters=k, random_state=0).fit(coords)
        cl = kmeans.labels_ ;   clcentre = kmeans.cluster_centers_
        mstd.append( np.sum( [ np.sum(np.power(np.sum( np.power( coords[cl == i, :] - clcentre[i], 2), axis = 1), 0.5)) for i in range(k) ] ) )
    kk = np.arange(1, limit);     i = 1; mstd = np.array(mstd)
    mstd = mstd/max(mstd); delmstd =np.gradient(mstd) ; del2mstd =np.gradient(delmstd)       
    while(del2mstd[i] > 0.0001): i+=1; 
#p1 = np.polyfit(kk[0:5], delmstd[0:5], 1);
#p2 = np.polyfit(kk[-10:], delmstd[-10:], 1)
##plt.plot(kk,mstd); 
#plt.plot(kk,delmstd ); plt.hold(True); plt.plot(kk[0:12], np.polyval(p1, kk[0:12]) );
#plt.plot(kk[-50:], np.polyval(p2, kk[-50:]) )
    return kk[i], [mstd, delmstd, del2mstd]


'K-means Clustering - not verified/completed'
## datasciencelab.wordpress.com/2013/12/12/clustering-with-k-means-in-python/  Lloyds algorithm for k-means
def cluster_points(X, mu):
    clusters  = {}
    for x in X:
        bestmukey = min([(i[0], np.linalg.norm(x-mu[i[0]])) \
                    for i in enumerate(mu)], key=lambda t:t[1])[0]
        try:
            clusters[bestmukey].append(x)
        except KeyError:
            clusters[bestmukey] = [x]
    return clusters
 
def reevaluate_centers(mu, clusters):
    newmu = []
    keys = sorted(clusters.keys())
    for k in keys:
        newmu.append(np.mean(clusters[k], axis = 0))
    return newmu
 
def has_converged(mu, oldmu):
    return (set([tuple(a) for a in mu]) == set([tuple(a) for a in oldmu]) )
 
def find_centers(X, K):
    # Initialize to K random centers
    oldmu = random.sample(X, K)
    mu = random.sample(X, K)
    while not has_converged(mu, oldmu):
        oldmu = mu
        # Assign all points in X to clusters
        clusters = cluster_points(X, mu)
        # Reevaluate centers
        mu = reevaluate_centers(oldmu, clusters)
    return(mu, clusters)
    
###########################################################################3333############    
###########################################################################3333############
## datasciencelab.wordpress.com/2013/12/12/clustering-with-k-means-in-python/  Lloyds algorithm for k-means
## customized for dihedralcoordinates rmsd matrix 
def cluster_points(X, mu, distmat):
    clusters  = {}
    for x in X:
        bestmukey = min([(i[0], np.linalg.norm(x-mu[i[0]])) \
                    for i in enumerate(mu)], key=lambda t:t[1])[0]
        try:
            clusters[bestmukey].append(x)
        except KeyError:
            clusters[bestmukey] = [x]
    return clusters
 
def reevaluate_centers(mu, clusters):
    newmu = []
    keys = sorted(clusters.keys())
    for k in keys:
        newmu.append(np.mean(clusters[k], axis = 0))
    return newmu
 
def has_converged(mu, oldmu):
    return (set([tuple(a) for a in mu]) == set([tuple(a) for a in oldmu]) )
 
def find_centers(X, K, distmat):
    # Initialize to K random centers
    oldmu = random.sample(X, K)
    mu = random.sample(X, K)
    while not has_converged(mu, oldmu):
        oldmu = mu
        # Assign all points in X to clusters
        clusters = cluster_points(X, mu)
        # Reevaluate centers
        mu = reevaluate_centers(oldmu, clusters)
    return(mu, clusters)
    
###########################################################################3333############      

### 'Multiple liner Regression' - sample code
### 'X and Y must be np.array'
alpha = 0.05 ; ndigits = 3
def mlregress(X, Y):
#	Y = np.matrix(data_array[:,0]).T
#	X = np.matrix(data_array[:, vvlist])
	n = np.shape(X) ; p = n[1]+1 ; n = n[0]
	## Multiple Regression
	X = np.append(X, np.ones(np.shape(Y)), axis=1)
	XT = X.T
	if np.linalg.det(XT * X) == 0:
		text = str(k) + '\t{' + ','.join([Header[x] for x in vvlist]) + '}\t' + '## Error: determinant Zero ##'  + '\t' +  ' '  + '\t' +  ' ' + '\t' +  ' '  + '\t' +   ' ' + '\t' +  ' '+ '\t' +  ' ' + '\n'
		print(text); f.write(text) ;
		return -1
		
	XX = (XT * X).I
	B = XX * XT * Y
	Yp = X * B
	## Metrics
	B = list(B.getA1())
	SSE = np.sum(np.power(Y - Yp, 2))
	sig2 = SSE/(n-p);		covB = sig2 * XX;
	SST = np.sum(np.power(Y - np.sum(Y)/n, 2));		SSR = SST - SSE;
	out.R2 = 1 - SSE/SST ;		out.R = pow(R2, 0.5)

	## Stats
	out.tcal = [B[i]/pow(sig2 * covB[i,i], 0.5)  for i in range(len(B))]
	out.talpha = scipy.stats.t.ppf(1-alpha/2.0, n-p)
	out.Brange = [ str( round(B[i] - out.talpha * pow(covB[i,i], 0.5), ndigits))+'-'+str(round(B[i] + out.talpha * pow(covB[i,i], 0.5),ndigits)) for i in range(len(B)) ] 
	#Or whatever you want your alpha to be.
	out.Fp_value = 2 * min(scipy.stats.f.cdf(SSR/SSE, p-1, n-p), scipy.stats.f.sf(SSR/SSE, p-1, n-p))
	return out

######################################################################################

### Segment overlap score - SOV 
### not generalized functions
### todo: add documentation and generalize it
def calSOV_corr(predicted, minmatch=[1,1]):
    #winj = 4
    segs_act_list = []
    segs_act_list1 = []
    for i in range(predicted.shape[0]):
        a = predicted.iloc[i,0]
        
        segs_act_ind = bool2ind(a)  
        segs_act = find_uninterruptedseq(segs_act_ind)            
        segs_act_list.append(segs_act)        
        
        segs_act_ind1 = bool2ind(a, 0)
        segs_act1 = find_uninterruptedseq(segs_act_ind1)            
        segs_act_list1.append(segs_act1)

        #segs_act_list.append(find_uninterruptedseq2(a, winj))    
    matNi = []; matsovo = []
    matova = []; matovb = []; matovbt = []
    matovar = []; matovbr = []

    matNi1 = []; matsovo1 = []
    matova1 = []; matovb1 = []; matovbt1 = []
    matovar1 = []; matovbr1 = []
    #pred_new1 = []
    #pred_new2 = []
    for i in range(predicted.shape[0]):    
        segs_act = segs_act_list[i]
        segs_act1 = segs_act_list1[i]
        Nilist = []; sovolist = []
        ovalist = []; ovblist = []; ovbtlist = []
        ovarlist = []; ovbrlist = [];

        Nilist1 = []; sovolist1 = []
        ovalist1 = []; ovblist1 = []; ovbtlist1 = []
        ovarlist1 = []; ovbrlist1 = [];

        for j in range(predicted.shape[1]):
            a = predicted.iloc[i,j]            
            b = bool2ind(a)
            segs_pred = find_uninterruptedseq(b, )             
            temp = []
            for x in segs_pred:
                temp.extend(x)                
            temp = np.array(temp)
#            if j == 1 :
#                pred_new1.append(ind2bool( temp, len(a)))
#            if j == 3 :
#                pred_new2.append(ind2bool( temp, len(a)))
            
            if len(segs_pred) > 0:
                segs_pred_bool = ind2bool(np.concatenate( segs_pred ).tolist(), len(a))
                segs_pred_ind1 = bool2ind(np.array(segs_pred_bool), 0)                            
            else:
                segs_pred_ind1 = bool2ind(a, 0)                                        
            segs_pred1 = find_uninterruptedseq(segs_pred_ind1)
            
            Ni, SOVo, ova, ovb, ovbt, ovar, ovbr = SOVraw(segs_act, segs_pred, minmatch[0])
            Ni1, SOVo1, ova1, ovb1, ovbt1, ovar1, ovbr1 = SOVraw(segs_act1, segs_pred1, minmatch[1])
            
            Nilist.append(Ni)
            sovolist.append(SOVo)
            ovalist.append(ova)
            ovblist.append(ovb)
            ovbtlist.append(ovbt)
            ovarlist.append(ovar)
            ovbrlist.append(ovbr)            
            
            Nilist1.append(Ni1)
            sovolist1.append(SOVo1)
            ovalist1.append(ova1)
            ovblist1.append(ovb1)
            ovbtlist1.append(ovbt1)
            ovarlist1.append(ovar1)
            ovbrlist1.append(ovbr1) 
            
        matNi.append(Nilist)
        matsovo.append(sovolist)
        matova.append(ovalist)
        matovb.append(ovblist)
        matovbt.append(ovbtlist)
        matovar.append(ovarlist)
        matovbr.append(ovbrlist)
        
        matNi1.append(Nilist1)
        matsovo1.append(sovolist1)
        matova1.append(ovalist1)
        matovb1.append(ovblist1)
        matovbt1.append(ovbtlist1)
        matovar1.append(ovarlist1)
        matovbr1.append(ovbrlist1)
    return np.array(matNi), np.array(matsovo), np.array(matova), np.array(matovb), np.array(matovbt), np.array(matovar), np.array(matovbr), np.array(matNi1), np.array(matsovo1), np.array(matova1), np.array(matovb1), np.array(matovbt1), np.array(matovar1), np.array(matovbr1)

def process_calSOV_corr(predicted, flag_mean=True, minmatch=[1,1]):
    Na, SOVao, ovaa, ovba, ovbta, ovara, ovbra, Nn, SOVno, ovan, ovbn, ovbtan, ovarn, ovbrn = calSOV_corr(predicted)
    #Nn, SOVno, ovan, ovbn, ovbtan, ovarn, ovbrn = calSOV(predicted_ind2, winj)
    No = SOVao + SOVno
    nSOVa = SOVao/Na
    nSOVn = SOVno/Nn
    nSOVo = (SOVao + SOVno)/(Na + Nn) 
            
    result = pd.DataFrame( index=predicted.columns[0:])
    
    if flag_mean:
        result['SOVa'] = np.nanmean(nSOVa, axis = 0).T * 100
        result['SOVn'] = np.nanmean(nSOVn, axis = 0).T * 100
        result['SOVo'] = np.nanmean(nSOVo, axis = 0).T * 100
    else:
        result['SOVa'] = np.nanmedian(nSOVa, axis = 0).T * 100
        result['SOVn'] = np.nanmedian(nSOVn, axis = 0).T * 100
        result['SOVo'] = np.nanmedian(nSOVo, axis = 0).T * 100
    
    result['SOVav'] = (result.SOVa + result.SOVn)/2
    result['TPR_seg'] = ovaa.sum(axis=0)/ovaa[:,0].sum(axis=0)  * 100
    result['TNR_seg'] = ovan.sum(axis=0)/ovan[:,0].sum(axis=0)  * 100
    #result['Preci_seg_a'] = ovba.sum(axis=0)/ovbta.sum(axis=0)  * 100
    #result['Preci_seg_b'] = ovbn.sum(axis=0)/ovbtan.sum(axis=0)  * 100
    P = ovara[:,0].sum(axis=0)
    N = ovarn[:,0].sum(axis=0)
    TP = np.nansum(ovara, axis=0)
    TN = np.nansum(ovarn, axis=0)
    FN = P - TP
    FP = N - TN
    tpr_prot = ovara.T/ovara[:,0].T
    tnr_prot = ovarn.T/ovarn[:,0].T
    tpr_prot[np.isnan(tpr_prot)] = 0
    tnr_prot[np.isnan(tnr_prot)] = 0
    ACC2 = (tpr_prot  + tnr_prot ).T/2
    ACC2 = pd.DataFrame(ACC2, columns=predicted.columns[0:])
    
    result['TPR'] = ovara.sum(axis=0)/P * 100
    result['TNR'] = ovarn.sum(axis=0)/N * 100
    result['Precision'] = ovara.sum(axis=0)/ovbra.sum(axis=0) * 100
    result['ACC'] = (TP + TN)/(TP + TN + FP + FN) * 100
    result['ACC2'] = (result['TPR'] + result['TNR'])/2 
    result['FPR'] = (ovbra.sum(axis=0) - ovara.sum(axis=0))/ovbra[:,0].sum(axis=0) * 100
    result = result.round(2)
    result['MCC'] = ((TP * TN) - (FP * FN))/np.sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)  )
    result['MCC'] = result['MCC'].round(3)
    result['TP'] = TP
    result['TN'] = TN
    result['FP'] = FP
    result['FN'] = FN
    result['total'] = TP + TN + FP + FN      
    #result['P'] = P
    #result['N'] = N
    return result, ACC2

def calSOV_corr1(predicted, winj):
    #winj = 4
    segs_act_list = []
    segs_act_list1 = []
    for i in range(predicted.shape[0]):
        a = predicted.iloc[i,0]
        
        segs_act_ind = bool2ind(a)  
        segs_act = find_uninterruptedseq(segs_act_ind)            
        segs_act_list.append(segs_act)        
        
        segs_act_ind1 = bool2ind(a, 0)
        segs_act1 = find_uninterruptedseq(segs_act_ind1)            
        segs_act_list1.append(segs_act1)

        #segs_act_list.append(find_uninterruptedseq2(a, winj))    
    matNi = []; matsovo = []
    matova = []; matovb = []; matovbt = []
    matovar = []; matovbr = []

    matNi1 = []; matsovo1 = []
    matova1 = []; matovb1 = []; matovbt1 = []
    matovar1 = []; matovbr1 = []
    
    for i in range(predicted.shape[0]):    
        segs_act = segs_act_list[i]
        segs_act1 = segs_act_list1[i]
        Nilist = []; sovolist = []
        ovalist = []; ovblist = []; ovbtlist = []
        ovarlist = []; ovbrlist = [];

        Nilist1 = []; sovolist1 = []
        ovalist1 = []; ovblist1 = []; ovbtlist1 = []
        ovarlist1 = []; ovbrlist1 = [];

        for j in range(predicted.shape[1]):
            a = predicted.iloc[i,j]
            b = bool2ind(a)  
            segs_pred = find_uninterruptedseq2(b, winj[j]) 
            if winj[j] > 0 and len(segs_pred) > 0:
                segs_pred_bool = ind2bool(np.concatenate( segs_pred ).tolist(), len(a))
                segs_pred_ind1 = bool2ind(np.array(segs_pred_bool), 0)                            
            else:
                segs_pred_ind1 = bool2ind(a, 0)                                        
            segs_pred1 = find_uninterruptedseq(segs_pred_ind1)
            
            Ni, SOVo, ova, ovb, ovbt, ovar, ovbr = SOVraw(segs_act, segs_pred)
            Ni1, SOVo1, ova1, ovb1, ovbt1, ovar1, ovbr1 = SOVraw(segs_act1, segs_pred1)
            
            Nilist.append(Ni)
            sovolist.append(SOVo)
            ovalist.append(ova)
            ovblist.append(ovb)
            ovbtlist.append(ovbt)
            ovarlist.append(ovar)
            ovbrlist.append(ovbr)            
            
            Nilist1.append(Ni1)
            sovolist1.append(SOVo1)
            ovalist1.append(ova1)
            ovblist1.append(ovb1)
            ovbtlist1.append(ovbt1)
            ovarlist1.append(ovar1)
            ovbrlist1.append(ovbr1) 
            
        matNi.append(Nilist)
        matsovo.append(sovolist)
        matova.append(ovalist)
        matovb.append(ovblist)
        matovbt.append(ovbtlist)
        matovar.append(ovarlist)
        matovbr.append(ovbrlist)
        
        matNi1.append(Nilist1)
        matsovo1.append(sovolist1)
        matova1.append(ovalist1)
        matovb1.append(ovblist1)
        matovbt1.append(ovbtlist1)
        matovar1.append(ovarlist1)
        matovbr1.append(ovbrlist1)
    return np.array(matNi), np.array(matsovo), np.array(matova), np.array(matovb), np.array(matovbt), np.array(matovar), np.array(matovbr), np.array(matNi1), np.array(matsovo1), np.array(matova1), np.array(matovb1), np.array(matovbt1), np.array(matovar1), np.array(matovbr1)

def SOVraw(segs_act, segs_pred, minmatch=1, delflag=True):
    SOVi = []
    SOVo = []
    N = 0
    overlapa = np.zeros(len(segs_act))
    overlapb = np.zeros(len(segs_pred))
    overlapar = np.zeros(len(segs_pred))
    overlapbr = np.zeros(len(segs_pred))
    for a in range(len(segs_act)):
        s1 = segs_act[a]
        l1 = len(s1)
       # if l1 > 20:
        #    continue
        for b in range(len(segs_pred)):
            s2 = segs_pred[b]
            l2 = len(s2)
            l3 = len(set(s1 + s2))
            l4 = l1 + l2 - l3
            if l4 > 0:
                N = N + l1
                if l4 >= minmatch:
                    overlapa[a] = 1
                    overlapb[b] = 1
                overlapar[b] = l4 
                overlapbr[b] = l2
                if delflag :
                    ldel = min( [ l3-l4, l4, int(l1/2), int(l2/2) ] )
                    SOVi.append( (l4 + ldel)/l3 )
                    SOVo.append( l1 * (l4 + ldel)/l3 )    
                else:
                    SOVi.append( (l4)/l3 )
                    SOVo.append( l1 * (l4 )/l3 )    
        if overlapa[a] == 0:
                N = N + l1
    return N, np.sum(SOVo), overlapa.sum(), overlapb.sum(), len(segs_pred), overlapar.sum(), overlapbr.sum()
