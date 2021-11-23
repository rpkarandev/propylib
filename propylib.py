#!/user/bin/python
## load core library
import glob
import math
import numpy as np
import pandas as pd

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
aapair = np.array([ [y+x for x in aa1code] for y in aa1code]).flatten()
aaA2AAA = {'A' :'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER', 'T':'THR', 'V': 'VAL', 'W' : 'TRP', 'Y' : 'TYR', 'X': 'UNK' } #, 'B': 'MSE', 'J' : 'PCA', 'O': 'SL5', 'J': 'SWG'
aaAAA2A = {'ALA':'A', 'CYS':'C', 'ASP' : 'D',  'GLU' : 'E' , 'PHE' : 'F', 'GLY' : 'G', 'HIS' : 'H', 'ILE' : 'I', 'LYS': 'K' , 'LEU':'L', 'MET':'M' , 'ASN':'N', 'PRO':'P' , 'GLN':'Q' , 'ARG':'R' , 'SER':'S' , 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y', 'UNK' :'X'} #, 'MSE':'B', 'PCA' : 'B', 'SL5' : 'B', 'SWG' : 'B' , 'TPO' : 'B', 'MIR' : 'B', 'PTR' : 'B', 'PIA' : 'B', 'CRF' : 'B', 'CZZ' : 'B'  

aa2chph74 = {'A' :0, 'C':0, 'D':-1, 'E':-1, 'F':0, 'G':0, 'H':0.5, 'I':0, 'K':1, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':1, 'S':0, 'T':0, 'V':0, 'W' : 0, 'Y' :0, 'X':0} 
aa2chph26 = {'A' :0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':1, 'I':0, 'K':1, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':1, 'S':0, 'T':0, 'V':0, 'W' : 0, 'Y' :0, 'X':0} 

dssp2ss1 = {'H':'H', 'G':'H', 'I':'H', 'E':'E', 'B':'E', 'S':'T', 'T':'T', '':'C'}
dssp2ss2 = {'H':'H', 'G':'H', 'I':'H', 'E':'E', 'B':'E', 'S':'C', 'T':'C', '':'C'}
aa3code = ['ALA', 'ASP', 'CYS', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
aa1code = ['A', 'D', 'C', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
aatASA = [110.2, 144.1, 140.4, 174.7, 200.7, 78.7, 181.9, 185, 205.7, 183.1, 200.1, 146.4, 141.9, 178.6, 229, 117.2, 138.7, 153.7, 240.5, 213.7];
aaShphobic = [0.81, 0.88, 0.72, 0.46, 0.95, 0.77, 0.54, 1, 0.26, 0.92, 0.81, 0.45, 0.68, 0.43, 0, 0.6, 0.63, 0.92, 0.85, 0.71] 
aaShphobic_c = [x-0.77 for x in aaShphobic]
aa3toaa1 = { a:[b,c,d] for a,b,c,d in zip(aa3code,aa1code,aatASA,aaShphobic)}
stride2ss = {'H':'H', 'E':'E', 'T':'C', 'C':'C', 'G':'H', 'I':'H', 'b':'E', 'B':'E'} 
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

###########################################################################################




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
