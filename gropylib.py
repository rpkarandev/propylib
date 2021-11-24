r#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 10:10:05 2019
functions to handle gromacs output files
@author: karan
"""

import re
from glob import glob
import numpy  as np
import pandas as pd

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
      
def read_xvg(filename):
    data = []
    cols = []
    comment = []
    ylabel = []
    with open(filename, 'r') as infile:
        line = infile.readline()        
        while line:            
            if line.startswith('@'):
                line_st = line.split()
                if line_st[1] == 'xaxis':
                    m = re.search('"(.*)"', line) 
                    if m is not None:
                        cols.append( m.groups()[0] )                
                elif line_st[1] == 'yaxis':
                    m = re.search('"(.*)"', line)
                    if m is not None:
                        ylabel.append( m.groups()[0] )                
                elif line_st[1][0] == 's' and line_st[2] == 'legend' :
                    m = re.search('legend\s+"(.*)"', line) 
                    if m is not None:
                        cols.append( m.groups()[0] ) 
            elif line.startswith('#'):
                comment.append( line.strip() )
            else:
            #elif line.startswith(' '):
                #print(line)
                data.append( line.split() )
            line = infile.readline()
    for i in range(len(data[0])-len(cols)):
        cols = cols + [str(i)]
    data = pd.DataFrame( data, columns=cols ).astype('float')
    return data, comment, ylabel

def read_xvg_notnum(filename):
    data = []
    cols = []
    comment = []
    ylabel = []
    with open(filename, 'r') as infile:
        line = infile.readline()        
        while line:            
            if line.startswith('@'):
                line_st = line.split()
                if line_st[1] == 'xaxis':
                    m = re.search('"(.*)"', line) 
                    if m is not None:
                        cols.append( m.groups()[0] )                
                elif line_st[1] == 'yaxis':
                    m = re.search('"(.*)"', line)
                    if m is not None:
                        ylabel.append( m.groups()[0] )                
                elif line_st[1][0] == 's' and line_st[2] == 'legend' :
                    m = re.search('legend\s+"(.*)"', line) 
                    if m is not None:
                        cols.append( m.groups()[0] ) 
            elif line.startswith('#'):
                comment.append( line.strip() )
            else:
            #elif line.startswith(' '):
                #print(line)
                data.append( line.split() )
            line = infile.readline()
    for i in range(len(data[0])-len(cols)):
        cols = cols + [str(i)]
    data = pd.DataFrame( data, columns=cols )#.astype('float')
    return data, comment, ylabel

def read_rawxpm(filename):
    data = []    
    comment = []
    xlabel = ''
    ylabel = ''
    xaxis = []
    yaxis=[]
    strings = []
    legenddict = {}
    valdict = {}   
    with open(filename, 'r') as infile:
        line = infile.readline()
        while line:
            #print(line)
            if xlabel == '':
                m = re.search('/\* x-label: "(.*)" \*/', line)
                if m != None :
                    xlabel = m.groups()[0] 
            elif ylabel == '':
                m = re.search('/\* y-label: "(.*)" \*/', line)
                if m != None :
                    ylabel = m.groups()[0] 
            elif line[0] == '"' and len(xaxis) == 0 :            
                m = re.search('"(.*?)\s+[c]*\s+(.*)\s+" /\*\s"(.*)"\s\*/', line)
                if m != None:
                    string = m.groups()[0] 
                    legenddict[ m.groups()[0] ] = m.groups()[2] 
                    valdict[ m.groups()[0] ] = m.groups()[1]
            elif line[:10] == "/* x-axis:":
                m = re.search('/\* x-axis:\s+(.*) \*/', line)
                xaxis.extend( m.groups()[0].strip().split() )
            elif line[:10] == "/* y-axis:":
                m = re.search('/\* y-axis:\s+(.*) \*/', line)
                yaxis.extend( m.groups()[0].strip().split() )
            elif line[0] == '"' and len(xaxis) > 0 :            
                m = re.search('"(.*)"', line)
                string = m.groups()[0] 
                strings.append( list(string) )
            line = infile.readline()
        data = pd.DataFrame(strings, ).T
        data.columns = [ylabel + '_' + y  for y in yaxis]        
        data[xlabel] = pd.Series(xaxis).astype('float')
        data = data.loc[:, [ data.columns[-1] ] + data.columns[:-1].values.tolist() ]
        return data, comment, legenddict, valdict, ylabel

def writeNDXfile(grpIndexDic, indexFileName):
    with open(indexFileName, 'w') as ofile:
        s = ''
        for key in grpIndexDic.keys():
            s += f'[ {key} ]\n'
            ld = 10
            for ai, a in enumerate( grpIndexDic[key] ):
                if ai%ld == 0:
                    s += '\n'
                s += f'{a:5d} '
            s += '\n\n'
        ofile.write(s)
    return None

def read_numxpm(filename): 
    data, comment, legenddict, valdict, ylabel = read_rawxpm(filename)
    if len( data ) > 0:
        data.iloc[:,1:] = data.iloc[:,1:].apply(lambda x: x.apply( lambda y: legenddict[y] ) )
        return data, comment, legenddict, valdict, ylabel
    else:
        print("error reading xpm")
        return None
        
def read_ssxpm(filename):
    data, comment, legenddict, valdict, ylabel = read_rawxpm(filename)
    if len( data ) > 0:
        cols = []
        chsep = []
        try:
            chsep = [x for x in legenddict.keys() if legenddict[x] == "Chain_Separator"]
        except:
            print("No chainseparators are found")
        chno = 1
        resno = 1    
        for x in data.iloc[0, 1:]: # first row, assign chain based on that
            if x in chsep:
                chno += 1
                resno = 1
                cols.append( "chsep" )
            else:
                cols.append( "ch"+str(chno)+"_res"+str(resno) )
                resno += 1
        data.columns = [data.columns[0]] + cols
        data = data.loc[:, data.columns != 'chsep']        
        return data, comment, legenddict, valdict, ylabel
    else:
        print("error reading xpm")
        return None

def read_gro(filename):
    try :
        with open(filename, 'r') as infile:
            data = []
            line = infile.readline()
            title = line
            m = re.search('t=([\d.]*)', line)
            if m != None:
                frtime = float(m.groups()[0])
            else:
                frtime = 0.0
            line = infile.readline()    
            m = re.search('t=([\d.]*)', line)
            tnmol = int(line.strip())
            line = infile.readline()
            line1 = infile.readline()
            while line1 :
                data.append( [ line[:5].strip(), line[5:10].strip(), line[10:15].strip(), line[15:20].strip(), line[20:28].strip(), line[28:36].strip(), line[36:44].strip(), line[44:52].strip(), line[52:60].strip(), line[60:68].strip() ] )
                line = line1
                line1 = infile.readline()
            box = line.strip().split()
            box = [float(x) for x in box]
            cols =['resno', 'resname', 'aname', 'ano', 'x', 'y', 'z', 'vx', 'vy', 'vz']
            dtypes = ['int64', 'object', 'object', 'int64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64']
            dtypes = ['int', 'object', 'object', 'int', 'float', 'float', 'float', 'float', 'float', 'float']
            dtypesdict = dict( zip( cols, dtypes) )
            data = pd.DataFrame( data, columns=cols) 
            data = data.astype( dtypesdict )
            return data, title, frtime, box
    except:
        print("read_gro: error -> unable to open the file")

class read_pdbtraj:
    'Class for storing PDB data'
    def __init__(self, filename, framestart=0, frameend=0, chainids=[], nresperch=0, verbose=1):
        self.fname = ((filename.split('/'))[-1].split('.'))[0]
        if type(chainids) != list :
            print("read_pdbtraj: Error chainids must be a list")
            return -100
        elif type(nresperch) != int or nresperch < 0 :
            print("read_pdbtraj: Enter valid nresperch ")
            return -100
        elif type(framestart) != int or framestart < 0 :
            print("read_pdbtraj: Enter valid framestart ")
            return -100
        elif type(frameend) != int or frameend < 0 :
            print("read_pdbtraj: Enter valid framestart ")
            return -100
        frameend = 10000000 if frameend < 1 else frameend
        nchs_to_get = len(chainids)            
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
        atom_recs_l = []
        atom_coords_l = []
        atom_coords_list = []
        hetatom_recs_l = []
        hetatom_coords_l = []
        hetatom_coords_list = []
        readmodelcounter = 0
        framelist = []
        timelist = []
        boxlist = []
        
        print( " framend ", frameend) if verbose > 0 else None
        with open(filename, 'r') as pdbfile :                        
            line = pdbfile.readline()			
            while line and modelcounter <= frameend:
                #print( line )
                if line[:5] == 'TITLE': 
                    #print(line)
                    #m = re.search('t=\s+([0-9]*\.[0-9]*)\s+', line)
                    m = re.search('t=\s+([0-9\.]*)\s+', line)
                    if m != None:
                        #print(m.groups()[0].strip())
                        frtime = float(m.groups()[0].strip())
                    else:
                        frtime = 0.5
                if line[:6] == 'CRYST1': 
                    #print(line)
                    m = re.search('CRYST1\s+([0-9\.]*)\s+([0-9\.]*)\s+([0-9\.]*)\s+', line)
                    if m != None:                        
                        #print(m.groups()[0].strip())
                        box = [ float(m.groups()[0].strip()), float(m.groups()[0].strip()), float(m.groups()[0].strip())]
                    else:
                        box = [0.0, 0.0, 0.0]
                    #print(box)
                if line[:5] == 'MODEL':                         
                    modelcounter += 1                    
                    if modelcounter >= framestart and modelcounter<= frameend:
                        modelflag = True
                        print("read_pdbtraj: reading model ", modelcounter)  if verbose > 0 else None
                        readmodelcounter += 1
                        framelist.append(modelcounter)
                        timelist.append( frtime )
                        boxlist.append( box )
                        if len( atom_coords_l ) > 0 or len( hetatom_coords_l ) > 0 :
                            atom_coords_list.append( atom_coords_l )
                            atom_coords_l = []
                            hetatom_coords_list.append( hetatom_coords_l )
                            hetatom_coords_l = []
                    else:
                        modelflag = False                        
                if readmodelcounter == 1:      
                    if modelflag and nchs_to_get == 0: ## Read atom coordinates
                        if line[:6] == 'ATOM  ' :
                            atom_recs_l.append( [int(line[6:11].strip()), line[12:16].strip(), line[17:20].strip(), line[21:22].strip(), int(line[22:26].strip()), float(line[54:60].strip()), float(line[60:66].strip()), line[76:78].strip() ]  )   
                            #atom_coords_l.append([modelcounter, frtime, float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] )                                  
                            atom_coords_l.append([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] )                                  
                        if line[:6] == 'HEATOM':
                            hetatom_recs_l.append( [int(line[6:11].strip()), line[12:16].strip(), line[17:20].strip(), line[21:22].strip(), int(line[22:26].strip()), float(line[54:60].strip()), float(line[60:66].strip()), line[76:78].strip() ]  )   
                            #hetatom_coords_l.append([modelcounter, frtime, float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] )                                  
                            hetatom_coords_l.append([ float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] )                                  
                        #line = pdbfile.readline()
                    if modelflag and nchs_to_get > 0: ## Read atom coordinates from chain
                        if line[:6] == 'ATOM  ' and line[20:22].strip() in chainids :                        
                            atom_recs_l.append( [int(line[6:11].strip()), line[12:16].strip(), line[17:20].strip(), line[21:22].strip(), int(line[22:26].strip()), float(line[54:60].strip()), float(line[60:66].strip()), line[76:78].strip() ]  )   
                            #atom_coords_l.append([modelcounter, frtime, float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] )
                            atom_coords_l.append([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] )                                                                    
                        if line[:6] == 'HEATOM':
                            hetatom_recs_l.append( [int(line[6:11].strip()), line[12:16].strip(), line[17:20].strip(), line[21:22].strip(), int(line[22:26].strip()), float(line[54:60].strip()), float(line[60:66].strip()), line[76:78].strip() ]  )   
                            #hetatom_coords_l.append([modelcounter, frtime, float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] )
                            hetatom_coords_l.append([ float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ])
                        #line = pdbfile.readline()
                else:
                    if modelflag and nchs_to_get == 0: ## Read atom coordinates
                        if line[:6] == 'ATOM  ' :
                            #atom_coords_l.append([modelcounter, frtime, float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] )                                  
                            atom_coords_l.append([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] )                                  
                        if line[:6] == 'HEATOM':
                            #hetatom_coords_l.append([modelcounter, frtime, float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] )                                  
                            hetatom_coords_l.append([ float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] )                                  
                        #line = pdbfile.readline()
                    if modelflag and nchs_to_get > 0: ## Read atom coordinates from chain
                        if line[:6] == 'ATOM  ' and line[20:22].strip() in chainids :                        
                            #atom_coords_l.append([modelcounter, frtime, float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] )
                            atom_coords_l.append([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] )                                                                    
                        if line[:6] == 'HEATOM':
                            #hetatom_coords_l.append([modelcounter, frtime, float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] )
                            hetatom_coords_l.append([ float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ])
                        #line = pdbfile.readline()
                line = pdbfile.readline()
                #print(line)
        cols1 = ['atomno', 'atomname', 'resname', 'chain', 'resno', 'occupancy', 'bfac', 'atomtype']
        #cols2 = ['frame', 'time', 'X', 'Y', 'Z']
        self.atom_recs = pd.DataFrame( atom_recs_l, columns=cols1 )
        self.hetatom_recs = pd.DataFrame( hetatom_recs_l, columns=cols1 )
        #self.atom_coords = pd.DataFrame( atom_coords_l, columns=cols2 )
        #self.hetatom_coords = pd.DataFrame( hetatom_coords_l, columns=cols2 )
        self.atom_coords = np.array( atom_coords_list )
        self.hetatom_coords = pd.DataFrame( hetatom_coords_list )
        self.framelist = pd.DataFrame(  )
        self.framelist['frameno'] = framelist
        self.framelist['time'] = timelist
        self.framelist['boxxyz'] = boxlist
        
        
        if nresperch > 0:            
            resnoc_list = [] # continuous resno starts from 0
            resnoc = -1
            resno_old = -100
            for resno in self.atom_recs.resno:
                if resno_old != resno:
                    resnoc += 1
                    resnoc_list.append( resnoc )
                    resno_old = resno
                else:
                    resnoc_list.append( resnoc )
            self.atom_recs['resnoc'] = resnoc_list
            self.atom_recs['chainidc'] = self.atom_recs.resnoc.apply(lambda x: int( np.floor(x/nresperch) ) )
            self.atom_recs['resnoch'] = self.atom_recs.resnoc.apply(lambda x: int( x%nresperch ) )
            self.atom_recs['m'] = self.atom_recs.atomtype.apply( assign_mass )
#        self.resnos = [ [ self.atom.resno[y] for y in x ] for x in self.indCalpha]
#        self.resnames = [ [ self.atom.resname[y] for y in x ] for x in self.indCalpha]
#        self.resnames1 = [ [ aminoAAA2A(y) for y in x ] for x in self.resnames]
#        self.indBB = [[x for x in range(len(self.atom.atomname)) if (self.atom.atomname[x] == 'N') or (self.atom.atomname[x] == 'CA') or (self.atom.atomname[x] == 'C')] for y in self.chains]    
#        self.indnH = [[x for x in range(len(self.atom.atomname)) if (self.atom.atomname[x] != 'H') or (self.atom.atomname[x] != 'D') ] for y in self.chains]    
def cal_principalaxes(chcoordsfr, chm):
    com = ( chcoordsfr * chm ).sum(axis=0)/chm.sum()
    chcoordsfr = chcoordsfr - com
    chcoordsfr2 = chcoordsfr * chcoordsfr
    I11 = (chm * (chcoordsfr2[:, 1] + chcoordsfr2[:, 2])).sum()
    I22 = (chm * (chcoordsfr2[:, 0] + chcoordsfr2[:, 2])).sum()
    I33 = (chm * (chcoordsfr2[:, 0] + chcoordsfr2[:, 1])).sum()
    I12 = -(chm * chcoordsfr[:, 0] * chcoordsfr[:, 1] ).sum()
    I13 = -(chm * chcoordsfr[:, 0] * chcoordsfr[:, 2] ).sum()
    I23 = -(chm * chcoordsfr[:, 1] * chcoordsfr[:, 2] ).sum()
    I = np.array( [I11, I12, I13, I12, I22, I23, I13, I23, I33]).reshape( (3,3) )
    eigval, eigvec = np.linalg.eig( I )
    return eigvec, com, eigval
def cal_largestaxis(chcoordsfr, chm):
    #n = 1000
    # chcoordsfr = np.random.normal(size=3*n).reshape( (n, 3) ) +  100
    #chm = np.ones( (n,1) )                                     
    eigvec, com, eigval = cal_principalaxes(chcoordsfr, chm)
    eigvalmaxind = np.argmax(np.abs(eigval))
    eigvecmax = eigvec[:, eigvalmaxind ]
    return eigvecmax, com
def get_neighbours(clusters_temp, slist):
        temp = []
        for s in slist:
            temp += clusters_temp[s]
        return list( set(temp) )
    
def get_clusters( distmat0, distcutoff ):    
    distmat = distmat0 <= distcutoff
    cluster_rows = [x for x in range(len(distmat))]
    clusters_temp = [np.where( distmat[ci, :] )[0].tolist() for ci in cluster_rows]    
    clusters_temp1 = []      
    while len(cluster_rows) > 0:
        ci = cluster_rows.pop( 0 )
        clust = clusters_temp[ ci ]
        nclust = get_neighbours(clusters_temp, clust)
        while len(clust) != len(nclust):
            clust = nclust.copy()
            nclust = get_neighbours(clusters_temp, clust)
        cluster_rows = [x for x in cluster_rows if x not in nclust]
        clusters_temp1.append( nclust )
    clus_centers = []
    for clus in clusters_temp1:
        ind_min = np.argmin(np.mean( distmat0[ clus, : ][:, clus], axis=1))
        clus_centers.append( clus[ ind_min ] )
    return clusters_temp1, clus_centers

def get_all_clusters(pairdist, distcutoff  ):
    nch1 = int( np.sqrt( pairdist.shape[1]-1 ) )
    distmats = pairdist.iloc[:,1:].values        
    clus_stat = []
    sizecounts = {}
    ncounts = {}
    for x in range( distmats.shape[0] ) :
        sizecount = {'1':0, '2':0, '3-4':0, '5-6':0, '7-8':0, '9+':0, '11+':0,'13+':0}
        ncount = {'1':0, '2':0, '3-4':0, '5-6':0, '7-8':0, '9+':0, '11+':0,'13+':0}
        distmat0 = distmats[x,:].reshape( (nch1, nch1) ) 
        clusters, clus_centers = get_clusters( distmat0, distcutoff )
        clus_lens = [ len(x) for x in clusters ]
        clus_stat.append( [len(clusters), np.max(clus_lens), np.mean(clus_lens) ])
        for cl in clus_lens:
            if cl == 1:
                sizecount['1'] += cl
                ncount['1'] += 1
            elif cl == 2:
                sizecount['2'] += cl
                ncount['2'] += 1
            elif cl in [3,4]:
                sizecount['3-4'] += cl
                ncount['3-4'] += 1
            elif cl in [5,6] :
                sizecount['5-6'] += cl
                ncount['5-6'] += 1
            elif cl in [7,8] :
                sizecount['7-8'] += cl
                ncount['7-8'] += 1
            elif cl > 8 :
                sizecount['9+'] += cl
                ncount['9+'] += 1
            if cl > 10 :
                sizecount['11+'] += cl
                ncount['11+'] += 1
            if cl > 12 :
                sizecount['13+'] += cl
                ncount['13+'] += 1
        sizecounts[x] = sizecount 
        ncounts[x] = ncount
    return np.array( clus_stat ),  pd.DataFrame( sizecounts ).T, pd.DataFrame( ncounts ).T   
def get_all_clusters_update1(pairdist, distcutoff  ):
    nch1 = int( np.sqrt( pairdist.shape[1]-1 ) )
    distmats = pairdist.iloc[:,1:].values        
    clus_stat = []
    sizecounts = {}
    ncounts = {}
    for x in range( distmats.shape[0] ) :
        sizecount = {'1':0, '2':0, '3-8':0, '9-16':0, '17-24':0, '24+':0}
        ncount = {'1':0, '2':0, '3-8':0, '9-16':0, '17-24':0, '24+':0}
        distmat0 = distmats[x,:].reshape( (nch1, nch1) ) 
        clusters, clus_centers = get_clusters( distmat0, distcutoff )
        clus_lens = [ len(x) for x in clusters ]
        clus_stat.append( [len(clusters), np.max(clus_lens), np.mean(clus_lens) ])
        for cl in clus_lens:
            if cl == 1:
                sizecount['1'] += cl
                ncount['1'] += 1
            elif cl == 2:
                sizecount['2'] += cl
                ncount['2'] += 1
            elif cl < 9:
                sizecount['3-8'] += cl
                ncount['3-8'] += 1
            elif cl < 17 :
                sizecount['9-16'] += cl
                ncount['9-16'] += 1
            elif cl < 25 :
                sizecount['17-24'] += cl
                ncount['17-24'] += 1
            elif cl >= 25 :
                sizecount['24+'] += cl
                ncount['24+'] += 1
        sizecounts[x] = sizecount 
        ncounts[x] = ncount
    return np.array( clus_stat ),  pd.DataFrame( sizecounts ).T, pd.DataFrame( ncounts ).T   

def get_all_clusters_detail(pairdist, distcutoff, frameno=-1  ):
    nch1 = int( np.sqrt( pairdist.shape[1]-1 ) )
    distmats = pairdist.iloc[:,1:].values        
    clus_stat = []
    sizecount = {'1':0, '2':0, '3-4':0, '5-6':0, '7-8':0, '9+':0, '11+':0,'13+':0}
    ncount = {'1':0, '2':0, '3-4':0, '5-6':0, '7-8':0, '9+':0, '11+':0,'13+':0}
    ## Take only the last frame
    distmat0 = distmats[frameno,:].reshape( (nch1, nch1) ) 
    clusters, clus_centers = get_clusters( distmat0, distcutoff )
    clus_lens = [ len(x) for x in clusters ]
    clus_stat.append( [len(clusters), np.max(clus_lens), np.mean(clus_lens) ])
    
    for cl in clus_lens:
        if cl == 1:
            sizecount['1'] += cl
            ncount['1'] += 1
        elif cl == 2:
            sizecount['2'] += cl
            ncount['2'] += 1
        elif cl in [3,4]:
            sizecount['3-4'] += cl
            ncount['3-4'] += 1
        elif cl in [5,6] :
            sizecount['5-6'] += cl
            ncount['5-6'] += 1
        elif cl in [7,8] :
            sizecount['7-8'] += cl
            ncount['7-8'] += 1
        elif cl > 8 :
            sizecount['9+'] += cl
            ncount['9+'] += 1
        if cl > 10 :
            sizecount['11+'] += cl
            ncount['11+'] += 1
        if cl > 12 :
            sizecount['13+'] += cl
            ncount['13+'] += 1                  
    return clusters, clus_centers, clus_lens, np.array( clus_stat ),  sizecount, ncount  

        
        
