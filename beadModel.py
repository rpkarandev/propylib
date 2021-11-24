import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy
import logging

import nglview as ngl
import networkx as nx


logging.basicConfig(format='%(asctime)s %(message)s')

DEBUG = True
if DEBUG == True:
    logging.basicConfig(level = logging.DEBUG )
    #logging.basicConfig( 'root':{'handlers':('console', 'file'), 'level':'DEBUG'} )
    
    
## bead properties

aaList = list( 'ACDEFGHIKLMNPQRSTVWY' ) 
aaList += ['CA', 'C', 'H', 'N', 'O', 'OH', 'S', 'NT', 'NE', 'NH', 'OT', 'CB', 'CG', 'CD', 'CE', 'CZ', 'CH'] 
#aaList = ['CA', 'C', 'H', 'N', 'O', 'S'] 
#print( len(aaList))

beadDictRadius = { a:0.77 for a in aaList }
beadDictProp1 = { a:0 for a in aaList }

ndim = 3

def getDistance( p ):
    return np.sqrt( np.sum( (p[0] - p[1]) ** 2 ) )

def getAngle( p ):
    p -= p[1]
    v1 = p[0] / np.linalg.norm(p[0])
    v2 = p[2] / np.linalg.norm(p[2])
    d = np.dot(v1, v2)
    angle = np.arccos(d)

    if np.isnan( angle ):
        return 0
    else:
        return np.degrees( angle )

def getTorsionAngle(p):
    ## https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    
    p = np.array( p )

    b0 = -1.0 * ( p[1] - p[0])
    b1 = p[2] - p[1]
    b2 = p[3] - p[2]

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    if ndim == 2:
        b1 = b1.tolist() + [0]
        v = v.tolist() + [0]
        y0 = np.cross( b1, v )[:-1]
        
    else:
        y0 = np.cross( b1, v )
    y = np.dot( y0, w)
    theta = np.degrees(np.arctan2(y, x))
    if np.isnan( theta ):
        return 0
    else:
        return theta 


class bead: # residue or atom or any smallest entity in the model
    def __init__( self, code='C', name='C', bno=0, resno=0, chain='A' ):
        self.loc = np.zeros( (ndim) ) # location or coordinates
        self.name = name # one-letter code
        self.code = code
        self.radius = beadDictRadius[ code ] # radius or coarse-grained radius
        self.prop1 = beadDictProp1[ code ] # some property - just an example
        self.interaction = 0 # count or some quantity
        self.no = bno # some interation no for referencing within chain
        self.resno = resno # some residue no for referencing within chain 
        self.chain = chain
        
    def move(self, dx ): ## update location by an increment
        self.loc += dx 
    
    def __str__(self):
        return  'bead {}({}):{} in res {} and chain {} '.format( self.name, self.code, self.no, self.resno, self.chain  ) 
        
    def __repr__(self):
        return self.__str__()

class chain: ## a molecule
    def __init__(self, seqCode=[], seq=[], chname='A', loc=np.zeros( (ndim) ), initConfor='zigzag' ):
        self.bindex = 0 # loop over beads like a iterator
        self.name = chname
        self.seqCode = []
        self.seq = []
        self.bbInd = []
        #self.seqCode = seqCode
        if len(seq) == 0:
            seq = seqCode
        else:
            if len(seq)  < len(seqCode):
                logging.warning(  'Warning bead Codes and Names differ in count {}:{}'.format( len(seqCode), len(seq) ) )
                
        self.residues = []
        self.bondPair = []
        self.initConfor = initConfor
        #self.beads = [bead(s, bi, bi, self.name) for bi, s in enumerate( seq ) ]
        if len( seqCode ) == 0:
            self.nbeads = 0
            self.beads = []
        else:
            self.nbeads = 1
            self.beads = [ bead( seqCode[0], seq[0], 0, 0, self.name) ]
            self.beads[0].loc = loc
        self.update()
        for bi in range( 1, len(seqCode) ):
            self.grow(bcode=seqCode[bi], bname=seq[bi], addtoBead=bi-1, resno=bi, conType=initConfor )
            #self.bondPair.append( (bi-1, bi) )
            #self.setConformation(bi-1, bi, conType=initConfor)
        #self.update()
        
    ## set loc/orientation second bead with respect previous
    ## based bead radiuss
    def setConformation(self, b1, b2, conType ='zigzag' ):
        if type(conType) == str:
            if  conType == 'linear':
                r1 = self.beads[b1].radius
                r2 = self.beads[b2].radius
                d = r1 + r2
                sign = 1
                self.beads[b2].loc = self.beads[b1].loc + [0, d, 0]
                
            elif  conType == 'zigzag':
                r1 = self.beads[b1].radius
                r2 = self.beads[b2].radius
                d = np.sqrt( r1 + r2 )
                sign = 1 if b2%2 == 0 else -1
                self.beads[b2].loc = self.beads[b1].loc + [d, d * sign, 0]
            
        else: ## conType must be a dict
            if 'd' in conType:
                r1 = 0
                r2 = 0
                d = np.sqrt( conType['d'] )
            else:
                r1 = self.beads[b1].radius
                r2 = self.beads[b2].radius
                d = np.sqrt( r1 + r2 )
                
            sign = 1 if b2%2 == 0 else -1
            self.beads[b2].loc = self.beads[b1].loc + [d, d * sign, 0]
            
            if False and 'ang' in conType : # and 'angleTriple' in conType:
                ang = conType['ang'][0]
                if b1 - 1 > 0 :
                    self.setAngle( [b1 - 1, b1, b2], ang  )
            
            ''' 
            if 'ang' in conType and 'angleTriple' in conType:
                for ang, angleTriple in zip( conType['ang'], conType['angleTriple'] ):
                    self.setAngle( angleTriple, ang  )
                
            if 'tor' in conType and 'torsionQuad' in conType:
                for tor, torsionQuad in zip( conType['tor'], conType['torsionQuad'] ):
                    self.setTorsionAngle( torsionQuad, tor  )
            '''  
        if  True in np.isnan( self.beads[b2].loc ):
            logging.error( ' NaN value in location, something went wrong while assiging coordinates {}{}:{}'.format(r1, r2, d) )
            logging.error(   '{}, {} '.format( d,  sign ) )
        return True
    
    def get_coords(self):
        return np.asarray( [b.loc for b in self.beads  ])
    
    ## assumming bondPair is in order
    def get_angleTriple(self):
        angleTriple = []
        n = len(self.bondPair)
        for i in range(n-1):
            a = self.bondPair[i]
            for j in range( i+1, n):
                b = self.bondPair[j]
                if a[1] == b[0]:
                    angleTriple.append( a + ( b[-1], )  )
        return angleTriple
    
    ## assuming angleTriple is in order
    def get_torsionQuad(self):
        torsionQuad = []
        n = len(self.angleTriple)
        for i in range(n-1):
            a = self.angleTriple[i]
            for j in range( i+1, n):
                b = self.angleTriple[j]
                if a[1:] == b[:-1] :
                    torsionQuad.append( a + ( b[-1], ) )
        return torsionQuad

    def update(self): ## topology change
        self.update0()
        self.update1()
        
    def update0(self): ## topology change
        self.angleTriple = self.get_angleTriple()
        self.torsionQuad = self.get_torsionQuad()
    
    def update1(self): ## conformation change
        self.coords = self.get_coords()
        self.torsions = [ getTorsionAngle( self.coords[q,:] ) for q in self.torsionQuad ]
        self.update2()
        
    def update2(self):
        self.bonds = [ getDistance( self.coords[bp,:]) for bp in self.bondPair ]
        self.angles = [ getAngle( self.coords[t,:]) for t in self.angleTriple ]
            
        
    def grow(self, bcode='C', bname='C', addtoBead=-1, resno=0, conType='zigzag' ):
        bi = self.nbeads
        if addtoBead == -1:
            addtoBead = bi-1
        b = bead(bcode, bname, bi, resno, self.name)
        self.beads.append( b )
        self.nbeads += 1
        self.seq.append( bname )
        self.seqCode.append( bcode )
        self.update()
        if bi != 0:
            self.setConformation(addtoBead, bi, conType)
            self.bondPair.append( (addtoBead, bi) )
        self.update()
        
    def __str__(self):
        return 'chain {} with {} bead(s)'.format(self.name, self.nbeads) 
        
    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        yield from self.beads

    def __next__(self):
        if self.bindex == 0:
            raise StopIteration
        self.bindex = self.bindex - 1
        return self.beads[self.bindex]
    
    def __getitem__(self, index):
        return self.beads[index]
    
    def copy(self):
        return copy.deepcopy( self )
    
    def draw(self, ax):
        pradius = 100
        coords = self.coords
        radiuss = [ b.radius * pradius for b in self ]
        for pair in self.bondPair:
            if pair[0] in self.bbInd and pair[1] in self.bbInd :
                plt.plot( coords[ pair, 0], coords[ pair, 1], 'k-'  )
            else:
                plt.plot( coords[ pair, 0], coords[ pair, 1], 'r-'  )
                
        ax.scatter( coords[:,0], coords[:, 1], radiuss, marker='o' )

    def rotateAngle(self, angleTriple, theta=10, endIndex=0, beadIndices=[] ): # one entry from angleTriple list
        startIndex = angleTriple[1] 
        if endIndex == 0:
            endIndex = self.nbeads - 1
        c = self.coords[ startIndex, :  ] # temp center
        if len(beadIndices) == 0:
            beadIndices = list( range(startIndex, endIndex) )
        coordsSub = self.coords[ beadIndices, :  ].copy() - c ## rest of the chain
        #axis = coordsSub[1]  # rotation axis
        v0 = self.coords[ angleTriple[0], :  ] - c
        v0 /= np.linalg.norm( v0 ) # unit vector
        v1 = self.coords[ angleTriple[2], :  ] - c
        v1 /= np.linalg.norm( v1 ) # unit vector
        axis = np.cross( v0, v1 ) ### considering an axis perpendicular to the two vectors
        axis /= np.linalg.norm( axis ) 
        ##
        ## alternate axis 
        v1s = [ [0, 0, 1], [0, 1, 0], [1, 0, 0], [1, 1, 1] ]
        v1i = 0
        while True in np.isnan(axis): ## two parallel lines in an axis
            axis = np.cross( v0, v1s[v1i] )
            axis /= np.linalg.norm( axis ) 
            v1i + 1
            
        rad = np.radians( theta )
        q0 = np.cos( rad/2 )
        q1 = np.sin( rad/2 ) * axis[0] 
        q2 = np.sin( rad/2 ) * axis[1]
        q3 = np.sin( rad/2 ) * axis[2]
        Q = [ [(q0**2 + q1**2 - q2**2 - q3**2),      2*(q1*q2 - q0*q3),          2*(q1*q3 + q0*q2)], 
             [ 2*(q2*q1 + q0*q3),     (q0**2 - q1**2 + q2**2 - q3**2),      2*(q2*q3 - q0*q1) ], 
             [ 2*(q3*q1 - q0*q2),          2*(q3*q2 + q0*q1),     (q0**2 - q1**2 - q2**2 + q3**2) ] ]
        Q = np.array( Q )
        newCoords = np.matmul( coordsSub, Q.T )  

        for i, bi in enumerate( beadIndices ):
            self.beads[bi].loc = newCoords[i] + c
        self.update1()
        
        if True in np.isnan( newCoords ):
            logging.error( 'Nan found in coordinates during rotation. Info below:' )
            logging.error( ['theta, rad, axis, v0, v1', theta, rad, axis, v0, v1 ] )
            logging.error( ['angleTriple, self.nbeads', angleTriple, self.nbeads ])
            logging.error( [ 'self.coords[ angleTriple ]', self.coords[ angleTriple ] ]  )
            logging.error(  ['Q', Q] )
            return False
        return True

    def rotateTorsion(self, torsionQuad, theta=10, endIndex=0, beadIndices=[] ): # one entry from torsionQuad list
        startIndex = torsionQuad[1] 
        if endIndex == 0:
            endIndex = self.nbeads - 1
        c = self.coords[ startIndex, :  ] # temp center
        if len(beadIndices) == 0:
            beadIndices = list( range(startIndex, endIndex) )
        coordsSub = self.coords[ beadIndices, :  ].copy() - c ## rest of the chain
        #axis = coordsSub[1]  # rotation axis
        axis = self.coords[ torsionQuad[2], :  ] - c
        axis /= np.linalg.norm( axis ) # unit vector
        
        rad = np.radians( theta )
        q0 = np.cos( rad/2 )
        q1 = np.sin( rad/2 ) * axis[0] 
        q2 = np.sin( rad/2 ) * axis[1]
        q3 = np.sin( rad/2 ) * axis[2]
        Q = [ [(q0**2 + q1**2 - q2**2 - q3**2),      2*(q1*q2 - q0*q3),          2*(q1*q3 + q0*q2)], 
             [ 2*(q2*q1 + q0*q3),     (q0**2 - q1**2 + q2**2 - q3**2),      2*(q2*q3 - q0*q1) ], 
             [ 2*(q3*q1 - q0*q2),          2*(q3*q2 + q0*q1),     (q0**2 - q1**2 - q2**2 + q3**2) ] ]
        Q = np.array( Q )
        newCoords = np.matmul( coordsSub, Q.T )  

        for i, bi in enumerate( beadIndices ):
            self.beads[bi].loc = newCoords[i] + c
        self.update1()
        
        return True
    
    def setAngle( self, angleTriple, theta=10, endIndex=0, beadIndices=[]  ):
        theta0 = getAngle( self.coords[angleTriple,:] ) 
        logging.debug( ['setangle', angleTriple, theta0] )
        s = self.rotateAngle(angleTriple, theta - theta0 )
        if DEBUG:
            theta0 = getAngle( self.coords[angleTriple,:] )
            logging.debug( ['setangle', angleTriple, theta0] )
        return s
 
    def setTorsionAngle( self, torsionQuad, theta=10, endIndex=0, beadIndices=[]  ):
        theta0 = getTorsionAngle( self.coords[torsionQuad,:] )
        logging.debug( ['settorsion', torsionQuad, theta0] )
        s = self.rotateTorsion(torsionQuad, theta - theta0 )
        if DEBUG:
            theta0 = getTorsionAngle( self.coords[torsionQuad,:] )
            logging.debug( ['settorsion', torsionQuad, theta0] )
        return s
    
def write_pdb(chain, filename, rows=[], mode='a'): #writes pdb based on atom records 
    #	try:
        with open(filename, mode) as f:
            if len(rows) == 0:
                rows = np.arange( chain.nbeads)
            if len( chain.residues ) == 0:
                residues = [ 'UNK' ] * chain.nbeads
            else:
                residues = chain.residues
            lineDict  = { }
            for i in rows:

                a = chain[i] #atom
                atomname = '{0:<4s}'.format( a.name )
                line = 'ATOM  {0:>5d} {1:<4s}{2:1s}{3:>3s} {4:1s}{5:>4d}{6:1s}   {7:>8.3f}{8:>8.3f}{9:>8.3f}{10:>6.2f}{11:>6.2f}{12:<10s}{13:>2s}{14:>2.1s}'
                line = line.format(a.no, atomname, '', residues[a.resno], a.chain, a.resno, '', a.loc[0], a.loc[1], a.loc[2], 1, 0, ' ', a.code[0], '')
                line += '\n'
                
                if a.resno in lineDict.keys():
                    lineDict[a.resno].append( line )
                else:
                    lineDict[a.resno] = [line]
                    
            resnos = list( set( list( lineDict.keys() ) ) )
            for i in resnos:
                for line in lineDict[i]:
                    f.write(line)
        return None
    
CT_N_Ca = { 'd': 1.47,  'bb' : True }
CT_Ca_C = { 'd': 1.53,  'bb' : True, 'ang' : [ 121 ] } # CA-C-N or CA-C-O
CT_C_N = { 'd': 1.32,  'bb' : True,  'ang' : [ 123 ] } # C-N-C
CT_C_O = {'d': 1.24,  'bb' : True, }

CT_C_Na = { 'd': 1.352,  'bb' : True, }
CT_C_Na_His = { 'd': 1.37,  'bb' : True, }
CT_C_Nd = { 'd': 1.32,  'bb' : True, }
CT_C_OH = {'d': 1.43,  'bb' : True, }
CT_C_Od = {'d': 1.23,  'bb' : True, }

CT_C_Na_Trp = { 'd': 1.384,  'bb' : True, }
CT_Car_Car_Trp = {'d':1.384, 'bb' : False, 'ang' : [121.3] }

CT_Car_Car = {'d':1.41, 'bb' : False, 'ang' : [120] }
CT_C_C = {'d':1.535, 'bb' : False, 'ang' : [111.17] }
CT_Cd_Cd = {'d':1.339, 'bb' : False, 'ang' : [121.3] }
CT_C_S = {'d':1.83, 'bb' : False, 'ang' : [111.17] }
#if initConfor == 'alpha':
    
        
def buildPeptide( resList, chname='A', loc=np.zeros( (ndim) ), initConfor0='alpha', Ncap='' ):
    ch1 = chain( [], [], chname, loc, initConfor0)
    ch1.residues = resList
    ch1.bbInd = []
    ch1.bbCode = []
    
    initConfor = 'linear' # sc
    initConfor1 = 'linear' # mc
    phi = -140; psi = 130
    
    phi = -57; psi=-47
    #phi = -80; psi=150
    ## add backbone first
    ## add first residue with N terminal atom/ Ncap
    resno = -1
    
    ## declare conType
    if Ncap == '':
        ch1.grow(bcode='NT', bname='NT', addtoBead=-1, resno=0, conType=CT_C_N )
        biN = 0 ## reference for backbone N atom
        #ch1.bbInd.extend( []  ) ## add to backbone except the last atom - which will be added at the end
    else:
        ch1.grow(bcode='NT', bname='NT', addtoBead=-1, resno=0, conType=CT_C_N )
        biN = 0 
        #ch1.bbInd.extend( []  )

    biCAo = '' ## old residue
    biCo = ''  ## old residue
    for resno, resName in enumerate( resList ):
        if resno > 0:
            biCAo = biCA ## old residue
            biCo = biC  ## old residue

        ## add mainchain
        if resno != 0: ## add spl Ncaps
            ch1.grow(bcode='N', bname='N', addtoBead=biC, resno=resno, conType=CT_N_Ca )
            biN = ch1[-1].no
        ch1.grow(bcode='CA', bname='CA', addtoBead=biN, resno=resno, conType=CT_N_Ca )
        biCA = biN + 1 ## index for CA atom too grow sidechain
        ### add sidechain - ACDEFGHIKLMNPQRSTVWY
        initConfor = 'linear'
        bi = biCA ## both same in this version
        if resName == 'GLY':
            None
        elif resName == 'ALA':
            ch1.grow(bcode='CB', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
        elif resName == 'VAL':
            ch1.grow(bcode='CB', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='CG', bname='CG1', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='CG', bname='CG2', addtoBead=bi+1, resno=resno, conType=CT_C_C )
        elif resName == 'LEU':
            ch1.grow(bcode='CB', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='CG', bname='CG', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='CD', bname='CD1', addtoBead=bi+2, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='CD', bname='CD2', addtoBead=bi+2, resno=resno, conType=CT_C_C )
        elif resName == 'ILE':
            ch1.grow(bcode='C', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CG', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CG', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CD', addtoBead=bi+2, resno=resno, conType=CT_C_C )

        elif resName == 'PRO':
            ch1.grow(bcode='C', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CG', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CD', addtoBead=bi+2, resno=resno, conType=CT_C_C )
            ch1.bondPair.append( [biCA-1, bi+3] ) ## add bond between N to Cdelta
        elif resName == 'CYS':
            ch1.grow(bcode='C', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='S', addtoBead=bi+1, resno=resno, conType=CT_C_C )
        elif resName == 'MET':
            ch1.grow(bcode='C', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CG', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='S', bname='S', addtoBead=bi+2, resno=resno, conType=CT_C_S )
            ch1.grow(bcode='C', bname='CE', addtoBead=bi+3, resno=resno, conType=CT_C_C )

        elif resName == 'THR':
            ch1.grow(bcode='CB', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='CG', bname='CG', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='OH', bname='OH', addtoBead=bi+1, resno=resno, conType=CT_C_OH )
 
        elif resName == 'SER':
            ch1.grow(bcode='C', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='OH', bname='OH', addtoBead=bi+1, resno=resno, conType=CT_C_OH )
        elif resName == 'ASN':
            ch1.grow(bcode='C', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='C', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='O', addtoBead=bi+2, resno=resno, conType=CT_C_Od )
            ch1.grow(bcode='C', bname='N', addtoBead=bi+2, resno=resno, conType=CT_C_Na )
        elif resName == 'GLN':
            ch1.grow(bcode='C', bname='CG', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='C', addtoBead=bi+2, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='O', addtoBead=bi+3, resno=resno, conType=CT_C_Od )
            ch1.grow(bcode='C', bname='N', addtoBead=bi+3, resno=resno, conType=CT_C_Na )
            ch1.grow(bcode='C', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
        elif resName == 'ASP':
            ch1.grow(bcode='C', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='C', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='O', bname='O', addtoBead=bi+2, resno=resno, conType=CT_C_O )
            ch1.grow(bcode='O', bname='O', addtoBead=bi+2, resno=resno, conType=CT_C_O )
        elif resName == 'GLU':
            ch1.grow(bcode='C', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CG', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='C', addtoBead=bi+2, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='O', bname='O', addtoBead=bi+3, resno=resno, conType=CT_C_O )
            ch1.grow(bcode='O', bname='O', addtoBead=bi+3, resno=resno, conType=CT_C_O )
        elif resName == 'HIS':
            ch1.grow(bcode='C', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CG', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CD2', addtoBead=bi+2, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='NE2', addtoBead=bi+3, resno=resno, conType=CT_C_Na_His )
            ch1.grow(bcode='C', bname='CE1', addtoBead=bi+4, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='ND1', addtoBead=bi+5, resno=resno, conType=CT_C_Nd )
            ch1.bondPair.append( [bi+2, bi+6] ) ## add bond between NH to Cgamma
        elif resName == 'LYS':
            ch1.grow(bcode='C', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CG', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CD', addtoBead=bi+2, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CE', addtoBead=bi+3, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='N', bname='N', addtoBead=bi+4, resno=resno, conType=CT_C_Na )
        elif resName == 'ARG':
            ch1.grow(bcode='C', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CG', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CD', addtoBead=bi+2, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='NE', addtoBead=bi+3, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CZ', addtoBead=bi+4, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='N', bname='NH', addtoBead=bi+5, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='N', bname='NH', addtoBead=bi+5, resno=resno, conType=CT_C_C )
            
        elif resName == 'TRP':
            ch1.grow(bcode='C', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CG', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CD', addtoBead=bi+2, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CE', addtoBead=bi+3, resno=resno, conType=CT_Car_Car_Trp )
            ch1.grow(bcode='N', bname='NE', addtoBead=bi+4, resno=resno, conType=CT_C_Na_Trp )
            ch1.grow(bcode='C', bname='CD', addtoBead=bi+5, resno=resno, conType=CT_C_C )
            ch1.bondPair.append( [bi+1, bi+6] )
            
            ch1.grow(bcode='C', bname='CZ', addtoBead=bi+4, resno=resno, conType=CT_Car_Car_Trp )
            ch1.grow(bcode='C', bname='CH', addtoBead=bi+7, resno=resno, conType=CT_Car_Car_Trp )
            ch1.grow(bcode='C', bname='CZ', addtoBead=bi+8, resno=resno, conType=CT_Car_Car_Trp )
            ch1.grow(bcode='C', bname='CE', addtoBead=bi+9, resno=resno, conType=CT_Car_Car_Trp)
            ch1.bondPair.append( [bi+6, bi+10] ) ## complete the ring
        elif resName == 'TYR':
            ch1.grow(bcode='C', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CG', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CD', addtoBead=bi+2, resno=resno, conType=CT_Car_Car )
            ch1.grow(bcode='C', bname='CE', addtoBead=bi+3, resno=resno, conType=CT_Car_Car )
            ch1.grow(bcode='C', bname='CZ', addtoBead=bi+4, resno=resno, conType=CT_Car_Car )
            ch1.grow(bcode='C', bname='CE', addtoBead=bi+5, resno=resno, conType=CT_Car_Car )
            ch1.grow(bcode='C', bname='CD', addtoBead=bi+6, resno=resno, conType=CT_Car_Car )
            ch1.grow(bcode='C', bname='OH', addtoBead=bi+5, resno=resno, conType=CT_C_OH )
            ch1.bondPair.append( [bi+2, bi+7] ) ## complete the ring
        elif resName == 'PHE':
            ch1.grow(bcode='C', bname='CB', addtoBead=biCA, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CG', addtoBead=bi+1, resno=resno, conType=CT_C_C )
            ch1.grow(bcode='C', bname='CD', addtoBead=bi+2, resno=resno, conType=CT_Car_Car )
            ch1.grow(bcode='C', bname='CE', addtoBead=bi+3, resno=resno, conType=CT_Car_Car )
            ch1.grow(bcode='C', bname='CZ', addtoBead=bi+4, resno=resno, conType=CT_Car_Car )
            ch1.grow(bcode='C', bname='CE', addtoBead=bi+5, resno=resno, conType=CT_Car_Car )
            ch1.grow(bcode='C', bname='CD', addtoBead=bi+6, resno=resno, conType=CT_Car_Car )
            ch1.bondPair.append( [bi+2, bi+7] ) ## complete the ring
            print( bi )
            tsets = [ [2, 3, 4], [3, 4, 5], [4, 5, 6], [5, 6, 7], [6, 7, 2], [7, 2, 3] ]
            tall = [ bi+x for x in range(2,7) ]
            for t in tsets:
                a = [ bi + x for x in t ]
                ch1.setAngle( a, 120, beadIndices= a )
            
        ### added side chain - add rest of backbone   
        ch1.grow(bcode='C', bname='C', addtoBead=biCA, resno=resno, conType=CT_Ca_C )
        biC = ch1[-1].no
        if resno < len(resList)-1:
            ch1.grow(bcode='O', bname='O', addtoBead=biC, resno=resno, conType=CT_C_O )
        else:
         ## add C-terminal oxygen  
            ch1.grow(bcode='OT', bname='OT', addtoBead=biC, resno=resno, conType=CT_C_O )
        ch1.bbInd.extend( [ biN, biCA, biC, biC+1 ] )
         # [ 114, 121 ] } # CA-C-N or CA-C-O
         
        
        ### adjust bond angle and torsion angles
        ## adjust bond angle
        sign = 1 #  (-1 if resno%2 != 0 else 1) * -1
        
        if False:
            if resno > 0:
                ch1.setAngle( [biCo, biCA-1, biCA ], -123 * sign ) # C-N-Ca
                ch1.setAngle( [ biCAo, biCo, biCA-1 ], 114 * sign )  #Ca-C-N
                ch1.setAngle( [ biCo+1, biCo, biCA-1 ], -125 * sign ) # O-C-N
            ch1.setAngle( [ biCA-1, biCA, biC ], 110 * sign )  # N-Ca-C
            ch1.setAngle( [ biCA, biC, biC+1 ], 121 * sign )  # Ca-C-O
            
            if resName != 'GLY':
                ch1.setAngle( [biCA-1, biCA, biCA+1 ], -111.17 * sign, beadIndices= list(range( biCA, biC )) ) 
                    
            ## check
            if False and resno > 0:
                for a in [[biCo, biCA-1, biCA ], [ biCo+1, biCo, biCA-1 ], [ biCAo, biCo, biCA-1 ], 
                          [ biCA-1, biCA, biC ], [ biCA, biC, biC+1 ] , [biCA-1, biCA, biCA+1 ] ]:
                    print( [resno, [-123, -125, -114, -119, -125, 120], getAngle( ch1.coords[ a ] ) ] )

        ## adust torsion angles -- bb
        if True and resno > 0:
            if resName != 'PRO':
                ch1.setTorsionAngle( [biCAo, biCo, biCA-1, biCA ], 180 ) # Ca-C-N-Ca
            else:
                ch1.setTorsionAngle( [biCAo, biCo, biCA-1, biCA ], 0 ) 

            ch1.setTorsionAngle( [biCo, biCA-1, biCA, biC ], phi * sign ) # C-N-Ca-C
            ch1.setTorsionAngle( [biCAo-1, biCAo, biCo, biCA-1 ], psi * sign) # N-Ca-C-N

        ## check
        if False and resno > 0:
            for a in [[biCo, biCA-1, biCA ], [ biCo+1, biCo, biCA-1 ], [ biCAo, biCo, biCA-1 ], 
                      [ biCA-1, biCA, biC ], [ biCA, biC, biC+1 ] , [biCA-1, biCA, biCA+1 ] ]:
                print( [resno, [-123, -125, -114, -119, -125, 120], getAngle( ch1.coords[ a ] ) ] )
                
        if True and resno > 0:
            for a in [ [biCAo, biCo, biCA-1, biCA ], [biCo, biCA-1, biCA, biC ], [biCAo-1, biCAo, biCo, biCA-1 ] ]:
                print( [resno, ['omega, phi, psi'], getTorsionAngle( ch1.coords[ a ] ) ] )

    for resno, res in enumerate( resList):           
        if True and resno > 0:
            for a in [ [biCAo, biCo, biCA-1, biCA ], [biCo, biCA-1, biCA, biC ], [biCAo-1, biCAo, biCo, biCA-1 ] ]:
                print( [resno, ['omega, phi, psi'], getTorsionAngle( ch1.coords[ a ] ) ] )

    ####
    ch1.bbCode = [ ch1[b].code for b in ch1.bbInd ]
    return ch1