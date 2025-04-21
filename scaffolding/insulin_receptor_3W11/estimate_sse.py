'''
Secondary Structure Elements (SSE) of a peptide chain
SSE refers to its local spatial arrangements, 
primarily α-helices and β-sheets, along with turns and coils.
'''
import numpy as np

_radians_to_angle = 2 * np.pi / 360

_r_helix = [i * _radians_to_angle for i in (89-12, 89+12)]
_a_helix = [i * _radians_to_angle for i in (50-20, 50+20)]
_d2_helix = (5.5-0.5, 5.5+0.5)
_d3_helix = (5.3-0.5, 5.3+0.5)
_d4_helix = (6.4-0.6, 6.4+0.6)

_r_strand = [i * _radians_to_angle for i in (124-14, 124+14)]
_a_strand = [i * _radians_to_angle for i in (-180, -125, 145, 180)]
_d2_strand = (6.7-0.6, 6.7+0.6)
_d3_strand = (9.9-0.9, 9.9+0.9)
_d4_strand = (12.4-1.1, 12.4+1.1)

class EstimateSSE:

    def __init__(self, xyz):
        # coordinates of atoms of a chain
        self.ca_coord = xyz
        self.ca_num = len(xyz)
        self.sse = ["L"] * self.ca_num

    def get_sse(self):
        '''
        calculates the SSE of a peptide chain based on the P-SEA algorithm (Labesse 1997)
        '''
        # print('###ca_coord:', ca_coord.shape, ca_coord)

        # Filter all CA atoms in the relevant chain.
        d2i = self.cal_d2i()
        d3i = self.cal_d3i()
        d4i = self.cal_d4i()
        ri = self.cal_ri()
        ai = self.cal_ai()
           
        # Annotate helices
        # Find CA that meet criteria for potential helices
        is_pot_helix = self.cal_is_pot_helix(d3i, d4i, ri, ai)
        # Real helices are 5 consecutive helix elements
        is_helix = self.cal_is_helix(is_pot_helix)
        # update self.sse
        self.sse_helix(d3i, ri, is_helix)
        
        # Annotate sheets
        # Find CA that meet criteria for potential strands
        is_pot_strand = self.cal_is_pot_strand(d2i,d3i, d4i, ri, ai)
        # Real strands are 5 consecutive strand elements,
        # or shorter fragments of at least 3 consecutive strand residues,
        # if they are in hydrogen bond proximity to 5 other residues
        is_strand = self.cal_is_strand(is_pot_strand)
        # update self.sse
        self.sse_strand(d3i, is_strand)
        return self.sse

###################
    def cal_is_pot_helix(self, d3i, d4i, ri, ai):
        '''
        Annotate helices: 
        Find CA that meet criteria for potential helices
        '''
        is_pot_helix = np.zeros(self.ca_num, dtype=bool)
        for i in range(self.ca_num):
            cond1 = _d3_helix[0] <= d3i[i] <= _d3_helix[1] \
                and _d4_helix[0] <= d4i[i] <= _d4_helix[1]
            cond2 = _r_helix[0] <= ri[i] <= _r_helix[1] \
                and _a_helix[0] <= ai[i] <= _a_helix[1]
            if cond1 or cond2:
                is_pot_helix[i] = True
        return is_pot_helix

    def cal_is_helix(self, is_pot_helix):
        '''
        Real helices are 5 consecutive helix elements
        '''
        is_helix = np.zeros(self.ca_num, dtype=bool)
        counter = 0
        for i in range(self.ca_num):
            if is_pot_helix[i]:
                counter += 1
            else:
                if counter >= 5:
                    is_helix[i-counter:i] = True
                counter = 0
        return is_helix

    def sse_helix(self, d3i, ri, is_helix):
        '''
        'H' represent helix
        Extend the helices by one at each end if CA meets extension criteria
        '''
        for i in np.where(is_helix == True)[0]:
            # current: i
            self.sse[i] = "H"
            for j in (i-1, i+1):
                if _d3_helix[0] <= d3i[j] <= _d3_helix[1] or \
                    _r_helix[0] <= ri[j] <= _r_helix[1]:
                    self.sse[j] = "H"
        return self.sse

    def cal_is_pot_strand(self, d2i, d3i, d4i, ri, ai):
        '''
        Annotate sheets
        Find CA that meet criteria for potential strands
        '''
        is_pot_strand = np.zeros(self.ca_num, dtype=bool)
        for i in range(self.ca_num):
            cond1 = _d2_strand[0] <= d2i[i] <= _d2_strand[1] \
                and _d3_strand[0] <= d3i[i] <= _d3_strand[1] \
                and _d4_strand[0] <= d4i[i] <= _d4_strand[1]
            cond2 = ri[i] >= _r_strand[0] and ri[i] <= _r_strand[1] \
                and (_a_strand[0] <= ai[i] <= _a_strand[1]) \
                    or (_a_strand[2] <= ai[i] <= _a_strand[3])
            if cond1 or cond2:
                is_pot_strand[i] = True
        return is_pot_strand



    # def get_pot_strand_coord(self, is_pot_strand):
    #     pot_strand_coord = self.ca_coord[is_pot_strand]
    #     return pot_strand_coord

    def cal_is_strand(self, is_pot_strand):
        '''
        Real strands are 5 consecutive strand elements,
        or shorter fragments of at least 3 consecutive strand residues,
        if they are in hydrogen bond proximity to 5 other residues
        '''
        is_strand = np.zeros(self.ca_num, dtype=bool)
        counter = 0
        contacts = 0
        for i in range(self.ca_num):
            if is_pot_strand[i]:
                counter += 1
                coord = self.ca_coord[i]
                for strand_coord in self.ca_coord:
                    dist = self.distance(coord, strand_coord)
                    if 4.2 <= dist <= 5.2:
                        contacts += 1
            else:
                if counter >= 4 or (counter == 3 and contacts >= 5):
                    is_strand[i-counter : i] = True
                counter = 0
                contacts = 0
        return is_strand

    def sse_strand(self, d3i, is_strand):
        '''
        'E' represent strand
        '''
        for i in np.where(is_strand==True)[0]:
            self.sse[i] = "E"
            for j in (i-1, i+1):
                if _d3_strand[0] <= d3i[j] <= _d3_strand[1]:
                    self.sse[j] = "E"
        return self.sse

    def cal_d2i(self):
        '''
        2D distance to index
        '''
        d2i_coord = np.full((self.ca_num, 2, 3 ), np.nan)
        for i in range(1, self.ca_num-1): 
            d2i_coord[i] = (self.ca_coord[i-1], self.ca_coord[i+1])
        d2i = self.distance(d2i_coord[:,0], d2i_coord[:,1])
        return d2i

    def cal_d3i(self):
        '''
        3D distance to index
        '''
        d3i_coord = np.full((self.ca_num, 2, 3 ), np.nan)
        for i in range(1, self.ca_num-2): 
            d3i_coord[i] = (self.ca_coord[i-1], self.ca_coord[i+2])
        d3i = self.distance(d3i_coord[:,0], d3i_coord[:,1])
        return d3i

    def cal_d4i(self):
        '''
        Distance-Based 4-Dimensional Indexing
        '''
        d4i_coord = np.full((self.ca_num, 2, 3 ), np.nan)
        for i in range(1, self.ca_num-3): 
            d4i_coord[i] = (self.ca_coord[i-1], self.ca_coord[i+3])
        d4i = self.distance(d4i_coord[:,0], d4i_coord[:,1])
        return d4i

    def cal_ri(self):
        '''
        radius to index
        '''
        ri_coord = np.full((self.ca_num, 3, 3 ), np.nan)
        for i in range(1, self.ca_num-1): 
            ri_coord[i]  = (self.ca_coord[i-1], self.ca_coord[i], self.ca_coord[i+1])
        ri = self.angle(ri_coord[:,0], ri_coord[:,1], ri_coord[:,2])
        return ri
    
    def cal_ai(self):
        '''
        '''
        ai_coord = np.full((self.ca_num, 4, 3 ), np.nan)
        for i in range(1, self.ca_num-2): 
            ai_coord[i]  = (
                self.ca_coord[i-1], self.ca_coord[i],
                self.ca_coord[i+1], self.ca_coord[i+2]
            )
        ai = self.dihedral(
            ai_coord[:,0], ai_coord[:,1],
            ai_coord[:,2], ai_coord[:,3]
        )
        return ai

###################
    def vector_dot(self, v1,v2):
        return (v1*v2).sum(-1)

    def norm_vector(self, v):
        scaler = np.linalg.norm(v, axis=-1, keepdims=True)
        return v / scaler

    def displacement(self, atoms1, atoms2):
        v1 = np.asarray(atoms1)
        v2 = np.asarray(atoms2)
        if len(v1.shape) <= len(v2.shape):
            diff = v2 - v1
        else:
            diff = -(v1 - v2)
        return diff

    def distance(self, atoms1, atoms2):
        diff = self.displacement(atoms1, atoms2)
        dist = self.vector_dot(diff, diff)
        return np.sqrt(dist)

    def angle(self, atoms1, atoms2, atoms3):
        v1 = self.norm_vector(self.displacement(atoms1, atoms2))
        v2 = self.norm_vector(self.displacement(atoms3, atoms2))
        v = self.vector_dot(v1,v2)
        return np.arccos(v)

    def dihedral(self, atoms1, atoms2, atoms3, atoms4):
        v1 = self.norm_vector(self.displacement(atoms1, atoms2))
        v2 = self.norm_vector(self.displacement(atoms2, atoms3))
        v3 = self.norm_vector(self.displacement(atoms3, atoms4))
        
        n1 = np.cross(v1, v2)
        n2 = np.cross(v2, v3)
        
        # Calculation using atan2, to ensure the correct sign of the angle 
        x = self.vector_dot(n1,n2)
        y = self.vector_dot(np.cross(n1,n2), v2)
        return np.arctan2(y, x)