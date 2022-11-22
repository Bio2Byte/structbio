"""
DISCLAIMER:
    IF YOU ARE LOOKING AT THIS CODE, AND YOU ARE NOT ME, 
    YOU PROBABLY SHOULD NOT
    
    THIS CODE WAS WRITTEN IN A QUICK AND DIRTY FASHION TO RUN ONCE,
    AND BE REFACTORED BEFORE EVER BEING RUN AGAIN.

Author: David Bickel
Date:   Nov 21, 2022
"""
import argparse
import os, sys

from Bio import SeqUtils, pairwise2
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.base import AnalysisBase
import numpy as np
from scipy.spatial.transform import Rotation

from .blosum62 import BLOSUM62
from .ndx_parser import GmxNDX


def shared_atoms(atomgroup_A: mda.AtomGroup, atomgroup_B: mda.AtomGroup) -> tuple[mda.AtomGroup]:
    """ Compare two atom groups and return corresponding CA, based on a pairwise sequence alignment """
    # Only use Calphas
    calphas_A, calphas_B = atomgroup_A.select_atoms("name CA"),  atomgroup_B.select_atoms("name CA")
    aln = pairwise2.align.globaldd(get_sequence(calphas_A), get_sequence(calphas_B), BLOSUM62, -10, -.5, -10, -.5).pop()
    for aaA, caA, aaB, caB in zip(aln.seqA, calphas_A, aln.seqB, calphas_B):
        if aaA == '-':
            calphas_B -= caB
        elif aaB == '-':
            calphas_A -= caA
    return calphas_A, calphas_B    
    
def get_sequence(ag: mda.AtomGroup) -> str:
    """ Returns a sequence for an AtomGroup """
    aa = [SeqUtils.IUPACData.protein_letters_3to1[r.title()] for r in ag.resnames]
    return "".join(aa)

def unify_orientation(reference_coords, mobile_coords):
    """ Unifies the orientation of a structure such, that x, y, z == PC1, PC2, PC3 """
    # Translation
    tram = np.mean(reference_coords, axis=0)   # translation matrix
    rcoord = reference_coords - tram
    mcoord = mobile_coords - tram

    # Calculate rotation matrix so that x == PC1, y == PC2, z == PC3
    covariance_matrix = np.cov(reference_coords.T)
    eigen_values, eigen_vectors = np.linalg.eig(covariance_matrix)
    order = np.argsort(eigen_values)[::-1]
    eigen_vectors, eigen_values = eigen_vectors[:, order], eigen_values[order]
    eigen_vectors[:,2] *= np.linalg.det(eigen_vectors)
    rotm = np.linalg.inv(eigen_vectors)
    
    rcoord, mcoord = np.matmul(rotm, rcoord.T).T, np.matmul(rotm, mcoord.T).T
    return rcoord, mcoord

class DimerAnalysis(AnalysisBase):
    
    REFERENCE_FILE = f"{os.path.dirname(__file__)}/1x75_superimposed.pdb"
    
    def __init__(self, static_group, mobile_group,  **kwargs):
        super(self.__class__, self).__init__(static_group.universe.trajectory, **kwargs)
        self._staticgroup = static_group
        self._mobilegroup = mobile_group
        self._reference = mda.Universe(self.__class__.REFERENCE_FILE)
        self._reference_group = self._reference.select_atoms("chainID A")
        
        chain_A, chain_B = shared_atoms(self._reference_group, self._reference.select_atoms("chainID B"))
        self._refR, _rms = align.rotation_matrix(chain_B.positions - chain_B.atoms.center_of_mass(), chain_A.positions - chain_A.atoms.center_of_mass())

    def _prepare(self):
        # OPTIONAL
        # Called before iteration on the trajectory has begun.
        # Data structures can be set up at this time
        # self.results.example_result = []
        self._staticgroup, self._mobilegroup = shared_atoms(self._staticgroup, self._mobilegroup)
        self._reference_group, self._staticgroup = shared_atoms(self._reference_group, self._staticgroup)
        self._reference_group, self._mobilegroup = shared_atoms(self._reference_group, self._mobilegroup)        
        self._reference_group.atoms.translate(-self._reference_group.atoms.center_of_mass())
        self.results.time = []
        self.results.rotation_angle = []

    def _single_frame(self):
        # REQUIRED
        # Called after the trajectory is moved onto each new frame.
        # store an example_result of `some_function` for a single frame
        stccrd, mobcrd = self._staticgroup.positions, self._mobilegroup.positions
        
        # Translate to center staticcrd on origin
        T = -self._staticgroup.center_of_mass()
        self._staticgroup.translate(T)
        self._mobilegroup.translate(T)
        
        # Rotate coordinates to align with reference structure
        R, _rms = align.rotation_matrix(self._staticgroup.positions, self._reference_group.positions)
        self._staticgroup.atoms.rotate(R)
        self._mobilegroup.atoms.rotate(R)
        
        # Actual calculation
        self._mobilegroup.atoms.translate(-self._mobilegroup.atoms.center_of_mass())
        R, _rms = align.rotation_matrix(self._mobilegroup.positions, self._staticgroup.positions)
        
        # Write output
        self.results.time.append(self._staticgroup.ts.time)
        self.results.rotation_angle.append(self.theta_rotation(R, self._refR))   
        
    def _conclude(self):
        # OPTIONAL
        # Called once iteration on the trajectory is finished.
        # Apply normalisation and averaging to results here.
        self.results.time = np.array(self.results.time)
        self.results.rotation_angle = np.array(self.results.rotation_angle)

    @staticmethod
    def theta_rotation(ra, rb=None):
        if rb is None: rb = np.identity(ra.shape[0])
        rab = np.matmul(ra, rb.T)
        theta = np.arccos((np.trace(rab) - 1) / 2.)
        return np.rad2deg(theta)
    
def parse_arguments(arguments):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--structure", help="A structure file [*.gro, *tpr, *pdb]")
    parser.add_argument("-f", "--trajectory", help="A trajectory file [*.xtc]")
    parser.add_argument("-n", "--index", help="An index file [*.ndx]")
    parser.add_argument("-o", "--outfile", help="The output file to which the results are written.")
    args = parser.parse_args(arguments)
    return args


def main(arguments):
    args = parse_arguments(arguments)
    
    # Load universe
    u = mda.Universe(args.structure, args.trajectory)
    ndx = GmxNDX(args.index)
    
    # Select dimers
    ag0 = ndx.select_atomgroup(u, "Select static group:")
    ag1 = ndx.select_atomgroup(u, "Select mobile group:")
    
    # Run analysis
    analy = DimerAnalysis(ag0, ag1)
    analy.run()
    
    with open(args.outfile, "w") as fhandle:
        fhandle.write("#Time     DeltaRotation\n")
        for time, angle in zip(analy.results.time, analy.results.rotation_angle):
            fhandle.write(f"{time:8.3f}  {angle:8.3f}\n")
            
    return 0
    
if __name__ == "__main__":
    print(sys.argv)
    sys.exit(main(sys.argv[1:]))
    