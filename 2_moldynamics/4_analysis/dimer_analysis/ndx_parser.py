import MDAnalysis as mda
import numpy as np

class GmxNDX:
    
    def __init__(self, filename: str):
        self._filename = filename
        self._groupnames = []
        self._groupindices = []
        self._parse(filename)

    def __getitem__(self, item):
        if isinstance(item, int):
            grpname = self._groupnames[item]
            grpidcs = self._groupindices[item]
        elif isinstance(item, str):
            item = self._groupnames.index(item)
            grpname = self._groupnames[item]
            grpidcs = self._groupindices[item]
        else:
            raise ValueError(f"Cannot handle input of type: {type(item)}")
        return grpname, grpidcs
    
    def display_groups(self):
        for i, (grpname, grpidcs) in enumerate(zip(self._groupnames, self._groupindices)):
            print(f"Group{i:6d} ({grpname:>15}) has {grpidcs.size:5d} elements")
        
    def select_atomgroup(self, universe: mda.Universe, prompt="Select a group:") -> mda.AtomGroup:
        self.display_groups()
        sele = input(prompt)
        if sele.isdigit():
            sele = int(sele)
        grpname, grpidcs = self[sele]
        print(f"Selected '{grpname}'")
        return universe.atoms[grpidcs - 1]
        
    def _parse(self, filename):
        grpname = None
        with open(filename, "r") as fhandler:
            for ln in fhandler:
                if ln.startswith("["):
                    if grpname is not None:
                        self._groupnames.append(grpname)
                        self._groupindices.append(np.array(grpidcs, dtype=int))
                    grpname = ln.split()[1]
                    grpidcs = []
                else:
                    grpidcs += ln.split()
            self._groupnames.append(grpname)
            self._groupindices.append(np.array(grpidcs, dtype=int))