from mpinterfaces.nanoparticle import Nanoparticle
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from mpinterfaces import get_struct_from_mp
from pymatgen.ext.matproj import MPRester
from pymatgen.io.ase import AseAtomsAdaptor

aaa = AseAtomsAdaptor()

mpid = 'mp-2'

with MPRester() as m:
    surf_data = m.get_surface_data(mpid)['surfaces']
    struct = m.get_structure_by_material_id(mpid, conventional_unit_cell=True)

surfen_data = {tuple(k['miller_index']): k['surface_energy'] for k in surf_data}

from ase.cluster import wulff_construction

surfaces = list(surfen_data.keys())
esurf = list(surfen_data.values())  # Surface energies.
lc = struct.lattice.a
size = 200  # Number of atoms
atoms = wulff_construction('Pd', surfaces, esurf,
                           size, 'fcc',
                           rounding='above')

nanoparticle = aaa.get_structure(atoms)
