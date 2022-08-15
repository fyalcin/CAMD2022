from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab, SlabGenerator
import numpy as np


# struct = Structure.from_file("POSCAR")
# asf = AdsorbateSiteFinder(struct)

def slab_from_structure(miller, structure):
    return Slab(lattice=structure.lattice,
                species=structure.species_and_occu,
                coords=structure.frac_coords,
                miller_index=miller,
                oriented_unit_cell=structure,
                shift=0,
                scale_factor=np.eye(3, dtype=np.int),
                site_properties=structure.site_properties)


from pymatgen.ext.matproj import MPRester

with MPRester() as m:
    struct = m.get_structure_by_material_id('mp-2', conventional_unit_cell=True)

miller_index = (1, 1, 1)
sg = SlabGenerator(initial_structure=struct,
                   miller_index=miller_index,
                   min_slab_size=10,
                   min_vacuum_size=10,
                   max_normal_search=max([abs(i) for i in miller_index]),
                   primitive=True,
                   lll_reduce=True,
                   in_unit_planes=True,
                   center_slab=True)

slab = sg.get_slabs()[0]
slab.make_supercell([2, 2, 1])


# struct = Structure.from_file('POSCAR_Pd_79.vasp')
# slab = slab_from_structure((1, 1, 1), struct)
asf = AdsorbateSiteFinder(slab)
sites = asf.find_adsorption_sites(symm_reduce=False)
from pymatgen.core.structure import Molecule
site_coords = list(sites.values())
site_coords = [a for b in site_coords for a in b]
mol = Molecule(['H'], [[0, 0, 0]])
# for site in site_coords:
#     asf = AdsorbateSiteFinder(slab)
#     slab = asf.add_adsorbate(molecule=mol, ads_coord=site)
# slab.to('poscar', 'pd_slab.vasp')
