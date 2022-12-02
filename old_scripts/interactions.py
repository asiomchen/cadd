import prolif as plf
import MDAnalysis as mda

#load protein
prot = mda.Universe("/home/anton/in_dev/Docking_tools/master/data/5kir_cleaned.pdb",guess_bonds=True, vdwradii={'CO': 1.5})
prot = plf.Molecule.from_mda(prot)
print(prot.n_residues)

# load ligands
path = str('/home/anton/in_dev/Docking_tools/master/smina_16-results/3-CHEMBL6_smina_docked.sdf')
lig_suppl = list(plf.sdf_supplier(path))
# generate fingerprint
fp = plf.Fingerprint()
fp.run_from_iterable(lig_suppl, prot)
results_df = fp.to_dataframe(return_atoms=True)
print(results_df)