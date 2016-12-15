from component_contribution.compound import Compound
from component_contribution.thermodynamic_constants import default_T

def formation_energy(listKeggMolecules, cc, pH=7, I=0.1, T=298.15):

    for cid in listKeggMolecules:
        comp = Compound.from_kegg(cid)
        print comp


        major_ms_dG0_f = cc.get_major_ms_dG0_f(cid)
    #    for d in comp.get_species(major_ms_dG0_f, default_T):
    #        print d['dG0_f']
        a = comp.get_species(major_ms_dG0_f, default_T)
        dG0_f = a.next()['dG0_f']

        #Compound.transform_pH7(comp, 7, .1, 298.15) # it gives the difference between the transformed molecule and the MS species at pH7
        #Compound.transform_neutral(comp, 7, .1, 298.15)  # ???
        dG0f_dG0tf = Compound.get_transform(comp, pH, I, T)

        return dG0_f + dG0f_dG0tf

    #    for i in np.arange(len(comp.nHs)):
    #        # DO THIS FOR EVERY SPECIES!
    #        print  comp.transform(i, 7.0, 0.1, 298)

            # print "%s: ddG0_prime = %.2f" % (cid, comp.transform(i,7.0, 0.5, 300))
            # print "%s: ddG0_prime(z=0) = %.2f" % (cid, comp.transform_neutral(7.0, 0.2, 300))