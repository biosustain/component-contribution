from component_contribution.compound import Compound
from component_contribution.thermodynamic_constants import default_T
from numpy import arange
import openbabel

def check_if_already_exists_inchi(inchi):

    with open('../data/compounds.tsv', 'r') as tsv:
        AoA = [line.strip().split('\t') for line in tsv]

    for i in arange(1,len(AoA)):
        if AoA[i].__len__()==3:
            found = AoA[i][2].find(inchi)
            if found != -1:
                #matched_ids[matched_ids.__len__()]=AoA[i][0]
                cid = AoA[i][0]
                name = AoA[i][1]
                return cid, name
    return None, None

def smile2standard_inchi(s):
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats('smiles', 'inchi')
    #conv.AddOption("F", conv.OUTOPTIONS)
    conv.AddOption("T", conv.OUTOPTIONS)
    conv.AddOption("x", conv.OUTOPTIONS, "noiso")
    conv.AddOption("w", conv.OUTOPTIONS)
    obmol = openbabel.OBMol()
    if not conv.ReadString(obmol, str(s)):
        return None
    inchi = conv.WriteString(obmol, True) # second argument is trimWhitespace
    if inchi == '':
        return None
    else:
        return inchi

def smile2standard_inchi2(s):
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats('smiles', 'inchi')
    #conv.AddOption("F", conv.OUTOPTIONS)
    #conv.AddOption("T", conv.OUTOPTIONS)
    conv.AddOption("x", conv.OUTOPTIONS, "noiso")
    conv.AddOption("w", conv.OUTOPTIONS)
    obmol = openbabel.OBMol()
    if not conv.ReadString(obmol, str(s)):
        return None
    inchi = conv.WriteString(obmol, True) # second argument is trimWhitespace
    if inchi == '':
        return None
    else:
        return inchi

def smile2standard_inchi3(s):
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats('smiles', 'inchi')
    #conv.AddOption("F", conv.OUTOPTIONS)
    conv.AddOption("T", conv.OUTOPTIONS)
    #conv.AddOption("x", conv.OUTOPTIONS, "noiso")
    conv.AddOption("w", conv.OUTOPTIONS)
    obmol = openbabel.OBMol()
    if not conv.ReadString(obmol, str(s)):
        return None
    inchi = conv.WriteString(obmol, True) # second argument is trimWhitespace
    if inchi == '':
        return None
    else:
        return inchi

def smile2standard_inchi4(s):
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats('smiles', 'inchi')
    #conv.AddOption("F", conv.OUTOPTIONS)
    conv.AddOption("T", conv.OUTOPTIONS)
    conv.AddOption("x", conv.OUTOPTIONS, "noiso")
    #conv.AddOption("w", conv.OUTOPTIONS)
    obmol = openbabel.OBMol()
    if not conv.ReadString(obmol, str(s)):
        return None
    inchi = conv.WriteString(obmol, True) # second argument is trimWhitespace
    if inchi == '':
        return None
    else:
        return inchi

def smile2standard_inchi5(s):
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats('smiles', 'inchi')
    #conv.AddOption("F", conv.OUTOPTIONS)
    #conv.AddOption("T", conv.OUTOPTIONS)
    #conv.AddOption("x", conv.OUTOPTIONS, "noiso")
    conv.AddOption("w", conv.OUTOPTIONS)
    obmol = openbabel.OBMol()
    if not conv.ReadString(obmol, str(s)):
        return None
    inchi = conv.WriteString(obmol, True) # second argument is trimWhitespace
    if inchi == '':
        return None
    else:
        return inchi

def smile2standard_inchi6(s):
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats('smiles', 'inchi')
    #conv.AddOption("F", conv.OUTOPTIONS)
    #conv.AddOption("T", conv.OUTOPTIONS)
    conv.AddOption("x", conv.OUTOPTIONS, "noiso")
    #conv.AddOption("w", conv.OUTOPTIONS)
    obmol = openbabel.OBMol()
    if not conv.ReadString(obmol, str(s)):
        return None
    inchi = conv.WriteString(obmol, True) # second argument is trimWhitespace
    if inchi == '':
        return None
    else:
        return inchi

def smile2standard_inchi7(s):
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats('smiles', 'inchi')
    #conv.AddOption("F", conv.OUTOPTIONS)
    conv.AddOption("T", conv.OUTOPTIONS)
    #conv.AddOption("x", conv.OUTOPTIONS, "noiso")
    #conv.AddOption("w", conv.OUTOPTIONS)
    obmol = openbabel.OBMol()
    if not conv.ReadString(obmol, str(s)):
        return None
    inchi = conv.WriteString(obmol, True) # second argument is trimWhitespace
    if inchi == '':
        return None
    else:
        return inchi

def mol2standard_inchi(s):
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats('mol', 'inchi')
    #conv.AddOption("F", conv.OUTOPTIONS)
    conv.AddOption("T", conv.OUTOPTIONS)
    conv.AddOption("x", conv.OUTOPTIONS, "noiso")
    conv.AddOption("w", conv.OUTOPTIONS)
    obmol = openbabel.OBMol()
    if not conv.ReadString(obmol, str(s)):
        return None
    inchi = conv.WriteString(obmol, True) # second argument is trimWhitespace
    if inchi == '':
        return None
    else:
        return inchi

def smile2mol(s):
    #openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats('smiles','mol')
    # conv.AddOption("F", conv.OUTOPTIONS)
    # conv.AddOption("T", conv.OUTOPTIONS)
    # conv.AddOption("x", conv.OUTOPTIONS, "noiso")
    # conv.AddOption("w", conv.OUTOPTIONS)
    obmol = openbabel.OBMol()
    if not conv.ReadString(obmol, str(s)):
        return None
    mol = conv.WriteString(obmol, True)  # second argument is trimWhitespace
    if mol == '':
        return None
    else:
        return mol



    #
    #
    #
    # import openbabel, pybel
    # conv = openbabel.OBConversion()
    # conv.ReadString()
    # a = pybel.readstring('smi', s)
    #
    # a.write

def mol_to_inchi_to_cid(mol):
    inchi = Compound.mol2inchi(mol)
    inchiS = mol2standard_inchi(mol)
    cid, Name = check_if_already_exists_inchi(inchiS)
    if cid == None:
        pass
        # gotta cash it!!!
    return cid, inchi

def all_same(items):
    return all(x == items[0] for x in items)

def consistent_smile_to_inchi(molecule_input):
    # this is to make sure that we get the same inchi regardless of the different parameters
    inchiS1 = []
    inchiS1.append(smile2standard_inchi(molecule_input))
    inchiS1.append(smile2standard_inchi2(molecule_input))
    inchiS1.append(smile2standard_inchi3(molecule_input))
    inchiS1.append(smile2standard_inchi4(molecule_input))
    inchiS1.append(smile2standard_inchi5(molecule_input))
    inchiS1.append(smile2standard_inchi6(molecule_input))
    inchiS1.append(smile2standard_inchi7(molecule_input))
    if all_same(inchiS1):
        return inchiS1[0]

def formation_energy(listMolecules, cc, pH=7, I=0.1, T=298.15, format_molecule='kegg'):

    gibbs_formation_transformed_list = []
    compound_list = []
    cid_list = []

    for molecule_input in listMolecules:

        # based on the input format (KEGG, MOL, or SMILEs) this section creates a compound object.
        if format_molecule == 'kegg':
            comp = Compound.from_kegg(molecule_input)

        elif format_molecule == 'smile':
            #corrected_smiles = Compound.smiles2smiles(molecule_input)

            # direct smile to inchi and database check
            inchi = consistent_smile_to_inchi(molecule_input)
            cid, Name = check_if_already_exists_inchi(inchi)
            if cid == None:
                # smile to mol to inchi and then check database
                mol = smile2mol(molecule_input)
                cid, inchi = mol_to_inchi_to_cid(mol)
            comp = Compound.from_inchi('KEGG', cid, inchi)


        elif format_molecule == 'mol':

            cid, inchi = mol_to_inchi_to_cid(molecule_input)
            comp = Compound.from_inchi('KEGG', cid, inchi)

        print comp

        major_ms_dG0_f = cc.get_major_ms_dG0_f(cid)

        # get the deltaGf transformed of the first species. uses the pKa's calculated previously. if interested in the
        # other species, uncomment below

        # for d in comp.get_species(major_ms_dG0_f, default_T):
        # print d['dG0_f']
        a = comp.get_species(major_ms_dG0_f, default_T)
        dG0_f = a.next()['dG0_f']

        # calculate the difference between the transformed formation energy of molecule and the first major species
        # Compound.transform_pH7(comp, 7, .1, 298.15)
        # Compound.transform_neutral(comp, 7, .1, 298.15)  # ???
        dG0f_dG0tf = Compound.get_transform(comp, pH, I, T)

        # total transformed gibbs formation energy
        Gibbs_trans_formation = dG0_f + dG0f_dG0tf

        # output everything
        gibbs_formation_transformed_list.append(Gibbs_trans_formation)
        compound_list.append(comp)
        cid_list.append(cid)

    return gibbs_formation_transformed_list, compound_list, cid_list