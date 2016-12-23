from component_contribution.compound import Compound
from component_contribution.thermodynamic_constants import default_T
from component_contribution import inchi2gv
import openbabel, numpy

def check_if_already_exists_inchi(inchi):

    with open('../data/compounds.tsv', 'r') as tsv:
        AoA = [line.strip().split('\t') for line in tsv]

    for i in numpy.arange(1,len(AoA)):
        if AoA[i].__len__()==3:
            found = AoA[i][2].find(inchi)
            if found != -1:
                cid = AoA[i][0]
                name = AoA[i][1]
                return cid, name
    return None, None

def anyFormat2standard_inchi(s,format):
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats(format, 'inchi')
    # conv.AddOption("F", conv.OUTOPTIONS)
    conv.AddOption("T", conv.OUTOPTIONS)
    conv.AddOption("x", conv.OUTOPTIONS, "noiso")
    conv.AddOption("w", conv.OUTOPTIONS)
    obmol = openbabel.OBMol()
    if not conv.ReadString(obmol, str(s)):
        return None
    inchi = conv.WriteString(obmol, True)  # second argument is trimWhitespace
    if inchi == '':
        return None
    else:
        return inchi

def anyFormat2standard_inchi2(s,format):
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats(format, 'inchi')
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

def anyFormat2standard_inchi3(s,format):
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats(format, 'inchi')
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

def anyFormat2standard_inchi4(s,format):
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats(format, 'inchi')
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

def anyFormat2standard_inchi5(s,format):
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats(format, 'inchi')
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

def anyFormat2standard_inchi6(s,format):
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats(format, 'inchi')
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

def anyFormat2standard_inchi7(s,format):
    openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats(format, 'inchi')
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

def all_same(items):
    return all(x == items[0] for x in items)

def consistent_to_inchi(molecule_input, format):
    # this is to make sure that we get the same inchi regardless of the different parameters
    inchiS1 = []
    inchiS1.append(anyFormat2standard_inchi(molecule_input, format))
    inchiS1.append(anyFormat2standard_inchi2(molecule_input, format))
    inchiS1.append(anyFormat2standard_inchi3(molecule_input, format))
    inchiS1.append(anyFormat2standard_inchi4(molecule_input, format))
    inchiS1.append(anyFormat2standard_inchi5(molecule_input, format))
    inchiS1.append(anyFormat2standard_inchi6(molecule_input, format))
    inchiS1.append(anyFormat2standard_inchi7(molecule_input, format))
    if all_same(inchiS1):
        return inchiS1[0]
    else:
        raise ValueError("inchi changes depending on parameters! need to pick one or try them all?")

def convert2any(s, formatin, formatout):
    #openbabel.obErrorLog.SetOutputLevel(-1)

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats(formatin, formatout)
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

def mol_to_inchi_to_cid(mol):
    inchi = Compound.mol2inchi(mol)
    inchiS = mol2standard_inchi(mol)
    cid, Name = check_if_already_exists_inchi(inchiS)
    if cid == None:
        pass
        # gotta cash it!!!
    return cid, inchi

def decompose_without_cache(cc,comp):

    try:
        group_vec = cc.decomposer.smiles_to_groupvec(comp.smiles_pH7)
        g = numpy.matrix(group_vec.ToArray())
        dG0_gc = cc.params['dG0_gc'][0:cc.Ng, :]
        return float(numpy.dot(g, dG0_gc))
    except inchi2gv.GroupDecompositionError:
        return numpy.nan

def formation_energy(listMolecules, cc, pH=7, I=0.1, T=298.15, format_molecule='kegg'):

    gibbs_formation_transformed_list = []
    compound_list = []
    cid_list = []

    for molecule_input in listMolecules:

        # based on the input format (KEGG, MOL, or SMILEs) this section creates a compound object.
        if format_molecule == 'kegg':
            inchiS = Compound.get_inchi(molecule_input)
            cid = molecule_input

        elif format_molecule == 'smiles':
            #corrected_smiles = Compound.smiles2smiles(molecule_input)

            # direct smile to inchi and database check
            inchiS = consistent_to_inchi(molecule_input, 'smiles')
            cid, Name = check_if_already_exists_inchi(inchiS)
            # if this consistently doesnt work try changing to mol first: mol=convert2any(molecule_input,'smiles','mol')

        elif format_molecule == 'mol':
            inchiS = consistent_to_inchi(molecule_input, 'mol')
            cid, Name = check_if_already_exists_inchi(inchiS)
        else:
            raise ValueError('you did not input one of the valid format for the metabolites: \'kegg\', \'smiles\', or \'mol\'')

        # this only gets the pKa information of the molecule
        comp = Compound.from_inchi('KEGG', cid, inchiS)
        print comp

        # get the thermodynamic information of the molecule
        if cid == None:
            major_ms_dG0_f = decompose_without_cache(cc, comp)
        else:
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