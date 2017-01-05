from component_contribution.compound import Compound
from component_contribution.thermodynamic_constants import default_T
from component_contribution import inchi2gv
import openbabel, numpy, csv, json, copy

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
            inchiS = Compound.get_inchi_from_kegg(molecule_input)
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
        comp = Compound.from_inchi_with_keggID('KEGG', cid, inchiS)
        print comp

        # get the thermodynamic information of the molecule # get the deltaGf transformed of the first species.
        if cid == None:
            major_ms_dG0_f = decompose_without_cache(cc, comp)
        else:
            major_ms_dG0_f = cc.get_major_ms_dG0_f(cid)

        # uses the pKa's calculated previously. if interested in the
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

def load_bigg_dict(filename = '../data/bigg_models_metabolites.tsv'):

    with open(filename, mode='r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        mydict = {rows[0]: rows[4] for rows in reader}
    return mydict

def _decompose_bigg_reaction(cc, reaction, bigg_dict):
    if cc.params is None:
        cc.train()

    # copy the reaction because we will be modifying it
    reaction_sparese_with_kegg =  copy.deepcopy(reaction.sparse)
    only_kegg = True

    # dict to replace the id's in the model
    model_ids_to_replace = {}

    # map all ids I can to kegg. bigg dict is the downloaded bigg metabolite database
    for compound_id, coeff in reaction.iteritems():
        all_references_json = bigg_dict[compound_id]
        all_references_readable = json.loads(all_references_json)
        kegg_reference = all_references_readable.get('KEGG Compound', None)

        if kegg_reference != None:
            kegg_id = kegg_reference[0]['id']
            reaction_sparese_with_kegg[kegg_id] = reaction_sparese_with_kegg.pop(compound_id)
            model_ids_to_replace[compound_id] = kegg_id

        else:
            # there were some non-kegg mets!
            only_kegg = False
            # what to do if there is at least one non-kegg item...?

            # model_ids_to_replace

            raise(ValueError)

    if only_kegg:
        reaction.sparse = copy.deepcopy(reaction_sparese_with_kegg)
        reaction.format = 'kegg'

    # load all cid's in the database, and G is all the groups in the cids_joined (reactants used to train)
    cids = list(cc.params['cids'])
    G = cc.params['G']

    # calculate the reaction stoichiometric vector and the group incidence
    # vector (x and g)
    x = numpy.matrix(numpy.zeros((cc.Nc, 1)))
    x_prime = []
    G_prime = []

    for compound_id, coeff in reaction.iteritems():
        if compound_id in cc.cids_joined: # cids_joined is the list of cid used in the training data
            i = cids.index(compound_id)
            x[i, 0] = coeff
        else:
            # Decompose the compound and calculate the 'formation energy'
            # using the group contributions.
            # Note that the length of the group contribution vector we get
            # from CC is longer than the number of groups in "groups_data"
            # since we artifically added fictive groups to represent all the
            # non-decomposable compounds. Therefore, we truncate the
            # dG0_gc vector since here we only use GC for compounds which
            # are not in cids_joined anyway.
            x_prime.append(coeff)
            comp = cc.ccache.get_compound(compound_id)
            group_vec = cc.decomposer.smiles_to_groupvec(comp.smiles_pH7)
            G_prime.append(group_vec.ToArray())

    if x_prime != []:
        g = numpy.matrix(x_prime) * numpy.vstack(G_prime)
    else:
        g = numpy.matrix(numpy.zeros((1, 1)))

    # g1=copy.deepcopy(g)
    g.resize((G.shape[1], 1), refcheck=False)# referencing is very weird in numpy?

    return x, g, model_ids_to_replace

def replace_ids_with_cids(model):

    new_cids = [model.model_ids_to_replace.get(item, item) for item in model.cids]
    model.cids = new_cids
    return model