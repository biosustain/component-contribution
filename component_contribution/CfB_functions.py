from component_contribution.compound import Compound
from component_contribution.thermodynamic_constants import default_T
from component_contribution import inchi2gv
from component_contribution.kegg_reaction import KeggReaction
import openbabel, numpy, csv, json, os.path
from component_contribution.Chebi_web import ChEBI


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

def load_mnx_file(filename = '../data/chem_prop'):

    with open(filename, mode='r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        mydict = {rows[0]: rows[5] for rows in reader}
    return mydict

def _decompose_bigg_reaction(cc, reaction, bigg_dict):
    if cc.params is None:
        cc.train()

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

def separate_S_matrix(S,input_ids):

    '''
    Component contribution as it is written gets all of its information from a cache file that matches up a "kegg id"
    which means that if we can conveniently use this ID then we should.

    Conversely, in reality all of the group decomposition is based on smiles. if we can not retrieve the kegg id then we
    shall do group decomposition with the smiles.

    :param S:
    :param input_ids:
    :return:
    '''

    # first load bigg database to translate. bigg dict is the downloaded bigg metabolite database
    bigg_dict = load_bigg_dict(filename='../data/bigg_models_metabolites.tsv')

    #create variables to save output
    to_kegg = []
    to_kegg_dict = {}
    non_kegg =[]
    non_kegg_refs ={}

    # map all ids I can to kegg.
    for compound_id in input_ids:
        # classifying met. first see if id matches the bigg database
        all_references_json = bigg_dict.get(compound_id, 'its not bigg')
        print compound_id

        # if it didn't map, it's not bigg, and we assume it's a mol file
        if all_references_json == 'its not bigg':

            if os.path.isfile('../examples/molfiles/' + compound_id):
                # read mol file
                test_file = open('../examples/molfiles/' + compound_id, 'r')
                mol = test_file.read()
                inchiS = consistent_to_inchi(mol, 'mol')
                cid, Name = check_if_already_exists_inchi(inchiS)
                if cid == None:
                    comp = Compound.from_inchi_with_keggID('MOL', cid, inchiS)
                    # save the references
                    non_kegg_refs[compound_id] = comp
                    non_kegg.append(compound_id)
                else:
                    kegg_id = cid
                    to_kegg_dict[compound_id] = kegg_id
                    to_kegg.append(compound_id)
            else: raise ValueError


        # we got a hit so it's a bigg metabolite, lets get either the kegg id or the structure
        else:
            all_references_readable = json.loads(all_references_json)

            if 'KEGG Compound' in all_references_readable: # we matched it to kegg id and will use component contribution's database for this guy
                kegg_reference = all_references_readable['KEGG Compound']
                kegg_id = kegg_reference[0]['id']
                to_kegg_dict[compound_id] = kegg_id
                to_kegg.append(compound_id)

            else: # bigg id with no kegg id!!!
                # get compound info... currently only checks if there is chebi info.
                comp_data = get_compound_info(all_references_readable)
                # save the references
                non_kegg_refs[compound_id] = comp_data
                non_kegg.append(compound_id)


    # # now begin to separate the S matrix according to category
    # S_kegg = S[[input_ids.index(b) for b in to_kegg], :]
    # S_bigg_non_kegg = S[[input_ids.index(b) for b in bigg_non_kegg],:]
    # S_non_bigg = S[[input_ids.index(b) for b in non_bigg],:]

    # new cids list
    output_ids = [to_kegg_dict.get(id, id) for id in input_ids]

    return { 'to_kegg':to_kegg,
             'to_kegg_dict':to_kegg_dict,
             #'S_kegg': S_kegg,

             'non_kegg':non_kegg,
             'non_kegg_refs':non_kegg_refs,
             #'S_bigg_non_kegg':S_bigg_non_kegg,

             'output_ids':output_ids}


def only_decompose(cc, reactions):
    """
        Arguments:
            reaction - a KeggReaction object

        Returns:
            the CC estimation for this reaction's untransformed dG0 (i.e.
            using the major MS at pH 7 for each of the reactants)

            X and G are the decompositions of the reaction into reactions and groups respectively
            eq.10 in the paper
    """
    X = []
    G = []

    for reaction in reactions:
        try:
            x, g = cc._decompose_reaction(reaction)
        except inchi2gv.GroupDecompositionError:
            x = numpy.zeros((cc.Nc, 1))
            g = numpy.zeros((cc.params['G'].shape[1], 1))
        X.append(list(x.flat))
        G.append(list(g.flat))
    X = numpy.matrix(X).T
    G = numpy.matrix(G).T
    return X, G

def get_compound_info(info):
    '''
    this function retrieves the structural information included in info
    :param info: this is a dictionary including alternative references for a bigg id
    :return:

    TODO:
        1 define an order of databases
        2 determine what to do when there are multiple compounds, for now taking first
        3 what to return when there is nothing
    '''
    # take ChEBI
    if 'CHEBI' in info:
        ok = ChEBI()
        res = ok.getCompleteEntity(info['CHEBI'][0]['id'])
        inchi = res.inchi
        compound_data = Compound.from_inchi_with_keggID('CHEBI', 'xly', inchi)

    else:
        compound_data = None

    return compound_data

def add_thermo_comp_info(self, cc):
    # check that all CIDs in the reaction are already cached by CC
    Nc, Nr = self.S.shape
    reactions = []
    for j in xrange(Nr):
        sparse = {self.cids[i]: self.S[i, j] for i in xrange(Nc)
                  if self.S[i, j] != 0}
        reaction = KeggReaction(sparse)
        reactions.append(reaction)

    self.dG0, self.cov_dG0 = cc.get_dG0_r_multi(reactions,self.separated_by_id_type)