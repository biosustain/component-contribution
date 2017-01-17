from component_contribution.compound import Compound
from component_contribution.thermodynamic_constants import default_T
from component_contribution import inchi2gv
from component_contribution.kegg_reaction import KeggReaction
from component_contribution.chemaxon import getNonPolarArea
import openbabel, numpy, csv, json, os.path, warnings, urllib
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
    conv.AddOption("--gen3d", conv.OUTOPTIONS)
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
    return cid, inchiS

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

def process_input_mets(input_metabolites, ccache):

    '''
    Component contribution as it is written gets all of its information from a cache file that matches up a "kegg met"
    which means that if we can conveniently use this ID then we should.

    Conversely, in reality all of the group decomposition is based on smiles. if we can not retrieve the kegg met then we
    shall do group decomposition with the smiles.

    :param input_metabolites:
    :return:
    '''

    # first load bigg database to translate. bigg dict is the downloaded bigg metabolite database
    bigg_dict = load_bigg_dict(filename='../data/bigg_models_metabolites.tsv')

    #create variables to save output
    met_to_comp_dict = {}
    non_kegg =[] # this one is used  /all_model.py", line 314, in _get_transform_ddG0
    all_comp_data ={}

    # map all ids I can to kegg (or whatever internal ID with Compound info).
    for met in input_metabolites:
        # classifying met. first see if met matches the bigg database
        all_references_json = bigg_dict.get(met, 'its not bigg')

        # it didn't match bigg, so we assume it's a mol file
        if all_references_json == 'its not bigg':

            if os.path.isfile('../examples/molfiles/' + met):
                # read mol file
                test_file = open('../examples/molfiles/' + met, 'r')
                mol = test_file.read()
                inchiS = consistent_to_inchi(mol, 'mol')
                cid, Name = check_if_already_exists_inchi(inchiS)
                if cid == None:
                    comp = Compound.from_inchi_with_keggID('MOL', met, inchiS)
                    comp.molfile = mol
                    # save the references
                    met_to_comp_dict[met] = met
                    all_comp_data[met] = comp
                    non_kegg.append(met)
                    print( met+ ' was inserted as mol file. it had no id mapping to the internal reference. will use input file name')
                else:
                    comp = ccache.get_compound(cid)
                    comp.molfile = mol
                    all_comp_data[cid] = comp
                    met_to_comp_dict[met] = cid
                    print( met+ ' was inserted as mol file. it mapped to the kegg id' + cid)
            else: raise ValueError


        # we got a hit so it's a bigg metabolite, lets get either the kegg met or the structure
        else:
            all_references_readable = json.loads(all_references_json)

            if 'KEGG Compound' in all_references_readable: # we matched it to kegg met and will use component contribution's database for this guy
                kegg_reference = all_references_readable['KEGG Compound']
                kegg_id = kegg_reference[0]['id']
                met_to_comp_dict[met] = kegg_id
                try:
                    comp = ccache.get_compound(kegg_id)
                    # s_mol = urllib.urlopen('http://rest.kegg.jp/get/cpd:%s/mol' % kegg_id).read()
                    # comp.molfile = s_mol
                    all_comp_data[kegg_id] = comp
                    print( met+ ' mapped to the kegg id' + kegg_id)

                except:
                    print(met + ' had kegg id in the bigg database but wasn\'t in component contribution database')
                    #comp = get_info_by_mapping(ccache, all_references_readable, met)
                    comp = None

            else: # bigg met with no kegg met!!!
                # get compound info... currently only checks if there is chebi info.

                # comp = get_info_by_mapping(ccache, all_references_readable, met)
                # all_comp_data[met] = comp
                # non_kegg.append(met)
                print(met + ' had no kegg id in the bigg database. it could potentially be found using other databases, but ignoring for now')
                comp = None

    # # now begin to separate the S matrix according to category
    # S_kegg = S[[input_metabolites.index(b) for b in to_kegg], :]
    # S_bigg_non_kegg = S[[input_metabolites.index(b) for b in bigg_non_kegg],:]
    # S_non_bigg = S[[input_metabolites.index(b) for b in non_bigg],:]

    # new cids list
    output_ids = [met_to_comp_dict.get(met, met) for met in input_metabolites]

    return {'met_to_comp_dict':met_to_comp_dict,
             'non_kegg':non_kegg,
             'all_comp_data':all_comp_data,
             'output_ids':output_ids}

def get_info_by_mapping(ccache, all_references_readable,met):
    comp = get_compound_info(ccache, all_references_readable)
    if comp == None:
        print(met + ' had no id mapping to the internal reference OR it mapped, but in has no cached comp info')
    else:
        print(met + 'had compound information found from mapping databases')
    return comp

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


def add_thermo_comp_info(self, cc):
    # check that all CIDs in the reaction are already cached by CC
    Nc, Nr = self.S.shape
    reactions = []
    for j in xrange(Nr):
        if j%50==0: print('creating reactions to decompose: '+ str(j))
        sparse = {self.cids[i]: self.S[i, j] for i in xrange(Nc)
                  if self.S[i, j] != 0}
        reaction = KeggReaction(sparse)
        reactions.append(reaction)

    self.dG0, self.cov_dG0 = cc.get_dG0_r_multi(reactions,self.comp_data)

def addHydrogens(input_mol):
    import openbabel

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "pdb")

    mol_object = openbabel.OBMol()
    obConversion.ReadString(mol_object, input_mol)
    print mol_object.NumAtoms()

    mol_object.AddHydrogens()
    print mol_object.NumAtoms()

    return obConversion.WriteString(mol_object)

def calc_concentration(calculate_conc, input_metabolites, comp_data=None, ccache=None):
    '''
    input is the metabolites in the stoichiometric matrix of the thermo model. These are already converted to keggIds
    :return:
    '''
    import math
    from component_contribution.molecule import Molecule
    from component_contribution.compound_cacher import KeggCompoundCacher
    from component_contribution.CfB_functions import process_input_mets


    default_C = 0.001
    conc_M =  numpy.matrix(numpy.full((len(input_metabolites), 1), default_C))

    input_metabolites = list(set(input_metabolites))
    # hydrogen = [i for i in input_metabolites if i.startswith('h_') or i.startswith('h2o_')]
    # water = [i for i in input_metabolites if i.startswith('h2o_')]
    # input_metabolites = [i for i in input_metabolites if not i.startswith('h_') and not i.startswith('h2o_')]

    if ccache == None:
        ccache = KeggCompoundCacher()

    if comp_data == None:
        comp_data = process_input_mets(input_metabolites, ccache)

    comp_data_mapper = comp_data['all_comp_data']

    # calc charge and NPSA, and [] for all comps
    countcharge = {}
    NPSA = {}
    concentration = {}

    for compound_id,comp in zip(comp_data_mapper.keys(), comp_data_mapper.values()):
        print compound_id

        if comp == None:
            countcharge[compound_id] = None
            NPSA[compound_id] = None
            concentration[compound_id] = default_C
            continue

        # water or H+: putting H as "1" because its log(1)=0 and it will not affect reaction energy
        if compound_id in ('C00001','C00080'):
            countcharge[compound_id] = None
            NPSA[compound_id] = None
            concentration[compound_id] = 1
            continue

        # glutamate or glutamine
        if compound_id in ('C00064', 'C00025'):
            countcharge[compound_id] = None
            NPSA[compound_id] = None
            concentration[compound_id] = 0.1
            continue

        if calculate_conc:
            # calculate charge ##

            # charge calculations work better from smiles
            if comp.smiles_pH7:
                smiles_pH7 = comp.smiles_pH7
                mol = Molecule.FromSmiles(smiles_pH7)
            else:
                countcharge[compound_id] = None
                NPSA[compound_id] = None
                concentration[compound_id] = default_C
                continue

            charge_list = mol.GetAtomCharges()
            charge_count = sum(x != 0 for x in charge_list)
            countcharge[compound_id] = charge_count

            # calc NPSA ##
            NPSA1 = getNonPolarArea(comp.inchi, pH=7)
            NPSA[compound_id] = NPSA1

            # old way to calc NPSA
            # # NPSA calculations work better from mol files
            # molfile = comp.molfile
            # mol = Molecule.FromMol(molfile)
            #
            # # create the pdb file. in the futre can cache this
            # test = Molecule._ToFormat(mol.obmol,'pdb')
            # file_name = "../examples/pdb_files/" + input_id +".pdb" # have to decide what we're going to save all the files as
            # with open(file_name, "w") as text_file:
            #     text_file.write(test)
            #
            # # load pdb create " soup".     A Soup is essentially a collection of atoms, which we can grab by:
            # soup = pdbatoms.Soup(file_name)
            # atoms = soup.atoms()
            #
            # # asa - calculates the accessible surface - area of every atom in list of atom, with respect to the other atoms. which  assigns the asa to to each atom.asa
            # pdbatoms.add_radii(atoms) # this calculates the radius of each atom.
            # areas = asa.calculate_asa(atoms, 1.4) # 1.4 is the probe i.e. the size of water, changing this can change results ~ 25%
            # total_area = sum(areas)
            #
            # # get atom neighbors
            # adj = asa.adjacency_list(atoms, 1.8)
            # adj = [[atoms[c].element for c in at] for at in adj]
            #
            # # get polar surface area, i.e. the area contributed by polar atoms only (oxygen, nitrogen and the hydrogen atoms attached to them
            # polar_area=0
            # for a, area, adj1 in zip(atoms, areas,adj):
            #     print a, area
            #     if a.element in ['O','N']:       # a.bfactor = area
            #         polar_area += area
            #     if a.element=='H' and  any([i in ['O','N'] for i in adj1]):
            #         polar_area += area
            #
            # NPSA1 = total_area - polar_area

            conc = math.exp(charge_count*1.0425-NPSA1*0.0272)/1000

            concentration[compound_id] = conc
        else:
            countcharge[compound_id] = None
            NPSA[compound_id] = None
            concentration[compound_id] = default_C

    # 1. compounds can have multiple input ids associated.
    # 2. some input ids have no compounds associated leave them with default concentrations
    for i, met in enumerate(input_metabolites):
        comp_id = comp_data['met_to_comp_dict'].get(met, None)
        if comp_id == None: continue
        conc_M[i] = concentration[comp_id]

    return conc_M, concentration, NPSA, countcharge

def remove_duplicate(input_list):
    repeat = []
    uniq = []
    for x in input_list:
        if x not in uniq:
            uniq.append(x)
            # seen.add(x)
        else:
            repeat.append(x)

def get_compound_info(ccahce,info):
    '''
    this function retrieves the structural information for the references included in the bigg database. this function
    should no longer be necessary once a proper internal database is generated.

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
        try:
            if res.ChemicalStructures[0].type == 'mol':
                molfile = res.ChemicalStructures[0].structure
                cid,inchi = mol_to_inchi_to_cid(molfile)
                try: compound = ccahce.get_compound(cid)
                except:
                    compound = Compound.from_inchi_with_keggID('CHEBI', met, inchi)
                compound.molfile = molfile
            else:
                compound = None

        except:
            try:
                inchi = res.inchi
                compound = Compound.from_inchi_with_keggID('CHEBI', met, inchi)
            except:
                compound = None
    else:
        compound = None
    return compound
