import logging
import csv
import re
import platform
import io
from subprocess import Popen, PIPE
import openbabel

from component_contribution.exceptions import ChemAxonError


logger = logging.getLogger(__name__)

CXCALC_BIN = 'cxcalc'
if platform.system() == 'Windows':
    use_shell_for_echo = True
else:
    use_shell_for_echo = False

MID_PH = 7.0
N_PKAS = 20
_obElements = openbabel.OBElementTable()


def run_cxcalc(molstring, args):
    with open(platform.DEV_NULL, 'w') as dev_null:
        try:
            logger.debug("INPUT: echo %s | %s" % (molstring, ' '.join([CXCALC_BIN] + args)))
            p1 = Popen(["echo", molstring], stdout=PIPE, shell=use_shell_for_echo)
            p2 = Popen([CXCALC_BIN] + args, stdin=p1.stdout,
                       executable=CXCALC_BIN, stdout=PIPE, stderr=dev_null, shell=False)
            res = p2.communicate()[0]
            if p2.returncode != 0:
                raise ChemAxonError(str(args))
            logger.debug("OUTPUT: %s" % res)
            return res.decode("utf-8")
        except OSError:
            raise Exception("Marvin (by ChemAxon) must be installed to calculate pKa data.")


def parse_pka_output(s, n_acidic, n_basic):
    """
    Returns
    -------
    atom2pka : dict
        A dictionary that maps the atom index to a list of pKas that are assigned to that atom.
    """
    atom2pka = {}

    pkaline = s.split('\n')[1]
    splitline = pkaline.split('\t')
    splitline.pop(0)
    
    if n_acidic + n_basic > 0:
        if len(splitline) != (n_acidic + n_basic + 2):
            raise ChemAxonError('ChemAxon failed to find any pKas')
        
        pka_list = []
        acid_or_base_list = []
        for i in range(n_acidic + n_basic):
            x = splitline.pop(0) 
            if x == '':
                continue
            
            pka_list.append(float(x))
            if i < n_acidic:
                acid_or_base_list.append('acid')
            else:
                acid_or_base_list.append('base')
        
        atom_list = splitline.pop(0)

        if len(atom_list) > 0:  # a comma separated list of the deprotonated atoms
            atom_numbers = [int(y)-1 for y in atom_list.split(',')]
            for i, j in enumerate(atom_numbers):
                atom2pka.setdefault(j, [])
                atom2pka[j].append((pka_list[i], acid_or_base_list[i]))
    
    smiles_list = splitline
    return atom2pka, smiles_list


def _dissociation_constants(molstring, n_acidic=N_PKAS, n_basic=N_PKAS, ph=MID_PH):
    """
    Returns
    -------
    pka_list : list
        A list is of the pKa values in ascending order.
    pseudoisomers : str
        The pseudoisomer SMILES string at the given pH.
    """
    args = []
    if n_acidic + n_basic > 0:
        args += ['pka', '-a', str(n_acidic), '-b', str(n_basic),
                 'majorms', '-M', 'true', '--pH', str(ph)]
    
    output = run_cxcalc(molstring, args)
    atom2pka, smiles_list = parse_pka_output(output, n_acidic, n_basic)
    
    all_pkas = []
    for pKa_list in atom2pka.values():
        all_pkas += [pKa for pKa, _ in pKa_list]
    
    return sorted(all_pkas), smiles_list


def dissociation_constants(molstring, n_acidic=N_PKAS, n_basic=N_PKAS, ph=MID_PH):
    """
    Arguments
    ---------
    molstring : str
        A text description of the molecule (SMILES or InChI).
    n_acidic : int
        The max number of acidic pKas to calculate.
    n_basic : int
        The max no. of basic pKas to calculate.
    ph : float
        The pH for which the major pseudoisomer is calculated.

    Returns
    -------
    all_pkas : list
        A list of floats (pKa values)
    major_ms : str
        SMILES string of the major pseudoisomer at pH
    """
    all_pkas, smiles_list = _dissociation_constants(molstring, n_acidic, n_basic, ph)
    major_ms = smiles_list[0]
    return all_pkas, major_ms


def formula_and_charge(molstring) -> tuple:
    """
    Arguments
    ---------
    molstring : str
        A text description of the molecule (SMILES or InChI).

    Returns
    -------
    formula : str
        Chemical formula of the molecule.
    formal_charge : int
        The molecule charge.
    """
    args = ['formula', 'formalcharge']
    output = run_cxcalc(molstring, args)
    # the output is a tab separated table whose columns are:
    # id, Formula, Formal charge
    f = io.StringIO(output)
    tsv_output = csv.reader(f, delimiter='\t')
    headers = next(tsv_output)
    if headers != ['id', 'Formula', 'Formal charge']:
        raise ChemAxonError('cannot get the formula and charge for: ' + molstring)
    _, formula, formal_charge = next(tsv_output)

    try:
        formal_charge = int(formal_charge)
    except ValueError:
        formal_charge = 0
    
    return formula, formal_charge


def atom_bag_and_charge(molstring):
    """
    Arguments
    ---------
    molstring : str
       A text description of the molecule (SMILES or InChI).
    """
    formula, formal_charge = formula_and_charge(molstring)

    atom_bag = {}
    for mol_formula_times in formula.split('.'):
        for times, mol_formula in re.findall('^(\d+)?(\w+)', mol_formula_times):
            if not times:
                times = 1
            else:
                times = int(times)
            for atom, count in re.findall("([A-Z][a-z]*)([0-9]*)", mol_formula):
                if count == '':
                    count = 1
                else:
                    count = int(count)
                atom_bag[atom] = atom_bag.get(atom, 0) + count * times
    
    n_protons = sum([c * _obElements.GetAtomicNum(str(elem)) for (elem, c) in atom_bag.items()])
    atom_bag['e-'] = n_protons - formal_charge

    return atom_bag, formal_charge


def non_polar_surface_area(molstring, pH=7):
    """

    Calculates the non-polar surface area.

    Arguments
    ---------
    molstring : str
        A text description of the molecule (SMILES or InChI).

    Returns
    -------
    float

    """
    args = ['msa', '--type', 'ASA_H', '--pH', str(pH)]

    output = run_cxcalc(molstring, args)
    npsa = float(re.split('\t|\n', output)[-2])

    return npsa


