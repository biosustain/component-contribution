import logging
import re

import openbabel
import numpy as np
import requests

from scipy.misc import logsumexp

import bioservices

from cobra.core.metabolite import Metabolite
from component_contribution.exceptions import OpenBabelError

from component_contribution import chemaxon
from component_contribution.chemaxon import non_polar_surface_area
from component_contribution.core.molecule import Molecule
from component_contribution.thermodynamics.constants import R
from component_contribution.thermodynamics.functions import debye_huckel


KEGG = bioservices.KEGG()
CHEBI = bioservices.ChEBI()
HMDB_URL = "http://www.hmdb.ca/structures/metabolites/%s.mol"

KEGG_CPD_RE = re.compile("^C\d{5}$")
CHEBI_CPD_RE = re.compile("^CHEBI:\d+$")
HMDB_CPD_RE = re.compile("^HMDB\d+$")

logger = logging.getLogger(__name__)


MIN_PH = 0.0
MAX_PH = 14.0


EXCEPTIONS = {
    # We add an exception for H+ (and put nH = 0) in order to eliminate its effect of the Legendre transform.
    'C00080': ({'H': 1}, [], None, 0, [0], [0]),
    # ChemAxon gets confused with the structure of sulfur # (returns a protonated form, [SH-], at pH 7).
    'C00087': ({'S': 1, 'e-': 16}, [], 'S', 0, [0], [0]),
    # ChemAxon gets confused with the structure of carbon monoxide # (returns a protonated form, [CH]#[O+], at pH 7).
    'C00237': ({'C': 1, 'O': 1, 'e-': 14}, [], '[C-]#[O+]', 0, [0], [0]),
    # ChemAxon gets confused with the structure of hydrogen
    'C00282': ({'H': 2, 'e-': 2}, [], None, 0, [2], [0]),
    # When given the structure of carbonic acid, ChemAxon returns the pKas for CO2(tot), i.e. it assumes the
    # non-hydrated CO2 species is one of the pseudoisomers, and the lower pKa value is 6.05 instead of 3.78.
    # Here, we introduce a new "KEGG" compound that will represent pure bicarbonate (without CO2(sp)) and therefore
    # plug in the pKa values from Alberty's book.
    'C01353': ({'C': 1, 'H': 1, 'O': 3, 'e-': 32}, [10.33, 3.43], 'OC(=O)[O-]', 1, [0, 1, 2], [-2, -1, 0]),
    # Metal Cations get multiple pKa values from ChemAxon, which is a bug. We override the important ones here:
    # Ca2+
    'C00076': ({'Ca': 1, 'e-': 18}, [], '[Ca++]', 0, [0], [2]),
    # K+
    'C00238': ({'K': 1, 'e-': 18}, [], '[K+]', 0, [0], [1]),
    # Mg2+
    'C00305': ({'Mg': 1, 'e-': 10}, [], '[Mg++]', 0, [0], [2]),
    # Fe2+
    'C14818': ({'Fe': 1, 'e-': 24}, [], '[Fe++]', 0, [0], [2]),
    # Fe3+
    'C14819': ({'Fe': 1, 'e-': 23}, [], '[Fe+++]', 0, [0], [3]),
    # ferredoxin(red)
    'C00138': ({'Fe': 1, 'e-': 26}, [], None, 0, [0], [0]),
    # ferredoxin(ox)
    'C00139': ({'Fe': 1, 'e-': 25}, [], None, 0, [0], [1])
}


class Compound(Metabolite):
    """
    Representation of a chemical compound.

    Attributes
    ----------
    database : str
        The database where the compound was obtained.
    compound_id : str
        The identifier of the compound in the database
    inchi : str
        The InChI representation of the chemical
    atom_bag : ?
        ?????
    p_kas : list

    smiles_ph_7 : str
        SMILES for the major microspicies at pH7
    major_ms_ph_7 : ?
        ????
    number_of_protons : int
        The number of protons
    charges : list
        Charges (zs?)

    """
    
    def __init__(self, database, compound_id, inchi, inchi_key, atom_bag, p_kas, smiles_ph_7,
                 major_ms_ph_7, number_of_protons, charges, npsa=None, molfile=None):
        super(Compound, self).__init__(inchi_key or compound_id)
        self.database = database
        self.compound_id = compound_id
        self.inchi = inchi
        self.inchi_key = inchi_key
        self.atom_bag = atom_bag
        self._npsa = npsa
        self.p_kas = p_kas
        self.smiles_ph_7 = smiles_ph_7
        self.major_ms_ph_7 = major_ms_ph_7
        self.number_of_protons = number_of_protons
        self.charges = charges
        self.molfile = molfile

    def __repr__(self):
        return "<Compound (%s:%s) %s>" % (self.database, self.compound_id, self.smiles_ph_7)

    @classmethod
    def from_database(cls, compound_id):
        if KEGG_CPD_RE.match(compound_id):
            return cls.from_kegg(compound_id)
        elif CHEBI_CPD_RE.match(compound_id):
            return cls.from_chebi(compound_id)
        elif HMDB_CPD_RE.match(compound_id):
            return cls.from_hmdb(compound_id)
        elif compound_id.startswith("InChI="):
            return cls.from_inchi(compound_id, "UNKNOWN", cls.inchi2inchi_key(compound_id))
        else:
            return cls("UNKNOWN", compound_id, None, None, {}, [], None, 0, [0], [0])

    @classmethod
    def from_kegg(cls, compound_id):
        inchi = cls.get_inchi_from_kegg(compound_id)
        return cls.from_inchi(inchi, "KEGG", compound_id)

    @classmethod
    def from_chebi(cls, compound_id):
        inchi = cls.get_inchi_from_chebi(compound_id)
        return cls.from_inchi(inchi, "CHEBI", compound_id)

    @classmethod
    def from_hmdb(cls, compound_id):
        inchi = cls.get_inchi_from_hmdb(compound_id)
        return cls.from_inchi(inchi, "HMDB", compound_id)

    @classmethod
    def from_inchi(cls, inchi, database, compound_id):
        if compound_id in EXCEPTIONS:
            inchi_key = cls.inchi2inchi_key(inchi)
            return cls(database, compound_id, inchi, inchi_key, *EXCEPTIONS[compound_id])
        # If the compound has no explicit structure, we assume that it has
        # no proton dissociation in the relevant pH range
        elif inchi is None:
            return cls(database, compound_id, inchi, None, {}, [], None, 0, [0], [0])

        # Otherwise, we use ChemAxon's software to get the pKas and the 
        # properties of all microspecies

        try:
            p_kas, major_ms_smiles = chemaxon.dissociation_constants(inchi)
            major_ms_smiles = Compound.smiles2smiles(major_ms_smiles)
            p_kas = sorted([pka for pka in p_kas if MIN_PH < pka < MAX_PH], reverse=True)
        except chemaxon.ChemAxonError:
            logging.warning('chemaxon failed to find pKas for this molecule: ' + inchi)
            # use the original InChI to get the parameters (i.e. assume it 
            # represents the major microspecies at pH 7)
            major_ms_smiles = Compound.inchi2smiles(inchi)
            p_kas = []
        
        if major_ms_smiles:
            atom_bag, major_ms_charge = chemaxon.atom_bag_and_charge(major_ms_smiles)
            major_ms_number_of_protons = atom_bag.get('H', 0)
        else:
            atom_bag = {}
            major_ms_charge = 0
            major_ms_number_of_protons = 0

        n_species = len(p_kas) + 1
        if len(p_kas) == 0:
            major_ms_ph7 = 0
        else:
            major_ms_ph7 = len([1 for pka in p_kas if pka > 7])
            
        number_of_protons = []
        zs = []

        for i in range(n_species):
            zs.append((i - major_ms_ph7) + major_ms_charge)
            number_of_protons.append((i - major_ms_ph7) + major_ms_number_of_protons)

        inchi_key = cls.inchi2inchi_key(inchi)
        return cls(database, compound_id, inchi, inchi_key, atom_bag, p_kas,
                   major_ms_smiles, major_ms_ph7, number_of_protons, zs)

    def to_json_dict(self):
        return {'database': self.database,
                'compound_id': self.compound_id,
                'inchi': self.inchi,
                'inchi_key': self.inchi_key,
                'atom_bag': self.atom_bag,
                'p_kas': self.p_kas,
                'smiles_ph_7': self.smiles_ph_7,
                'major_ms_ph_7': self.major_ms_ph_7,
                'number_of_protons': self.number_of_protons,
                'zs': self.charges,
                'npsa': self._npsa}
    
    @classmethod
    def from_json_dict(cls, d):
        inchi = d['inchi']
        inchi_key = d.get('inchi_key')
        if inchi is not None and inchi_key is None:  # Legacy data has no inchi_key
            inchi_key = cls.inchi2inchi_key(inchi)

        return cls(d['database'], d['compound_id'], inchi, inchi_key, d['atom_bag'], d['p_kas'],
                   d['smiles_ph_7'], d['major_ms_ph_7'], d['number_of_protons'], d['zs'], npsa=d.get('npsa', None))

    @classmethod
    def get_inchi_from_kegg(cls, compound_id):
        s_mol = KEGG.get('cpd:%s' % compound_id, 'mol')
        return cls.mol2inchi(s_mol)

    @staticmethod
    def get_inchi_from_chebi(compound_id):
        entity = CHEBI.getCompleteEntity(compound_id)
        if hasattr(entity, 'inchi'):
            return entity.inchi

    @classmethod
    def get_inchi_from_hmdb(cls, compound_id):
        s_mol = requests.get(HMDB_URL % compound_id)
        return cls.mol2inchi(s_mol)

    @staticmethod
    def inchi2inchi_key(inchi):
        openbabel.obErrorLog.SetOutputLevel(-1)

        converter = openbabel.OBConversion()
        converter.SetInAndOutFormats('inchi', 'inchikey')
        converter.AddOption("x", converter.OUTOPTIONS, "noiso")
        converter.AddOption("w", converter.OUTOPTIONS)
        obmol = openbabel.OBMol()
        if not converter.ReadString(obmol, str(inchi)):
            return None
        inchikey = converter.WriteString(obmol, True)  # second argument is trimWhitespace
        if inchikey == '':
            return None
        else:
            return inchikey

    @staticmethod
    def mol2inchi(s):
        openbabel.obErrorLog.SetOutputLevel(-1)

        converter = openbabel.OBConversion()
        converter.SetInAndOutFormats('mol', 'inchi')
        converter.AddOption("x", converter.OUTOPTIONS, "noiso")
        converter.AddOption("w", converter.OUTOPTIONS)
        obmol = openbabel.OBMol()
        if not converter.ReadString(obmol, str(s)):
            return None
        inchi = converter.WriteString(obmol, True) # second argument is trimWhitespace
        if inchi == '':
            return None
        else:
            return inchi

    @staticmethod
    def inchi2smiles(inchi):
        openbabel.obErrorLog.SetOutputLevel(-1)
        
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats('inchi', 'smiles')
        obmol = openbabel.OBMol()
        conv.ReadString(obmol, str(inchi))
        smiles = conv.WriteString(obmol, True) # second argument is trimWhitespace
        if smiles == '':
            return None
        else:
            return smiles
            
    @staticmethod
    def smiles2smiles(smiles_in):
        openbabel.obErrorLog.SetOutputLevel(-1)
        
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats('smiles', 'smiles')
        obmol = openbabel.OBMol()
        conv.ReadString(obmol, str(smiles_in))
        smiles_out = conv.WriteString(obmol, True) # second argument is trimWhitespace
        if smiles_out == '':
            return None
        else:
            return smiles_out

    @staticmethod
    def smiles2inchi(smiles):
        openbabel.obErrorLog.SetOutputLevel(-1)
        
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats('smiles', 'inchi')
        conv.AddOption("F", conv.OUTOPTIONS)
        conv.AddOption("T", conv.OUTOPTIONS)
        conv.AddOption("x", conv.OUTOPTIONS, "noiso")
        conv.AddOption("w", conv.OUTOPTIONS)
        obmol = openbabel.OBMol()
        conv.ReadString(obmol, str(smiles))
        inchi = conv.WriteString(obmol, True) # second argument is trimWhitespace
        if inchi == '':
            return None
        else:
            return inchi

    def __str__(self):
        return "%s\nInChI: %s\npKas: %s\nmajor MS: nH = %d, charge = %d" % \
            (self.compound_id, self.inchi, ', '.join(['%.2f' % p for p in self.p_kas]),
             self.number_of_protons[self.major_ms_ph_7], self.charges[self.major_ms_ph_7])

    @property
    def molecule(self):
        try:
            return Molecule.from_smiles(self.smiles_ph_7)
        except OpenBabelError:
            try:
                return Molecule.from_inchi(self.inchi)
            except OpenBabelError:
                return None

    @property
    def non_polar_surface_area(self):
        if self._npsa is None:
            self._npsa = non_polar_surface_area(self.inchi)
        return self._npsa

    @property
    def concentration(self):
        if self.inchi.startswith("InChI=1S/H2O/h1H2"):  # Water
            return 55
        molecule = self.molecule
        if molecule is None:
            return 1/1000
        charged_atoms = len([c for c in molecule.atom_charges if c != 0])
        return np.math.exp(charged_atoms * 1.0425 - self.non_polar_surface_area * 0.0272) / 1000

    def _dG0_prime_vector(self, ph, ionic_strength, temperature):
        """
        Calculates the difference in kJ/mol between dG'0 and
        the dG0 of the MS with the least hydrogens (dG0[0])

        Parameters
        ----------
        ph : float
            The pH of the environment.
        ionic_strength : float
            The ionic strength of the environment.
        temperature : float
            The temperature of the environment.

        Returns
        -------
        dG0_prime_vector : ndarray
            dG'0 - dG0[0]
        """
        if self.inchi is None:
            return 0
        elif len(self.p_kas) == 0:
            dG0s = np.zeros((1, 1))
        else:
            dG0s = -np.cumsum([0] + self.p_kas) * R * temperature * np.log(10)
            dG0s = dG0s
        dh = debye_huckel(ionic_strength, temperature)
        
        # dG0' = dG0 + nH * (R T ln(10) pH + DH) - charge^2 * DH
        pseudoisomers = np.vstack([dG0s, np.array(self.number_of_protons), np.array(self.charges)]).T
        dG0_prime_vector = pseudoisomers[:, 0] + \
                           pseudoisomers[:, 1] * (R * temperature * np.log(10) * ph + dh) - \
                           pseudoisomers[:, 2]**2 * dh
        return dG0_prime_vector
        
    def _transform(self, ph, ionic_strength, temperature):
        dg0_prime_vector = self._dG0_prime_vector(ph, ionic_strength, temperature)
        return -R * temperature * logsumexp(dg0_prime_vector / (-R * temperature))

    def _ddG(self, i_from, i_to, temperature):
        """
        Calculates the difference in kJ/mol between two MSs.
            
        Returns:
            dG0[i_to] - dG0[i_from]
        """
        if not (0 <= i_from <= len(self.p_kas)):
            raise ValueError('MS index is out of bounds: 0 <= %d <= %d' % (i_from, len(self.p_kas)))

        if not (0 <= i_to <= len(self.p_kas)):
            raise ValueError('MS index is out of bounds: 0 <= %d <= %d' % (i_to, len(self.p_kas)))

        if i_from == i_to:
            return 0
        elif i_from < i_to:
            return sum(self.p_kas[i_from:i_to]) * R * temperature * np.log(10)
        else:
            return -sum(self.p_kas[i_to:i_from]) * R * temperature * np.log(10)

    def get_transform(self, ph, ionic_strength, temperature):
        """
        Returns the difference in kJ/mol between dG'0 and the dG0 of the MS with index 'i'.

        Returns
        -------
        (dG'0 - sp. dG0[0])
        """
        return self._transform(ph, ionic_strength, temperature)

    def transform(self, index, ph, ionic_strength, temperature):
        """
        Returns the difference in kJ/mol between dG'0 and the dG0 of the MS with index 'i'.

        Returns
        -------
        (dG'0 - sp. dG0[0]) + (sp. dG0[0] - sp. dG0[i]) = dG'0 - sp. dG0[i]
        """
        return self._transform(ph, ionic_strength, temperature) + self._ddG(0, index, temperature)

    def transform_ph7(self, ph, ionic_strength, temperature):
        """
        Returns the transform for the major MS in pH 7
        """
        try:
            return self.transform(self.major_ms_ph_7, ph, ionic_strength, temperature)
        except Exception as e:
            logger.error("%s: %s" % (self.compound_id, str(e)))
            raise e

    def transform_neutral(self, ph, ionic_strength, temperature):
        """
        Returns the transform for the MS with no charge
        """
        try:
            return self.transform(self.charges.index(0), ph, ionic_strength, temperature)
        except ValueError:
            raise ValueError("The compound (%s) does not have a microspecies with 0 charge" % self.compound_id)

    def get_species(self, major_ms_dG0_f, temperature):
        """
        Given the chemical formation energy of the major microspecies,
        uses the pKa values to calculate the chemical formation energies
        of all other species, and returns a list of dictionaries with
        all the relevant data: dG0_f, nH, nMg, z (charge)
        """
        for i, (nH, z) in enumerate(zip(self.number_of_protons, self.charges)):
            dG0_f = major_ms_dG0_f + self._ddG(i, self.major_ms_ph_7, temperature)
            d = {'phase': 'aqueous', 'dG0_f': np.round(dG0_f, 2), 'nH': nH, 'z': z, 'nMg': 0}
            yield d
