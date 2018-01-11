import re
import openbabel
from component_contribution import chemaxon
from component_contribution.exceptions import OpenBabelError


class Molecule(object):
    # for more rendering options visit:
    # http://www.ggasoftware.com/opensource/indigo/api/options#rendering
    _ob_elements = openbabel.OBElementTable()
    _ob_smarts = openbabel.OBSmartsPattern()

    @staticmethod
    def _convert_molecule(obmol, fmt='inchi'):
        converter = openbabel.OBConversion()
        converter.AddOption("w", converter.OUTOPTIONS)
        converter.SetOutFormat(fmt)
        res = converter.WriteString(obmol)
        if not res:
            raise OpenBabelError("Cannot convert OBMol to %s" % fmt)

        elif fmt == 'smiles' or fmt == 'smi':
            res = res.split()
            if len(res) == 0:
                raise OpenBabelError("Cannot convert OBMol to %s" % fmt)
            else:
                return res[0]

        elif fmt == 'inchi':
            return res.strip()

        else:
            return res

    def __init__(self, title: str, obmol: openbabel.OBMol, smiles: str, inchi: str):
        self._title = title
        self._obmol = obmol
        self.smiles = smiles
        self.inchi = inchi

    @classmethod
    def get_number_of_elements(cls):
        return cls._ob_elements.GetNumberOfElements()
    
    @classmethod
    def get_all_elements(cls):
        return [cls._ob_elements.GetSymbol(i) for i in range(cls.get_number_of_elements())]

    @classmethod
    def get_symbol(cls, atomic_num):
        return cls._ob_elements.GetSymbol(atomic_num)
            
    @classmethod
    def get_atomic_mum(cls, elem):
        return cls._ob_elements.GetAtomicNum(elem)
    
    @classmethod
    def verify_smarts(cls, smarts):
        return cls._ob_smarts.Init(smarts)

    @classmethod
    def from_smiles(cls, smiles):
        if smiles is None:
            raise OpenBabelError("Cannot read the SMILES: None")
        obmol = openbabel.OBMol()
        converter = openbabel.OBConversion()
        converter.AddOption("w", converter.OUTOPTIONS)
        converter.SetInFormat("smiles")
        if not converter.ReadString(obmol, smiles):
            raise OpenBabelError("Cannot read the SMILES string: " + smiles)
        else:
            inchi = cls._convert_molecule(obmol, 'inchi')
            return cls(smiles, obmol, smiles, inchi)

    @classmethod
    def from_inchi(cls, inchi):
        obmol = openbabel.OBMol()
        converter = openbabel.OBConversion()
        converter.AddOption("w", converter.OUTOPTIONS)
        converter.SetInFormat("inchi")
        if not converter.ReadString(obmol, inchi):
            raise OpenBabelError("Failed to create Molecule from InChI: " + inchi)
        else:
            smiles = cls._convert_molecule(obmol, 'smi')
            return cls(inchi, obmol, smiles, inchi)

    @classmethod
    def from_mol(cls, mol):
        obmol = openbabel.OBMol
        converter = openbabel.OBConversion()
        converter.AddOption("w", converter.OUTOPTIONS)
        converter.SetInFormat("mol")
        if not converter.ReadString(obmol, mol):
            raise OpenBabelError("Failed to create Molecule from MOL file:\n" + mol)
        else:
            inchi = cls._convert_molecule(obmol, 'inchi')
            smiles = cls._convert_molecule(obmol, 'smi')
            return Molecule("", obmol, smiles, inchi)

    @classmethod
    def from_ob_mol(cls, obmol):
        try:
            inchi = cls._convert_molecule(obmol, 'inchi')
            smiles = cls._convert_molecule(obmol, 'smi')
            return Molecule("", obmol, smiles, inchi)
        except OpenBabelError:
            raise OpenBabelError("Failed to create Molecule from OBMol")

    @classmethod
    def _from_format(cls, string, fmt='inchi'):
        if fmt == 'smiles' or fmt == 'smi':
            return cls.from_smiles(string)
        if fmt == 'inchi':
            return cls.from_inchi(string)
        if fmt == 'mol':
            return cls.from_mol(string)
        if fmt == 'obmol':
            return cls.from_ob_mol(string)

    @staticmethod
    def smiles2inchi(smiles):
        converter = openbabel.OBConversion()
        converter.AddOption("w", converter.OUTOPTIONS)
        converter.SetInAndOutFormats("smiles", "inchi")
        obmol = openbabel.OBMol()
        if not converter.ReadString(obmol, smiles):
            raise OpenBabelError("Cannot read the SMILES string: " + smiles)
        return converter.WriteString(obmol).strip()

    @staticmethod
    def inchi2smiles(inchi):
        converter = openbabel.OBConversion()
        converter.AddOption("w", converter.OUTOPTIONS)
        converter.SetInAndOutFormats("inchi", "smiles")
        obmol = openbabel.OBMol()
        if not converter.ReadString(obmol, inchi):
            raise OpenBabelError("Cannot read the InChI string: " + inchi)
        return converter.WriteString(obmol).split()[0]

    def __str__(self):
        return self.title or self.smiles or self.inchi or ""
        
    def __len__(self):
        return self.number_of_atoms
    
    def clone(self):
        return self.__class__(self.title, openbabel.OBMol(self.obmol), self.smiles, self.inchi)

    @property
    def title(self):
        return self._title

    @title.setter
    def title(self, title):
        self.title = title 

    def remove_hydrogens(self):
        self.obmol.DeleteHydrogens()
    
    def remove_atoms(self, indices):
        self.obmol.BeginModify()
        for i in sorted(indices, reverse=True):
            self.obmol.DeleteAtom(self.obmol.GetAtom(i+1))
        self.obmol.EndModify()
        self.smiles = self.to_string('smi')
        self.inchi = self.to_string('inchi')
        
    def set_atomic_num(self, index, new_atomic_num):
        self.obmol.GetAtom(index+1).SetAtomicNum(new_atomic_num)
        self.smiles = self.to_string('smi')
        self.inchi = self.to_string('inchi')

    @property
    def obmol(self):
        return self._obmol
    
    def to_string(self, fmt='inchi'):
        return self._convert_molecule(self.obmol, fmt=fmt)
    
    def to_mol(self):
        return self.to_string('mol')

    def to_inchi(self):
        return self.to_string(fmt='inchi')

    @property
    def formula(self):
        tokens = re.findall('InChI=1S?/([0-9A-Za-z\.]+)', self.inchi)
        if len(tokens) == 1:
            return tokens[0]
        elif len(tokens) > 1:
            raise ValueError('Bad InChI: ' + self.inchi)
        else:
            return ''

    @property
    def exact_mass(self):
        return self.obmol.GetExactMass()

    @property
    def atom_bag_and_charge(self):
        atom_bag, major_ms_charge = chemaxon.atom_bag_and_charge(self.inchi)
        return atom_bag, major_ms_charge

    @property
    def hydrogens_and_charge(self):
        atom_bag, charge = self.atom_bag_and_charge
        return atom_bag.get('H', 0), charge

    @property
    def number_of_electrons(self):
        """
        Calculates the number of electrons in a given molecule.
        """
        atom_bag, fixed_charge = self.atom_bag_and_charge
        return atom_bag.get('e-', 0)

    @property
    def number_of_atoms(self):
        return self.obmol.NumAtoms()

    @property
    def atoms(self):
        return [self.obmol.GetAtom(i+1) for i in range(self.obmol.NumAtoms())]
    
    def find_smarts(self, smarts):
        """
        Corrects the pyBel version of Smarts.findall() which returns results as tuples,
        with 1-based indices even though Molecule.atoms is 0-based.

        Parameters
        ----------
        smarts: str
            The SMARTS query to search for.
        
        Returns
        -------
        matches : list
            The re-mapped list of SMARTS matches.
        """
        def shift_left(m):
            return [(n - 1) for n in m]

        self._ob_smarts.Init(smarts)
        if self._ob_smarts.Match(self.obmol):
            match_list = self._ob_smarts.GetMapList()

            return map(shift_left, match_list)
        else:
            return []

    @property
    def atom_charges(self):
        """
        Returns
        -------
        atom_charges : list
            A list of charges, according to the number of atoms in the molecule
        """
        return [atom.GetFormalCharge() for atom in self.atoms]


if __name__ == '__main__':
    mol = Molecule.from_inchi('InChI=1S/H2/h1H')
    print(mol.exact_mass)
