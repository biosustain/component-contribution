
from component_contribution.component_contribution_trainer import ComponentContribution
from component_contribution.CfB_functions import formation_energy
import os

cc = ComponentContribution.init()

listKeggMolecules = ['C01598','C00002']


#atpmol='\n \n \n 31 33  0  0  1  0  0  0  0  0999 V2000\n   29.4250  -14.6015    0.0000 N   0  0  3  0  0  0  0  0  0  0  0  0\n   30.4825  -15.3378    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   28.1393  -15.0165    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n   29.9040  -13.2050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   31.7155  -14.4321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   30.6227  -16.7578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   26.9999  -14.1633    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   27.7071  -16.2669    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n   31.3532  -13.2108    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   32.9718  -14.9814    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   31.9434  -17.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   25.8897  -14.9580    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n   26.3221  -16.2669    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n   28.5133  -17.4063    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   33.1353  -16.4655    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   34.0819  -14.1575    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   24.5691  -14.5315    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   25.8780  -17.5874    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   23.5349  -15.4664    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   22.1325  -15.4664    0.0000 P   0  0  3  0  0  0  0  0  0  0  0  0\n   20.7361  -15.4664    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   22.1268  -16.8629    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   22.1268  -14.0698    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   19.3394  -15.4605    0.0000 P   0  0  3  0  0  0  0  0  0  0  0  0\n   17.9430  -15.4605    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   19.3337  -16.8571    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   19.3337  -14.0581    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   16.5406  -15.4547    0.0000 P   0  0  3  0  0  0  0  0  0  0  0  0\n   16.5348  -16.8512    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   15.1441  -15.4430    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   16.5348  -14.0523    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0     0  0\n  3  1  1  1     0  0\n  1  4  1  0     0  0\n  2  5  2  0     0  0\n  2  6  1  0     0  0\n  3  7  1  0     0  0\n  3  8  1  0     0  0\n  4  9  2  0     0  0\n  5 10  1  0     0  0\n  6 11  2  0     0  0\n  7 12  1  0     0  0\n  8 13  1  0     0  0\n  8 14  1  6     0  0\n 10 15  2  0     0  0\n 10 16  1  0     0  0\n 12 17  1  1     0  0\n 13 18  1  6     0  0\n 17 19  1  0     0  0\n 19 20  1  0     0  0\n 20 21  1  0     0  0\n 20 22  1  0     0  0\n 20 23  2  0     0  0\n 21 24  1  0     0  0\n 24 25  1  0     0  0\n 24 26  1  0     0  0\n 24 27  2  0     0  0\n 25 28  1  0     0  0\n 28 29  1  0     0  0\n 28 30  1  0     0  0\n 28 31  2  0     0  0\n  5  9  1  0     0  0\n 11 15  1  0     0  0\n 12 13  1  0     0  0\nM  END\n\n> <ENTRY>\ncpd:C00002\n\n$$$$\n'
with open(os.path.join('../..','del_atpMolFound_string'), 'r') as myfile:
    data = myfile.read()
    data = [data]


# ATP in two different formats and this molecule: https://pubchem.ncbi.nlm.nih.gov/compound/443290#section=Top
smileMolecules = [ 'C1=NC2=C(C(=N1)N)N=CN2[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O', 'C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O','C1=CC(C(C(=C1)CCC(=O)O)O)O'] #

# testing melatonin and other
test_file = open("../../melatonin.sdf",'r')
atpmol = test_file.read()
data = [atpmol]

test_file = open("../../Structure2D_CID_115348.sdf",'r')
atpmol = test_file.read()
data.append(atpmol)


print formation_energy(listKeggMolecules,cc, 7, 0.1)

print formation_energy(data,cc, 7, 0.1,format_molecule='mol')

print formation_energy(smileMolecules,cc, 7, 0.1,format_molecule='smiles')



