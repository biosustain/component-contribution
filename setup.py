from setuptools import setup, find_packages

requirements = ['scipy>=0.14.0',
                'numpy>=1.6.2',
                'pandas>=0.21',
                'bioservices>=1.5',
                'requests>=2.18']

setup(
    name='component_contribution',
    version='2.0',
    author='Elad Noor, Joao Cardoso',
    author_email='noor@imsb.biol.ethz.ch, joaca@biosustain.dtu.dk',
    description='Standard reaction Gibbs energy estimation for biochemical reactions',
    license='MIT',
    packages=find_packages(),
    url='https://github.com/eladnoor/component-contribution',
    install_requires=requirements,
    data_files=[('data', ['data/TECRDB.tsv', 
                          'data/redox.tsv',
                          'data/formation_energies_transformed.tsv',
                          'data/equilibrator_compounds.json.gz',
                          'data/kegg_additions.tsv',
                          'data/kegg_compounds.json.gz',
                          'data/cached_compounds.json.gz']),
                ('cache', ['cache/component_contribution_python.mat'])],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5'
    ],
)

