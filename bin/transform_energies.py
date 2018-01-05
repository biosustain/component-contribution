import csv

from component_contribution.compound_cache import CompoundCache


def main(file_name, ph, ionic_strength, temperature):
    ccache = CompoundCache()
    with open(file_name, 'r') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        csv_reader.next()
        for row in csv_reader:
            compound_id = row[0]  # compound_id = re.findall('(C[0-9]+)_10', row[0])[0]
            dg0 = float(row[3])
            comp = ccache.get_compound(compound_id)
            dg0_prime = dg0 + comp.transform_neutral(ph, ionic_strength, temperature)
            print('%s\t%f\t%f' % (compound_id, dg0, dg0_prime))
        ccache.dump()


if __name__ == '__main__':
    pH = 7
    I = 0.2
    T = 298.15
    main('cc_compounds.csv', pH, I, T)
