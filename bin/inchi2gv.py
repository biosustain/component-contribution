import sys
from component_contribution.core.group_vector import GroupDecomposer, GroupDecompositionError
from component_contribution.core.group_vector import init_groups_data, inchi_to_group_vector

from optparse import OptionParser


def make_opts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-s", "--silent",
                          dest="silent",
                          default=False,
                          action="store_true",
                          help="suppress all error and warning messages")
    opt_parser.add_option("-i", "--inchi",
                          dest="inchi",
                          default=None,
                          help="input InChI string")
    opt_parser.add_option("-l", "--list_groups",
                          dest="list_groups",
                          default=False,
                          action="store_true",
                          help="list all group names")
    return opt_parser


if __name__ == "__main__":

    parser = make_opts()
    options, _ = parser.parse_args(sys.argv)
    if options.inchi is None and not options.list_groups:
        sys.stderr.write(parser.get_usage())
        sys.exit(-1)

    groups_data = init_groups_data()
    if options.inchi:
        group_decomposer = GroupDecomposer(groups_data)
        try:
            groupvec = inchi_to_group_vector(group_decomposer, options.inchi)
        except GroupDecompositionError as e:
            if not options.silent:
                sys.stderr.write(str(e) + '\n')
                sys.stderr.write(e.debug_table())
            sys.exit(-1)
        if options.list_groups:
            sys.stdout.write(', '.join(groups_data.group_names) + '\n')
        sys.stdout.write(', '.join("%g" % i for i in groupvec.flatten()) + '\n')
    else:
        sys.stdout.write('\n'.join(groups_data.group_names))
