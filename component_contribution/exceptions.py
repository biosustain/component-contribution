class ParseException(Exception):
    pass


class ChemAxonError(Exception):
    pass


class KEGGNonCompoundException(Exception):
    pass


class KEGGReactionNotBalancedException(Exception):
    def __init__(self, msg, atom_bag=None):
        Exception.__init__(self, msg)
        self.atom_bag = atom_bag or {}


class KEGGMissingModuleException(Exception):
    pass


class OpenBabelError(Exception):
    pass


class GroupsDataError(Exception):
    pass


class MalformedGroupDefinitionError(GroupsDataError):
    pass


class GroupDecompositionError(Exception):
    def __init__(self, msg, decomposition=None):
        Exception.__init__(self, msg)
        self.decomposition = decomposition

    def __str__(self):
        return Exception.__str__(self)

    @property
    def debug_table(self):
        if self.decomposition is not None:
            return self.decomposition.ToTableString()
        else:
            return ''


class KEGGParsingException(Exception):
    pass
