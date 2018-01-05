#!/usr/bin/python

import logging
import re

import gzip

from component_contribution.exceptions import KEGGParsingException

KEGG_REACTION_REGEX = re.compile("R\d{5}")
KEGG_ORTHOLOGY_REGEX = re.compile('^(K\d{5})  (.*)$')


def normalize_names(name_str):
    """
    Normalize a KEGG-style list of names.
    """
    all_names = name_str.replace('\t', ' ').split(';')
    return [n.strip() for n in all_names]


def normalize_reactions(reactions_str):
    """

    Normalize a KEGG-style list of reaction IDs.
    
    NOTE(flamholz): Some enzymes have lists of reactions as such:
        "RXXXXX > RXXXXY RYYYZZ"
    where RXXXXX is a general reaction and the others have specific
    substrates. We may want special parsing for this, but for now 
    we don't have it. 
    
    Parameters
    ----------
    reactions_str : str
        The string containing a list of reactions.
        
    Returns
    -------
    reaction_ids : list
        A list of KEGG reaction IDs
    """
    if not reactions_str:
        return []
    
    reaction_ids = KEGG_REACTION_REGEX.findall(reactions_str)
    return reaction_ids


def normalize_organisms(organisms_str):
    """
    Normalize a KEGG-style list of organism names.

    Parameters
    ----------
    organisms_str : str
        A string with organisms

    """
    return organisms_str.split('\t')


def parse_orthology_mapping(orthology_str):
    """Parses the orthology string to a mapping.
    
    Parameters
    ----------
    orthology_str : str
        The orthology string in the KEGG file.
    
    Returns
    -------
    orthology_mapping : dict
        A mapping from orthology IDs to names.
    """
    entries = orthology_str.split('\t')

    orthology_mapping = {}

    for match in map(KEGG_ORTHOLOGY_REGEX.match, entries):
        if not match:
            continue
        
        groups = match.groups()
        orthology_mapping[groups[0]] = groups[1]
    return orthology_mapping


def parse_organism_to_gene_mapping(genes_str):
    """Parses the genes string to a mapping.
    
    TODO(flamholz): Keep open reading frame data as well.
    
    Parameters
    ----------
    genes_str : str
        The orthology string in the KEGG file.
    
    Returns
    -------
    organism_gene_mapping : dict
        A mapping from organisms to gene names.
    """
    splitted = genes_str.split('\t')
    pattern = re.compile('^([A-Z]{3}): (.*)$')

    organism_gene_mapping = {}
    for match in map(pattern.match, splitted):
        if not match:
            continue
        
        groups = match.groups()

        gene_ids_w_names = re.split('\s', groups[1])
        gene_ids = [s.split('(')[0] for s in gene_ids_w_names]
        if gene_ids:
            organism_gene_mapping[groups[0]] = gene_ids
    
    return organism_gene_mapping


class EntryDictWrapper(dict):
    
    def get_string_field(self, field_name, default_value=None):
        if field_name not in self:
            if default_value is not None:
                return default_value
            raise Exception("Missing obligatory string field: " + field_name)
            
        return self[field_name]
    
    def get_string_list_field(self, field_name, default_value=None):
        val = self.get_string_field(field_name, default_value=False)
        
        if val is False:
            if default_value is None:
                raise Exception("Missing obligatory string-list field: " + field_name)
            return default_value
        return val.split()
        
    def get_bool_field(self, field_name, default_value=True):
        val = self.get_string_field(field_name, default_value=False)
        
        if val is False:
            if default_value is None:
                raise Exception("Missing obligatory boolean field: " + field_name)
            return default_value
        elif val.upper() == 'TRUE':
            return True
        elif val.upper() == 'FALSE':
            return False
    
    def get_float_field(self, field_name, default_value=None):
        val = self.get_string_field(field_name, default_value=False)
        
        if val is False:
            if default_value is None:
                raise Exception("Missing obligatory float field: " + field_name)
            return default_value
        return float(val)
    
    def get_v_float_field(self, field_name, default_value=()):
        val = self.get_string_field(field_name, default_value=False)
        
        if val is False:
            if default_value is None:
                raise Exception("Missing obligatory vector-float field: " + field_name)
            return default_value
        return [float(x) for x in val.split()]


class KEGGData(dict):
    """A class encapsulating a parsed KEGG file."""

    def __init__(self):
        """Initialize the ParsedKeggFile object."""
        super(KEGGData, self).__init__(iterable=[])
        self.ordered_entries = []
        
    def _add_entry(self, entry, fields):
        """Protected helper for adding an entry from the file.
        
        Args:
            entry: the entry key.
            fields: the fields for the entry.
        """
        if entry in self:
            logging.warning('Overwriting existing entry for %s', entry)
        else:
            self.ordered_entries.append(entry)
        self[entry] = EntryDictWrapper(fields)

    def entries(self):
        return self.ordered_entries

    @classmethod
    def from_kegg_file(cls, file_name):
        """Parses a file from KEGG.
    
        Parameters
        ----------
        file_name : str
            The file handle or name of the file to parse.
        
        Returns
        -------
        data : KEGGData
            A dictionary mapping entry names to fields.
        """
        if isinstance(file_name, str):
            if file_name[-3:] == '.gz':
                kegg_file = gzip.open(file_name)
            else:
                kegg_file = open(file_name, 'r')
        else:
            kegg_file = file_name
        return cls._from_kegg_file_handle(kegg_file)

    @classmethod
    def _from_kegg_file_handle(cls, kegg_file):
        """Parses a file from KEGG. Uses a file handle directly.
        
        For testing.
    
        Parameters
        ----------
        kegg_file : file
            The file handle.
        
        Returns
        -------
        data : KEGGData
            A dictionary mapping entry names to fields.
        """
        parsed_file = cls()
    
        line_counter = 0
        line = kegg_file.readline()
        field_map = {}
        field = None
    
        while line:
            if line[0:3] == '///':
                if field_map:
                    entry = re.split('\s\s+', field_map['ENTRY'])[0].strip()
                    parsed_file._add_entry(entry, field_map)
                field = None
                field_map = {}
            elif line[0] in [' ', '\t']:
                if field is None:
                    raise KEGGParsingException('First line starts with a whitespace (space/tab)')
                value = line.strip()
                field_map[field] = field_map[field] + "\t" + value
            else:
                try:
                    field, value = line.split(None, 1)
                except ValueError:
                    raise KEGGParsingException('ERROR: line %d cannot be split: %s' % (line_counter, line))
                field_map[field] = value
    
            line = kegg_file.readline()
            line_counter += 1
        if 'ENTRY' in field_map:
            entry = re.split('\s\s+', field_map['ENTRY'])[0].strip()
            parsed_file._add_entry(entry, field_map)
        kegg_file.close()
        return parsed_file
    
    @classmethod
    def from_kegg_api(cls, response):
        """Parses a file from KEGG. The result string from the KEGG API.
        
        For testing.
    
        Parameters
        ----------
        response : str
            The string that is the result of serv.bget(...) using the KEGG API
        
        Returns
        -------
        data :  KEGGData
            A dictionary mapping entry names to fields.
        """
        parsed_file = KEGGData()
    
        curr_field = ""
        field_map = {}
    
        for line in response.split('\n'):
            field = line[0:12].strip()
            value = line[12:].strip()
    
            if field[:3] == "///":
                entry = re.split('\s\s+', field_map['ENTRY'])[0]
                parsed_file._add_entry(entry, field_map)
                field_map = {}
            else:
                if field != "":
                    curr_field = field
                if curr_field in field_map:
                    field_map[curr_field] = field_map[curr_field] + "\t" + value
                else:
                    field_map[curr_field] = value
    
        if 'ENTRY' in field_map:
            entry = re.split('\s\s+', field_map['ENTRY'])[0]
            parsed_file._add_entry(entry, field_map)
        return parsed_file

    @staticmethod
    def parse_kegg_reaction_line(line):
        rexp = '([a-zA-Z0-9,_]+)\s+([C\s\+\d\.]+)\s+(<?[-=]>?)\s+([C\s\+\d\.]+)(.*)'
        try:
            rid, left_clause, dir_clause, right_clause, remainder = next(re.finditer(rexp, line))
        except Exception as e:
            raise Exception(str(e) + ': ' + line)
        
        if dir_clause in ['=>', '->', '<=>', '<->', '=', '-']:
            reaction = left_clause + " <=> " + right_clause
        elif dir_clause in ['<=', '<-']:
            reaction = right_clause + " <=> " + left_clause
        else:
            raise ValueError("unclear reaction direction symbol: " + dir_clause)
    
        flux = 1
        if remainder != "":
            for (f) in re.findall('\(x([0-9\.\-\s]+)\)', remainder):
                flux = float(f)
        
        return reaction, rid, flux
    
    @staticmethod
    def parse_reaction_module(field_map):
        rids = []
        fluxes = []
        reactions = []
        for line in field_map["REACTION"].split('\t'):
            if line.strip() == '':
                continue
            reaction, rid, flux = KEGGData.parse_kegg_reaction_line(line)
            reactions.append(reaction)
            rids.append(rid)
            fluxes.append(flux)
        
        return rids, fluxes, reactions

    @staticmethod
    def parse_bound_module(field_map):
        bounds = {}  # a dictionary from KEGG IDs to a tuple of (low,up) bounds
        rexp = '(C[0-9]+)\s+([0-9e\-\+]+)\s*(.*)'
        for line in field_map["BOUND"].split('\t'):
            try:
                cid, low, up = re.findall(rexp, line)[0]
            except Exception as e:
                raise Exception(str(e) + ': ' + line)
            up = up or low
            low = float(low.strip())
            up = float(up.strip())
            bounds[cid] = (low, up)
        return bounds
