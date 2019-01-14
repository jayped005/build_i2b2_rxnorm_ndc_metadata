import io, os
from collections import defaultdict

'''
Module: modifier_metadata_builder.py

Author: Jay Pedersen, UNMC, Pathology/Microbiology Department, Feb 9, 2017

Purpose:
    Build metadata rows for modifiers specified in an input text file.
'''

class modifier_metadata_builder():
# --------------------------------------------------------------------------------------|
#          Creates metadata rows for information in MED_MODIFIERS_ING.txt                   |
#          Assumes the rows can be added as-is, no modification of paths or level.      |
# --------------------------------------------------------------------------------------|

    '''
    MOTE: MED_MODIFIERS_ING.txt content example
    "C_HLEVEL"|"C_FULLNAME"|"C_NAME"|"C_SYNONYM_CD"|"C_VISUALATTRIBUTES"|"C_TOTALNUM"|"C_BASECODE"|"C_METADATAXML"|"C_FACTTABLECOLUMN"|"C_TABLENAME"|"C_COLUMNNAME"|"C_COLUMNDATATYPE"|"C_OPERATOR"|"C_DIMCODE"|"C_COMMENT"|"C_TOOLTIP"|"M_APPLIED_PATH"|"UPDATE_DATE"|"DOWNLOAD_DATE"|"IMPORT_DATE"|"SOURCESYSTEM_CD"|"VALUETYPE_CD"|"M_EXCLUSION_CD"|"C_PATH"|"C_SYMBOL"|"PCORI_BASECODE"|"PCORI_CUI"|"PCORI_NDC"
    1|"\PCORI_MOD\RX_ING_BASIS\"|"* RX Basis"|"N"|"DAE"|5|||"modifier_cd"|"MODIFIER_DIMENSION"|"modifier_path"|"T"|"like"|"\PCORI_MOD\RX_ING_BASIS\"||"Basis of Medication Event (Dispense vs. Prescribe). ** This modifier is used by SCILHS and must be present or the ETL script will not choose the appropriate PopMedNet table."|"\PCORI\MEDICATION\INGREDIENT\%"|2015/07/06 12:17:00 PM|2015/07/06 12:18:00 PM|2015/07/06 12:18:00 PM|"PCORNET_CDM"|||"\PCORI_MOD\"|"RX_BASIS"|""||
    '''

    def __init__(self):
        self.metadata_rows = {}
    # end constructor

    def build(self, path_prefix, path_prefix_level, modifiers_fn):

        def chomp(s): return s.rstrip('\n').rstrip('\r')

        # Read MED_MODIFIERS_ING.txt, same folder as the script
        fn = os.path.dirname(os.path.abspath(__file__))+os.path.sep+modifiers_fn
        f = io.open(fn, 'r', encoding='utf-8')
        header = chomp(f.readline())
        field_name_list = [x.strip('"').strip() for x in header.split('|')]
        fields_d = { nm: idx for idx,nm in enumerate(field_name_list) }
        name_count = len(field_name_list)
        delim = '|'
        line_num = 1 # read header already
        while True:
            # read and concatenate lines until the number of fields >= expected (support multi-line cells)
            line = ''; field_values = []
            while len(field_values) < name_count:
                rawline = f.readline()
                if not rawline: break  # EOF
                line += chomp(rawline)
                field_values = [x.strip('"') for x in line.split(delim)]
                line_num += 1
            # end while len
            if len(line)==0: break # EOF with no additional data is okay
            if len(field_values) != name_count:
                raise ValueError('Cant parse line %d in [%s], field count wrong (%d vs %d)'
                                 % (line_num, fn, len(field_values), name_count))
            # create metadata row
            fields = field_values
            c_fullname = fields[fields_d['C_FULLNAME']]
            self.metadata_rows[c_fullname] = \
                { 'c_fullname': c_fullname,
                  'c_hlevel':   fields[fields_d['C_HLEVEL']],
                  'c_name':     fields[fields_d['C_NAME']],
                  'c_basecode': fields[fields_d['C_BASECODE']],
                  'm_applied_path': fields[fields_d['M_APPLIED_PATH']],
                  'c_tooltip':  fields[fields_d['C_TOOLTIP']],
                  'children': 0,
                  }
            for tup in [('C_VISUALATTRIBUTES', 'c_visualattributes'),
                        ('C_OPERATOR', 'c_operator'),
                        ('C_METADATAXML', 'c_metadataxml'),
                        ('C_TABLENAME', 'c_tablename'),
                        ('C_FACTTABLECOLUMN', 'c_facttablecolumn'),
                        ('C_COLUMNNAME', 'c_columnname'),
                        ('C_COLUMNDATATYPE', 'c_columndatatype'),
                        ('C_DIMCODE', 'c_dimcode')]:
                field_name, dictionary_name = tup
                if field_name in field_name_list:
                    self.metadata_rows[c_fullname][dictionary_name] = fields[fields_d[field_name]]
            # end
        # end reading MED_MODIFIERS
        f.close()
        return self.metadata_rows # Done
    # end build

# end class modifier_metadata_builder
