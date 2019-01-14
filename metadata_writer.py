from __future__ import print_function
import time, io, sys

'''
Module: metadata_writer.py

Author:
    Jay Pedersen, UNMC, Pathology/Microbiology Department, September 2018

Purpose:
    Define metadata_writer class, which is given python metadata structures, and writes
    equivalent i2b2-format metadata to a specified text file.
'''

class metadata_writer():
# Metadata file writer class
    def __init__(self, fn, append=False, include_dates=False, include_tty=False, \
                 csv_fields_file=None, csv_fields_file_delim=None, strftime_format=None):
        self.fout = open(fn,'a' if append else 'w') # create or append
        if csv_fields_file:
            with io.open(csv_fields_file, 'r', encoding='utf-8') as f:
                header = f.readline().rstrip('\n').rstrip('\r')
                self.fieldnames = [x.strip('"').strip().upper() for x in header.split(csv_fields_file_delim)]
                # NOTE: "c_name" ==> C_NAME, etc
        else:
            self.fieldnames = \
                ('C_FULLNAME|C_HLEVEL|C_NAME|C_BASECODE|C_VISUALATTRIBUTES|M_APPLIED_PATH|C_SYNONYM_CD|' \
                 + 'C_TABLENAME|C_COLUMNNAME|C_COLUMNDATATYPE|C_OPERATOR|C_DIMCODE|C_FACTTABLECOLUMN|' \
                 + 'C_TOOLTIP|C_TOTALNUM|FACT_COUNT|SOURCESYSTEM_CD|C_METADATAXML').split('|')
            if include_dates: self.fieldnames.extend(['UPDATE_DATE','DOWNLOAD_DATE','IMPORT_DATE'])
            if include_tty: self.fieldnames.append('TTY')
        self.fieldname_set = set(self.fieldnames)
        self.date_str = time.strftime('%Y/%m/%d %I:%M:%S %p' if not strftime_format else strftime_format)
        if not append: # created (not appending) ==> write header line
            print('|'.join(self.fieldnames),file=self.fout)

    def write_metadata_rows(self, metadata_paths):
        def clean_str(s):
            ''' clean_str: for SQL*LOADER files.  enclose whole string with double quotes (embedded double-quotes replaced by \" ) '''
            return make_utf8('"' + s.strip().replace('"', r'\"') + '"')  # CSV-safe string, safe for Oracle import

        def make_utf8(s):
            ''' make_unicode: for python 2, return unicode(s) otherwise return s '''
            return unicode(s) if sys.version_info[0] < 3 else s

        fieldnames, fieldname_set, date_str = self.fieldnames, self.fieldname_set, self.date_str
        field = {fieldname: '' for fieldname in fieldnames}  # dictionary of field names, '' for unrecognized
        for path in sorted(metadata_paths.keys()):
            field['C_FULLNAME'] = clean_str(path)
            field['C_HLEVEL'] = str(metadata_paths[path]['c_hlevel'])  # prefix is constant, don't concern with it
            field['C_NAME'] = clean_str(metadata_paths[path]['c_name'])
            if 'C_NAME_ORIG' in fieldname_set: field['C_NAME_ORIG'] = clean_str(metadata_paths[path]['c_name'])

            # start basecode processing
            c_basecode = metadata_paths[path]['c_basecode']
            field['C_BASECODE'] = clean_str(c_basecode) if c_basecode else '' # leave empty string alone, dont quote it
            colon_pos = c_basecode.find(':')
            if 'PCORI_BASECODE' in fieldname_set: # SHRINE/HARVARD
                field['PCORI_BASECODE'] = '' if colon_pos < 0 else clean_str(c_basecode[colon_pos + 1:])
            if 'PCORI_NDC' in fieldname_set:
                field['PCORI_NDC'] = '' if (not c_basecode.startswith('NDC:')) else clean_str(c_basecode[colon_pos + 1:])
            if 'PCORI_CUI' in fieldname_set: # SHRINE/HARVARD
                if c_basecode.startswith('RXNORM:'):
                    field['PCORI_CUI'] = clean_str(c_basecode[colon_pos + 1:])
                elif c_basecode.startswith('NDC:'):
                    # parse C_FULLNAME and look for RXNORM code in previous position
                    c_fullname_elements = metadata_paths[path]['c_fullname'].split('\\')
                    field['PCORI_CUI'] = '' if len(c_fullname_elements) < 2 else clean_str(c_fullname_elements[-2])
            # end basecode processing

            field['C_DIMCODE'] = clean_str(metadata_paths[path].get('c_dimcode', path)) # JGP 2017-07-02, works for ITCP and path-based
            field['C_TOOLTIP'] = clean_str(metadata_paths[path].get('c_tooltip', path[1:-1]))
            field['C_SYNONYM_CD'] = clean_str('N')
            if 'C_TOTALNUM' in fieldname_set: field['C_TOTALNUM'] = '0'
            if 'FACT_COUNT' in fieldname_set: field['FACT_COUNT'] = '0'
            field['C_VISUALATTRIBUTES'] = clean_str(metadata_paths[path].get('c_visualattributes',
                                                                             ('FA ' if metadata_paths[path]['children'] > 0 else 'LA '))) # folder?
            field['C_FACTTABLECOLUMN'] = clean_str(metadata_paths[path].get('c_facttablecolumn', 'concept_cd'))
            field['C_TABLENAME'] = clean_str(metadata_paths[path].get('c_tablename', 'concept_dimension'))  # NOT the ITCP table
            field['C_OPERATOR'] = clean_str(metadata_paths[path].get('c_operator', 'LIKE'))  # path-based transitive-closure computation
            field['C_COLUMNNAME'] = clean_str(metadata_paths[path].get('c_columnname', 'concept_path'))
            field['C_COLUMNDATATYPE'] = clean_str(metadata_paths[path].get('c_columnndatatype', 'T'))  # text
            field['M_APPLIED_PATH'] = clean_str(metadata_paths[path].get('m_applied_path', '@'))
            field['SOURCESYSTEM_CD'] = clean_str(metadata_paths[path].get('sourcesystem_cd', 'rxnav.nlm.nih.gov'))
            for datefield in ['UPDATE_DATE', 'DOWNLOAD_DATE', 'IMPORT_DATE']:
                if datefield in fieldname_set: field[datefield] = clean_str(date_str)
            c_metadataxml = metadata_paths[path].get('c_metadataxml', '')
            field['C_METADATAXML'] = clean_str(c_metadataxml) if c_metadataxml else ''
            if 'TTY' in fieldname_set: field['TTY'] = clean_str(metadata_paths[path].get('tty',''))  # TTY
            # Deal with fields, requested by SHRINE metadata (etc), not 'standard' but need to be in CSV.
            # ==> For now, initialized as an empty string in field and will show up as such in the result. e.g. |||
            '''
            if 'C_PATH' in fieldname_set:
                field['C_PATH'] = clean_str(opts.pathprefix + path_sep \
                                                         + path_sep.join(reversed(path[1:])) + path_sep)
            if 'C_SYMBOL' in fieldname_set:
                field['C_SYMBOL'] = clean_str(path[0])
            '''
            print(make_utf8('|'.join(field[x] for x in fieldnames)), file=self.fout)
        # end for path
        return

    def close(self):
        self.fout.close()

# end class metadata_writer