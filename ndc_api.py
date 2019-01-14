import zipfile

'''
Module: ndc_api.py

Author:
    Jay Pedersen, UNMC, Pathology/Microbiology Department, September 19, 2018

Purpose:
    Define ndc_api class, whose current sole purpose is to provide the get_ndc_name method
    for determining the name given by the FDA for an NDC code.  This is determined from a file
    provided by the FDA.
'''

class ndc_api:

    def __init__(self):

        def read_name_mapping_from_fda_file():
            ndc_zip_filename = 'ndc_names_from_fda.zip'
            zf = zipfile.ZipFile(ndc_zip_filename)
            ndc_csv_filename = 'ndc_names_from_fda.csv'
            ndc_to_name = {}
            ndc_code_set = set()
            fin = None
            try:
                fin = zf.open(ndc_csv_filename)
            except Exception as e:
                print('Could not open <<<%s>>> in <<<%s>>>' % (ndc_csv_filename, ndc_zip_filename))
                print(str(e))  # error message
                # terminate program by not passing the exception
            first_line = True
            while True:
                rawline = fin.readline()
                if not rawline: break  # EOF
                if first_line: first_line = False; continue  # skip header, (NDC_CODE NDC_NAME)
                line = rawline.decode('latin-1').rstrip('\n').rstrip('\r')  # decode (because of zipfile) and chomp
                # eg: NDC:00003085722 SPRYCEL (dasatinib) 1 BOTTLE in 1 CARTON (0003-0857-22)  > 30 TABLET in 1 BOTTLE
                fields = line.split('\t')  # tab-separated
                ndc_code = fields[0][4:]  # skip past 'NDC:'
                if len(ndc_code) != 11: raise ValueError('NDC code [%s] is not length 11' % ndc_code)
                ndc_name = fields[1].strip('"')  # drop optional surrounding quotes
                ndc_to_name[ndc_code] = (ndc_name, 'NLM:RXNORM:RESTAPI')
                ndc_code_set.add(ndc_code)
            fin.close()
            # done building ndc_to_name dictionary
            return ndc_to_name
        # end read_name_mapping_from_fda_file

        self.ndc_to_name = read_name_mapping_from_fda_file()

    # end constructor

    def get_ndc_name(self, ndc_code):
        return self.ndc_to_name[ndc_code][0] if ndc_code in self.ndc_to_name else ndc_code
    # end get_ndc_name

# end class ndc_api

