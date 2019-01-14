from ndc_api import ndc_api
import time
from collections import defaultdict

'''
Module: ingredient_metadata_builder.py

Author: Jay Pedersen, UNMC, Pathology/Microbiology Department, Feb 9, 2017

Purpose:
    (a) Define ingredient_metadata_builder class to build the by-ingredient hierarchy
        for RxNorm drugs for use in SCILHS Medications metadata.
'''

class ingredient_metadata_builder():
# --------------------------------------------------------------------------------------|
#                                 Drugs by Ingredient metadata                          |
# --------------------------------------------------------------------------------------|

    def __init__(self, rxnorm_code_maps, rxnav):
        self.rxnorm_coding = rxnorm_code_maps
        self.rxnav = rxnav
        self.start_t = None # start time, set by build method
        self.ndc_api = ndc_api()
        self.ingredient_rxcui_set_raw = self.rxnorm_coding.get_ingredient_rxcui_set()
        # Treat "unmapped MIN codes" the same as IN codes, they are not tied to any ingredient
        # and will be treated as a stand-alone ingredient.
        self.ingredient_rxcui_set = self.ingredient_rxcui_set_raw | self.rxnorm_coding.get_unmapped_min_code_set()
        self.multi_ingredient_rxcui_set = self.rxnorm_coding.get_multi_ingredient_rxcui_set()
        self.min_set = self.multi_ingredient_rxcui_set
        self.drug_rxcui_set = self.rxnorm_coding.get_drug_rxcui_set()
        self.ingredient_to_drug_set_map = self.rxnorm_coding.get_ingredient_to_drug_set_map()
        self.in_to_min_map = self.rxnorm_coding.get_ingredient_to_min_set_map()
        self.metadata_rows = {}
        self.rxcui_map = {} # track name and type (eg. SCD vs GPCK,etc)
        self.ingredient_rxcui_to_paths = {}
        self.undetermined_ingredient_rxcui_map = self.rxnorm_coding.get_undetermined_ingredient_rxcui_map() # A-Z codes
        self.undetermined_ingredient_rxcuis = self.undetermined_ingredient_rxcui_map.values()
        self.path_prefix = None # set by build
    # end constructor

    def build(self, path_prefix, path_prefix_level):
        self.start_t = time.time()

        # Ingredient level
        print('Computing ingredient-level hierarchy metadata')
        root_path = '\\'+path_prefix+'\\'
        ingredient_path_prefix = root_path+'INGREDIENT'+'\\'
        self.path_prefix = ingredient_path_prefix # NOT path_prefix
        # create "by ingredients" folder, to sit under "medications" folder
        self.metadata_rows[ingredient_path_prefix] =\
            { 'c_fullname': ingredient_path_prefix,
              'c_hlevel':   path_prefix_level+1,
              'c_name':     'Drugs by Ingredient',
              'c_basecode': 'SCILHS:MED:INGREDIENT',
              'c_tooltip' : 'Medications ordered by ingredient name',
              'children':   0,
              'tty':        None,
              'rxcui':      None }
        self.build_letter_metadata(ingredient_path_prefix) # 'A'..'Z', with undetermined ingredient subfolder
        self.build_ingredient_metadata(ingredient_path_prefix) # folders for each letter, and ingredients for each
        self.build_drug_level_metadata()
        self.build_ndc_level_metadata()
        return self.metadata_rows # Done
    # end build

    def check_for_improper_path(self, build_step):
        if '\\PCORI\\MEDICATION\\INGREDIENT\\U\\100000004\\' in self.metadata_rows:
            raise ValueError('Found invalid row in [%s]' % build_step)

    def show_hlevels(self):
        hlevels = defaultdict(int) # default value is zero
        for x in self.metadata_rows.keys(): hlevels[self.metadata_rows[x]['c_hlevel']] += 1 # determine per-c_hlevel counts
        for x in sorted(hlevels.keys()): print('Level %2d ==> %5d' % (x, hlevels[x])) # display level counts

    def rxcui_type_name(self, tty):
        return 'Ingredient' if tty in ['IN','MIN','PIN'] else 'Orderable Drug'

    @staticmethod
    def get_first_letter(name):  # class level function, no instance (no self argument)
        # return first letter from name argument (in uppercase), None if no alphabetic character
        first_letter = None
        for idx, char in enumerate(name):
            next_char = name.strip()[idx].upper()
            if next_char.isalpha(): first_letter = next_char; break
        return first_letter

    def build_letter_metadata(self, path_prefix):
        ''' create metadata for letters A-Z and for 'undetermined ingredients' A-Z '''
        for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ': # 'A', 'B', ..., 'Z'
            my_parent_path = path_prefix
            self.metadata_rows[my_parent_path]['children'] += 1
            my_path = my_parent_path + letter + '\\'
            self.metadata_rows[my_path] = \
                {'c_fullname': my_path,
                 'c_hlevel': self.metadata_rows[my_parent_path]['c_hlevel'] + 1,
                 'c_name': letter,
                 'c_basecode': 'SCILHS:MED:INGREDIENT',
                 'c_tooltip': 'Medications (%s)' % letter,
                 'children': 0,
                 'tty': None,
                 'rxcui': None}
            # create (undetermined ingredient) row, child of letter (each letter will have this folder)
            my_parent_path = path_prefix + letter + '\\'  # same as my_path, at this moment
            self.metadata_rows[my_parent_path]['children'] += 1
            undetermined_ingredient_rxcui = self.undetermined_ingredient_rxcui_map[letter]
            my_path = self.get_undetermined_ingredient_path_for_letter(letter)
            c_basecode = ''  # (fabricated code, dont put in metadata), 'RXNORM:%d' % undetermined_ingredient_rxcui
            c_fullname = my_path
            tty = 'IN'
            self.metadata_rows[my_path] = \
                {'c_fullname': c_fullname,
                 'c_hlevel': self.metadata_rows[my_parent_path]['c_hlevel'] + 1,
                 'c_name': '(undetermined ingredient RxCUI [letter %s])' % letter,
                 'c_basecode': '',
                 'c_tooltip': '%s (RxNAV type %s)' % (self.rxcui_type_name(tty), tty),
                 'children': 0,
                 'tty': tty,
                 'rxcui': undetermined_ingredient_rxcui}
            # print(letter, my_path, my_parent_path)
        # end for letter in 'A'..'Z'
        return

    def build_ingredient_metadata(self, path_prefix):
        ''' build_ingredient_metadata '''
        fabricated_ingredient_rxcuis = self.rxnorm_coding.get_fabricated_rxcuis_for_undetermined_ingredients()
        letters = defaultdict(int)
        for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ': letters[letter] = 0 # force all keys to exist
        ingredient_rxcuis = self.ingredient_rxcui_set - fabricated_ingredient_rxcuis # exclude undetermined ingreds A-Z
        for ingredient_rxcui in ingredient_rxcuis: # IN, PIN ingredients ... MIN (mostly) excluded
            # NOTE: there are some MIN codes that are not excluded, if could not map to ingredient(s)
            name = self.rxnorm_coding.get_rxcui_name(ingredient_rxcui)
            letter = self.get_first_letter(name) # first letter in name gives lexical order
            if letter==None:
                raise ValueError('*** Invalid ingredient name for rxcui [%d], [%s], no alphabetic characters ***'
                                 % (ingredient_rxcui, name))
            letters[letter] += 1
            # a parent metadata row is guaranteed to exist, be subfolder of it
            my_parent_path = path_prefix + letter + '\\'
            self.metadata_rows[my_parent_path]['children'] += 1
            my_path = my_parent_path + str(ingredient_rxcui) + '\\'
            c_basecode = 'RXNORM:%d' % ingredient_rxcui
            c_fullname = my_path
            tty = self.rxnorm_coding.get_rxcui_tty(ingredient_rxcui) # IN/PIN
            self.metadata_rows[my_path] =\
                { 'c_fullname': c_fullname,
                  'c_hlevel':   self.metadata_rows[my_parent_path]['c_hlevel']+1,
                  'c_name':     name,
                  'c_basecode': c_basecode,
                  'c_tooltip': '%s (RxNAV type %s)' % (self.rxcui_type_name(tty), tty),
                  'children':   0,
                  'tty':        tty,
                  'rxcui':      ingredient_rxcui }
            if ingredient_rxcui in self.in_to_min_map: # has associated MIN ingredients
                # Make child rows for each associated MIN
                my_parent_path = my_path
                for min_rxcui in self.in_to_min_map[ingredient_rxcui]:
                    my_path = my_parent_path+str(min_rxcui) + '\\'
                    c_basecode = 'RXNORM:%d' % min_rxcui
                    c_fullname = my_path
                    tty = 'MIN'
                    name = self.rxnorm_coding.get_rxcui_name(min_rxcui) # should exist
                    self.metadata_rows[my_path] =\
                        {'c_fullname': c_fullname,
                         'c_hlevel': self.metadata_rows[my_parent_path]['c_hlevel'] + 1,
                         'c_name': name,
                         'c_basecode': c_basecode,
                         'c_tooltip': '%s (RxNAV type %s)' % (self.rxcui_type_name(tty), tty),
                         'children': 0,
                         'tty': tty,
                         'rxcui': ingredient_rxcui}
                # end for min_rxcui
            # end if ingredient_rxcui
        # end for ingredient_rxcui

        for letter in sorted(letters.keys()):
            print('%s : %d' % (letter, letters[letter]))
        print('After build_ingredient_metadata: %d REST API requests,\nSeconds since start: %s' \
              % (self.rxnav.get_request_count(),str(time.time()-self.start_t)))
        print('Cache hits: %s' % self.rxnav.get_cache_usage())
        self.show_hlevels()

    def get_undetermined_ingredient_path_for_letter(self, letter):
        path = self.path_prefix + letter + '\\' + str(self.undetermined_ingredient_rxcui_map[letter]) + '\\'
        return path

    def get_undetermined_ingredient_path_for_drug(self, drug_rxcui):
        # undetermined ingredient, put under (undetermined ingredient) folder for letter of ingredient
        rawname = self.rxnorm_coding.get_rxcui_rawname(drug_rxcui)
        first_letter_in_name = self.get_first_letter(rawname)
        if first_letter_in_name == None:
            raise ValueError('*** Invalid drug name for rxcui [%d], [%s], no alphabetic character ***' \
                             % (drug_rxcui, rawname))
        return self.get_undetermined_ingredient_path_for_letter(first_letter_in_name)

    def build_drug_level_metadata(self):
        print('Obtaining generic and branded drugs for ingredients from NLM (RxClass)')
        # Find generic drugs
        ingredient_paths = [x for x in self.metadata_rows.keys() if (self.metadata_rows[x]['tty'] in ['IN','PIN']) ]
        for ingredient_path in ingredient_paths:
            ingredient_rxcui = self.metadata_rows[ingredient_path]['rxcui']
            drug_rxcuis_for_ingredient = self.ingredient_to_drug_set_map[ingredient_rxcui]
            for drug_rxcui in drug_rxcuis_for_ingredient: # list of tuples of (rxcui,name,tty)
                name,tty = self.rxnorm_coding.get_rxcui_name(drug_rxcui), self.rxnorm_coding.get_rxcui_tty(drug_rxcui)
                my_parent_path = ingredient_path if ingredient_rxcui not in self.undetermined_ingredient_rxcuis else \
                                 self.get_undetermined_ingredient_path_for_drug(drug_rxcui) # handle undetermined IN
                self.metadata_rows[my_parent_path]['children'] += 1
                my_path = my_parent_path + str(drug_rxcui) + '\\'
                c_basecode = 'RXNORM:%d' % drug_rxcui
                c_fullname = my_path
                self.metadata_rows[my_path] =\
                    { 'c_fullname': c_fullname,
                      'c_hlevel':   self.metadata_rows[my_parent_path]['c_hlevel']+1,
                      'c_name':     name,
                      'c_basecode': c_basecode,
                      'c_tooltip': '%s (RxNAV tty:%s)' % (self.rxcui_type_name(tty), tty),
                      'children':   0,
                      'tty':        tty,
                      'rxcui':      drug_rxcui }

        print('After build_drug_level_metadata: %d REST API requests,\nSeconds since start: %s' %
              (self.rxnav.get_request_count(),str(time.time()-self.start_t)))
        print('Cache hits: %s' % self.rxnav.get_cache_usage())
        self.show_hlevels()
    # end build_drug_level_metadata

    def build_ndc_level_metadata(self):
        print('Computing NDCs for branded and generic drugs')
        drug_paths = [x for x in self.metadata_rows.keys()
                      if (self.metadata_rows[x]['tty'] in ['SBD','BPCK','SCD','GPCK']) ]
        for drug_path in drug_paths:
            my_parent_path, sbd_rxcui = drug_path, self.metadata_rows[drug_path]['rxcui'] # was -- = tup
            ndcs_for_sbd_rxcui = self.rxnav.get_ndcs_for_drug(sbd_rxcui) # list of ndc codes
            for ndc in ndcs_for_sbd_rxcui:
                child_path = my_parent_path+ndc+'\\'
                child_name = self.ndc_api.get_ndc_name(ndc)
                c_basecode = 'NDC:%s' % ndc
                c_fullname = child_path
                if child_name == str(ndc): child_name = '(%s) %s' % (ndc, self.metadata_rows[my_parent_path]['c_name'])
                self.metadata_rows[my_parent_path]['children'] += 1
                self.metadata_rows[child_path] =\
                    { 'c_fullname': c_fullname,
                      'c_hlevel': self.metadata_rows[my_parent_path]['c_hlevel']+1,
                      'c_name': child_name,  # descriptive name for NDC when available
                      'c_basecode': c_basecode,
                      'c_tooltip': 'Package for Orderable Drug %s' % self.metadata_rows[my_parent_path]['c_basecode'],
                      'children': 0,
                      'tty': 'NDC',
                      'rxcui': None }
        print('After NDCs: %d REST API requests,\nSeconds since start: %s' %
              (self.rxnav.get_request_count(),str(time.time()-self.start_t)))
        print('Cache hits: %s' % self.rxnav.get_cache_usage())
        self.show_hlevels()
    # end build_ndc_level_metadata

# end class ingredient_metadata_builder

