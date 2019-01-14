from __future__ import print_function
from collections import defaultdict

'''
Author:
    Jay Pedersen, UNMC Pathology/Microbiology, March 10, 2018
    
Abstract:
    Define rxnorm_code_maps class to ...

    Define the information model used by the ingredient metadata building, including
    ALL historical ingredients and drugs.  This information can not currently be
    determined by the RxNorm REST API provided by the NLM.
    
    Use the current RxCUI information returned form the RxNorm REST API, but also
    include each and every ingredient and drug coded defined in the RxCUI history
    file received at UNMC from the NLM on Nov 1, 2017.  THis file included the
    comprehensive set of RxCUI codes as of the end of Oct 2017.  This will later
    be replaced by a RxNorm REST API call that the NLM has agreed to provide.
    
    Use the RxNorm REST API to determine the set of active drugs associated
    with each active ingredient.  This is bolstered by adding historical
    drug codes which match exactly to the ingredient by name.
    
    For any drug code that can not be matched to an existing or historical ingredient code,
    attach it to 'ingredient undetermined' code, given RXCUI 1 billion.

Methods:
    get_ingredient_set()
    get_drug_set()
    get_multi_ingredient_set()
    ingredient_to_drug_set(ingredient_rxcui)
    drug_to_ingredient_set(drug_rxcui)
    show_unmapped_drugs() -- drugs not associated with ingredients
'''

class rxnorm_code_maps:

    '''
    "getters" which expose the information model of RxNorm ingredients and drugs
    '''
    def get_ingredient_rxcui_set(self):       return self.ingredient_rxcui_set
    def get_multi_ingredient_rxcui_set(self): return self.min_ingredient_rxcui_set
    def get_drug_rxcui_set(self):             return self.drug_rxcui_set

    def get_ingredient_to_drug_set_map(self): return self.ingredient_to_drug_set
    def get_drug_to_ingredient_set_map(self): return self.drug_to_ingredient_set
    def get_ingredient_to_min_set_map(self):  return self.ingredient_to_min_set

    def get_rxcui_name(self,rxcui):           return self.rxcui_attr[rxcui][self.attr_tups_d['NAME']]
    def get_rxcui_rawname(self,rxcui):        return self.rxcui_attr[rxcui][self.attr_tups_d['RAWNAME']]
    def get_rxcui_tty(self,rxcui):            return self.rxcui_attr[rxcui][self.attr_tups_d['TTY']]
    def get_rxcui_status(self,rxcui):         return self.rxcui_attr[rxcui][self.attr_tups_d['STATUS']]
    def get_rxcui_start_date(self,rxcui):     return self.rxcui_attr[rxcui][self.attr_tups_d['START']]
    def get_rxcui_end_date(self,rxcui):       return self.rxcui_attr[rxcui][self.attr_tups_d['END']]

    def get_undetermined_ingredient_rxcui_map(self):    return self.undetermined_ingredient_rxcui
    def get_undetermined_ingredient_rxcuis(self):       return self.undetermined_ingredient_rxcui.values()
    def get_undetermined_ingredient_rxcui(self,letter): return self.undetermined_ingredient_rxcui[letter]
    def get_unmapped_min_code_set(self): return self.unmapped_min_code_set
    def get_fabricated_rxcuis_for_undetermined_ingredients(self): return self.fabricated_rxcuis_for_undetermined_ingredient

    '''
    CONSTRUCTOR.  Interfaces (indirectly) with NLM's RxNorm REST API, by processing informatino
                  in sets and dictionaries created from rxnav_rest_api.
    '''
    def __init__(self, rxnav, ingredient_rxcui_set, drug_rxcui_set, be_verbose=False):
        if be_verbose: print('--- rxnorm_code_maps constructor ---')
        self.rxnav = rxnav # data already in cache
        self.rxcui_attr = {} # track TTY,RAWNAME,NAME,STATUS,START,END,SCDRXCUI
        self.hist_drug_to_ingredient_attr = {}
        self.drug_to_ingredient_tups_d = {} # track PROVENANCE,SBDRXCUI,IN_RXCUI
        self.attr_tups_order = ['PROVENANCE', 'STATUS', 'RAWNAME', 'NAME', 'TTY', 'START', 'END']
        self.attr_tups_d = {nm: idx for idx, nm in enumerate(self.attr_tups_order)}
        self.drug_to_ingredient_tups_order = ['PROVENANCE','BOSSRXCUIS','SCDRXCUI','IN_RXCUIS']
        self.drug_to_ingredient_tups_d = { nm : idx for idx,nm in enumerate(self.drug_to_ingredient_tups_order) }
        self.ingredient_tty_list = ['IN', 'MIN', 'PIN']
        self.drug_tty_list = ['SCD','GPCK','SBD','BPCK']
        self.valid_tty_list = self.ingredient_tty_list+self.drug_tty_list
        self.ingredient_rxcui_set = None
        self.min_ingredient_rxcui_set = None
        self.drug_rxcui_set = None
        self.drug_to_ingredient_set = {} # empty unless build_ingredients_to_drug_set_map called
        self.ingredient_to_drug_set = {} # empty unless build_ingredients_to_drug_set_map called
        self.ingredient_to_min_set = defaultdict(set)  # empty unless build_ingredients_to_drug_set_map called
        self.ingredient_to_drug_provenance = {} # provenance for ingredient-to-drug map, set, HISTAPI/NAMEPARSING/UNDET
        self.drug_to_ingredient_provenance = {} # provenance for drug-to-ingredient map
        self.undetermined_ingredient_rxcui = {} # on a per-letter basis, e.g. undetermined but starting with letter D
        self.unmapped_min_code_set = set()
        self.fabricated_rxcuis_for_undetermined_ingredient = set()

        # determine the sets of ingredient and drug RxCUI codes
        ingredient_rxcuis = ingredient_rxcui_set # set([x for x in rxcui_d if rxcui_d[x]['type']] in ['IN','MIN','PIN'])
        drug_rxcuis = drug_rxcui_set # set([x for x in rxcui_d if rxcui_d[x]['type']] in ['SCD','SBD','GPCK','BPCK'])
        hist_rxcuis = ingredient_rxcuis | drug_rxcuis # union

        # ==> create self.rxcui_tups_d with per-code information about TTY, NAME, RAWNAME, STATUS, START, END
        print('--- Get Historical RxCUI attributes ---')
        for rxcui in hist_rxcuis: # all historical codes of interest
            # information already in the rxnav cache
            tup = self.rxnav.get_historical_rxcui_attributes\
                (rxcui,['NAME','TTY','STATUS','START','END','SCDRXCUI','BOSSRXCUIS'])
            name_str, tty_str, status_str, start_str, end_str, scdrxcui, bossrxcuis = tup
            rawname_str = name_str
            if status_str == 'Retired':
                name_str = '(retired %s) %s' % (end_str,name_str)
            elif status_str == 'Never Active':
                name_str = '(never active) %s' % name_str
            self.rxcui_attr[rxcui] = ('HISTAPI', status_str, rawname_str, name_str, tty_str, start_str, end_str)
            # Set IN_RXCUI to None for Drug codes, rest are ingredient codes ... use their own rxcui
            if tty_str in self.drug_tty_list:
                self.hist_drug_to_ingredient_attr[rxcui] = ('HISTAPI', bossrxcuis, scdrxcui, bossrxcuis)

        # Add catchall ingredients for historic drugs whose ingredients can not be determined (e.g. Never-Active)
        # use RXNORM:<100-million> as the code, far in the future for RxCUIs, TEMP HACK
        self.undetermined_ingredient_rxcui_set = set()
        next_undetermined_ingredient_rxcui = 100000000 - 1 # first will be 100M, loop increments by 1
        for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ': # 'A', then 'B', ...
            next_undetermined_ingredient_rxcui += 1
            self.fabricated_rxcuis_for_undetermined_ingredient.add(next_undetermined_ingredient_rxcui)
            self.undetermined_ingredient_rxcui_set.add(next_undetermined_ingredient_rxcui)
            self.undetermined_ingredient_rxcui[letter] = next_undetermined_ingredient_rxcui
            self.rxcui_attr[next_undetermined_ingredient_rxcui] =\
                ('HIST', 'Retired', '(undetermined ingredients)', '(undetermined ingredients)', 'IN', '', '')

        # determine ingredients, min_ingredients and drugs -- ALL must have metadata
        # NOTE: UNDETERMINED INGREDIENT is in set
        self.ingredient_rxcui_set = set(rxcui
                                        for rxcui in self.rxcui_attr
                                        if self.rxcui_attr[rxcui][self.attr_tups_d['TTY']] in ['IN', 'PIN'])
        self.min_ingredient_rxcui_set = set(rxcui
                                            for rxcui in self.rxcui_attr
                                            if self.rxcui_attr[rxcui][self.attr_tups_d['TTY']] in ['MIN'])
        self.drug_rxcui_set = set(rxcui
                                  for rxcui in self.rxcui_attr
                                  if self.rxcui_attr[rxcui][self.attr_tups_d['TTY']] in self.drug_tty_list)

        # determine the ingredient to drug mappings
        print('--- Determine ingredient to drug mappings ---')
        self.build_ingredients_to_drug_set_map()

        # done
        print('--- rxnorm_code_maps constructor (END) ---')
    # end of constructor

    '''
    Methods used by the constructor to build the information model of RxNorm ingredients and drugs
    '''
    def find_ingredient_from_drug_name(self,name_in):
        # e.g. Brompheniramine 2 MG/ML
        name = name_in
        for x in ['24 HR ', '12 HR ']:
            if name.startswith(x):
                name = name[len(x):].strip()
        parsed_name = name.split(' ')
        # find position of first numeric
        # e.g. lamotrigine 100 MG Extended Release Tablet
        #      ['lamotrigine', '100', 'MG', 'Extended', 'Release', 'Tablet']
        found = False;
        pos = 0
        for idx, s in enumerate(parsed_name):
            if s.replace('.', '').isdigit():
                found = True; pos = idx; break
        ingredient_name = ' '.join(parsed_name[:pos]) if found else ' '.join(parsed_name)
        return ingredient_name
    # end find_ingredient_from_drug_name

    def build_ingredients_to_drug_set_map(self):

        def get_first_letter(name):  # class level function, no instance (no self argument)
            # return first letter from name argument, None if no alphabetic character
            first_letter = None
            for idx, char in enumerate(name):
                next_char = name.strip()[idx].upper()
                if next_char.isalpha():
                    first_letter = next_char; break
            return first_letter

        # Prep for parsing -- get canonical ingredient names (uppercase), exact match checking for drug name parsing
        self.ingredient_name_to_rxcui_d = {}
        for rxcui in self.ingredient_rxcui_set: # IN, PIN codes
            provenance, name = [self.rxcui_attr[rxcui][self.attr_tups_d[x]]
                                for x in ['PROVENANCE', 'RAWNAME']]
            name = name.upper() # canonical form is uppercase
            if name not in self.ingredient_name_to_rxcui_d:
                self.ingredient_name_to_rxcui_d[name] = rxcui
            elif provenance=='API':
                print('## API ingredient RxCUI %d overriding in ingredient name to RxCUI map ##' % rxcui)
                self.ingredient_name_to_rxcui_d[name] = rxcui

        all_drugs_set = self.drug_rxcui_set
        all_ingredients = self.ingredient_rxcui_set # IN/PIN, not MIN
        never_active_drugs_set = set([x
                                      for x in self.rxcui_attr
                                      if (self.rxcui_attr[x][self.attr_tups_d['STATUS']] == 'Never Active'
                                          and self.rxcui_attr[x][self.attr_tups_d['TTY']] in self.drug_tty_list)])
        min_ingredient_set = self.min_ingredient_rxcui_set # MIN codes, multi-ingredient

        print('... Building ingredient to drug set map ...')
        print('Found %d MIN ingredient codes (multi ingredient)' % len(min_ingredient_set))
        print('Found %d drug codes' % len(all_drugs_set))
        print('Found %d Never Active drug codes' % len(never_active_drugs_set))

        # Initialize mapping sets and provenances
        for drug_rxcui in all_drugs_set:
            self.drug_to_ingredient_set[drug_rxcui] = set()
            self.drug_to_ingredient_provenance[drug_rxcui] = set()
        for ingredient_rxcui in all_ingredients:
            self.ingredient_to_drug_set[ingredient_rxcui] = set()
            self.ingredient_to_drug_provenance[ingredient_rxcui] = set()

        # Create IN to MIN map, from REST API and parsing the name and mapping it to ingredients directly
        empty_related_IN_sets = 0
        min_codes_mapped_to_ingredients = set()
        for min_ingredient_rxcui in min_ingredient_set:
            # Note: ingredient_to_min_set is a defaultdict(set), don't need to check if key exists defaults to empty set
            related_ingredient_rxcuis = self.rxnav.get_related_ingredients_for_multi_ingredient(min_ingredient_rxcui)
            if len(related_ingredient_rxcuis) == 0: empty_related_IN_sets += 1
            for ingredient_rxcui in related_ingredient_rxcuis:
                self.ingredient_to_min_set[ingredient_rxcui].add(min_ingredient_rxcui)  # add MIN code for this IN code
                min_codes_mapped_to_ingredients.add(min_ingredient_rxcui)
            # This is a MIN rxcui -- parse NAME to determine ingredients, e.g. Streptodornase / Streptokinase
            name = self.rxcui_attr[min_ingredient_rxcui][self.attr_tups_d['RAWNAME']]
            result_names = name.split(' / ')
            for name in result_names:
                ingredient_name = self.find_ingredient_from_drug_name(name)
                if ingredient_name in self.ingredient_name_to_rxcui_d:
                    ingredient_rxcui = self.ingredient_name_to_rxcui_d[ingredient_name]
                    self.ingredient_to_min_set[ingredient_rxcui].add(min_ingredient_rxcui)  # add MIN code for this IN code
                    min_codes_mapped_to_ingredients.add(min_ingredient_rxcui)
            # end for name
        # end for min_ingredient_rxcui

        # The unmapped MIN codes will be treated like IN ingredients, since there are no ingredient
        # codes to tie them to.  They are stand-alone.  This set is small, a set of 233 as of May 2018. (JGP)
        # There are over 4,000 MIN codes which are mapped to one or more ingredients (presumably at least 2)
        # and are not treated as stand-alone.
        self.unmapped_min_code_set = min_ingredient_set - min_codes_mapped_to_ingredients
        print('[MIN codes: %d, mapped to IN/PIN: %d, unmapped: %d, zero-related-IN+PIN cases: %d'
              % (len(min_ingredient_set), len(min_codes_mapped_to_ingredients),
                 len(self.unmapped_min_code_set), empty_related_IN_sets))

        # Use the historical drug to ingredient mapping
        for drug_rxcui in all_drugs_set:
            if drug_rxcui in self.hist_drug_to_ingredient_attr: # might be API-determined drug
                bossrxcuis = self.hist_drug_to_ingredient_attr[drug_rxcui][self.drug_to_ingredient_tups_d['BOSSRXCUIS']]
                #scdrxcui = self.hist_drug_to_ingredient_tups_d[rxcui][self.drug_to_ingredient_tups_d['SCDRXCUI']]
                for ingredient_rxcui in bossrxcuis:
                    if ingredient_rxcui not in self.ingredient_rxcui_set:
                        the_tty = self.rxcui_attr[ingredient_rxcui][self.attr_tups_d['TTY']] \
                                if ingredient_rxcui in self.rxcui_attr else '<unknown-RxCUI>'
                        print(('*** ISSUE -- for drug RxCUI %d, tty[%s] -- ' \
                               + 'history BOSSRXCUI %d -- NOT a known ingredient (ignore) ***') \
                              % (drug_rxcui,the_tty,ingredient_rxcui))
                    else:
                        self.drug_to_ingredient_set[drug_rxcui].add(ingredient_rxcui)
                        self.drug_to_ingredient_provenance[drug_rxcui].add('HISTAPI') # on per-drug RXCUI basis
                        self.ingredient_to_drug_set[ingredient_rxcui].add(drug_rxcui)
                        self.ingredient_to_drug_provenance[ingredient_rxcui].add('HISTAPI')

        # Parse these drugs to determine the ingredient linkages ==> hist_ingredient_to_drug_set_map
        for drug_rxcui in all_drugs_set: # never_active_drugs_set:
            self.parse_drug_name_and_map_to_ingredients(drug_rxcui)

        # Find ingredients to drug mapping by REST API
        all_ingredients = self.ingredient_rxcui_set
        for ingredient_rxcui in [rxcui for rxcui in all_ingredients
                                 if rxcui not in self.undetermined_ingredient_rxcui_set]: # RxCUI ingredient
            drug_rxcuis = self.rxnav.get_drug_rxcuis_for_ingredient(ingredient_rxcui)
            for drug_rxcui in drug_rxcuis:
                if not (drug_rxcui in self.rxcui_attr \
                        and self.rxcui_attr[drug_rxcui][self.attr_tups_d['TTY']] in self.drug_tty_list):
                    the_tty = self.rxcui_attr[drug_rxcui][self.attr_tups_d['TTY']] \
                              if drug_rxcui in self.rxcui_attr else '<unknown-RxCUI>'
                    print(('*** ISSUE -- for ingredient RxCUI %d -- get_drug_rxcuis returned ' \
                           +' %d, tty:[%s] -- NOT a known drug (ignore) ***') \
                          % (ingredient_rxcui, drug_rxcui, the_tty))
                else:
                    self.ingredient_to_drug_set[ingredient_rxcui].add(drug_rxcui)
                    self.ingredient_to_drug_provenance[ingredient_rxcui].add('API')
                    self.drug_to_ingredient_set[drug_rxcui].add(ingredient_rxcui)
                    self.drug_to_ingredient_provenance[drug_rxcui].add('API')
            # end for drug_rxcui
        # end for ingredient_rxcui

        # deal with no-result ==> add to undetermined-ingredient for the first letter of the drug name
        for drug_rxcui in all_drugs_set:
            if len(self.drug_to_ingredient_set[drug_rxcui]) == 0:
                rawname = self.rxcui_attr[drug_rxcui][self.attr_tups_d['RAWNAME']]
                first_letter = get_first_letter(rawname)
                undetermined_ingredient_rxcui = self.undetermined_ingredient_rxcui[first_letter]
                self.ingredient_to_drug_set[undetermined_ingredient_rxcui].add(drug_rxcui)
                self.ingredient_to_drug_provenance[undetermined_ingredient_rxcui].add('UNDETERMINED_INGREDIENT')
                self.drug_to_ingredient_set[drug_rxcui].add(undetermined_ingredient_rxcui)
                self.drug_to_ingredient_provenance[drug_rxcui].add('UNDETERMINED_INGREDIENT')
        return

    def parse_drug_name_and_map_to_ingredients(self,drug_rxcui):

        tty, name = [self.rxcui_attr[drug_rxcui][self.attr_tups_d[x]].upper().strip()
                     for x in ['TTY', 'RAWNAME']]
        # By-ingredient parsing
        # NOTE: Use RAWNAME, NAME can be modified, e.g. '(RETIRED 2013-02) CHLORPROMAZINE HYDROCHLORIDE 10 MG ...'
        # whereas RAWNAME is 'CHLORPROMAZINE HYDROCHLORIDE 10 MG ...'.
        result_names = []
        if tty in ['GPCK','BPCK']:
            # syntax - e.g. '{7 (lamotrigine 100 MG Extended Release Tablet) / .. / .. ) }
            if name[0]=='{':
                brace_idx = name.find('}')
                if brace_idx > 0: # if cant find this, poorly formatted and cant parse it
                    result_names = []
                    name = name[1:brace_idx]
                    names1 = [x for x in name.split(')') if len(x.strip()) > 0]
                    for n1 in names1:
                        # e.g. 7 (24 HR PRAMIPEXOLE DIHYDROCHLORIDE 0.375 MG EXTENDED RELEASE TABLET [MIRAPEX] / ..
                        # drop leading number, e.g. 7 (...
                        if n1.startswith(' / '): n1 = n1[len(' / '):]
                        temp = n1.split(' ')
                        if len(temp) > 0 and temp[0].isdigit(): n1 = ' '.join(temp[1:])
                        n1 = n1.lstrip('(') # drop leading parenthesis
                        # process the rest,  e.g. 24 HR PRAMIPEXOLE DIHYDROCHLORIDE 0.375 MG / ...
                        names2 = [x for x in n1.split(' / ') if len(x.strip()) > 0]
                        if len(names2) > 0: result_names.extend(names2)
            else:
                result_names.append(name)
        else: # ['SCD','SBD']
            # syntax - e.g. '24 HR lamotrigine 200 MG Extended Release Capsule'
            #          e.g. 'Nitrofurantoin 100 MG Oral Capsule [Nitro Macro]'
            #          e.g. 'Brompheniramine 2 MG/ML / Phenylephrine 0.25 MG/ML Oral Solution'
            result_names.extend(name.split(' / '))
        # At this point we have: rxcui of drug and a ingredient names, map to ingredient RXCUI values
        for name in result_names:
            ingredient_name = self.find_ingredient_from_drug_name(name)
            if ingredient_name in self.ingredient_name_to_rxcui_d:
                ingredient_rxcui = self.ingredient_name_to_rxcui_d[ingredient_name]
                self.drug_to_ingredient_set[drug_rxcui].add(ingredient_rxcui)
                self.drug_to_ingredient_provenance[drug_rxcui].add('NAMEPARSING')
                self.ingredient_to_drug_set[ingredient_rxcui].add(drug_rxcui)
                self.ingredient_to_drug_provenance[ingredient_rxcui].add('NAMEPARSING')
        # end for name
    # end parse_drug_name_and_map_to_ingredients

    '''
    utility methods
    '''
    def format_hist_date(self,datestr): # deal with date strings from the RxNav REST API
        ''' e.g. '22015' (len 5) ==> 2015-02, '102015' (len 6) ==> 2015-10 '''
        return datestr[-4:]+'-0'+datestr[0] if len(datestr)==5 \
          else datestr[-4:]+'-'+datestr[:2] if len(datestr)==6 \
          else datestr # unknown format, return as-is
    # end format_hist_date

    '''
    The following methods are "getters"
    '''

    def get_multi_ingredient_rxcuis_for_ingredient(self, ingredient_rxcui):
        '''
        Match signature of corresponding routine from rxnav_rest_api class.
        Return tuples with (rxcui,name,tty), rxcui as an integer
        '''
        return [] if ingredient_rxcui not in self.ingredient_to_min_set \
          else [tuple([drug_rxcui] + [self.rxcui_attr[drug_rxcui][self.attr_tups_d[x]]
                                      for x in ['NAME', 'TTY']])
                for drug_rxcui in self.ingredient_to_min_set[ingredient_rxcui]]

    def get_min_rxcuis_for_ingredient(self, ingredient_rxcui):
        '''
        Match signature of corresponding routine from rxnav_rest_api class.
        Return list of rxcuis of MIN ingredients associated with this IN ingredient.
        '''
        return [] if ingredient_rxcui not in self.ingredient_to_min_set \
          else self.ingredient_to_min_set[ingredient_rxcui]

    def get_drug_rxcuis_for_ingredient(self, ingredient_rxcui):
        '''
        Match signature of corresponding routine from rxnav_rest_api class.
        Return tuples with (rxcui,name,tty), rxcui as an integer
        '''
        return ([] if ingredient_rxcui not in self.ingredient_to_drug_set
                else [tuple([drug_rxcui] + [self.rxcui_attr[drug_rxcui][self.attr_tups_d[x]]
                                                    for x in ['NAME', 'TTY']])
                      for drug_rxcui in self.ingredient_to_drug_set[ingredient_rxcui]])

    '''
    The following methods are reporting methods, "dumpers", not used in metadata building
    '''

    def show_undetermined_ingredient_drugs(self,status='Retired'):
        print('--- UNDETERMINED INGREDIENT ---')
        for ingredient_rxcui in self.undetermined_ingredient_rxcui.values(): # rxcuis for letter A,B,...,Z
            ingred_name, ingred_tty, ingred_status = [self.rxcui_attr[ingredient_rxcui][self.attr_tups_d[x]]
                                                      for x in ['NAME','TTY','STATUS']]
            provenance = ','.join(x for x in self.ingredient_to_drug_provenance[ingredient_rxcui])
            for idx,drug_rxcui in enumerate(self.ingredient_to_drug_set[ingredient_rxcui]):
                drug_name, drug_tty, drug_status = [self.rxcui_attr[drug_rxcui][self.attr_tups_d[x]]
                                                    for x in ['NAME', 'TTY','STATUS']]
                bossrxcuis,scdrxcui = [ self.hist_drug_to_ingredient_attr\
                                            [drug_rxcui][self.drug_to_ingredient_tups_d[x]]
                                        for x in ['BOSSRXCUIS','SCDRXCUI'] ]

                print('%s : "%s" (%s) (%s) (%d)' \
                      % ((('"%s" (%s) (%s) (%s) (%d)' % (ingred_name,ingred_tty,ingred_status,provenance,ingredient_rxcui))
                          if idx == 0 else \
                          ' '*(len(ingred_name)+9+len('%d' % ingredient_rxcui)+len(ingred_tty)+len(ingred_status) \
                                                       +len(provenance))),
                          drug_name,drug_tty,drug_status,drug_rxcui))
                print(self.rxcui_attr[drug_rxcui])
                print(bossrxcuis)
                print(scdrxcui)
            # end drug iteration
        # end ingredient iteration
        print('--- UNDETERMINED INGREDIENT (END) ---')
    # end show_undetermined_ingredient_drugs

    def show_ingredient_to_drug_map(self):
        '''
        Example output:
        1921069 : [1921082, 1921083]
        abaloparatide : abaloparatide 2 MG/ML Pen Injector (SCD)
                      : abaloparatide 2 MG/ML Pen Injector [Tymlos] (SBD)
        '''
        print('--- INGREDIENT TO DRUG ---')
        for ingredient_rxcui in self.ingredient_to_drug_set:
            print('%s : %s' % (ingredient_rxcui, str(list(self.ingredient_to_drug_set[ingredient_rxcui]))))
            ingred_name, ingred_tty, ingred_status = [self.rxcui_attr[ingredient_rxcui][self.attr_tups_d[x]]
                                                      for x in ['NAME','TTY','STATUS']]
            provenance = ','.join(x for x in self.ingredient_to_drug_provenance[ingredient_rxcui])
            for idx,drug_rxcui in enumerate(self.ingredient_to_drug_set[ingredient_rxcui]):
                drug_name, drug_tty = [self.rxcui_attr[drug_rxcui][self.attr_tups_d[x]]
                                       for x in ['NAME', 'TTY']]
                print('%s : %s (%s) (%d)' \
                      % ((('%s (%s) (%s) (%s) (%d)' % (ingred_name,ingred_tty,ingred_status,provenance,ingredient_rxcui))
                          if idx == 0 else \
                          ' '*(len(ingred_name)+9+len('%d' % ingredient_rxcui)+len(ingred_tty)+len(ingred_status) \
                                                       +len(provenance))),
                         drug_name,drug_tty,drug_rxcui))
        print('--- INGREDIENT TO DRUG (END) ---')
    # end show_ingredient_to_drug_map

    def show_drug_to_ingredient_map(self):
        print('--- DRUG TO INGREDIENT ---')
        for drug_rxcui in self.drug_rxcui_set:
            drug_name, drug_tty, drug_status = [self.rxcui_attr[drug_rxcui][self.attr_tups_d[x]]
                                                for x in ['NAME','TTY','STATUS']]
            provenance = ','.join(x for x in self.drug_to_ingredient_provenance[drug_rxcui]) # set to string
            ingredient_rxcuis = self.drug_to_ingredient_set[drug_rxcui]
            for idx, ingredient_rxcui in enumerate(ingredient_rxcuis):
                ingredient_name, ingredient_tty = [self.rxcui_attr[ingredient_rxcui][self.attr_tups_d[x]]
                                                   for x in ['NAME','TTY']]
                print('%s : %s (%s) (%d)' % ((('%s (%s) (%s) (%s) (%d)' \
                                               % (drug_name, drug_tty, drug_status, provenance, drug_rxcui))
                                               if idx == 0 else \
                                               ' ' * (len(drug_name) + 9 + len('%d' % ingredient_rxcui) \
                                                      + len(drug_tty) + len(drug_status) + len(provenance))),
                                             ingredient_name, ingredient_tty, ingredient_rxcui))
        print('--- DRUG TO INGREDIENT (END) ---')
    # end show_drug_to_ingredient_map

    def show_ingredient_to_min_map(self):
        print('--- INGREDIENT TO MIN ---')
        for ingredient_rxcui in self.ingredient_to_min_set:
            print('%s : %s' % (ingredient_rxcui, str(list(self.ingredient_to_min_set[ingredient_rxcui]))))
            ingred_name, ingred_tty, ingred_status = [self.rxcui_attr[ingredient_rxcui][self.attr_tups_d[x]]
                                                      for x in ['NAME','TTY','STATUS']]
            for idx,min_rxcui in enumerate(self.ingredient_to_min_set[ingredient_rxcui]):
                min_name, min_tty, min_status = [self.rxcui_attr[min_rxcui][self.attr_tups_d[x]]
                                                 for x in ['NAME', 'TTY','STATUS']]
                print('%s : %s (%s) (%s) (%d)' \
                      % ((('%s (%s) (%s) (%d)' % (ingred_name,ingred_tty,ingred_status,ingredient_rxcui))
                          if idx == 0 else \
                          ' '*(len(ingred_name)+7+len('%d' % ingredient_rxcui)+len(ingred_tty)+len(ingred_status))),
                         min_name,min_tty,min_status,min_rxcui))
        print('--- INGREDIENT TO MIN (END) ---')
    # end show_ingredient_to_min_map

    def show_never_active_drug_coverage(self):
        print('--- NEVER ACTIVE DRUGS COVERAGE ---')
        never_active_drugs_set = set([x
                                      for x in self.rxcui_attr
                                      if (self.rxcui_attr[x][self.attr_tups_d['STATUS']] == 'Never Active'
                                          and self.rxcui_attr[x][self.attr_tups_d['TTY']] in self.drug_tty_list)])
        never_active_drugs_coverage = 0
        for rxcui in never_active_drugs_set:
            if rxcui in self.drug_to_ingredient_set:
                ingredient_rxcuis = self.drug_to_ingredient_set[rxcui]
                if len(ingredient_rxcuis) > 0:
                    never_active_drugs_coverage += 1
        print('Never active drugs code coverage: %d of %d are linked to known ingredients' %
              (never_active_drugs_coverage, len(never_active_drugs_set)))
        print('--- NEVER ACTIVE DRUGS COVERAGE (END) ---')
    # end show_never_active_drug_coverage

# end class