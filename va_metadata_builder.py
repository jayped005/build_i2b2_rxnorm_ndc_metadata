import time
from collections import defaultdict
from ndc_api import ndc_api

'''
Module: va_metadata_builder.py

Author: Jay Pedersen, UNMC, Pathology/Microbiology Department, Feb 9, 2017

Purpose:
    (a) Define va_metadata_builder class to build the VA hierarchy drug class
        for use in SCILHS metadata.
'''

class va_metadata_builder():
# --------------------------------------------------------------------------------------|
#         VA Drug Classification -- Metadata Build  (Olivier Bodenrieder input)         |
# --------------------------------------------------------------------------------------|

    def __init__(self, rxnav, rxnorm_coding):
        self.rxnav = rxnav
        self.rxnorm_coding = rxnorm_coding
        self.start_t = None # set by build method
        self.ndc_api = ndc_api()
        self.metadata_rows = {}
        self.classid_paths = {}
        self.scd_for_va_classid = {} # generic drugs associated with each VA class
        self.ingred_for_va_classid = {} # ingredient rxcui associated with each VA class
        self.scd_for_ingred_for_va_classid = {} # generic drugs associated with each ingredient, specific for VA class
        self.ingredient_rxcui_to_paths = {}
        self.rootpath = None # set in build_va_folders
        self.ndfrt_name = {} # given NDRFT code, return associated name, set in build_va_folders

    def show_hlevels(self):
        hlevels = defaultdict(int) # default value is zero
        for x in self.metadata_rows.keys(): hlevels[self.metadata_rows[x]['c_hlevel']] += 1
        for x in sorted(hlevels.keys()): print('Level %2d ==> %5d' % (x, hlevels[x]))

    def determine_va_readable_path(self, mypath):
        readable_path = ''
        # determine the c_name for the first element in the VA drug class path
        # For example, the path
        #   \PCORI\MEDICATION\RXNORM_CUI\N0000029074\
        #     Is for "Antimicrobials
        #   \PCORI\MEDICATION\RXNORM_CUI\N0000029074\N0000029080\"
        #     Is "Chloramphenicol" child of "Antimicrobials"
        # Find the c_name for the first two NDFRT codes in the path
        pathlist = mypath.split('\\') # eg: ['PCORI','MEDICATION','RXNORM_CUI','N0000029074']
        va_classes = [x for x in pathlist if x.startswith('N') and (x[1:]).isdigit()] # N<digits> are NDFRT codes
        ndfrt_list = va_classes[:3] # get up to three NDFRT codes
        readable_path = ''
        for list_idx,ndfrt_code in enumerate(ndfrt_list):
            if len(readable_path)==0: readable_path = ', VA: '
            readable_path += (', ' if list_idx > 0 else '') + self.ndfrt_name[ndfrt_code]
        return readable_path

    def build_va_folders(self, tlist, level, parent_path, rootpath):
        self.rootpath = rootpath
        for x in tlist:
            myclassid = x['rxclassMinConceptItem']['classId'] # eg: N0000010574, which are NDFRT codes, per Jeff Klann ==> NDFRT:<code>
            if parent_path!=rootpath:
                self.metadata_rows[parent_path]['children'] += 1
                mypath = parent_path+myclassid+'\\'
                name = x['rxclassMinConceptItem']['className'].lower().capitalize()
                self.ndfrt_name[myclassid] = name
            else:
                mypath = parent_path+'RXNORM_CUI'+'\\'  # Match Jeff Klann, want \PCORI\MEDICATION\RXNORM_CUI\% paths
                name = 'VA Drug Classes' # Match Jeff Klann, replace 'VA Classes (VA)'
            # x['rxclassMinConceptItem'] example ==> {u'classId': u'N0000029360', u'className': u'MULTIVITAMINS', u'classType': u'VA'}
            # ==> Use className -- but lowercase and then capitalize (only uppercase first letter)
            c_basecode = 'VACLASS:%s' % myclassid # no longer NDFRT, per Lee Peters, Aug 2018
            c_fullname = mypath
            readable_path = self.determine_va_readable_path(c_fullname)
            self.metadata_rows[mypath] = \
                { 'c_fullname': c_fullname,
                  'c_hlevel': level,
                  'c_name': name,
                  'c_basecode': c_basecode,
                  'c_tooltip': 'VA drug class%s' % readable_path,
                  'classId': myclassid,
                  'children': 0,
                  'tty': 'VAclass',
                  'rxcui': None }
            self.classid_paths[myclassid] = mypath
            if 'rxclassTree' in x: # build next level of tree (recursive, depth-first)
                self.build_va_folders(x['rxclassTree'], level+1, mypath, rootpath)
        return

    def rxcui_type_name(self, tty):
        return 'Ingredient' if tty=='IN' else 'Orderable Drug'

    def add_rxcui_metadata_child(self, child_rxcui, my_parent_path, child_path, child_name, tty):
        self.metadata_rows[my_parent_path]['children'] += 1
        c_basecode = 'RXNORM:%d' % child_rxcui
        c_fullname = child_path
        hlevel = self.metadata_rows[my_parent_path]['c_hlevel']+1
        readable_path = self.determine_va_readable_path(c_fullname)
        if hlevel>12: print('[HIST3] hlevel=%d, parent [%s]\n                    child [%s]' %
                            (hlevel,my_parent_path,child_path)) # DEBUG
        self.metadata_rows[child_path] = \
            { 'c_fullname': c_fullname,
              'c_hlevel':   hlevel,
              'c_name':     child_name,
              'c_basecode': c_basecode,
              'c_tooltip': '%s (RxNAV tty:%s)%s' % (self.rxcui_type_name(tty), tty, readable_path),
              'children':   0,
              'tty':        tty,
              'rxcui':      child_rxcui }

    def find_generic_drugs_for_VA_classes(self):
        print('Obtaining generic drugs for VA classes from NLM (RxClass)')
        for va_classid in self.classid_paths.keys():
            self.scd_for_va_classid[va_classid] = set() # only care at leaf level
            if self.metadata_rows[ self.classid_paths[va_classid] ]['children'] == 0: # leaf
                result_rxcuis = self.rxnav.get_generic_drugs_for_VA_class(va_classid)
                for rxcui in result_rxcuis:
                    self.scd_for_va_classid[va_classid].add(rxcui)
        # end for va_classid
    # end find_generic_drugs_for_VA_classes

    def find_ingredients_for_generic_drugs(self):
        print('Obtaining ingredients for generic drugs returned for VA classes from NLM (RxClass)')
        for va_classid in sorted(self.classid_paths.keys()):
            self.ingred_for_va_classid[va_classid] = set()
            for scd_rxcui in self.scd_for_va_classid[va_classid]: # next generic drug associated with class
                ingredient_rxcuis = self.rxnav.get_ingredients_for_generic_drug(scd_rxcui)
                for ingred_rxcui in ingredient_rxcuis: # (rxcui,name,tty) of next ingredient for this generic drug
                    self.ingred_for_va_classid[va_classid].add(ingred_rxcui) # ingredient for this VA class
                    # Want to know the class-specific set of generic drugs associated with the ingredient.
                    # ==> we are going through the list of generic drugs for the class one by one
                    # ==> the generic drug defined by 'scd' is a generic drug associated with the ingredient
                    # ==> we are processing for this class (and associated at the class level).
                    if va_classid not in self.scd_for_ingred_for_va_classid:
                        self.scd_for_ingred_for_va_classid[va_classid] = {}
                    if ingred_rxcui not in self.scd_for_ingred_for_va_classid[va_classid]:
                        self.scd_for_ingred_for_va_classid[va_classid][ingred_rxcui] = set()
                    self.scd_for_ingred_for_va_classid[va_classid][ingred_rxcui].add(scd_rxcui)
                    # self.ingredient_rxcui_to_paths -- needed for historical code processing
                    if ingred_rxcui not in self.ingredient_rxcui_to_paths:
                        self.ingredient_rxcui_to_paths[ingred_rxcui] = set()
                    self.ingredient_rxcui_to_paths[ingred_rxcui]\
                        .add(self.classid_paths[va_classid]+str(ingred_rxcui)+'\\')
        return

    def find_associated_branded_drugs(self):
        print('Determine associated branded drugs for generic drugs returned for VA classes from NLM (RxClass)')
        # Find generic drugs
        generic_drug_paths = [x for x in self.metadata_rows.keys() if (self.metadata_rows[x]['tty'] in ['SCD','GPCK']) ]
        for generic_drug_path in generic_drug_paths:
            scd_rxcui = self.metadata_rows[generic_drug_path]['rxcui']
            # See if there is a branded drug associated with the generic drug whose RXCUI is held in scd
            branded_drug_rxcuis = self.rxnav.get_branded_drugs_for_generic_drug(scd_rxcui)
            for sbd_rxcui in branded_drug_rxcuis:
                name = self.rxnorm_coding.get_rxcui_name(sbd_rxcui)
                tty = self.rxnorm_coding.get_rxcui_tty(sbd_rxcui)
                my_parent_path = '\\'.join(['']+[x for x in generic_drug_path.split('\\')
                                                 if len(x) > 0][:-1]+['']) # ingredient path, was -- generic_drug_path
                child_name = name
                child_path = my_parent_path+str(sbd_rxcui)+'\\'
                if child_path not in self.metadata_rows: # may have already been seen
                    self.add_rxcui_metadata_child(sbd_rxcui, my_parent_path, child_path, child_name, tty)
        # end for generic_drug_path
    # end find_associated_branded_drugs

    def find_ndc_codes_for_drugs(self):
        print('Computing NDCs for branded and generic drugs')
        drug_paths = [x for x in self.metadata_rows.keys()
                      if (self.metadata_rows[x]['tty'] in ['SBD','BPCK','SCD','GPCK']) ]
        for drug_path in drug_paths:
            my_parent_path, sbd_rxcui = drug_path, self.metadata_rows[drug_path]['rxcui'] # was -- = tup
            ndcs_for_sbd_rxcui = self.rxnav.get_ndcs_for_drug(sbd_rxcui) # list of ndc codes
            for ndc in ndcs_for_sbd_rxcui:
                child_path = my_parent_path+ndc+'\\'
                child_name = self.ndc_api.get_ndc_name(ndc)
                if child_name == str(ndc): child_name = '(%s) %s' % (ndc, self.metadata_rows[my_parent_path]['c_name'])
                self.metadata_rows[my_parent_path]['children'] += 1
                c_basecode = 'NDC:%s' % ndc
                c_fullname = child_path
                readable_path = self.determine_va_readable_path(c_fullname)
                self.metadata_rows[child_path] = \
                    { 'c_fullname': c_fullname,
                      'c_hlevel': self.metadata_rows[my_parent_path]['c_hlevel']+1,
                      'c_name': child_name,  # descriptive name for NDC when available
                      'c_basecode': c_basecode,
                      'c_tooltip': 'Package for Orderable Drug %s%s' %
                                   (self.metadata_rows[my_parent_path]['c_basecode'], readable_path),
                      'children': 0,
                      'tty': 'NDC',
                      'rxcui': None }
        print('After NDCs: %d REST API requests,\nSeconds since start: %s' %
              (self.rxnav.get_request_count(),str(time.time()-self.start_t)))
        print('Cache hits: %s' % self.rxnav.get_cache_usage())
        self.show_hlevels()

    def build(self, path_prefix, metadata_root_level):

        self.start_t = time.time()

        # (Levels 1 to 5) Determine VA drug hierarchy
        print('Computing VA drug hierarchy metadata from classTree obtained from RxNav')
        VA_classId = 'VA000' # Per Lee Peters -- as of Aug 6, 2018 -- NLM supports VA class ids and not NDFRT ids
        d = self.rxnav.get_class_tree(VA_classId) # obtain VA hierarchy tree from NLM, VA hierarchy root 'VA000'
        root_path = '\\'+path_prefix+'\\'
        self.build_va_folders(d['rxclassTree'], metadata_root_level, root_path, root_path)
        print('hlevel values after compute VA hierarchy'); self.show_hlevels()
        print('After VA hierarchy: %d REST API requests,\nSeconds since start: %s' %
              (self.rxnav.get_request_count(),str(time.time()-self.start_t)))
        print('Cache hits: %s' % self.rxnav.get_cache_usage())

        # find path to the MISCELLANEOUS AGENTS folder, classId 'N0000029353'
        misc_agents_classId = 'XX000' # hardcode of MISCELLANEOUS AGENTS class (was NDFRT 'N0000029353')
        misc_agents_path = self.classid_paths[misc_agents_classId]
        print('c_name for MISCELLANEOUS classid [%s] is [%s]' %
              (misc_agents_classId,self.metadata_rows[misc_agents_path]['c_name']))

        # Determine generic drugs for the VA classes
        self.find_generic_drugs_for_VA_classes()
        print('After find generic drugs for VA: %d REST API requests,\nSeconds since start: %s' %
              (self.rxnav.get_request_count(),str(time.time()-self.start_t)))
        print('Cache hits: %s' % self.rxnav.get_cache_usage())

        # Determine ingredients associated with generic drugs
        self.find_ingredients_for_generic_drugs()
        print('After find ingredients for generic drugs for VA: %d REST API requests,\nSeconds since start: %s' %
              (self.rxnav.get_request_count(),str(time.time()-self.start_t)))
        print('Cache hits: %s' % self.rxnav.get_cache_usage())

        # Create metadata rows for ingredients and generic drugs
        for va_classid in sorted(self.classid_paths.keys()):
            va_classpath = self.classid_paths[va_classid]
            for ingred_rxcui in self.ingred_for_va_classid[va_classid]:
                my_parent_path = va_classpath
                child_name = self.rxnorm_coding.get_rxcui_name(ingred_rxcui)
                child_path = my_parent_path+str(ingred_rxcui)+'\\'
                self.add_rxcui_metadata_child(ingred_rxcui, my_parent_path, child_path, child_name, 'IN')
                for scd_rxcui in self.scd_for_ingred_for_va_classid[va_classid][ingred_rxcui]:
                    my_parent_path = va_classpath+str(ingred_rxcui)+'\\'
                    child_name = self.rxnorm_coding.get_rxcui_name(scd_rxcui)
                    child_path = my_parent_path+str(scd_rxcui)+'\\'
                    self.add_rxcui_metadata_child(scd_rxcui, my_parent_path, child_path, child_name, 'SCD')
        print('Created metadata rows for ingredients and generic drugs'); self.show_hlevels()

        # (Level 6) add in historical RxNorm generic drug codes, based on ingredients
        # JGP -- 2018-03-01, suspending this for now.
        # self.place_historical_rxnorm_codes(misc_agents_path)
        # We are NOT missing any RXCUI codes by not doing this, as our historical comprehensiveness
        # comes from the by-ingredient hierarchy and not the VA hierarchy.
        # This method not currently sustainable (no way to determine these replacements versus
        # obsoletes going forward, and knowoledge of which drug replaced which drug.
        # There is no way to know if a historical code is one that should be associated with the VA class --
        # as that is curated.  We don't have a history of that curation.

        # Find branded drugs associated with generic drugs
        self.find_associated_branded_drugs()
        print('Created metadata rows for branded drugs'); self.show_hlevels()

        # (Level 8) Determine NDC codes for branded and generic drugs
        self.find_ndc_codes_for_drugs()

        return self.metadata_rows # Done

    def get_ingredient_rxcui_set(self):
        return set(self.ingredient_rxcui_to_paths.keys())

# end class va_metadata_builder
