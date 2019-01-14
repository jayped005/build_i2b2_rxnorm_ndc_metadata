from __future__ import print_function
import multiprocessing
import sys, io, time, datetime, requests, json, signal
from collections import defaultdict
from collections import deque

'''
Module: rxnav_rest_api_mp.py

Author:
    Jay Pedersen, UNMC, Pathology/Microbiology Department, August 16, 2018

Purpose:
    Define rxnav_rest_api_rxcui_tables_mp class, which interfaces with the NLM's REST API for RxNorm,
    for the purpose of determining the 'allrelated' results for all RxCUI codes.
        
    Results are maintained in a cache file, which only the Cache Writer writes.
        
    The wrinkle with this library, is that it is meant to be used simultaneoulsy by multiple processes,
    but only one of the processes is writing to the cache file, but all can read from the cache file.

    Do worker processes know the current state of the data in the cache? They don't.
    The design is that they should not find any results in the cache and always make the
    REST API request, and then always send the result to the Cache Writer.
    
    The workers will initialize with is_worker=True, causing rxnav_get_data to
    forward results to the Cache Writer process.
    
    The workers will read the initial contents of the cache file, so that they don't
    request information that was in the cache at start time.  This is to allow for
    restart situations where the cache is not empty.
    
    This this library does not automatically write results to the cache file when a new result
    is returned from the REST API.  That result is returned to the caller who forwards it to the
    Cache Writer process to write to the cache file.
    
    The get_nxnav_data routine's interface also sligbhtly differs, instead of writing to the
    cache file when a new result is determined ... it sends the result to the Cache Writer.
  
    The cache file can then be post-processed to create the RXCUI_RELATED table for the RxCUI tables.
'''

def chomp(s): return s.rstrip('\n').rstrip('\r')

def make_utf8(s):
    ''' make_utf8: for python 2, return unicode(s) otherwise return s '''
    return unicode(s) if sys.version_info[0] < 3 else s

class DelayedKeyboardInterrupt(object):
    '''
    Citation:
    https://stackoverflow.com
      /questions/842557/how-to-prevent-a-block-of-code-from-being-interrupted-by-keyboardinterrupt-in-py
    '''
    def __enter__(self):
        self.signal_received = False
        self.old_handler = signal.signal(signal.SIGINT, self.handler)

    def handler(self, sig, frame):
        self.signal_received = (sig, frame)
        #logging.debug('SIGINT received. Delaying KeyboardInterrupt.')

    def __exit__(self, type, value, traceback):
        signal.signal(signal.SIGINT, self.old_handler)
        if self.signal_received:
            self.old_handler(*self.signal_received)
# end DelayedKeyboardInterrupt class

class rxnav_rest_api_mp():
# RxNav REST API interface class (NLM service interface)

    def __init__(self,
                 cache_filename,
                 local_functions,
                 logfile,
                 cache_writer_queue,
                 readonly_access_to_cache=False,
                 forward_result_to_cache_writer=False,
                 fail_if_not_in_cache=False): # constructor
        ''' Constructor '''

        def load_existing_cache(logfile):
            # Process all information already in the cache file ==> fill self.nxnav_cache
            # Each cache entry consists of 3 lines
            #   1. URL (REST API url whose result is cached)
            #   2. Date (date of information, YYYYMMDD)
            #   3. JSON result
            print('Reading existing cache', file=logfile); logfile.flush()
            cache_entries = 0
            line_number = 0
            f = self.rxnav_cache_file  # f refers to cache file in this method
            f.seek(0)  # top-of-file
            while True:
                file_position_of_cache_entry = f.tell() # file pos of cache entry ==> (URL,file-pos) in dictionary
                # Read next cache entry (3 lines)
                rawlines = [ x for x in [ f.readline() for idx_temp in range(3) ] if len(x) > 0 ]
                if len(rawlines)==0: break # acceptable ending to cache file, multiple of 3 lines
                if len(rawlines)!=3: raise ValueError('*** Cache file format error, not in groups of 3 lines ***')
                line_number += 3
                # Extract data ==> L1 is REST API url, L2 is data date (YYYYMMDD), L3 is JSON result
                rest_api_url, data_date = [ chomp(rawlines[i]) for i in [0,1] ]
                if len(data_date) != 8:
                    raise ValueError('*** Cache file format, date not YYYYMMDD, line %d ***' % (line_number -1))
                # Cache this entry ==> NOT the JSON result, keep memory usage down (would be massive)
                # ==> if cache hit later, read JSON result from the file, based on file position stored in cache
                # ==> keyed access (fast)
                # ==> That is why we dont bother to chomp(rawlines[2]), but need to read past it to prep for next entry.
                self.rxnav_cache[rest_api_url] = (file_position_of_cache_entry, data_date)
                cache_entries += 1
                if cache_entries % 10000 == 0:
                    print('..Read %d entries from cache' % cache_entries, file=logfile)
                    logfile.flush()
            # done reading, at EOF, we can write new entries here
            print('[Done reading cache, found %d existing entries]' % cache_entries, file=logfile)
            logfile.flush()
        # end load_existing_cache

        self.readonly_access_to_cache = readonly_access_to_cache
        self.forward_result_to_cache_writer = forward_result_to_cache_writer
        self.fail_if_not_in_cache = fail_if_not_in_cache
        self.use_caching = cache_filename is not None # Worker process uses static cache, cache contents at start
        self.local_functions = local_functions
        self.logfile = logfile
        self.cache_writer_queue = cache_writer_queue
        self.start_t = time.time()
        self.last_request_displayed_time = time.time() # used in first display, after 500 requests
        self.last_request_displayed_count = 0
        self.rest_api_timing = deque([], maxlen=500)
        self.rxnav_cache = {} # cache of REST API url => (data date, file-position)
        if cache_filename is not None: # even if is_worker
            self.rxnav_cache_file = io.open(cache_filename,
                                            ('r' if self.readonly_access_to_cache else 'a+'),
                                            encoding='utf-8') # worker only has read-access, all have tell(), seek()
            load_existing_cache(self.logfile) # read any existing information from the cache file
        else:
            self.rxnav_cache_file = None
        self.rxnav_cache_file_at_EOF = True
        self.rxnav_cache_hits = 0
        self.request_count = 0 # total requests, some/all can be resolved from cache
        self.rxnav_request_count = 0 # REST API requests
        self.rxnorm_to_name = {} # rxnorm to name mappings determined
        self.rxnorm_to_ndcs = {} # rxnorm to NDC code mappings
        self.ingred_rxcui_to_prescript_rxcuis = {}
        self.ingred_rxcui_to_drug_rxcuis = {}
        self.min_rxcui_to_ingredients = {}
        self.cache_hits = hlevels = defaultdict(int)

    # end constructor

    def get_timestamp_string(self):
        return time.strftime("%Y-%m-%d %H:%M:%S")
    # end timestamp_string

    def get_rxnav_data(self, request_url):
        self.request_count += 1 # total requests, not only those that go to REST API

        # try to get data from cache
        if self.use_caching and (request_url in self.rxnav_cache):
            self.rxnav_cache_hits += 1
            file_position_of_cache_entry = self.rxnav_cache[request_url][0] # tuple in cache (file_position, data_date)
            self.rxnav_cache_file.seek(file_position_of_cache_entry)
            self.rxnav_cache_file_at_EOF = False # Track that position is moved
            # At this point we are positioned to the 3 line cache entry -- the 3rd line is the JSON result
            for idx in [1,2]: self.rxnav_cache_file.readline() # skip past 1st 2 lines in cache entry
            json_text = chomp(self.rxnav_cache_file.readline()) # JSON text, r.text value from REST API call
            return json.loads(json_text) # convert JSON text to python structure

        # Data is NOT in cache, must request from NLM's REST API
        if self.fail_if_not_in_cache: # if dont want to access REST API ==> fail
            raise ValueError('RxNAV data NOT In cache for [%s]' % request_url)

        # attempt up to 40 times, sleeping 15 seconds between retries (10 minutes max)
        rest_api_request_start_time = time.time()
        retry_limit = 40
        r = None
        for idx in range(retry_limit):
            try:
                r = requests.get(request_url) # JSON text
                break # request completed ==> stop retrying
            except requests.exceptions.RequestException as e:
                print('[Communications error with RxNav REST API ... attempt %d of 3%s]' \
                      % (idx+1,', retrying in 30 seconds' if idx<retry_limit else ''),
                      file=self.logfile)
                print('Communication Error: [%s]' % str(e), file=self.logfile)
                print('NOTE: RxNav requests: %d, seconds since start: %s' %
                      (self.rxnav_request_count,str(time.time()-self.start_t)),
                      file=self.logfile)
                self.logfile.flush()
                time.sleep(15) # wait 15 seconds before retry
                pass
        # end comm retry loop
        if r == None:
            raise ConnectionError # no response after max retries
        rest_api_request_duration = time.time() - rest_api_request_start_time
        self.rest_api_timing.append(rest_api_request_duration)
        # Request completed successfully, break from retry loop occurred
        result = json.loads(r.text) # convert JSON text into python structure
        self.rxnav_request_count += 1 # variable in procedure that we are nested in
        #if self.rxnav_request_count % 39 == 0: time.sleep(1) # limit around 20 requests/second
        if self.rxnav_request_count % 500 == 0:
            batch_size = self.rxnav_request_count-self.last_request_displayed_count
            if self.rxnav_request_count == 500: # DEBUG
                print('[Timings of first 500 requests]', file=self.logfile)
                timings = list(self.rest_api_timing)
                for idx in range(0, len(timings), 20):
                    print(str([round(x, 3) for x in timings[idx:idx + 20]]), file=self.logfile)
                print('[End timings]', file=self.logfile)
                self.logfile.flush()
            print('[%s] Sum of request timings of last batch of %d ==> %s seconds' \
                  % (self.get_timestamp_string(), batch_size, str(sum(self.rest_api_timing))),
                  file=self.logfile)
            seconds_since_last_display = time.time()-self.last_request_displayed_time # floating point result
            print('[%s] RxNav requests: %d, REST API calls: %d, seconds: %s, rate/sec: %s, cache (size: %d, hits: %d)' \
                  % (self.get_timestamp_string(),
                     self.request_count,
                     self.rxnav_request_count,
                     str(seconds_since_last_display),
                     str(round(batch_size/seconds_since_last_display, 3)),
                     len(self.rxnav_cache),
                     self.rxnav_cache_hits),
                  file=self.logfile)
            self.logfile.flush()
            self.last_request_displayed_time = time.time()
            self.last_request_displayed_count = self.rxnav_request_count
        if self.forward_result_to_cache_writer: # Send result to Cache Writer
            self.cache_writer_queue.put([request_url, r.text]) # send list with url and json result strings
        # return result
        return result # python structure generated from JSON result
    # end get_rxnav_data

    def write_cache_entry(self, request_url, json_string):
        # write 3 lines to end of cache file -- URL, data_date, JSON result
        if not self.rxnav_cache_file_at_EOF:
            self.rxnav_cache_file.seek(0, 2)  # seek to EOF
            self.rxnav_cache_file_at_EOF = True
        # NOTE: rxnav_cache_file_at_EOF is guaranteed to be True here, no need to reset
        # JGP 2018-05-22 - write all 3 lines at once, the position of the first line is the cache position
        file_position = self.rxnav_cache_file.tell()
        data_date = '{:%Y%m%d}'.format(datetime.datetime.now())
        print(make_utf8(request_url+'\n'+data_date+'\n'+json_string), file=self.rxnav_cache_file)
        self.rxnav_cache[request_url] = (file_position, data_date) # track position of JSON in file in the cache
    # end write_cache_entry

    def get_class_tree(self, classId):  # eg: VA root is VA000 as of Aug 6, 2018 (per Lee Peters) (was N0000010574)
        ''' Main use is determining VA drug class hierarchy (aka NDFRT). '''
        d = self.get_rxnav_data('https://rxnav.nlm.nih.gov/REST/rxclass/classTree/json?classId=%s' % classId)
        return d
    # end get_class_tree

    def get_generic_drugs_for_VA_class(self, va_classid):
        # relaSource changed to VA from NDFRT, per Lee Peters, after August 6, 2018
        # also -- https://rxnav.nlm.nih.gov/RxClassAPIREST.html#uLink=RxClass_REST_getClassMembers
        # NOTE: can't simply change this to get_allrelated -- special curation not contained in allrelated
        query_url = ('https://rxnav.nlm.nih.gov/REST/rxclass/classMembers.json?classId=%s'+\
                     '&relaSource=VA&rela=has_VAClass&ttys=SCD+GPCK') % va_classid
        d = self.get_rxnav_data(query_url)
        return ([] if not ('drugMemberGroup' in d and 'drugMember' in d['drugMemberGroup'])
                else [int(x['minConcept']['rxcui'])
                      for x in d['drugMemberGroup']['drugMember']])
    # end get_generic_drugs_for_VA_class

    def get_ndcs_for_drug(self, drug_rxcui): # See "Part 4" above for JSON spec
        d = self.get_rxnav_data('https://rxnav.nlm.nih.gov/REST/rxcui/%d/allhistoricalndcs/json' % drug_rxcui)
        # Comprehension below generates a list of list of NDC codes (due to JSON structure) ==> flatten to list
        return ([] if not ('historicalNdcConcept' in d
                           and d['historicalNdcConcept'] \
                           and 'historicalNdcTime' in d['historicalNdcConcept'])
                else
                    self.local_functions.flatten(
                        [x['ndc']
                         for x in self.local_functions.flatten(
                            [x['ndcTime']
                             for x in d['historicalNdcConcept']['historicalNdcTime']])]))
    # end get_ndcs_for_drug

    def get_allrelated(self, rxcui):
        '''
        Return the 'allrelated' data structure from the NLM
        '''
        query_url = '''https://rxnav.nlm.nih.gov/REST/rxcui/%s/allrelated.json''' % rxcui
        d = self.get_rxnav_data(query_url)
        # Example JSON result returned from NLM, returned to the caller as a python structure
        '''
        {"allRelatedGroup":
         {"rxcui":"1049214",
          "conceptGroup":
           [{"tty":"BN",
             "conceptProperties":
              [{"rxcui":"216903",
                "name":"Endocet",
                "synonym":"",
                "tty":"BN",
                "language":"ENG",
                "suppress":"N",
                "umlscui":"C0720206"},
               ]
            }
           ]
         }
        }
        '''

        return d
    # end get_allrelated

    def get_historical_rxcui(self, rxcui):
        '''
        Return historical information about the given RxCUI from the NLM, using their rxcuihistory API.
        '''
        query_url = '''https://rxnav.nlm.nih.gov/REST/rxcuihistory/concept.json?rxcui=%s''' % rxcui
        d = self.get_rxnav_data(query_url)
        # Example JSON result returned from NLM, returned to the caller as a python structure
        '''
        {"rxcuiHistoryConcept":
         {"rxcuiConcept":
            {"status":"Retired",
             "rxcui":"991041",
             "tty":"SBD",
             "str":"Chlorpromazine hydrochloride 10 MG Oral Tablet [Thorazine]",
             "sab":"RXNORM",
             "doseform":"Oral Tablet",
             "doseformRxcui":"317541",
             "mid":"0",
             "branded":"1",
             "qf":"",
             "qd":"",
             "startDate":"062010",
             "endDate":"022013",
             "isCurrent":"0",
             "currentRxcui":"",
             "packAlias":"",
             "scdName":"Chlorpromazine hydrochloride 10 MG Oral Tablet",
             "scdRxcui":"991039"
            },
          "bossConcept":
            [
             {"baseRxcui":"2403",
              "baseName":"Chlorpromazine",
              "bossRxcui":"104728",
              "bossName":"Chlorpromazine hydrochloride",
              "actIngredRxcui":"",
              "actIngredName":"",
              "moietyRxcui":"",
              "moietyName":"",
              "numeratorValue":"10",
              "numeratorUnit":"MG",
              "denominatorValue":"",
              "denominatorUnit":""
             }
            ]
          }
        }
        '''
        return d
    # end get_historical_rxcui

    def get_historical_rxcui_attributes(self, rxcui, attributes):
        # support 'NAME','TTY','STATUS','START','END','SCDRXCUI','BOSSRXCUIS'
        d = self.get_historical_rxcui(rxcui)
        if not ('rxcuiHistoryConcept' in d and 'rxcuiConcept' in d['rxcuiHistoryConcept']):
            return None
        d1 = d['rxcuiHistoryConcept']['rxcuiConcept']
        result = []
        for attr_name in attributes:
            if attr_name == 'NAME':
                result.append(d1['str'])
            elif attr_name == 'TTY':
                result.append(d1['tty'])
            elif attr_name == 'STATUS':
                result.append(d1['status'])
            elif attr_name == 'START':
                result.append(d1['startDate'])
            elif attr_name == 'END':
                result.append(d1['endDate'])
            elif attr_name == 'SCDRXCUI':
                result.append(int(d1['scdRxcui']) if (d1['scdRxcui'] and len(d1['scdRxcui']) > 0) else None)
            elif attr_name=='BOSSRXCUIS':
                bossrxcuis = []
                if 'bossConcept' in d['rxcuiHistoryConcept']:
                    for b in d['rxcuiHistoryConcept']['bossConcept']:
                        if 'bossRxcui' in b and len(b['bossRxcui']) > 0:
                            bossrxcuis.append( int(b['bossRxcui']))
                result.append(bossrxcuis)
            else:
                raise ValueError('Inavlid attribute [%s] in get_historical_rxcui_attributes' % attr_name)
        return tuple(result)
    # end get_historical_rxcui_attributes

    def determine_category_from_tty(self, tty):
        # two categories of interest for now -- INGREDIENT and DRUG
        category =       'INGREDIENT' if tty in ['IN', 'MIN', 'PIN'] \
                    else 'DRUG'       if tty in ['SCD', 'SBD', 'GPCK', 'BPCK'] \
                    else 'OTHER'
        return category
    # end determine_category_from_tty

    def get_historical_rxcuis(self,
                              target_status_values=["ACTIVE", "RETIRED", "NON-RXNORM", "NEVER%20ACTIVE"],
                              verbose=False):
        ''' Use NLM's rxcuihistory API to obtain RxCUI values for various classes of RxNorm codes (active, et al) '''
        target_status_values_set = set(target_status_values)
        valid_status_values = ["ACTIVE", "RETIRED", "NON-RXNORM", "NEVER%20ACTIVE"]
        valid_status_values_set = set(valid_status_values)
        intersection = target_status_values_set.intersection(valid_status_values_set)
        if len(intersection) != len(target_status_values_set):
            raise ValueError('get_historical_rxcuis -- invalid target_status_values, valid values: %s' % str(valid_status_values))
        rxcui_set = set()
        rxcuis_status_d = {}
        prev_size = 0
        for status_value in target_status_values:
            query_url = 'https://rxnav.nlm.nih.gov/REST/rxcuihistory/status.json?type=%s' % status_value
            d = self.get_rxnav_data(query_url) # {"rxcuiList": {"rxcuis": ["211", "292", ...] } }
            code_set = set([ int(x) for x in d['rxcuiList']['rxcuis'] ])
            rxcui_set.update(code_set)
            for x in code_set: rxcuis_status_d[x] = status_value
            if verbose:
                print(query_url, file=self.logfile)
                print('get_historical_rxcuis: [%s] : %d, delta:%d, tot: %d' \
                      % (status_value, len(code_set), len(rxcui_set)-prev_size, len(rxcui_set)),
                      file=self.logfile)
                self.logfile.flush()
            prev_size = len(rxcui_set)
        return rxcui_set, rxcuis_status_d # set
    # end get_historical_rxcuis

    def get_allrelated_rxcuis_for_tty_list(self, d, tty_list):
        result = []
        for idx in range(len(d['allRelatedGroup']['conceptGroup'])):
            d1 = d['allRelatedGroup']['conceptGroup'][idx]
            if 'conceptProperties' in d1:
                level2 = d1['conceptProperties']
                if d1['tty'] in tty_list:
                    result.extend([int(x['rxcui']) for x in level2])
        return result
    # end get_allrelated_rxcuis_for_tty_list

    def get_next_allrelated_tty_and_list(self, d):
        # process allrelated result ... return next segment of codes -- return tty and list
        # ==> iterator (uses yield to return next result
        # ==> for tty, rxcui_list in get_next_related_tty_and_list(d): <do-something-with-tty-rxcui_list>
        for idx in range(len(d['allRelatedGroup']['conceptGroup'])):
            d1 = d['allRelatedGroup']['conceptGroup'][idx]
            if 'conceptProperties' in d1:
                level2 = d1['conceptProperties']
                if isinstance(level2, list):
                    yield d1['tty'], [int(x['rxcui']) for x in level2]
                else:
                    raise ValueError('*** Non-List in allrelated result [%s] ***' % str(level2))
    # end get_next_allrelated_tty_and_list

    def get_related_ingredients_for_multi_ingredient(self, min_rxcui):
        # JGP - 2018-08-17, make this work with the 'allrelated' result instead of 'related'
        # ==> find the related 'IN' or 'PIN' codes associated with the specified 'MIN' code.
        # ==> we have 'allrelated' for all RxNORM RxCUI codes, it will be in the cache.

        d = self.get_allrelated(min_rxcui) # get 'allrelated'
        return self.get_allrelated_rxcuis_for_tty_list(d, ['IN', 'PIN'])
    # end get_related_ingredients_for_multi_ingredient

    def get_drug_rxcuis_for_ingredient(self, ingredient_rxcui):
        d = self.get_allrelated(ingredient_rxcui)
        return self.get_allrelated_rxcuis_for_tty_list(d, ['SBD', 'SCD', 'GPCK', 'BPCK'])
    # end get_drug_rxcuis_for_ingredient

    def get_ingredients_for_generic_drug(self, scd_rxcui):
        d = self.get_allrelated(scd_rxcui)
        return self.get_allrelated_rxcuis_for_tty_list(d, ['IN'])
    # end get_ingredients_for_generic_drug

    def get_branded_drugs_for_generic_drug(self, scd_rxcui):
        d = self.get_allrelated(scd_rxcui)
        return self.get_allrelated_rxcuis_for_tty_list(d, ['SBD','BPCK'])
    # end get_branded_drugs_for_generic_drug

    def get_ndc_codes_for_drug(self, drug_rxcui):
        d = self.get_rxnav_data('https://rxnav.nlm.nih.gov/REST/rxcui/%d/allhistoricalndcs/json' % drug_rxcui)
        # Comprehension below generates a list of list of NDC codes (due to JSON structure) ==> flatten to list
        if not ('historicalNdcConcept' in d \
                and d['historicalNdcConcept'] \
                and 'historicalNdcTime' in d['historicalNdcConcept']):
            result_ndcs = []
        else:
            result_ndcs = self.local_functions.flatten(\
                [x['ndc']
                 for x in self.local_functions.flatten([x['ndcTime']
                                                        for x in d['historicalNdcConcept']['historicalNdcTime']])])
        return set(result_ndcs)
    # end get_ndc_codes_for_drug

    def get_request_count(self):
        return self.rxnav_request_count
    # end get_request_count

    def get_cache_usage(self):
        return str(self.rxnav_cache_hits) #','.join(['%s:%d' % (nm,self.cache_hits[nm]) for nm in self.cache_hits])
    # end get_cache_usage

# end class
