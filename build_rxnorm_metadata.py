#!/usr/bin/python

from __future__ import print_function
import multiprocessing as mp
import sys, io, math, time

# Classes specific to this project, in separate python scripts in the same folder as this script.
from utility_functions import utility_functions
from rxnav_rest_api_mp import rxnav_rest_api_mp

'''
Module: build_rxnorm_metadata.py

Author:
    Jay Pedersen, UNMC, Pathology/Microbiology Department, September 2018

Purpose:
    Build i2b2 medications metadata including RxNORM and NDC code information
    and the VA drug class hierarchy.

    The information for these tables is extracted from the NLM's RxNORM REST API.
    
    Simultaneous processing of REST API requests is performed using python's
    multiprocessing library.  This reduces the issue created by the significant
    latency of obtaining a result from the REST API by performing multiple
    simultaneous requests from different processes.  This allows us to process
    more RxCUI codes per unit of time.

Architecture:    
    Create cache file with 'historic' and 'allrelated' information for every RxCUI code,
    which is the bulk of the information needed for the RxCUI tables.

    Extract the set of NDC codes for drug RxCUIs using the 'allhistoricalndcs' endpoint
    from the NLM's REST API.

    Extract the VA drug classification hierarchy, needed for building 'VA class' metadata.
    That information is also extracted into a CSV file.

    Use multiple processes simultaneously requesting disjoint pieces of the
    information from REST API (via python's multiprocessing library).
    This overcomes some of the latency problem associated with the REST API where
    each single piece of information is obtained in isolation, requring a round-trip
    request and response ... which takes a mean of 0.25 seconds for the request+response.
    For example, using 4 simultaneously executing processes, we can obtain about 16 pieces of
    information per second instead of 4.

Processing algorithm:

    Manager process:

      Create manager.log

      Main process determines the set of RxCUI codes from NLM Rest API,
      which is around 650,000 codes as of Augest 2018.  It then whittles
      this list down to about 350,000 codes of interest by excluding
      the 'Non-RxNORM' RxCUI codes.  There are about 300,000 of those.

      Write these codes to the cache file.

      Take the remaining set of approximately 350,000 and create four equal segments
      of codes (around 90,000 per segment).  Assign each segment to a worker.

      Create (--workers <count>) worker processes, passing the associated code segment to each.     
      Create 1 Cache Writer process to handle the write requests to the cache file.
      Wait for processes to terminate.
      Stop Cache Writer

      Done.

    RxCUI Worker process: (any of a set of workers), passed a segment of the RxCUI codes

      For each RxCUI code in given segment of RxCUI codes:
          Request the 'allrelated' information from the NLM's REST API.
          Request the 'historicalrxcui' information form the NLM's REST API.
          Send the result to the Cache writer process queue.          
      Done.

    Cache Writer process: passed a queue where messgaes from workers received

      Open cache file for writing, position to end of file.     
      Forever:
          Read next message from queue
          if 'kill' message from manager then stop
          Otherwise message is from a worker => Write to cache.
       Done.      
'''

def get_timestamp_string():
    return time.strftime('%Y-%m-%d %H:%M:%S')
# end get_timestamp_string

# ------------------- BasicTree class -------------------

class BasicTree:
    # embedded Node class
    class Node:
        def __init__(self, data):
            self.data = data
            self.children = []
            self.parent_node = None

        def add_child(self, child_node):
            self.children.append(child_node)
            child_node.parent_node = self

    # end Node

    # BasicTree:

    def __init__(self):
        self.root = None  # root is a Node, root node of tree

    def add_node(self, data, parent=None):  # identifier need not be unique
        if self.root and parent == None:
            raise ValueError('BasicTree add_node -- no parent specified and Root established')
        if parent and self.root == None:
            raise ValueError('BasicTree add_node -- Parent specified and root not set')
        node = self.Node(data)  # create node, set its data
        if self.root == None:  # if no Root, this is the root, not a child
            self.root = node
        else:  # must be a child, root already established
            parent.add_child(node)
        return node

    def traverse(self, node, depth_first=True):
        if node == None:
            yield None  # empty tree, yield immediately, traversal is finished
        else:
            yield node  # yield root node, then follow DFS/BFS rules of traversal of rest of tree
            q = node.children  # children of root node
            while q:
                yield q[0]  # yield node at head of queue
                # proceed with traversal in DFS/BFS manner, queue order updated in DFS/BFS manner
                child_nodes = q[0].children
                q = (child_nodes + q[1:]) if depth_first \
                    else (q[1:] + child_nodes)  # BFS expansion, prev line is DFS expansion
# end BasicTree

# --------------------- multiprocessing process definitions ==> workers and cache writer ---------------

# -------------------         Worker task           ----------------------------

# Worker -- process a segment of RxCUI codes, find their 'allrelated' definitions
#           cause the information to be written to the cache file.

def worker_task_rxcuis(mp_queue, worker_number, cache_filename, barrier, rxcui_code_list, log_filename, utility_fns):
    ''' RxCUI Worker -- determine 'allrelated' and 'rxcuihistory' for RxCUI codes, send results to Cache Writer '''

    logfile = io.open(log_filename, 'w', encoding='utf-8')
    print('[%s] Worker %d -- processing %d RxCUI codes from [%s] to [%s]'
          % (get_timestamp_string(), worker_number, len(rxcui_code_list), rxcui_code_list[0], rxcui_code_list[-1]),
          file=logfile)
    logfile.flush()
    rxnav = rxnav_rest_api_mp(cache_filename, utility_fns, logfile, mp_queue,
                              readonly_access_to_cache=True,
                              forward_result_to_cache_writer=True)

    # At this point, we have initialized with the REST API, which means that the
    # initial contents of the cache file have been read.
    # Wait for all other workers to arrive at barrier.
    print('[%s] Waiting at the barrier' % get_timestamp_string(), file=logfile);
    logfile.flush()
    barrier.wait()
    print('[%s] Passing the barrier' % get_timestamp_string(), file=logfile);
    logfile.flush()

    for idx, rxcui in enumerate(rxcui_code_list):
        rxnav.get_allrelated(rxcui)
        # NOTE: the request caused the result to be queued to the Cache Writer
        # ==> nothing more needs to be done.
        rxnav.get_historical_rxcui(rxcui)
        # As with 'allrelated', the act of making this call sends data to the Cache Writer
        if idx % 1000 == 0:
            print('[%s] Processed %d/%d RxCUIs, last was [%s]'
                  % (get_timestamp_string(), (idx + 1), len(rxcui_code_list), rxcui_code_list[idx]),
                  file=logfile)
            logfile.flush()
    # end processing rxcui_codes

    print('[%s] Finished processing %d RxCUI codes ... terminating.' % (get_timestamp_string(), len(rxcui_code_list)),
          file=logfile)
    logfile.flush()
    logfile.close()

# end worker_task_rxcuis


def worker_task_ndcs(mp_queue, worker_number, cache_filename, barrier, drug_rxcui_list, log_filename, utility_fns):
    ''' NDC Worker -- determine NDC codes for drug RxCUI codes, results to Cache Writer '''

    logfile = io.open(log_filename, 'w', encoding='utf-8')
    print('[%s] Worker %d -- processing %d RxCUI codes from [%s] to [%s]'
          % (get_timestamp_string(), worker_number, len(drug_rxcui_list), drug_rxcui_list[0], drug_rxcui_list[-1]),
          file=logfile)
    logfile.flush()
    rxnav = rxnav_rest_api_mp(cache_filename, utility_fns, logfile, mp_queue,
                              readonly_access_to_cache=True,
                              forward_result_to_cache_writer=True)

    # At this point, we have initialized with the REST API, which means that the
    # initial contents of the cache file have been read.
    # Wait for all other workers to arrive at barrier.
    print('[%s] Waiting at barrier 2' % get_timestamp_string(), file=logfile);
    logfile.flush()
    barrier.wait()
    print('[%s] Passing barrier 2' % get_timestamp_string(), file=logfile);
    logfile.flush()

    for idx, rxcui in enumerate(drug_rxcui_list):
        rxnav.get_ndc_codes_for_drug(rxcui)
        # NOTE: the request caused the result to be queued to the Cache Writer
        # ==> nothing more needs to be done.
        # As with 'allrelated', the act of making this call sends data to the Cache Writer
        if idx % 1000 == 0:
            print('[%s] Processed %d/%d drug RxCUIs, last was [%s]'
                  % (get_timestamp_string(), (idx + 1), len(drug_rxcui_list), drug_rxcui_list[idx]),
                  file=logfile)
            logfile.flush()
    # end processing rxcui_codes

    print('[%s] Finished processing %d drug RxCUI codes ... terminating.' % (
    get_timestamp_string(), len(drug_rxcui_list)),
          file=logfile)
    logfile.flush()
    logfile.close()

# end worker_task_ndcs


# ------------------            Cache Writer            ----------------------------


def cache_writer_task(mp_queue, cache_filename, log_filename, utility_fns):
    ''' Cache Writer -- waits for cache messages on mp_queue, writes to cache file. '''

    logfile = io.open(log_filename, 'w', encoding='utf-8')
    print('[%s] Starting Cache Writer' % get_timestamp_string(), file=logfile);
    logfile.flush()
    rxnav = rxnav_rest_api_mp(cache_filename, utility_fns, logfile, mp_queue,
                              readonly_access_to_cache=False,
                              forward_result_to_cache_writer=False)

    count = 0
    while True:
        m = mp_queue.get()
        if isinstance(m, list):  # message from Worker -- [request_url, json_result]
            request_url, json_result_string = m
            rxnav.write_cache_entry(request_url, json_result_string)
            count += 1
            if count % 1000 == 0:
                print('Processed %d messages' % count, file=logfile);
                logfile.flush()
        elif m == 'kill':  # manager is saying to stop, all workers have completed
            break
    # end processing mp_queue

    print('[%s] Processed %d cache messages' % (get_timestamp_string(), count), file=logfile)
    print('[%s] Received kill message from manager process ... terminating.'
          % (get_timestamp_string(),), file=logfile)
    logfile.flush()
    logfile.close()

# end cache_writer_task

# Phase 0 task -- per RxCUI status values ==> ACTIVE/RETIRED/etc

def phase0_task(mp_queue,
                log_filename,
                cache_filename,
                utility_fns):

    logfile = io.open(log_filename, 'w', encoding='utf-8')
    print('[%s] Starting' % get_timestamp_string(), file=logfile);
    logfile.flush()
    rxnav = rxnav_rest_api_mp(cache_filename, utility_fns, logfile, mp_queue,
                              readonly_access_to_cache=True,
                              forward_result_to_cache_writer=True)
    rxnav.get_historical_rxcuis(target_status_values=["ACTIVE", "RETIRED", "NEVER%20ACTIVE","NON-RXNORM"], verbose=True)

# end phase0_task

# Phase 3 task -- VA Drug Class information

def phase3_task(mp_queue,
                log_filename,
                cache_filename,
                utility_fns):

    def get_next_id():
        ''' Use the outer function next_id[0] integer value for the "next id", and then update it '''
        next_available_id = next_id[0]  # "next id", list from outer context (build_va_drug_class_csv context)
        next_id[0] += 1  # prep sequence for next allocation, this will be the id returned on the next call
        return next_available_id
    # end get_next_id

    def build_va_folders(t, tlist, parent_node=None):
        ''' Process nested structure of JSON with class hierarchy, recursively process nested levels '''
        for x in tlist:
            va_classid = x['rxclassMinConceptItem']['classId']  # eg: VA000, VA classid codes
            # NOTE: per Lee Peters, as of Aug 6, 2018 in NLM result, NDFRT no longer supported or used by NLM
            va_classid_set.add(va_classid)  # inherit vbl from environment
            raw_va_drug_class_name = x['rxclassMinConceptItem']['className'].lower().capitalize()
            if va_classid not in va_classid_d:
                va_classid_d[va_classid] = raw_va_drug_class_name
            va_drug_class_name = 'VA Drug Classes' if not parent_node else raw_va_drug_class_name
            if parent_node: parent_node.data['children'] += 1
            new_node = t.add_node({'level_num': 0 if not parent_node else parent_node.data['level_num'] + 1,
                                   'code': va_classid,
                                   'id': get_next_id(),
                                   'parent_id': None if not parent_node else parent_node.data['id'],
                                   'children': 0,
                                   'type': 'VACLASS',
                                   'name': va_drug_class_name},
                                  parent_node)
            va_node_d[va_classid] = new_node  # inherit from environment
            if 'rxclassTree' in x:  # build next level of tree (recursive, depth-first)
                build_va_folders(t, x['rxclassTree'], new_node)
        # end for x in tlist
    # end build_va_folders

    def cache_generic_drugs_for_VA_classes(rxnav, va_classids, va_node_d):
        print('Obtaining generic drugs for VA classes from NLM (RxClass)', file=logfile)
        scd_rxcuis_for_va_classid_d = {}  # SCD RxCUIs are generic drugs, SBD's are branded
        for va_classid in va_classids:
            node = va_node_d[va_classid]
            if node.data['children'] == 0:  # leaf VA classid code
                rxnav.get_generic_drugs_for_VA_class(va_classid)
        # end for va_classid loop
        return scd_rxcuis_for_va_classid_d
    # end cache_generic_drugs_for_VA_classes

    # phase3_task <start>:
    logfile = io.open(log_filename, 'w', encoding='utf-8')
    print('[%s] Starting' % get_timestamp_string(), file=logfile);
    logfile.flush()
    rxnav = rxnav_rest_api_mp(cache_filename, utility_fns, logfile, mp_queue,
                              readonly_access_to_cache=True,
                              forward_result_to_cache_writer=True)
    # get VA drug classes
    print('Determine VA drug class hierarchy', file=logfile)
    VA_classId = 'VA000'  # Aug 6, 2018 change, per Lee Peters, root code for VA classes
    d = rxnav.get_class_tree(VA_classId)  # obtain VA hierarchy tree from NLM, VA hierarchy root 'VA000'
    # Determine VA classid set
    next_id = [1]
    t = BasicTree()
    va_classid_set = set()
    va_classid_d = {}
    va_node_d = {}
    build_va_folders(t, d['rxclassTree'])
    # Find generic drugs for VA class information and load into cache
    cache_generic_drugs_for_VA_classes(rxnav, va_classid_set, va_node_d)

# end phase3_task

# -------------------- Part 1.  Building Cache file using REST API and Multiprocessing ------------

def build_cache_file_using_multiprocessing(opts):

    # Local methods.
    # Implicitly depend on variables defined at the outer layer
    #  -- opts, logfile, rxnav, mp_queue, utility_fns

    # PHASE 1 -- determine 'allrelated' and 'historicalrxcui' values
    def phase1__get_allrelated_plus_historicalrxcui(rxcui_list):

        # Start workers to process RxCUI code segments
        # Wait for all to reach barrier ==> they have read initial contents of cache fle
        worker_process_count = opts.workers  # defaults to 4
        print('[%s] Creating %d worker processes'
              % (get_timestamp_string(), worker_process_count), file=logfile)
        logfile.flush()
        barrier = mp.Barrier(worker_process_count + 1)  # workers and manager wait at barrier
        rxcui_segsize = math.ceil(len(rxcui_list) / float(worker_process_count))  # RxCUIs per worker
        worker_processes = []
        for idx in range(worker_process_count):
            worker_number = idx + 1
            worker_process = mp.Process(target=worker_task_rxcuis,
                                        args=(mp_queue,
                                              worker_number,
                                              opts.cache,
                                              barrier,
                                              rxcui_list[(rxcui_segsize * idx):(rxcui_segsize * (idx + 1))],
                                              opts.log_dir + ('rxcui_worker_%d.log' % worker_number),
                                              utility_fns))
            worker_process.start()
            worker_processes.append(worker_process)
        # end worker process creation loop

        # wait for workers to initialize, wait at same barrier that workers wait at.
        # At this point, all have read initial state of cache.  Go ahead and start adding to cache.
        print('[%s] Waiting at barrier 1' % get_timestamp_string(), file=logfile);
        logfile.flush()
        barrier.wait()
        print('[%s] Passing barrier 1' % get_timestamp_string(), file=logfile);
        logfile.flush()

        # wait for workers to finish
        for worker in worker_processes:
            worker.join()  # wait for termination
        # workers have completed phase 1

        print('[%s] Done with Phase 1' % get_timestamp_string(), file=logfile);
        logfile.flush()
        return  # Done

    # end phase1

    # PHASE 2 -- determine NDC codes for drugs
    def phase2__get_ndc_for_drugs(drug_rxcui_list):

        # Start workers to determine NDC codes associated with RxCUI drug codes
        worker_process_count = opts.workers  # defaults to 4
        print('[%s] Creating %d NDC worker processes'
              % (get_timestamp_string(), worker_process_count), file=logfile)
        logfile.flush()
        barrier = mp.Barrier(worker_process_count + 1)  # workers and manager wait at barrier
        ndc_segsize = math.ceil(len(drug_rxcui_list) / float(worker_process_count))  # drug RxCUI count per worker
        worker_processes = []
        for idx in range(worker_process_count):
            worker_number = idx + 1
            worker_process = mp.Process(target=worker_task_ndcs,
                                        args=(mp_queue,
                                              worker_number,
                                              opts.cache,
                                              barrier,
                                              drug_rxcui_list[(ndc_segsize * idx):(ndc_segsize * (idx + 1))],
                                              opts.log_dir + ('ndc_worker_%d.log' % worker_number),
                                              utility_fns))
            worker_process.start()
            worker_processes.append(worker_process)
        # end worker process creation loop

        # wait for workers to initialize, wait at same barrier that workers wait at.
        # At this point, all have read initial state of cache.  Go ahead and start adding to cache.
        print('[%s] Waiting at barrier 2' % get_timestamp_string(), file=logfile);
        logfile.flush()
        barrier.wait()
        print('[%s] Passing barrier 2' % get_timestamp_string(), file=logfile);
        logfile.flush()

        # wait for workers to finish
        for worker in worker_processes:
            worker.join()  # wait for termination
        # workers have completed phase 2

        print('[%s] Done with Phase 2' % get_timestamp_string(), file=logfile);
        logfile.flush()
        return  # Done

    # end phase2

    # PHASE0 -- determine ACTIVE/RETIRED/etc status for each RxCUI
    def phase0():
        worker_process = mp.Process(target=phase0_task,
                                    args=(mp_queue,
                                          opts.log_dir + 'phase0.log',
                                          opts.cache,
                                          utility_fns,))
        worker_process.start()
        worker_process.join()  # wait for termination

    # end phase0

    # PHASE3 -- anything else needed in cache for i2b2 metadata
    def phase3():
        worker_process = mp.Process(target=phase3_task,
                                    args=(mp_queue,
                                          opts.log_dir + 'phase3.log',
                                          opts.cache,
                                          utility_fns,))
        worker_process.start()
        worker_process.join()  # wait for termination
    # end phase3

    def start_cache_writer():  # start Cache Writer
        cache_writer_process = mp.Process(target=cache_writer_task,
                                          args=(mp_queue,
                                                opts.cache,
                                                opts.log_dir + 'cache_writer.log',
                                                utility_fns))
        cache_writer_process.start()
        return  # Done

    # end start_cache_writer

    def stop_cache_writer():  # stop Cache Writer
        mp_queue.put('kill')
        return  # Done

    # end stop_cache_writer

    def create_rxnav_object():
        rxnav = rxnav_rest_api_mp(opts.cache, utility_fns, logfile, mp_queue,
                                  readonly_access_to_cache=True,
                                  forward_result_to_cache_writer=False)
        return rxnav  # Done

    # end create_rxnav_object

    def get_rxcuis_and_rxcuis_status_d():
        rxcui_set, rxcuis_status_d = \
            rxnav.get_historical_rxcuis(target_status_values=["ACTIVE", "RETIRED", "NEVER%20ACTIVE"],
                                        verbose=True)  # verbose
        return rxcui_set, list(rxcui_set), rxcuis_status_d  # Done

    # end get_rxcuis_and_rxcuis_status_d

    def determine_drug_rxcui_set(rxnav, rxcuis, rxcuis_status_d):
        drug_rxcui_set = set()
        target_attributes = ['TTY']
        target_attributes_d = {nm: idx for idx, nm in enumerate(target_attributes)}
        for idx, rxcui in enumerate(rxcuis):
            if rxcuis_status_d[rxcui] == 'NON-RXNORM': continue  # not of interest, no tty, etc
            ''' It is now established that this is an RxNorm code '''
            rxcui_attributes = rxnav.get_historical_rxcui_attributes(rxcui, target_attributes)  # tuple
            tty = rxcui_attributes[target_attributes_d['TTY']]
            if tty in ['SCD', 'SBD', 'GPCK', 'BPCK']:
                drug_rxcui_set.add(rxcui)
        return drug_rxcui_set  # Done

    # end determine_drug_rxcui_set

    def determine_rxcui_set():
        # Determine historically comprehensive set of RxNorm codes from NLM (rxcuihistory api)
        # NOTE: attributes not returned, only a set of codes
        # NOTE: skipping "NON-RXNORM" codes, which accounts for over 300,000 codes as of May 2018,
        #       which are all unrelated to the ingredient and drug RXCUI codes
        rxcui_set, rxcuis_status_d = \
            rxnav.get_historical_rxcuis(target_status_values=["ACTIVE", "RETIRED", "NEVER%20ACTIVE"], verbose=True)
        nonrxnorm_rxcui_set, nonrxnorm_status_d = \
            rxnav.get_historical_rxcuis(target_status_values=["NON-RXNORM"], verbose=True)
        # Sanity check the returned values -- rxcuis and nonrxnorm_rxcuis should be mutually exclusive
        overlap_of_rxnorm_nonrxnorm = rxcui_set & nonrxnorm_rxcui_set
        if len(overlap_of_rxnorm_nonrxnorm) > 0:
            print('*** [Sanity checking] failure: found %d NON-RXNORM RXCUIs which overlapped with RXNORM RXCUIS'
                  % len(overlap_of_rxnorm_nonrxnorm), file=logfile)
            print(str(list(overlap_of_rxnorm_nonrxnorm)), file=logfile)
            print('[End of list]', file=logfile)
            nonrxnorm_rxcui_set = nonrxnorm_rxcui_set - rxcui_set
            print('Size of non-overlapping set of NON-RXNORM RXCUIS is %d' % len(nonrxnorm_rxcui_set), file=logfile)
        else:
            print('[Sanity checking] Passed: No overlap of rxcuis and nonrxnorm_rxcuis, as expected.',
                  file=logfile)
        logfile.flush()
        return rxcui_set  # Done

    # end determine_rxcui_set

    # build_cache_file_using_multiprocessing:
    logfile = io.open(opts.log_dir + 'manager.log', 'w', encoding='utf-8')
    print('Starting', file=logfile);
    logfile.flush()
    f = io.open(opts.cache, 'a+', encoding='utf-8')  # create cache file if it does not exist, else open append
    f.close()  # the only point was to ensure the file exists
    manager = mp.Manager()  # multiprocessing initialization
    mp_queue = manager.Queue()  # create shared Queue, used by manager, workers, Cache Writer
    utility_fns = utility_functions()  # needed by rxnav interface -- e.g. flatten fn
    rxnav = create_rxnav_object()  # Initialize with NLM's REST API interface class
    rxcui_set = determine_rxcui_set()  # Determine all RxNORM RxCUI codes, excludes NON-RXNORM
    rxcui_list = sorted(list(rxcui_set))  # values are distributed to workers in sorted order, easy to follow

    # MAIN processing logic:
    start_cache_writer()  # start the Cache Writer
    phase0() # get per RxCUI status -- ACTIVE/RETIRED/etc
    phase1__get_allrelated_plus_historicalrxcui(rxcui_list)  # get 'allrelated' and 'historicalrxcui' results
    # Get a new rxnav object, which will re-read cache file, now with 'historicalrxcuis' and 'allrelated' elements.
    print('[%s] Creating new rxnav object, read the updated cache file' % get_timestamp_string(), file=logfile)
    logfile.flush()
    rxnav = create_rxnav_object()  # generate a new rxnav object, re-reads the updated cahe file
    rxcui_set, rxcui_list, rxcuis_status_d = get_rxcuis_and_rxcuis_status_d()
    drug_rxcui_set = determine_drug_rxcui_set(rxnav, rxcui_list, rxcuis_status_d)
    drug_rxcui_list = list(drug_rxcui_set)
    phase2__get_ndc_for_drugs(drug_rxcui_list)  # get NDC codes associated with drug RxCUIs
    phase3() # anything else needed in cache for i2b2 metadata
    stop_cache_writer()  # stop the Cache Writer
    print('[%s] Terminating' % get_timestamp_string(), file=logfile);
    logfile.flush()

# end build_cache_file_using_multiprocessing


# ------------------- Part 2.  Build i2b2 metadata ------------------------------

# ---- \PROVENANCE\ folder

class provenance_metadata_builder():
    ''' create PROVENANCE metadata folder and children '''

    def __init__(self):
        self.start_t = None # set by build method
        self.metadata_rows = {}
    # end (constructor)

    def build_provenance_folders(self, level, parent_path, rootpath):
        self.rootpath = rootpath
        # Create PROVENANCE folder
        mypath = parent_path + 'PROVENANCE' + '\\'
        c_basecode = ''
        c_fullname = mypath
        c_hlevel = level
        c_name = 'Provenance'
        c_tooltip = 'metadata provenance'
        self.metadata_rows[mypath] = \
            {'c_fullname': c_fullname,
             'c_hlevel': c_hlevel,
             'c_name': c_name,
             'c_basecode': c_basecode,
             'c_visualattributes': 'FH ', # specifically control this
             'c_tooltip': c_tooltip,
             'children': 0,
             'rxcui': None}
        # Create children of PROVENANCE folder -- VERSION, SOURCE, AUTHOR
        parent_path = c_fullname
        parent_hlevel = c_hlevel
        for child_tup in [('VERSION','RXNORM_20181101'),
                          ('SOURCE', 'NLM'),
                          ('AUTHOR', 'jay.pedersen@unmc.edu'),
                          ('BUILD_DATE', time.strftime('%Y-%m-%d'))]:
            tagname, tooltip_text = child_tup
            mypath = parent_path + tagname + '\\'
            c_basecode = ''
            c_fullname = mypath
            c_hlevel = parent_hlevel + 1
            c_name = tagname
            c_tooltip = tooltip_text
            self.metadata_rows[parent_path]['children'] += 1
            self.metadata_rows[mypath] = \
                {'c_fullname': c_fullname,
                 'c_hlevel': c_hlevel,
                 'c_name': c_name,
                 'c_basecode': c_basecode,
                 'c_visualattributes': 'LH ',
                 'c_tooltip': c_tooltip,
                 'children': 0,
                 'rxcui': None}
        # end for child_tup
    # end build_provenance_folders

    def build(self, path_prefix, metadata_root_level):
        self.start_t = time.time()

        # (Levels 1 to 5) Determine VA drug hierarchy
        print('Adding PROVENANCE metadata folder')
        root_path = '\\'+path_prefix+'\\'
        self.build_provenance_folders(metadata_root_level + 1, root_path, root_path)
        return self.metadata_rows # Done
    # end build

# end class provenance_metadata_builder


# ----- Parent folder for MEDICATION and VA_CLASS folders

class pcori_meds_first_level_metadata_builder():
    ''' Needed for creating a single metadata row, the root level row '''

    def build(self, path_prefix, prefix_level):
        # create ATC drug class hierarchy for RxNorm
        path = '\\'+path_prefix+'\\'
        return { path : { 'c_fullname': path, 'c_hlevel': prefix_level, 'c_name': 'Medications',
                          'c_basecode': 'RXNORM_ROOT', 'children': 2, 'c_visualattributes': 'CA ' } }

# end class pcori_meds_first_level_metadata_builder

# ------------------            Cache Writer            ----------------------------


def metadata_writer_task(opts): # this is NOT a class, even though large, no __init__
    ''' Metadata Writer -- reads from NLM REST API cache and creates i2b2 metadata from that information '''

    from rxnorm_code_maps import rxnorm_code_maps
    from ingredient_metadata_builder import ingredient_metadata_builder
    from modifier_metadata_builder import modifier_metadata_builder
    from va_metadata_builder import va_metadata_builder
    from metadata_writer import metadata_writer

    # Local methods.
    # Implicitly depend on variables defined at the outer layer
    #  -- opts

    def determine_ingredient_rxcui_set(rxnav, rxcuis, rxcuis_status_d):
        ingredient_rxcui_set = set()
        target_attributes = ['TTY']
        target_attributes_d = {nm: idx for idx, nm in enumerate(target_attributes)}
        for idx, rxcui in enumerate(rxcuis):
            if rxcuis_status_d[rxcui] == 'NON-RXNORM': continue  # not of interest, no tty, etc
            ''' It is now established that this is an RxNorm code '''
            rxcui_attributes = rxnav.get_historical_rxcui_attributes(rxcui, target_attributes)  # tuple
            tty = rxcui_attributes[target_attributes_d['TTY']]
            if tty in ['IN', 'MIN', 'PIN']:
                ingredient_rxcui_set.add(rxcui)
        return ingredient_rxcui_set

    # end determine_ingredient_rxcui_set

    def determine_drug_rxcui_set(rxnav, rxcuis, rxcuis_status_d):
        drug_rxcui_set = set()
        target_attributes = ['TTY']
        target_attributes_d = {nm: idx for idx, nm in enumerate(target_attributes)}
        for idx, rxcui in enumerate(rxcuis):
            if rxcuis_status_d[rxcui] == 'NON-RXNORM': continue  # not of interest, no tty, etc
            ''' It is now established that this is an RxNorm code '''
            rxcui_attributes = rxnav.get_historical_rxcui_attributes(rxcui, target_attributes)  # tuple
            tty = rxcui_attributes[target_attributes_d['TTY']]
            if tty in ['SCD', 'SBD', 'GPCK', 'BPCK']:
                drug_rxcui_set.add(rxcui)
        return drug_rxcui_set

    # end determine_drug_rxcui_set

    def post_process(va_metadata_rows, ingred_metadata_rows):
        # PostProcess VA metadata, find ingredients whose RXCUI sets differ in VA versus 'by ingredient'
        #     ==> add '(VA subset)' to c_name of ingredient (directly modify metadata row)
        #     NOTE: need both VA and by-ingredient metadata rows to perform this comparison ==> "post-processing"
        #     NOTE: use sorted lists rather than processing .keys(), faster

        # embedded methods

        def find_positions_of_subpaths(sorted_paths, path):
            startpos = sorted_paths.index(path) + 1
            listsize = len(sorted_paths)
            if startpos >= listsize or (not sorted_paths[startpos].startswith(path)):
                startpos, endpos = (listsize, listsize)
            else:  # we know at least one subpath exists
                endpos = startpos
                idx = startpos + 1  # we know [startpos] starts with path, move on
                while idx < listsize:
                    if sorted_paths[idx].startswith(path):
                        endpos += 1  # this element starts with path
                        idx += 1  # move on
                    else:
                        break  # doesnt start with path, not subpath -- done
            return (startpos, endpos)
        # end find_positions_of_subpaths

        # post_process:
        va_paths = sorted(va_metadata_rows.keys())  # ALL paths in VA classes folders (sorted) => VA+RXNORM+NDC codes
        by_paths = sorted(ingred_metadata_rows.keys())  # ALL paths in ingredients (alphabetized) folders (sorted)
        va_ingred_paths = {va_metadata_rows[x]['rxcui']: x for x in va_metadata_rows.keys()
                           if (va_metadata_rows[x]['tty'] in ['IN'])}
        by_ingred_paths = {ingred_metadata_rows[x]['rxcui']: x for x in ingred_metadata_rows.keys()
                           if (ingred_metadata_rows[x]['tty'] in ['IN'])}
        va_name_change_count = 0
        for ingred_rxcui in sorted(va_ingred_paths.keys()):
            va_path = va_ingred_paths[ingred_rxcui]
            va_paths_startpos, va_paths_endpos = find_positions_of_subpaths(va_paths, va_path)
            va_rxnorm_codes = set(
                va_metadata_rows[x]['c_basecode'] for x in va_paths[va_paths_startpos:va_paths_endpos + 1]
                if (va_metadata_rows[x]['c_basecode'].startswith('RXNORM:')))
            by_path = by_ingred_paths[ingred_rxcui]
            by_paths_startpos, by_paths_endpos = find_positions_of_subpaths(by_paths, by_path)
            by_rxnorm_codes = set(
                ingred_metadata_rows[x]['c_basecode'] for x in by_paths[by_paths_startpos:by_paths_endpos + 1]
                if (ingred_metadata_rows[x]['c_basecode'].startswith('RXNORM:')))
            if va_rxnorm_codes != by_rxnorm_codes:  # see if difference in RxCUI code sets
                va_metadata_rows[va_path]['c_name'] += ' (VA subset)'
                va_name_change_count += 1
        # end for ingred_rxcui
        print('[[NOTE: changed %d VA metadata rows to (VA subset), differences with by-ingredient RXCUI sets]]' %
              va_name_change_count)

    # end post_process

    # ------ metadata_writer_task <start>: ---------
    logfile = io.open(opts.log_dir + 'metadata_writer.log', 'w', encoding='utf-8')
    # create global objects which interface to NLM, etc
    utility_fns = utility_functions()
    rxnav = rxnav_rest_api_mp(opts.cache, utility_fns, logfile, None,
                              readonly_access_to_cache=True,
                              forward_result_to_cache_writer=False,
                              fail_if_not_in_cache=True)

    # Determine historically comprehensive set of RxNorm codes from NLM (rxcuihistory api)
    # NOTE: attributes not returned, only a set of codes
    # NOTE: skipping "NON-RXNORM" codes, which accounts for over 300,000 codes as of May 2018,
    #       which are all unrelated to the ingredient and drug RXCUI codes
    rxcui_set, rxcuis_status_d = \
        rxnav.get_historical_rxcuis(target_status_values=["ACTIVE", "RETIRED", "NEVER%20ACTIVE"],
                                    verbose=True)  # verbose
    nonrxnorm_rxcui_set, nonrxnorm_status_d = \
        rxnav.get_historical_rxcuis(target_status_values=["NON-RXNORM"],
                                    verbose=True)  # verbose
    # Sanity check the returned values -- rxcuis and nonrxnorm_rxcuis should be mutually exclusive
    overlap_of_rxnorm_nonrxnorm = rxcui_set & nonrxnorm_rxcui_set
    if len(overlap_of_rxnorm_nonrxnorm) > 0:
        print('*** [Sanity check] failure: found %d NON-RXNORM RXCUIs which overlapped with RXNORM RXCUIS'
              % len(overlap_of_rxnorm_nonrxnorm), file=logfile)
        print(str(list(overlap_of_rxnorm_nonrxnorm)), file=logfile)
        print('[End of list]', file=logfile)
        nonrxnorm_rxcui_set = nonrxnorm_rxcui_set - rxcui_set
        print('Size of non-overlapping set of NON-RXNORM RXCUIS is %d' % len(nonrxnorm_rxcui_set), file=logfile)
    else:
        print('[Sanity check] Passed: No overlap of rxcuis and nonrxnorm_rxcuis, as expected.', file=logfile)
    # Obtain historical information for each RXCUI, needed to build RXNORM.csv ==>
    #    name, type, dates, establish category INGREDIENT, DRUG, OTHER
    rxcui_count = len(rxcui_set)
    print('Extracting attributes definitions for the %d historical RxCUI codes from the NLM' % rxcui_count,
          file=logfile)
    logfile.flush()
    rxcui_d = {}
    for idx, rxcui in enumerate(rxcui_set):  # integer rxcui value
        if (idx + 1) % 10000 == 0:  # Track progress
            print('[%s] Obtaining definition for rxcui %d of %d, %s%%' \
                  % (time.strftime("%Y-%m-%d %H:%M:%S"),
                     idx + 1, rxcui_count, str(round(((idx + 1) / rxcui_count) * 100, 3))),
                  file=logfile)
            logfile.flush()
        d = rxnav.get_historical_rxcui(rxcui)
        rxcui_d[rxcui] = d  # save NLM's definition of RxCUI, from rxcuihistory api
    print('[Done extracting RxCUI attribute definitions, found: %d]' % len(rxcui_d), file=logfile)
    logfile.flush()

    ''' Establish set of ingredient and drug codes '''
    ingredient_rxcui_set = determine_ingredient_rxcui_set(rxnav, rxcui_set, rxcuis_status_d)
    drug_rxcui_set = determine_drug_rxcui_set(rxnav, rxcui_set, rxcuis_status_d)
    ''' build maps of ingredients and drugs '''
    rxnorm_coding = rxnorm_code_maps(rxnav, ingredient_rxcui_set, drug_rxcui_set)
    ingredient_to_drug_set_map = rxnorm_coding.get_ingredient_to_drug_set_map()

    # build metadata - information is gathered
    path_prefix = 'i2b2_RXNORM_NDC' if not opts.prefix else opts.prefix
    # NOTE: path_prefix default from JRC's metadata interoperability guidelines, Nov 2018
    #       previously 'PCORI\\MEDICATION', which was Harvard SCILHS
    print('Metadata path prefix [%s]' % path_prefix)
    prefix_level = opts.prefix_level
    if opts.add_provenance: # create PROVENANCE portion of metadata if requested
        provenance_metadata_rows = provenance_metadata_builder().build(path_prefix, prefix_level)
    # Add MODIFIERS
    if not opts.no_modifiers:
        modifier_rows = modifier_metadata_builder().build(path_prefix, prefix_level, 'interop_MED_MODIFIERS.txt')
    # -------- va_meds_row -- \<prefix>\MEDICATION, row above VA and BY-INGREDIENT
    if not opts.only_by_ingredient:
        va_meds_row_builder = pcori_meds_first_level_metadata_builder()
        va_meds_row = va_meds_row_builder.build(path_prefix, prefix_level)
    # ------- VA/NDF-RT drug classification metadata build ---------------
    if not opts.only_by_ingredient:
        va_metadata_builder = va_metadata_builder(rxnav, rxnorm_coding)
        va_metadata_rows = va_metadata_builder.build(path_prefix, prefix_level + 1)  # lowest level is prefix level + 1
    # ------- BY-INGREDIENT metadata building (alphabetic/lexicographic order) build ---------
    ingred_metadata_builder = ingredient_metadata_builder(rxnorm_coding, rxnav)
    ingred_metadata_rows = ingred_metadata_builder.build(path_prefix, prefix_level)
    # -------- POST-PROCESS -----------------------------
    if not opts.only_by_ingredient:
        post_process(va_metadata_rows, ingred_metadata_rows)  # VA metadata update
    # ---- write metadata ----
    metadata_fn = opts.output  # 'scilhs_med.txt'
    print('[[ %s metadata file [%s] ]]' % (('Creating' if not opts.append else 'Appending to'), metadata_fn))
    metadata_writer = metadata_writer(metadata_fn,
                                      append=opts.append,
                                      include_dates=(not opts.dont_include_dates),
                                      include_tty=opts.include_tty,
                                      csv_fields_file=opts.csv_fields_file,
                                      csv_fields_file_delim=opts.csv_fields_file_delim,
                                      strftime_format=opts.strftime_format)
    if opts.add_provenance: # PROVENANCE folder
        metadata_writer.write_metadata_rows(provenance_metadata_rows)
    metadata_writer.write_metadata_rows(modifier_rows) # MODIFIER metadata
    if not opts.only_by_ingredient:
        metadata_writer.write_metadata_rows(va_meds_row)
        metadata_writer.write_metadata_rows(va_metadata_rows)
    metadata_writer.write_metadata_rows(ingred_metadata_rows)  # always write by-ingredient meds

    print('[%s] Finished writing i2b2 metadata ... terminating.'
          % (get_timestamp_string(),), file=logfile)
    logfile.flush()
    logfile.close()

# end metadata_writer_task


def build_i2b2_metadata_from_cache_file(opts):

    # build_i2b2_metadata
    #
    # Algorithm:
    #    Use a subprocess to perform the metadata file creation, which will resolve
    #    memory usage issues that would otherwise occur within a single process
    #    which is processing the cache file multiple times and creating massive
    #    dictionaries.

    worker_process = mp.Process(target=metadata_writer_task,
                                args=(opts,))
    worker_process.start()
    worker_process.join() # wait for termination

# end build_i2b2_metadata_from_cache_file


# -------------------- MAIN -------------------


def main():  # need separate main() for use of subprocess and task definitions in same file

    # Local methods to main(), to separate the logic into manageable pieces
    # They implicitly depend on varaibles created at the begining of main()
    #  -- opts, logfile, rxnav, mp_queue, utility_fns

    def parse_command():
        import optparse
        opt = optparse.OptionParser()
        opt.add_option('--cache', action='store', default='rxcui.cache')
        opt.add_option('--output_dir', action='store', default='./')
        opt.add_option('--output_filename', action='store', default='rxnorm_ndc.txt')
        opt.add_option('--append', action='store_true')
        opt.add_option('--log_dir', action='store', default='./')
        opt.add_option('--rxcui_relationships_csv', action='store')
        opt.add_option('--workers', action='store', type=int, default=4)
        opt.add_option('--verbose', action='store_true')
        opt.add_option('--logfile', action='store', default='-')
        opt.add_option('--only_by_ingredient', action='store_true')
        opt.add_option('--prefix', action='store')
        opt.add_option('--prefix_level', action='store', type=int, default=1)
        opt.add_option('--no_modifiers', action='store_true')
        opt.add_option('--dont_include_dates', action='store_true')
        opt.add_option('--include_tty', action='store_true')
        opt.add_option('--csv_fields_file', action='store')
        opt.add_option('--csv_fields_file_delim', action='store', default='|')
        opt.add_option('--strftime_format', action='store', default='%Y%m%d')
        opt.add_option('--add_provenance', action='store_true')
        opts, args = opt.parse_args()
        for folder_name in [opts.output_dir, opts.log_dir]:
            if folder_name[-1] != '/': folder_name += '/'  # folder spec, so we can append filename always
        opts.output = opts.output_dir + opts.output_filename # create folder+filename
        if opts.verbose:
            print('opts are %s' % str(opts))
        return opts, args  # Done

    # START OF MAIN

    # MAIN initialization:
    opts, args = parse_command()
    logfile = io.open(opts.log_dir + 'build.log', 'w', encoding='utf-8')
    print('[%s] Step 1 ==> Creating cache of needed REST API call values' % get_timestamp_string(), file=logfile)
    logfile.flush()
    build_cache_file_using_multiprocessing(opts) # historicalrxcui and getall related
    print('[%s] Step 2 ==> Creating i2b2 metadata from the cached information' % get_timestamp_string(), file=logfile)
    logfile.flush()
    build_i2b2_metadata_from_cache_file(opts)
    print('[%s] Completion.  Metadata created ... terminating.' % get_timestamp_string(), file=logfile)
    logfile.flush()
    print('Finished creating metadata.')

# end main()


# MAIN stub
if __name__ == '__main__':  # needed do to multiprocessing lib requirements
    main()
    sys.exit(0)
