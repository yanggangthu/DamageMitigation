'''
This module contains a two-byte repair function, a generalization of
the implementation of the  network repair methodology
introduced in Campbell and Albert (2014), BMC Syst. Biol.


Gang Yang
Contact: gzy105@psu.edu
Python Version: 2.7.x
Date: Aug 2014
'''


import networkx as nx
import numpy as np
from itertools import product
from random import choice, randrange, sample



def form_network(rules):
    '''
    Takes as input a list of rules in the format of sample_network.txt.

    Outputs a networkx DiGraph with node properties:
        'update_nodes': a list of regulating nodes
        'update_rules': a dictionary with binary strings as keys, corresponding
                        to the possible states of the update_nodes, and integers
                        as values, corresponding to the state of the node given
                        that input.

    Note that nodes are identified by their position in a sorted list of
    user-provided node names.
    '''
    def clean_states(x):
        #cleans binary representation of node input states
        out=x[2:]                                                               # Strip leading 0b
        return '0'*(len(inf)-len(out))+out                                      # Append leading 0's as needed

    stream = [x.rstrip('\n') for x in rules if x != '\n' and x[0]!='#']         # Remove comments and blank lines
    nodes = sorted([x.split(' ',1)[0][:-1] for x in stream])                    # Generate a sorted list of node names

    g = nx.DiGraph()
    g.graph['knockout'] = None                                                  # At creation, no node is flagged for knockout or overexpression
    g.graph['express'] = None

    for n in xrange(len(stream)):
        rule = stream[n].split('* = ')[1]
        rule = rule.replace(' AND ',' and ')                                    # Force decap of logical operators so as to work with eval()
        rule = rule.replace(' OR ',' or ')
        rule = rule.replace(' NOT ',' not ')
        if stream[n].find('True') >= 0 or stream[n].find('False') >= 0:         # For always ON or always OFF nodes
            g.add_node(n)                                                       # We refer to nodes by their location in a sorted list of the user-provided node names
            g.node[n]['update_nodes'] = []
            g.node[n]['update_rules'] = {'':str(int(eval(rule)))}
            continue

        inf = rule.split(' ')                                                   # Strip down to just a list of influencing nodes
        inf = [x.lstrip('(') for x in inf]
        inf = [x.rstrip(')') for x in inf]
        #The sort ensures that when we do text replacement (<node string>->'True' or 'False') below in this fn, we avoid problems where node 1 is a substring of node 2 (e.g. NODE1_phosph and NODE1)
        inf = sorted([x for x in inf if x not in ['','and','or','not']],key=len,reverse=True)

        for i in inf: g.add_edge(nodes.index(i),n)                              # Add edges from all influencing nodes to target node
        g.node[n]['update_nodes'] = [nodes.index(i) for i in inf]
        g.node[n]['update_rules'] = {}

        bool_states = map(bin,range(2**len(inf)))
        bool_states = map(clean_states,bool_states)
        for j in bool_states:
            rule_mod = rule[:]
            for k in range(len(j)):
                if j[k]=='0':
                    rule_mod=rule_mod.replace(nodes[g.node[n]['update_nodes'][k]],'False')      # Modify the rule to have every combination of True, False for input nodes
                else: rule_mod=rule_mod.replace(nodes[g.node[n]['update_nodes'][k]],'True')
            g.node[n]['update_rules'][j] = int(eval(rule_mod))                                  # Store outcome for every possible input

    return g,nodes

def find_attractor(graph,state=False):
    '''
    Takes a graph as formatted in form_network() as input.

    Chooses a random starting state (if state==False); synchronously advances
    until a SS or LC is found. Can accommodate a node knockout or overexpression
    as formed in damage_network().

    if state is not False, it must be a string of 0s and 1s, which specifies the
    starting state for the update iteration.

    The string bits must be are arranged in alphabetical (sorted) order,
    according to node names.

    Outputs a list of [next_state,attractor]. For a SS, the 'attractor' entries
    are '0's and '1's, representing the SS. For a LC, the 'attractor' is a list
    of state strings, representing every state in the LC.

    If no state is found after 1000 interations, the function returns False and
    prints a warning to the console.
    '''
    def update_state(x):
        #x is the node's index
        input_state = ''
        for i in graph.node[nodes[x]]['update_nodes']: input_state += str(trajectory[-1][nodes.index(i)])    # Acquire string of current states of node x's input nodes
        if input_state=='':
            return str(trajectory[-1][x])                                                                    #Gang changes the code here in case that this is a source node, input_state gonna be '' 07/23/2014
        else:
            return str(graph.node[nodes[x]]['update_rules'][input_state])

    nodes = sorted(graph.nodes())
    if not state: trajectory = [list(np.random.random_integers(0,1,nx.number_of_nodes(graph)))]         # Random starting state
    else: trajectory = [state]                                                                          # Provided starting state

    while True:
        trajectory += [map(update_state,xrange(len(nodes)))]
        if graph.graph['knockout'] != None:
            for node in graph.graph['knockout']:
                trajectory[-1][nodes.index(node)] = '0'  # If a node has been knocked out, it must be 0 even if it would normally be active
        elif graph.graph['express'] != None:
            for node in graph.graph['express']:
                trajectory[-1][nodes.index(node)] = '1'  # "  " "    "   "    expressed,   "  "    "  1 "    "  "  "     "        "  inactive
        if trajectory[-1] in trajectory[:-1]:                                                           # Return a list of [next state, attractor], once attractor is found (attractor list length 1 in case of SS)
            return [trajectory[1],trajectory[trajectory.index(trajectory[-1]):-1]]
        if len(trajectory) == 10000:
            print 'find_attractor() was unable to find an attractor in 1000 interations; returning False.'
            return False

def superset(a):
    '''
    Takes a limit cycle (list of lists of binary strings), and collapses it into
    one state with a 1 if the node is 1 in any states in 'a', and 0 otherwise.
    '''
    transpose = [[row[i] for row in a[1]] for i in range(len(a[1][0]))]         # Turn list of states into list of time sequences for each node [[n1s1,n2s1],[n1s2,n2s2]] -> [[n1s1,n1s2],[n2s1,n2s2]]
    superset = ['1' if '1' in x else '0' for x in transpose]                    # Evaluate each entry for any 1s; generate superset
    return superset


def damage_network(graph,a=None,force=None,force_type_list=['knockout','knockout'],num_damage=2):
    '''
    Takes a graph and superset attractor as input.
    Two nodes will be forced into another state, which can be knocted out or expressed.
    force=None leads either a transiently
    active or permanently inactive node to be knocked out or
    overexpressed, respectively.
    num_damage is the number of modifications in the network, which will be used to generate
    random positions in the network. Set default as 2.

    Alternatively, if force != False, it must be a list of integers of a node's position
    in a sorted list of graph nodes. That node is then set to 'knockout' or
    'express' according to force_type_list.
    force_type_list contains a list of force_type, corresponding to the list of force.

    Returns the modified graph.
    '''
    g_copy = graph.copy()                                                       # Don't modify the input graph
    g_copy.graph['knockout'] = []                                               # Initialization of these two properties of the graph
    g_copy.graph['express'] = []                                                # Notice the two properties exist as a list

    if force == None:
        nodelist = sample(range(len(a)),num_damage)                             # Generate a list of random positions to be modified, with a length as num_damage
        for node in nodelist:                                                   # Treat each node individually
            if a[node] == 0: g_copy.graph['express'].append(node)                         # Assign the index of an overexpressed (0 forced into 1) node
            else: g_copy.graph['knockout'].append(node)                                   # "      "   "     "  a knocked out (1 forced into 0)    "
    else:
        for i in range(len(force_type_list)):
            if force_type_list[i] == 'knockout': g_copy.graph['knockout'].append(force[i])
            elif force_type_list[i] == 'express': g_copy.graph['express'].append(force[i])

    return g_copy

def damage_state(graph,a):
    '''
    Form a damaged state (in superset form or regular form, depending on input)
    that mirrors the original attractor except for the knockout/damage.

    Called internally in other functions; likely not of direct interest to the
    end user.
    '''
    if graph.graph['knockout'] != None:
        for node in graph.graph['knockout']:
            a = a[:node]+['0']+a[node+1:]
    elif graph.graph['express'] != None:
        for node in graph.graph['express']:
            a = a[:node]+['1']+a[node+1:]
    else:
        pass

    return a



def compare_attractors(graph,a):
    '''
    Takes as input a damaged graph and the original, undamaged LC attractor.

    Returns (as the second value) a list where the first entry is the first
    entry of 'a' (preserving format) and the second is the largest component of
    the attractor that survives when duplicate states are collapsed (e.g. 011
    and 010 when node 3 is fixed).

    Returns (as the first value) True of False, respectively corresponding to
    if a pair of state collapsed or not.
    '''
    if len(a[1]) == 1:
        return False, [a[0],[damage_state(graph,a[1][0])]]

    new_attractors = []
    valid,invalid = [],[]
    for state in a[1]:
        if graph.graph['knockout'] != None:
            for node in graph.graph['knockout']:
                if state[node] == '1':
                    switch_state = state[:node] + ['0'] + state[node+1:]
                    if switch_state in a[1]:
                        invalid += [state]
                        valid += [switch_state]                                     # Store the 'invalid' states; inputs to them reroute to the 'valid' states
        elif graph.graph['express'] != None:
            for node in graph.graph['express']:
                if state[node] == '0':
                    switch_state = state[:node] + ['1'] + state[node+1:]
                    if switch_state in a[1]:
                        invalid += [state]
                        valid += [switch_state]






    if len(valid)==0: return False, [a[0],[damage_state(graph,x) for x in a[1]]]# If there are no state collapses, we straightforwardly damage every state in a according to graph damage.

    positions = range(1,len(a[1]))
    cur_pos = 0
    route=[]  #route = [a[1][0]]
    while True:                                                                 # We walk across the state transition map, sensitive to invalid states due to state collapse
        if a[1][cur_pos] in invalid:
            while a[1][cur_pos] in invalid:
                try: positions.remove(cur_pos)                                  # We can sometimes loop back to the same state multiple times; no need to remove it from the "points to walk across" list again
                except ValueError: pass
                cur_pos = a[1].index(valid[invalid.index(a[1][cur_pos])])       # Jump from invalid position to corresponding valid position when moving onto invalid state
            route+=[a[1][cur_pos]]
        else:
            route+=[a[1][cur_pos]]
        falsified = [x for x in route if route.count(x) > 1]
        if len(falsified) > 0:                                                  # Check to see if we've hit a node twice (found a LC)
            start = route.index(route[-1])
            new_attractors += [route[start+1:]]                                 # The last entry is always one of the repeats, so the LC runs from the first instance of the repeat to the end

            try: positions.remove(cur_pos)
            except ValueError: pass                                             # We can sometimes loop back to the same state multiple times; no need to remove it from the "points to walk across" list again
            if len(positions) == 0: break

            cur_pos = choice(positions)                                         # If we isolate an attractor and other states have yet to be walked, pick one randomly
            route = []
            continue
        else:                                                                   # If we are still progressing, simply iterate our position forward
            try: positions.remove(cur_pos)
            except ValueError: pass
            cur_pos += 1
            if cur_pos == len(a[1]):                                            # We can also walk back to 0
                cur_pos = 0

    if len(new_attractors) == 0:
        raise RuntimeError("No new attractors identified.")
    elif len(new_attractors) == 1:
        return True, [a[0],[damage_state(graph,x) for x in new_attractors[0]]]  # We now have the largest salvageable section of the LC
    else:
        new_attractor_lengths = map(len,new_attractors)
        max_index = new_attractor_lengths.index(max(new_attractor_lengths))     # We choose to look at the longest "new" attractor
        return True, [a[0],[damage_state(graph,x) for x in new_attractors[max_index]]]  # We now have the largest salvageable section of the LC


def check_stability(graph,a):
    '''
    Takes as input the damaged graph and damaged attractor.

    If the attractor reached from every state in the damaged attractor is
    identical to the damaged attractor, the damaged attractor is stable, and
    the function returns True. Otherwise, the function returns False.
    '''
    def shift(seq, n):
        n = n % len(seq)
        return seq[n:] + seq[:n]

    A_n = [find_attractor(graph,state=x) for x in a[1]]
    if [type(x) for x in A_n].count(bool) > 0:
        raise RuntimeError("Unable to find attractor in check_stability()")
    if len(A_n) == len(a[1]) == 1:                                              # If the damaged attractor is a SS and evaluates to a SS, we can directly compare them
        if A_n[0][1] == a[1]: return True
        else: return False

    for attr in A_n:                                                            # If one or both is a LC, We check to make sure that the procession of states matches
        if len(attr[1]) != len(a[1]): return False                              # Quickest check is to see if the lengths are the same
        else:
            try: asi = a[1].index(attr[1][0])                                   # See the 'starting index' of A_d relative to this attractor
            except Exception: return False                                      # If it isn't in this attractor, they obviously don't match
            rotated = shift(a[1],asi)                                           # Then see if the aligned processions are equivalent at each position (same procession)
                                                                                # Gang has changed shift(a,asi) into shift(a[1],asi) 08/14/2014
            if [i==j for i,j in zip(rotated,attr[1])].count(False) > 0: return False
    return True



def evaluate_repair(graph,a,a_s=None,method='fix_to_SS'):
    '''
    Takes the damaged graph, damaged attractor, and original superset as input.
    method == 'LC_repair' or method == 'fix_to_SS', to determine whether or not
    we attempt to preserve LC transitions or fix its superset to a SS.

    Returns the repaired graph, or if not possible, a string explaining the
    cause of failure.
    '''
    # Define internal functions ------------------------------------------------
    def disjoint(x,y):
        '''
        Takes two lists and outputs the positions where the two don't have the
        same value.
        '''
        return [v for v in xrange(len(x)) if x[v]!=y[v]]

    def examine_modifications(g_in,a_t):
        '''
        Takes damaged network g and target attractor a_t. Looks at all possible
        edge modifications (as enumerated in the report) to force every node to
        be in its desired state.

        results are stored in g.graph['modifications'][<node>] as list of
        tuples:
        (<approach #>,<fix to 1/0>,<interacting node 1>,<interacting node 2>)

        '''
        def approach_1(x,pn,an,fix_to=1):
            '''
            adds OR <present_new> to an inactive node to make it active, or
            an AND <absent_new> to an active node to make it inactive.
            here and below:
                fix_to = 0 - we are fixing an active node to be inactive
                fix_to = 1 - we are fixing an inactive node to be active
                pn - viable_present_new_nodes
                p  - viable_present_nodes
                an - viable_absent_new_nodes
                a  - viable_absent_nodes
            '''
            if fix_to == 1:
                for j in pn:
                    g.graph['modifications'][x] += [(1,1,j)]                          #output meaning: method, fix_to, interacting node #1 (similar for other approach_# functions)
            else:
                for j in an:
                    g.graph['modifications'][x] += [(1,0,j)]
        def approach_2(x,pn,an,fix_to=1):
            '''
            adds an OR NOT <absent_new> to an inactive node to make it active,
            or an AND NOT <present_new> to an active node to make it inactive.
            '''
            if fix_to == 1:
                for j in an:
                    g.graph['modifications'][x] += [(2,1,j)]
            else:
                for j in pn:
                    g.graph['modifications'][x] += [(2,0,j)]
        def approach_3(x,p,pn,a,an,fix_to=1):
            '''
            adds an OR <pres and pres_new> to an inactive node to make it
            active, or an AND <abs or abs_new> to an active node to make it
            inactive.
            '''
            if fix_to == 1:
                for j,k in product(p,pn):
                    g.graph['modifications'][x] += [(3,1,j,k)]
            else:
                for j,k in product(a,an):
                    g.graph['modifications'][x] += [(3,0,j,k)]
        def approach_4(x,p,pn,a,an,fix_to=1):
            '''
            adds an OR <pres and not abs_new> to an inactive node to make it
            active, or an AND <abs OR NOT pres_new> to an active node to make it
            inactive.
            '''
            if fix_to == 1:
                for j,k in product(p,an):
                    g.graph['modifications'][x] += [(4,1,j,k)]
            else:
                for j,k in product(a,pn):
                    g.graph['modifications'][x] += [(4,0,j,k)]
        def approach_5(x,p,pn,a,an,fix_to=1):
            '''
            adds an OR <NOT abs AND pres_new> to an inactive node to make it
            active, or an AND <NOT pres OR abs_new> to an active node to make it
            inactive.
            '''
            if fix_to == 1:
                for j,k in product(a,pn):
                    g.graph['modifications'][x] += [(5,1,j,k)]
            else:
                for j,k in product(p,an):
                    g.graph['modifications'][x] += [(5,0,j,k)]
        def approach_6(x,p,pn,a,an,fix_to=1):
            '''
            adds an OR <NOT abs AND NOT abs_new> to an inactive node to make it
            active, or an AND <NOT pres OR NOT pres_new> to an active node to
            make it inactive.
            '''
            if fix_to == 1:
                for j,k in product(a,an):
                    g.graph['modifications'][x] += [(6,1,j,k)]
            else:
                for j,k in product(p,pn):
                    g.graph['modifications'][x] += [(6,0,j,k)]

        g=g_in.copy()
        nodes = sorted(g.nodes())
        g.graph['modifications'] = {}
        for i in g.nodes_iter():
            if (i in g.graph['express']) or (i in g.graph['knockout']): continue     # retest whether i exist as a number or a character                                                                                 # Don't attempt to modify the knocked out/overexpressed node
            g.graph['modifications'][i] = []                                                                                                                    # Set container for all viable modifications to this node (0->1 and 1->0)
            viable_present_new_nodes = [y for y in g.nodes_iter() if int(a_t[nodes.index(y)]) == 1 and y not in g.node[i]['update_nodes'] and y not in g.graph['express'] and y != i]     # Possible new regulators that are present nodes (that have not been overexpressed)
            viable_absent_new_nodes = [y for y in g.nodes_iter() if int(a_t[nodes.index(y)]) == 0 and y not in g.node[i]['update_nodes'] and y not in g.graph['knockout']  and y != i]     # Possible new regulators that are absent nodes (that have not been knocked out)
            viable_present_nodes = [y for y in g.node[i]['update_nodes'] if int(a_t[nodes.index(y)]) == 1 and y not in g.graph['express']  and y != i]                                     # Possible existing regulators that are present nodes (that have not been overexpressed)
            viable_absent_nodes = [y for y in g.node[i]['update_nodes'] if int(a_t[nodes.index(y)]) == 0 and y not in g.graph['knockout'] and y != i]                                     # Possible existing regulators that are absent nodes (that have not been overexpressed)
            if set(viable_present_new_nodes)&set(viable_absent_new_nodes)&set(viable_present_nodes)&set(viable_absent_nodes) != set(): raise RuntimeError("Algorithm halted - duplicate node assignment.")
            #go through all 6 approaches, store possible combinations of nodes as a graph property
            approach_1(i,viable_present_new_nodes,viable_absent_new_nodes,fix_to=int(a_t[nodes.index(i)]))
            approach_2(i,viable_present_new_nodes,viable_absent_new_nodes,fix_to=int(a_t[nodes.index(i)]))
            approach_3(i,viable_present_nodes,viable_present_new_nodes,viable_absent_nodes,viable_absent_new_nodes,fix_to=int(a_t[nodes.index(i)]))
            approach_4(i,viable_present_nodes,viable_present_new_nodes,viable_absent_nodes,viable_absent_new_nodes,fix_to=int(a_t[nodes.index(i)]))
            approach_5(i,viable_present_nodes,viable_present_new_nodes,viable_absent_nodes,viable_absent_new_nodes,fix_to=int(a_t[nodes.index(i)]))
            approach_6(i,viable_present_nodes,viable_present_new_nodes,viable_absent_nodes,viable_absent_new_nodes,fix_to=int(a_t[nodes.index(i)]))
            approach_1(i,viable_present_new_nodes,viable_absent_new_nodes,fix_to=(int(a_t[nodes.index(i)])+1)%2)                                                #look at both the 0->1 and 1->0 possibilities
            approach_2(i,viable_present_new_nodes,viable_absent_new_nodes,fix_to=(int(a_t[nodes.index(i)])+1)%2)
            approach_3(i,viable_present_nodes,viable_present_new_nodes,viable_absent_nodes,viable_absent_new_nodes,fix_to=(int(a_t[nodes.index(i)])+1)%2)
            approach_4(i,viable_present_nodes,viable_present_new_nodes,viable_absent_nodes,viable_absent_new_nodes,fix_to=(int(a_t[nodes.index(i)])+1)%2)
            approach_5(i,viable_present_nodes,viable_present_new_nodes,viable_absent_nodes,viable_absent_new_nodes,fix_to=(int(a_t[nodes.index(i)])+1)%2)
            approach_6(i,viable_present_nodes,viable_present_new_nodes,viable_absent_nodes,viable_absent_new_nodes,fix_to=(int(a_t[nodes.index(i)])+1)%2)
        return g

    def fix_to_SS(graph,a,node_set):
        '''
        Takes a damaged graph (that has been run through examine_modifications)
        and a desired attractor, and makes that attractor a SS of an
        edge-modified version of the network. Returns the repaired graph.

        Randomly selects a viable edge modification for every node.

        node_set is a list of the nodes with rules to be modified.

        Returns repaired graph, steady state attractor, and dictionary of
        viable modifications for each node.
        '''
        g_r = graph.copy()
        nodes = sorted(g_r.nodes())
        choice_dict = {}
        for node in node_set:
            if g_r.graph['express'] == node or g_r.graph['knockout'] == node: continue                              # Don't attempt to modify the knocked out/overexpressed node
            choice_dict[node] = [x for x in g_r.graph['modifications'][node] if x[1] == int(a[nodes.index(node)])]  # Choose only from approaches that fix it to the appropriate choice of 0 or 1
            if choice_dict[node] != []:                                                                             # Gang inserted code here to take care of SS failure case 08/31/2014
                modification = choice(choice_dict[node])                                                                # Format: (approach from the 6 listed in examine_modifications(), method 0 or 1 for the block within the method, interacting node 1[, interacting node 2])
            elif choice_dict[node]==[] and set(nodes) <= (set(g_r.node[node]['update_nodes'])|set([node])|set([g_r.graph['express']])|set([g_r.graph['knockout']])):
                return False,a,False
            else:
                raise RuntimeError('Unexpected case in SS repair')

            g_r.node[node]['update_nodes'] += [modification[-1]]                                                    # Append new interacting nodes to previously selected interacting nodes
            new_rules = {}
            bool_suffixes = ['0','1']                                                                           # We always only ever add 1 new node
            for state in g_r.node[node]['update_rules']:
                for bs in bool_suffixes:
                    new_rules[state+bs]=g_r.node[node]['update_rules'][state]       # Append the suffixes but initially with the same outputs as before; then with complete set of input possibilities, change outcomes below
            g_r.node[node]['update_rules'] = new_rules.copy()

            if len(modification) == 4: ex_index = g_r.node[node]['update_nodes'].index(modification[-2])        # Slot in 'update_rule' keys that corresponds to the existing node whose edge is being modified
            #adding an "... OR p_n". So all rules where the final entry is a 1 must have an output of 1
            if modification[0] == 1 and modification[1] == 1: g_r.node[node]['update_rules'] = {key:(val if key[-1] == '0' else 1) for key,val in g_r.node[node]['update_rules'].iteritems()}

            #adding an "... AND a_n". So all rules where the final entry is a 0 must have an output of 0
            elif modification[0] == 1 and modification[1] == 0: g_r.node[node]['update_rules'] = {key:(val if key[-1] == '1' else 0) for key,val in g_r.node[node]['update_rules'].iteritems()}

            #adding an "... OR NOT a_n". So all rules where the final entry is a 0 must have an output of 1
            elif modification[0] == 2 and modification[1] == 1: g_r.node[node]['update_rules'] = {key:(val if key[-1] == '1' else 1) for key,val in g_r.node[node]['update_rules'].iteritems()}

            #adding an "... AND NOT p_n". So all rules where the final entry is a 1 must have an output of 0
            elif modification[0] == 2 and modification[1] == 0: g_r.node[node]['update_rules']  = {key:(val if key[-1] == '0' else 0) for key,val in g_r.node[node]['update_rules'].iteritems()}

            #adding an "... OR (p AND p_n)". So all rules where the existing and new node entries are both 1 must have an output of 1
            elif modification[0] == 3 and modification[1] == 1: g_r.node[node]['update_rules'] = {key:(val if key[ex_index] == '0' or key[-1] == '0' else 1) for key,val in g_r.node[node]['update_rules'].iteritems()}

            #adding an "... AND (a OR a_n)". So all rules where the existing and new node entries are both 0 must have an output of 0
            elif modification[0] == 3 and modification[1] == 0: g_r.node[node]['update_rules'] = {key:(val if key[ex_index] == '1' or key[-1] == '1' else 0) for key,val in g_r.node[node]['update_rules'].iteritems()}

            #adding an "... OR (p AND NOT a_n)". So all rules where the existing and new node entries are '10' must have an output of 1
            elif modification[0] == 4 and modification[1] == 1: g_r.node[node]['update_rules'] = {key:(val if key[ex_index] == '0' or key[-1] == '1' else 1) for key,val in g_r.node[node]['update_rules'].iteritems()}

            #adding an "... AND (a OR NOT p_n)". So all rules where the existing and new node entries are '01' must have an output of 0
            elif modification[0] == 4 and modification[1] == 0: g_r.node[node]['update_rules'] = {key:(val if key[ex_index] == '1' or key[-1] == '0' else 0) for key,val in g_r.node[node]['update_rules'].iteritems()}

            #adding an "... OR (NOT a AND p_n)". So all rules where the existing and new node entries are '01' must have an output of 1
            elif modification[0] == 5 and modification[1] == 1: g_r.node[node]['update_rules'] = {key:(val if key[ex_index] == '1' or key[-1] == '0' else 1) for key,val in g_r.node[node]['update_rules'].iteritems()}

            #adding an "... AND (NOT p OR a_n)". So all rules where the existing and new node entries are '10' must have an output of 0
            elif modification[0] == 5 and modification[1] == 0: g_r.node[node]['update_rules'] = {key:(val if key[ex_index] == '0' or key[-1] == '1' else 0) for key,val in g_r.node[node]['update_rules'].iteritems()}

            #adding an "... OR (NOT a AND NOT a_n)". So all rules where the existing and new node entries are both 0 must have an output of 1
            elif modification[0] == 6 and modification[1] == 1: g_r.node[node]['update_rules'] = {key:(val if key[ex_index] == '1' or key[-1] == '1' else 1) for key,val in g_r.node[node]['update_rules'].iteritems()}

            #adding an "... AND (NOT p OR NOT p_n)". So all rules where the existing and new node entries are both 1 must have an output of 0
            elif modification[0] == 6 and modification[1] == 0: g_r.node[node]['update_rules'] = {key:(val if key[ex_index] == '0' or key[-1] == '0' else 0) for key,val in g_r.node[node]['update_rules'].iteritems()}

        a_r = find_attractor(g_r,a)
        if len(a_r[1]) > 1: raise RuntimeError('LC found in fix_to_SS()')
        return g_r,a_r[1][0],choice_dict




        # End internal functions ---------------------------------------------------
    if method == 'fix_to_SS' or (method == 'LC_repair' and len(a[1]) == 1):     # If we pass a steady state, we evaluate it as a steady state even if we call 'LC_repair'
        if len(a[1]) > 1: a_t = damage_state(graph,a_s)                         # If we do pass a LC, we are concerned with its superset
        else: a_t = a[1][0][:]
        graph = examine_modifications(graph,a_t)                                # Look at the enumerated methods of forcing each node to be a specified state after update. Output is a graph with: g.graph['modifications'][node] properties
        node_set = find_attractor(graph,state=a_t)
        node_set = disjoint(node_set[0],a_t)                                    # Those nodes that change states in the first update from A_t
        g_r,a_r,choice_dict = fix_to_SS(graph,a_t,node_set)                     # Modify edges so the target attractor is a SS of the network
        if choice_dict == False:                                                # Gang inserted code here to deal with the failure case of SS repair 08/31/2014
            return 'SS repair impossible',False
        else:
            pass
        if a_r == a_t: return g_r,choice_dict                                   # The repair succeeded
        else: raise RuntimeError("evaluate_repair() failed.")                   # The repair failed (i.e. there is a bug)


def write_dict_to_file(d,n,fname=False,console_dump=False):
    '''
    The output of evaluate_repair(), if successful, includes a dictionary of all
    viable edge modifications. This function writes this dictionary (d) in a
    human-readable format to a .txt file at fname (if fname != False), and
    dumps the data to the console (if console_dump = True).

    Note that when fixing to a SS, the "adjustments" are permanent, whereas when
    fixing to a LC, the adjustments may not permanently adjust the state of the
    node, but rather ensures that is follows the prescribed oscillations.
    '''
    out = ''
    for node in d.iterkeys():
        out+='Modifications for node '+str(n[node])+':\n'
        for entry in d[node]:
            if entry[0] == 1 and entry[1] == 0: out+= '\tAdjust to %s via: ... AND %s\t(absent_new)\n'%(entry[1],n[entry[-1]])
            elif entry[0] == 1 and entry[1] == 1: out+= '\tAdjust to %s via: ... OR %s\t(present_new)\n'%(entry[1],n[entry[-1]])
            elif entry[0] == 2 and entry[1] == 0: out+= '\tAdjust to %s via: ... AND NOT %s\t(present_new)\n'%(entry[1],n[entry[-1]])
            elif entry[0] == 2 and entry[1] == 1: out+= '\tAdjust to %s via: ... OR NOT %s\t(absent_new)\n'%(entry[1],n[entry[-1]])
            elif entry[0] == 3 and entry[1] == 0: out+= '\tAdjust to %s via: ... AND (%s OR %s)\t(absent, absent_new)\n'%(entry[1],n[entry[-2]],n[entry[-1]])
            elif entry[0] == 3 and entry[1] == 1: out+= '\tAdjust to %s via: ... OR (%s AND %s)\t(present, present_new)\n'%(entry[1],n[entry[-2]],n[entry[-1]])
            elif entry[0] == 4 and entry[1] == 0: out+= '\tAdjust to %s via: ... AND (%s OR NOT %s)\t(absent, present_new)\n'%(entry[1],n[entry[-2]],n[entry[-1]])
            elif entry[0] == 4 and entry[1] == 1: out+= '\tAdjust to %s via: ... OR (%s AND NOT %s)\t(present, absent_new)\n'%(entry[1],n[entry[-2]],n[entry[-1]])
            elif entry[0] == 5 and entry[1] == 0: out+= '\tAdjust to %s via: ... AND (NOT %s OR %s)\t(present, absent_new)\n'%(entry[1],n[entry[-2]],n[entry[-1]])
            elif entry[0] == 5 and entry[1] == 1: out+= '\tAdjust to %s via: ... OR (NOT %s AND %s)\t(absent, present_new)\n'%(entry[1],n[entry[-2]],n[entry[-1]])
            elif entry[0] == 6 and entry[1] == 0: out+= '\tAdjust to %s via: ... AND (NOT %s OR NOT %s)\t(present, present_new)\n'%(entry[1],n[entry[-2]],n[entry[-1]])
            elif entry[0] == 6 and entry[1] == 1: out+= '\tAdjust to %s via: ... OR (NOT %s AND NOT %s)\t(absent, absent_new)\n'%(entry[1],n[entry[-2]],n[entry[-1]])
        out+='\n\n'

    if fname:
        f = open(fname,'wt')
        f.write(out)
        f.close()

    if console_dump:
        print out


