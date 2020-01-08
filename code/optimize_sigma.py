
from probmodels import log_prob_absent, log_prob_present, log_prob_mixed
import pandas as pd

###
#   Calculates the optimal sigma assignments for all mutations given state tree S
###
def get_optimal_sigma(S, BC, L):
    """
    For all mutations calculates the optimal sigma assignment
    S -- edgelist representation of copy-number state tree

    Note that there's an inherent assumption about C that the copy-number states are 
    continuous integers starting at 0

    BC -- pandas dataframe input matrix

    returns dataframe with columns corresponding to mutations and rows corresponding to
    copy-number states, and entries in {'Absent', 'Present', 'Mixed'}

    """
    mutations = set(sorted([v.split('_')[0] for v in BC.columns[1:]]))
    C = list(BC['c'])
    num_states = len(set(C))
    subtrees = enum_all_subtrees(num_states, S)

    sigmas = {}
    print L

    deletions = []
    for mutation in mutations:
        max_value = float('-inf')
        max_sigma = None
        max_deletions = None
        V = list(BC['{}_v'.format(mutation)])
        T = list(BC['{}_t'.format(mutation)])
        for subtree in subtrees:

            status1 = ['Mixed' if i == subtree[0] else 'Present' if i in subtree[1:] else 'Absent' for i in range(num_states)]
            
            # S is the edge list
            # For every edge, I need to check that if I go from present to absent
            # then there is an allowed loss
            valid_tree = True
            tree_deletions = []
            for edge in S:
                s,t = edge
                status_s = status1[s]
                status_t = status1[t]
                if status_s == 'Present' and status_t == 'Absent':
                    if mutation not in L[tuple(edge)]: 
                        valid_tree = False
                    else:
                        tree_deletions.append(('ANC:{}'.format(t), mutation))
                    
            if not valid_tree: continue
                
            p1 = log_prob_sigma(V,T,C, status1)

            if p1 > max_value:
                max_value = p1
                max_sigma = status1
                max_deletions = tree_deletions
    
        deletions += max_deletions
        sigmas[mutation] = max_sigma

    sigma_new = pd.DataFrame(sigmas)
    print "SIGMA"
    print sigma_new

    return sigma_new, deletions

###
#   Calculates the probability of a given sigma for one mutation
###
def log_prob_sigma(V,T,C,sigma):
    """
    V,T,C -- lists of length n (num cells) whose values correspond to the number of
                variant reads, total reads, and copy-number state assignment 
    sigma -- list of length k (number of copy-number states) whose entries are in 
                {'Absent', 'Present', 'Mixed'}

    Note that there's an inherent assumption about C that the copy-number states are 
    continuous integers starting at 0

    returns log probability for sigma
    """

    log_prob = 0
    for i,R in enumerate(zip(V,T,C)):
        # print '-------'
        v,t,c = R

        status = sigma[c]
        if status == 'Absent':
            log_prob += log_prob_absent(v,t)
        elif status == 'Present':
            log_prob += log_prob_present(v,t)
        elif status == 'Mixed':
            log_prob += log_prob_mixed(v,t)
        else:
            raise Exception('No such status: {}'.format(status))
    return log_prob


###
#   Enumerate possible colorings
###
def enum_all_subtrees(num_states, S):
    edgelist = S[:]
    all_subtrees = []
    for i in range(num_states):
        subtree = [i]
        frontier = [(s,t) for s,t in edgelist if s in subtree and t not in subtree]
        subtrees = enum_rooted_subtrees(subtree, edgelist, frontier)
        all_subtrees += subtrees

    return all_subtrees

def enum_rooted_subtrees(subtree, edgelist, frontier, treelist = []):
    if len(treelist) == 0: treelist = [subtree]
    while len(frontier) > 0:
        edge = frontier.pop()
        frontier_new = frontier + [(s,t) for s,t in edgelist[:] if s == edge[1]]
        subtree_new = subtree+[edge[1]]
        treelist.append(subtree_new)
        treelist = enum_rooted_subtrees(subtree_new, edgelist, frontier_new, treelist)
    return treelist   



if __name__ == "__main__": 
    
    import pandas as pd

    # TEST

    print "Testing Opitimize Sigma"
    input_file = "test_data/BC.csv"
    state_tree = "test_data/S.csv"
    BC = pd.read_csv(input_file)
    S = []
    L = {}



    with open(state_tree) as f:
        for line in f:
            line_split = line.strip().split(',')
            edge = map(int, line_split[:2])
            S.append(edge)
            try: 
                L[tuple(edge)] = line_split[2:]
            except:
                L[tuple(edge)] = []
    print "Input data"
    print "B,C"
    print BC
    print "State tree S"
    print S
    print "Result"
    print get_optimal_sigma(S,BC)
