from gurobipy import *
import pandas as pd

from optimize_sigma import get_optimal_sigma
from probmodels import log_prob_present, log_prob_absent

EPSILON = -0.00001

def solve_model(C):

    """
    C -- pandas dataframe with columns correspondng to mutations and rows corresponding to cells
    """
    try:
        vals = []
        # Create a new model
        m = Model("mip1")
        # Add Solver Time Limit
        m.setParam('TimeLimit', 5*60*60)

        Bs = {}
        # Create mutation matrix
        for p in C.index:
            for a in C.columns:
                Bs[(p,a)] = m.addVar(vtype=GRB.BINARY, name="b_{}_{}".format(p,a))
                vals.append((p,a))
        
        for i,a in enumerate(C.columns):
            for b in C.columns[i+1:]:

                x = m.addVar(vtype=GRB.BINARY, name="x_{}_{}".format(a,b))
                y = m.addVar(vtype=GRB.BINARY, name="y_{}_{}".format(a,b))
                z = m.addVar(vtype=GRB.BINARY, name="z_{}_{}".format(a,b))

                m.addConstr(x+y+z <= 2)

                Xps = []
                Yps = []
                Zps = []

                for p in C.index:
                    x_p = m.addVar(vtype=GRB.BINARY, name = "x_{}_{}_{}".format(p,a,b))
                    y_p = m.addVar(vtype=GRB.BINARY, name = "y_{}_{}_{}".format(p,a,b))
                    z_p = m.addVar(vtype=GRB.BINARY, name = "z_{}_{}_{}".format(p,a,b))

                    m.addConstr(x_p >= Bs[(p,a)] + Bs[(p,b)] -1)
                    m.addConstr(y_p >= (1-Bs[(p,a)]) + Bs[(p,b)] -1)
                    m.addConstr(z_p >= Bs[(p,a)] + (1-Bs[(p,b)]) -1)

                    m.addConstr(x_p <= Bs[(p,a)])
                    m.addConstr(x_p <= Bs[(p,b)])

                    m.addConstr(y_p <= 1-Bs[(p,a)])
                    m.addConstr(y_p <= Bs[(p,b)])

                    m.addConstr(z_p <= Bs[(p,a)])
                    m.addConstr(z_p <= 1-Bs[(p,b)])

                    m.addConstr(x >= x_p)
                    m.addConstr(y >= y_p)
                    m.addConstr(z >= z_p)

                    Xps.append(x_p)
                    Yps.append(y_p)
                    Zps.append(z_p)

                m.addConstr(x <= sum(Xps))
                m.addConstr(y <= sum(Yps))
                m.addConstr(z <= sum(Zps))



        # Set objective
        objn = sum(Bs[(p,a)]*C.loc[p][a] for p,a in vals)
        m.setObjective(objn, GRB.MAXIMIZE)

        m.optimize()

        B = C.copy()

        for v in m.getVars():
            if v.varName.startswith('b'):
                try:
                    ''' 
                    if v.varName.startswith('b_ANC:'):                                                                                                   
                        p = '_'.join(v.varName.split('_')[1:3])                                                                                          
                        a = '_'.join(v.varName.split('_')[3:])                                                                                           
                     else:                                                                                                                                
                        p = '_'.join(v.varName.split('_')[1:4])                                                                                          
                        a = '_'.join(v.varName.split('_')[4:])       
                    '''
                    if v.varName.startswith('b_ANC:'):
                        filter_var = v.varName.split(':')[1]
                        p = '_'.join(filter_var.split('_')[0:1])
                        a = '_'.join(filter_var.split('_')[1:])
                        if p.isdigit():
                            p = B.index[int(p)]
                    else:
                        filter_var = v.varName.split('_',1)[1]
                        p = '_'.join(filter_var.split('_')[0:3])
                        a = '_'.join(filter_var.split('_')[3:]) 
                        #p,a = v.varName.split('_')[1:]
                except:
                    print(v.varName)
                    raise
                B.loc[p][a] = v.x

        
        print("Optimized B --------------------")
        print(B)

        deletions = output_with_deletions(C,m)
        return B, deletions
        #print('Obj:', m.objVal)

    except GurobiError:
        print('Error reported')
        raise


def output_with_deletions(C,m):
    deletions = []
    print("C COLUMNS", C.columns)
    print(C)
    for v in C.index:
        if v.startswith('ANC:'):
            for c in C.columns:
                if C.loc[v][c] == EPSILON: 
                    value =  m.getVarByName("b_{}_{}".format(v, c)).x
                    if value == 1: 
                        print("DELETION DETECTED", v,c)
                        deletions.append((v,c))
                    print("---------------------------------------------------------------- EPSILON")
    return deletions
    #raise Exception()



def calculate_C(c, sigmas, DPs, BC):
    # For each copy-number state, we consider take the subet of mutations 
    # that are mixed in that copy-number state
    
    mixed_muts = [a for a in sigmas if sigmas[a][c]=='Mixed']

    C = pd.DataFrame()
    for a in mixed_muts:
        C[a] = BC[BC['c']==c].apply(lambda x: calc_c_observed_cell(x['{}_v'.format(a)], \
               x['{}_t'.format(a)]), axis=1)
    

    print("MIXED MUTS FOR STATE {}: {}".format(c, mixed_muts))
    def desc_scores(v):
        print(v)
        if v == 1: return 100000
        elif v == 0: return -100000
        else: 
            global EPSILON
            print("------------------------------", EPSILON)
            return EPSILON
    

    try:
        descs = DPs[c]
        print("----------------- DESCS", descs)
    except KeyError:
        descs=[]
    
    for i,D in enumerate(descs):
        child, desc = D
        #print(child)
        d = pd.DataFrame([[desc_scores(desc[a]) for a in mixed_muts]], columns = mixed_muts,\
                index = ['ANC:{}'.format(child)])
        C = C.append(d)

    #C = C.reset_index(drop = True)
    return C
        
def calc_c_observed_cell(v,t):
    return log_prob_present(v,t) - log_prob_absent(v,t)

def get_descendent_profiles(sigmas, mutations, S, L):
    DPs = {}
    print(L)
    for edge in S:
        parent, child = edge
        print(edge, L[tuple(edge)])
        child_status = {m:sigmas[m][child] for m in mutations}

        parent_status = {m:sigmas[m][parent] for m in mutations}

        print("child status", child_status)
        print("parent status", parent_status)
        child_founder_profile = {m:1 if child_status[m] == 'Present' and parent_status[m] in ['Mixed', 'Present'] \
                                 else 0 for m in mutations}
        descendent_profile = {m:1 if child_founder_profile[m] == 1 else '?' if m in L[tuple(edge)] \
                              and parent_status[m] in ['Mixed', 'Present'] else 0 for m in mutations}

        print([descendent_profile[v] for v in descendent_profile])
        if parent not in DPs:
            DPs[parent] = []
        DPs[parent].append((child, descendent_profile)) 

    return DPs 

def assemble_mutation_matrix(Bs,sigmas, BC, mutations):

    # Bs give partial information and sigmas assemble the rest
    B_tot = pd.DataFrame(columns = mutations, index = BC.index)
    #print(mutations)
    #print(Bs)

    for i in Bs:

        B = Bs[i]
        for p in B.index:
            v = B.loc[p]
            if p.startswith('ANC:'): continue
            for a in B.columns:
                B_tot.loc[p][a] = int(v[a])

        status = sigmas.loc[i]
        cells = BC[BC['c'] == i].index

        for cell in cells:
            for a in status.index:
                if status[a] == 'Absent':
                    B_tot.loc[cell][a] = 0
                elif status[a] == 'Present':
                    B_tot.loc[cell][a] = 1

    return B_tot

def assemble_mutation_matrix_with_ancestors(Bs,sigmas, BC, mutations):

    
    index = []
    for i in Bs:
            index+=Bs[i].index.tolist()
            print(Bs[i].index)


    # Bs give partial information and sigmas assemble the rest
    B_tot = pd.DataFrame(columns = ["CN"]+ mutations, index = index)
    #B_tot = pd.DataFrame(columns = mutations, index = index)

    for i in Bs:

        B = Bs[i]
        for p in B.index:
            v = B.loc[p]
            #if p.startswith('d'): continue
            for a in B.columns:
                B_tot.loc[p][a] = int(v[a])

        status = sigmas.loc[i]
        #cells = BC[BC['c'] == i].index
        cells = B.index

        for cell in cells:
            #B_tot["CN"]=i
            B_tot.loc[cell]['CN'] = i
            for a in status.index:
                if status[a] == 'Absent':
                    B_tot.loc[cell][a] = 0
                elif status[a] == 'Present':
                    B_tot.loc[cell][a] = 1
    return B_tot



  
if __name__ == '__main__':

    import pandas as pd


    input_file = "test_data/BC.csv"
    state_tree = "test_data/S.csv"
    BC, S, L = read_in_files(input_file, state_tree)
    mutations = set(sorted([v.rsplit('_', 1)[0] for v in BC.columns[1:]]))
    sigmas = get_optimal_sigma(S,BC)

    DPs = get_descendent_profiles(sigmas, mutations)

  
    C = calculate_C(1, sigmas, DPs)

    B = solve_model(C)
    
    #print(B)
