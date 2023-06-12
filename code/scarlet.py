from fileio import read_in_files, write_out_files
from optimize_sigma import get_optimal_sigma
from optimize_mutation_matrix import get_descendent_profiles, calculate_C, solve_model, assemble_mutation_matrix, assemble_mutation_matrix_with_ancestors
import pandas as pd
from probmodels import compute_LL_solution 

import sys
def correct_ternary_matrix(B, S, BC, deletions):
    T = B.copy()
    
    for deletion in deletions:
        print("  ----", deletion)
        child, mut = deletion
        #print('DELETION: ', deletion, child, mut)
        child = int(child.replace('ANC:', ''))
        # calculate the set of copy-number states affected by the deletion
        # not the most efficient way to do this but S should be small
        set_states = [child]
        while True:
            set_states_new = set_states + [t for s,t in S if s in set_states ]
            if set(set_states) == set(set_states_new): break
            set_states = set_states_new


        # This is the set of states where the mutation has been deleted. Thus it gets set to be
        # a two in the output data
        
        cells = BC[BC['c'].apply(lambda x: x in set_states)].index
        T[mut] = T.apply(lambda x: x if x.name not in cells else 2 , axis=1)[mut]
    return T

def main():
    BC_file = sys.argv[1]
    SL_file = sys.argv[2]
    
    #TODO: Add input file containing list of germline mutations (SNPs)
    SNP_file = sys.argv[3]
    
    output_file = sys.argv[4]

    BC, S, L, SNP = read_in_files(BC_file, SL_file, SNP_file)
    
    mutations = sorted(set([v.rsplit('_', 1)[0] for v in BC.columns[1:]]))
    sigmas, dels = get_optimal_sigma(S,BC,L,SNP)
    DPs = get_descendent_profiles(sigmas, mutations, S, L)

    all_deletions = dels
    cn_states = BC['c'].unique()
    Bs = {}
    
    for i in cn_states:
        
        C= calculate_C(i, sigmas, DPs, BC)
        
        B, deletions = solve_model(C)

        Bs[i]=B
        all_deletions += deletions

    print("All DELETIONS", all_deletions)

    result = assemble_mutation_matrix(Bs, sigmas, BC, mutations)
    result_with_ancestors = assemble_mutation_matrix_with_ancestors(Bs, sigmas, BC, mutations)

    ternary = result.copy()
    for deletion in all_deletions:
        ternary = correct_ternary_matrix(ternary, S, BC, all_deletions)


    solutionLL = compute_LL_solution(BC, result, mutations)
    write_out_files(result, result_with_ancestors, ternary, output_file, solutionLL)





if __name__ == '__main__':

    main()


