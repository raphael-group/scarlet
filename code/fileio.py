import pandas as pd

def read_in_files(input_file, state_tree_file):
    BC = pd.read_csv(input_file, index_col=0)
    S = []
    L = {}
    mutations = set(sorted([v.split('_')[0] for v in BC.columns[1:]]))
    BC.index = map(str, BC.index)

    with open(state_tree_file) as f:
        for line in f:
            line_split = line.strip().split(',')
            edge = list(map(int, line_split[:2]))
            S.append(edge)
            try: 
                L[tuple(edge)] = line_split[2:]
            except:
                L[tuple(edge)] = []

    return BC, S, L

def write_out_files(B, B_with_ancestors, T, filename, totalLL):
    B.to_csv('{}.B'.format(filename))
    B_with_ancestors.to_csv('{}.B_ancestor'.format(filename))
    T.to_csv('{}.T'.format(filename))

    with open('{}.LL'.format(filename), 'w') as out:
        out.write('{}\n'.format(totalLL))
