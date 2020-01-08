# read in results file
import sys
import pandas as pd

###
# Sorts columns such that if mutation i is ancestral to mutation j, then i < j
# Used by construct_perfect_phylogeny
###
def sort_columns_by_containment(matrix, cn_state, result):
    # mutations are ignored if none of the cells in the subtree contain the mutation
    # or if the parent of the root contains the mutation
    # because then it is not gained in this subtree
    def contains(v1, v2):
        return all([x >= y for x, y in zip(v1, v2)])
    values = []
    counts = []
    for column1 in matrix.columns[1:]:
        if column1 == 'CN': continue
        if all([x == 0 for x in matrix[column1]]): continue

        # need to check the parent. If the parent of the root does not contain the mutation,
        # then it is gained here
        try:
            parent_state = result[column1]['ANC:{}'.format(cn_state)]
            if all([x == 1 for x in matrix[column1]]) and parent_state != 0: continue
        except KeyError:
            if all([x == 1 for x in matrix[column1]]) and cn_state != 0: continue
        values.append(column1)
        count = 0
        for column2 in matrix.columns[1:]:
            if column2 == 'CN': continue
            if contains(matrix[column1], matrix[column2]):
                count += 1
        counts.append(len(matrix.columns) - count)

    columns_sorted = [y for x, y in sorted(zip(counts, values))]
    matrix_sorted = matrix[columns_sorted]
    return matrix_sorted

def construct_perfect_phylogeny(result, cn_state, vertex_state):
    """
    INPUT: perfect phylogeny matrix (pd.DataFrame)
    OUTPUT: perfect phylogeny tree (dict, parent -> [list of children])

    Vertices will be prepended by either
    - "CELL:" indicates that the vertex is a leaf of the subtree -- either an observed cell or
        parent of descendant
    - "MUT:" indicates that the vertex is an internal vertex
    - "ROOT:" indicates that the vertex is the root of the subtree
    """
    matrix = result[result['CN'] == cn_state]
    matrix_sorted = sort_columns_by_containment(matrix, cn_state, result)

    root_name = 'ROOT:{}'.format(cn_state)
    ###
    #   Assign copy-number state to all vertices. The set of vertices corresponds to the leaves
    #   (the rows of the matrix) and the internal vertices (the columns of the matrix) as well as the root
    ###
    for v in matrix_sorted.columns: vertex_state[v] = cn_state
    for v in matrix_sorted.index: vertex_state[v] = cn_state
    vertex_state[root_name] = cn_state

    ###
    # Construct perfect phylogeny tree from sorted matrix
    # Using the prefix method. Every observed cell corresponds to a path from the root
    # to the leaves. As the columns are sorted by containment, the sequence of mutations
    # contained in the cell in the matrix is the order in the tree
    ###
    tree = {}
    root = 'ROOT:{}'.format(cn_state)
    tree[root] = []
    for cell_id in matrix_sorted.index:
        cell = matrix_sorted.loc[cell_id]
        prefix = [v for v in cell.index if cell[v] == 1]
        try:
            if prefix[0] not in tree[root]: tree[root].append(prefix[0])
        except IndexError:
            tree[root].append(cell_id)
            continue

        for j, value in enumerate(prefix[1:]):
            if prefix[j] not in tree:
                tree[prefix[j]] = []
            if value not in tree[prefix[j]]:
                tree[prefix[j]].append(value)

        if prefix[-1] not in tree: tree[prefix[-1]] = []
        tree[prefix[-1]].append(cell_id)
    return tree, vertex_state

def identify_mutation_losses(result):
    ## A mutation is lost if it is present in the ancestor and not in the root
    ## It is not in the root if it is not in any of the parents
    mutation_losses = {}
    for cn_state in cn_states:
        try:
            parent_state = result.loc['ANC:{}'.format(cn_state)]
        except KeyError:
            continue
        matrix = result[result['CN'] == cn_state]
        mutation_losses[str(cn_state)] = []
        for column in matrix.columns:
            if column == 'CN': continue
            absent_in_child = all([x == 0 for x in matrix[column]])
            present_in_parent = parent_state[column] == 1
            if absent_in_child and present_in_parent:
                mutation_losses[str(cn_state)].append(column.split(':')[-1])

    return mutation_losses


def construct_full_tree(result, cn_tree):

    subtrees = {}
    vertex_state = {}
    for i, cn_state in enumerate(result['CN'].unique()):
        tree, vertex_state = construct_perfect_phylogeny(result, cn_state, vertex_state)
        subtrees[i] = tree

    # join subtrees:
    global_tree = {}
    for tree in subtrees:
        global_tree.update(subtrees[tree])

    # for every edge (i,j) in the CN tree, join dj in tree[i] with ROOT:j in tree[j]
    for i, j in cn_tree:
        # get rid of dj, have parent of dj connect directly to root
        vertex = 'ANC:{}'.format(j)
        for key in global_tree:
            if vertex in global_tree[key]:
                children = global_tree[key]
                children.remove(vertex)
                children.append('ROOT:{}'.format(j))
                global_tree[key] = children
                break
        #global_tree['CELL:d{}'.format(j)] = ['ROOT:{}'.format(j)]
    return global_tree, subtrees, vertex_state

def write_out_tree(global_tree, tree_file):
    with open(tree_file, 'w') as out:
        for key in global_tree:
            children = global_tree[key]
            for value in children:
                parent = "{} ({})".format(key, vertex_colors[key])
                color_child = vertex_colors[value]
                child = "{} ({})".format(value, color_child)
                out.write("{},{}\n".format(parent, child))

def read_inputs():
    if len(sys.argv) < 5:
        print "USAGE: plot_tree.py output_file.B_ancestor [CN Tree file] [plotting style] [output prefix]"
    output_file = sys.argv[1]
    tree_file = sys.argv[2]
    plotting_style = sys.argv[3].upper()
    output_prefix = sys.argv[4]

    if plotting_style not in ['ALL', 'COUNT', 'NONE']:
        raise ValueError
    result = pd.read_table(output_file, sep=',', index_col=0)
    cn_tree = []
    with open(tree_file) as f:
        for line in f:
            cn_tree.append(line.strip().split(',')[:2])

    return result, cn_tree, output_prefix, plotting_style

def add_all_leaves(result, out, vertex_colors):
    for v in result.index:
        if v.startswith('ANC:'): continue
        label = v.split(':')[-1]
        out.write(
            '\"{}\" [label = \"{}\", shape=circle, height = 0.3, width = 0.3, margin=0, color = {}, fontsize = 20]; \n'.format(
                v, label, vertex_colors[v]))


def output_dot_file(output_prefix, result, global_tree, vertex_colors, cn_states, draw_leaves, mutation_losses):

    ## Draw leaves is in {'ALL', 'COUNT', 'NONE'}
    if draw_leaves not in ['ALL', 'COUNT', 'NONE']:
        raise ValueError
    with open('{}.dot'.format(output_prefix), 'w') as out:
        out.write('digraph g{\n')
        if draw_leaves == 'ALL' : out.write('ratio=0.5\n')
        else: out.write('ratio=1.5\n')
        out.write('nodesep=0\n')
        out.write('graph [fontname = \"helvetica\", colorscheme=set19];\n')
        out.write('node [fontname = \"helvetica\", colorscheme=set19];\n')
        out.write('edge [fontname = \"helvetica\", colorscheme=set19];\n')

        if draw_leaves == 'ALL':
            add_all_leaves(result, out, vertex_colors)

        ## Add nodes for internal vertices. These are represented as a small point to save space
        for v in result.columns:
            if v != 'CN':
                out.write('\"{}\" [label = \"\", color = {}, height = 0.3, width = 0.3, shape=point]; \n'.format(v,vertex_colors[v]))
        for v in cn_states:
            out.write('\"ROOT:{}\" [label = \"\", color = {}, height = 0.3, width = 0.3, shape=point]; \n'.format(v, v))

        ## Add edges
        for parent in global_tree:
            children = global_tree[parent]
            for child in children:
                ### Dashed, narrow line to observed cells to de-emphasize
                if child.startswith('CELL:'):
                    if draw_leaves == 'ALL':
                        out.write("\"{}\" -> \"{}\" [penwidth=1, style=dashed, color = \"{};0.5:{}\"];".format(parent,\
                                                    child, vertex_colors[parent], vertex_colors[child]) + '\n')
                ### This is between two different copy-number states and thus needs two different colors for the edges
                ### TODO: mutation losses
                elif child.startswith('ROOT:'):
                    state = child.split(':')[-1]
                    if len(mutation_losses[state]) > 0: losses = 'LOSS: '+'\nLOSS: '.join(mutation_losses[state])
                    else: losses = ""
                    out.write("\"{}\" -> \"{}\" [label = \"{}\", fontsize=30,fontcolor=\"#8b0000\", penwidth=3, color = \"{};0.5:{}\"];".format(parent, child, \
                                                    losses, vertex_colors[parent], vertex_colors[child]) + '\n')
                ### Edges within a copy-number state. All labeled by mutations
                else:
                    label = child.split(':')[-1]
                    out.write("\"{}\" -> \"{}\" [penwidth=3, label = \"{}\", fontsize = 30, color=\"{};0.5:{}\"];".format(parent,\
                                                    child, label, vertex_colors[parent], vertex_colors[child]) + '\n')
            if draw_leaves == 'COUNT':
                leaves = [child for child in children if child.startswith('CELL')]
                if len(leaves) > 0:
                    out.write("\"{}\" [label = \"{}\", color = {}, fontsize=30]".format('CHILDREN:{}'.format(parent), len(leaves),
                                                                           vertex_colors[parent]))

                    out.write("\"{}\" -> \"{}\" [penwidth=1, style=dashed, fontsize=30, color = \"{};0.5:{}\"];".format(parent, \
                                                                                                           "CHILDREN:{}".format(
                                                                                                               parent),
                                                                                                           vertex_colors[
                                                                                                               parent],
                                                                                                           vertex_colors[
                                                                                                               parent]) + '\n')

        out.write('}\n')


if __name__ == "__main__":


    result_matrix, cn_tree, output_prefix, draw_leaves = read_inputs()
    result_matrix.columns = ['MUT:{}'.format(v) if v != 'CN' else v for v in result_matrix.columns ]
    result_matrix['CELL_ID'] = ['CELL:{}'.format(v) if 'ANC:' not in v else v for v in result_matrix.index]
    result_matrix = result_matrix.set_index('CELL_ID')

    global_tree, subtrees, vertex_colors = construct_full_tree(result_matrix, cn_tree)

    tree_filename = "{}.edgelist".format(output_prefix)
    print "Outputting edgelist to {}".format(tree_filename)
    write_out_tree(global_tree, tree_filename)
    cn_states = result_matrix['CN'].unique()
    mutation_losses = identify_mutation_losses(result_matrix)
    dot_filename = "{}.dot".format(output_prefix)
    print "Outputting DOT file to {}".format(dot_filename)
    output_dot_file(output_prefix, result_matrix, global_tree, vertex_colors, cn_states, draw_leaves, mutation_losses)


