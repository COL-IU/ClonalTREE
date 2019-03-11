from sys import *
from decimal import *
from itertools import *
from copy import *
from math import *
from time import *
from datetime import *
import numpy as np

one = Decimal('1.0')
zero = Decimal('0.0')
low = Decimal('0.00001')
FAIL_THRESHOLD = Decimal('0.0')


def remove_first_rowcol(mat):
    out = []
    for i in range(1, len(mat)):
        out.append(mat[i][1:])
    return out


def add_founder(F):
    my_F = [[zero for _ in range(len(F[0]) + 1)]]
    my_F[0][0] = one
    for row in F:
        my_F.append([one] + row)
    return my_F


def ancestry_graph(F):
    m = len(F)
    n = len(F[0])
    G = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(0, n):
        edges = set(list(range(n)))
        for j in range(0, m):
            temp = set()
            for k in range(0, n):
                if F[j][k] <= F[j][i]:
                    temp.add(k)
            edges = edges & temp
        for e in edges:
            G[e][i] = 1
    return G


def enum_spanning_trees(G_i):
    G = deepcopy(G_i)
    n = len(G)
    for i in range(0, n):
        G[i][i] = 0
    ones = []
    for i in range(1, n):
        indices = [j for j, x in enumerate(G[i]) if x == 1]
        ones.append(indices)
    num_trees = 1
    for item in ones:
        num_trees = num_trees * len(item)
    # stdout.write(" " + str(num_trees) + " trees\n")
    # stdout.flush()
    sTrees = []
    i = 0
    for element in product(*ones):
        # if i % 1000000 == 0:
        #     stdout.write(" i1 = " + str(i/1000000) + " M; " + str(datetime.now()) + "\n")
        #     stdout.flush()
        # i += 1
        sTrees.append(list_to_tree([0]+list(element)))
    return sTrees


def calc_dist(l):
    dist = []
    s = Decimal(sum(l))
    for i in l:
        dist.append(Decimal(i) / s)
    return dist


def list_to_tree(parents):
    T = [[0]*len(parents)]
    T[0][0] = 1
    for i in range(1, len(parents)):
        parent = parents[i]
        T.append(T[parent][:])
        T[i][i] = 1
    return T


def list_to_clones(parents, variants):
    clones = [[]]
    for i in range(1, len(parents)):
        parent = parents[i]
        clones.append(clones[parent][:])
        clones[i].append(variants[i-1])
    return clones[1:]


def tree_to_list(T):
    n = len(T)
    parents = [0] * n
    children = [[] for _ in range(n)]
    for i in range(1, n):
        row = T[i][:]
        row[i] = 0
        for j in range(0, n):
            if T[j] == row:
                parents[i] = j
                children[j].append(i)
    return parents, children


def get_c(F, T, fail_threshold):
    # Assumes square F
    n = len(T)
    parents, children = tree_to_list(T)
    C = [[zero for _ in range(n)] for _ in range(n)]
    C[0][0] = one
    fail = False
    for i in range(1, n):
        for j in range(0, n):
            c = children[j]
            cs = 0
            for child in c:
                cs += F[i][child]
            C[i][j] = F[i][j] - cs
            if C[i][j] < fail_threshold:
                fail = True
    if fail:
        p = -1
    else:
        p = one
        for i in range(1, n):
            p = p * C[i-1][parents[i]]
    return C, p


def get_p(F, parents, fail_threshold):
    tree = list_to_tree(parents)
    _, p = get_c(F, tree, fail_threshold)
    return p


def get_f(C, T):
    m = len(C)
    if m == 0:
        return []
    n = len(C[0])
    F = [[Decimal('0.0') for _ in range(n)] for _ in range(m)]
    for i in range(0, m):
        for j in range(0, n):
            temp = Decimal('0.0')
            for k in range(0, n):
                temp += (C[i][k]*T[k][j])
            F[i][j] = temp
    return F


def write_im(mat, out):
    m = len(mat)
    n = len(mat[0])
    for i in range(0, m):
        for j in range(0, n):
            out.write(str(mat[i][j]) + "\t")
        out.write("\n")


def write_dm(mat, out):
    m = len(mat)
    if m == 0:
        return
    n = len(mat[0])
    for i in range(0, m):
        for j in range(0, n):
            out.write("{0:.5f}".format(mat[i][j]) + "\t")
        out.write("\n")


def read_F(infile):
    F = []
    lines = infile.readlines()
    var = lines[0].split(None)
    for i in range(1, len(lines)):
        line = lines[i]
        words = line.split(None)
        F.append(list(map(Decimal, words)))
    return F, var


def read_dm(infile):
    out = []
    lines = infile.readlines()
    for line in lines:
        words = line.split(None)
        out.append(list(map(Decimal, words)))
    return out


def print_dm(mat):
    # Assumes square mat
    m = len(mat)
    n = len(mat[0])
    for i in range(0, m):
        for j in range(0, n):
            stdout.write("{0:.5f}".format(mat[i][j]) + "\t")
        stdout.write("\n")


def print_dl(l):
    for item in l:
        stdout.write("{0:.3f}".format(item) + "\t")
    stdout.write("\n")


def dl_to_str(l):
    out = []
    for item in l:
        out.append("{0:.3f}".format(item))
    return "\t".join(out)


def print_tree(T):
    n = len(T)
    for i in range(0, n):
        for j in range(0, n):
            stdout.write(str(T[i][j]) + " ")
        stdout.write("\n")


def get_children(parents, parent):
    children = []
    for i in range(1, len(parents)):
        if parents[i] == parent:
            children.append(i)
    return children


def read_variants_set(f):
    out = set()
    for line in f:
        out.add(int(line.strip()))
    return out


def read_variants_list(f):
    out = []
    for line in f:
        out.append(int(line.strip()))
    return out


def reorder(l, order):
    return [l[i] for i in order]


def rearrange_rows(mat1, variants1, variants2):
    mat2 = []
    for var in variants2:
        mat2.append(mat1[variants1.index(var)])
    return mat2


def rearrange(F, variants, remove_redundancy=False):
    F_out = deepcopy(F)
    m = len(F_out)
    if m == 0:
        return F, variants
    if len(F[0]) == 0:
        return F, variants
    # Make F_out a diagonal (step) matrix
    F_out = list(map(list, zip(*F_out)))
    order = []
    i = 0
    for row in F_out:
        for i in range(0, len(row)):
            if row[i] > zero:
                break
        order.append(i)
    indices = list(range(0, len(F_out)))
    [order, indices, F_out] = list(zip(*sorted(zip(order, indices, F_out), key=lambda x: x[0])))
    variants_out = reorder(variants, indices)

    '''
    j = 0
    # Fill possibly non-zero frequencies with low frequency
    for i in range(0, len(F_out)):
        for j in range(len(F_out[0]) - 1, -1, -1):
            if F_out[i][j] > zero:
                break
        for k in range(order[i], j):
            if F_out[i][k] == zero:
                F_out[i][k] = low
    '''

    F_out = list(map(list, zip(*F_out)))

    if remove_redundancy:
        # Remove unnecessary time points
        z = 0
        z_prev = -1
        new_F = []
        for j in range(0, len(F_out)):
            row = F_out[j]
            for z in range(len(row)-1, -2, -1):
                if row[z] > zero:
                    break
            if z - z_prev > 0:
                new_F.append(F_out[j])
            z_prev = z
        F_out = new_F

    return F_out, variants_out


def get_step_structure(F):
    # Assumes step structured F
    # Also assumes that redundant time points are removed
    m = len(F)
    if m == 0:
        return []
    n = len(F[0])
    multi_spawns = []
    z_prev = -1
    for i in range(m):
        z = n - 1
        for z in range(n - 1, -1, -1):
            if F[i][z] > 0:
                break
        multi_spawns.append(range(z_prev + 1, z + 1))
        z_prev = z
    return list(map(list, multi_spawns))


def squarify(F, step_structure):
    # Assumes F's step structure matches with step_structure
    if not F:
        return []
    m = len(F)
    n = len(F[0])
    new_order = list(chain.from_iterable(step_structure))
    r = len(new_order)
    if r < n:
        new_order = new_order + list(range(r, n))
    F_out = []
    k = 0
    for i in range(0, len(step_structure)):
        tup = step_structure[i]
        reordered = reorder(F[i], new_order)
        for j in range(0, len(tup)):
            F_out.append(reordered[:k] + [F[i][x] for x in tup[:(j + 1)]] + [zero] * (n - j - 1 - k))
        k += len(tup)
    if r < n:
        for i in range(len(step_structure), m):
            F_out.append(F[i])
    return F_out


def square_forms(F, variants):
    F1, variants1 = rearrange(F, variants, True)
    m = len(F1)
    if m == 0:
        return [], []

    # Getting the step structure
    step_structure = get_step_structure(F1)
    num_matrices = 1
    for item in step_structure:
        num_matrices = num_matrices * factorial(len(item))

    perm = map(permutations, step_structure)
    prod = product(*perm)
    stdout.write(str(num_matrices) + " matrices\n")
    stdout.flush()
    i = 0
    for order in prod:
        if i % 100000 == 0:
            stdout.write("i = " + str(i) + "; " + str(datetime.now()) + "\n")
            stdout.flush()
        i += 1
        F_yield = squarify(F1, order)
        yield F_yield, reorder(variants1, list(chain.from_iterable(order)))


def sub_f(in_F_file, in_var_file, sub_var_file, out_F_file, out_var_file):
    f1 = open(in_F_file)
    f2 = open(in_var_file)
    f3 = open(sub_var_file)
    f4 = open(out_F_file, "w")
    f5 = open(out_var_file, "w")

    in_F = read_dm(f1)
    in_var = read_variants_list(f2)
    sub_var = read_variants_list(f3)

    in_F = list(map(list, zip(*in_F)))
    out_F = []
    out_var = []
    for i in range(0, len(in_var)):
        var = in_var[i]
        if var in sub_var:
            out_var.append(var)
            out_F.append(in_F[i])
    out_F = list(map(list, zip(*out_F)))
    out_F, out_var = rearrange(out_F, out_var, True)
    write_dm(out_F, f4)
    for var in out_var:
        f5.write(str(var) + "\n")
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()


def remap_parents(indices, parents):
    remapped = [0]*len(parents)
    for i in range(1, len(parents)):
        if parents[i] != 0:
            remapped[indices[i-1]] = indices[parents[i]-1]
    return remapped


def read_results(f_path):
    f = open(f_path)
    out = list(map(lambda x: list(map(int, x.split(None))), f.readlines()))
    f.close()
    return out


def read_clones(f_path):
    f = open(f_path)
    out = []
    lines = f.readlines()
    for line in lines:
        out.append(list(map(lambda x: [0] + list(map(int, x.split("-"))), line.split(";"))))
    f.close()
    return out


def print_parents(var, parents, out):
    out.write("Clone\tParent\n")
    for i in range(0, len(var)):
        out.write(str(var[i]) + "\t")
        if parents[i+1] == 0:
            out.write("0\n")
        else:
            out.write(str(var[parents[i+1]-1]) + "\n")


def white_list(variants, clones):
    out = {}
    for variant in variants:
        temp = set()
        for clone in map(set, clones):
            if variant in clone:
                temp.update(clone)
        if len(temp) == 0:
            temp = set(variants)
        temp.remove(variant)
        out[variant] = temp
    return out


def black_list(variants, clones):
    out = {}
    for variant in variants:
        out[variant] = set()
    clones_sets = list(map(set, clones))
    for i in range(0, len(clones_sets)-1):
        clone1 = clones_sets[i]
        for j in range(i+1, len(clones_sets)):
            clone2 = clones_sets[j]
            inter = clone1.intersection(clone2)
            if len(inter) > 0:
                sym_diff = clone1.symmetric_difference(clone2)
                for variant in sym_diff:
                    out[variant].update(inter)
    return out


def translate_clones(variants, clones):
    out = []
    for clone in clones:
        temp = [0]
        for v in clone:
            if v in variants:
                temp.append(variants.index(v)+1)
        out.append(temp)
    return out


def sample_clones(truth, n):
    nodes = list(np.random.choice(range(1, len(truth)), replace=False, size=n))
    clones = []
    for node in nodes:
        clone = []
        temp = node
        while temp != 0:
            clone.append(temp)
            temp = truth[temp]
        clones.append(clone)
    return clones


def write_clones(clones, out_file):
    reps = []
    for clone in clones:
        clone_string = "-".join(list(map(str, clone)))
        reps.append(clone_string)
    out_file.write(";".join(reps) + "\n")


def parents_to_clones(parents):
    clones = set()
    for i in range(0, len(parents)):
        clone = set()
        clone.add(i)
        temp = parents[i]
        while temp != 0:
            clone.add(temp)
            temp = parents[temp]
        clone.add(0)
        clone = frozenset(clone)
        clones.add(clone)
    return clones
