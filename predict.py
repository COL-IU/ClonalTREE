from myutils import *
import numpy as np
from time import *


def penalty(parents, F):
    # F may not be square
    m = len(F)
    pen = 0
    cur_parent = parents[-1]
    children = get_children(parents, cur_parent)
    for t in range(cur_parent, m):
        children_sum = zero
        for child in children:
            children_sum += F[t][child]
        temp = F[t][cur_parent] - children_sum
        if temp < 0:
            pen = pen + (10 * temp)
    return pen


def valid_parent_value(parents, F, fail_threshold):
    # F may not be square
    m = len(F)
    vp = True
    cur_parent = parents[-1]
    children = get_children(parents, cur_parent)
    for t in range(cur_parent, m):
        children_sum = zero
        for child in children:
            children_sum += F[t][child]
        if F[t][cur_parent] - children_sum < fail_threshold:
            vp = False
            return vp
    return vp


def valid_parent_order(order_validity, parents):
    cur_variant = len(parents)-1
    cur_path = deepcopy(order_validity[cur_variant])
    temp = parents[cur_variant]
    while temp != 0:
        cur_path.discard(temp)
        temp = parents[temp]
    cur_path.discard(0)
    return len(cur_path) == 0


def c_row(f_row, parents):
    m = len(f_row)
    C = []
    for parent in range(0, m):
        children = get_children(parents, parent)
        children_sum = zero
        for child in children:
            children_sum += f_row[child]
        C.append(f_row[parent] - children_sum)
    return C


def smart_predict_original(F, m, fail_threshold, clones=[]):
    # Assumes square F
    if len(F) == 0:
        return [], [], [], -1

    # Add founder
    my_F = add_founder(F)
    # m = len(my_F)
    choices_stack = [[0]]
    chosen_parents = [0]

    valid_parents = white_list(list(range(m)), clones)
    order_validity = black_list(list(range(m)), clones)

    success = False
    while choices_stack:
        # print(choices_stack)
        chosen_parents = choices_stack.pop()
        i = len(chosen_parents) - 1
        if not valid_parent_value(chosen_parents, my_F, fail_threshold):
            continue
        if clones:
            if not valid_parent_order(order_validity, chosen_parents):
                continue
        if i == (m - 1):
            success = True
            break
        C_row = c_row(my_F[i], chosen_parents)
        next_choices = list((np.array(C_row)).argsort())
        for next_choice in next_choices:
            if C_row[next_choice] > fail_threshold and next_choice <= i and next_choice in valid_parents[i+1]:
                temp = chosen_parents[:]
                temp.append(next_choice)
                choices_stack.append(temp)
    if not success:
        return [], my_F, [], -1

    C = []
    for i in range(0, m):
        C.append(c_row(my_F[i], chosen_parents))

    p = one
    for i in range(1, m):
        p = p * C[i-1][chosen_parents[i]]
    # choices_stack.append(chosen_parents)
    return C, my_F, chosen_parents, p


def scores(F, parents):
    i = len(parents) - 1
    C_row = c_row(F[i], parents)
    for j in range(0, len(C_row)):
        pen = penalty(parents + [j], F)
        C_row[j] = C_row[j] + pen
    return C_row


def smart_predict_penalty(F, m, fail_threshold, clones=[]):
    # Assumes square F
    if len(F) == 0:
        return [], [], [], -1

    # Add founder
    my_F = add_founder(F)
    # m = len(my_F)
    choices_stack = [[0]]
    chosen_parents = [0]

    valid_parents = white_list(list(range(m)), clones)
    order_validity = black_list(list(range(m)), clones)

    success = False
    while choices_stack:
        # print(choices_stack)
        chosen_parents = choices_stack.pop()
        i = len(chosen_parents) - 1
        if clones:
            if not valid_parent_order(order_validity, chosen_parents):
                continue
        if i == (m - 1):
            success = True
            break
        sc = scores(my_F, chosen_parents)
        next_choices = list((np.array(sc)).argsort())
        for next_choice in next_choices:
            if next_choice <= i and next_choice in valid_parents[i+1]:
                temp = chosen_parents[:]
                temp.append(next_choice)
                choices_stack.append(temp)
    if not success:
        return [], my_F, [], -1

    C = []
    for i in range(0, m):
        C.append(c_row(my_F[i], chosen_parents))

    p = one
    for i in range(1, m):
        p = p * C[i-1][chosen_parents[i]]
    # choices_stack.append(chosen_parents)
    return C, my_F, chosen_parents, p


def slow_predict(F, fail_threshold):
    # Assumes square F

    if len(F) == 0:
        return [], [], [], -1

    # Add founder
    my_F = add_founder(F)
    G = ancestry_graph(my_F)
    trees = enum_spanning_trees(G)
    max_p = -1
    max_tree = []
    max_C = []
    i = 0
    for tree in trees:
        if i % 1000000 == 0:
            stdout.write(" i2 = " + str(i/1000000) + " M; " + str(datetime.now()) + "\n")
            stdout.flush()
        i += 1
        C, p = get_c(my_F, tree, fail_threshold)
        if p > max_p:
            max_p = p
            max_tree = tree
            max_C = C
    par, _ = tree_to_list(max_tree)
    return max_C, my_F, par, max_p


def slow_predict_no_founder(F):
    # Assumes square F

    if len(F) == 0:
        return [], [], [], -1

    my_F = deepcopy(F)
    G = ancestry_graph(my_F)
    trees = enum_spanning_trees(G)
    max_p = -1
    max_tree = []
    max_C = []
    # stdout.write(" " + str(len(trees)) + " trees\n")
    # stdout.flush()
    for tree in trees:
        C, p = get_c(my_F, tree, zero)
        if p > max_p:
            max_p = p
            max_tree = tree
            max_C = C
    par, _ = tree_to_list(max_tree)
    return max_C, my_F, par, max_p


def extend_matrix(vafs):
    out = []
    n = len(vafs)
    vafs = list(vafs)
    for i in range(0, n):
        out.append(vafs[:(i+1)] + [zero]*(n-i-1))
    return out


def best_permutation(indices, vafs):
    out = []
    max_p = -1
    mapping = list(zip(indices, vafs))
    for permutation in permutations(mapping):
        [cur_indices, cur_vafs] = list(zip(*permutation))
        F = extend_matrix(cur_vafs)
        _, _, par, p = slow_predict_no_founder(F)
        if p > max_p:
            max_p = p
            out = cur_indices
    return out, max_p


def valid_permutations(indices, vafs, fail_threshold):
    perms = []
    ps = []
    mapping = list(zip(indices, vafs))
    stdout.write(str(len(list(permutations(mapping)))) + " permutations\n")
    stdout.flush()
    i = 0
    for permutation in permutations(mapping):
        if i % 10000 == 0:
            stdout.write(" Perm " + str(i+1) + "\n")
            stdout.flush()
        i += 1
        [cur_indices, cur_vafs] = list(zip(*permutation))
        F = extend_matrix(cur_vafs)
        _, _, par, p = slow_predict_no_founder(F)
        if p > fail_threshold:
            ps.append(p)
        else:
            ps.append(fail_threshold)
        perms.append(cur_indices)
    return perms, ps


def non_square_predict_random(F, variants):
    sf = square_forms(F, variants)
    F_iter, variants_iter = next(sf)
    my_F = add_founder(F_iter)
    parents_out = [0, 0]
    for i in range(1, len(variants_iter)):
        cur_parent = np.random.choice(list(range(0, i+1)))
        parents_out.append(cur_parent)
    p_out = get_p(my_F, parents_out, -100000)
    return variants_iter, parents_out, p_out


def non_square_predict_slow(F, variants, predictor, fail_threshold):
    parents_out = []
    variants_out = []
    p_out = -1
    for F_iter, variants_iter in square_forms(F, variants):
        c_out, _, parents, p = predictor(F_iter, fail_threshold)
        if p > p_out:
            p_out = p
            parents_out = parents
            variants_out = variants_iter

    if p_out == -1:
        return [], [], -1, []

    return variants_out, parents_out, p_out, c_out


def non_square_predict_smart(F, variants, smart_predict_algo, fail_threshold=FAIL_THRESHOLD, clones=[]):
    F_iter, variants_iter = rearrange(F, variants, True)

    step_structure = get_step_structure(F_iter)
    choices_stack = [(0, [0], [])]
    (cur_score, cur_T, cur_P) = (0, [0], [])

    success = False
    while choices_stack:
        (cur_score, cur_T, cur_P) = choices_stack.pop()
        # print(str(cur_T) + "\t" + str(cur_P))
        cur_m = 1
        for item in cur_P:
            cur_m += len(item)
        i = len(cur_P)
        if i == len(step_structure):
            success = True
            break
        cur_step = step_structure[i]
        temp = []
        j = 0
        for permutation in permutations(cur_step):
            # if j % 10000 == 0:
            #     stdout.write(" Perm " + str(j + 1) + "\n")
            #     stdout.flush()
            # j += 1
            perm_ss = deepcopy(cur_P)
            perm_ss.append(list(permutation))
            perm_F = squarify(F_iter, perm_ss)
            if len(perm_ss) == len(step_structure):
                new_clones = translate_clones(reorder(variants_iter, chain.from_iterable(perm_ss)), clones)
            else:
                new_clones = []
            _, _, perm_T, perm_score = smart_predict_algo(perm_F, cur_m + len(permutation), fail_threshold, new_clones)
            if perm_score != -1:
                temp.append((perm_score, perm_T, perm_ss))
        if temp:
            temp_sorted = sorted(temp, key=lambda x: x[0])
            choices_stack = choices_stack + temp_sorted

    if success:
        variants_out = reorder(variants_iter, chain.from_iterable(cur_P))
        F_temp = squarify(F_iter, cur_P)
        F_out = add_founder(F_temp)
        C_out = []
        for i in range(0, len(F_out)):
            C_out.append(c_row(F_out[i], cur_T))
        return variants_out, cur_T, cur_score, C_out
    else:
        return [], [], 0, []


def non_square_predict_fast(F, variants, predictor, fail_threshold, clones=[]):
    F_iter, variants_iter = rearrange(F, variants, True)
    step_structure = get_step_structure(F_iter)
    temp1 = []
    temp2 = []
    for i in range(0, len(step_structure)):
        indices = step_structure[i]
        if len(indices) > 1:
            perms, ps = valid_permutations(indices, F_iter[i][indices[0]:indices[-1] + 1], zero)
            if not perms:    # No valid solution!
                return [], [], -1
        else:
            perms = [indices]
            ps = [one]
        temp1.append(perms)
        temp2.append(ps)
    prod1 = list(product(*temp1))
    prod2 = list(map(np.prod, product(*temp2)))
    temp3 = zip(prod1, prod2)
    temp4 = sorted(temp3, key=lambda x: x[1], reverse=True)
    options, _ = zip(*temp4)

    stdout.write("Number of options to explore: " + str(len(options))+"\n")
    stdout.flush()

    i = 0
    for step_structure in options:
        if i % 10000 == 0:
            stdout.write(" Option " + str(i+1) + "\n")
            stdout.flush()
        i += 1
        F_sq = squarify(F_iter, step_structure)
        new_clones = translate_clones(reorder(variants_iter, chain.from_iterable(step_structure)), clones)
        c_out, _, parents_out, p_out = predictor(F_sq, len(F_sq)+1, fail_threshold, new_clones)
        if p_out == -1:
            continue
        variants_out = reorder(variants_iter, chain.from_iterable(step_structure))
        return variants_out, parents_out, p_out, c_out

    return [], [], 0, []


def non_square_predict(F, variants, algo, fail_threshold=FAIL_THRESHOLD, clones=[]):
    if algo == 0:
        return non_square_predict_random(F, variants)
    elif algo == 1:
        return non_square_predict_slow(F, variants, slow_predict, fail_threshold)
    elif algo == 2:
        return non_square_predict_slow(F, variants, smart_predict_original, fail_threshold)
    # elif algo == 3:
    #     return non_square_predict_fast(F, variants, smart_predict_original, fail_threshold, clones)
    elif algo == 3:
        return non_square_predict_smart(F, variants, smart_predict_original, fail_threshold, clones)
    elif algo == 4:
        return non_square_predict_smart(F, variants, smart_predict_penalty, fail_threshold, clones)
    else:
        exit("Invalid parameter for algo.")


def predict_batch(batch, out_file1, out_file2, out_file3, algo, truths=[], clone_count=0):
    for i in range(0, len(batch)):
        F = batch[i]
        clones = []
        if truths:
            truth = truths[i]
            clones = sample_clones(truth, clone_count)
            write_clones(clones, out_file3)
        stdout.write("Sim " + str(i) + "; " + str(datetime.now()) + "\n")
        stdout.flush()
        if not F:
            out_file1.write("\n")
            out_file2.write("0\t-1")
        else:
            var = list(range(1, len(F[0])+1))
            start = process_time()
            var1, parents1, p1, c1 = non_square_predict(F, var, algo, clones=clones)
            end = process_time()
            parents1 = remap_parents(var1, parents1)
            out_file1.write(" ".join(map(str, parents1)) + "\n")
            out_file2.write(str(end-start) + "\t" + str(p1) + "\n")


def get_match_histo(results1, results2):
    if len(results1) != len(results2) or len(results1) == 0:
        return []
    match_counts = []
    for i in range(0, len(results1)):
        r1 = results1[i][2:]
        r2 = results2[i][2:]
        if len(r1) != len(r2):
            match_counts.append(0)
        else:
            match_counts.append(len([j for j, k in zip(r1, r2) if j == k]))
    match_histo = [0] * (len(results1[0])-1)
    for c in match_counts:
        match_histo[c] += 1
    return match_histo


def get_match(results1, results2, clone_count=0):
    if len(results1) != len(results2) or len(results1) == 0:
        return []
    match_percs = []
    for i in range(0, len(results1)):
        r1 = results1[i][1:]
        r2 = results2[i][1:]
        if len(r1) != len(r2):
            match_percs.append(0.0)
        else:
            total = len(r1)
            match_count = len([j for j, k in zip(r1, r2) if j == k])
            match_percs.append((match_count + 0.0 - clone_count) / (total - clone_count))
            # match_percs.append(match_count)
    return match_percs


def get_clonal_match(results1, results2, clone_sets):
    if len(results1) != len(results2) or len(results1) == 0:
        return []
    match_percs = []
    for i in range(0, len(results1)):
        c1 = parents_to_clones(results1[i])
        c2 = parents_to_clones(results2[i])
        if clone_sets:
            clone_set = set(map(frozenset, clone_sets[i]))
            clone_set.add(frozenset([0]))
        else:
            clone_set = {frozenset([0])}
        c1 = c1.difference(clone_set)
        c2 = c2.difference(clone_set)
        match_count = len(c1.intersection(c2))
        total = len(c1)
        match_percs.append((match_count + 0.0) / total)
    return match_percs


def evaluate(truths, predictionss, clone_setss, out_f):
    out_f_file = open(out_f, "w")
    out = []
    for i in range(0, len(predictionss)):
        predictions = predictionss[i]
        if clone_setss:
            clone_sets = clone_setss[i]
        else:
            clone_sets = []
        out.append(get_clonal_match(truths, predictions, clone_sets))
    out = list(map(list, zip(*out)))
    for blah in out:
        out_f_file.write(" ".join(map(str, blah)) + "\n")


def get_clonal_mismatch(results1, results2, clone_sets):
    if len(results1) != len(results2) or len(results1) == 0:
        return []
    match_percs = []
    for i in range(0, len(results1)):
        c1 = parents_to_clones(results1[i])
        c2 = parents_to_clones(results2[i])
        if clone_sets:
            clone_set = set(map(frozenset, clone_sets[i]))
            clone_set.add(frozenset([0]))
        else:
            clone_set = {frozenset([0])}
        c1 = c1.difference(clone_set)
        c2 = c2.difference(clone_set)
        mismatch_count = len(c2.difference(c1))
        total = len(c1)
        match_percs.append(mismatch_count)
    return match_percs


def evaluate_mismatch(truths, predictionss, clone_setss, out_f):
    out_f_file = open(out_f, "w")
    out = []
    for i in range(0, len(predictionss)):
        predictions = predictionss[i]
        if clone_setss:
            clone_sets = clone_setss[i]
        else:
            clone_sets = []
        out.append(get_clonal_mismatch(truths, predictions, clone_sets))
    out = list(map(list, zip(*out)))
    for blah in out:
        out_f_file.write(" ".join(map(str, blah)) + "\n")


def batch_predict_test(batch, out_prefix, algo, truths=[], clone_count=0):
    stdout.write("\n" + str(datetime.now()) + "\n")
    stdout.write("Algo " + str(algo) + "\n")
    stdout.flush()
    file1 = open(out_prefix + "." + str(algo) + ".pred", "w")
    file2 = open(out_prefix + "." + str(algo) + ".time", "w")
    file3 = open(out_prefix + "." + str(algo) + ".clones", "w")
    predict_batch(batch, file1, file2, file3, algo, truths, clone_count)
    file1.close()
    file2.close()
    file3.close()
