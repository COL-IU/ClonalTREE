import numpy as np
from myutils import *
from math import *
from predict import *
from time import *
import sys


def split_one(n):
    out = list(map(Decimal, np.random.rand(n)))
    s = sum(out)
    for i in range(n):
        out[i] = out[i]/s
    return out


def simple_simulate_square(n):
    parents = [0] * n
    C = [[0] * n]
    C[0][0] = 1
    indices = list(range(n))
    likelihood = one
    for i in range(1, n):
        parent = np.random.choice(indices, p=C[i-1])
        parents[i] = parent
        likelihood = likelihood * C[i-1][parent]
        new_row = split_one(i+1) + ([zero] * (n - i - 1))
        C.append(new_row)
    return C, parents, likelihood


def simple_simulator(params):
    [n, p] = params
    while True:
        C, parents, likelihood = simple_simulate_square(n)
        n_time_points = np.random.binomial(n, p)
        chosen_time_points = list(np.random.choice(range(n-1), replace=False, size=n_time_points-1))
        chosen_time_points.append(n-1)
        C_chosen = []
        for i in sorted(chosen_time_points):
            C_chosen.append(C[i])
        C_new, parents_new = remove_zero_clones(C_chosen, parents)
        if not parents_new:
            continue
        T = list_to_tree(parents_new)
        F = get_f(C_new, T)
        for row in F:
            del row[0]

        var = list(range(1, len(parents_new)))
        F, var_r = rearrange(F, var, True)

        if not all([v == 0 for v in parents_new]):
            step_structure = get_step_structure(F)
            if max(map(len, step_structure)) <= 6:
                yield F, remap_parents(var_r, parents_new), var, likelihood


def simple_simulator_fixed(params):
    [n, p] = params
    while True:
        C, parents, likelihood = simple_simulate_square(n)
        n_time_points = p
        chosen_time_points = list(np.random.choice(range(n-1), replace=False, size=n_time_points-1))
        chosen_time_points.append(n-1)
        C_chosen = []
        for i in sorted(chosen_time_points):
            C_chosen.append(C[i])
        C_new, parents_new = remove_zero_clones(C_chosen, parents)
        if not parents_new:
            continue
        T = list_to_tree(parents_new)
        F = get_f(C_new, T)
        for row in F:
            del row[0]

        var = list(range(1, len(parents_new)))
        F, var_r = rearrange(F, var, True)

        if not all([v == 0 for v in parents_new]):
            step_structure = get_step_structure(F)
            if max(map(len, step_structure)) <= 6:
                yield F, remap_parents(var_r, parents_new), var, likelihood


def clonal_composition(ts):
    C = []
    for row in ts:
        dist = calc_dist(row)
        C.append(dist)
    return C


def bottleneck(C, b):
    dist = calc_dist(C)
    out = list(map(lambda x: Decimal(x.item()), list(np.random.multinomial(b, dist))))
    return out


def spawn(C, mu, d, parents, fitnesses, fitness_low, fitness_high):
    C_out = C[:]
    parents_out = parents[:]
    fitnesses_out = fitnesses[:]
    dist = calc_dist(C_out)
    d_res = d - floor(d)
    nspawns = np.random.binomial(floor(d), mu)
    newparents = np.random.choice(range(len(dist)), size=nspawns, p=dist)
    for parent in newparents:
        if C_out[parent] > 1:
            C_out[parent] -= 1
            C_out.append(one)
            parents_out.append(parent)
            fitnesses_out.append(Decimal(np.random.uniform(fitness_low, fitness_high)).quantize(Decimal('.1')))
    return C_out, d_res, parents_out, fitnesses_out


def simulate(params):
    [n, parents, fitnesses, lamb, mu, T_max, b, cc, fitness_low, fitness_high] = params
    parents_out = parents[:]
    fitnesses_out = fitnesses[:]
    T = 1
    t = 1
    timeseries = [n]
    N = [sum(n)]
    d_res = 0
    while T < T_max:
        nclones = len(fitnesses_out)
        # print("T"+str(T) + "\t" + str(nclones))
        C = []
        for i in range(0, nclones):
            C.append(n[i] * (2 ** (lamb * fitnesses_out[i] * t)))
        timeseries.append(C)
        N.append(sum(C))
        d = N[-1] - N[-2] + d_res
        t += 1
        if d > 1:
            (C, d_res, parents_out, fitnesses_out) = spawn(C, mu, d, parents_out, fitnesses_out, fitness_low,
                                                           fitness_high)
            n = C[:]
            timeseries[-1] = C
            t = 1
        if N[-1] > cc:
            n = bottleneck(C, b)
            timeseries[-1] = n
            N[-1] = sum(n)
            t = 1
        T += 1
    return timeseries, parents_out, fitnesses_out


def fill_timeseries(timeseries, nclones):
    ts = timeseries[:]
    out = []
    for t in range(0, len(ts)):
        row = []
        C = ts[t]
        m = len(C)
        for i in range(0, nclones):
            if i < m:
                c = C[i]
            else:
                c = zero
            row.append(c)
        out.append(row)
    return out


def sample_timeseries(timeseries, cutoff, period):
    ts_sample = []
    for t in range(0, len(timeseries)):
        if t < cutoff or t % period != 0:
            continue
        else:
            ts_sample.append(timeseries[t])
    return ts_sample


def write_timeseries(timeseries_file, parents_file, timeseries, parents, fitnesses, cutoff, period):
    n = len(timeseries)
    m = len(timeseries[0])
    timeseries_file.write("Generation,Identity,Population,Fitness\n")
    for t in range(0, n):
        if t < cutoff or t % period != 0:
            continue
        for i in range(0, m):
            timeseries_file.write(str(t) + ",c" + str(i) + "," + "{0:.4f}".format(timeseries[t][i])
                                  + "," + "{0:.1f}".format(fitnesses[i]) + "\n")
    parents_file.write("Parent,Identity\n")
    for i in range(1, m):
        parents_file.write("c" + str(parents[i]) + ",c" + str(i) + "\n")


def remove_zero_clones(ts_sample, parents):
    ts_sample_in = list(map(list, zip(*ts_sample)))
    out = []
    mapping = [0] * len(parents)
    new_parents = []
    removed = []
    j = 0
    for i in range(0, len(ts_sample_in)):
        row = ts_sample_in[i]
        if not all([v == zero for v in row]):
            out.append(row)
            mapping[i] = j
            parent = parents[i]
            while parent in removed:
                if parent == 0:
                    parent = -1
                    break
                parent = parents[parent]
            new_parents.append(mapping[parent])
            j += 1
        else:
            removed.append(i)
    out = list(map(list, zip(*out)))
    return out, new_parents


def timeseries_to_f(timeseries, parents, params):
    [cutoff, period] = params
    filled = fill_timeseries(timeseries, len(parents))
    sample = sample_timeseries(filled, cutoff, period)
    timeseries_new, parents_new = remove_zero_clones(sample, parents)
    C = clonal_composition(timeseries_new)
    T = list_to_tree(parents_new)
    F = get_f(C, T)
    likelihood = one
    for i in range(1, len(T)):
        likelihood = likelihood * C[i - 1][parents_new[i]]
    for row in F:
        del row[0]
    # del parents_new[0]
    return F, parents_new, likelihood


def exponential_simulator(params):
    while True:
        (timeseries, parents, fitnesses) = simulate(params[0:10])
        F, parents, likelihood = timeseries_to_f(timeseries, parents, params[10:12])
        var = list(range(1, len(parents)))
        F, var_r = rearrange(F, var, True)

        if not all([v == 0 for v in parents]):
            step_structure = get_step_structure(F)
            if max(map(len, step_structure)) <= 6 and len(parents) <= 10:
                yield F, remap_parents(var_r, parents), var, likelihood


def evaluate_non_square(simulator, params, N, t, predictor_pairs):
    fails = [0.0] * len(predictor_pairs)
    tps = [0.0] * len(predictor_pairs)
    close_calls = [0.0] * len(predictor_pairs)
    times = [0.0] * len(predictor_pairs)

    sim = simulator(params)
    for i in range(N):
        print(i)
        (F, parents, var) = sim.__next__()

        print_dm(F)
        print(parents)

        for j in range(0, len(predictor_pairs)):
            (inner_predictor, outer_predictor) = predictor_pairs[j]
            start = time()
            var1, parents1, p1 = non_square_predict_fast(F, var, inner_predictor, outer_predictor)
            end = time()
            parents1 = remap_parents(var1, parents1)
            num_matches = sum([v1 == v2 for (v1, v2) in zip(parents1[1:], parents[1:])])
            freq_matches = (num_matches + 0.0) / (len(parents) - 1)
            print(parents1)

            if not parents1:
                fails[j] += 1
            if freq_matches == 1:
                tps[j] += 1
            if freq_matches > t:
                close_calls[j] += 1
            times[j] += (end - start)

    for i in range(0, len(times)):
        fails[i] = fails[i]/N
        tps[i] = tps[i]/N
        close_calls[i] = close_calls[i]/N
        times[i] = times[i]/N
    return fails, tps, close_calls, times


def write_simulations(simulator, params, N, f_out, parents_out, like_out):
    f_out_file = open(f_out, "w")
    parents_out_file = open(parents_out, "w")
    like_out_file = open(like_out, "w")
    sim = simulator(params)
    for i in range(N):
        (F, parents, var, likelihood) = sim.__next__()
        f_out_file.write("#\n")
        write_dm(F, f_out_file)
        parents_out_file.write(" ".join(map(str, parents)) + "\n")
        like_out_file.write(str(likelihood) + "\n")
    f_out_file.close()
    parents_out_file.close()
    like_out_file.close()


def read_simulations(f_in, parents_in):
    f_in_file = open(f_in)
    parents_in_file = open(parents_in)
    f_out = []
    parents_out = []
    lines = f_in_file.readlines()
    if lines[0][0] != '#':
        sys.exit("Invalid input file!")
    curr_f = []
    for line in lines:
        if len(line) > 0:
            if line[0] == '#':
                if curr_f:
                    f_out.append(curr_f)
                    curr_f = []
            else:
                curr_f.append(list(map(Decimal, line.split(None))))
    if curr_f:
        f_out.append(curr_f)
    for line in parents_in_file:
        parents_out.append(list(map(int, line.split(None))))
    f_in_file.close()
    parents_in_file.close()
    return f_out, parents_out


def evaluate_batch(batch_prefix, algos):
    batch, truths = read_simulations(batch_prefix + ".sims", batch_prefix + ".truth")
    results = []
    for i in algos:
        results.append(read_results(batch_prefix + "." + str(i) + ".pred"))
    evaluate(truths, results, [], batch_prefix + ".histo")
    likess = []
    timess = []
    f = open(batch_prefix + ".like")
    true_likes = list(map(float, f.readlines()))
    f.close()
    for i in algos:
        f = open(batch_prefix + "." + str(i) + ".time")
        [times, likes] = list(map(list, zip(*(list(map(lambda x: list(map(float, x.split(None))), f.readlines()))))))
        times = list(map(lambda x: x if x != 0 else 0.0001, times))
        likes = list(map(lambda x, y: np.log10(abs(x)/y), likes, true_likes))
        times = list(map(np.log10, times))
        likess.append(likes)
        timess.append(times)
        f.close()
    likess = list(map(list, zip(*likess)))
    timess = list(map(list, zip(*timess)))
    f = open(batch_prefix + ".llike", "w")
    for blah in likess:
        f.write("\t".join(list(map(str, blah))) + "\n")
    f.close()
    f = open(batch_prefix + ".ltime", "w")
    for blah in timess:
        f.write("\t".join(list(map(str, blah))) + "\n")
    f.close()
