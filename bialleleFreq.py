import sys

if len(sys.argv) < 3:
    sys.exit("Usage: bialleleFreq.py <printBases> <time_tag>")

f = open(sys.argv[1])
w = open(sys.argv[1]+".allele","w")
freqs = open(sys.argv[1]+".freqs","w")


for line in f:
    words = line.split(None)
    contig = words[0]
    pos = words[1]
    ref = words[2]
    depth = int(words[3])
    A = int(words[4])
    a = int(words[5])
    C = int(words[6])
    c = int(words[7])
    G = int(words[8])
    g = int(words[9])
    T = int(words[10])
    t = int(words[11])
    dict = {}
    dict['A'] = A+a+0.0
    dict['C'] = C+c+0.0
    dict['G'] = G+g+0.0
    dict['T'] = T+t+0.0

    ab = 0 if (A+a) == 0 else min((A+0.0)/(A+a), (a+0.0)/(A+a))
    cb = 0 if (C+c) == 0 else min((C+0.0)/(C+c), (c+0.0)/(C+c))
    gb = 0 if (G+g) == 0 else min((G+0.0)/(G+g), (g+0.0)/(G+g))
    tb = 0 if (T+t) == 0 else min((T+0.0)/(T+t), (t+0.0)/(T+t))

    if ab < 0.25:
        dict['A'] = 0
        A = 0
        a = 0
    if cb < 0.25:
        dict['C'] = 0
        C = 0
        c = 0
    if gb < 0.25:
        dict['G'] = 0
        G = 0
        g = 0
    if tb < 0.25:
        dict['T'] = 0
        T = 0
        t = 0

    calls = []
    if dict['A'] > 0:
        calls.append('A')
    if dict['C'] > 0:
        calls.append('C')
    if dict['G'] > 0:
        calls.append('G')
    if dict['T'] > 0:
        calls.append('T')


    nrefCalls = calls[:]
    if ref in nrefCalls:
        nrefCalls.remove(ref)
    if len(nrefCalls) == 0:
        continue
    maxvaf = 0.0
    maxallele = ''
    for nonRefAllele in nrefCalls:
        vaf = (dict[nonRefAllele]+0.0)/(dict[nonRefAllele] + dict[ref])
        if vaf > maxvaf:
            maxvaf = vaf
            maxallele = nonRefAllele
    nraf = maxvaf
    if nraf < 0.05 or (dict[nonRefAllele] + dict[ref]) < 10 or dict[nonRefAllele] < 6:
        continue
    nonRefAllele = maxallele
    raf = (dict[ref]+0.0)/(dict[nonRefAllele] + dict[ref])
    if raf < nraf:
        maf = raf
    else:
        maf = nraf

    freqs.write(sys.argv[2] + "\t" + contig + "_" + pos + "\t" + ref + "\t" + nonRefAllele + "\t" + "{0:.4f}".format(raf) + "\t" + "{0:.4f}".format(nraf) + "\t" + "{0:.4f}".format(maf) + "\n")

    callsString = calls[0]
    if len(calls) > 1:
        for i in range(1,len(calls)):
            callsString = callsString + "/" + calls[i]
    w.write(contig + "\t" + pos + "\t" + ref + "\t" + callsString + "\t" + str(dict[nonRefAllele] + dict[ref]) + "\t" + str(A) + "\t" + str(a) + "\t" + str(C) + "\t" + str(c) + "\t" + str(G) + "\t" + str(g) + "\t" + str(T) + "\t" + str(t) + "\n")

f.close()
w.close()
freqs.close()

