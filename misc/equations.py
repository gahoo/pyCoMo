#! /usr/bin/env python
#coding=utf-8
#from mpmath import loggamma
#from math import exp
import math,mpmath
import pp
import sys


def logchoose(ni, ki):
    try:
        lgn1 = mpmath.loggamma(ni + 1)
        lgk1 = mpmath.loggamma(ki + 1)
        lgnk1 = mpmath.loggamma(ni - ki + 1)
    except ValueError:
        #print ni,ki
        raise ValueError
    return lgn1 - (lgnk1 + lgk1)


def gauss_hypergeom(X, n, m, N):
    """Returns the probability of drawing X successes of m marked items
     in n draws from a bin of N total items."""

    assert N >= m, 'Number of items %i must be larger than the number of marked items %i' % (N, m)
    assert m >= X, 'Number of marked items %i must be larger than the number of sucesses %i' % (m, X)
    assert n >= X, 'Number of draws %i must be larger than the number of sucesses %i' % (n, X)
    assert N >= n, 'Number of draws %i must be smaller than the total number of items %i' % (n, N)

    r1 = logchoose(m, X)
    try:
        r2 = logchoose(N - m, n - X)
    except ValueError:
        return 0
    r3 = logchoose(N, n)

    return math.exp(r1 + r2 - r3)


def hypergeo_cdf(X, n, m, N):
    '''
    用来计算共现的p
    '''

    assert N >= m, 'Number of items %i must be larger than the number of marked items %i' % (N, m)
    assert m >= X, 'Number of marked items %i must be larger than the number of sucesses %i' % (m, X)
    assert n >= X, 'Number of draws %i must be larger than the number of sucesses %i' % (n, X)
    assert N >= n, 'Number of draws %i must be smaller than the total number of items %i' % (n, N)
    assert N - m >= n - X, 'There are more failures %i than unmarked items %i' % (N - m, n - X)

    s = 0
    for i in xrange(X, min(m, n) + 1):
        s += max(gauss_hypergeom(i, n, m, N), 0.0)
    return min(max(s, 0.0), 1)



def hypergeo_cdf_PAN(X, n, m, N):
    '''
    用来计算富集的p
    N：数据库中基因总数
    m：数据库中特定motif数
    n：列表中基因总数
    x：列表中特定motif数
    '''
    assert N >= m, 'Number of items %i must be larger than the number of marked items %i' % (N, m)
    assert m >= X, 'Number of marked items %i must be larger than the number of sucesses %i' % (m, X)
    assert n >= X, 'Number of draws %i must be larger than the number of sucesses %i' % (n, X)
    assert N >= n, 'Number of draws %i must be smaller than the total number of items %i' % (n, N)
    assert N - m >= n - X, 'There are more failures %i than unmarked items %i' % (N - m, n - X)

    s = 0
    for i in xrange(X, n + 1):
        s += max(gauss_hypergeom(i, n, m, N), 0.0)
    return min(max(s, 0.0), 1)


def hypergeo_cdf_enrich(X, n, m, N):
    '''
    用来计算富集的p
    N：数据库中基因总数
    m：数据库中特定motif数
    n：列表中基因总数
    x：列表中特定motif数
    '''
    assert N >= m, 'Number of items %i must be larger than the number of marked items %i' % (N, m)
    assert m >= X, 'Number of marked items %i must be larger than the number of sucesses %i' % (m, X)
    assert n >= X, 'Number of draws %i must be larger than the number of sucesses %i' % (n, X)
    assert N >= n, 'Number of draws %i must be smaller than the total number of items %i' % (n, N)
    assert N - m >= n - X, 'There are more failures %i than unmarked items %i' % (N - m, n - X)

    s = 0
    for i in range(0, X + 1):
        s += max(gauss_hypergeom(i, n, m, N), 0.0)
    return min(max(1 - s, 0.0), 1)


def enrichmen_ratio(x, list_size, m, db_size):
    return (float(x) / list_size) / (float(m) / db_size)


def bonferroni_adjustment(p_values):
    """ Correct a list of p-values using the Bonferroni adjustment

    Return a list of corrected p-values; null values are ignored.
    cf. http://en.wikipedia.org/wiki/Bonferroni_correction
    """
    n = len(filter(lambda x: x != None, p_values)) + 0.0

    adjusted_p_values = []
    for p_value in p_values:
        if (p_value == None):
            adjusted_p_values.append(None)
        else:
            adjusted_p_values.append(min(p_value * n, 1))

    return adjusted_p_values


def holm_adjustment(p_values):
    """ Correct a list of p-values using the Holm-Bonferroni adjustment

    Return a list of corrected p-values; null values are ignored.
    cf. http://en.wikipedia.org/wiki/Holm-Bonferroni_method
    """
    # multiply p-values by a corrective factor, ignoring null entries
    n, c = len(filter(lambda x: x != None, p_values)), 0
    m = []

    adjusted_and_ranked_p_values = []
    for i, (i_, p_value) in enumerate(sorted(enumerate(p_values), lambda x, y: cmp(x[1], y[1]))):
        m.append(i_)
        if (p_value == None):
            adjusted_and_ranked_p_values.append(None)
        else:
            adjusted_and_ranked_p_values.append(min(p_value * (n - c), 1))
        c += 1

    # correct the p-values out of their proper order
    adjusted_p_values = [0 for c in range(i + 1)]
    for i, p_value in enumerate(adjusted_and_ranked_p_values):
        if(p_value != None):
            adjusted_and_ranked_p_values[i] = max(filter(lambda x: x != None, adjusted_and_ranked_p_values[:i + 1]))

        adjusted_p_values[m[i]] = adjusted_and_ranked_p_values[i]

    return adjusted_p_values


def ppCacl(job_server, inputs, equation, funs=(gauss_hypergeom, logchoose), packages=("math", "mpmath"), Progress=True):
    num_inputs = len(inputs) / 100 + 1
    #ppservers = ()
    #job_server = pp.Server(ppservers=ppservers)
    #inputs = [(10, 20, 30, 50)] * 1000
    jobs = [job_server.submit(equation, pars, funs, packages) for pars in inputs]
    #ps = [job() for job in jobs]
    ps = []
    for i, job in enumerate(jobs):
        if Progress:
            sys.stderr.write("Progress:%d%%\r" % (i / num_inputs))
        ps.append(job())
    return ps


def ssCacl(inputs, equation, Progress=True):
    num_inputs = len(inputs) / 100 + 1
    ps = []
    for i, (X, n, m, N) in enumerate(inputs):
        if Progress:
            sys.stderr.write("Progress:%d%%\r" % (i / num_inputs))
        ps.append(equation(X, n, m, N))
    return ps


if __name__ == '__main__':
    (X, n, m, N) = (10, 20, 30, 50)
    (X, n, m, N) = (18, 65, 3781, 41953)
    print hypergeo_cdf(X, n, m, N)
    print hypergeo_cdf_enrich(X, n, m, N)
    print hypergeo_cdf_PAN(X, n, m, N)
    '''
    ppservers = ()
    job_server = pp.Server(ppservers=ppservers)
    inputs = [(10, 20, 30, 50)] * 1000
    jobs = [((X, n, m, N), job_server.submit(hypergeo_cdf,(X, n, m, N, ), (gauss_hypergeom, logchoose), ("math","mpmath",))) for (X, n, m, N) in inputs]
    for (X, n, m, N), job in jobs:
        print "Sum of primes below", (X, n, m, N), "is", job()
    '''
    print hypergeo_cdf_pp([(10, 20, 30, 50)] * 500, hypergeo_cdf_enrich)
