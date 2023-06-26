# This is a test of the formula
# L_m^{(n)}(t_0+\ldots+t_{n-1})
# = \sum_{k_0+\ldots+k_{n-1}\le m} \prod_{j=0}^{n-1} L_{k_j}(t_j),
# where L denotes the generalized Laguerre(-Sonin) polynomial.


def lists_with_bounded_sum(n, s):
    # lists [k_0, k_1, ..., k_{n-1}] with k_0 + k_1 + ... + k_{n-1} \le s
    r = [[]]
    for j in range(n):
        rprev = r
        r = []
        for x in rprev:
            sprev = sum(x)
            for y in range(s - sprev + 1):
                z = x + [y]
                r += [z]
    return r


def lists_with_given_sum(n, s):
    # lists [k_0, k_1, ..., k_{n-1}] with k_0 + k_1 + ... + k_{n-1} = s
    r = [[]] if (n == 0 and s == 0) else []
    if n > 0:
        rprev = lists_with_bounded_sum(n - 1, s)
        for x in rprev:
            y = s - sum(x)
            z = x + [y]
            r += [z]
    return r


def product_laguerre(k, t):
    # k is a multi-index
    n = len(t)
    PR = parent(t[0])
    result = PR.one()
    for j in range(n):
        result *= laguerre(k[j], t[j])
    return result


def sum_products_laguerre(m, t):
    n = len(t)
    PR = parent(t[0])
    ks = lists_with_bounded_sum(n, m)
    result = PR.zero()
    for k in ks:
        result += product_laguerre(k, t)
    return result


def num_summands(n, m):
    result0 = len(lists_with_bounded_sum(n, m))
    result1 = binomial(n + m, n)
    return (result0, result1)


def test_num_summands(n, m):
    (result0, result1) = num_summands(n, m)
    return result0 == result1


def big_test_num_summands(nmax, mmax):
    result = True
    for n in range(1, nmax + 1):
        for m in range(mmax + 1):
            r0 = test_num_summands(n, m)
            result = result and r0
    return result


def symbolic_test_laguerre_addition_formula(n, m, verb):
    varnames = ['t' + str(k) for k in range(n)]
    PR = PolynomialRing(QQ, varnames)
    t = list(PR.gens())
    result0 = PR(gen_laguerre(m, n, sum(t)))
    result1 = PR(sum_products_laguerre(m, t))
    if verb:
        print('result0 = ', result0)
        print('result1 = ', result1)
    return result0 == result1


def random_symbolic_test_laguerre_addition_formula():
    n = ZZ.random_element(1, 6)
    m = ZZ.random_element(0, 6)
    print('n = ', n)
    print('m = ', m)
    return symbolic_test_laguerre_addition_formula(n, m, True)


def big_symbolic_test_laguerre_addition_formula(nmax, mmax):
    result = True
    for n in range(1, nmax + 1):
        for m in range(mmax + 1):
            r0 = symbolic_test_laguerre_addition_formula(n, m, False)
            print(n, m, r0)
            result = result and r0
    return result


#print(symbolic_test_laguerre_addition_formula(2, 2, True))
#print(random_symbolic_test_laguerre_addition_formula())
#print(num_summands(2, 6))

print(big_test_num_summands(5, 5))
print(big_symbolic_test_laguerre_addition_formula(5, 5))

