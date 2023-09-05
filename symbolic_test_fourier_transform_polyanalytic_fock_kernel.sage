# Here we test the formula "(Fourier \otimes Identity) K = L" from the paper
# Horizontal Fourier transform of the reproducing kernel of the polyanalytic Bargmann--Segal--Fock space


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


def myvars(n):
    unames = ['u' + str(j) for j in range(n)]
    ynames = ['y' + str(j) for j in range(n)]
    vnames = ['v' + str(j) for j in range(n)]
    xinames = ['xi' + str(j) for j in range(n)]
    allvars = var(unames + ynames + vnames + xinames)
    u = allvars[0 : n]
    y = allvars[n : 2 * n]
    v = allvars[2 * n : 3 * n]
    xi = allvars[3 * n : 4 * n]
    return u, y, v, xi


def mynorm(a):
    return sum([elem ** 2 for elem in a])


def mydotproduct(a, b):
    return sum([a[j] * b[j] for j in range(len(a))])


def sumvec(a, b):
    return list([a[j] + b[j] for j in range(len(a))])


def difvec(a, b):
    return list([a[j] - b[j] for j in range(len(a))])


def mykernel(m, u, y, v):
    # flattened poly-Fock kernel, space Hm from our article
    n = len(u)
    coef = (2 / pi) ** (n / 2)
    v_minus_y = difvec(v, y)
    v_plus_y = sumvec(v, y)
    u_norm2 = mynorm(u)
    v_minus_y_norm2 = mynorm(v_minus_y)
    u_product_v_plus_y = mydotproduct(u, v_plus_y)
    expr_power = - (u_norm2 + v_minus_y_norm2) / 2 - I * u_product_v_plus_y
    exp_factor = exp(expr_power)
    laguerre_arg = u_norm2 + v_minus_y_norm2
    laguerre_factor = gen_laguerre(m - 1, n, laguerre_arg)
    return coef * exp_factor * laguerre_factor


def fourier_mykernel(m, xi, u, y, v):
    n = len(u)
    K = mykernel(m, u, y, v)
    character_factor = exp(- I * mydotproduct(u, xi))
    result = K * character_factor
    for j in range(n):
        old_result = result
        result = integral(old_result, u[j], -Infinity, Infinity)
    result = ((2 * pi) ** (- n / 2)) * result
    return result


def hermite_function(n, t):
    coef = ((2 ** n) * factorial(n) * sqrt(pi)) ** (- 1 / 2)
    exp_factor = exp(- t ** 2 / 2)
    hermite_factor = hermite(n, t)
    return coef * exp_factor * hermite_factor


def q_factor(k, xi, v):
    n = len(xi)
    coef = 2 ** (n / 4)
    return coef * prod([hermite_function(k[j], (xi[j] + 2 * v[j]) / sqrt(2)) for j in range(n)])


def L_formula(m, xi, u, y, v):
    n = len(xi)
    ks = lists_with_bounded_sum(n, m - 1)
    result = 0
    for k in ks:
        result += q_factor(k, xi, y) * q_factor(k, xi, v)
    return result


def test_fourier_mykernel(n, m, verbose):
    u, y, v, xi = myvars(n)
    L0 = fourier_mykernel(m, xi, u, y, v)
    L0 = L0.full_simplify()
    L1 = L_formula(m, xi, u, y, v)
    L1 = L1.full_simplify()
    er = L0 - L1
    er = er.full_simplify()
    result = bool(er.is_zero())
    if verbose >= 1:
        print('n = %d, m = %d, result = %s' % (n, m, str(result)))
    if verbose >= 2:
        print('\nFourier transform of K:\n', L0)
        print('\nL by formula:\n', L1, '\n')
    return result


def small_test():
    nmax = 2
    mmax = 2
    verbose = 2
    samples = [(n, m) for n in range(1, nmax + 1) for m in range(1, mmax + 1)]
    result = True
    for n, m in samples:
        result = result and test_fourier_mykernel(n, m, verbose)
    print('global result = ' + str(result))
    return result


def big_test():
    nmax = 4
    mmax = 4
    print('big test with nmax = %d and mmax = %d\n' % (nmax, mmax))
    verbose = 1
    samples = [(n, m) for n in range(1, nmax + 1) for m in range(1, mmax + 1)]
    result = True
    for n, m in samples:
        result = result and test_fourier_mykernel(n, m, verbose)
    print('global result = ' + str(result))
    return result


test_fourier_mykernel(2, 2, 2)
print('')
print(big_test())

