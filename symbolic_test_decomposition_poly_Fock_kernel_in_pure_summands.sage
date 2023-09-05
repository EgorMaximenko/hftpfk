import time


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


def my_binomial(a, m):
    F = parent(a)
    result = F(1)
    for k in range(m):
        result *= F(a - k) / F(k + 1)
    return F(result)


def my_laguerre(n, al, t):
    F = parent(t)
    result = F(0)
    for k in range(n + 1):
        factor0 = (-1) ** k
        factor1 = my_binomial(n + al, n - k)
        factor2 = F(t ** k)
        denom = factorial(k)
        result += F(factor0 * factor1 * factor2 / denom)
    return result


def my_inner_product(w, z):
    n = len(w)
    return sum(w[j] * conjugate(z[j]) for j in range(n))


def my_norm2(z):
    n = len(z)
    return sum(z[j] * conjugate(z[j]) for j in range(n))


def my_diff(w, z):
    n = len(w)
    return [w[j] - z[j] for j in range(n)]


def kernel_poly_fock(m, w, z, alpha):
    n = len(z)
    CF = parent(z[0])
    factor0 = exp(alpha * my_inner_product(w, z))
    diff0 = my_diff(w, z)
    factor1 = my_laguerre(m - 1, n, alpha * my_norm2(diff0))
    result = CF(factor0 * factor1)
    result = result.full_simplify()
    return result


def kernel_pure_poly_fock_one_dim(m, w, z, alpha):
    CF = parent(z)
    w_minus_z = w - z
    factor0 = exp(alpha * w * conjugate(z))
    factor1 = my_laguerre(m - 1, 0, alpha * w_minus_z * conjugate(w_minus_z))
    return CF(factor0 * factor1)


def kernel_pure_poly_fock(beta, w, z, alpha):
    n = len(z)
    CF = parent(z[0])
    result = CF(1)
    for r in range(n):
        result *= kernel_pure_poly_fock_one_dim(beta[r], w[r], z[r], alpha)
    return result


def kernel_poly_fock_decomposed_in_pure_summands(m, w, z, alpha):
    n = len(z)
    CF = parent(z[0])
    ks = lists_with_bounded_sum(n, m - 1)
    result = CF(0)
    for k in ks:
        beta = [k[r] + 1 for r in range(n)]
        result += kernel_pure_poly_fock(beta, w, z, alpha)
    result = result.full_simplify()
    return result


def myvars(letter, n):
    var_names = [letter + str(j) for j in range(n)]
    return var(var_names)


def symbolic_test_decomposition_poly_fock_kernel_in_pure_summands(n, m):
    xs = myvars('x', n)
    ys = myvars('y', n)
    us = myvars('u', n)
    vs = myvars('v', n)
    alpha = var('a')
    z = [xs[j] + I * ys[j] for j in range(n)]
    w = [us[j] + I * vs[j] for j in range(n)]
    result0 = kernel_poly_fock(m, w, z, alpha)
    result1 = kernel_poly_fock_decomposed_in_pure_summands(m, w, z, alpha)
    result0_minus_result1 = (result0 - result1).full_simplify()
    r = result0_minus_result1.is_zero()
    print('symbolic_test_decomposition_poly_fock_kernel_in_pure_summands')
    print('n = ' + str(n) + ', m = ' + str(m))
    print('result0 = ' + str(result0))
    print('result1 = ' + str(result1))
    print(r, '\n')
    return r


def big_symbolic_test_decomposition_poly_fock_kernel_in_pure_summands(nmax, mmax):
    print('big_symbolic_test_decomposition_poly_fock_kernel_in_pure_summands')
    print('nmax = ' + str(nmax) + ', mmax = ' + str(mmax))
    ns = list(range(1, nmax + 1))
    ms = list(range(1, mmax + 1))
    samples = [(n, m) for n in ns for m in ms]
    print('number of samples = ' + str(len(samples)) + '\n')
    big_result = True
    for (n, m) in samples:
        r = symbolic_test_decomposition_poly_fock_kernel_in_pure_summands(n, m)
        big_result = big_result and r
    print('big_result = ', big_result)
    return big_result


#print(symbolic_test_decomposition_poly_fock_kernel_in_pure_summands(2, 2))


print(big_symbolic_test_decomposition_poly_fock_kernel_in_pure_summands(5, 5))

