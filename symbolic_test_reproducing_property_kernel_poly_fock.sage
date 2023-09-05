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


def kernel_poly_Fock(m, w, z, alpha):
    n = len(z)
    CF = parent(z[0])
    factor0 = exp(alpha * my_inner_product(w, z))
    diff0 = my_diff(w, z)
    factor1 = my_laguerre(m - 1, n, alpha * my_norm2(diff0))
    return CF(factor0 * factor1)


def my_power(a, m):
    F = parent(a)
    result = F(1)
    for k in range(m):
        result = F(result * a)
    return result


def monomial(p, q, z):
    n = len(z)
    F = parent(z[0])
    result = F(1)
    for j in range(n):
        result *= z[j] ** p[j]
        result *= conjugate(z[j]) ** q[j]
    return result


def myvars(letter, n):
    var_names = [letter + str(j) for j in range(n)]
    return var(var_names)


def symbolic_test_reproducing_property_kernel_poly_Fock_monomial(n, m, p, q):
    xs = myvars('x', n)
    ys = myvars('y', n)
    us = myvars('u', n)
    vs = myvars('v', n)
    alpha = var('a')
    for j in range(n):
        assume(xs[j] > 0)
        assume(ys[j] > 0)
        assume(us[j] > 0)
        assume(vs[j] > 0)
    assume(alpha > 0)
    zs = [xs[j] + I * ys[j] for j in range(n)]
    ws = [us[j] + I * vs[j] for j in range(n)]
    monz = monomial(p, q, zs).full_simplify()
    monw = monomial(p, q, ws).full_simplify()
    kernel = kernel_poly_Fock(m, ws, zs, alpha)
    weight = ((alpha / pi) ** n) * exp(- alpha * my_norm2(ws))
    f0 = weight * monw * conjugate(kernel)
    f = f0
    for j in range(n):
        f1 = f
        f2 = integral(f1, us[j], -Infinity, Infinity)
        f = integral(f2, vs[j], -Infinity, Infinity)
    f = f.full_simplify()
    diff0 = f - monz
    diff1 = diff0.full_simplify()
    r = diff1.is_zero()
    print('symbolic_test_repr_kernel_poly_Fock_monomial, n = %d, m = %d, p = %s, q = %s' % (n, m, str(p), str(q))) 
    #print('f0 = ' + str(f0))
    print('integral = ' + str(f))
    print('monomial = ' + str(monz))
    print('result = ' + str(r) + '\n')
    return r


def big_symbolic_test_reproducing_property_kernel_poly_Fock_monomial(n, m, psum_max):
    ps = lists_with_bounded_sum(n, psum_max)
    qs = lists_with_bounded_sum(n, m - 1)
    samples = [(p, q) for p in ps for q in qs]
    big_result = True
    for (p, q) in samples:
        r = symbolic_test_reproducing_property_kernel_poly_Fock_monomial(n, m, p, q)
        big_result = big_result and r
    return big_result


def huge_symbolic_test_reproducing_property_kernel_poly_Fock_monomial(nmax, mmax, psum_max):
    t0 = time.time()
    ns = list(range(1, nmax + 1))
    ms = list(range(1, mmax + 1))
    samples = [(n, m) for n in ns for m in ms]
    print('number of pairs (n, m):' + str(len(samples)))
    big_result = True
    for (n, m) in samples:
        r = big_symbolic_test_reproducing_property_kernel_poly_Fock_monomial(n, m, psum_max)
        big_result = big_result and r
    t1 = time.time()
    print('total time = ' + str(t1 - t0))
    return big_result


#print(symbolic_test_reproducing_property_kernel_poly_Fock_monomial(1, 2, [1], [1]))
#print(big_symbolic_test_reproducing_property_kernel_poly_Fock_monomial(2, 2, 2))
print(huge_symbolic_test_reproducing_property_kernel_poly_Fock_monomial(3, 3, 5))

