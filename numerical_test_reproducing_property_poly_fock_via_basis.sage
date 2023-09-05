# This program is a part of tests for the paper
# "Horizontal Fourier transform of the reproducing kernel of the polyanalytic Fock space"
# by Erick Lee-Guzm\'an, Egor A. Maximenko, Gerardo Ramos-Vazquez, Armando S\'anchez-Nungaray.


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
    return sum(abs(z[j]) ** 2 for j in range(n))


def kernel_poly_Fock(al, m, w, z):
    n = len(z)
    CF = parent(z[0])
    factor0 = exp(al * my_inner_product(w, z))
    factor1 = my_laguerre(m - 1, n, al * my_norm2(w - z))
    return CF(factor0 * factor1)


def my_power(a, m):
    F = parent(a)
    result = F(1)
    for k in range(m):
        result = F(result * a)
    return result


def basic_polynomial_onedim(p, q, z):
    CF = parent(z)
    RF = RealField(CF.precision())
    result = 0
    zabs2 = RF(abs(z) * abs(z))
    if p >= q:
        factor0 = sqrt(RF(factorial(q) / factorial(p)))
        factor1 = (-1) ** q
        factor2 = CF(my_power(z, p - q))
        factor3 = RF(my_laguerre(q, p - q, zabs2))
    else:
        factor0 = sqrt(RF(factorial(p) / factorial(q)))
        factor1 = (-1) ** p
        zc = conjugate(z)
        factor2 = CF(my_power(zc, q - p))
        factor3 = RF(my_laguerre(p, q - p, zabs2))
    return CF(factor0 * factor1 * factor2 * factor3)


def basic_polynomial(al, p, q, z):
    n = len(p)
    return prod(basic_polynomial_onedim(p[j], q[j], sqrt(al) * z[j]) for j in range(n))


def fast_basic_polynomial(p, q, values_onedim_z):
    n = len(p)
    return prod(values_onedim_z[j][p[j], q[j]] for j in range(n))


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


def kernel_poly_Fock_via_basis(al, m, w, z, maxdeg):
    n = len(w)
    CF = parent(w[0])
    ps = lists_with_bounded_sum(n, maxdeg)
    qs = lists_with_bounded_sum(n, m - 1)
    pqs = list([(p, q) for p in ps for q in qs])
    print('number of basic polynomials = %d' % len(pqs))
    result = CF(0)
    for p, q in pqs:
        factorw = CF(basic_polynomial(al, p, q, w))
        factorz = CF(conjugate(basic_polynomial(al, p, q, z)))
        result += factorw * factorz
    return result


def fast_kernel_poly_Fock_via_basis(al, m, w, z, maxdeg):
    n = len(w)
    CF = parent(w[0])
    ps = lists_with_bounded_sum(n, maxdeg)
    qs = lists_with_bounded_sum(n, m - 1)
    pqs = list([(p, q) for p in ps for q in qs])
    print('number of basic polynomials = %d' % len(pqs))
    values_onedim_z = [matrix(CF, maxdeg + 1, m) for j in range(n)]
    values_onedim_w = [matrix(CF, maxdeg + 1, m) for j in range(n)]
    for d in range(n):
        for j in range(maxdeg + 1):
            for k in range(m):
                values_onedim_w[d][j, k] = CF(basic_polynomial_onedim(j, k, sqrt(al) * w[d]))
                values_onedim_z[d][j, k] = CF(basic_polynomial_onedim(j, k, sqrt(al) * z[d]))
    result = CF(0)
    for p, q in pqs:
        factorw = CF(fast_basic_polynomial(p, q, values_onedim_w))
        factorz = CF(fast_basic_polynomial(p, q, values_onedim_z))
        result += factorw * conjugate(factorz)
    return CF(result)


def slow_random_test_kernel_poly_Fock_via_basis(n, m, maxdeg, prec, verbose_level):
    if verbose_level >= 1:
        print('slow random test kernel for n = %d, m = %d, maxdeg = %d, prec = %d' % (n, m, maxdeg, prec))
    CF = ComplexField(prec)
    RF = RealField(prec)
    al = RF.random_element(1, 2)
    w0 = random_vector(CF, n)
    z0 = random_vector(CF, n)
    w = (1 / (norm(w0) + 2)) * w0
    z = (1 / (norm(z0) + 2)) * z0
    r0 = kernel_poly_Fock(al, m, w, z)
    r1 = kernel_poly_Fock_via_basis(al, m, w, z, maxdeg)
    if verbose_level >= 2:
        print('r0 = ', r0)
        print('r1 = ', r1)
    if verbose_level >= 1:
        print('error = ', N(abs(r0 - r1), 32))
    return abs(r0 - r1)


def fast_random_test_kernel_poly_Fock_via_basis(n, m, maxdeg, prec, verbose_level):
    if verbose_level >= 1:
        print('fast random test kernel for n = %d, m = %d, maxdeg = %d, prec = %d' % (n, m, maxdeg, prec))
    CF = ComplexField(prec)
    RF = RealField(prec)
    al = RF.random_element(1, 2)   
    w0 = random_vector(CF, n)
    z0 = random_vector(CF, n)
    w = (1 / (norm(w0) + 2)) * w0
    z = (1 / (norm(z0) + 2)) * z0
    r0 = kernel_poly_Fock(al, m, w, z)
    r1 = fast_kernel_poly_Fock_via_basis(al, m, w, z, maxdeg)
    if verbose_level >= 2:
        print('r0 = ', r0)
        print('r1 = ', r1)
    if verbose_level >= 1:
        print('error = ', N(abs(r0 - r1), 32))
    return abs(r0 - r1)


def big_random_test_kernel_poly_Fock_via_basis(nmax, mmax, maxdeg, nreps, prec):
    print('big random test kernel poly Fock via basis,')
    print('nmax = %d, mmax = %d, maxdeg = %d, prec = %d' % (nmax, mmax, maxdeg, prec))
    RF = RealField(prec)
    max_error = RF(0)
    for n in range(1, nmax + 1):
        for m in range(1, mmax + 1):
            for r in range(nreps):
                er = fast_random_test_kernel_poly_Fock_via_basis(n, m, maxdeg, prec, 0)
                max_error = max(max_error, er)
                print('n = %d, m = %d, er = %s' % (n, m, str(N(er, 32))))
    print('maximal error = ' + str(N(max_error, 32)))
    return max_error


def some_random_tests_kernel_poly_Fock_via_basis():
    print(slow_random_test_kernel_poly_Fock_via_basis(1, 4, 128, 128, 2), '\n')
    print(fast_random_test_kernel_poly_Fock_via_basis(1, 4, 128, 128, 2), '\n')
    print(slow_random_test_kernel_poly_Fock_via_basis(2, 4, 128, 128, 2), '\n')
    print(fast_random_test_kernel_poly_Fock_via_basis(2, 4, 128, 128, 2), '\n')
    print(slow_random_test_kernel_poly_Fock_via_basis(3, 4, 128, 128, 2), '\n')
    print(fast_random_test_kernel_poly_Fock_via_basis(3, 4, 128, 128, 2), '\n')


#some_random_tests_kernel_poly_Fock_via_basis()


big_random_test_kernel_poly_Fock_via_basis(3, 4, 128, 5, 128)

