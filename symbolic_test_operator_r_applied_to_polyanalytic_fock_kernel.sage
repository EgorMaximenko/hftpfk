# Here we test the formula for R K_{iy} from our paper
# Horizontal Fourier transform of the reproducing Fock kernel,
# by Erick Lee-Guzm\'an, Egor A. Maximenko, Gerardo Ramos-Vazquez, Armando S\'anchez-Nungaray


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


def myvars(letter, m):
    return var([letter + str(j) for j in range(m)])


def vars_aluvyxi(n):
    al = var('al')
    u = myvars('u', n)
    v = myvars('v', n)
    y = myvars('y', n)
    xi = myvars('xi', n)
    return al, u, v, y, xi


def my_binomial(a, m):
    F = parent(a)
    result = F(1)
    for k in range(m):
        result *= F(a - k) / F(k + 1)
    return F(result)


def laguerre_coef(n, param, k):
    factor0 = (-1) ** k
    factor1 = my_binomial(n + param, n - k)
    denom = factorial(k)
    return factor0 * factor1 / denom


def laguerre_pol(n, param, t):
    F = parent(t)
    result = F(0)
    for k in range(n + 1):
        coef = laguerre_coef(n, param, k)
        t_power = F(t ** k)
        term = F(coef * t_power)
        result += term
    return result


def hermite_coef(n, k):
    numer0 = factorial(n)
    numer1 = 2 ** (n - 2 * k)
    numer2 = (-1) ** k
    denom0 = factorial(k)
    denom1 = factorial(n - 2 * k)
    return (numer0 * numer1 * numer2) / (denom0 * denom1)


def hermite_pol(n, t):
    # physicist's Hermite polynomial
    F = parent(t)
    result = F(0)
    n2 = n // 2
    for k in range(n2 + 1):
        coef = hermite_coef(n, k)
        tpower = t ** (n - 2 * k)
        term = F(coef * tpower)
        result += term
    return result


def hermite_function(n, t):
    f0 = hermite_pol(n, t)
    f1 = exp(- t * t / 2)
    coef0 = pi ** (-1 / 4)
    coef1 = 2 ** (- n / 2)
    coef2 = 1 / sqrt(factorial(n))
    coef = coef0 * coef1 * coef2
    return coef * f0 * f1


def real_dot_product(a, b):
    return sum([a[j] * b[j] for j in range(len(a))])


def real_norm2(a):
    return real_dot_product(a, a)


def complex_dot_product(w, z):
    n = len(w)
    return sum(w[j] * conjugate(z[j]) for j in range(n))


def complex_norm2(z):
    return complex_dot_product(z, z)


def product_scalar_by_vec(la, a):
    return list([la * a[j] for j in range(len(a))])


def sum_vec(a, b):
    return list([a[j] + b[j] for j in range(len(a))])


def dif_vec(a, b):
    return list([a[j] - b[j] for j in range(len(a))])


def kernel_poly_fock(m, w, z, alpha):
    n = len(z)
    CurrentRing = parent(z[0])
    factor0 = exp(alpha * complex_dot_product(w, z))
    diff0 = dif_vec(w, z)
    factor1 = laguerre_pol(m - 1, n, alpha * complex_norm2(diff0))
    result = CurrentRing(factor0 * factor1)
    result = result.full_simplify()
    return result


def list_of_multiindices(n, m):
    # J_{n,m} from the article
    result0 = lists_with_bounded_sum(n, m - 1)
    result1 = sorted(result0)
    return list(result1)


def q(k, xi, v):
    n = len(k)
    factor0 = 2 ** (n / 2)
    factor1 = pi ** (n / 4)
    result = factor0 * factor1
    for r in range(n):
        arg0 = (xi[r] + 2 * v[r]) / sqrt(2)
        f = hermite_function(k[r], arg0)
        result *= f
    return result


def operator_r(n, m, f):
    # f must be a function
    al, u, v, y, xi = vars_aluvyxi(n)
    assume(al > 0)
    z = [I * y[j] for j in range(n)]
    w = [u[j] + I * v[j] for j in range(n)]
    argument_of_f = product_scalar_by_vec(1 / sqrt(al), w)
    factor0 = f(argument_of_f)
    u2 = real_norm2(u)
    v2 = real_norm2(v)
    # v_minus_xi = dif_vec(v, xi)
    v_plus_xi = sum_vec(v, xi)
    # u_by_v_minus_xi = real_dot_product(u, v_minus_xi)
    u_by_v_plus_xi = real_dot_product(u, v_plus_xi)
    factor1 = exp(- (u2 / 2) - (v2 / 2) - I * u_by_v_plus_xi)
    J = list_of_multiindices(n, m)
    d = len(J)
    result = vector(SR, d)
    for j in range(d):
        k = J[j]
        factor2 = SR(1)
        for r in range(n):
            arg0 = (xi[r] + 2 * v[r]) / sqrt(2)
            factor2 *= hermite_function(k[r], arg0)
        g = factor0 * factor1 * factor2
        g = g.full_simplify()
        for r in range(n):
            g = g.full_simplify()
            gprev = g
            g = integral(gprev, u[r], -Infinity, Infinity, algorithm='giac')
            g = g.full_simplify()
            gprev = g
            g = integral(gprev, v[r], -Infinity, Infinity, algorithm='giac')
        g = (pi ** (- 3 * n / 4)) * g
        g = g.full_simplify()
        result[j] = g
    return result


def operator_r_applied_to_kernel(n, m):
    # we use the point i y / sqrt{\al}
    al, u, v, y, xi = vars_aluvyxi(n)
    assume(al > 0)
    z = [I * y[j] / sqrt(al) for j in range(n)]
    w = [u[j] + I * v[j] for j in range(n)]
    f = lambda argf : kernel_poly_fock(m, argf, z, al)
    result = operator_r(n, m, f)
    return result


def operator_r_applied_to_kernel_formula(n, m):
    al, u, v, y, xi = vars_aluvyxi(n)
    y2 = real_norm2(y)
    z = [I * y[j] for j in range(n)]
    w = [u[j] + I * v[j] for j in range(n)]
    factor0 = 2 ** (- n / 2)
    factor1 = exp(y2 / 2)
    J = list_of_multiindices(n, m)
    d = len(J)
    result = vector(SR, d)
    for j in range(d):
        k = J[j]
        factor2 = q(k, xi, y)
        comp = factor0 * factor1 * factor2
        result[j] = comp.full_simplify()
    return result


def symbolic_vectors_are_equal(vector0, vector1):
    d = len(vector0)
    result = True
    for j in range(d):
        comp = vector0[j] - vector1[j]
        comp = comp.full_simplify()
        result = result and comp.is_zero()
    return result


def test_operator_r_applied_to_kernel(n, m, verbose_level):
    if verbose_level >= 1:
        print('test_operator_r_applied_to_kernel')
        print('n = ' + str(n) + ', m = ' + str(m))
    result0 = operator_r_applied_to_kernel(n, m)
    result1 = operator_r_applied_to_kernel_formula(n, m)
    is_correct = symbolic_vectors_are_equal(result0, result1)
    if verbose_level >= 2:
        print('result using R:')
        print(result0, '\n')
        print('result that we have written in the example:')
        print(result1, '\n')
    if verbose_level >= 1:
        print('are equal? ' + str(is_correct) + '\n')
    return is_correct


def big_test_operator_r_applied_to_kernel(nmax, mmax, verbose_level):
    if verbose_level >= 1:
        print('big_test_operator_r_applied_to_kernel')
        print('nmax = ' + str(nmax) + ', mmax = ' + str(mmax) + '\n')
    big_result = True
    for n in range(1, nmax + 1):
        for m in range(1, mmax + 1):
            r = test_operator_r_applied_to_kernel(n, m, verbose_level)
            big_result = big_result and r
    if verbose_level >= 1:
        print('big_result = ' + str(big_result))
    return big_result


#print(test_operator_r_applied_to_kernel(1, 1, 2))

print(big_test_operator_r_applied_to_kernel(3, 3, 1))

