# Here we test the formula for R_{\cH_m} K_{x + iy} from our paper
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


def vars_uvxyxi(n):
    u = myvars('u', n)
    v = myvars('v', n)
    x = myvars('x', n)
    y = myvars('y', n)
    xi = myvars('xi', n)
    return u, v, x, y, xi


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


def kernel_space_h(m, u, v, x, y):
    n = len(u)
    CurrentRing = parent(u[0])
    u_minus_x = dif_vec(u, x)
    v_minus_y = dif_vec(v, y)
    v_plus_y = sum_vec(v, y)
    term0 = real_norm2(u_minus_x) + real_norm2(v_minus_y)
    term1 = real_dot_product(u_minus_x, v_plus_y)
    factor0 = 2 ** n
    factor1 = exp((- 1 / 2) * term0 - I * term1)
    factor2 = laguerre_pol(m - 1, n, term0)
    result = CurrentRing(factor0 * factor1 * factor2)
    result = result.full_simplify()
    return result


def kernel_space_h_example(m, u, v, x, y):
    n = len(u)
    CurrentRing = parent(u[0])
    factor0 = 2 ** (- n / 2)
    x2 = real_norm2(x)
    y2 = real_norm2(y)
    xy = real_dot_product(x, y)
    factor1 = exp((1 / 2) * (x2 + y2) - I * xy)
    factor2 = kernel_space_h(m, u, v, x, y)
    result = CurrentRing(factor0 * factor1 * factor2)
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


def operator_rh(n, m, g):
    # operator R^H from the article
    # g must be a function
    u, v, x, y, xi = vars_uvxyxi(n)
    u_xi = real_dot_product(u, xi)
    factor0 = exp(- I * u_xi)
    factor1 = g(u, v)
    J = list_of_multiindices(n, m)
    d = len(J)
    result = vector(SR, d)
    for j in range(d):
        k = J[j]
        factor2 = q(k, xi, v)
        g = factor0 * factor1 * factor2
        g = g.full_simplify()
        for r in range(n):
            g = g.expand()
            gprev = g
            g = integral(gprev, u[r], -Infinity, Infinity, algorithm='giac')
            g = g.expand()
            gprev = g
            g = integral(gprev, v[r], -Infinity, Infinity, algorithm='giac')
        coef = ((2 * pi) ** (-n))
        g = coef * g
        g = g.full_simplify()
        result[j] = g
    return result


def operator_rh_explicit(n, m, g):
    # operator R^H from the article, in a more explicit form
    # g must be a function
    u, v, x, y, xi = vars_uvxyxi(n)
    u_xi = real_dot_product(u, xi)
    factor0 = exp(- I * u_xi)
    factor1 = g(u, v)
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
            g = g.expand()
            gprev = g
            g = integral(gprev, u[r], -Infinity, Infinity, algorithm='giac')
            g = g.expand()
            gprev = g
            g = integral(gprev, v[r], -Infinity, Infinity, algorithm='giac')
        coef = (2 ** (- n / 2)) * (pi ** (- 3 * n / 4))
        g = coef * g
        g = g.full_simplify()
        result[j] = g
    return result


def operator_rh_applied_to_kernel(n, m):
    u, v, x, y, xi = vars_uvxyxi(n)
    f = lambda s, t: kernel_space_h_example(m, s, t, x, y)
    result = operator_rh(n, m, f)
    return result


def operator_rh_explicit_applied_to_kernel(n, m):
    u, v, x, y, xi = vars_uvxyxi(n)
    f = lambda s, t: kernel_space_h_example(m, s, t, x, y)
    result = operator_rh_explicit(n, m, f)
    return result


def example_vector_function_expr(n, m):
    # x, y are parameters, xi is a variable
    u, v, x, y, xi = vars_uvxyxi(n)
    x2 = real_norm2(x)
    y2 = real_norm2(y)
    y_plus_xi = sum_vec(y, xi)
    x_by_y_plus_xi = real_dot_product(x, y_plus_xi)
    factor0 = 2 ** (- n / 2)
    factor1 = exp((x2 + y2) / 2 - I * x_by_y_plus_xi)
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


def test_operator_rh_applied_to_kernel(n, m, verbose_level):
    if verbose_level >= 1:
        print('test_operator_rh_applied_to_kernel')
        print('n = ' + str(n) + ', m = ' + str(m))
    result0 = operator_rh_applied_to_kernel(n, m)
    result1 = operator_rh_explicit_applied_to_kernel(n, m)
    result2 = example_vector_function_expr(n, m)
    is_correct0 = symbolic_vectors_are_equal(result0, result2)
    is_correct1 = symbolic_vectors_are_equal(result1, result2)
    is_correct = is_correct0 and is_correct1
    if verbose_level >= 2:
        print('result using RH:')
        print(result0, '\n')
        print('result using RH in a more explicit form:')
        print(result1, '\n')
        print('result that we have written in the example:')
        print(result2, '\n')
    if verbose_level >= 1:
        print('are equal? ' + str(is_correct) + '\n')
    return is_correct


def big_test_operator_rh_applied_to_kernel(nmax, mmax, verbose_level):
    if verbose_level >= 1:
        print('big_test_operator_rh_applied_to_kernel')
        print('nmax = ' + str(nmax) + ', mmax = ' + str(mmax) + '\n')
    big_result = True
    for n in range(1, nmax + 1):
        for m in range(1, mmax + 1):
            r = test_operator_rh_applied_to_kernel(n, m, verbose_level)
            big_result = big_result and r
    if verbose_level >= 1:
        print('big_result = ' + str(big_result))
    return big_result


print(test_operator_rh_applied_to_kernel(1, 1, 2))

print(big_test_operator_rh_applied_to_kernel(3, 3, 1))

