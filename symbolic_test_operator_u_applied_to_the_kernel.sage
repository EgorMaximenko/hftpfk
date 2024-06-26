# Here we test the formula for U K_{(x + iy) / sqrt{a}} from our paper
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


def vars_aluvxy(n):
    al = var('al')
    u = myvars('u', n)
    v = myvars('v', n)
    x = myvars('x', n)
    y = myvars('y', n)
    return al, u, v, x, y


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


def kernel_poly_fock(m, alpha, w, z):
    n = len(z)
    CurrentRing = parent(z[0])
    factor0 = exp(alpha * complex_dot_product(w, z))
    diff0 = dif_vec(w, z)
    factor1 = laguerre_pol(m - 1, n, alpha * complex_norm2(diff0))
    result = CurrentRing(factor0 * factor1)
    result = result.full_simplify()
    return result


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


def operator_u_applied_to_kf_example_formula(m, u, v, x, y):
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


def operator_u(al1, n, m, f, u1, v1):
    w = [(u1[j] + I * v1[j]) / sqrt(al1) for j in range(n)]
    factor0 = 2 ** (n / 2)
    u2 = real_norm2(u1)
    v2 = real_norm2(v1)
    uv = real_dot_product(u1, v1)
    factor1 = exp(- (1 / 2) * (u2 + v2) - I * uv)
    factor2 = f(w)
    g = factor0 * factor1 * factor2
    g = g.full_simplify()
    return g


def complex_vector_from_real_vectors(u1, v1):
    n = len(u1)
    return [u1[j] + I * v[j] for j in range(n)]


def test_operator_u_applied_to_kf_example(n, m, verbose_level):
    if verbose_level >= 1:
        print('test_operator_u_applied_to_kf_example')
        print('n = ' + str(n) + ', m = ' + str(m))
    al, u, v, x, y = vars_aluvxy(n)
    assume(al > 0)
    z = [(x[j] + I * y[j]) / sqrt(al) for j in range(n)]
    f1 = lambda w: kernel_poly_fock(m, al, w, z)
    result0 = operator_u(al, n, m, f1, u, v)
    result1 = operator_u_applied_to_kf_example_formula(m, u, v, x, y)
    diff_result0_result1 = (result0 - result1).full_simplify()    
    is_correct = diff_result0_result1.is_zero()
    if verbose_level >= 2:
        print('operator_u_applied_to_kf_example:')
        print(result0, '\n')
        print('operator_u_applied_to_kf_example_formula:')
        print(result1, '\n')
    if verbose_level >= 1:
        print('are equal? ' + str(is_correct) + '\n')
    return is_correct


def big_test_operator_u_applied_to_kf_example(nmax, mmax, verbose_level):
    if verbose_level >= 1:
        print('big_test_operator_u_applied_to_kernel')
        print('nmax = ' + str(nmax) + ', mmax = ' + str(mmax) + '\n')
    big_result = True
    for n in range(1, nmax + 1):
        for m in range(1, mmax + 1):
            r = test_operator_u_applied_to_kf_example(n, m, verbose_level)
            big_result = big_result and r
    if verbose_level >= 1:
        print('big_result = ' + str(big_result))
    return big_result


print(test_operator_u_applied_to_kf_example(1, 1, 2))

print(big_test_operator_u_applied_to_kf_example(4, 4, 1))

