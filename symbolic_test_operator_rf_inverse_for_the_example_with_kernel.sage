# Here we test the formula for R_{\cF_{m,\alpha}}^\ast from our paper
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


def vars_aluvxyxi(n):
    al = var('al')
    u = myvars('u', n)
    v = myvars('v', n)
    x = myvars('x', n)
    y = myvars('y', n)
    xi = myvars('xi', n)
    return al, u, v, x, y, xi


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


def kernel_poly_fock_example_expr(n, m):
    # Here we return a symbolic expression.
    # We use the point (x + i y) / sqrt{\al}
    al, u, v, x, y, xi = vars_aluvxyxi(n)
    assume(al > 0)
    coef = 1 / sqrt(al)
    z = [coef * (x[j] + I * y[j]) for j in range(n)]
    w = [u[j] + I * v[j] for j in range(n)]
    f = kernel_poly_fock(m, w, z, al)
    return f


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


def operator_rf_inverse(n, m, h):
    # operator R_{F_{m,\alpha}}^\ast from the article
    # h must be a vector-function
    al, u, v, x, y, xi = vars_aluvxyxi(n)
    assume(al > 0)
    factor0 = (2 ** (-n)) * (pi ** (- n / 2))
    u2 = real_norm2(u)
    v2 = real_norm2(v)
    uv = real_dot_product(u, v)
    u_xi = real_dot_product(u, xi)
    factor1 = exp((al * u2 / 2) + (al * v2 / 2) + (I * al * uv))
    J = list_of_multiindices(n, m)
    d = len(J)
    mysum = SR(0)
    for j in range(d):
        k = J[j]
        factor2 = exp(sqrt(al) * I * u_xi)
        argument_of_q = product_scalar_by_vec(sqrt(al), v)
        factor3 = q(k, xi, argument_of_q)
        factor4 = h[j]
        fun = factor2 * factor3 * factor4
        fun = fun.full_simplify()
        for r in range(n):
            fun_prev = fun
            fun_prev = fun_prev.expand()
            fun = integral(fun_prev, xi[r], -Infinity, Infinity, algorithm='giac')
            fun = fun.full_simplify()
        mysum += fun
    mysum = mysum.full_simplify()
    result = factor0 * factor1 * mysum
    result = result.full_simplify()
    return result


def operator_rf_inverse_explicit(n, m, h):
    # operator R_{F_{m,\alpha}}^\ast from the article, in a more explicit form
    # h must be a vector whose components are expressions
    al, u, v, x, y, xi = vars_aluvxyxi(n)
    assume(al > 0)
    factor0 = (2 ** (- n / 2)) * (pi ** (- n / 4))
    u2 = real_norm2(u)
    v2 = real_norm2(v)
    uv = real_dot_product(u, v)
    u_xi = real_dot_product(u, xi)
    factor1 = exp((al * u2 / 2) + (al * v2 / 2) + (I * al * uv))
    J = list_of_multiindices(n, m)
    d = len(J)
    mysum = SR(0)
    for j in range(d):
        factor2 = exp(sqrt(al) * I * u_xi)
        k = J[j]
        factor3 = SR(1)
        for r in range(n):
            arg0 = (xi[r] + 2 * sqrt(al) * v[r]) / sqrt(2)
            factor3 *= hermite_function(k[r], arg0)
        factor4 = h[j]
        fun = factor2 * factor3 * factor4
        fun = fun.full_simplify()
        for r in range(n):
            fun_prev = fun
            fun_prev = fun_prev.expand()
            fun = integral(fun_prev, xi[r], -Infinity, Infinity, algorithm='giac')
            fun = fun.full_simplify()
        mysum += fun
    mysum = mysum.expand()
    mysum = mysum.full_simplify()
    result = factor0 * factor1 * mysum
    mysum = mysum.expand()
    result = result.full_simplify()
    return result


def example_vector_function_expr(n, m):
    # x, y are parameters, xi is a variable
    al, u, v, x, y, xi = vars_aluvxyxi(n)
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


def operator_rf_inverse_applied_to_example_vector_function(n, m):
    h = example_vector_function_expr(n, m)
    f = operator_rf_inverse(n, m, h)
    return f


def operator_rf_inverse_explicit_applied_to_example_vector_function(n, m):
    h = example_vector_function_expr(n, m)
    f = operator_rf_inverse_explicit(n, m, h)
    return f


def test_operator_rf_inverse_applied_to_example(n, m, verbose_level):
    if verbose_level >= 1:
        print('test_operator_rf_inverse_applied_to_example')
        print('n = ' + str(n) + ', m = ' + str(m))
    result0 = operator_rf_inverse_applied_to_example_vector_function(n, m)
    result1 = operator_rf_inverse_explicit_applied_to_example_vector_function(n, m)
    result2 = kernel_poly_fock_example_expr(n, m)
    diff0 = result0 - result2
    diff0 = diff0.full_simplify()
    is_correct0 = diff0.is_zero()
    diff1 = result1 - result2
    diff1 = diff1.full_simplify()
    is_correct1 = diff1.is_zero()
    is_correct1 = True
    is_correct = is_correct0 and is_correct1
    if verbose_level >= 2:
        print('result using RFinverse:')
        print(result0, '\n')
        print('result using RFinverse, in a more explicit form:')
        print(result1, '\n')
        print('result corresponding to our computations in the example:')
        print(result2, '\n')
    if verbose_level >= 1:
        print('are equal? ' + str(is_correct) + '\n')
    return is_correct


def big_test_operator_rf_inverse_applied_to_example(nmax, mmax, verbose_level):
    if verbose_level >= 1:
        print('big_test_operator_rf_inverse_applied_to_example')
        print('nmax = ' + str(nmax) + ', mmax = ' + str(mmax) + '\n')
    big_result = True
    for n in range(1, nmax + 1):
        for m in range(1, mmax + 1):
            r = test_operator_rf_inverse_applied_to_example(n, m, verbose_level)
            big_result = big_result and r
    if verbose_level >= 1:
        print('big_result = ' + str(big_result))
    return big_result


print(test_operator_rf_inverse_applied_to_example(2, 2, 2))

print(big_test_operator_rf_inverse_applied_to_example(3, 3, 1))

