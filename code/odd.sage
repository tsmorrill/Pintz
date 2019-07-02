print("Siegel zeroes for primitive odd characters.")

var('t')

def error(A, tau, x, C1, C2):
    """Calculate E(q, tau, x) for precomputed A, tau, C1, C2."""

    alpha = 1 - tau

    def H6(d, tau):
        return log(d)/d^(1-tau)
    def H7(d, tau):
        return 1/d^(1-tau)
    def H8(d, tau):
        return 1/d

    # h6(z) = |d/dz H6(z)| &c

    def h6(d, tau):
        return ((1-tau)*log(d) - 1)/d^(2-tau)
    def h7(d, tau):
        return (1 - tau)/d^(2-tau)
    def h8(d, tau):
        return d**(-2)

    K1 = abs(C1)
    K2 = abs(C2)
    K3 = abs(x^tau/tau*(1/tau - log(x)))

    three_six = three_seven = three_eight = Infinity  # Initialize

    if x <= A:                   # trivial bound
        three_six = K1*(2*x*H6(x, tau) + numerical_integral(h6(t, tau), x, A)[0])
    three_six = min(three_six, K1*A*log(x)/x^(1-tau))

    if x <= A:                   # trivial bound
        three_seven = K2*(2*x*H7(x, tau) + numerical_integral(h7(t, tau), x, A)[0])
    three_seven = min(three_seven, K2*A/x^(1-tau))

    if x <= A:                   # trivial bound
        three_eight = K3*(2*x*H8(x, tau) + numerical_integral(h8(t, tau), x, A)[0])
    three_eight = min(three_eight, K3*A/x)

    # Spitting sum
    D1 = x**tau*log(x)/2/x     # lead term of W
    D2 = A*x**tau*log(x)/tau   # upper_sum
    z = sqrt(D2/D1)            # minimize D1*z + D2/z
    W = (x**tau*log(x)/2*z/x
         + x**tau*(1 + alpha*log(x))/24*z/x*(z/x + 1/x)
         + x**tau*(alpha*(alpha + 1)*log(x) + (3*alpha**2 + 6*alpha + 2)/(alpha + 2))
         /216/sqrt(3)*z/x*(z/x + 1/x)*(2*z/x + 1/x))
    upper_sum = A/z*x**tau*log(x)/tau

    # Not splitting sum
    # W = (x**tau*log(x)*(1/2 + (1-tau)/12 + (1-tau)*(2-tau)/36/sqrt(3))
    #    + x**tau/36/sqrt(3)*((1-tau)*(2-tau)*log(x) + (3*(1-tau)**2 + 6*(1-tau) + 2)/(3-tau)))

    # the_rest = minimize(W + upper_sum, [sqrt(A*x)])[0]

    # print(three_six.n(), three_seven.n(), three_eight.n(), W.n(), upper_sum.n())

    number = (three_six + three_seven + three_eight + W + upper_sum)
    return number

def F(c, q0, q1, x):
    """Calculate an upper bound for F on the interval [q0, q1] for fixed c and x."""

    A = sqrt(q1) * log(q1) * (pi**-2 + 0.5/log(q0))    # chi is even
    # A = sqrt(q1) * log(q1) * (0.5/pi + 1/log(q0))    # chi is odd
    tau = c/log(q1)

    alpha = 1 - tau
    C1 = (0.5 + alpha/12 - 1/(1-alpha) - alpha*(alpha + 1)*(alpha + 2)/6
         *numerical_integral((frac(t)**3 - 1.5*frac(t)**2 + 0.5*frac(t))
                             /t**(alpha + 3), 1, Infinity)[0])
    C2 = (11/12 + (1 - alpha)**-2 + 1/6*numerical_integral(
          (2 + 6*alpha + 3*alpha^2 - alpha*(alpha + 1)*(alpha + 2)*log(t))
          *(frac(t)**3 - 1.5*frac(t)**2 + 0.5*frac(t))
          /t**(alpha+3), 1, Infinity)[0])

    if q0 <= 10**7:
        Bennet = 79.2
    else:
        Bennet = 12

    L1 = (1/tau - log(x))*Bennet*x^tau/tau/sqrt(q1)
    lower_bound = max(-2*zetaderiv(1, 2-2*tau)
                      -2*(1 + (0.5 - tau)*log(x))/x^(0.5 - tau)/(1 - 2*tau)^2,
                      log(4)/4)
    return (float(error(A, tau, x, C1, C2) - L1 - lower_bound))

def search(c, q0, q1, x0, x1, step):
    """Find an x on the interval [x0, x1] which minimizes F for c and q in [q0, q01],
    searching in increments of 10**step. Then, return the tuple (bool(F < 0), log(x, 10)).
    Returns None if the specified range and step size have no test values.
    """

    tau = c/log(q1)
    A = sqrt(q1) * log(q1) * (pi**-2 + 0.5/log(q0))    # chi is even
    # A = sqrt(q1) * log(q1) * (0.5/pi + 1/log(q0))    # chi is odd

    x_min = max(x0, log((exp(1/(2*tau - 1)) + 1)**2, 10))    # lower bound from Lemma 4
    x_max = min(x1, log(q0^(1/c), 10))    # upper bound from (10)
    x_range = [x_min + i*step for i in range(floor((x_max-x_min)/step) + 1)]
    values = []

    if not x_range:
        return None

    alpha = 1 - tau
    C1 = (0.5 + alpha/12 - 1/(1-alpha) - alpha*(alpha + 1)*(alpha + 2)/6
          *numerical_integral((frac(t)**3 - 1.5*frac(t)**2 + 0.5*frac(t))
                              /t**(alpha + 3), 1, Infinity)[0])
    C2 = 11/12 + (1/(1 - alpha)^2 + 1/6*numerical_integral(
                       (2 + 6*alpha + 3*alpha^2 - alpha*(alpha + 1)*(alpha +  2)
                        *log(t))*(frac(t)**3 - 1.5*frac(t)**2 + 0.5*frac(t))
                        /t**(alpha+3), 1, Infinity)[0])

    for log_x in x_range:
        x = 10^log_x
        if q0 <= 10**7:
            Bennet = 79.2
        else:
            Bennet = 12

        E = error(A, tau, x, C1, C2)
        L1 = (1/tau - log(x))*Bennet*x^tau/tau/sqrt(q1)
        lower_bound = max(-2*zetaderiv(1, 2-2*tau)
                              -2*(1 + (0.5 - tau)*log(x))/x^(0.5 - tau)/(1 - 2*tau)^2,
                              log(4)/4)

        # print('E = {}'.format(E.n()))
        # print('L1 = {}'.format(-L1.n()))
        # print('lower_bound = {}'.format(-lower_bound.n()))

        values.append((float(E - L1 - lower_bound), log_x))
    F, x = (float(i) for i in min(values))
    output = (F < 0)

    return(output, x)

def best_c(q0, q1, significant_figures=2):
    """Calculate the maximal c so that F < 0 on the interval [q0, q1]."""

    c_true, c_step = 0, 1.0
    record_significant_figures, current_figures = False, 0
    done = False

    while not done:
        print('Trying increments of {}.'.format(c_step))
        c, it_works = c_true + c_step, True
        while it_works:
            x0, x1, step = 0, log(q0^(1/c), 10), 1
            for i in range(3):
                result, x = search(c, q0, q1, x0, x1, step)
                x0, x1, step = max(0, x - step), min(x1, x + step), step/10
                if result:
                    it_works = True
                    record_significant_figures = True
                    x_true, c_true = x, c
                    print("q in [{}, {}], c >= {}.".format(q0, q1, c_true))
                    c += c_step
                    break
                else:
                    it_works = False
        if record_significant_figures == True:
            current_figures += 1
            if current_figures == significant_figures:
                done = True
        c_step /= 10

    q0_magnitude = int(log(q0, 10))
    q0_lead = round(q0/10**q0_magnitude, 2)
    q1_magnitude = int(log(q1, 10))
    q1_lead = round(q1/10**q1_magnitude, 2)

    str1 = "${{{}}} \cdot 10^{{{}}}$ & ${{{}}} \cdot 10^{{{}}}$ & \\num{{{}}} & $10^{{{}}}$ \\\\".format(q0_lead, q0_magnitude, q1_lead, q1_magnitude, round(c_true,5), x_true)
    str2 = 'F({}, {}, {}, 10^{})'.format(c_true, q0, q1, x_true)

    print('')
    print('c = {}.'.format(c_true))
    print('')

    return(str1, str2)

def best_q1(q0, c):
    """Calculate the maximal q1 so that F < 0 on the interval [q0, q1] for fixed c."""
    q1_true, q_step = q0, 10**int(log(q0, 10) + 4)
    done = False
    while not done:
        q1, it_works = q1_true + q_step, True
        q1 = int(q1/10**int(log(q1, 10)-1))*10**int(log(q1, 10))
        while it_works:
            x0, x1, step = 0, log(q0^(1/c), 10), 1
            for i in range(4):
                result, x = search(c, q0, q1, x0, x1, step)
                x0, x1, step = max(0, x - step), min(x1, x + step), step/10
                if result:
                    it_works = True
                    record_significant_figures = True
                    x_true, q1_true = x, q1
                    print("q1 >= {}.".format(q1_true))
                    q1 += q_step
                else:
                    it_works = False
        q_step /= 10
        done = (int(log(q_step, 10)) <= int(log(q1_true, 10) - 2))
    str = 'F({}, {}, {}, 10^{})'.format(c, q0, q1_true, x_true)
    print(str)
    return(q1_true)

def q_list(q0, c, stop):
    """Calculate a list of q values starting at q0 so that F < 0 for fixed c."""
    q_list, q1 = [q0], None
    while q1 < stop:
        print(q0)
        q1 = best_q1(q0, c)
        q_list.append(q1)
        q0 = q1
    return(q_list)

def cq_table(q_list, significant_figures=4):
    """Generate a table of c values corresponding to q_list formatted for LaTeX, then generate Sage
    commands to verify these calculations.
    """
    TeX_list, Sage_list = [], []
    for q0, q1 in zip(q_list[:-1], q_list[1:]):
        print(q0, q1)
        TeX_str, Sage_str = best_c(q0, q1, significant_figures=significant_figures)
        print('====')
        TeX_list.append(TeX_str)
        Sage_list.append(Sage_str)
    for item in TeX_list:
        print(item)
    for item in Sage_list:
        print(item)
