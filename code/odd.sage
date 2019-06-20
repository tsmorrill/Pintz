print("Siegel zeroes for primitive odd characters.")

var('t')

def error(A, tau, T, C1, C2):
    """Calculate E(q, tau, T + 1/2) for precomputed A, tau, C1, C2."""

    x = int(T) + 1/2

    def h6(d, tau):
        return log(d)/d^(1-tau)
    def h7(d, tau):
        return 1/d^(1-tau)
    def h8(d, tau):
        return 1/d

    K1 = C1
    K2 = C2
    K3 = x^tau/tau*(log(x) - 1/tau)

    if x <= A:                   # trivial bound
        three_six = -K1*(2*x*h6(x, tau) + numerical_integral(h6(t, tau), x, A)[0])
    else:                        # PV inequality
        three_six = -K1*A*log(x)/x^(1-tau)

    if x <= A:                   # trivial bound
        three_seven = -K2*(2*x*h7(x, tau) + numerical_integral(h7(t, tau), x, A)[0])
    else:                        # PV inequality
        three_seven = -K2*A/x^(1-tau)

    if x <= A:                   # trivial bound
        three_eight = -K3*(2*x*h8(x, tau) + numerical_integral(h8(t, tau), x, A)[0])
    else:                        # PV inequality
        three_eight = -K3*A/x

    W = x^tau*((1-tau)*log(x) - 1)*(1 + 1/x)/24

    number = (three_six + three_seven + three_eight + W)
    return number

def constants(c, q0, q1):
    """Calculate C_1, C_2 on the interval [q0, q1] for fixed c."""

    A = sqrt(q1) * log(q1) * (1/2/pi + 1/log(q0))    # chi is odd
    tau = c/log(q1)
    alpha = 1 - tau
    C1 = 0.5 + alpha/12 - 1/(1-alpha) - alpha*(alpha + 1)*(alpha + 2)/6*numerical_integral((frac(t)**3 - 1.5*frac(t)**2 + 0.5*frac(t))/t**(alpha + 3), 1, Infinity)[0]
    C2 = 11/12 + (1/(1 - alpha)^2 + 1/6*numerical_integral(
          (2 + 6*alpha + 3*alpha^2 - alpha*(alpha + 1)*(alpha +  2)
          *log(t))*(frac(t)**3 - 1.5*frac(t)**2 + 0.5*frac(t))
          /t**(alpha+3), 1, Infinity)[0])

    return(C1, C2)

def F(c, q0, q1, T):
    """Calculate an upper bound for F on the interval [q0, q1] for fixed c and x."""

    x = T + 1/2

    A = sqrt(q1) * log(q1) * (1/2/pi + 1/log(q0))    # chi is odd
    tau = c/log(q1)
    alpha = 1 - tau
    C1 = 0.5 + alpha/12 - 1/(1-alpha) - alpha*(alpha + 1)*(alpha + 2)/6*numerical_integral((frac(t)**3 - 1.5*frac(t)**2 + 0.5*frac(t))/t**(alpha + 3), 1, Infinity)[0]
    C2 = 11/12 + (1/(1 - alpha)^2 + 1/6*numerical_integral(
          (2 + 6*alpha + 3*alpha^2 - alpha*(alpha + 1)*(alpha +  2)
          *log(t))*(frac(t)**3 - 1.5*frac(t)**2 + 0.5*frac(t))
          /t**(alpha+3), 1, Infinity)[0])
    if q0 <= 10**7:
        Bennet = 79.2
    else:
        Bennet = 12
    lower_bound = (-(1/tau-log(x))*Bennet*x^tau/tau/sqrt(q1)
                    + 2*zetaderiv(1, 2-2*tau)
                    + 2*(1 + (0.5 - tau)*log(x))/x^(0.5 - tau)/(1 - 2*tau)^2)
    return (float(error(A, tau, x, C1, C2) + lower_bound))

def search(c, q0, q1, x0, x1, step):
    """Find an x on the interval [x0, x1] which minimizes F for c and q in [q0, q01],
    searching in increments of 10**step. Then, return the tuple (bool(F < 0), log(x, 10)).
    """

    tau = c/log(q1)
    A = sqrt(q1) * log(q1) * (1/2/pi + 1/log(q0))    # chi is odd
    x_min = max(x0, log((exp(1/(2*tau - 1)) + 1)**2, 10))    # lower bound from Lemma 4
    x_max = min(x1, log(q0^(1/c), 10))    # upper bound from (10)
    x_range = [x_min + i*step for i in range(floor((x_max-x_min)/step) + 1)]
    values = []

    for log_x in x_range:
        T = int(10^log_x)
        x = T + 1/2
        alpha = 1 - tau
        C1 = 0.5 + alpha/12 - 1/(1-alpha) - alpha*(alpha + 1)*(alpha + 2)/6*numerical_integral((frac(t)**3 - 1.5*frac(t)**2 + 0.5*frac(t))/t**(alpha + 3), 1, Infinity)[0]
        C2 = 11/12 + (1/(1 - alpha)^2 + 1/6*numerical_integral(
              (2 + 6*alpha + 3*alpha^2 - alpha*(alpha + 1)*(alpha +  2)
              *log(t))*(frac(t)**3 - 1.5*frac(t)**2 + 0.5*frac(t))
              /t**(alpha+3), 1, Infinity)[0])
        if q0 <= 10**7:
            Bennet = 79.2
        else:
            Bennet = 12
        lower_bound = (-(1/tau-log(x))*Bennet*x^tau/tau/sqrt(q1)
                       + 2*zetaderiv(1, 2-2*tau)
                       + 2*(1 + (0.5 - tau)*log(x))/x^(0.5 - tau)/(1 - 2*tau)^2)
        values.append((float(error(A, tau, x, C1, C2)
                      + lower_bound), log_x))
    F, x = (float(i) for i in min(values))
    if F < 0:
        output = True
    else:
        output = False
    return(output, x)

def best_c(q0, q1):
    """Calculate the maximal c so that F < 0 on the interval [q0, q1]."""

    c_true, c_step = 0, 1.0
    record_significant_figures, significant_figures = False, 0
    done = False

    while not done:
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
            significant_figures += 1
            if significant_figures == 4:
                done = True
        c_step /= 10
    q0_magnitude = int(log(q0, 10))
    q0_lead = round(q0/10**q0_magnitude, 2)
    q1_magnitude = int(log(q1, 10))
    q1_lead = round(q1/10**q1_magnitude, 2)
    str1 = "${{{}}} \cdot 10^{{{}}}$ & ${{{}}} \cdot 10^{{{}}}$ & \\num{{{}}} & $10^{{{}}}$ \\\\".format(q0_lead, q0_magnitude, q1_lead, q1_magnitude, round(c_true,5), x_true)
    str2 = 'F({}, {}, {}, 10^{})'.format(c_true, q0, q1, x_true)
    return(str1, str2)

def best_q1(q0, c):
    """Calculate the maximal q1 so that F < 0 on the interval [q0, q1] for fixed c."""
    q1_true, q_step = q0, 10**int(log(q0, 10) + 4)
    significant_figures, record_significant_figures = 0, False
    done = False
    while not done:
        q1, it_works = q1_true + q_step, True
        while it_works:
            x0, x1, step = 0, log(q0^(1/c), 10), 1
            for i in range(3):
                result, x = search(c, q0, q1, x0, x1, step)
                x0, x1, step = max(0, x - step), min(x1, x + step), step/10
                if result:
                    it_works = True
                    record_significant_figures = True
                    x_true, q1_true = x, q1
                    print("q1 >= {}.".format(q1_true))
                    q1 += q_step
                    break
                else:
                    it_works = False
        if record_significant_figures == True:
            significant_figures += 1
            if significant_figures == 2:
                done = True
        q_step /= 10
    str = 'F({}, {}, {}, 10^{})'.format(c, q0, q1_true, x_true)
    return(str)

def q_list(q0, c):
    """Calculate a list of q values starting at q0 so that F < 0 for fixed c."""
    q_list, q1 = [q0], None


def cq_table(q_list):
    """Generate a table of c values corresponding to q_list formatted for LaTeX, then generate Sage
    commands to verify these calculations.
    """
    TeX_list, Sage_list = [], []
    for q0, q1 in zip(q_list[:-1], q_list[1:]):
        TeX_str, Sage_str = best_c(q0, q1)
        TeX_list.append(TeX_str)
        Sage_list.append(Sage_str)
    for item in TeX_list:
        print(item)
    for item in Sage_list:
        print(item)
