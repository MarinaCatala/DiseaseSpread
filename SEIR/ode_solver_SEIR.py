from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

'''
SEIR model with different R0 values
'''

# System constants
# R0 < 1
"""
S_0 = 7000
E_0 = 3000
I_0 = 0
R_0 = 0
N = S_0 + E_0 + I_0 + R_0
gamma = 0.6
beta = 0.00003
epsilon = 0.1
mu = 0.001
Lambda = mu * N
startTime = 0.0
stopTime = 10000.0
delta_t = 0.1
"""
# R0 = 1
#"""
S_0 = 7000
E_0 = 3000
I_0 = 0
R_0 = 0
N = S_0 + E_0 + I_0 + R_0
gamma = 0.5
beta = 0.00005
epsilon = 0.5
mu = 0.001
Lambda = mu * N
startTime = 0.0
stopTime = 10000.0
delta_t = 0.1
#"""
#R0 > 1
"""
S_0 = 9000
E_0 = 1000
I_0 = 0
R_0 = 0
N = S_0 + E_0 + I_0 + R_0
gamma = 0.008
beta = 0.00001
epsilon = 0.01
mu = 0.001
Lambda = mu * N
startTime = 0.0
stopTime = 10000.0
delta_t = 0.1
"""



R0 = (epsilon*beta*Lambda) / (mu*(epsilon + mu)*(gamma + mu))


# The ODEs to be solved
def odes(t, y):
    S,E,I,R=y
    dSdt = Lambda - beta*S*I - mu*S
    dEdt = beta*S*I - epsilon*E - mu*E
    dIdt = epsilon*E - gamma*I - mu*I
    dRdt = gamma*I - mu*R
    return [dSdt, dEdt, dIdt, dRdt]

if __name__ == '__main__':
    # initial values for S, E, I, R
    ic = [S_0, E_0, I_0, R_0]
    print(R0)
    # times where the solution is calculated
    t = np.arange(startTime, stopTime, delta_t)

    sol = solve_ivp(odes,
            t_span=[startTime, stopTime],
            y0=ic,
            t_eval=t)

    np.savetxt('odeEq.csv', np.c_[sol.t, sol.y.transpose()], fmt='%.4g', delimiter=',', header='t, S, E, I, R')

    S,E,I,R=sol.y
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
    ax.plot(t, S/N, 'b', alpha=0.5, lw=2, label='Susceptible')
    ax.plot(t, E/N, 'orange', alpha=0.5, lw=2, label='Exposed')
    ax.plot(t, I/N, 'r', alpha=0.5, lw=2, label='Infectious')
    ax.plot(t, R/N, 'g', alpha=0.5, lw=2, label='Removed')
    ax.set_xlabel('Time')
    ax.set_ylabel('Population')
    ax.set_ylim(0,1.1)
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(visible=True, which='major', c='w', lw=2, ls='-')
    ax.set_title("$N$ = {0} $R_0$ = {1:.2f}".format(str(N).zfill(4), R0))
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
    plt.show()