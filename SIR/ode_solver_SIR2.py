from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

'''
SIR model on networks
'''

# System constants
# R0 < 1
"""
S_0 = 7000
I_0 = 3000
R_0 = 0
gamma = 0.6
beta = 0.2
"""
# R0 = 1
"""
S_0 = 7000
I_0 = 3000
R_0 = 0
gamma = 0.2
beta = 0.2
"""
# R0 > 1
#""""
S_0 = 9000
I_0 = 1000
R_0 = 0
gamma = 0.03
beta = 0.2
#"""


N_0 = S_0 + I_0 + R_0
startTime = 0.0
stopTime = 175.0
delta_t = 0.1
R0 = (beta / gamma) 


# The ODEs to be solved
def odes(t, y):
    S,I,R=y
    dSdt = -beta*S*I/N_0
    dIdt = beta*S*I/N_0 - gamma*I 
    dRdt = gamma*I
    return [dSdt, dIdt, dRdt]

if __name__ == '__main__':
    # initial values for S, I, R
    ic = [S_0, I_0, R_0]
    print(R0)
    # times where the solution is calculated
    t = np.arange(startTime, stopTime, delta_t)

    sol = solve_ivp(odes,
            t_span=[startTime, stopTime],
            y0=ic,
            t_eval=t)

    np.savetxt('odeGreater.csv', np.c_[sol.t, sol.y.transpose()], fmt='%.4g', delimiter=',', header='t, S, I, R')

    # We plot the data on three separate curves for S(t), I(t) and R(t)
    S,I,R=sol.y
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
    ax.plot(t, S/N_0, 'b', alpha=0.5, lw=2, label='Susceptible')
    ax.plot(t, I/N_0, 'r', alpha=0.5, lw=2, label='Infectious')
    ax.plot(t, R/N_0, 'g', alpha=0.5, lw=2, label='Removed')
    ax.set_xlabel('Time')
    ax.set_ylabel('Population')
    ax.set_ylim(0,1.1)
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(visible=True, which='major', c='w', lw=2, ls='-')
    ax.set_title("$N$ = {0} $R_0$ = {1:.2f}".format(str(N_0).zfill(4), R0))
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
    plt.show()
