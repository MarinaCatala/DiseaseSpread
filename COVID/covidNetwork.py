import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import csv

'''
COVID-19 Small World Network Simulation
'''

MAX_NUM_GENERATIONS = 1500
NUM_ITERATIONS = 20
IS_LOCKDOWN = False
IS_VACCINATION = False
S_0 = 9900
E_0 = 100
N = S_0 + E_0
S = [[] for i in range(NUM_ITERATIONS)]
V = [[] for i in range(NUM_ITERATIONS)]
E = [[] for i in range(NUM_ITERATIONS)]
EV = [[] for i in range(NUM_ITERATIONS)]
I = [[] for i in range(NUM_ITERATIONS)]
IV = [[] for i in range(NUM_ITERATIONS)]
H = [[] for i in range(NUM_ITERATIONS)]
HV = [[] for i in range(NUM_ITERATIONS)]
R = [[] for i in range(NUM_ITERATIONS)]
RV = [[] for i in range(NUM_ITERATIONS)]
D = [[] for i in range(NUM_ITERATIONS)]

beta = 0.4
sigma = 1/7
eta = 0.0055
gamma = 0.0083
gammap = 2.7E-4
delta1 = 1/7
delta1p = 1/7
delta2 = 3.5E-4
delta2p = 2.8E-4
delta3 = 0.01
delta3p = 0.005
zeta2 = 0.018
zeta2p = 0.018
zeta1 = 1/14
zeta1p = 1/14

if(IS_VACCINATION):
    p = 0.05
else:
    p = 0

if(IS_LOCKDOWN):
    lockdownStart = 15
else:
    lockdownStart = -1

lockdownPeriod = 60
lockdownEnd = lockdownStart + lockdownPeriod
lockdownEssentialContactRate = 0.5
lockdownNonEssentialContactRate = 0.1
removedEdges = []


rateEssential = 0.15

def initialise_graph(G):
    # initialise the graph

    for n in G.nodes():
        G.nodes[n]['state'] = 'S'
        if(np.random.random() < rateEssential):
            G.nodes[n]['essential'] = True
        else:
            G.nodes[n]['essential'] = False

    # randomly pick E_0 nodes to be exposed
    initial_infected_node = np.random.choice(G.nodes(),E_0)

    for i in initial_infected_node:
        G.nodes[i]['state'] = 'E'


def plot_data():

    plt.figure()
    plt.plot(S, color = '#0000FF', linestyle = 'solid', label = 'S')
    plt.plot(V, color = '#00FFFF', linestyle = 'dashed', label = 'V')
    plt.plot(E+EV, color = '#FF8000', linestyle = 'solid', label = 'E+EV')
    plt.plot(I+IV, color = '#FF0000', linestyle = 'solid', label = 'I+IV')
    plt.plot(H+HV, color = '#00bfaf', linestyle = 'solid', label = 'H+HV')
    plt.plot(R+RV, color = '#00FF00', linestyle = 'solid', label = 'R+RV')
    plt.plot(D, color = '#000000', linestyle = 'solid', label = 'D')
    plt.axvline(x = 15, color = '#FFD9EE',linestyle="dashed")
    plt.axvline(x = 75, color = '#FFD9EE', linestyle="dashed")
    

    plt.xlabel("Days")
    plt.ylabel("Proportion of people")
    plt.xlim(0, 365)

    plt.title("COVID-19 Simulation")
    plt.legend()
    plt.savefig("covida.png")

    plt.figure()
    plt.plot(S, color = '#0000FF', linestyle = 'solid', label = 'S')
    plt.plot(V, color = '#00FFFF', linestyle = 'dashed', label = 'V')
    plt.plot(E+EV, color = '#FF8000', linestyle = 'solid', label = 'E+EV')
    plt.plot(I+IV, color = '#FF0000', linestyle = 'solid', label = 'I+IV')
    plt.plot(H+HV, color = '#00bfaf', linestyle = 'solid', label = 'H+HV')
    plt.plot(R+RV, color = '#00FF00', linestyle = 'solid', label = 'R+RV')
    plt.plot(D, color = '#000000', linestyle = 'solid', label = 'D')
    plt.axvline(x = 15, color = '#FFD9EE', linestyle="dashed")
    plt.axvline(x = 75, color = '#FFD9EE', linestyle="dashed")

    plt.xlabel("Days")
    plt.ylabel("Proportion of people")

    plt.xlim(0, 150)
    plt.ylim(0, 0.3)

    plt.title("COVID-19 Simulation")
    plt.legend()
    plt.savefig("covidb.png")


    plt.close()




def get_infected_nodes(G):
    infected_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'I']
    return infected_nodes

def get_infected_vaccinated_nodes(G):
    infected_vaccinated_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'IV']
    return infected_vaccinated_nodes

def get_all_infected_nodes(G):
    infected_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'I' or G.nodes()[n]['state'] == 'IV']
    return infected_nodes

def get_susceptible_nodes(G):
    susceptible_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'S']
    return susceptible_nodes

def get_exposed_nodes(G):
    exposed_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'E']
    return exposed_nodes

def get_exposed_vaccinated_nodes(G):
    exposed_vaccinated_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'EV']
    return exposed_vaccinated_nodes

def get_vaccinated_nodes(G):
    vaccinated_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'V']
    return vaccinated_nodes

def get_hospitalized_nodes(G):
    hospitalized_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'H']
    return hospitalized_nodes

def get_hospitalized_vaccinated_nodes(G):
    hospitalized_vaccinated_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'HV']
    return hospitalized_vaccinated_nodes

def get_recovered_nodes(G):
    recovered_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'R']
    return recovered_nodes

def get_recovered_vaccinated_nodes(G):
    recovered_vaccinated_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'RV']
    return recovered_vaccinated_nodes

def get_dead_nodes(G):
    dead_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'D']
    return dead_nodes

def calculate_vaccine_efficiency(day):
    return 0.7*(1-day/3650)

def simulate_epidemic(G,j,k):

    for day in range(MAX_NUM_GENERATIONS):

        allInfected = get_all_infected_nodes(G)
        infectious = get_infected_nodes(G)
        infected_vaccinated = get_infected_vaccinated_nodes(G)
        exposed = get_exposed_nodes(G)
        exposed_vaccinated = get_exposed_vaccinated_nodes(G)
        susceptible = get_susceptible_nodes(G)
        vaccinated = get_vaccinated_nodes(G)
        hospitalized = get_hospitalized_nodes(G)
        hospitalized_vaccinated = get_hospitalized_vaccinated_nodes(G)
        recovered = get_recovered_nodes(G)
        recovered_vaccinated = get_recovered_vaccinated_nodes(G)

        for i in allInfected:
            neighbours = nx.neighbors(G, i)

            for n in neighbours:
                # V -> EV
                if G.nodes()[n]['state'] == 'V':
                    if np.random.uniform(0,1) < beta * (1 - calculate_vaccine_efficiency(day)):
                        G.nodes[n]['state'] = 'EV'
                # S -> E
                elif G.nodes()[n]['state'] == 'S':
                    if np.random.uniform(0,1) < beta:
                        G.nodes[n]['state'] = 'E'

        # E -> I
        for i in exposed:
            if np.random.uniform(0,1) < sigma:
                G.nodes[i]['state'] = 'I'

        # EV -> IV
        for i in exposed_vaccinated:
            if np.random.uniform(0,1) < sigma:
                G.nodes[i]['state'] = 'IV'
        # I -> H
        for i in infectious:
            if np.random.uniform(0,1) < delta3:
                G.nodes[i]['state'] = 'H'
        # IV -> HV
        for i in infected_vaccinated:
            if np.random.uniform(0,1) < delta3p:
                G.nodes[i]['state'] = 'HV'
        # I -> R
        for i in infectious:
            if np.random.uniform(0,1)< delta1:
                G.nodes[i]['state'] = 'R'
        # IV -> RV
        for i in infected_vaccinated:
            if np.random.uniform(0,1)< delta1p:
                G.nodes[i]['state'] = 'RV' 
        # I -> D
        for i in infectious:
            if np.random.uniform(0,1) < delta2:
                G.nodes[i]['state'] = 'D'
        # IV -> D
        for i in infected_vaccinated:
            if np.random.uniform(0,1) < delta2p:
                G.nodes[i]['state'] = 'D'
        # H -> R
        for i in hospitalized:
            if np.random.uniform(0,1) < zeta1:
                G.nodes[i]['state'] = 'R'
        # HV -> RV
        for i in hospitalized_vaccinated:
            if np.random.uniform(0,1) < zeta1p:
                G.nodes[i]['state'] = 'RV'
        # H -> D      
        for i in hospitalized:
            if np.random.uniform(0,1) < zeta2:
                G.nodes[i]['state'] = 'D'
        # HV -> D
        for i in hospitalized_vaccinated:
            if np.random.uniform(0,1) < zeta2p:
                G.nodes[i]['state'] = 'D'
        # R -> S
        for i in recovered:
            if np.random.uniform(0,1) < gamma:
                G.nodes[i]['state'] = 'S'
        # RV -> S
        for i in recovered_vaccinated:
            if np.random.uniform(0,1) < gammap:
                G.nodes[i]['state'] = 'S'
        # S -> V   
        for i in susceptible:
            if np.random.uniform(0,1) < p:
                G.nodes[i]['state'] = 'V'
        # V -> S      
        for i in vaccinated:
            if np.random.uniform(0,1) < eta:
                G.nodes[i]['state'] = 'S'
        
        #calculate the number of nodes in each state
        S[j].append(len([n for n in G.nodes() if G.nodes()[n]['state'] == 'S'])/N)
        V[j].append(len([n for n in G.nodes() if G.nodes()[n]['state'] == 'V'])/N)
        E[j].append(len([n for n in G.nodes() if G.nodes()[n]['state'] == 'E'])/N)
        EV[j].append(len([n for n in G.nodes() if G.nodes()[n]['state'] == 'EV'])/N)
        I[j].append(len([n for n in G.nodes() if G.nodes()[n]['state'] == 'I'])/N)
        IV[j].append(len([n for n in G.nodes() if G.nodes()[n]['state'] == 'IV'])/N)
        H[j].append(len([n for n in G.nodes() if G.nodes()[n]['state'] == 'H'])/N)
        HV[j].append(len([n for n in G.nodes() if G.nodes()[n]['state'] == 'HV'])/N)
        R[j].append(len([n for n in G.nodes() if G.nodes()[n]['state'] == 'R'])/N)
        RV[j].append(len([n for n in G.nodes() if G.nodes()[n]['state'] == 'RV'])/N)
        D[j].append(len([n for n in G.nodes() if G.nodes()[n]['state'] == 'D'])/N)

        #remove edges if lockdown is in effect
        if day == lockdownStart:
            for e in G.edges():
                if(G.nodes[e[0]]['essential'] == False and G.nodes[e[1]]['essential'] == False):
                    if(np.random.uniform(0,1) < (1 - lockdownNonEssentialContactRate)):
                        G.remove_edge(e[0], e[1])
                        removedEdges.append(e)
                else:
                    if(np.random.uniform(0,1) < (1 - lockdownEssentialContactRate)):
                        G.remove_edge(e[0], e[1])
                        removedEdges.append(e)

        #add edges back if lockdown is over
        if day == lockdownEnd:
            for e in removedEdges:
                G.add_edge(e[0], e[1])


if __name__== "__main__":

    np.random.seed(49)

    for j in range(NUM_ITERATIONS):
        G = nx.watts_strogatz_graph(N, 17, 0.01)
        initialise_graph(G)
        simulate_epidemic(G,j,0)
        print(j)
    
    S = np.mean(S, axis = 0)
    V = np.mean(V, axis = 0)
    E = np.mean(E, axis = 0)
    EV = np.mean(EV, axis = 0)
    I = np.mean(I, axis = 0)
    IV = np.mean(IV, axis = 0)
    H = np.mean(H, axis = 0)
    HV = np.mean(HV, axis = 0)
    R = np.mean(R, axis = 0)
    RV = np.mean(RV, axis = 0)
    D = np.mean(D, axis = 0)

    plot_data()

    with open('covid.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['S', 'V', 'E', 'EV', 'I', 'IV', 'H', 'HV', 'R', 'RV', 'D'])
        for i in range(MAX_NUM_GENERATIONS):
            writer.writerow([S[i], V[i], E[i], EV[i], I[i], IV[i], H[i], HV[i], R[i], RV[i], D[i]])