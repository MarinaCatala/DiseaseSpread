import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import csv

'''
COVID-19 Total Infected Simulation
'''

MAX_NUM_GENERATIONS = 1500
NUM_ITERATIONS = 20
S_0 = 9900
E_0 = 100
N = S_0 + E_0
I = [[] for i in range(NUM_ITERATIONS)]
I2 = [[] for i in range(NUM_ITERATIONS)]
I3 = [[] for i in range(NUM_ITERATIONS)]
I4 = [[] for i in range(NUM_ITERATIONS)]


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
zeta1 = 0.018
zeta1p = 0.018
zeta2 = 1/14
zeta2p = 1/14


lockdownEssentialContactRate = 0.5
lockdownNonEssentialContactRate = 0.1
removedEdges = []


rateEssential = 0.15


def initialise_graph(G):
    # initialise the graph

    for n in G.nodes():
        G.nodes[n]['state'] = 'S'
        G.nodes[n]['alreadyInfected'] = False
        if (np.random.random() < rateEssential):
            G.nodes[n]['essential'] = True
        else:
            G.nodes[n]['essential'] = False

    # randomly pick a node as an initial infection (or always pick node 1 !)
    initial_infected_node = np.random.choice(G.nodes(), E_0)

    for i in initial_infected_node:
        G.nodes[i]['state'] = 'E'


def plot_graph(G, step):
    """
    Read a graph from a file, each node is labelled according to its state.
    """
    my_labels = {}
    for n in G.nodes():
        my_labels[n] = n

    # TODO: may want to colour each noded according to its state, e.g S=green, I=red, R=grey.
    state_color_map = {'S': 'green', 'E': 'orange', 'I': 'red', 'R': 'blue'}
    node_colors = [state_color_map[G.nodes[n]['state']] for n in G.nodes()]
    nx.draw(G, with_labels=True, labels=my_labels, node_color=node_colors)

    # plt.show()
    plt.savefig("diseases_on_networkѕ_step_"+str(step)+".png")
    plt.close()


def plot_data():

    plt.figure()

    plt.plot(I, color='#0000FF', linestyle='solid',
             label='Without Vaccination and Lockdown')
    plt.plot(I2, color='#00FFFF', linestyle='solid', label='With Vaccination')
    plt.plot(I3, color='#FF8000', linestyle='solid', label='With Lockdown')
    plt.plot(I4, color='#FF0000', linestyle='solid',
             label='With Vaccination and Lockdown')

    plt.xlabel("Days")
    plt.ylabel("Proportion of people")
    plt.xlim(0, 365)
    plt.ylim(0, 1.2)

    plt.title("COVID-19 Simulation")
    plt.legend()
    plt.savefig("covidTotalInfected.png")

    plt.close()


def get_infected_nodes(G):
    infected_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'I']
    return infected_nodes


def get_infected_vaccinated_nodes(G):
    infected_vaccinated_nodes = [n for n in G.nodes() if G.nodes()[
        n]['state'] == 'IV']
    return infected_vaccinated_nodes


def get_all_infected_nodes(G):
    infected_nodes = [n for n in G.nodes() if G.nodes(
    )[n]['state'] == 'I' or G.nodes()[n]['state'] == 'IV']
    return infected_nodes


def get_susceptible_nodes(G):
    susceptible_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'S']
    return susceptible_nodes


def get_exposed_nodes(G):
    exposed_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'E']
    return exposed_nodes


def get_exposed_vaccinated_nodes(G):
    exposed_vaccinated_nodes = [n for n in G.nodes() if G.nodes()[
        n]['state'] == 'EV']
    return exposed_vaccinated_nodes


def get_vaccinated_nodes(G):
    vaccinated_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'V']
    return vaccinated_nodes


def get_hospitalized_nodes(G):
    hospitalized_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'H']
    return hospitalized_nodes


def get_hospitalized_vaccinated_nodes(G):
    hospitalized_vaccinated_nodes = [
        n for n in G.nodes() if G.nodes()[n]['state'] == 'HV']
    return hospitalized_vaccinated_nodes


def get_recovered_nodes(G):
    recovered_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'R']
    return recovered_nodes


def get_recovered_vaccinated_nodes(G):
    recovered_vaccinated_nodes = [n for n in G.nodes() if G.nodes()[
        n]['state'] == 'RV']
    return recovered_vaccinated_nodes


def get_dead_nodes(G):
    dead_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'D']
    return dead_nodes


def calculate_vaccine_efficiency(day):
    return 0.7*(1-day/3650)


def simulate_epidemic(G, j, k):

    numTotalInfected = 0

    if k == 0:
        p = 0
        lockdownPeriod = -1
        lockdownStart = -1
    elif k == 1:
        p = 0.05
        lockdownPeriod = -1
        lockdownStart = -1
    elif k == 2:
        p = 0
        lockdownPeriod = 60
        lockdownStart = 15
    elif k == 3:
        p = 0.05
        lockdownPeriod = 60
        lockdownStart = 15

    lockdownEnd = lockdownStart + lockdownPeriod

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
                    if np.random.uniform(0, 1) < beta * (1 - calculate_vaccine_efficiency(day)):
                        G.nodes[n]['state'] = 'EV'
                # S -> E
                elif G.nodes()[n]['state'] == 'S':
                    if np.random.uniform(0, 1) < beta:
                        G.nodes[n]['state'] = 'E'

        # E -> I
        for i in exposed:
            if np.random.uniform(0, 1) < sigma:
                G.nodes[i]['state'] = 'I'
                if (G.nodes[i]['alreadyInfected'] == False):
                    numTotalInfected += 1
                    G.nodes[i]['alreadyInfected'] = True

        for i in exposed_vaccinated:
            if np.random.uniform(0, 1) < sigma:
                G.nodes[i]['state'] = 'IV'
                if (G.nodes[i]['alreadyInfected'] == False):
                    numTotalInfected += 1
                    G.nodes[i]['alreadyInfected'] = True

        for i in infectious:
            if np.random.uniform(0, 1) < delta3:
                G.nodes[i]['state'] = 'H'

        for i in infected_vaccinated:
            if np.random.uniform(0, 1) < delta3p:
                G.nodes[i]['state'] = 'HV'

        for i in infectious:
            if np.random.uniform(0, 1) < delta1:
                G.nodes[i]['state'] = 'R'

        for i in infected_vaccinated:
            if np.random.uniform(0, 1) < delta1p:
                G.nodes[i]['state'] = 'RV'

        for i in infectious:
            if np.random.uniform(0, 1) < delta2:
                G.nodes[i]['state'] = 'D'

        for i in infected_vaccinated:
            if np.random.uniform(0, 1) < delta2p:
                G.nodes[i]['state'] = 'D'

        for i in hospitalized:
            if np.random.uniform(0, 1) < zeta2:
                G.nodes[i]['state'] = 'R'

        for i in hospitalized_vaccinated:
            if np.random.uniform(0, 1) < zeta2p:
                G.nodes[i]['state'] = 'RV'

        for i in hospitalized:
            if np.random.uniform(0, 1) < zeta1:
                G.nodes[i]['state'] = 'D'

        for i in hospitalized_vaccinated:
            if np.random.uniform(0, 1) < zeta1p:
                G.nodes[i]['state'] = 'D'

        for i in recovered:
            if np.random.uniform(0, 1) < gamma:
                G.nodes[i]['state'] = 'S'

        for i in recovered_vaccinated:
            if np.random.uniform(0, 1) < gammap:
                G.nodes[i]['state'] = 'S'

        for i in susceptible:
            if np.random.uniform(0, 1) < p:
                G.nodes[i]['state'] = 'V'

        for i in vaccinated:
            if np.random.uniform(0, 1) < eta:
                G.nodes[i]['state'] = 'S'

        if k == 0:
            I[j].append(numTotalInfected/N)
        elif k == 1:
            I2[j].append(numTotalInfected/N)
        elif k == 2:
            I3[j].append(numTotalInfected/N)
        elif k == 3:
            I4[j].append(numTotalInfected/N)

        if day == lockdownStart:
            for e in G.edges():
                if (G.nodes[e[0]]['essential'] == False and G.nodes[e[1]]['essential'] == False):
                    if (np.random.uniform(0, 1) < (1 - lockdownNonEssentialContactRate)):
                        G.remove_edge(e[0], e[1])
                        removedEdges.append(e)
                else:
                    if (np.random.uniform(0, 1) < (1 - lockdownEssentialContactRate)):
                        G.remove_edge(e[0], e[1])
                        removedEdges.append(e)

        if day == lockdownEnd:
            for e in removedEdges:
                G.add_edge(e[0], e[1])


if __name__ == "__main__":

    nodeDegrees = []

    for j in range(NUM_ITERATIONS):
        G = nx.watts_strogatz_graph(N, 17, 0.01)
        initialise_graph(G)
        print(nx.is_connected(G))
        # for node in G.nodes():
        #     nodeDegrees.append(G.degree(node))
        simulate_epidemic(G, j, 0)
        print(j)

    # print(np.mean(nodeDegrees))
    # print(np.min(nodeDegrees))
    # print(np.max(nodeDegrees))

    nodeDegrees = []


    I = np.mean(I, axis=0)

    for j in range(NUM_ITERATIONS):
        G = nx.watts_strogatz_graph(N, 17, 0.01)
        initialise_graph(G)
        print(nx.is_connected(G))
        # for node in G.nodes():
        #     nodeDegrees.append(G.degree(node))
        simulate_epidemic(G, j, 1)
        print(j)

    # print(np.mean(nodeDegrees))
    # print(np.min(nodeDegrees))
    # print(np.max(nodeDegrees))

    nodeDegrees = []


    I2 = np.mean(I2, axis=0)

    for j in range(NUM_ITERATIONS):
        G = nx.watts_strogatz_graph(N, 17, 0.01)
        initialise_graph(G)
        print(nx.is_connected(G))
        # for node in G.nodes():
        #     nodeDegrees.append(G.degree(node))
        simulate_epidemic(G, j, 2)
        print(j)

    # print(np.mean(nodeDegrees))
    # print(np.min(nodeDegrees))
    # print(np.max(nodeDegrees))

    nodeDegrees = []

    I3 = np.mean(I3, axis=0)

    for j in range(NUM_ITERATIONS):
        G = nx.watts_strogatz_graph(N, 17, 0.01)
        initialise_graph(G)
        print(nx.is_connected(G))
        for node in G.nodes():
            nodeDegrees.append(G.degree(node))
        simulate_epidemic(G, j, 3)
        print(j)

    # print(np.mean(nodeDegrees))
    # print(np.min(nodeDegrees))
    # print(np.max(nodeDegrees))

    I4 = np.mean(I4, axis=0)

    plot_data()
#
    with open('dataTI.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['I', 'I2', 'I3', 'I4'])
        for i in range(MAX_NUM_GENERATIONS):
            writer.writerow([I[i], I2[i], I3[i], I4[i]])
