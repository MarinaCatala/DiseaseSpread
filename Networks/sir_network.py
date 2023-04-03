import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import csv

'''
SIR model on networks
'''

MAX_NUM_GENERATIONS = 4000
NUM_ITERATIONS = 20
S_0 = 990
I_0 = 10
R_0 = 0
N = S_0 + I_0 + R_0
S = [[[] for i in range(NUM_ITERATIONS)] for i in range(3)]
I = [[[] for i in range(NUM_ITERATIONS)] for i in range(3)]
R = [[[] for i in range(NUM_ITERATIONS)] for i in range(3)]
gamma = 0.006
beta = 0.0002


def initialise_graph(G):
    # initialise the graph
    for n in G.nodes():
        G.nodes[n]['state'] = 'S'

    # randomly pick a node as an initial infection (or always pick node 1 !)
    initial_infected_node = np.random.choice(G.nodes(), I_0)
    for i in initial_infected_node:
        G.nodes[i]['state'] = 'I'


def plot_graph(G, step):
    """
    Read a graph from a file, each node is labelled according to its state.
    """
    my_labels = {}
    for n in G.nodes():
        my_labels[n] = G.nodes()[n]['state']

    # TODO: may want to colour each noded according to its state, e.g S=green, I=red, R=grey.
    state_color_map = {'S': 'green', 'I': 'red', 'R': 'blue'}
    node_colors = [state_color_map[G.nodes[n]['state']] for n in G.nodes()]
    nx.draw(G, with_labels=True, labels=my_labels, node_color=node_colors)

    # plt.show()
    plt.savefig("diseases_on_networkѕ_step_"+str(step)+".png")
    plt.close()


def plot_data():
    # Erdo-Renyi Network
    #plt.plot(S[0], color='blue', linestyle='solid', label='p = 0.2')
    #plt.plot(I[0], color='red', linestyle='solid')
    #plt.plot(R[0], color='green', linestyle='solid')
    #plt.plot(S[1], color='blue', linestyle='dashed', label='p = 0.5')
    #plt.plot(I[1], color='red', linestyle='dashed')
    #plt.plot(R[1], color='green', linestyle='dashed')
    #plt.plot(S[2], color='blue', linestyle='dotted', label='p = 0.8')
    #plt.plot(I[2], color='red', linestyle='dotted')
    #plt.plot(R[2], color='green', linestyle='dotted')
##
    #plt.xlim(0, 1000)
    #plt.title("Erdos-Renyi Network")
    #plt.ylabel('Population')
    #plt.xlabel('Time')
    #plt.legend()
    #plt.savefig("ER1.png")
    #plt.close()

    # Watts-Strogatz Network
    plt.plot(S[0], color = 'blue', linestyle = 'solid', label = 'k = 50')
    plt.plot(I[0], color = 'red', linestyle = 'solid')
    plt.plot(R[0], color = 'green', linestyle = 'solid')
    plt.plot(S[1], color = 'blue', linestyle = 'dashed', label = 'k = 100')
    plt.plot(I[1], color = 'red', linestyle = 'dashed')
    plt.plot(R[1], color = 'green', linestyle = 'dashed')
    plt.plot(S[2], color = 'blue', linestyle = 'dotted', label = 'k = 300')
    plt.plot(I[2], color = 'red', linestyle = 'dotted')
    plt.plot(R[2], color = 'green', linestyle = 'dotted')
##
    plt.xlim(0, 1000)
    plt.title("Watts-Strogatz Network")
    plt.ylabel('Population')
    plt.xlabel('Time')
    plt.legend()
    plt.savefig("WS.png")
    plt.close()
#
    # Barabasi-Albert Network
    # plt.plot(S[0], color = 'blue', linestyle = 'solid', label = 'm = 50')
    # plt.plot(I[0], color = 'red', linestyle = 'solid')
    # plt.plot(R[0], color = 'green', linestyle = 'solid')
    # plt.plot(S[1], color = 'blue', linestyle = 'dashed', label = 'm = 100')
    # plt.plot(I[1], color = 'red', linestyle = 'dashed')
    # plt.plot(R[1], color = 'green', linestyle = 'dashed')
    # plt.plot(S[2], color = 'blue', linestyle = 'dotted', label = 'm = 300')
    # plt.plot(I[2], color = 'red', linestyle = 'dotted')
    # plt.plot(R[2], color = 'green', linestyle = 'dotted')
##
    # plt.xlim(0, 1000)
    # plt.title("Barabasi-Albert Network")
    # plt.ylabel('Population')
    # plt.xlabel('Time')
    # plt.legend()
    # plt.savefig("BA.png")
    # plt.close()


def get_infected_nodes(G):
    """
    Useful function to get a list of all the infected nodes, i.e. all those whos state is 'I'
    """
    infected_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'I']
    return infected_nodes


def simulate_epidemic(G, j, k):

    for step in range(MAX_NUM_GENERATIONS):

        infectious = get_infected_nodes(G)

        # 1. Determine the S->I transitions
        for i in infectious:
            # let's get the neighbours and filter out those who are susceptible
            neighbours = nx.neighbors(G, i)
            susceptibles = [n for n in neighbours if G.nodes()[n]
                            ['state'] == 'S']

            for s in susceptibles:
                # TODO do the s->i infection, marking s now as state 'I'
                if np.random.uniform(0, 1) < beta:
                    G.nodes[s]['state'] = 'I'

        # 2. Determine the I->R transitions
        for i in infectious:
            # TODO do the i->r infection, marking s now as state 'R'
            if np.random.uniform(0, 1) < gamma:
                G.nodes[i]['state'] = 'R'


        # calculate the number of S, I, R
        S[j][k].append(
            len([n for n in G.nodes() if G.nodes()[n]['state'] == 'S'])/N)
        I[j][k].append(
            len([n for n in G.nodes() if G.nodes()[n]['state'] == 'I'])/N)
        R[j][k].append(
            len([n for n in G.nodes() if G.nodes()[n]['state'] == 'R'])/N)


if __name__ == "__main__":

    probErdos = [0.2, 0.5, 0.8]
    watts = [50, 100, 300]
    mBarabasi = [50, 100, 300]

    for j, p in enumerate(probErdos):
        for k in range(NUM_ITERATIONS):
            G = nx.erdos_renyi_graph(N, p)
            initialise_graph(G)
            simulate_epidemic(G, j, k)
            print(j, k)
        S[j] = np.mean(S[j], axis=0)
        I[j] = np.mean(I[j], axis=0)
        R[j] = np.mean(R[j], axis=0)

    for j, w in enumerate(watts):
       for k in range(NUM_ITERATIONS):
           G = nx.watts_strogatz_graph(N,w,0.01)
           initialise_graph(G)
           simulate_epidemic(G,j,k)
           print(j,k)
       S[j] = np.mean(S[j], axis = 0)
       I[j] = np.mean(I[j], axis = 0)
       R[j] = np.mean(R[j], axis = 0)


    for j, m in enumerate(mBarabasi):
       for k in range(NUM_ITERATIONS):
           G = nx.barabasi_albert_graph(N,m)
           initialise_graph(G)
           simulate_epidemic(G,j,k)
       S[j] = np.mean(S[j], axis = 0)
       I[j] = np.mean(I[j], axis = 0)
       R[j] = np.mean(R[j], axis = 0)

    plot_data()

    with open('data.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['S', 'I', 'R'])
        for i in range(MAX_NUM_GENERATIONS):
            writer.writerow([S[0][i], I[0][i], R[0][i]])

    with open('data1.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['S', 'I', 'R'])
        for i in range(MAX_NUM_GENERATIONS):
            writer.writerow([S[1][i], I[1][i], R[1][i]])

    with open('data2.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['S', 'I', 'R'])
        for i in range(MAX_NUM_GENERATIONS):
            writer.writerow([S[2][i], I[2][i], R[2][i]])
