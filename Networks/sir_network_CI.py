import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import csv

'''
SIR model on networks confidence interval
'''

MAX_NUM_GENERATIONS = 1000
NUM_ITERATIONS = 20
S_0 = 990
I_0 = 10
R_0 = 0
N = S_0 + I_0 + R_0
I = [[[] for i in range(NUM_ITERATIONS)] for i in range(3)] 
CI = [[] for i in range(3)]
gamma = 0.006
beta = 0.0002

def initialise_graph(G):
    # initialise the graph
    for n in G.nodes():
        G.nodes[n]['state'] = 'S'

    initial_infected_node = np.random.choice(G.nodes(),I_0)
    for i in initial_infected_node:
        G.nodes[i]['state'] = 'I'


def plot_graph(G,step):
    """
    Read a graph from a file, each node is labelled according to its state.
    """
    my_labels = {}
    for n in G.nodes():
        my_labels[n] = G.nodes()[n]['state']

    # TODO: may want to colour each noded according to its state, e.g S=green, I=red, R=grey.
    state_color_map = {'S': 'green', 'I': 'red', 'R': 'blue'}
    node_colors = [state_color_map[G.nodes[n]['state']] for n in G.nodes()]
    nx.draw(G, with_labels = True, labels=my_labels, node_color=node_colors)

    #plt.show()
    plt.savefig("diseases_on_networkѕ_step_"+str(step)+".png")
    plt.close()

def plot_data():
    # Erdo-Renyi Network
    #fig, ax = plt.subplots()
    #ax.plot(I[0], color = 'red', linestyle = 'solid', label = 'p = 0.2')
    #ax.fill_between(range(1000), (I[0]-CI[0]), (I[0]+CI[0]), color = 'red', alpha = 0.2)
    #ax.plot(I[1], color = 'blue', linestyle = 'solid', label = 'p = 0.5')
    #ax.fill_between(range(1000), (I[1]-CI[1]), (I[1]+CI[1]), color = 'blue', alpha = 0.2)
    #ax.plot(I[2], color = 'green', linestyle = 'solid', label = 'p = 0.8')
    #ax.fill_between(range(1000), (I[2]-CI[2]), (I[2]+CI[2]), color = 'green', alpha = 0.2)
    #ax.set_ylim(0,1.1)
    #plt.title("Erdos-Renyi Network")
    #plt.ylabel('Population')
    #plt.xlabel('Time')
    #plt.legend()
    #plt.savefig("ERCI.png")
    #plt.close()
    
    # Watts-strogatz Network
    fig , ax = plt.subplots()
    ax.plot(I[0], color = 'red', linestyle = 'solid', label = 'k = 50')
    ax.fill_between(range(1000), (I[0]-CI[0]), (I[0]+CI[0]), color = 'red', alpha = 0.2)
    ax.plot(I[1], color = 'blue', linestyle = 'solid', label = 'k = 100')
    ax.fill_between(range(1000), (I[1]-CI[1]), (I[1]+CI[1]), color = 'blue', alpha = 0.2)
    ax.plot(I[2], color = 'green', linestyle = 'solid', label = 'k = 300')
    ax.fill_between(range(1000), (I[2]-CI[2]), (I[2]+CI[2]), color = 'green', alpha = 0.2)
    ax.set_ylim(0,1.1)
    plt.title("Watts-Strogatz Network")
    plt.ylabel('Population')
    plt.xlabel('Time')
    plt.legend()
    plt.savefig("WSCI.png")
    plt.close()
#
    ## Barabasi-Albert Network
    #fig , ax = plt.subplots()
    #ax.plot(I[0], color = 'red', linestyle = 'solid', label = 'm = 50')
    #ax.fill_between(range(1000), (I[0]-CI[0]), (I[0]+CI[0]), color = 'red', alpha = 0.2)
    #ax.plot(I[1], color = 'blue', linestyle = 'solid', label = 'm = 100')
    #ax.fill_between(range(1000), (I[1]-CI[1]), (I[1]+CI[1]), color = 'blue', alpha = 0.2)
    #ax.plot(I[2], color = 'green', linestyle = 'solid', label = 'm = 300')
    #ax.fill_between(range(1000), (I[2]-CI[2]), (I[2]+CI[2]), color = 'green', alpha = 0.2)
    #ax.set_ylim(0,1.1)
    #plt.title("Barabasi-Albert Network")
    #plt.legend()
    #plt.ylabel('Population')
    #plt.xlabel('Time')
    #plt.savefig("BACI.png")
    #plt.close()


def get_infected_nodes(G):
    """
    Useful function to get a list of all the infected nodes, i.e. all those whos state is 'I'
    """
    infected_nodes = [n for n in G.nodes() if G.nodes()[n]['state'] == 'I']
    return infected_nodes



def simulate_epidemic(G,j,k):

    for step in range(MAX_NUM_GENERATIONS):

        infectious = get_infected_nodes(G)

        for i in infectious:
            neighbours = nx.neighbors(G, i)
            susceptibles = [n for n in neighbours if G.nodes()[n]['state'] == 'S']

            for s in susceptibles:
                if np.random.uniform(0,1) < beta:
                    G.nodes[s]['state'] = 'I'

        for i in infectious:
            if np.random.uniform(0,1)< gamma:
                G.nodes[i]['state'] = 'R'

        
        #calculate the number of infected nodes
        I[j][k].append(len([n for n in G.nodes() if G.nodes()[n]['state'] == 'I'])/N)




if __name__== "__main__":

    probErdos = [0.2, 0.5, 0.8]
    watts = [50,100,300]
    mBarabasi = [50,100,300]

    for j, p in enumerate(probErdos):
        for k in range(NUM_ITERATIONS):
            G = nx.erdos_renyi_graph(N, p)    
            initialise_graph(G)
            simulate_epidemic(G,j,k)
        CI[j]=1.96*np.std(I[j])/np.sqrt(NUM_ITERATIONS)
        I[j] = np.mean(I[j], axis = 0)

    for j, w in enumerate(watts):
        for k in range(NUM_ITERATIONS):
            G = nx.watts_strogatz_graph(N, w, 0.01)
            initialise_graph(G)
            simulate_epidemic(G,j,k)
        CI[j]=1.96*np.std(I[j])/np.sqrt(NUM_ITERATIONS)
        I[j] = np.mean(I[j], axis = 0)

    for j, m in enumerate(mBarabasi):
        for k in range(NUM_ITERATIONS):
            G = nx.barabasi_albert_graph(N,m)
            initialise_graph(G)
            simulate_epidemic(G,j,k)
        CI[j]=1.96*np.std(I[j])/np.sqrt(NUM_ITERATIONS)
        I[j] = np.mean(I[j], axis = 0)

    plot_data()

    with open('data.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['I'])
        for i in range(MAX_NUM_GENERATIONS):
            writer.writerow([I[0][i]])

    with open('dataCI.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['I'])
        for i in range(MAX_NUM_GENERATIONS):
            writer.writerow([CI[0]])
    
    with open('data1.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['I'])
        for i in range(MAX_NUM_GENERATIONS):
            writer.writerow([I[1][i]])

    with open('data1CI.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['I'])
        for i in range(MAX_NUM_GENERATIONS):
            writer.writerow([CI[1]])

    with open('data2.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['I'])
        for i in range(MAX_NUM_GENERATIONS):
            writer.writerow([I[2][i]])

    with open('data2CI.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['I'])
        for i in range(MAX_NUM_GENERATIONS):
            writer.writerow([CI[2]])