import sys
import argparse
import time
from Functions.construct_mutation_matrix import *
from Functions.construct_bipartite_graph import *
from Functions.LBSA import *
from Functions.results import *
from Functions.ARO import *



def main():
    #Here we read the necessary parameters required for running the algorithm
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, \
                                    usage='\n\npython main.py <-n network file> <-d dataset name> <-a algorithm name> [-c cooling factor] [-i number of iterations]\n', \
                                    description='', epilog='Implementation of MYTHESIS.\n')

    parser.add_argument('-n', '--network', required=True, type=str, help='mutations file.', metavar='', dest="network")
    parser.add_argument('-d', '--dataset', required=True, type=str, help='dataset name.', metavar='', dest="dataset_name")
    parser.add_argument('-a', '--algorithm', required=True, default='ARO', type=str, help='algorithm name.', metavar='', dest="algorithm_name")
    parser.add_argument('-c', '--cooling', default=(1 - (10 ** -2)), type=int, help='Cooling factor', metavar='', dest="cooling")
    parser.add_argument('-i', '--iterations', default=(10 ** 5), type=int, help='Number of iterations', metavar='', dest="iterations")

    args = parser.parse_args() #all input arguments are parsed into args variable

    #get input data
    maf_file_name = args.network
    dataset_name = args.dataset_name
    algorithm_name = args.algorithm_name
    cooling_fact = args.cooling
    iterations = args.iterations

    #read data
    print ("START RUN...\n Reading input data...")
    mutation_matrix = get_mutation_matrix_from_maf("Data/TCGA/" + maf_file_name)
    influence_graph = nx.read_edgelist("Data/Reactome_FIsInGene_2021.txt", delimiter='\t')
    gene_expression_matrix = pd.read_csv('Data/GBMPatientOutlierMatrix.csv', index_col=0)
    print ("Finished Reading input.")

    #create bipartite graph
    print ("Creating Bipartite Graph. Start...")
    bipartite_graph, green_nodes = create_bipartite_graph(mutation_matrix, influence_graph, gene_expression_matrix)
    print ("Bipartite Graph Created Successfully.")

    drivernet_starttime = time.time()
    if (algorithm_name == 'LBSA'):
        # implement List Based Simulated Annealing
        print("Apply List Based Simulated Annealing")
        drivers_list, drivers_order_list2 = List_Based_Simulated_Annaeling(dataset_name, green_nodes,
                                                                                bipartite_graph, cooling_fact,
                                                                                iterations)
        print("Algorithm successful...")
    elif (algorithm_name == 'ARO'):
        # implement Artificial Rabbits Optimization (ARO)
        # fun_index = 1: Sphere
        # fun_index = 2: Schwefel 2.22
        # fun_index = 3: Schwefel 1.2
        # fun_index = 4: Schwefel 2.21
        # fun_index = 5: Rosenbrock

        fun_index = 2
        max_it = iterations
        print("Apply Artificial Rabbits Optimization")
        drivers_list, drivers_order_list2 = ARO(fun_index, max_it, bipartite_graph, green_nodes, dataset_name)
        print("Algorithm successful...")
    else:
        print("Enter valid algorithm name.")
        return

    drivernet_endtime = time.time()
    #write output
    print ("Writing Data to output...")
    report(drivers_list, drivers_order_list2, dataset_name)
    print ("Output file -- ranked_driver_genes_" + dataset_name + ".txt -- created.\n Check Output Folder")
    print ("Cheers!!!")

    total_algorithm_runtime = drivernet_endtime - drivernet_starttime
    total_algorithm_runtime_mins = total_algorithm_runtime / 60
    print("---Total program runtime ----%s minutes ----" % (total_algorithm_runtime_mins))

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print ("Error Encountered. EXITING....")
        sys.exit(0)

