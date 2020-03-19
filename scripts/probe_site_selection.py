import numpy as np
import itertools

def conseq_point_dist(array):
    """Finds the difference (distance) between consequtive points in an array"""
    first_probe_site = []
    dist_between_probe_sites = np.append(first_probe_site, [t - s for s, t in zip(array[:], array[1:])])
    return dist_between_probe_sites

def func_ideal_length(target_position, no_probe_sites):
    """Finds the optimum length between insertion sites, based on the nucleotide length"""
    length = target_position[-1] - target_position[0] #length of sequence available for probe sites
    no_intervals = no_probe_sites - 1 #number of intervals between probe sites
    ideal_length = length / no_intervals #if positions of the sites available for probes are continuous
    return ideal_length

def my_func_dist(array, ideal_length):
    """Finds the euclidean distance between the optimum length and the insertion site length"""
    square_of_dist_from_ideal_per_point = [(ideal_length - t)*(ideal_length - t) for t in array] #1st line of crit
    square_of_dist_from_ideal_SUM = np.sum(square_of_dist_from_ideal_per_point, axis=1) #2nd line of criterion
    return square_of_dist_from_ideal_SUM

def choose_combination(target_position, no_probe_sites):
    # all possible combinations for probe site insertions
    all_combinations = list(itertools.combinations(target_position, no_probe_sites))

    # use this to specify sites of probe insertions
    # e.g. if you want a probe site at the first available candidate and the last available candidate
    selected_combinations =[]
    #last_element = len(target_position) - 1
    for x in all_combinations:
        if x[0] == target_position[0] and x[-1] == target_position[-1]: #I specify the first and the last probe sites
            selected_combinations.append(x)

    # make sure selected_combinations is an array
    selected_combinations = np.asarray(selected_combinations)


    # Applies the above function to every row of the selected combinations array
    # Finds the distance between selected insertion sites for every row
    dist_between_probe_sites_array = np.apply_along_axis(conseq_point_dist, 1, selected_combinations)

    ideal_length = func_ideal_length(target_position, no_probe_sites)
    min_distance = np.min(my_func_dist(dist_between_probe_sites_array, ideal_length)) #Finding the distance that minimises criterion
    loc= np.where(my_func_dist(dist_between_probe_sites_array, ideal_length) == np.min(min_distance)) #Index of minimum
    return selected_combinations[loc[0][0]]
    #print(selected_combinations[loc[0][0]])

def test():

    test_targets = [
    (3, 7, 9.5, 13, 19, 22, 23, 24, 25, 27, 30, 34),
    (1, 8, 10, 13, 23, 25, 82, 92, 96, 97, 98)
    ]
    #target_position = (3, 7, 9.5, 13, 19, 22, 23, 24, 25, 27, 30, 34) #list for coordinates of target positions
    no_probe_sites = 5
    for t in test_targets:
        print('options', t)
        comb = choose_combination(t, no_probe_sites)
        print('chosen', comb)

if __name__ == '__main__':
    test()
