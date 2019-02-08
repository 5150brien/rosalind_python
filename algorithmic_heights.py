
def fibonacci_numbers(n):
    """
    Given an integer n, returns the value of the nth (0-based) Fibonacci number
    """
    fibo_sequence = []
    for x in range(0, n + 1):
        if x == 0:
            fibo_sequence.append(0)
        elif x == 1:
            fibo_sequence.append(1)
        else:
            fibo_sequence.append(
                fibo_sequence[-1] + fibo_sequence[-2]
            )

    return fibo_sequence[-1]

def binary_search(reference_list, value_list):
    """
    Finds the index of value_list entries in reference list

    Values that are not found in reference_list get an index of -1

    :param reference_list: a sorted list of integers
    :type reference_list: list
    :param value_list: a list of any integers in any order
    :type value_list: list
    :rtype: list
    :return: the indexes where value_list entries are found in reference_list
    """
    indexes = []
    n = len(reference_list)

    for item in value_list:
        lower_boundary = 0
        upper_boundary = n - 1
        index = -1

        while lower_boundary <= upper_boundary and index < 0:
            middle = round((lower_boundary + upper_boundary) / 2)
            if reference_list[middle] < item:
                lower_boundary = middle + 1
            elif reference_list[middle] > item:
                upper_boundary = middle - 1
            else:
                index = middle + 1

        indexes.append(index)

    return indexes

def degree_array(edge_file_path):
    """
    Returns a sorted list of the degrees for each vertex in an edge file

    :param edge_file_path: the path to a file of data in edge list format
    :type edge_file_path: str
    :rtype: list
    :return: a list of degrees for each vertex specified in edge_file_path
    """
    degree_counts = []

    with open(edge_file_path, "r") as edge_data:
        for i, line in enumerate(edge_data):
            edge = line.strip().split(" ")
            if i == 0:
                # First line is vertex & edge counts (not an edge)
                total_vertices = int(edge[0])
                degree_counts = [0 for x in range(total_vertices)]
            else:
                for vertex in edge:
                    degree_counts[int(vertex) - 1] += 1

    return degree_counts

def insertion_sort(length, sort_list):
    """
    Finds the total number of swaps required to complete an insertion sort

    :param length: the number of items in sort_list
    :type length: int
    :param sort_list: an unsorted list of numbers
    :type sort_list: list
    :rtype: int
    :return: the number of swaps required to complete an insertion sort
    """
    swaps = 0
    for i in range(1, length):
        k = i
        while k > 0 and sort_list[k] < sort_list[k - 1]:
            sort_list[k - 1], sort_list[k] = sort_list[k], sort_list[k - 1]
            swaps += 1
            k -= 1

    print(sort_list)
    return swaps

def double_degree_array(edge_file_path):
    """
    Calculates the degrees for all neighbors of each vertex in an edge file

    :param edge_file_path: the path to a file of data in edge list format
    :type edge_file_path: str
    :rtype: list
    :return: a list of total degrees for all neighbors of each vertex
    """
    neighbor_degrees = []

    with open(edge_file_path, "r") as edge_data:
        # Count the degrees and neighbors for each vertex
        for i, line in enumerate(edge_data):
            edge = line.strip().split(" ")
            if i == 0:
                length = int(edge[0])
                degrees = [0 for x in range(length)]
                neighbors = [[] for x in range(length)]
            else:
                a, b = int(edge[0]), int(edge[1])
                degrees[a - 1] += 1
                degrees[b - 1] += 1
                neighbors[a - 1].append(int(edge[1]))
                neighbors[b - 1].append(int(edge[0]))

    # Sum the degrees of all neighbors for each vertex
    for i, neighbor_group in enumerate(neighbors):
        total = 0
        for neighbor in neighbor_group:
            total += degrees[neighbor - 1]
        neighbor_degrees.append(total)

    return neighbor_degrees
