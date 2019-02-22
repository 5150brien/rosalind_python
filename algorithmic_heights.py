
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

def majority_element(value_list):
    """
    Returns the majority element for value_list or -1 if there is none.

    Majority element in a list of length n is a value that occurs
    more than n/2 times. Implements Moore's Voting Algorithm.

    :param value_list: a list of positive integers
    :type value_list: list
    :rtype: int
    :return: the majority element or -1 if none exists
    """
    candidate = value_list[0]
    count = 1

    for val in value_list:
        if candidate == val:
            count += 1
        else:
            count -= 1

        # majority element occurrences minus all others must be > 0
        if count == 0:
            candidate, count = val, 1

    if value_list.count(candidate) > len(value_list) / 2:
        return candidate
    else:
        return -1

def merge_two_sorted_arrays(a, b):
    """
    Combines sorted lists a and b into a single, sorted list
    """
    c = []
    i = 0
    j = 0

    while i < len(a) and j < len(b):
        if a[i] < b[j]:
            c.append(a[i])
            i += 1
        else:
            c.append(b[j])
            j += 1

    if i < len(a):
        c = c + a[i:]

    if j < len(b):
        c = c + b[j:]

    return c

def load_arrays(file_path):
    """
    Loads data for 'merge two sorted arrays' problem into two lists
    """
    with open(file_path, "r") as data:
        for i, line in enumerate(data):
            if i == 1:
                a = [int(x) for x in line.strip().split(" ")]
            if i == 3:
                b = [int(x) for x in line.strip().split(" ")]

        return a, b

def two_sum(file_path):
    """
    Prints p, q such that A[p] = -A[q] in each of k lists of length n

    Only one match of A[p] = -A[q] should be printed per list/line. If no
    match is found, -1 is printed for the line.
    
    The data file format will contain two integers k and n, the total
    number of lists and the length of each list, on the first line and a
    list of n integers on each subsequent line.

    :param file_path: the path to a file of two_sum problem data
    :type file_path: str
    :rtype: None
    :return: None
    """
    with open(file_path, "r") as data:
        for i, line in enumerate(data):
            values = line.strip().split()
            if i == 0:
                # List length, n
                n = int(values[1])
            else:
                inventory = {}
                solved = False
                for x in range(0, n):
                    val = int(values[x])
                    if -val in inventory and not solved:
                        p = inventory[-val] + 1
                        q = x + 1
                        solved = True
                        print(p, q)
                    else:
                        inventory[val] = x
                if not solved:
                    print("-1")

def parse_edge_file(file_path):
    """
    Returns graph data, given a file in edge file format

    :param file_path: the path to an edge file
    :type file_path: str
    :rtype: tuple
    :return: a 3-tuple with (1) the total number of vertices, (2) the total
             number of edges, and (3) a list of edges, each of which is also a
             (2-item) list
    """
    total_vertices = 0
    total_edges = 0
    edge_list = []

    with open(file_path, "r") as data:
        for i, line in enumerate(data):
            values = line.strip().split(" ")
            if i == 0:
                total_vertices = int(values[0])
                total_edges = int(values[1])
            else:
                u, v = int(values[0]), int(values[1])
                edge_list.append([u, v])

    return total_vertices, total_edges, edge_list

def breadth_first_search(file_path, starting_point=1):
    """
    Returns a list of shortest path distances for a graph defined by an edge file
    """
    queue = []
    vertices, _, edges = parse_edge_file(file_path)
    distances = [-1 for x in range(0, vertices)]

    distances[0] = 0
    queue.append(starting_point)

    while queue:
        u = queue.pop(0)
        for edge in edges:
            if int(edge[0]) == u:
                v = int(edge[1])
                if distances[v - 1] == -1:
                    queue.append(v)
                    distances[v - 1] = distances[u - 1] + 1

    #print(" ".join([str(x) for x in distances]))
    return distances

def depth_first_search(edges, start):
    """
    Finds all vertices in an undirected graph that are accessible from a vertex

    :param edges: a list of 2-member lists defining a graph's edges
    :type edges: list
    :param start: the name of vertex from which to begin exploring paths
    :type start: int
    :rtype: list
    :return: a list of the vertices accessible from start
    """
    stack, explored = [start], []

    while stack:
        u = stack.pop()
        if u not in explored:
            explored.append(u)
        for edge in edges:
            if u == edge[0] and edge[1] not in explored:
                # This is a new, unexplored neighbor
                stack.append(edge[1])
            elif u == edge[1] and edge[0] not in explored:
                # This condition makes DFS work for UNdirected graphs
                stack.append(edge[0])

    return explored

def connected_components(edge_file):
    """
    Returns the number of connected components in a graph

    :param edge_file: the path to a file in edge format that defines a graph
    :type edge_file: str
    :rtype: int
    :return: the number of connected components in the graph
    """
    explored = []
    components = 0
    total_vertices, _, edges = parse_edge_file(edge_file)

    for x in range(1, total_vertices + 1):
        if x not in explored:
            explored += depth_first_search(edges, x)
            components += 1

    return components
