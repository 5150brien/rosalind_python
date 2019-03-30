import re
from decimal import Decimal
from utilities import UniprotClient, factorial, list_permutations
from exceptions import (InvalidSequenceError, SequenceLengthMismatchError,
                        PopulationError)

def count_dna_nucleotides(dna_string):
    """
    Counts the occurrence of each nucleotide bases in a DNA sequence.

    :param dna_string: A string of nucleotide bases
    :type dna_string: String
    :rtype: Dict
    :return: A dictionary of the 4 nucleotide bases and their occurrences
    """
    if valid_dna(dna_string):
        sequence = dna_string.upper()
        base_counts = {
            'Adenine': sequence.count("A"),
            'Cytosine': sequence.count("C"),
            'Guanine': sequence.count("G"),
            'Thymine': sequence.count("T"),
        }

        print("{a} {c} {g} {t}".format(
            a=base_counts['Adenine'],
            c=base_counts['Cytosine'],
            g=base_counts['Guanine'],
            t=base_counts['Thymine'],
        ))

        return base_counts
    else:
        raise InvalidSequenceError()

def transcribe_dna_to_rna(dna_string):
    """
    Converts a DNA sequence into the corresponding RNA sequence.

    :param dna_string: A string of nucleotide bases
    :type dna_string: String
    :rtype: String
    :return: A string of RNA nucleotide bases
    """
    if valid_dna(dna_string):
        return dna_string.upper().replace("T", "U")
    else:
        raise InvalidSequenceError()

def reverse_complement(dna_string):
    """
    Returns the reverse complement of a DNA sequence.

    :param dna_string: A string of nucleotide bases
    :type dna_string: String
    :rtype: String
    :return: The reverse compelement of dna_string
    """
    if valid_dna(dna_string):
        reverse_complement = ''
        for base in dna_string[::-1]:
            if base.upper() == "A":
                reverse_complement += "T"
            elif base.upper() == "T":
                reverse_complement += "A"
            elif base.upper() == "C":
                reverse_complement += "G"
            elif base.upper() == "G":
                reverse_complement += "C"
        return reverse_complement
    else:
        raise InvalidSequenceError()

def rabbit_recurrence(months=0, pairs_per_litter=1):
    """
    Calculates rabbit population as a Fibonacci recurrence relation.

    Given a number of months and a number of rabbit pairs per litter,
    calculates the total number of rabbit pairs that will exist at the end
    of the specified number of months. **Note that the sequence ALWAYS 
    starts with a litter of one pair, as specified by the problem parameters

    :param months: The total number of months to be considered
    :type months: Integer
    :param pairs_per_litter: the number of rabbit pairs born in each litter
    :type pairs_per_litter: Integer
    :rtype: Integer
    :return: The total number of rabbits after the specified number of months
    """
    baby_pairs = 1
    adult_pairs = 0
    for month in range(1, months):
        baby_pairs_last_month = baby_pairs
        baby_pairs = adult_pairs * pairs_per_litter
        adult_pairs += baby_pairs_last_month

    return adult_pairs + baby_pairs

def compute_gc_content(dna_string):
    """
    Calculates the percentage of a DNA sequence made of Cytosine or Guanine.

    GC content is useful for quickly differentiating between prokaryotic 
    and eukaryotic DNA sequences. 

    :param dna_string: A string of nucleotide bases
    :type dna_string: String
    :rtype: Decimal
    :return: GC content as a percentage
    """
    if valid_dna(dna_string):
        length = len(dna_string)
        sequence = dna_string.upper()
        gc_count = sequence.count("G") + sequence.count("C")
        gc_content = round(Decimal((gc_count)/length),6)*100
        return gc_content
    else:
        raise InvalidSequenceError()

def load_fasta(file_path):
    """
    A parser for FASTA datasets. 
    
    Reads a text file of DNA strings in FASTA format and returns them as a
    dictionary (ex {'Rosalind_4995': 'ATTGCTTGACCG'}).

    :param file_path: The path to the FASTA text file
    :type file_path: String
    :rtype: Dictionary
    :return: A dictionary of DNA strings with their IDs
    """
    lines = None
    sequence_id = None
    sequence = ""
    dna_dict = {}

    with open(file_path, "r") as fasta_file:
        for line in fasta_file:
            data = line.strip()
            if data[0] == '>':
                if sequence:
                    # A complete DNA sequence is ready to add
                    dna_dict[sequence_id] = sequence
                    sequence = ""
                sequence_id = data.replace('>', '')
            else:
                # This is sequence data and may be multiple lines
                sequence += data

        # FASTA data is delimited in a way that orphans the last sequence
        if sequence and sequence_id:
            dna_dict[sequence_id] = sequence

    return dna_dict

def get_highest_gc_content(dna_string_dict):
    """
    Finds the sequence with the highest GC content in a dictionary of DNA strings

    :param dna_string_dict: A dict of DNA strings (ex: {"Rosalind_0808": "ACATG"})
    :type dna_string_dict: Dictionary
    :rtype: String, Decimal
    :return: ID of dna string with highest GC content, GC content of that string
    """
    gc_leader = None
    highest_gc = 0.0

    for dna_id, dna_string in dna_string_dict.items():
        current_gc = compute_gc_content(dna_string)
        if current_gc > highest_gc:
            gc_leader = dna_id
            highest_gc = current_gc

    return gc_leader, highest_gc

def hamming_distance(dna_string_1, dna_string_2):
    """
    Calculates the Hamming distance between two strings of the same length.

    Hamming distance counts the number of point mutations separating two DNA
    strings. 

    :param dna_string_1: A string of nucleotide bases
    :type dna_string_1: String
    :param dna_string_2: A string of nucleotide bases
    :type dna_string_2: String
    :rtype: Integer
    :return: The Hamming distance between dna_string_1 and dna_string_2
    """
    if valid_dna(dna_string_1) and valid_dna(dna_string_2):
        if len(dna_string_1) == len(dna_string_2):
            point_mutations = 0
            for i, base in enumerate(dna_string_1):
                if base != dna_string_2[i]:
                    point_mutations += 1

            return point_mutations
        else:
            raise SequenceLengthMismatchError()
    else:
        raise InvalidSequenceError()

def valid_dna(dna_string):
    """
    Determines whether dna_string is a valid DNA sequence
    """
    if type(dna_string) is str:
        if all(base in 'ATCG' for base in dna_string.upper()):
            return True
    return False

def valid_rna(rna_string):
    """
    Determines whether rna_string is a valid RNA sequence
    """
    if type(rna_string) is str:
        if all(base in 'AUCG' for base in rna_string.upper()):
            return True
    return False

def valid_amino_acid(aa):
    """
    Returns True if aa is a valid amino acid, otherwise returns False
    """
    if isinstance(aa, str) and len(aa) == 1:
        valid_aas = [
            'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
        ]

        if aa.upper() in valid_aas:
            return True
    return False

def mendelian_inheritance(k, m, n):
    """
    Determines the probability of offspring expressing a dominant factor.

    In a population of k+m+n members, the probability that two
    randomly-selected mating organisms will produce offspring
    possessing a dominant allele is calculated. k represents the 
    number of homozygous dominant members of the population, m 
    represents the number that are heterozygous, and n represents the
    number that are homozygous recessive.

    Punnett Squares for mating in this population:
    k+k -> AA, AA, AA, AA : 4/4 (100%) dominant
    k+m -> AA, AA, Aa, Aa : 4/4 (100%) dominant
    k+n -> Aa, Aa, Aa, Aa : 4/4 (100%) dominant 
    n+m -> aA, aA, aa, aa : 2/4 (50%) dominant
    n+n -> aa, aa, aa, aa : 0/4 (0%) dominant
    m+m -> AA, Aa, Aa, aa : 1/4 (75%) dominant

    :param k: population count that are homozygous dominant
    :type k: Integer
    :param m: population members that are heterozygous
    :type m: Integer
    :param n: population members that are homozygous recessive
    :type n : Integer
    :rtype: Float
    :return: The probability of offspring with a dominant allele
    """
    population = k + m + n
    outcomes = []

    # K parent 1 + K parent 2 (100% chance of dominant offspring)
    outcomes.append((k/population) * ((k-1)/(population-1)) * 1.0)

    # Parent 1 is K, parent 2 is M (100% chance of dominant offspring)
    outcomes.append((k/population) * (m/(population-1)) * 1.0)

    # Parent 1 is K, parent 2 is N (100% chance of dominant offspring)
    outcomes.append((k/population) * (n/(population-1)) * 1.0)

    # Parent 1 is m, parent 2 is k (100% chance of dominant offspring)
    outcomes.append((m/population) * (k/(population-1)) * 1.0)

    # Parent 1 is m, parent 2 is m (75% chance of dominant offspring)
    outcomes.append((m/population) * ((m-1)/(population-1)) * 0.75)

    # parent 1 is m, parent 2 is n (50% chance of dominant offspring)
    outcomes.append((m/population) * (n/(population-1)) * 0.5)

    # Parent 1 is n, parent 2 is k (100% chance of dominant offspring)
    outcomes.append((n/population) * (k/(population-1)) * 1.0)

    # Parent 1 is n, parent 2 is m (50% chance of dominant offspring)
    outcomes.append((n/population) * (m/(population-1)) * 0.5)

    # Parent 1 is n, parent 2 is n (0% chance of dominant offspring)
    outcomes.append(0.0)

    return round(sum(outcomes), 6)

def translate_rna(rna_string):
    """
    Translates an RNA sequence into an amino acid sequence
    
    :param rna_string: an RNA sequence
    :type rna_string: String
    :rtype: String
    :return: A sequence of amino acids coded by rna_string
    """

    rna_codons = {
        'auu': 'i', 'auc': 'i', 'aua': 'i', 'aug': 'm', 'acu': 't', 
        'acc': 't', 'aca': 't', 'acg': 't', 'aau': 'n', 'aac': 'n',
        'aaa': 'k', 'aag': 'k', 'agu': 's', 'agc': 's', 'aga': 'r',
        'agg': 'r', 'cuu': 'l', 'cuc': 'l', 'cua': 'l', 'cug': 'l',
        'ccu': 'p', 'ccc': 'p', 'cca': 'p', 'ccg': 'p', 'cau': 'h',
        'cac': 'h', 'caa': 'q', 'cag': 'q', 'cgu': 'r', 'cgc': 'r',
        'cga': 'r', 'cgg': 'r', 'guu': 'v', 'guc': 'v', 'gua': 'v',
        'gug': 'v', 'gcu': 'a', 'gcc': 'a', 'gca': 'a', 'gcg': 'a',
        'gau': 'd', 'gac': 'd', 'gaa': 'e', 'gag': 'e', 'ggu': 'g',
        'ggc': 'g', 'gga': 'g', 'ggg': 'g', 'uuu': 'f', 'uuc': 'f',
        'uua': 'l', 'uug': 'l', 'ucu': 's', 'ucc': 's', 'uca': 's',
        'ucg': 's', 'uau': 'y', 'uac': 'y', 'uaa': 'stop', 
        'uag': 'stop', 'ugu': 'c', 'ugc': 'c', 'uga': 'stop', 
        'ugg': 'w',
    }

    if valid_rna(rna_string):
        amino_acid_sequence = []
        start = False
        stop = False

        for i in range(0, len(rna_string)+1, 3):
            codon = rna_string[i:i+3]
            if len(codon) == 3:
                try:
                    codon_value = rna_codons[codon.lower()]
                    if codon_value == 'm':
                        # AUG (Methionine) is the start codon
                        start = True
                    elif codon_value == 'stop':
                        stop = True

                    if start == True and stop == False:
                        amino_acid_sequence.append(codon_value)
                except KeyError:
                    # Not a valid codon sequence
                    pass

        return ''.join(amino_acid_sequence).upper()
    else:
        message = "rna_string was not valid RNA"
        raise InvalidSequenceError(message=message)

def find_motif(motif=None, sequence=None):
    """
    Finds a motif of nucleotides or amino acids in a larger sequence.

    :param motif: the nucleotides/amino acids pattern to be found in sequence
    :type motif: String
    :param sequence: the nucleotide/amino acid sequence to search for motif
    :type sequence: String
    :rtype: Integer, List
    :return: The number of times motif occurs in sequence and a list of the
             starting indexes (using 1-based numbering) of each occurrence
    """
    motif_length = len(motif)
    motif_locations = []

    for i in range(0, len(sequence)+1):
        if sequence[i:i+motif_length] == motif:
            motif_locations.append(i+1)

    occurrences = len(motif_locations)

    return occurrences, motif_locations

def consensus_and_profile(dna_dict):
    """
    Returns a consensus and profile for a dict of same-length DNA sequences.

    Profile here is a matrix of the counts for each nucleotide base in each
    position of a group of DNA sequences of the same length. Consensus is the
    DNA sequence that is the likely ancester of all the sequences in the dict,
    based on the frequency a particular base occurs in each position of that
    sequence.

    **Note problem parameters specify that profile rows are ordered A, C, G, T
    Ex: A: 0 0 0 0 0 0
        C: 0 0 0 0 0 0
        G: 0 0 0 0 0 0
        T: 0 0 0 0 0 0

    :param dna_dict: A dict of dna sequences
    :type dna_dict: dict
    :rtype: list, list
    :return: a profile matrix for the DNA sequences in dna_dict (returned as
             a list of lists), a list of ALL POSSIBLE consensus strings
    """
    if len(dna_dict):
        sequence_length = len(next(iter(dna_dict.values())))
    else:
        sequence_length = 0

    profile_matrix = [[0]*sequence_length for x in range(0, 4)]
    row_positions = { 0: 'A', 1: 'C', 2: 'G', 3:'T'}
    consensus_list = ['',]

    # Fill in the profile matrix
    for key, value in dna_dict.items():
        for position, base in enumerate(value):
            if base.upper() == 'A':
                profile_matrix[0][position] += 1
            elif base.upper() == 'C':
                profile_matrix[1][position] += 1
            elif base.upper() == 'G':
                profile_matrix[2][position] += 1
            elif base.upper() == 'T':
                profile_matrix[3][position] += 1
            else:
                m = "{0} is an invalid DNA sequence".format(key)
                raise InvalidSequenceError(message=m)

    # Write all possible consensus strings
    for i in range(0, sequence_length):
        base_counts = [profile_matrix[x][i] for x in range(0,4)]
        max_count = max(base_counts)
        base_indexes = [j for j, x in enumerate(base_counts) if x == max_count]

        for i, base_index in enumerate(base_indexes):
            consensus_list_copy = consensus_list
            if i == 0:
                # Only 1 consensus base to be added to existing strings
                for j, consensus_string in enumerate(consensus_list):
                    new_base = row_positions[base_index]
                    consensus_list[j] = consensus_string + new_base
            else:
                # Multiple consensus bases to be added (requires branching)
                new_branches = consensus_list_copy
                for j, consensus_string in new_branches:
                    new_base = row_positions[base_index]
                    new_branches[j] = consensus_string + new_base
                consensus_list += new_branches

    return consensus_list, profile_matrix

def mortal_fibonacci_rabbits(n, m):
    """
    Returns the rabbit population at month n for rabbits with m-month lifespan

    The ages of rabbit pairs are tracked in a list of length m where each
    element in the list represents the number of pairs aged list[i]. All adult
    pairs mate every month (producing a 1-pair litter) and adults are shifted
    over one position for each month they age until they fall off the list.

    :param n: the number of months for which the population is considered
    :type n: int
    :param m: the lifespan of a rabbit, in months
    :type m: int
    :rtype: int
    :return: the number of rabbits alive after n months
    """
    age_counts_by_month = [1] + ([0] * (m - 1))

    for month in range(n - 1):
        new_pairs = sum(age_counts_by_month[1:])
        surviving_pairs = age_counts_by_month[:-1]
        age_counts_by_month = [new_pairs] + surviving_pairs

    return sum(age_counts_by_month)

def overlap_graphs(dna_dict):
    """
    Returns an adjacency list to describe an overlap graph for dna sequences

    The overlap space, k, is defined to be exactly 3 for this problem.
    
    :param dna_dict: a dictionary of DNA sequences
    :type dna_dict: dict
    :rtype: list
    :return: an adjacency list; each item is a list of 2 edges with overlap
    """
    k = 3
    adjacency_list = []
    for key1, sequence1 in dna_dict.items():
        for key2, sequence2 in dna_dict.items():
            if key1 != key2:
                if sequence1[-k:] == sequence2[:k]:
                    adjacency_list.append([key1, key2])

    return adjacency_list

def calculate_expected_offspring(population):
    """
    Returns total offspring expected with dominant phenotype in a population

    The population will be expressed as a list of six positive integers,
    corresponding to the number of couples with the following genotypes:

    1. AA + AA
    2. AA + Aa
    3. AA + aa
    4. Aa + Aa
    5. Aa + aa
    6. aa + aa
    
    Each couple will produce exactly two offspring.

    :param population_list: a list of the total couples with each genotype
    :type population_list: list
    :rtype: float
    :return: the number of offspring expected to display the dominant phenotype
    """
    if len(population) == 6 and all(isinstance(x, int) for x in population):
        # Couples in groups 1-3 will always have dominant phenotype offspring
        a = 2 * (population[0] + population[1] + population[2])

        # Couples in group 4 will have ~75% dominant phenotype offspring
        b = 2 * population[3] * (3/4)

        # Couples in group 5 will have ~50% dominant phenotype offspring
        c = 2 * population[4] * (1/2)

        return a + b + c
    else:
        raise PopulationError("population must consist of six integers")

def find_shared_motif(dna_dict):
    """
    Returns a longest common substring from a dictionary of DNA sequences
    """
    sequences = sorted(dna_dict.values(), key=len)
    shortest = sequences.pop(0)
    motif_length = len(shortest)
    lcs = None

    while motif_length > 0 and not lcs:
        for start in range(len(shortest) - motif_length + 1):
            motif = shortest[start:start + motif_length]
            if all(find_motif(motif, s)[0] > 0 for s in sequences):
                lcs = motif

        motif_length -= 1
    return lcs

def independent_alleles(k, n):
    """
    Finds the probability that at least n children in generation k are Aa/Bb

    Uses binomial cumulative distribution function with a 'success rate' (the
    probability of Aa allele * Bb allele in any offspring) over 2^k trials

    :param k: the number of generations of reproduction to consider
    :type k: int
    :param n: the minimum number of successes observed (Aa/Bb offspring)
    :type n: int
    :rtype: float
    :return: the probability of at least n Aa/Bb offspring in generation k
    """
    trials = 2**k
    p = 0.0
    for x in range(n, trials + 1):
        p += (factorial(trials) / (factorial(x) * factorial(trials - x))) \
             * (0.25**x) * (0.75**(trials - x))
    return p

def find_protein_motif(protein_ids):
    """
    Searches uniprot for the N-glycosylation motif in a list of proteins

    **Nota bene: results dictionary gives 0-based motif locations, but the 
    print output provides 1-based formatting, as specified by the problem 
    parameters.

    :param protein_ids: a list of uniprot access id strings
    :type protein_ids: list
    :rtype: dict
    :return: a dictionary of motif locations for each protein
    """
    motif = re.compile('(?=(N[ACDEFGHIKLMNQRSTVWY][ST][ACDEFGHIKLMNQRSTVWY]))')
    results = {}

    u = UniprotClient()
    for _id in protein_ids:
        results[_id] = []
        data = u.get_protein(_id)
        results[_id] = [x.start() for x in re.finditer(motif, data['sequence'])]
    
    # Problem-specific print formatting
    for key, val in results.items():
        if val:
            print(key)
            print(' '.join(str(v + 1) for v in val))

    return results

def enumerate_gene_orders(n):
    """
    Lists all permutations length n

    :param n: The permutation length
    :type n: int
    :rtype: NoneType
    :return: No return; problem parameters specify print formatting
    """
    if isinstance(n, int):
        A = [x for x in range(1, n+1)]
        p = list_permutations(A)

        # Problem-specific print formatting
        print(len(p))
        for permutation in p:
            print(' '.join(str(x) for x in permutation))

        return p
    else:
        msg = "n must be an integer"
        raise TypeError(msg)

def calculate_protein_mass(protein_sequence):
    """
    Calculates the weight of a protein as the sum of amino acid masses

    Monoisotopic masses are given in Daltons/atomic mass units.

    :param protein_sequence: the AA sequence defining the protein
    :type protein_sequence: str
    :rtype: float
    :return: the total weight of protein_sequence in Daltons (Da)
    """
    monoisotopic_masses = {
        "A": 71.03711, "C": 103.00919, "D": 115.02694, "E": 129.04259,
        "F": 147.06841, "G": 57.02146, "H": 137.05891, "I": 113.08406,
        "K": 128.09496, "L": 113.08406, "M": 131.04049, "N": 114.04293,
        "P": 97.05276, "Q": 128.05858, "R": 156.10111, "S": 87.03203,
        "T": 101.04768, "V": 99.06841, "W": 186.07931, "Y": 163.06333,
    }

    if isinstance(protein_sequence, str):
        total_weight = 0.0
        for aa in protein_sequence:
            if valid_amino_acid(aa):
                aa_weight = monoisotopic_masses[aa]
                total_weight += aa_weight
            else:
                aa_msg = "protein_sequence contained an invalid amino acid"
                raise InvalidSequenceError(aa_msg)

        return round(total_weight, 3)
    else:
        type_msg = "protein_sequence must be a string"
        raise TypeError(type_msg)
