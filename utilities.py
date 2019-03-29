import requests
from exceptions import ProteinNotFoundError


def factorial(val):
    """
    Returns val! for any positive integer val
    """
    if val >= 0:
        answer = 1
        for x in range(0, val):
            answer *= val - x

        return answer
    else:
        msg = "Factorials of negative numbers equate to division by zero."
        raise ZeroDivisionError(msg)

def list_permutations(A):
    """
    Implements Heap's algorithm to list all permutations of values in list A

    :param A: A list of values to be considered
    :type A: list
    :rtype: list
    :return: a list of all possible permutations of A
    """
    n = len(A)
    c = [0 for x in range(1, n+1)]  # Stores stack state across iterations
    permutations = []
    i = 0

    # The first permutation is the original ordering of values
    permutations.append(A)

    while i < n:
        if c[i] < i:
            if i % 2 == 0:
                A[0], A[i] = A[i], A[0]
            else:
                A[c[i]], A[i] = A[i], A[c[i]]

            # Add the latest permutation to the list
            B = A.copy()
            permutations.append(B)

            c[i] += 1   # Increment the state placeholder
            i = 0       # reset pointer to the beginning of A
        else:
            # Reset state and increment the pointer
            c[i] = 0
            i += 1

    return permutations


class UniprotClient(object):
    """
    A simple library for working with Uniprot data

    Note: data can be retrieved in XML format with https://www.uniprot.org/uniprot/protein_id.xml
    but significantly more data is returned. For just base info and sequence data, FASTA seems
    to be a better choice for client and server.
    """
    def __init__(self):
        self.fasta_url = "https://www.uniprot.org/uniprot/{0}.fasta"
        self.xml_url = "https://www.uniprot.org/uniprot/{0}."

    def get_protein(self, protein_id):
        """
        Returns a protein sequence given it's uniprot ID
        """
        protein_url = self.fasta_url.format(protein_id)
        response = requests.get(protein_url)
        
        if response.encoding == 'ISO-8859-1':
            # Valid response content with protein data
            fasta_text = response.text.strip()
            protein_dict = self.__parse_fasta(fasta_text)

            return protein_dict
        else:
            # 'Not Found' response is a UTF-8 encoded HTML page
            msg = "Protein ID {0} not found".format(protein_id)
            raise ProteinNotFoundError(msg)

    def __parse_fasta(self, fasta_text):
        """
        Returns an amino acid sequence given a uniprot FASTA response
        """
        protein = {}
        protein['sequence'] = ""

        if fasta_text:
            lines = fasta_text.split('\n')
            for line in lines:
                if line[0] == '>':
                    # Protein Meta info
                    meta = line.split(" ")
                    short_labels = meta[0].split("|")
                    protein['primary_accession_number'] = short_labels[1]
                    protein['name'] = short_labels[2]
                    protein['gene'] = self.__split_tag(meta[-3])
                    protein['taxonomic_id'] = self.__split_tag(meta[-4])
                    del meta[0]
                    del meta[-4:]

                    org_index = -1
                    for i, chunk in enumerate(meta):
                        if chunk[0:2] == 'OS':
                            org_index = i
                    
                    protein['organism'] = self.__split_tag(' '.join(meta[org_index:]))
                    protein['recommended_name'] = ' '.join(meta[0:org_index])

                else:
                    # Lines of sequence data
                    protein['sequence'] += line

        return protein

    def __split_tag(self, tagged_value):
        """
        Returns just the cleaned value for a tagged piece of meta info
        """
        return tagged_value.split("=")[-1]
