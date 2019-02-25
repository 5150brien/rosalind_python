import requests
from exceptions import ProteinNotFoundError


class UniprotClient(object):
    """
    A simple library for working with Uniprot data
    """
    def __init__(self):
        self.uniprot_url = "https://www.uniprot.org/uniprot/"

    def get_protein(self, protein_id):
        """
        Returns a protein sequence given it's uniprot ID
        """
        protein_url = self.uniprot_url + protein_id + '.fasta'
        response = requests.get(protein_url)
        
        if response.encoding == 'ISO-8859-1':
            # Valid response content with protein data
            fasta_text = response.text.strip()
            sequence = self.__parse_fasta(fasta_text)

            return sequence
        else:
            # 'Not Found' response is a UTF-8 encoded HTML page
            msg = "Protein ID {0} not found".format(protein_id)
            raise ProteinNotFoundError(msg)

    def __parse_fasta(self, fasta_text):
        """
        Returns an amino acid sequence given a uniprot FASTA response
        """
        sequence = ""

        if fasta_text:
            lines = fasta_text.split('\n')
            for line in lines:
                if line[0] == '>':
                    # Name, ID, Meta info... not needed yet
                    pass
                else:
                    # Lines of sequence data
                    sequence += line

        return sequence
                
