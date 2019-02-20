class InvalidSequenceError(Exception):
    def __init__(self, message=None):
        if message:
            self.message = message
        else:
            self.message = "Nucleotide sequence was not valid"

    def __str__(self):
        return self.message


class SequenceLengthMismatchError(Exception):
    def __init__(self, message=None):
        if message:
            self.message = message
        else:
            self.message = "Input sequences were not the same length"

    def __str__(self):
        return self.message


class PopulationError(Exception):
    def __init__(self, message=None):
        if message:
            self.message = message
        else:
            self.message = "Invalid population"

    def __str__(self):
        return self.message
