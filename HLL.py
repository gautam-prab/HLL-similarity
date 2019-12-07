import Cardinality
import Hash

from Random_Generators import Rangen_jaccard
# HyperLogLog implementation for our project
# Using cardinality estimation algorithms from:
# Ertl, O. (2017). New cardinality estimation algorithms for HyperLogLog sketches. ArXiv.

class HLL:
    """HyperLogLog implementation"""

    def __init__(self, p, registers=None):
        # Initialize HLL with 2^p registers
        self.p = p
        self.q = 64 - p  # assuming 64-bit hash value
        self.m = 2 ** p

        if registers != None:
            assert len(registers) == self.m
            self.registers = registers
        else:
            self.registers = [0 for i in range(self.m)]

    """
    HLL.insert():
    Inserts reads into the HLL, utilizes a uniform hashing function and methods
    Comes from Algorithm 1 in Ertl paper
    
    Input: The HLL, the read
    """
    def insert(self, read):
        # Insert read into the HLL, convert read to base 4 -> base 10 -> base 2
        read.upper()  # Capitalization ensures proper conversion
        b_four_read = ""
        # Base 4 conversion
        for b in read:
            if b == 'A':
                b_four_read += "0"
            elif b == 'C':
                b_four_read += "1"
            elif b == 'G':
                b_four_read += "2"
            elif b == 'T':
                b_four_read += "3"

        b_ten_read = 0
        read_len = len(b_four_read)

        # Base 10 conversion
        for c in b_four_read:
            b_ten_read = b_ten_read + (int(c) * (4 ** (read_len - 1)))
            read_len = read_len - 1

        # Call to hash function, returns a 64-bit hash value
        read_hash = Hash.hash64shift(b_ten_read)

        read_hash = "{0:b}".format(read_hash)  # Converts into binary

        len_read_hash = len(read_hash)

        # Ensures proper length of the hashed value in case the value is too small for extend to 64 bits
        if len_read_hash < 64:
            zeroes = "".join(["0"] * (64 - len_read_hash))
            read_hash = zeroes + read_hash

        a = read_hash[0:self.p]  # First p bits of (p + q)-bit hash value of read
        b = read_hash[self.p:]  # Following q bits of (p + q)-bit hash value of read

        k = b.find('1') + 1  # be sure to add 1 to match match from Ertl paper
        if k == -1:
            k = self.q + 1

        i = int(a, 2)  # Indexing into HLL

        # Insert into HLL
        if k > self.registers[i]:
            self.registers[i] = k
        return

    """
    HLL.cardinality():
    Calls the getMultiplicity() function to get the count vector from HLL registers
    Calls the estimateCardinality() function to calculate the maximum likelihood estimator for the cardinality of a HLL
    
    Input: The HLL
    Output: The cardinality estimate
    """
    def cardinality(self):
        C = Cardinality.getMultiplicity(self)
        return Cardinality.estimateCardinality(C)

    """
    HLL.simple_cardinality():
    Calls the simpleCardinality() function to calculate cardinality of HLL by bias-corrected harmonic mean
    
    Input: The HLL
    Output: The cardinality estimate
    """
    def simple_cardinality(self):
        return Cardinality.simpleCardinality(self.getRegisters())

    """
    HLL.getRegisters():
    Returns the registers of the HLL, a list of values
    
    Input: The HLL
    Output: The registers
    """
    def getRegisters(self):
        return self.registers

# Test basic HLL functionality
def main():
    h = HLL(8)
    for i in range(10000):
        h.insert(Rangen_jaccard.generate_random_string(40)) # probabilistically these are all distinct
    print('MLE: {}, Raw: {}'.format(h.cardinality(), h.simple_cardinality()))

if __name__ == "__main__":
    main()
