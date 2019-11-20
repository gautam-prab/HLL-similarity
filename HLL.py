import Cardinality
import Hash

import Random_gen
# HyperLogLog implementation for our project
# Using cardinality estimation algorithms from:
# Ertl, O. (2017). New cardinality estimation algorithms for HyperLogLog sketches. ArXiv.

class HLL:
    """HyperLogLog implementation"""

    def __init__(self, p, registers=None):
        """Initialize HLL with 2^p registers"""
        self.p = p
        self.q = 64 - p # assuming 64-bit hash value
        self.m = 2 ** p

        if registers != None:
            assert len(registers) == self.m
            self.registers = registers
        else:
            self.registers = [0 for i in range(self.m)]

    def insert(self, read):
        """Insert read into the HLL"""
        # read is the sequence, convert read to base 4 -> base 10 -> base 2
        read.upper()
        b_four_read = ""
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

        for c in b_four_read:
            b_ten_read = b_ten_read + (int(c) * (4 ** (read_len - 1)))
            read_len = read_len - 1

        read_hash = Hash.hash64shift(b_ten_read)

        read_hash = "{0:b}".format(read_hash)  # Converts into binary

        len_read_hash = len(read_hash)

        if len_read_hash < 64:
            zeroes = "".join(["0"] * (64 - len_read_hash))
            read_hash = zeroes + read_hash

        a = read_hash[0:self.p]  # First p bits of (p + q)-bit hash value of read
        b = read_hash[self.p:]  # Following q bits of (p + q)-bit hash value of read

        k = b.find('1')+1 # be sure to add 1 to match match from Ertl paper
        if k == -1:
            k = self.q + 1

        i = int(a, 2)

        if k > self.registers[i]:
            self.registers[i] = k
        return

    def cardinality(self):
        C = Cardinality.getMultiplicity(self)
        print(C)
        return Cardinality.estimateCardinality(C, self.m), Cardinality.simpleCardinality(self.getRegisters())
        #return Cardinality.getMultiplicity(self)

    def getRegisters(self):
        return self.registers

# Test basic HLL functionality
def main():
    h = HLL(8)
    for i in range(1000):
        h.insert(Random_gen.generate_random_string(40)) # probabilistically these are all distinct
    print(h.cardinality())
    #print(h.getRegisters())

if __name__ == "__main__":
    main()
