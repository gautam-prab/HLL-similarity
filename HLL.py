import Cardinality
import Hash
# HyperLogLog implementation for our project
# Using cardinality estimation algorithms from:
    # Ertl, O. (2017). New cardinality estimation algorithms for HyperLogLog sketches. ArXiv.

class HLL:
    """HyperLogLog implementation"""

    def __init__(self, p, registers=None):
        """Initialize HLL with 2^p registers"""
        self.p = p
        self.q = 64 - p # assuming 32-bit hash value
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

        read_hash = Hash.hash64shift()
        read_hash = "{0:b}".format(read_hash)  # Converts into binary
        a = read_hash[0:self.p]  # First p bits of (p + q)-bit hash value of read
        b = read_hash[self.p:]  # Following q bits of (p + q)-bit hash value of read

        k = b.find('1')
        if k == -1:
            k = self.q

        i = int(a) + 1

        if k > self.registers[i]:
            self.registers[i] = k
        return

    def cardinality(self):
        C = Cardinality.getMultiplicity(self)
        Cardinality.estimateCardinality(C)
        return

    def getRegisters(self):
        return self.registers
