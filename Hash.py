
"""
Hash.rshift32()
Account for integer overflow, replaces the ">>>" operator found in
Thomas Wang's implementation since Python does not do overflow for integers
Input: The value to be hashed, the number to shift right by
Output: The shifted value
"""
def rshift32(val, n): return (val % 0x100000000) >> n

"""
Hash.rshift64()
Account for integer overflow, replaces the ">>>" operator found in 
Thomas Wang's implementation since Python does not do overflow for integers
Input: The value to be hashed, the number to shift right by
Output: The shifted value
"""
def rshift64(val, n): return (val % 0x10000000000000000) >> n

# Source: Thomas Wang, March 2007
# https://gist.github.com/badboy/6267743

"""
Hash.hash64shift()
Creates a 64-bit hash value, uniform hashing function
Input: The value to be hashed
"""
def hash64shift(key):
    key = (~key) + (key << 21)  # key = (key << 21) - key - 1
    key = key ^ (rshift64(key, 24))
    key = (key + (key << 3)) + (key << 8)  # key * 265
    key = key ^ (rshift64(key, 14))
    key = (key + (key << 2)) + (key << 4)  # key * 21
    key = key ^ (rshift64(key, 28))
    key = (key + (key << 31)) % 0x10000000000000000
    return key

"""
Hash.hash32shift()
Creates a 32-bit hash value, uniform hashing function
Input: The value to be hashed
"""

def hash32shift(key):
    key = ~key + (key << 15)  # key = (key << 15) - key - 1;
    key = key ^ (rshift32(key, 12))
    key = key + (key << 2)
    key = key ^ (rshift32(key, 4))
    key = ((key + (key << 3)) + (key << 11)) % 0x100000000
    key = key ^ (rshift32(key, 16))
    return key
