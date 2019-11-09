def rshift32(val, n): return (val % 0x100000000) >> n

def rshift64(val, n): return (val % 0x10000000000000000) >> n

def hash64shift(key):
    key = (~key) + (key << 21)  # key = (key << 21) - key - 1
    key = key ^ (rshift64(key, 24))
    key = (key + (key << 3)) + (key << 8)  # key * 265
    key = key ^ (rshift64(key, 14))
    key = (key + (key << 2)) + (key << 4)  # key * 21
    key = key ^ (rshift64(key, 28))
    key = (key + (key << 31)) % 0x10000000000000000
    return key

def hash32shift(key):
    key = ~key + (key << 15)  # key = (key << 15) - key - 1;
    key = key ^ (rshift32(key, 12))
    key = key + (key << 2)
    key = key ^ (rshift32(key, 4))
    key = ((key + (key << 3)) + (key << 11)) % 0x100000000
    key = key ^ (rshift32(key, 16))
    return key
