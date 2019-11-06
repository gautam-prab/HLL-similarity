import sys
import math
import random

def variable_assertions(argv):
    assert(len(argv) == 4)
    jaccard_val = float(argv[0])
    num_reads_a = int(argv[1])
    num_reads_b = int(argv[2])
    read_length = int(argv[3])
    assert(jaccard_val >= 0 and jaccard_val <= 1)
    assert(num_reads_a > 0)
    assert(num_reads_b > 0)
    assert(read_length > 0)
    return jaccard_val, num_reads_a, num_reads_b, read_length

def generate_random_string(read_length):
    ACGT = "ACGT"
    string = ""

    for i in range(read_length):
        string += random.choice(ACGT)
        
    return string

def generate_reads(jaccard_val, num_reads_a, num_reads_b, read_length):
    # Instantiate key variables
    total_reads = num_reads_a + num_reads_b
    num_overlapped = math.ceil(total_reads * jaccard_val)
    a = []
    b = []

    # If the num_overlapped amount is greater than either read
    # amount, then tell user and print the maximum jaccard value
    # in this case
    if num_overlapped > num_reads_a or num_overlapped > num_reads_b:
        print("Given the number of reads for both, the Jaccard value must be less than:")
        print(str(min(num_reads_a, num_reads_b) / total_reads))
        return None, None

    # Handle overlaps
    for i in range(num_overlapped):
        string = generate_random_string(read_length)
        a.append(string)
        b.append(string)
    
    for i in range(num_reads_a - num_overlapped):
        string = generate_random_string(read_length)
        a.append(string)

    for i in range(num_reads_b - num_overlapped):
        string = generate_random_string(read_length)
        b.append(string)
    
    random.shuffle(a)
    random.shuffle(b)

    return a, b

def main(argv):
    jaccard_val, num_reads_a, num_reads_b, read_length = variable_assertions(argv)
    a, b = generate_reads(jaccard_val, num_reads_a, num_reads_b, read_length)
    if a is None or b is None:
        return 

    experiment = input("Experiment Name: ")
    filename_a = experiment + "_" + str(num_reads_a) + "_" + str(jaccard_val) + "_a.txt"
    filename_b = experiment + "_" + str(num_reads_b) + "_" + str(jaccard_val) + "_b.txt"

    file_a = open(filename_a, "w")
    for i, line in enumerate(a):
        file_a.write("> " + str(i) + "\n")
        file_a.write(line)
        file_a.write("\n")
    file_a.close()

    file_b = open(filename_b, "w")
    for i, line in enumerate(b):
        file_b.write("> " + str(i) + "\n")
        file_b.write(line)
        file_b.write("\n")
    file_b.close()


if __name__ == "__main__":
    main(sys.argv[1:])