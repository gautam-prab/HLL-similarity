import sys
import math
import random

def variable_assertions(argv):
    assert(len(argv) == 4)
    forbes_val = float(argv[0])
    num_reads_a = int(argv[1])
    num_reads_b = int(argv[2])
    read_length = int(argv[3])
    assert(forbes_val >= 0)
    assert(num_reads_a > 0)
    assert(num_reads_b > 0)
    assert(read_length > 0)
    return forbes_val, num_reads_a, num_reads_b, read_length

def generate_random_string(read_length):
    ACGT = "ACGT"
    string = ""

    for i in range(read_length):
        string += random.choice(ACGT)
        
    return string

def generate_reads(forbes_val, num_reads_a, num_reads_b, read_length):
    # Instantiate key variables
    num_overlapped = math.ceil(forbes_val * num_reads_a * num_reads_b)
    a = []
    b = []
    # If the num_overlapped amount is greater than either read
    # amount, then tell user and print the maximum forbes value
    # in this case
    if num_overlapped > num_reads_a or num_overlapped > num_reads_b:
        print("Given the number of reads for both, the Forbes value must be less than:")
        print(str(min(num_overlapped, num_reads_a, num_reads_b) / (num_reads_a * num_reads_b)))
        return None, None, None

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

    # Calculate new jaccard value - if changed at all
    forbes = num_overlapped / (num_reads_a * num_reads_b)

    random.shuffle(a)
    random.shuffle(b)

    return a, b, forbes

def main(argv):
    forbes_val, num_reads_a, num_reads_b, read_length = variable_assertions(argv)
    a, b, forbes = generate_reads(forbes_val, num_reads_a, num_reads_b, read_length)
    if a is None or b is None:
        return 

    experiment = input("Experiment Name: ")
    filename_a = experiment + "_" + str(num_reads_a) + "_" + str(forbes_val) + "_a.txt"
    filename_b = experiment + "_" + str(num_reads_b) + "_" + str(forbes_val) + "_b.txt"

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

    print("The Forbes value is: " + str(forbes))

if __name__ == "__main__":
    main(sys.argv[1:])