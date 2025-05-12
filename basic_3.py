import sys
import time
import psutil

def process_memory():
    process = psutil.Process()
    memory_info = process.memory_info()
    memory_consumed = int(memory_info.rss / 1024)
    return memory_consumed

# hardcoded penalty matrix
gap_penalty = 30


def mismatch_penalty(base1, base2):
    base = 'ACGT'
    i1 = base.find(base1)
    i2 = base.find(base2)
    penalty_matrix = [[0, 110, 48, 94],
                      [110, 0, 118, 48],
                      [48, 118, 0, 110],
                      [94, 48, 110, 0]]
    return penalty_matrix[i1][i2]


# read input from txt files
def parse_input():
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    f = open(input_file, "r")
    input_strings = []
    num1 = []
    num2 = []
    is_second = True
    for line in f.readlines():
        line = line.strip()
        if not line.isnumeric():
            input_strings.append(line)
            if len(input_strings) == 2:
                is_second = False
        else:
            if is_second:
                num1.append(int(line))
            else:
                num2.append(int(line))
    f.close()

    # return a list of two generated strings
    i = 0
    seq1 = input_strings[0]
    while i < len(num1):
        seq1 = seq1[0:num1[i] + 1] + seq1 + seq1[num1[i] + 1:]
        i = i + 1

    i = 0
    seq2 = input_strings[1]
    while i < len(num2):
        seq2 = seq2[0:num2[i] + 1] + seq2 + seq2[num2[i] + 1:]
        i = i + 1
    return seq1, seq2, output_file


# dynamic programming
def seq_alignment(seq1, seq2):
    s1_length = len(seq1)
    s2_length = len(seq2)

    # fill the opt matrix by dynamic programming
    DP = [[0 for col in range(s2_length + 1)] for row in range(s1_length + 1)]
    for i in range(s1_length + 1):
        DP[i][0] = gap_penalty * i

    for j in range(s2_length + 1):
        DP[0][j] = gap_penalty * j

    for i in range(1, s1_length + 1):
        for j in range(1, s2_length + 1):
            if seq1[i - 1] == seq2[j - 1]:
                min_value = min(DP[i - 1][j - 1],
                                DP[i - 1][j] + gap_penalty,
                                DP[i][j - 1] + gap_penalty)
                DP[i][j] = min_value
            else:
                min_value = min(DP[i - 1][j - 1] + mismatch_penalty(seq1[i - 1], seq2[j - 1]),
                                DP[i - 1][j] + gap_penalty,
                                DP[i][j - 1] + gap_penalty)
                DP[i][j] = min_value

    cost = DP[s1_length][s2_length]

    # find the alignment in the matrix
    i = s1_length
    j = s2_length
    path = [[0 for col in range(s2_length + 1)] for row in range(s1_length + 1)]
    path[i][j] = 1
    path[0][0] = 1
    while i != 0 or j != 0:
        if DP[i][j] == DP[i - 1][j] + gap_penalty:
            path[i - 1][j] = 1
            i = i - 1
        elif DP[i][j] == DP[i][j - 1] + gap_penalty:
            path[i][j - 1] = 1
            j = j - 1
        else:
            path[i - 1][j - 1] = 1
            i = i - 1
            j = j - 1

    s1_aligned = ''
    s2_aligned = ''
    i = 0
    j = 0

    while i != s1_length and j != s1_length:
        if i != s1_length and path[i + 1][j] == 1:
            s1_aligned = s1_aligned + seq1[i]
            s2_aligned = s2_aligned + '_'
            i = i + 1
        elif j != s2_length and path[i][j + 1] == 1:
            s1_aligned = s1_aligned + '_'
            s2_aligned = s2_aligned + seq2[j]
            j = j + 1
        else:
            s1_aligned = s1_aligned + seq1[i]
            s2_aligned = s2_aligned + seq2[j]
            i = i + 1
            j = j + 1

    return s1_aligned, s2_aligned, cost


## RUNNER CODE
start_time = time.time()
start_memory = process_memory()

s1, s2, output_file = parse_input()
s1_aligned, s2_aligned, cost = seq_alignment(s1, s2)

end_memory = process_memory()
end_time = time.time()

used_memory = end_memory - start_memory
runtime = (end_time - start_time) * 1000

# write the output
output_data = [str(cost) + '\n',
               s1_aligned + '\n',
               s2_aligned + '\n',
               str(runtime) + '\n',
               str(used_memory)]
f = open(output_file, 'w')
f.writelines(output_data)
f.close()