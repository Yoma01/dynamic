import sys
import time
import psutil
import tracemalloc


def generate_strings(input_path):
    with open(input_path, 'r') as f:
        lines = f.read().splitlines()

    idx = 0
    s1 = lines[idx]
    idx += 1
    #counting how many insertion indices for s1
    j = 0
    while idx < len(lines) and lines[idx].isdigit():
        j += 1
        idx += 1
    s1_indices = list(map(int, lines[1:1+j]))
    s2 = lines[idx]
    s2_indices = list(map(int, lines[idx+1:]))

    #creating complete s1...expanding s1 by inserting it at each index
    complete_s1 = s1
    for pos in s1_indices:
        complete_s1 = complete_s1[:pos+1] + complete_s1 + complete_s1[pos+1:]

    #creating s2
    complete_s2 = s2
    for pos in s2_indices:
        complete_s2 = complete_s2[:pos+1] + complete_s2 + complete_s2[pos+1:]

    # validation check
    if len(complete_s1) != (2 ** len(s1_indices)) * len(s1):
        raise ValueError(f"Length of string s1 ({len(complete_s1)}) doesnt match the length expected: ({ (2 ** len(s1_indices)) * len(s1) })")
    if len(complete_s2) != (2 ** len(s2_indices)) * len(s2):
        raise ValueError(f"Length of string s2 ({len(complete_s2)}) doesnt match the length expected: ({ (2 ** len(s2_indices)) * len(s2) })")
    

    if len(complete_s1) > 20000:
        raise ValueError(f"Expanded s1 is bigger than 20000")
    if len(complete_s2) > 20000:
        raise ValueError(f"Expanded s2 is bigger than 20000")
        

    return complete_s1, complete_s2

#alignment constants
GAP = 30
MISMATCHCOST = {
    ('A', 'A'): 0, ('A', 'C'): 110, ('A', 'G'): 48, ('A', 'T'): 94,
    ('C', 'A'): 110, ('C', 'C'): 0, ('C', 'G'): 118, ('C', 'T'): 48,
    ('G', 'A'): 48, ('G', 'C'): 118, ('G', 'G'): 0, ('G', 'T'): 110,
    ('T', 'A'): 94, ('T', 'C'): 48, ('T', 'G'): 110, ('T', 'T'): 0,
}

#computing space efficient forward alignment cost
def efficientAlignment(X, Y):
    m, n = len(X), len(Y)
    prev = [i * GAP for i in range(n+1)]

    for i in range(1, m+1):
        curr = [i * GAP] + [0] * n
        for j in range(1, n+1):
            match = prev[j-1] + MISMATCHCOST[(X[i-1], Y[j-1])]
            delete = prev[j] + GAP
            insert = curr[j-1] + GAP
            curr[j] = min(match, delete, insert)
        prev = curr
    return prev

#computing cost from reversed substrings
def efficientAlignment_rev(X, Y):
    X = X[::-1]
    Y = Y[::-1]
    return efficientAlignment(X, Y)

#the full complete DP alignment with traceback
def basic_alignment(X, Y):
    m, n = len(X), len(Y)
    dp = []
    for i in range(m+1):
        row = [0] * (n+1)
        dp.append(row)

    #init DP table with gap penalties
    for i in range(m+1):
        dp[i][0] = i * GAP
    for j in range(n+1):
        dp[0][j] = j * GAP

    #filling DP table with min alignment costs
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = dp[i-1][j-1] + MISMATCHCOST[(X[i-1], Y[j-1])]
            delete = dp[i-1][j] + GAP
            insert = dp[i][j-1] + GAP
            dp[i][j] = min(match, delete, insert)

    #traceback for reconstruction
    alignedX = []
    alignedY = []
    i, j = m, n
    while i > 0 and j > 0:
        if dp[i][j] == dp[i-1][j-1] + MISMATCHCOST[(X[i-1], Y[j-1])]:
            alignedX.append(X[i-1])
            alignedY.append(Y[j-1])
            i -= 1
            j -= 1
        elif dp[i][j] == dp[i-1][j] + GAP:
            alignedX.append(X[i-1])
            alignedY.append('_')
            i -= 1
        else:
            alignedX.append('_')
            alignedY.append(Y[j-1])
            j -= 1
    #handling remaining characters
    while i > 0:
        alignedX.append(X[i-1])
        alignedY.append('_')
        i -= 1
    while j > 0:
        alignedX.append('_')
        alignedY.append(Y[j-1])
        j -= 1
    return ''.join(reversed(alignedX)), ''.join(reversed(alignedY))


def hirschberg(X, Y):
    if len(X) == 0:
        return '_' * len(Y), Y
    if len(Y) == 0:
        return X, '_' * len(X)
    if len(X) == 1 or len(Y) == 1:
        return basic_alignment(X, Y)

    mid = len(X) // 2
    
    #forward and reverse scores
    score_l = efficientAlignment(X[:mid], Y)
    score_r = efficientAlignment_rev(X[mid:], Y)
    
    #finding the split point in Y using forward and reverse DP
    total = []
    rev_score_r = list(reversed(score_r))
    for i in range(len(score_l)):
        total.append(score_l[i] + rev_score_r[i])
    split = total.index(min(total))

    XL, YL = hirschberg(X[:mid], Y[:split])
    XR, YR = hirschberg(X[mid:], Y[split:])
    return XL + XR, YL + YR

def main():
    inputfile = sys.argv[1]
    outputfile = sys.argv[2]

    X, Y = generate_strings(inputfile)

    tracemalloc.start()
    starttime = time.time()

    alignedX, alignedY = hirschberg(X, Y)

    endtime = time.time()
    _, peakmemory = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    cost = 0
    for a, b in zip(alignedX, alignedY):
        if a == '_' or b == '_':
            cost += GAP
        else:
            cost += MISMATCHCOST[(a, b)]

    timeMS = (endtime - starttime) * 1000
    memory = peakmemory / 1024

    with open(outputfile, 'w') as f:
        f.write(str(cost) + '\n')
        f.write(alignedX + '\n')
        f.write(alignedY + '\n')
        f.write(f"{timeMS}" + '\n')
        f.write(f"{memory}" + '\n')



if __name__ == "__main__":
    main()
