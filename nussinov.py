## Nussinov RNA folding algorithm + recursive backtrack. Implemented by Carlos G. Oliver ##

import sys

import numpy as np

min_loop_length = 4

def pair_check(tup):
    if tup in [('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')]:
        return True
    return False

def OPT(i,j, sequence):
    """ 
    returns the score of the optimal pairing between indices i and j
    """
    #base case: no pairs allowed when i and j are less than 4 bases apart
    if i >= j-min_loop_length:
        return 0
    else:
        #i and j can either be paired or not be paired, if not paired then the optimal score is OPT(i,j-1)
        unpaired = OPT(i, j-1, sequence)

        #check if j can be involved in a pairing with a position t
        pairing = [1 + OPT(i, t-1, sequence) + OPT(t+1, j-1, sequence) for t in range(i, j-min_loop_length)\
                   if pair_check((sequence[t], sequence[j]))]
        if not pairing:
            pairing = [0]
        paired = max(pairing)


        return max(unpaired, paired)


def traceback(i, j, structure, DP, sequence):
    #in this case we've gone through the whole sequence. Nothing to do.
    if j <= i:
        return
    #if j is unpaired, there will be no change in score when we take it out, so we just recurse to the next index
    elif DP[i][j] == DP[i][j-1]:
        traceback(i, j-1, structure, DP, sequence)
    #consider cases where j forms a pair.
    else:
        #try pairing j with a matching index k to its left.
        for k in [b for b in range(i, j-min_loop_length) if pair_check((sequence[b], sequence[j]))]:
            #if the score at i,j is the result of adding 1 from pairing (j,k) and whatever score
            #comes from the substructure to its left (i, k-1) and to its right (k+1, j-1)
            if k-1 < 0:
                if DP[i][j] == DP[k+1][j-1] + 1:
                    structure.append((k,j))
                    traceback(k+1, j-1, structure, DP, sequence)
                    break
            elif DP[i][j] == DP[i][k-1] + DP[k+1][j-1] + 1:
                #add the pair (j,k) to our list of pairs
                structure.append((k,j))
                #move the recursion to the two substructures formed by this pairing
                traceback(i, k-1, structure, DP, sequence)
                traceback(k+1, j-1, structure, DP, sequence)
                break

def write_structure(sequence, structure):
    dot_bracket = ["." for _ in range(len(sequence))]
    for s in structure:
        dot_bracket[min(s)] = "("
        dot_bracket[max(s)] = ")"
    return "".join(dot_bracket)


#initialize matrix with zeros where can't have pairings
def initialize(N):
    #NxN matrix that stores the scores of the optimal pairings.
    DP = np.empty((N,N))
    DP[:] = np.NAN
    for k in range(0, min_loop_length):
        for i in range(N-k):
            j = i + k
            DP[i][j] = 0
    return DP

def nussinov(sequence):
    N = len(sequence)
    DP = initialize(N)
    structure = []

    #fill the DP matrix diagonally
    for k in range(min_loop_length, N):
        for i in range(N-k):
            j = i + k
            DP[i][j] = OPT(i,j, sequence)

    #copy values to lower triangle to avoid null references
    for i in range(N):
        for j in range(0, i):
            DP[i][j] = DP[j][i]


    traceback(0,N-1, structure, DP, sequence)
    return (write_structure(sequence, structure), structure)
if __name__ == "__main__":
	print(nussinov(sys.argv[1]))
