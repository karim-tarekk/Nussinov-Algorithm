# Developer:
# Karim Tarek Emam
# ---------------------------------------------------------

def CreateMatrix(seq):
    # This Function will create a 2D array of size seqlen * seqlen
    seqlen = len(seq)
    matrix = [[0 for x in range(seqlen)] for y in range(seqlen)] 
    return matrix

def couple(pair):
    # This Function will check if the two AA are couplable or not
    pairs = {"A": "U", "U": "A", "G": "C", "C": "G"} 
    # check if pair is couplable
    if pair in pairs.items():
        return True
    return False


def fillMatrix(mat, seq):
    # This Funtion will fill the matrix diagonally to find the optimal solution
    seqlen = len(seq)
    for g in range(1,seqlen):
        i = 0
        j = g
        while j < seqlen :
            daig = mat[i + 1][j - 1] + couple((seq[i], seq[j]))
            down = mat[i+1][j]
            left = mat[i][j-1] 
            value = max([mat[i][k] + mat[k + 1][j] for k in range(i, j)])
            mat[i][j] = max(left, down, value, daig)
            i += 1
            j += 1
    return mat


def trace(mat, seq, fold, i, L):
    # This Function Traceback recursively through complete matrix to find secondary structure according to the four rules
    j = L
    if i < j:
        if mat[i][j] == mat[i][j - 1]: 
            trace(mat, seq, fold, i, j - 1)
        elif mat[i][j] == mat[i + 1][j - 1] + couple((seq[i], seq[j])): 
            fold.append((i, j))
            trace(mat, seq, fold, i + 1, j - 1)
        elif mat[i][j] == mat[i + 1][j]: 
            trace(mat, seq, fold, i + 1, j)
        else:
            for k in range(i + 1, j - 1):
                if mat[i][j] ==  mat[i][k] + mat[k + 1][j]:
                    trace(mat, seq, fold, i, k)
                    trace(mat, seq, fold, k + 1, j)
                    break
    return fold

def dot_write(seq, fold):
    # This Function will form our dotbracket string that converts the tuble in each index in Fold list and add '(' to the smallest number
    # which represents the first AA in the dot list which is equal to seq len
    dot = ["." for i in range(len(seq))]
    for s in fold:
        dot[min(s)] = "("
        dot[max(s)] = ")"
    return "".join(dot)

if __name__ == '__main__':
    seq = "GGGAAAUCC"
    seqlen = len(seq)
    matrix = CreateMatrix(seq)
    matrix = fillMatrix(matrix,seq)
    folds = []
    dotBracket = trace(matrix,seq,folds,0,seqlen-1)
    res = dot_write(seq, folds)
    print(res)
