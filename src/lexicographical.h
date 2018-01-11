#ifndef H_LEXICO
#define H_LEXICO

#include <stdint.h>

// optimal amino acid ordering based on BLOSUM62
// ACMLJIVTSKRQZEDBNHYFWXP*G
// A  C  M  L  J  I  V  T  S  K  R  Q  Z  E  D  B  N  H  Y  F  W  X  P  *  G
// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
// A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y


// substitution matrix order
//  A    K    Q    O    B    L    N    Y    R    F    D    J    C    T    W    I    H    U    S    G    P    E    M    V    X
//  A,   R,   N,   D,   C,   Q,   E,   G,   H,   I,   L,   K,   M,   F,   P,   S,   T,   W,   Y,   V,   B,   J,   Z,   X,   *
//  0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24
//  A    K    Q    O    B    L    N    Y    R    F    D    J    C    T    W    I    H    U    S    G    P    E    M    V    X 

static int8_t aa_table2[128] = {
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//                                           *
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
//       A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
    23,  0,  4, 12, 10, 21,  9, 19, 16, 15, 11,  1,  5, 22,  6,  3,
//   P   Q   R   S   T   U   V   W   X   Y   Z
    20,  2,  8, 18, 13, 17, 23, 14, 24,  7, 23, 23, 23, 23, 23, 23,
//       A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
    23,  0,  4, 12, 10, 21,  9, 19, 16, 15, 11,  1,  5, 22,  6,  3,
//   P   Q   R   S   T   U   V   W   X   Y   Z
    20,  2,  8, 18, 13, 17, 23, 14, 24,  7, 23, 23, 23, 23, 23, 23
  };

static inline char upper_aa(char c) {
    switch(c) {
        case 'A': case 'B': case 'C': case 'D':
        case 'E': case 'F': case 'G': case 'H': 
        case 'I': case 'J': case 'K': case 'L': 
        case 'M': case 'N': case 'O': case 'P': 
        case 'Q': case 'R': case 'S': case 'T': 
        case 'U': case 'V': case 'W': case 'X': 
        case 'Y': case '*': case '\n': 
            return c;

        case 'a': case 'b': case 'c': case 'd':
        case 'e': case 'f': case 'g': case 'h':
        case 'i': case 'j': case 'k': case 'l':
        case 'm': case 'n': case 'o': case 'p': 
        case 'q': case 'r': case 's': case 't':
        case 'u': case 'v': case 'w': case 'x':
        case 'y': 
            return c - 32;        

        default:
            return 'X';
    }
}

static inline char internal2aa(char c) {
    switch(c) {
        case 'A': return 'A'; // 0
        case 'B': return 'C'; // 1
        case 'C': return 'M'; // 2
        case 'D': return 'L'; // 3
        case 'E': return 'J'; // 4
        case 'F': return 'I'; // 5
        case 'G': return 'V'; // 6
        case 'H': return 'T'; // 7
        case 'I': return 'S'; // 8
        case 'J': return 'K'; // 9
        case 'K': return 'R'; // 10
        case 'L': return 'Q'; // 11
        case 'M': return 'Z'; // 12
        case 'N': return 'E'; // 13
        case 'O': return 'D'; // 14
        case 'P': return 'B'; // 15
        case 'Q': return 'N'; // 16
        case 'R': return 'H'; // 17
        case 'S': return 'Y'; // 18
        case 'T': return 'F'; // 19
        case 'U': return 'W'; // 20
        case 'V': return 'X'; // 21
        case 'W': return 'P'; // 22
        case 'X': return '*'; // 23
        case 'Y': return 'G'; // 24
        case '\n': return '\n';
        default : return 'X';
    }
}

static inline char aa2internal(char c) {
    switch(c) {
        case 'A': return 'A';
        case 'C': return 'B';
        case 'M': return 'C';
        case 'L': return 'D';
        case 'J': return 'E';
        case 'I': return 'F';
        case 'V': return 'G';
        case 'T': return 'H';
        case 'S': return 'I';
        case 'K': return 'J';
        case 'R': return 'K';
        case 'Q': return 'L';
        case 'Z': return 'M';
        case 'E': return 'N';
        case 'D': return 'O';
        case 'B': return 'P';
        case 'N': return 'Q';
        case 'H': return 'R';
        case 'Y': return 'S';
        case 'F': return 'T';
        case 'W': return 'U';
        case 'X': return 'V';
        case 'P': return 'W';
        case '*': return 'X';
        case 'G': return 'Y';
        case '\n': return '\n';
        default : return 'V'; // i.e. 'X'
    }
}

#endif

