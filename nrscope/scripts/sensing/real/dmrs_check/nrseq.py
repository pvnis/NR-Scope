"""NR pseudo-random sequences, written from TS 38.211 (5.2.1, 7.4.1.1.1), independent of NR-Scope.
  gold(c_init, n)              c(n), Nc = 1600
  dmrs_cinit(slot, l, n_id, n_scid)   PDSCH DM-RS, 14 symbols per slot (lambda = 0)
  qpsk(c_init, n)              r(m) = (1 - 2 c(2m) + j (1 - 2 c(2m+1))) / sqrt(2), m < n"""
import numpy as np

NC = 1600


def gold(c_init, n):
    L = n + NC + 31; x1 = np.zeros(L, np.uint8); x2 = np.zeros(L, np.uint8); x1[0] = 1
    for i in range(31):
        x2[i] = (c_init >> i) & 1
    for i in range(L - 31):
        x1[i + 31] = x1[i + 3] ^ x1[i]
        x2[i + 31] = x2[i + 3] ^ x2[i + 2] ^ x2[i + 1] ^ x2[i]
    return x1[NC:NC + n] ^ x2[NC:NC + n]


def dmrs_cinit(slot, l, n_id, n_scid=0):
    return ((1 << 17) * (14 * slot + l + 1) * (2 * n_id + 1) + 2 * n_id + n_scid) % (1 << 31)


def qpsk(c_init, n):
    c = gold(c_init, 2 * n).astype(float)
    return ((1 - 2 * c[0::2]) + 1j * (1 - 2 * c[1::2])) / np.sqrt(2)


def qpsk_many(c_inits, n):
    """qpsk for many c_init at once -> (len(c_inits), n)"""
    c_inits = np.asarray(c_inits, np.int64); L = 2 * n + NC + 31; K = len(c_inits)
    x1 = np.zeros(L, np.uint8); x1[0] = 1
    for i in range(L - 31):
        x1[i + 31] = x1[i + 3] ^ x1[i]
    x2 = np.zeros((K, L), np.uint8)
    for i in range(31):
        x2[:, i] = (c_inits >> i) & 1
    for i in range(L - 31):
        x2[:, i + 31] = x2[:, i + 3] ^ x2[:, i + 2] ^ x2[:, i + 1] ^ x2[:, i]
    c = (x1[None, NC:NC + 2 * n] ^ x2[:, NC:NC + 2 * n]).astype(float)
    return ((1 - 2 * c[:, 0::2]) + 1j * (1 - 2 * c[:, 1::2])) / np.sqrt(2)
