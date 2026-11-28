"""Resource grid of one slot from the raw IQ, aligned to NR-Scope's slot numbering by matching its raw TRS REs.
  Grid(iq_path, chunk_start_s, slot0_sample, a_offset): .slot(a, ants) -> (n_ant, 14, 3276) complex (carrier REs)"""
import numpy as np
import ofdm

NSC, D = 3276, (4096 - 3276) // 2   # carrier centred on the FFT (found by matching)


class Grid:
    def __init__(self, paths, chunk_start_s, slot0, a_off):
        self.paths, self.s0, self.slot0, self.a_off = paths, int(round(chunk_start_s * ofdm.FS)), slot0, a_off

    def slot(self, a, ants=(0,)):
        n = a - self.a_off; st = self.slot0 + n * ofdm.SLOT; out = np.zeros((len(ants), 14, NSC), complex)
        for i, c in enumerate(ants):
            x = ofdm.read(self.paths[c], self.s0 + st, ofdm.SLOT + 512).astype(np.complex128)
            for l in range(14):
                out[i, l] = ofdm.symbol(x, ofdm.START[l], ofdm.CP[l])[D:D + NSC]
        return out
