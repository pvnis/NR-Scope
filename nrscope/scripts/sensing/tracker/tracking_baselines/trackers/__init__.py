"""Tracker registry: name -> function(frames, **params) returning per frame [(id, x, y, vx, vy)].
None = our tracker (sim/track.py run_ekf_tracker, default settings, with its site knowledge).
Names: OURS_*, C* = classic trackers from their original papers, P* = trackers of the 5G/ISAC
papers we compare to (first author + year + what)."""
from . import c1_gnn_ekf, c2_gm_phd, c3_jipda, c4_mht, c5_bp_meyer, c6_ggiw_phd, p1_bauhofer2026, p2_bsense2026, p3_cellsense2026, p4_gunes2025

REGISTRY = {
    "OURS_full": None,
    "C1_GNN-EKF": c1_gnn_ekf.run,
    "C2_GM-PHD-VoMa2006": c2_gm_phd.run,
    "C3_JIPDA-Musicki2004": c3_jipda.run,
    "C4_MHT-CoxHingorani1996": c4_mht.run,
    "C5_BP-Meyer2017": c5_bp_meyer.run,
    "C6_GGIW-PHD-Granstrom2012": c6_ggiw_phd.run,
    "P1_Bauhofer2026_GM-PHD-range-doppler": p1_bauhofer2026.run,
    "P2_BSense2026_core": p2_bsense2026.run,
    "P3_CellSense2026_window-Hungarian": p3_cellsense2026.run,
    "P4_Gunes2025_plotfusion-GNN": p4_gunes2025.run,
}
SINGLE_SENSOR = {"P1_Bauhofer2026_GM-PHD-range-doppler"}   # no fused setup
