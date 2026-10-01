"""
Sensing observation-window bounds: visualising the feasible t_span.

The range-Doppler map is built over an observation window of duration t_span.
That duration is squeezed from two sides (see nr_ue_map.h, "OBSERVATION WINDOW"):

  LOWER bound  (want a LONG window)
    4. velocity resolution:  dv = lambda / (2 t)  ->  t >= lambda / (2 dv)

  UPPER bounds (want a SHORT window), smallest wins:
    1. grid budget:      4 f_max t <= N_freq - 1,  f_max = 2 v / lambda + 1/T_tdd
    2. range migration:  target stays in one range cell, t <= cell_m / v
    3. Doppler/accel:    target stays in one Doppler cell while accelerating

The window is FEASIBLE only where the lower bound is below the smallest upper
bound. This script plots that band and the crossings.

Ground truth is the C pipeline (nrscope/src/libs/sensing/nr_ue_map.c). Two
formulas there differ from the first draft of this file, and both choices move
the feasibility line, so both variants are plotted (solid = C code, dashed =
alternative):

  * migration:  C uses cell_m / v ; alt uses cell_m / (2 v)
  * accel:      C code uses cbrt(lambda/2a) = ^(1/3), but its own derivation
                comment (nr_ue_map.c:176) and nothing physical gives ^(1/2).
                Toggle ACCEL_EXP below.
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Radio / cell constants (config_x410_benetel.yaml + sensing C defines)
# ---------------------------------------------------------------------------
CENTER_FREQUENCY = 3.45e9        # Hz, rx_center_freq
SCS              = 30_000        # Hz, scs_index 1
SPEED_OF_LIGHT   = 299_792_458.0 # m/s
TDD_PERIOD_SLOTS = 10            # nr_sensing_tdd_period_slots (7D 1S 2U Benetel cell)
SLOT_DUR_S       = 1e-3 / (SCS / 15_000.0)          # 0.5 ms at 30 kHz
TDD_PERIOD_S     = TDD_PERIOD_SLOTS * SLOT_DUR_S     # 5 ms -> 1/T = 200 Hz

# Doppler axis point budget, NR_SENSING_MAP_MAX_BINS_FREQ
DOPPLER_POINTS   = 1024

# Widest expected grant: full 273-PRB carrier, DM-RS type 1 (half the carriers).
N_RB      = 273
N_PILOTS  = N_RB * 6             # 6 type-1 pilots per PRB per symbol = 1638
K_STEP    = 2                    # TEST_K_STEP, DM-RS type 1 rank 1

# Target scene assumption, NR_SENSING_TARGET_MAX_ACCEL_MS2
MAX_ACCELERATION = 1.0           # m/s^2

# Resolution targets to compare. C default NR_SENSING_TARGET_DV_MS = 0.30.
DV_TARGETS = [0.1, 0.2, 0.3]     # m/s

# ^(1/2) is the physically-derived exponent (nr_ue_map.c:176); the C code's cbrt
# (^1/3) is treated here as a bug to fix.
ACCEL_EXP = 1.0 / 2.0


def lambda_():
    return SPEED_OF_LIGHT / CENTER_FREQUENCY


# ---------------------------------------------------------------------------
# The four criteria, each a function of velocity coverage v (m/s).
# All return seconds. These mirror nr_ue_sensing_span_* in nr_ue_map.c.
# ---------------------------------------------------------------------------
def f_max(v):
    """Half-width of the Doppler axis: fastest target + one replica period."""
    return 2.0 * v / lambda_() + 1.0 / TDD_PERIOD_S


def span_grid(v, doppler_points=DOPPLER_POINTS):
    """Criterion 1 (upper): 4 f_max t <= N_freq - 1."""
    return (doppler_points - 1) / (4.0 * f_max(v))


def span_migration(v, n_pilots=N_PILOTS, k_step=K_STEP):
    """Criterion 2 (upper): t <= cell_m / 2v."""
    cell_m = SPEED_OF_LIGHT / (n_pilots * k_step * SCS)
    denom = 2.0 * v
    return cell_m / denom


def span_accel(exp=ACCEL_EXP):
    """Criterion 3 (upper): (lambda / 2a)^exp. Independent of v."""
    return (lambda_() / (2.0 * MAX_ACCELERATION)) ** exp


def span_min(dv):
    """Criterion 4 (lower): t >= lambda / (2 dv). Independent of v."""
    return lambda_() / (2.0 * dv)


def t_max(v, accel_exp=ACCEL_EXP, n_pilots=N_PILOTS):
    """Effective upper bound: the smallest of the three upper criteria."""
    return np.minimum.reduce([
        span_grid(v),
        span_migration(v, n_pilots=n_pilots),
        np.full_like(np.asarray(v, dtype=float), span_accel(accel_exp)),
    ])


# ---------------------------------------------------------------------------
# Figure 1: t_span vs velocity coverage
# ---------------------------------------------------------------------------
def plot_window_vs_velocity(accel_exp=ACCEL_EXP):
    v = np.linspace(0.5, 30.0, 600)

    grid = span_grid(v)
    mig  = span_migration(v)
    acc  = np.full_like(v, span_accel(accel_exp))
    tmax = np.minimum.reduce([grid, mig, acc])

    fig, ax = plt.subplots(figsize=(9, 6))

    # Upper bounds
    ax.plot(v, grid, color="#4C78A8", lw=1.6, label="t_span^max (grid budget)")
    ax.plot(v, mig,  color="#F58518", lw=1.6,
            label="t_span^max (migration)")
    ax.plot(v, acc,  color="#54A24B", lw=1.6,
            label=f"t_span^max (acceleration)")
    ax.plot(v, tmax, color="black", lw=2.4, label="min(t_span^max)")

    # Lower bounds (resolution targets) + feasible shading per target
    reds = ["#E45756", "#B279A2", "#9D755D"]
    for dv, c in zip(DV_TARGETS, reds):
        tmin = span_min(dv)
        ax.axhline(tmin, color=c, ls="--", lw=1.4,
                   label=f"dv = {dv} m/s  (t_min = {tmin*1e3:.0f} ms)")

    ax.set_xlabel("velocity coverage  v  [m/s]")
    ax.set_ylabel("observation window  t_span  [s]")
    ax.set_title("Feasible observation window vs max target velocity\n"
                 f"{CENTER_FREQUENCY/1e9:.2f} GHz, {SCS//1000} kHz SCS, "
                 f"TDD Period={TDD_PERIOD_S*1e3:.0f} ms, N_pilots={N_PILOTS}, K_step={K_STEP}, max accel={MAX_ACCELERATION} m/s²")
    ax.set_ylim(0, 0.6)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure 2: feasibility map over (velocity coverage, resolution target)
#   This is the map_bound* idea: sweep the bounds over a parameter grid.
#   Colour = feasible window width (t_max - t_min), masked where infeasible.
# ---------------------------------------------------------------------------
def plot_feasibility_map(accel_exp=ACCEL_EXP):
    v  = np.linspace(0.5, 30.0, 400)
    dv = np.linspace(0.05, 0.5, 400)
    V, DV = np.meshgrid(v, dv)

    TMAX = t_max(V, accel_exp)
    TMIN = span_min(DV)
    width = TMAX - TMIN               # >0 feasible, <=0 infeasible
    width_ms = np.where(width > 0, width * 1e3, np.nan)

    fig, ax = plt.subplots(figsize=(9, 6))
    pcm = ax.pcolormesh(V, DV, width_ms, shading="auto", cmap="viridis")
    cb = fig.colorbar(pcm, ax=ax)
    cb.set_label("feasible window width  t_max - t_min  [ms]")

    # Feasibility boundary (width = 0)
    ax.contour(V, DV, width, levels=[0.0], colors="red", linewidths=2.0)

    ax.set_xlabel("velocity coverage  v  [m/s]")
    ax.set_ylabel("velocity resolution  dv  [m/s]")
    ax.set_title("Feasibility map")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure 3: grant-width (N_pilots) sweep.
#   Only range migration depends on N_pilots: cell_m = c/(N_pilots k_step scs),
#   so a narrower grant widens the cell and LOOSENS the migration bound -- at the
#   cost of range resolution and SNR. This shows how much headroom that buys.
# ---------------------------------------------------------------------------
def plot_npilots_sweep(dv=0.3, accel_exp=ACCEL_EXP):
    v = np.linspace(0.5, 30.0, 600)
    # Grants from a few PRB up to the full 273-PRB carrier (6 type-1 pilots/PRB).
    prb_list = [12, 51, 106, 273]
    cmap = plt.cm.plasma(np.linspace(0.15, 0.85, len(prb_list)))

    fig, ax = plt.subplots(figsize=(9, 6))
    tmin = span_min(dv)

    for prb, c in zip(prb_list, cmap):
        n_pilots = prb * 6
        tmax = t_max(v, accel_exp=accel_exp,
                     n_pilots=n_pilots)
        ax.plot(v, tmax, color=c, lw=1.8, label=f"N_pilots={n_pilots}")
        feas = tmax > tmin
        if feas.any():
            ax.plot(v[feas].max(), tmin, "o", color=c, ms=7)

    ax.axhline(tmin, color="black", ls="--", lw=1.4,
               label=f"dv={dv} m/s (t_min={tmin*1e3:.0f} ms)")

    ax.set_xlabel(f"velocity coverage  v  [m/s]")
    ax.set_ylabel(f"min(t_span^max)  [s]")
    ax.set_title(f"Grant width vs feasibility (dv={dv} m/s)")
    ax.set_ylim(0, 0.35)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    return fig


def plot_simple(dv=0.3):
    """The whole story in two lines: window you NEED vs window you're ALLOWED."""
    v = np.linspace(0.0, 20.0, 400)

    need    = np.full_like(v, span_min(dv) * 1e3)   # ms, flat
    allowed = t_max(np.maximum(v, 1e-6)) * 1e3       # ms, falls with speed

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(v, need, color="#D62728", lw=2.5,
            label=f"window size needed for {dv} m/s resolution")
    ax.plot(v, allowed, color="#1F77B4", lw=2.5,
            label="longest window allowed wrt criterions")

    ok = allowed >= need
    ax.fill_between(v, need, allowed, where=ok, color="#2CA02C", alpha=0.18)
    ax.fill_between(v, allowed, need, where=~ok, color="#D62728", alpha=0.12)

    # The crossing: where "allowed" falls below "needed"
    if ok.any() and (~ok).any():
        v_cross = v[ok].max()
        ax.axvline(v_cross, color="black", ls=":", lw=1.2)
        ax.text(v_cross / 2, need[0] * 0.5, "ALLOWED", ha="center",
                fontsize=13, color="#2CA02C", weight="bold")
        ax.text((v_cross + 20) / 2, need[0] * 0.5, "IMPOSSIBLE", ha="center",
                fontsize=13, color="#D62728", weight="bold")

    ax.set_xlabel("maximum velocity [m/s]")
    ax.set_ylabel("window size [ms]")
    ax.set_title(f"Can we resolve {dv} m/s?  Only below {v_cross:.0f} m/s.")
    ax.set_xlim(0, 20)
    ax.set_ylim(0, max(need[0], allowed.max()) * 1.4)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    # f0 = plot_simple(dv=0.3)
    # f0.savefig("bounds_simple.png", dpi=130)
    # f1 = plot_window_vs_velocity()
    # f2 = plot_feasibility_map()
    f3 = plot_npilots_sweep(dv=0.3)
    # f1.savefig("bounds_window_vs_velocity.png", dpi=130)
    # f2.savefig("bounds_feasibility_map.png", dpi=130)
    # f3.savefig("bounds_npilots_sweep.png", dpi=130)
    print("wrote bounds_window_vs_velocity.png, bounds_feasibility_map.png, "
          "bounds_npilots_sweep.png")
    plt.show()
