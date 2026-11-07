"""Main figures for the RESI design analysis, from saved results only.

Usage: python figures.py        (writes figures/fig1-4 as PDF and PNG)
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from design import uniform_distinct

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'
FIG = HERE / 'figures'

INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#8a8984', '#e4e3df'
BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
MAGENTA, VIOLET = '#e87ba4', '#6250d6'
LABEL_COLORS = [BLUE, ORANGE, AQUA, YELLOW]
UNTREATED = ['fig4untreated', 'edf10rep1', 'edf10rep2']
TREATED = ['rtx_rep1', 'rtx_rep2', 'rtx_rep3']
SHORT = {'fig4untreated': 'U1', 'edf10rep1': 'U2', 'edf10rep2': 'U3',
         'rtx_rep1': 'R1', 'rtx_rep2': 'R2', 'rtx_rep3': 'R3'}
MM = 1 / 25.4

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7, 'axes.titlesize': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.edgecolor': MUTED,
    'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'lines.linewidth': 1.2,
    'lines.solid_capstyle': 'round', 'axes.spines.top': False, 'axes.spines.right': False,
    'legend.frameon': False, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
})


def load(name):
    return json.loads((RES / name).read_text())


def grid(ax, axis='y'):
    ax.grid(True, axis=axis, color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)


def panel(ax, letter, title):
    ax.set_title(title, loc='left', fontweight='normal', color=INK, pad=4)
    ax.text(-0.18, 1.04, letter, transform=ax.transAxes, fontweight='bold', fontsize=8,
            ha='left', va='bottom')


def save(fig, name):
    FIG.mkdir(exist_ok=True)
    for ext, meta in (('pdf', {'CreationDate': None}), ('png', None)):
        fig.savefig(FIG / f'{name}.{ext}', dpi=300, bbox_inches='tight', facecolor='white', metadata=meta)
    plt.close(fig)


def dots_ci(ax, x, val, lo, hi, color, size=4.5, **kw):
    ax.vlines(x, lo, hi, color=color, linewidth=1.2)
    ax.plot(x, val, 'o', color=color, markersize=size, markeredgecolor='white',
            markeredgewidth=0.8, zorder=3, **kw)


# ---------------------------------------------------------------- Figure 1
def fig1():
    fig = plt.figure(figsize=(180 * MM, 58 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1, 1], wspace=0.45)

    ax = fig.add_subplot(gs[0])
    ax.set_title('Random permanent labels', loc='left', color=INK, pad=4)
    ax.text(-0.06, 1.0, 'a', transform=ax.transAxes, fontweight='bold', fontsize=8, ha='right', va='bottom')
    rng = np.random.default_rng(11)
    centres = np.array([[1.2, 3.3], [3.4, 3.8], [2.4, 1.6], [4.6, 1.5], [0.9, 0.8], [4.4, 3.0]])
    kinds = [2, 1, 3, 2, 1, 1]           # complex sizes
    lab_seq = iter(rng.choice(4, 20))
    for (cx, cy), m in zip(centres, kinds):
        th = rng.uniform(0, np.pi)
        pts = [(cx + (j - (m - 1) / 2) * 0.55 * np.cos(th), cy + (j - (m - 1) / 2) * 0.55 * np.sin(th))
               for j in range(m)]
        labs = [next(lab_seq) for _ in pts]
        for i in range(m - 1):
            ax.plot(*zip(pts[i], pts[i + 1]), color=MUTED, linewidth=0.8, zorder=1)
        for (x, y), l in zip(pts, labs):
            ax.plot(x, y, 'o', markersize=7, color=LABEL_COLORS[l], markeredgecolor='white',
                    markeredgewidth=1.0, zorder=2)
    ax.set_xlim(0.2, 5.4); ax.set_ylim(0.2, 4.6); ax.set_aspect('equal'); ax.set_anchor('N'); ax.axis('off')
    for i, c in enumerate(LABEL_COLORS):
        ax.plot(0.45 + i * 1.2, -0.25, 'o', color=c, markersize=5, markeredgecolor='white', clip_on=False)
        ax.text(0.62 + i * 1.2, -0.25, f'R{i + 1}', va='center', color=INK2, fontsize=6.5)
    ax.text(0.2, -0.75, 'Different labels imply different molecules', color=INK2, fontsize=6.5)

    ax = fig.add_subplot(gs[1])
    panel(ax, 'b', 'Labels distinct across m molecules')
    K = np.arange(2, 17)
    for m, c in zip((2, 3, 4), (BLUE, ORANGE, AQUA)):
        y = [uniform_distinct(m, k) for k in K]
        ax.plot(K, y, color=c, linewidth=1.5, label=f'm = {m}')
        ax.text(16.3, y[-1], f'm = {m}', va='center', fontsize=6.5, color=INK2)
    for m in (2, 3, 4):
        ax.plot(4, uniform_distinct(m, 4), 'o', color=INK, markersize=3.5, zorder=3)
    ax.annotate('RESI, K = 4: 75%, 38%, 9%', xy=(4, uniform_distinct(4, 4)), xytext=(4.6, 0.03),
                fontsize=6.5, color=INK2, arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.6))
    ax.set_xlabel('Labels K'); ax.set_ylabel('Probability'); ax.set_xlim(2, 16); ax.set_ylim(0, 1)
    ax.set_xticks([2, 4, 8, 12, 16]); grid(ax, 'both')
    ax.legend(loc='center right', bbox_to_anchor=(1.0, 0.36), handlelength=1.2)

    ax = fig.add_subplot(gs[2])
    panel(ax, 'c', 'Exact identities, simulated label draws')
    v = load('verify_exact.json')
    zs = []
    for case in v['cases']:
        for key in ('z_same', 'z_cross', 'z_S', 'z_W'):
            zs.extend([z for z in case[key] if z is not None and np.isfinite(z)])
        zs.append(case['z_star'])
    zs = np.array(zs)
    bins = np.arange(-3, 3.01, 0.5)
    ax.hist(zs, bins=bins, color=BLUE, edgecolor='white', linewidth=0.8)
    xx = np.linspace(-3, 3, 200)
    ax.plot(xx, len(zs) * 0.5 * np.exp(-xx ** 2 / 2) / np.sqrt(2 * np.pi), color=INK2, linewidth=1.0)
    top = ax.get_ylim()[1] * 1.45
    ax.set_ylim(0, top)
    ax.text(1.0, 0.62 * top, 'N(0, 1)', color=INK2, fontsize=6.5)
    s = v['summary']
    ax.text(0.02, 0.98, f"{s['n_z']} comparisons\nmax |z| = {s['max_abs_z']:.2f}\n"
            f"variance ratio {s['var_ratio_min']:.2f}–{s['var_ratio_max']:.2f}",
            transform=ax.transAxes, va='top', fontsize=6.5, color=INK2)
    ax.set_xlabel('(Monte Carlo − exact) / s.e.'); ax.set_ylabel('Comparisons'); ax.set_xlim(-3, 3)
    grid(ax)
    save(fig, 'fig1')


# ---------------------------------------------------------------- Figure 2
def fig2():
    fig, axes = plt.subplots(1, 2, figsize=(140 * MM, 58 * MM), gridspec_kw=dict(wspace=0.4))
    ax = axes[0]
    panel(ax, 'a', 'Excess neighbours, simulated fields')
    rows = load('verify_model.json')
    jitter = {}
    for r in rows:
        k = r['true_f']; jitter[k] = jitter.get(k, -1) + 1
        x = r['realized_f'] + (jitter[k] - 3.5) * 0.006
        dots_ci(ax, x, r['excess'], *r['excess_ci'], BLUE, size=3.5)
    ax.plot([0, 0.7], [0, 0.7], color=MUTED, linewidth=0.8, zorder=0)
    ax.text(0.17, 0.07, 'identity', color=INK2, fontsize=6.5, rotation=36)
    ax.set_xlabel('True fraction of subunits in dimers'); ax.set_ylabel('Estimate (95% CI)')
    ax.set_xlim(-0.05, 0.72); ax.set_ylim(-0.1, 0.92); grid(ax, 'both')

    ax = axes[1]
    panel(ax, 'b', 'Connected triplets, simulated fields')
    rows = load('verify_triplets.json')
    order = ['monomer only, heterogeneous', 'monomer/dimer', 'monomer/dimer, heterogeneous',
             'with trimers', 'chains up to 5, heterogeneous']
    names = ['Mono*', 'Dimer', 'Dimer*', 'Trimer', 'Chain*']
    for i, case in enumerate(order):
        rr = [r for r in rows if r['case'] == case]
        for j, r in enumerate(rr):
            x = i + (j - 1) * 0.18
            dots_ci(ax, x, r['m2'], *r['ci'], ORANGE, size=3.5)
            ax.plot([x - 0.08, x + 0.08], [r['oracle']] * 2, color=INK, linewidth=1.0, zorder=4)
    ax.axhline(0, color=MUTED, linewidth=0.6, zorder=0)
    ax.set_xticks(range(len(order))); ax.set_xticklabels(names, fontsize=6)
    ax.set_ylabel(r'$m_2$ (95% CI)'); grid(ax)
    ax.plot([], [], 'o', color=ORANGE, markersize=3.5, label='estimate')
    ax.plot([], [], '-', color=INK, linewidth=1.0, label='true value')
    ax.legend(loc='upper left', handlelength=1.0)
    ax.text(1.0, -0.2, '* uneven density', transform=ax.transAxes, ha='right', va='top', fontsize=6, color=INK2)
    save(fig, 'fig2')


# ---------------------------------------------------------------- Figure 3
def fig3():
    fig = plt.figure(figsize=(180 * MM, 58 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1], wspace=0.42)
    ax = fig.add_subplot(gs[0])
    panel(ax, 'a', 'Cross-label correlation, RESI molecules')
    for views, color, name in ((UNTREATED, BLUE, 'untreated'), (TREATED, ORANGE, 'rituximab')):
        for k, v in enumerate(views):
            d = np.load(RES / f'cd20-molecules-{v}.npz')
            e = d['edges']; c = (e[:-1] + e[1:]) / 2
            lo, hi = np.percentile(d['gboot'], [2.5, 97.5], axis=0)
            ax.fill_between(c, lo, hi, color=color, alpha=0.1, linewidth=0)
            ax.plot(c, d['g'], color=color, linewidth=1.2, label=name if k == 0 else None)
    ax.axhline(1, color=MUTED, linewidth=0.6, zorder=0)
    ax.set_xlim(0, 100); ax.set_yscale('log'); ax.set_ylim(0.8, 30)
    ax.set_yticks([1, 2, 5, 10, 20]); ax.set_yticklabels(['1', '2', '5', '10', '20'])
    ax.set_xlabel('Distance between molecules with different labels (nm)'); ax.set_ylabel('g(r)')
    grid(ax, 'both'); ax.legend(loc='upper right', handlelength=1.2)

    mol = load('cd20-molecules.json')
    trip = load('cd20-triplets.json')
    views = UNTREATED + TREATED
    ax = fig.add_subplot(gs[1])
    panel(ax, 'b', 'Excess labelled neighbours ≤ 30 nm')
    for i, v in enumerate(views):
        val, (lo, hi) = mol[v]['excess_neighbours_30nm']
        dots_ci(ax, i, val, lo, hi, BLUE if v in UNTREATED else ORANGE)
    for y, x, txt in ((0.5, 2.55, 'dimer maximum,\n50% labelling'), (1.0, -0.45, 'dimer maximum, any labelling')):
        ax.axhline(y, color=MUTED, linewidth=0.6, zorder=0)
        ax.text(x, y - 0.03, txt, fontsize=6, color=INK2, va='top')
    ax.set_xticks(range(6)); ax.set_xticklabels([SHORT[v] for v in views])
    ax.set_xlim(-0.6, 5.6); ax.set_ylim(0, 2.45); ax.set_ylabel('Excess neighbours (95% CI)'); grid(ax)

    ax = fig.add_subplot(gs[2])
    panel(ax, 'c', 'Connected triplets ≤ 30 nm')
    for i, v in enumerate(views):
        dots_ci(ax, i, trip[v]['m2'], *trip[v]['m2_ci'], BLUE if v in UNTREATED else ORANGE)
    ax.axhline(0, color=MUTED, linewidth=0.6, zorder=0)
    ax.text(2.6, 0.12, 'zero for any\nmonomer/dimer\nmixture', fontsize=6, color=INK2, va='bottom')
    ax.set_xticks(range(6)); ax.set_xticklabels([SHORT[v] for v in views])
    ax.set_xlim(-0.6, 5.6); ax.set_ylim(-0.1, 3.3); ax.set_ylabel(r'Partner pairs per molecule, $m_2$ (95% CI)')
    grid(ax)
    save(fig, 'fig3')


# ---------------------------------------------------------------- Figure 4
def fig4():
    fig = plt.figure(figsize=(180 * MM, 112 * MM))
    gs = fig.add_gridspec(2, 2, hspace=0.6, wspace=0.32)
    views = UNTREATED + TREATED

    ax = fig.add_subplot(gs[0, 0])
    panel(ax, 'a', 'Label permanence: same molecule in two rounds?')
    perm = load('cd20-permanence.json')
    mol = load('cd20-molecules.json')
    for i, v in enumerate(views):
        r = perm['views'][v]['translation_corrected']
        obs, null = np.array(r['observed_by_nm'], float), np.array(r['null_by_nm'], float)
        ref = obs[2:8].sum() / null[2:8].sum()           # organization-driven enrichment at 2-8 nm
        queries = 3 * sum(mol[v]['molecules'])           # nearest neighbour in each other round
        o2, e2 = obs[:2].sum(), ref * null[:2].sum()
        val = 100 * (o2 - e2) / queries
        half = 100 * 1.96 * np.sqrt(o2 + ref ** 2 * null[:2].sum()) / queries
        dots_ci(ax, i, val, val - half, val + half, BLUE if v in UNTREATED else ORANGE)
    ax.axhline(0, color=MUTED, linewidth=0.6, zorder=0)
    ax.set_xticks(range(6)); ax.set_xticklabels([SHORT[v] for v in views]); ax.set_xlim(-0.6, 5.6)
    ax.set_ylim(-0.3, 0.45)
    ax.set_ylabel('Excess cross-round neighbours\nwithin 2 nm (% of queries)'); grid(ax)

    ax = fig.add_subplot(gs[0, 1])
    panel(ax, 'b', 'Per-round registration shifts (no fiducials)')
    der = load('cd20-derived.json')
    for i, v in enumerate(views):
        t = np.array(der[v]['round_translation_nm'])
        for r in range(4):
            ax.plot(i + (r - 1.5) * 0.14, np.hypot(*t[r]), 'o', color=LABEL_COLORS[r], markersize=4.5,
                    markeredgecolor='white', markeredgewidth=0.8)
    for r, col in enumerate(LABEL_COLORS):
        ax.plot([], [], 'o', color=col, markersize=4, label=f'R{r + 1}')
    ax.set_xticks(range(6)); ax.set_xticklabels([SHORT[v] for v in views]); ax.set_xlim(-0.6, 5.6)
    ax.set_ylim(0, 1.8); ax.set_ylabel('Shift magnitude (nm)'); grid(ax)
    ax.legend(loc='upper left', ncol=4, handletextpad=0.2, columnspacing=0.8)

    bud = load('cd20-budget.json')
    phis = np.array(bud['phis'])
    th_colors = {'min16': BLUE, 'min8': ORANGE, 'min4': AQUA}
    th_names = {'min16': '≥ 16', 'min8': '≥ 8', 'min4': '≥ 4'}
    ax = fig.add_subplot(gs[1, 0])
    panel(ax, 'c', 'Detection versus round length')
    for key, col in th_colors.items():
        P = np.array([[bud['views'][v][key]['p_det'][str(p)] for p in phis] for v in views])
        med = np.median(P, 0)
        ax.fill_between(phis, P.min(0), P.max(0), color=col, alpha=0.1, linewidth=0)
        ax.plot(phis, med, color=col, linewidth=1.5, marker='o', markersize=3, label=th_names[key])
    ax.set_xscale('log', base=2); ax.set_xticks(phis); ax.set_xticklabels(['1', '1/2', '1/4', '1/8', '1/16'])
    ax.set_xlim(1.1, 0.055); ax.set_ylim(0, 1.02)
    ax.set_xlabel('Round length relative to published rounds'); ax.set_ylabel('Molecules detected')
    grid(ax, 'both'); ax.legend(loc='lower left', handlelength=1.2, title='Localizations\nper molecule',
                                title_fontsize=6.5)

    ax = fig.add_subplot(gs[1, 1])
    panel(ax, 'd', 'Triplet information at fixed imaging time')
    Ks = bud['K']
    for key, col in th_colors.items():
        I = np.array([[bud['views'][v][key]['order3'][str(k)] for k in Ks] for v in views])
        med = np.median(I, 0)
        ax.fill_between(Ks, I.min(0), I.max(0), color=col, alpha=0.1, linewidth=0)
        ax.plot(Ks, med, color=col, linewidth=1.5, marker='o', markersize=3, label=th_names[key])
    ax.axhline(1, color=MUTED, linewidth=0.6, zorder=0)
    ax.set_xlim(2, 16); ax.set_xticks([2, 4, 6, 8, 12, 16]); ax.set_ylim(0, 1.75)
    ax.set_xlabel('Rounds K sharing the imaging time of four full rounds')
    ax.set_ylabel('Triplet information\n(relative to four full rounds)'); grid(ax, 'both')
    ax.text(0.98, 0.97, 'colours as in c', transform=ax.transAxes, ha='right', va='top', fontsize=6, color=INK2)
    save(fig, 'fig4')


# ---------------------------------------------------------------- Figure 5
def fig5():
    v = load('npc-validation.json')
    fig = plt.figure(figsize=(180 * MM, 58 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.05, 0.95], wspace=0.75)

    ax = fig.add_subplot(gs[0])
    panel(ax, 'a', 'Nup96 corners, cross-label neighbours')
    e = np.array(v['nn_edges']); c = (e[:-1] + e[1:]) / 2
    ax.bar(c, v['observed_nn_hist'], width=3.4, color=BLUE, edgecolor='white', linewidth=0.8, label='data')
    ax.axvline(v['radius_nm'], color=MUTED, linewidth=0.6)
    top = 1.3 * max(v['observed_nn_hist'])
    ax.set_ylim(0, top)
    ax.text(v['radius_nm'], 0.99 * top, f"R = {v['radius_nm']:g} nm", fontsize=6, color=INK2, va='top', ha='center',
            bbox=dict(facecolor='white', edgecolor='none', pad=1))
    ax.text(12, 0.82 * top, 'same\ncorner', fontsize=6, color=INK2, ha='center', va='bottom')
    ax.text(42, 0.82 * top, 'adjacent\ncorner', fontsize=6, color=INK2, ha='center', va='bottom')
    ax.set_xlabel('Nearest other-label molecule, same pore (nm)'); ax.set_ylabel('Fraction of molecules')
    ax.set_xlim(0, 62); grid(ax)

    ax = fig.add_subplot(gs[1])
    panel(ax, 'b', 'Efficiency, 1,218 NPCs')
    d = v['data']
    routes = [('p_counts', 'counts / 32', INK2),
              ('p_pairs', 'pairs', BLUE),
              ('p_triples', 'triples', BLUE),
              ('p_joint', 'pairs + triples,\nn free', BLUE),
              ('p_corner_pairs', 'corner partners', ORANGE),
              ('p_corner_triplets', 'corner pairs', ORANGE),
              ('p_corner_ratio', 'corner ratio', ORANGE)]
    for i, (k, name, col) in enumerate(routes):
        lo, hi = d['ci'][k]
        ax.hlines(i, lo, hi, color=col, linewidth=1.2)
        ax.plot(d[k], i, 'o', color=col, markersize=4.5, markeredgecolor='white', markeredgewidth=0.8, zorder=3)
    ax.axvline(d['p_counts'], color=GRID, linewidth=0.8, zorder=0)
    ax.set_yticks(range(len(routes))); ax.set_yticklabels([r[1] for r in routes], fontsize=6)
    ax.invert_yaxis(); ax.set_xlim(0.2, 0.45); ax.set_xlabel('Efficiency p (95% CI)'); grid(ax, 'x')
    ax.text(0.02, 0.98, f"joint copies per pore\nn = {d['n_joint']:.0f} "
            f"({d['ci']['n_joint'][0]:.0f}–{d['ci']['n_joint'][1]:.0f});\nstructure: 32",
            transform=ax.transAxes, ha='left', va='top', fontsize=6, color=INK2)

    ax = fig.add_subplot(gs[2])
    panel(ax, 'c', 'Recovery in simulated NPCs')
    sims = v['simulations']
    keys = [('p_pairs', BLUE, 'o'), ('p_triples', BLUE, 's'), ('p_corner_ratio', ORANGE, 'D')]
    for s_ in sims:
        off = -0.008 if s_['merge_nm'] == 0 else 0.008
        for k, col, mk in keys:
            merged = s_['merge_nm'] > 0
            ax.plot(s_['p_true'] + off, s_[k], mk, color=col, markersize=3.5,
                    markeredgecolor=col if merged else 'white', markeredgewidth=0.8 if merged else 0.6,
                    fillstyle='none' if merged else 'full')
    ax.plot([0.25, 0.55], [0.25, 0.55], color=MUTED, linewidth=0.8, zorder=0)
    for k, col, mk in keys:
        ax.plot([], [], mk, color=col, markersize=3.5, label={'p_pairs': 'pairs', 'p_triples': 'triples',
                                                               'p_corner_ratio': 'corner ratio'}[k])
    ax.legend(loc='upper left', handlelength=1.0)
    ax.text(0.98, 0.02, 'filled: no merging\nopen: 8-nm same-label\nmerging', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=6, color=INK2)
    ax.set_xlim(0.25, 0.55); ax.set_ylim(0.25, 0.55); ax.set_xticks([0.3, 0.4, 0.5]); ax.set_yticks([0.3, 0.4, 0.5])
    ax.set_xlabel('True efficiency'); ax.set_ylabel('Estimate'); grid(ax, 'both')
    save(fig, 'fig5')


# ---------------------------------------------------------------- Figure 6
def fig6():
    v = load('cd20-shape.json')
    e = np.rad2deg(np.array(v['angle_edges'])); c = (e[:-1] + e[1:]) / 2
    fig = plt.figure(figsize=(180 * MM, 58 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.9], wspace=0.42)
    shown = [('straight chains', BLUE), ('flexible chains, 60 deg', ORANGE), ('compact clusters', AQUA),
             ('hexameric rings', VIOLET)]
    names = {'straight chains': 'straight chains', 'flexible chains, 60 deg': 'flexible chains (60°)',
             'compact clusters': 'compact clusters', 'hexameric rings': 'hexameric rings'}

    ax = fig.add_subplot(gs[0])
    panel(ax, 'a', 'Reference shapes (simulated)')
    for k, col in shown:
        ax.plot(c, v['models'][k], color=col, linewidth=1.5, label=names[k])
    ax.set_xlim(0, 180); ax.set_xticks([0, 60, 120, 180]); ax.set_ylim(0, 0.2)
    ax.set_xlabel('Angle between two partners (°)'); ax.set_ylabel('Share of connected triplets')
    grid(ax, 'both'); ax.legend(loc='upper center', handlelength=1.2, fontsize=6)

    ax = fig.add_subplot(gs[1])
    panel(ax, 'b', 'Rituximab-treated CD20')
    for k, col in shown[:3]:
        ax.plot(c, v['models'][k], color=col, linewidth=1.0, alpha=0.5)
    for i, view in enumerate(TREATED):
        o = v['views'][view]
        y, se = np.array(o['normalized']), np.array(o['normalized_se'])
        ax.errorbar(c + (i - 1) * 2.5, y, yerr=1.96 * se, fmt='o', color=INK, markersize=2.8, linewidth=0.8,
                    elinewidth=0.8, capsize=0, label=SHORT[view])
    ax.set_xlim(0, 180); ax.set_xticks([0, 60, 120, 180]); ax.set_ylim(0, 0.2)
    ax.set_xlabel('Angle between two partners (°)'); grid(ax, 'both')
    ax.text(0.98, 0.97, 'black: R1-R3 (95% CI)\nlines: models as in a', transform=ax.transAxes, ha='right',
            va='top', fontsize=6, color=INK2)

    ax = fig.add_subplot(gs[2])
    panel(ax, 'c', 'Linearity index')
    views = UNTREATED[:2] + TREATED
    for i, view in enumerate(views):
        o = v['views'][view]
        dots_ci(ax, i, o['linearity'], *o['linearity_ci'], BLUE if view in UNTREATED else ORANGE)
    for k, col in shown:
        y = v['model_linearity'][k]
        ax.axhline(y, color=col, linewidth=0.8, zorder=0)
        va = {'flexible chains, 60 deg': 'bottom', 'compact clusters': 'top'}.get(k, 'center')
        ax.text(4.45, y, names[k], fontsize=5.5, color=INK2, va=va, ha='left')
    ax.set_xticks(range(len(views))); ax.set_xticklabels([SHORT[x] for x in views]); ax.set_xlim(-0.5, 4.4)
    ax.set_ylim(0, 0.55); ax.set_ylabel('Share of triplets above 150°'); grid(ax)
    save(fig, 'fig6')


# ---------------------------------------------------------------- Supplementary robustness
def figS1():
    r = load('robustness.json')
    fig = plt.figure(figsize=(180 * MM, 58 * MM))
    gs = fig.add_gridspec(1, 3, wspace=0.42)
    views = UNTREATED + TREATED
    for j, (stat, title, ylab, kfmt, opts) in enumerate((
            ('excess', 'Excess neighbours vs radius', 'Excess labelled neighbours', 'r{:g}_b60-100', (20, 25, 30, 40)),
            ('m2', 'Connected triplets vs radius', r'$m_2$', 'r{:g}_a100-250', (20, 25, 30, 40)))):
        ax = fig.add_subplot(gs[j])
        panel(ax, 'ab'[j], title)
        for v in views:
            col = BLUE if v in UNTREATED else ORANGE
            y = [r['cd20'][v][stat][kfmt.format(x)][0] for x in opts]
            lo = [r['cd20'][v][stat][kfmt.format(x)][1][0] for x in opts]
            hi = [r['cd20'][v][stat][kfmt.format(x)][1][1] for x in opts]
            ax.fill_between(opts, lo, hi, color=col, alpha=0.1, linewidth=0)
            ax.plot(opts, y, color=col, linewidth=1.2, marker='o', markersize=3)
            if v in TREATED:
                ax.text(41, y[-1], SHORT[v], fontsize=6, color=INK2, va='center')
        u_end = max(r['cd20'][v][stat][kfmt.format(opts[-1])][0] for v in UNTREATED)
        ax.text(41, u_end, 'U1–U3', fontsize=6, color=INK2, va='bottom')
        ax.set_xlabel('Radius (nm)'); ax.set_ylabel(ylab); ax.set_xticks(opts); ax.set_xlim(18, 44); grid(ax, 'both')
        if stat == 'excess':
            ax.axhline(0.5, color=MUTED, linewidth=0.6, zorder=0)
            ax.text(31, 0.53, 'dimer maximum\n(50% labelling)', fontsize=5.5, color=INK2, va='bottom')
    ax = fig.add_subplot(gs[2])
    panel(ax, 'c', 'NPC efficiency vs analysis choices')
    a = r['npc']['assignment']
    for i, key in enumerate(('p_counts', 'p_pairs', 'p_triples')):
        ax.plot([float(x) for x in a], [a[x][key] for x in a], marker='o', markersize=3, linewidth=1.2,
                color=[INK2, BLUE, AQUA][i], label=['counts / 32', 'pairs', 'triples'][i])
    cr = r['npc']['corner_range']
    ax.axhspan(min(v[0] for v in cr.values()), max(v[1] for v in cr.values()), color=ORANGE, alpha=0.12, linewidth=0)
    ax.text(72, 0.445, 'corner-level range over\ngeometry assumptions', fontsize=6, color=INK2, va='top')
    ax.set_xlabel('Pore assignment radius (nm)'); ax.set_ylabel('Efficiency p'); ax.set_ylim(0.28, 0.48)
    ax.set_xticks([70, 90, 110]); grid(ax, 'both'); ax.legend(loc='lower right', handlelength=1.2)
    save(fig, 'figS1')


if __name__ == '__main__':
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6(); figS1()
    print('wrote', sorted(p.name for p in FIG.glob('*.png')))
