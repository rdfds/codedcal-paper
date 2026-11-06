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
    for ext in ('pdf', 'png'):
        fig.savefig(FIG / f'{name}.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
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


if __name__ == '__main__':
    fig1(); fig2(); fig3(); fig4()
    print('wrote', sorted(p.name for p in FIG.glob('*.png')))
