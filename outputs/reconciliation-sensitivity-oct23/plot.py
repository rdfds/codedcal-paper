"""Standalone scientific sensitivity plot; does not replace a main figure."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from sensitivity import *


def run():
    results = [json.loads((P/f"result-{r['tag']}.json").read_text()) for r in registry()]
    gain = np.array([[a['full_to_diagonal_fit_gap']-a['selected_to_diagonal_fit_gap']
        if a['selected_to_diagonal_fit_gap'] is not None else np.nan for a in r['arms']] for r in results])
    coverage = np.array([[100*a['localization_coverage'] for a in r['arms']] for r in results])
    labels = [a['name'].replace('-', '\n') for a in results[0]['arms']]
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.2), constrained_layout=True)
    extent = max(float(np.nanmax(abs(gain))), 1e-6)
    image = axes[0].imshow(gain, cmap='RdBu', norm=TwoSlopeNorm(0, -extent, extent))
    fig.colorbar(image, ax=axes[0], shrink=.7, label='Full gap minus matched gap (nm²)')
    image2 = axes[1].imshow(coverage, cmap='viridis', vmin=0, vmax=100)
    fig.colorbar(image2, ax=axes[1], shrink=.7, label='Selected code-input localization rows (%)')
    axes[0].set_title('Population matching can narrow or widen gaps')
    axes[1].set_title('Reference selection changes across the grid')
    for ax in axes:
        ax.set_xticks(range(9), labels, fontsize=8)
        ax.set_yticks(range(9), [r['tag'] for r in results], fontsize=9)
        ax.axvline(2.5, color='white', linewidth=1); ax.axvline(5.5, color='white', linewidth=1)
        ax.set_xlabel('Radius factor × eligibility; original setting is r1-original')
    for y in range(9):
        for x in range(9):
            arm = results[y]['arms'][x]
            supported = all(arm[k]['uncertainty']['all_frequencies_usable']
                for k in ('selected_code', 'diagonal'))
            label = (f'{gain[y,x]:.2f}'+('' if supported else '*')
                if np.isfinite(gain[y,x]) else 'NA')
            axes[0].text(x, y, label,
                ha='center', va='center', fontsize=7,
                color='white' if abs(gain[y,x]) > extent*.55 else 'black')
            axes[1].text(x, y, f'{coverage[y,x]:.1f}', ha='center', va='center',
                fontsize=7, color='white' if coverage[y,x] < 45 else 'black')
    fig.suptitle('Eight acquisitions, nine views, all 81 specified settings; post-hoc conditional sensitivity', fontsize=12)
    fig.supxlabel('Descriptive point-fit gaps are not accuracy or equivalence; * marks unsupported individual Fieller sets.', fontsize=9)
    fig.savefig(P/'sensitivity.png', dpi=180)
    fig.savefig(P/'sensitivity.pdf')
    plt.close(fig)


if __name__ == '__main__':
    run()
