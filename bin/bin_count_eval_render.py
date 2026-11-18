"""Render frozen ordination coordinates and original-space distance scores."""
from pathlib import Path
import os


def pyplot():
    os.environ.setdefault('MPLCONFIGDIR', str(Path(os.environ.get('TMPDIR', '/tmp'))/'coverage-bin-eval-matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans', 'pdf.fonttype':42, 'ps.fonttype':42,
                         'axes.spines.top':False, 'axes.spines.right':False, 'font.size':10})
    return plt


def scatter(coordinates, rows, folder, title, method):
    import numpy as np
    from bin.bin_count_eval_io import write_tsv
    plt = pyplot()
    groups = sorted({r['gcf'] for r in rows})
    palette = plt.get_cmap('turbo')
    colors = {g: palette(i/max(1,len(groups)-1)) for i,g in enumerate(groups)}
    fig, ax = plt.subplots(figsize=(9,7), layout='constrained')
    for group in groups:
        mask = np.array([r['gcf']==group for r in rows])
        ax.scatter(coordinates[mask,0], coordinates[mask,1], s=13, alpha=.65,
                   color=colors[group], linewidths=0, label=group, rasterized=len(rows)>5000)
    ax.set(title=title, xlabel=method.upper()+' 1', ylabel=method.upper()+' 2')
    if len(groups)<=20:
        ax.legend(loc='center left', bbox_to_anchor=(1, .5), fontsize=7, frameon=False)
    ax.text(0, -.13, f'{len(rows):,} sample-GCF records; {len(groups):,} GCFs. Coordinates are for display; scores use coverage space.',
            transform=ax.transAxes, fontsize=8)
    for suffix in ('pdf','png'):
        fig.savefig(folder/f'scatter.{suffix}', dpi=200, bbox_inches='tight')
    plt.close(fig)
    write_tsv(folder/'gcf_colors.tsv', ['gcf','rgba'],
              [dict(gcf=g,rgba=','.join(f'{v:.8g}' for v in colors[g])) for g in groups])


def curves(scores, output):
    plt = pyplot()
    for stratum in sorted({r['stratum'] for r in scores}):
        for metric in sorted({r['distance'] for r in scores if r['stratum']==stratum}):
            rows = sorted([r for r in scores if r['stratum']==stratum and r['distance']==metric],key=lambda r:r['bin_count'])
            good = [r for r in rows if r['ratio'] is not None]
            fig, ax = plt.subplots(figsize=(6.8,4.6), layout='constrained')
            if good:
                ax.plot([r['bin_count'] for r in good], [r['ratio'] for r in good], marker='o', color='#167F82')
                if max(r['bin_count'] for r in good)/min(r['bin_count'] for r in good)>10:
                    ax.set_xscale('log')
            else:
                ax.text(.5,.5,'No estimable between/within ratio',ha='center',transform=ax.transAxes)
            ax.set(xlabel='Bin Count',ylabel='Mean between / mean within distance',title=f'{stratum} | {metric}')
            ax.grid(alpha=.15)
            ax.text(0,-.17,'Scores in normalized coverage space; equal weight per GCF within, per center pair between.',
                    transform=ax.transAxes,fontsize=7)
            for suffix in ('pdf','png'):
                fig.savefig(output/f'{stratum}_{metric}_score_curve.{suffix}',dpi=200,bbox_inches='tight')
            plt.close(fig)
