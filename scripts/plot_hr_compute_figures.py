"""Draw separate HR/source-compute and hydrogen figures from frozen records."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import RegularPolygon
from matplotlib.ticker import FuncFormatter
import numpy as np


def read(path):
    return json.loads(Path(path).read_text())


def plot(directory, published):
    primary = [json.loads(l) for l in gzip.decompress((directory/'metal_history.jsonl.gz').read_bytes()).decode().splitlines()]
    fixed = [json.loads(l) for l in gzip.decompress((directory/'consistent_cn_history.jsonl.gz').read_bytes()).decode().splitlines()]
    pp_path = directory/'pp_comparison_history.csv'
    with pp_path.open() as stream:
        pp = [{k:float(v) for k,v in row.items() if v!=''} for row in csv.DictReader(stream)]
    data = read(directory/'source_compute.json')
    assert hashlib.sha256((directory/'metal_history.jsonl.gz').read_bytes()).hexdigest() == data['history_sha256']
    d = read(published)
    p = np.array(d['pixels']['hr']); a = d['calibration']['hr']
    lt = a['left_Teff_K']+(p[:,0]-a['left_px'])/(a['right_px']-a['left_px'])*(a['right_Teff_K']-a['left_Teff_K'])
    ll = a['top_log10_L_Lsun']+(p[:,1]-a['top_px'])/(a['bottom_px']-a['top_px'])*(a['bottom_log10_L_Lsun']-a['top_log10_L_Lsun'])
    p = np.array(d['pixels']['hydrogen']); a = d['calibration']['composition']
    la = (p[:,0]-a['zero_age_px'])/(a['age_6T_px']-a['zero_age_px'])*6
    lx = (a['zero_fraction_px']-p[:,1])/(a['zero_fraction_px']-a['unit_fraction_px'])
    c = lambda k: np.array([r[k] for r in primary])
    o = lambda k: np.array([r[k] for r in fixed])
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False, 'pdf.fonttype':42})
    orange, pale, gray = '#ad4d28', '#3f7887', '#777777'

    def style(ax):
        ax.grid(alpha=.12, zorder=0)
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.4g}'))

    xy, values = [], []
    channels = data['channels']
    for job in data['jobs']:
        hours = job['estimated_worker_seconds']/3600/len(job['points'])
        for point in job['points']:
            xy.append(point)
            row = [0., 0., 0.]
            row[channels.index(job['channel'])] = hours
            values.append(row)
    xy, values = np.asarray(xy), np.asarray(values)
    assert np.allclose(values.sum(axis=0), [data['channel_worker_hours'][c] for c in channels],
                       rtol=2e-12, atol=1e-10)
    extent = [2500, 6100, -3.2, -1.9]
    inside = ((xy[:,0]>=extent[0]) & (xy[:,0]<=extent[1]) &
              (xy[:,1]>=extent[2]) & (xy[:,1]<=extent[3]))
    # Clip explicitly: partial off-frame bins must not import off-frame cost.
    px, vv = xy[inside], values[inside]
    fig = plt.figure(figsize=(7.5,5.2))
    ax = fig.add_axes([.115,.14,.71,.70])
    sets = [ax.hexbin(px[:,0], px[:,1], C=vv[:,i], gridsize=(30,12),
                     extent=extent, reduce_C_function=np.sum, mincnt=0,
                     linewidths=.35, edgecolors='white', zorder=1) for i in range(3)]
    offsets = sets[0].get_offsets().copy()
    assert all(np.array_equal(s.get_offsets(), offsets) for s in sets)
    cost = np.column_stack([np.ma.filled(s.get_array(), 0) for s in sets])
    assert np.allclose(cost.sum(axis=0), vv.sum(axis=0), rtol=2e-12, atol=1e-10)
    total = cost.sum(axis=1)
    visible_worker_hours = float(total.sum())
    sequence = read(directory/'sequence_compute.json')
    assert sequence['sequence_states']['metal'] == len(primary)
    assert sequence['sequence_states']['fixed'] == len(fixed)
    flop_source = read(directory/'flop_estimate.json')
    rates = flop_source['adopted_flops_per_worker_second']
    operations = {key: sum(h*3600*rates[channel][key] for channel,h in
                          zip(channels,cost.sum(axis=0)))
                  for key in ['central','low','high']}
    flop_exponent = int(np.floor(np.log10(operations['central'])))
    flop_mantissa = operations['central']/10**flop_exponent
    # One common absolute scale: channel shares set mixed RGB color, and total
    # Estimated CPU-hours set opacity. Historical worker durations are the
    # source-time estimator, not a second unit on the displayed scale.
    basis = np.array([[.88,.15,.12], [.13,.63,.26], [.15,.32,.87]])
    maximum = float(total.max())
    maximum_alpha = .38
    alpha = lambda h: maximum_alpha*np.log1p(np.asarray(h)/.1)/np.log1p(maximum/.1)
    fractions = np.divide(cost, total[:,None], out=np.zeros_like(cost), where=total[:,None]>0)
    rgb = fractions @ basis
    rgba = np.column_stack([rgb, alpha(total)])
    for s in sets[1:]:
        s.remove()
    sets[0].set_array(None)
    sets[0].set_facecolors(rgba)
    ax.plot(lt,ll,color=pale,alpha=.38,lw=1.1,zorder=3,label='LBA97: published track')
    pp_line, = ax.plot([r['Teff_K'] for r in pp],np.log10([r['L_Lsun'] for r in pp]),
                      color='#656b75',alpha=.42,lw=1.1,ls='--',zorder=3,
                      label='pp comparison: no element separation')
    pp_line.set_path_effects([pe.Stroke(linewidth=2.4,foreground='white',alpha=.65),pe.Normal()])
    ax.plot(o('Teff'),np.log10(o('luminosity')/3.828e33),color=gray,alpha=.48,lw=1.1,zorder=3,label='Fixed-metal comparison')
    line, = ax.plot(c('Teff'),np.log10(c('luminosity')/3.828e33),color=orange,lw=1.8,zorder=5,label='Ember: settling metals')
    line.set_path_effects([pe.Stroke(linewidth=3.3,foreground='white',alpha=.9),pe.Normal()])
    ax.plot(c('Teff')[-1],np.log10(c('luminosity')[-1]/3.828e33),'o',color=orange,ms=4,zorder=6)
    ax.set(xlim=extent[:2][::-1],ylim=extent[2:],xlabel='Effective surface temperature (K)',ylabel=r'$\log_{10}(L/L_\odot)$')
    style(ax)
    fig.legend(*ax.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.48,.985),
               ncol=2,frameon=False,fontsize=9)
    legend = fig.add_axes([.845,.16,.15,.65]); legend.set_axis_off()
    for i, (label,color) in enumerate(zip(['Opacity','EOS','Atmosphere'],basis)):
        y = .96 - i*.12
        legend.plot(.07,y,marker='s',ms=8,color=color,transform=legend.transAxes)
        legend.text(.21,y,label,fontsize=8.5,va='center',transform=legend.transAxes)
    legend.text(0,.53,'Estimated\nCPU-hours\nper cell',fontsize=8.5,linespacing=1.35,transform=legend.transAxes)
    levels = [v for v in [.1,1,5] if v <= maximum]
    for i,h in enumerate(levels):
        y = .34-i*.13
        legend.plot(.08,y,marker='h',ms=14,color=basis.mean(axis=0),alpha=float(alpha(h)),transform=legend.transAxes)
        legend.text(.27,y,f'{h:g}',fontsize=9,va='center',transform=legend.transAxes)
    ax.text(.035,.06,
            rf'Tables: $\sim 10^{{{flop_exponent}}}$ FLOP'
            '\n'+rf"Evolution: $\sim {sequence['estimated_operations']['central']/1e14:.0f}\times10^{{14}}$ FLOP",
            fontsize=7.2,color='black',ha='left',va='bottom',
            transform=ax.transAxes,zorder=7)
    fig.savefig(directory/'consistent_cn_comparison.pdf',metadata={'CreationDate':None})
    fig.savefig(directory/'hr_compute.png',dpi=180)
    plt.close(fig)
    # Both comparison panels use the same current primary history as Fig. 1.
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.65))
    fig.subplots_adjust(left=.105,right=.985,bottom=.17,top=.76,wspace=.34)
    for comparison in axes:
        comparison.plot(lt,ll,color=pale,alpha=.35,lw=1.1,label='LBA97: published track')
        comparison.plot([r['Teff_K'] for r in pp],np.log10([r['L_Lsun'] for r in pp]),
                        color=gray,alpha=.35,lw=1,ls='--',label='pp burning; no microscopic separation')
        comparison.plot(o('Teff'),np.log10(o('luminosity')/3.828e33),color=gray,
                        alpha=.45,lw=1.1,label='CN burning; fixed-metal material')
        comparison.plot(c('Teff'),np.log10(c('luminosity')/3.828e33),color=orange,
                        lw=1.8,label='Ember: settling metals')
        comparison.plot(c('Teff')[-1],np.log10(c('luminosity')[-1]/3.828e33),
                        'o',color=orange,ms=4)
        comparison.set(xlabel='Effective surface temperature (K)',ylabel=r'$\log_{10}(L/L_\odot)$')
        style(comparison)
    axes[0].set(xlim=(6100,1500),ylim=(-5.4,-1.9))
    axes[1].set(xlim=(5900,2500),ylim=(-3.2,-1.9))
    fig.legend(*axes[0].get_legend_handles_labels(),loc='upper center',ncol=2,
               frameon=False,fontsize=7.5)
    fig.savefig(directory/'lba97_hr.pdf',metadata={'CreationDate':None})
    plt.close(fig)

    fig,ax = plt.subplots(figsize=(7.5,4.1))
    fig.subplots_adjust(left=.115,right=.97,bottom=.17,top=.75)
    ax.plot(la,lx,color=pale,alpha=.38,lw=1.1)
    for key,ls in [('central_H','-'),('surface_H','--')]:
        ax.plot(o('age_years')/1e12,o(key),color=gray,alpha=.4,lw=1.1,ls=ls)
    ax.plot(c('age_years')/1e12,c('central_H'),color=orange,lw=1.8)
    ax.plot(c('age_years')/1e12,c('surface_H'),color=orange,lw=1.5,ls='--')
    # The arrow follows actual points on the rising dashed surface trajectory.
    i0 = int(np.argmin(abs(c('surface_H')-.73)))
    i1 = int(np.argmin(abs(c('surface_H')-.84)))
    assert i1 > i0 and c('surface_H')[i1] > c('surface_H')[i0]
    arrow = [[float(c('age_years')[i]/1e12),float(c('surface_H')[i])] for i in (i0,i1)]
    ax.annotate('',xy=arrow[1],xytext=arrow[0],arrowprops=dict(arrowstyle='-|>',color=orange,lw=1.5,
                mutation_scale=13,shrinkA=0,shrinkB=0),zorder=5)
    ax.set(xlim=(0,6.4),ylim=(0,1),xlabel='Age (Tyr)',ylabel='Hydrogen mass fraction')
    style(ax)
    handles = [Line2D([],[],color=orange,lw=1.8,label='Ember center'),
               Line2D([],[],color=orange,lw=1.5,ls='--',label='Ember surface'),
               Line2D([],[],color=gray,alpha=.4,lw=1.1,label='Fixed-metal comparison'),
               Line2D([],[],color=pale,alpha=.38,lw=1.1,label='LBA97: published curve')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.54,.97),ncol=2,frameon=False,fontsize=9)
    fig.savefig(directory/'hydrogen_evolution.pdf',metadata={'CreationDate':None})
    fig.savefig(directory/'hydrogen_evolution.png',dpi=180)
    plt.close(fig)
    result = dict(history_states=len(primary),endpoint=primary[-1],
                  displayed_time_unit='CPU-hour',
                  source_cpu_hours_estimated=cost.sum(axis=0).tolist(),
                  source_time_method='Serial worker durations approximate CPU time; EOS batch-worker durations are upper estimates. Waiting and unused batch capacity are not fully separated in historical records.',
                  sequence_time_method='Measured native process user plus system CPU time, shared histories counted once',
                  channels=channels,extent=extent,gridsize=[30,12],
                  bins=[dict(teff_K=float(p[0]),log_L=float(p[1]),worker_hours=v.tolist())
                        for p,v in zip(offsets,cost) if v.sum()>0],
                  visible_worker_hours=cost.sum(axis=0).tolist(),
                  compute_annotation=dict(
                      scope='Separate totals for source work in the displayed cells and the retained plotted model sequences',
                      time_unit='CPU-hour',
                      estimated_source_cpu_hours=visible_worker_hours,
                      estimated_worker_hours=visible_worker_hours,
                      measured_operation_count=False,
                      note_position_axes=[.035,.06],note_color='black',
                      basis='Order-of-magnitude estimate from retained worker times and workload throughput assumptions informed by measured atmosphere kernels',
                      estimated_operations=operations,
                      adopted_flops_per_worker_second=rates,
                      uncertainty='Approximately a factor of three; not a statistical confidence interval',
                      estimate_source_sha256=hashlib.sha256((directory/'flop_estimate.json').read_bytes()).hexdigest(),
                      sequence_estimated_operations=sequence['estimated_operations'],
                      sequence_cpu_hours=sequence['total_cpu_hours'],
                      sequence_compute_sha256=hashlib.sha256((directory/'sequence_compute.json').read_bytes()).hexdigest()),
                  allocated_worker_hours=values.sum(axis=0).tolist(),
                  color_basis=basis.tolist(),maximum_hexagon_alpha=maximum_alpha,
                  intensity='0.38*log1p(hours/0.1)/log1p(maximum_cell_hours/0.1)',
                  maximum_cell_hours=maximum,hydrogen_arrow=arrow,
                  pp_comparison_sha256=hashlib.sha256(pp_path.read_bytes()).hexdigest(),
                  source_compute_sha256=hashlib.sha256((directory/'source_compute.json').read_bytes()).hexdigest(),
                  plot_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (directory/'compute_hexagons.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(cells=len(cost),visible_hours=result['visible_worker_hours'],maximum_cell_hours=maximum),indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('published_digitization',type=Path)
    a = parser.parse_args()
    plot(a.directory,a.published_digitization)
