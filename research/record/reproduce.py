#!/usr/bin/env python3
"""Verify the retained record, redraw PNGs, or print a checked replay command.

This program never launches stellar evolution or source-table calculations.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shlex
import tempfile

HERE = Path(__file__).resolve().parent


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def verify():
    manifest = json.loads((HERE/'manifest.json').read_text())
    for entry in manifest['files']:
        path = HERE/entry['path']
        if digest(path) != entry['sha256']:
            raise ValueError(f'research record checksum mismatch: {path}')
        content = gzip.decompress(path.read_bytes()) if entry['encoding']=='gzip' else path.read_bytes()
        if hashlib.sha256(content).hexdigest() != entry['original_sha256']:
            raise ValueError(f'original content mismatch: {path}')
    return len(manifest['files'])


def read_tracks(study):
    with gzip.open(HERE/study/'tracks.json.gz', 'rt') as stream:
        return json.load(stream)['tracks']


def figures(output):
    verify()
    os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'ember-matplotlib'))
    os.environ.setdefault('XDG_CACHE_HOME', str(Path(tempfile.gettempdir())/'ember-font-cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LogNorm
    from matplotlib.ticker import FixedLocator, FormatStrFormatter, NullFormatter
    import numpy as np

    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                        'axes.spines.top': False, 'axes.spines.right': False,
                        'savefig.facecolor': 'white'})
    rendered = []

    def save(fig, name):
        path = output/name
        fig.savefig(path, dpi=180, bbox_inches='tight')
        plt.close(fig)
        rendered.append({'path': str(path), 'sha256': digest(path)})

    tracks = read_tracks('exploratory')
    fig, ax = plt.subplots(figsize=(8,6), layout='constrained')
    norm = LogNorm(.005,.235)
    for track in sorted(tracks, key=lambda t:t['mass'], reverse=True):
        rows, mass = track['rows'], track['mass']
        ax.plot([r['Teff_K'] for r in rows], [r['luminosity_Lsun'] for r in rows],
                color=plt.cm.viridis(norm(mass)), linewidth=2.3 if mass==.1 else 1.1)
    ax.set(xscale='log', yscale='log', xlabel=r'$T_{\rm eff}$ (K)', ylabel=r'$L/L_\odot$')
    ax.invert_xaxis()
    ax.xaxis.set_major_locator(FixedLocator([4000,2500,1500,800,400,200]))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_formatter(NullFormatter())
    fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='viridis'), ax=ax,
                 label=r'$M/M_\odot$', shrink=.7)
    save(fig,'exploratory-hr.png')

    tracks = read_tracks('boundary')
    fig, axes = plt.subplots(2,3,figsize=(12,7),sharex=True,layout='constrained')
    groups = [('lowz','Primordial / trace Z'),('z001',r'$Z=0.001$'),
              ('z005',r'$Z=0.005$'),('half',r'$0.5\,Z_\odot$'),
              ('solar',r'$Z_\odot$'),('enriched',r'$2$–$8\,Z_\odot$')]
    for ax,(group,title) in zip(axes.flat,groups):
        selected = [t for t in tracks if t['group']==group]
        for i,track in enumerate(selected):
            rows = [r for r in track['rows'] if r['years']>0]
            ax.plot([r['years']/1e9 for r in rows], [r['Teff_K'] for r in rows],
                    color=plt.cm.viridis(.1+.8*i/max(1,len(selected)-1)), linewidth=1.2)
        ax.set(xscale='log',title=title,xlabel='Elapsed time (Gyr)',ylabel=r'$T_{\rm eff}$ (K)')
    save(fig,'boundary-temperature.png')

    tracks = read_tracks('low-metallicity')
    fig,axes = plt.subplots(1,2,figsize=(10,4),sharex=True,layout='constrained')
    norm = LogNorm(1e-5,2e-3)
    for ax,group,title in zip(axes,['primordial','z005'],['Pop III',r'$Z=0.005$']):
        for track in tracks:
            if track['group'] != group:
                continue
            values=np.asarray(track['data'])
            values=values[values[:,0]>0]
            xy=values[:,:2]
            line=LineCollection(np.stack([xy[:-1],xy[1:]],axis=1),cmap='magma',norm=norm,
                                linewidths=2 if track.get('new') else 1.1)
            line.set_array((values[:-1,2]+values[1:,2])/2)
            ax.add_collection(line,autolim=False)
            ax.update_datalim(xy[(xy[:,0]>=10)&(xy[:,0]<=400)])
        ax.autoscale_view()
        ax.set(xscale='log',xlim=(10,400),title=title,xlabel='Elapsed time (Gyr)',ylabel=r'$T_{\rm eff}$ (K)')
    fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='magma'),ax=axes,
                 label=r'$L/L_\odot$',shrink=.8)
    save(fig,'low-metallicity.png')
    (output/'receipt.json').write_text(json.dumps({'figures':rendered,
        'new_stellar_runs':0,'new_source_calculations':0},indent=2)+'\n')
    print(f'Rendered {len(rendered)} PNGs in {output}')


def replay_command(launch_path, output):
    launch=json.loads(launch_path.read_text())
    command=list(launch['command'])
    if command[1] != '--lifetime':
        raise ValueError('this helper supports recorded --lifetime invocations')
    for path,expected in [(command[0],launch['executable_sha256']),
                          (command[2],launch['configuration_sha256'])]:
        if digest(Path(path)) != expected:
            raise ValueError(f'replay identity mismatch: {path}')
    if '--restart' in command:
        path=command[command.index('--restart')+1]
        if digest(Path(path)) != launch['restart_sha256']:
            raise ValueError(f'restart identity mismatch: {path}')
    if output.exists():
        raise ValueError('replay needs a new output directory')
    command[4]=str(output.resolve())
    print(shlex.join(command))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='action',required=True)
    commands.add_parser('verify')
    plot=commands.add_parser('figures')
    plot.add_argument('--output',type=Path,default=HERE/'figures')
    replay=commands.add_parser('command')
    replay.add_argument('launch',type=Path)
    replay.add_argument('output',type=Path)
    args=parser.parse_args()
    if args.action=='verify':
        print(f'Verified {verify()} retained files against their original contents.')
    elif args.action=='figures':
        figures(args.output)
    else:
        replay_command(args.launch,args.output)


if __name__=='__main__':
    main()
