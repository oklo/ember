"""Combine cluster catalogues, retaining each survey's own selection evidence.

This is an observed candidate union, not a completeness-corrected birth census.
No published evolutionary masses are used. Gaia DR2 IDs are not treated as DR3.
"""
from pathlib import Path
import json
import math
import re
import sqlite3
import csv
import io

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u

W = Path(__file__).resolve().parent
C = W / "catalogues"
objects, records, matches = {}, [], []
counter = {}


def number(value):
    try:
        f = float(value)
        return f if math.isfinite(f) else None
    except (ValueError, TypeError):
        return None


def clean(row):
    return {str(k): (v.item() if isinstance(v, np.generic) else v)
            for k, v in row.items() if not pd.isna(v)}


def read(name, ids=()):
    return pd.read_csv(C / name, dtype={key: str for key in ids})


def pos2016(ra, dec, epoch=2016, pmra=0, pmdec=0):
    # pmra includes cos(dec); first-order propagation is adequate over decades
    # at these proper motions. Gaia catalogue crosschecks establish epoch 2016.
    return (ra + (2016-epoch)*(pmra or 0)/3.6e6/math.cos(math.radians(dec)),
            dec + (2016-epoch)*(pmdec or 0)/3.6e6)


def add(cluster, source, rid, row, ra, dec, *, gaia=None, epoch=2016,
        pmra=0, pmdec=0, probability=None, selected=False, kind="cluster_or_tail",
        radius=0.7, values=None, priority=0):
    ra, dec = pos2016(ra, dec, epoch, pmra, pmdec)
    oid = f"{cluster}:gaia:{gaia}" if gaia else None
    sep = None
    if not oid:
        candidates = [o for o in objects.values() if o["cluster"] == cluster
                      and abs(o["dec_2016_deg"]-dec) < radius/3600]
        if candidates:
            sky = SkyCoord([o["ra_2016_deg"] for o in candidates]*u.deg,
                           [o["dec_2016_deg"] for o in candidates]*u.deg)
            distances = SkyCoord(ra*u.deg, dec*u.deg).separation(sky).arcsec
            order = np.argsort(distances)
            # An ambiguous match stays separate for review.
            if distances[order[0]] <= radius and (len(order) == 1 or distances[order[1]] > radius):
                target = candidates[order[0]]
                oid, sep = target["object_id"], float(distances[order[0]])
    if not oid:
        counter[cluster] = counter.get(cluster, 0)+1
        oid = f"{cluster}:source:{source}:{rid}"
    if oid not in objects:
        objects[oid] = dict(object_id=oid, cluster=cluster, gaia_dr3_id=gaia,
                            ra_2016_deg=ra, dec_2016_deg=dec,
                            selected_by_any=False, source_names=set(), kinds=set(),
                            values={}, priorities={})
    obj = objects[oid]
    obj["source_names"].add(source)
    obj["kinds"].add(kind)
    obj["selected_by_any"] |= bool(selected)
    for key, value in (values or {}).items():
        if value is not None and not pd.isna(value) and priority > obj["priorities"].get(key,-1):
            obj["values"][key] = value
            obj["values"][key+"_source"] = source
            obj["priorities"][key] = priority
    records.append(dict(object_id=oid, cluster=cluster, source=source,
                        source_id=str(rid), membership_probability=probability,
                        catalogue_grade=row.get('Grade'),original_catalogue=row.get('Ref'),
                        selected=int(selected), association=kind,
                        match_separation_arcsec=sep,
                        original=json.dumps(clean(row), allow_nan=False)))
    if sep is not None:
        matches.append(dict(object_id=oid,source=source,source_id=str(rid),separation_arcsec=sep))


# All Gaia records are merged by the exact DR3 source ID, regardless of whether
# they passed a catalogue's recommended membership threshold.
for cl in ["hyades", "pleiades"]:
    for _, s in read(f"{cl}_kos2024.csv", ["GaiaDR3"]).iterrows():
        row = clean(s); plx=number(row.get("plx")); d=1000/plx if plx and plx>0 else None
        g=number(row.get("Gmag")); bp=number(row.get("BPmag")); rp=number(row.get("RPmag"))
        ag=number(row.get("AG")); ebp=number(row.get("E(BP-RP)"))
        vals=dict(distance_pc=d,parallax_mas=plx,Gmag=g,
                  BP_RP=(bp-rp if bp is not None and rp is not None else None),
                  G_RP=(g-rp if g is not None and rp is not None else None),
                  MG_observed=(g-5*math.log10(d/10) if d and g is not None else None),
                  A_G=ag,E_BP_RP=ebp,RUWE=number(row.get("RUWE")))
        add(cl,"kos2024",row["GaiaDR3"],row,row["RA_ICRS"],row["DE_ICRS"],
            gaia=row["GaiaDR3"],probability=row["pbint"],selected=row["pbint"]>.5,
            values=vals,priority=10)
    for _, s in read(f"{cl}_liu2025_members.csv", ["GaiaDR3"]).iterrows():
        row=clean(s);d=row["rpgeo"];mg=row["GMag"]
        add(cl,"liu2025",row["GaiaDR3"],row,row["RAJ2000"],row["DEJ2000"],
            gaia=row["GaiaDR3"],probability=row["Prob"],selected=row["Prob"]>=.5,
            kind=("outside_tidal_radius" if row["Memb"]==2 else "inside_tidal_radius"),
            values=dict(distance_pc=d,parallax_mas=row["plx"],MG_observed=mg,
                        MG_dereddened=row["G0Mag"],BP_RP=row["BP-RP"],
                        BP_RP_dereddened=row["(BP-RP)0"],Gmag=mg+5*math.log10(d/10)),
            priority=30)

# Tail catalogues use different selection geometries. Retain all published
# candidates, including the longer tail; FCorrSamp is a comparison cut, not a
# membership verdict. Source classifications and flags remain queryable.
lines=[s for s in (W/'raw/risbud2025_members.tsv').read_text().splitlines()
       if s.strip() and not s.startswith('#')]
tail=pd.read_csv(io.StringIO('\n'.join([lines[0]]+lines[3:])),sep='\t',dtype=str)
tail=tail.map(lambda x:x.strip() if isinstance(x,str) else x)
for _,s in tail[tail.Cluster.isin(['Melotte_22','Melotte_25'])].iterrows():
    row=clean(s);cl='hyades' if row['Cluster']=='Melotte_25' else 'pleiades'
    plx=float(row['plx']);kind=('leading_tail' if row['FleadingTail']=='1' else
        'trailing_tail' if row['FtrailingTail']=='1' else 'inside_tidal_radius')
    add(cl,'risbud2025',row['GaiaDR3'],row,float(row['RA_ICRS']),float(row['DE_ICRS']),
        gaia=row['GaiaDR3'],selected=True,kind=kind,priority=12,
        values=dict(distance_pc=1000/plx if plx>0 else None,Gmag=number(row['Gmag']),
                    MG_observed=number(row['GMAG']),BP_RP=number(row['BP-RP'])))

# EDR3 source identifiers are shared with DR3. DR2 identifiers are different;
# merge those positions after proper-motion propagation, retaining unmatched
# DR2 candidates without falsely labelling their IDs as DR3.
for cl in ['hyades','pleiades']:
    lines=[s for s in (W/f'raw/{cl}_jadhav2025.tsv').read_text().splitlines()
           if s.strip() and not s.startswith('#')]
    data=pd.read_csv(io.StringIO('\n'.join([lines[0]]+lines[3:])),sep='\t',dtype=str)
    data=data.map(lambda x:x.strip() if isinstance(x,str) else x)
    for _,s in data.iterrows():
        row=clean(s);kind={'C':'inside_tidal_radius','L':'leading_tail','T':'trailing_tail'}[row['Class']]
        # A catalogue grade assesses the recovered tail as a whole. It is not
        # an individual membership probability and is not substituted for one.
        add(cl,'jadhav2025_'+row['Ref'],row['GaiaDR3'],row,float(row['RA_ICRS']),float(row['DE_ICRS']),
            gaia=row['GaiaDR3'],selected=True,kind=kind,priority=5,
            values=dict(distance_pc=number(row['rmedgeo']),Gmag=number(row['Gmag']),
                        MG_observed=number(row['GMAG']),BP_RP=number(row['BP-RP'])))

for release,filename,specs in [
    ('edr3','hyades_jerabkova2021_b2.dat',[(4,23),(24,44),(45,66),(67,88),(89,111)]),
    ('dr2','hyades_jerabkova2021_b1.dat',[(5,24),(25,45),(46,67),(74,95),(96,118)])]:
    rows=pd.read_fwf(W/'raw'/filename,colspecs=specs,names=['GaiaID','ra','dec','pmra','pmdec'],dtype={'GaiaID':str})
    assert len(rows)==(862 if release=='edr3' else 1109)
    for _,s in rows.iterrows():
        row=clean(s)
        add('hyades',f'jerabkova2021_{release}',row['GaiaID'],row,row['ra'],row['dec'],
            gaia=row['GaiaID'] if release=='edr3' else None,
            epoch=2016 if release=='edr3' else 2015.5,pmra=row['pmra'],pmdec=row['pmdec'],
            selected=True,kind='model_selected_extended_tail',radius=.7)

for _, s in read("pleiades_boyle2025.csv", ["GaiaID"]).iterrows():
    row=clean(s)
    add("pleiades","boyle2025",row["GaiaID"],row,row["RAdeg"],row["DEdeg"],
        gaia=row["GaiaID"],probability=row["Prob"],selected=row["Prob"]>=.5,
        kind="greater_pleiades_complex",values=dict(MG_dereddened=row["GMag"],
        BP_RP_dereddened=row["Bp-Rp"],binary_flag=row["FlagB"],cmd_consistent=row["FlagP"]),priority=20)

for i,s in read("pleiades_hobart2026.csv").iterrows():
    row=clean(s);plx=number(row.get("par"));d=1000/plx if plx and plx>0 else None
    add("pleiades","hobart2026",i,row,row["ra"],row["dec"],selected=True,
        kind="cluster_selection",radius=.3,priority=40,
        values=dict(distance_pc=d,parallax_mas=plx,e_parallax_mas=number(row.get("epar")),
                    Gmag=row.get("gmag"),MG_observed=row.get("G_ABS"),
                    BP_RP=(row["bmag"]-row["rmag"] if "bmag" in row and "rmag" in row else None),
                    G_RP=row.get("g_r"),Jmag=row.get("jmag"),Kmag=row.get("kmag")))

for _,s in read("pleiades_dance_resolved.csv", ["SarroID"]).iterrows():
    row=clean(s)
    add("pleiades","olivares2018_bouy2015",row["SarroID"],row,row["RAJ2000"],row["DEJ2000"],
        epoch=row["Epoch"],pmra=row["pmRA"],pmdec=row["pmDE"],
        probability=row["pc"],selected=row["pc"]>.84,kind="dance_footprint",priority=15,
        values=dict(Jmag=row.get("Jmag"),Kmag=row.get("Kmag"),
                    Lbol_Lsun_136p2pc=row.get("Lbol"),e_Lbol_Lsun_136p2pc=row.get("e_Lbol"),
                    equal_mass_binary_probability=row.get("pemb")))


def cwise_name(text):
    m=re.search(r'(J\d{6}\.\d+)\s*\++\s*\+?\s*(\d{6}\.\d+)',text)
    assert m,text
    return m[1]+"+"+m[2]


def cwise_position(name):
    m=re.fullmatch(r'J(\d\d)(\d\d)(\d\d\.\d+)([+-])(\d\d)(\d\d)(\d\d\.\d+)',name)
    h,mi,se,sign,de,dm,ds=m.groups()
    return 15*(int(h)+int(mi)/60+float(se)/3600), (1 if sign=='+' else -1)*(int(de)+int(dm)/60+float(ds)/3600)


def pair(text):
    nums=re.findall(r'[-+]?\d+(?:\.\d+)?',text)
    return [float(x) for x in nums[:2]]


# Schneider 2022 Table 6: original observed distances and membership evidence,
# excluding that paper's evolutionary mass estimates. Bracketed distances are
# photometric; all others come from parallaxes. Photometry is parsed separately
# from Tables 1 and 2. The source line is preserved for checking.
hy22 = """J031042.59+204629.3|L5|32||photometric|47.7
J031627.87+265027.2|T2.5|29.9|2.8|parallax|85.6
J032905.95+132231.5|T3.5|22.6|1.5|parallax|92.5
J033817.87+171744.1|L7|35||photometric|83.6
J035246.42+211232.7|L1.5|56.5|6.4|parallax|24.4
J035304.34+041820.0|L6pec|31||photometric|14.7
J035542.11+225700.9|L3|25.6|10.7|parallax|82.0
J041024.02+145910.1|L0.5|55.9|3.7|parallax|72.6
J041232.79+104408.0|L5pec|46||photometric|96.2
J041733.97+143015.2|L2|35.1|4.8|parallax|99.2
J041835.00+213126.6|L5|38.8|4.4|parallax|98.3
J041953.55+203628.0|T4|37||photometric|93.3
J042418.72+063745.5|L4|64||photometric|93.4
J042731.38+074344.9|L7|37||photometric|90.0
J042922.88+153529.4|T1|29.9|11.3|parallax|74.7
J043018.70+105857.1|T3|49||photometric|94.4
J043038.87+130956.7|T2.5|32.5|1.6|parallax|98.7
J043543.04+132344.8|L6pec|41.5|3.6|parallax|99.5
J043642.79+190134.6|L6|35||photometric|96.2
J043855.29+042300.6|T2|27.3|4.3|parallax|86.0
J043941.41+202514.8|T3|47||photometric|62.0
J044105.60+213001.3|L5pec|45||photometric|90.5
J044635.44+145125.7|L3.5|48.5|5.9|parallax|98.5
J045845.76+121234.1|L0.5|44.0|1.7|parallax|99.0
J053204.60+111955.1|L7|20||photometric|7.1"""
hytext=(W/'literature/hyades_schneider2022.txt').read_text()
phot22={}
for lo,hi,columns in [(hytext.index('Table 1. Recovered'),hytext.index('Table 2. New'),4),
                      (hytext.index('Table 2. New'),hytext.index('a Based on the UKIRT',hytext.index('Table 2. New')),9)]:
    for line in hytext[lo:hi].splitlines():
        m=re.match(r'\s*(J\d{6}\.\d+\+\d{6}\.\d+)[a-z]?\s+',line)
        if not m:continue
        # Missing bands remain explicit tokens, so J/K indices stay fixed.
        vals=re.findall(r'([-+]?\d+\.\d+)±(\d+\.\d+)|(\.\.\.)',line[m.end():])
        if len(vals)!=columns:continue
        jk=vals[2:4] if columns==4 else vals[5:7]
        phot22[m[1]]=dict(Jmag=float(jk[0][0]),e_Jmag=float(jk[0][1]),
                         Kmag=float(jk[1][0]),e_Kmag=float(jk[1][1]))
for line in hy22.splitlines():
    name,sp,d,ed,method,p=line.split('|');d=float(d);p=float(p)/100
    raw=dict(name=name,spectral_type=sp,distance_pc=d,e_distance_pc=number(ed),
             distance_method=method,membership_probability=p,**phot22.get(name,{}))
    source_line=next(l for l in hytext[hytext.index('Table 6. Substellar'):].splitlines() if l.startswith(name))
    proper=re.findall(r'([-+]?\d+\.\d+)±(\d+\.\d+)',source_line)[:2]
    raw['source_line']=source_line
    ra,dec=cwise_position(name)
    vals=dict(spectral_type=sp,**phot22.get(name,{}))
    vals['distance_pc' if method=='parallax' else 'photometric_distance_pc']=d
    if ed:vals['e_distance_pc']=float(ed)
    add('hyades','schneider2022',name,raw,ra,dec,epoch=2014,
        pmra=float(proper[0][0]),pmdec=float(proper[1][0]),probability=p,
        selected=p>.5,kind='substellar_candidate',radius=2,priority=25,values=vals)


# Small 2024 Hyades table: retain both distance-dependent and purely kinematic
# membership estimates. Include all eight published candidates; the paper's
# distance-conditioned selection includes one that has a lower probability
# from proper motion alone. Preserve both, rather than silently dropping it.
with (W/'raw/hyades_schneider2024_table1.csv').open() as f: phot=list(csv.reader(f))[2:]
with (W/'raw/hyades_schneider2024_table3.csv').open() as f: kin=list(csv.reader(f))[2:]
phot={cwise_name(r[0]):r for r in phot if r and r[0].startswith('J')}
for row in kin:
    if not row or not row[0].startswith('J'):continue
    name=cwise_name(row[0]);p=phot[name];j,ej=pair(p[5]);k,ek=pair(p[6]);dp,edp=pair(row[2])
    raw=dict(name=name,spectral_type=row[1],photometric_distance_pc=dp,
             e_photometric_distance_pc=edp,kinematic_distance_pc=float(row[3]),
             banyan_distance_pc=float(row[4]),p_kinematic=float(row[7])/100,
             p_with_photometric_distance=float(row[8])/100,Jmag=j,e_Jmag=ej,Kmag=k,e_Kmag=ek)
    ra,dec=cwise_position(name)
    add('hyades','schneider2024',name,raw,ra,dec,epoch=2014,
        pmra=pair(p[1])[0],pmdec=pair(p[2])[0],probability=raw['p_kinematic'],
        selected=True,kind='substellar_candidate',radius=2,
        priority=30,values=dict(Jmag=j,Kmag=k,photometric_distance_pc=dp,
            e_photometric_distance_pc=edp,kinematic_distance_pc=float(row[3]),spectral_type=row[1]))

# Calar 2014 unresolved candidates; the faintest confirmed objects additionally
# have measured SED luminosities in the 2018 follow-up. Other faint candidates
# retain their uncertainty instead of being promoted to confirmed detections.
text=(W/'literature/pleiades_deep2014.txt').read_text()
for line in text.splitlines():
    m=re.search(r'Calar Pleiades\s+(\d+)e?\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+\+(\d+)\s+(\d+)\s+(\d+\.\d+)',line)
    if not m:continue
    num,h,mi,se,d,dm,ds=m.groups();num=int(num)
    ra=15*(int(h)+int(mi)/60+float(se)/3600);dec=int(d)+int(dm)/60+float(ds)/3600
    rest=line[m.end():].replace('−','-')
    pairs=re.findall(r'([-+]?\d+\.\d+)\s*±\s*(\d+\.\d+)',rest)
    if len(pairs)<4:continue
    pmra,pmdec=[float(p[0]) for p in pairs[:2]]
    j,ej=map(float,pairs[2]);hv,eh=map(float,pairs[3]);kv,ek=(map(float,pairs[4]) if len(pairs)>4 else (None,None))
    # Only the top block (8--26) passes the original photometric selection;
    # 23 was rejected, and 26 is questioned by the 2018 follow-up.
    selected=(8<=num<=26 and num not in (23,26))
    status='confirmed_2018' if num in (21,22,25) else ('rejected_2018' if num==23 else 'candidate')
    raw=dict(calar=num,Jmag=j,e_Jmag=ej,Hmag=hv,Kmag=kv,e_Kmag=ek,
             followup_status=status,paper_line=line.strip())
    vals=dict(Jmag=j,Kmag=kv)
    if num in (21,22,25):
        vals.update(logL_SED_2018={21:-4.29,22:-4.33,25:-4.66}[num],e_logL_SED_2018=.10)
    add('pleiades','zapatero2014_2018',str(num),raw,ra,dec,epoch=2007,
        pmra=pmra,pmdec=pmdec,selected=selected,kind=status,radius=1.5,
        priority=50,values=vals)

# Save a queryable database with full source evidence and a compact object view.
flat=[]
for o in objects.values():
    d={k:v for k,v in o.items() if k not in ('values','priorities','source_names','kinds')}
    d.update(o['values']);d['source_names']=';'.join(sorted(o['source_names']))
    d['associations']=';'.join(sorted(o['kinds']))
    flat.append(d)
frame=pd.DataFrame(flat)
db=W/'cluster_catalogues.sqlite'
with sqlite3.connect(db) as con:
    frame.to_sql('objects',con,if_exists='replace',index=False)
    pd.DataFrame(records).to_sql('source_records',con,if_exists='replace',index=False)
    con.execute('CREATE UNIQUE INDEX objects_key ON objects(object_id)')
    con.execute('CREATE INDEX records_object ON source_records(object_id)')
    con.execute('CREATE INDEX objects_cluster ON objects(cluster)')
frame.to_csv(C/'combined_objects.csv',index=False)
pd.DataFrame(matches).to_csv(C/'coordinate_matches.csv',index=False)
summary={}
for cl in ('hyades','pleiades'):
    f=frame[frame.cluster==cl];s=f[f.selected_by_any]
    summary[cl]=dict(candidate_objects=len(f),selected_by_at_least_one_catalogue=len(s),
                    with_direct_SED_luminosity=int(s.get('Lbol_Lsun_136p2pc',pd.Series()).notna().sum()),
                    source_counts={source:sum(r['cluster']==cl and r['source']==source and r['selected'] for r in records)
                                   for source in sorted({r['source'] for r in records if r['cluster']==cl})})
(C/'membership_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
