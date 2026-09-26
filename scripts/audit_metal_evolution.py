"""Independent mass, nuclear-energy and first-law audits for moving metals.

The seven physical inventories are H1, He3, He4, C12, C13, N14 and inert
metals. The inert group has no nuclear source; transport conserves its mass.
"""
import numpy as np
from audit_cn_thermal import ATOMIC, LIGHT, MASS, YEAR, weights


def physical(model):
    cn = model[:, 7:10]*MASS[3:]
    return np.column_stack((model[:, 5:7], 1-model[:, 5:7].sum(axis=1)-model[:, 10],
                            cn, model[:, 10]-cn.sum(axis=1)))


def audit(old, result, age, dt):
    if not result.get('converged'):
        raise ValueError(result.get('message', result.get('error', 'missing native result')))
    new = np.array(result['model']); source = np.array(result['physical_sources'])
    energy = np.array(result['energy_cells']); w = weights(old[:, 0]); mass = old[-1, 0]
    before, after = physical(old), physical(new)
    if (new.shape != old.shape or source.shape != after.shape or energy.shape != (len(old), 6)
            or not np.isfinite(new).all() or not np.isfinite(source).all()
            or not np.isfinite(energy).all() or after.min() < 0
            or not np.array_equal(new[:, 0], old[:, 0])):
        raise ValueError('invalid physical state, source, energy or mesh')
    if not result['input_preserved'] or result['age_seconds'] != age+dt:
        raise ValueError('input state or absolute age changed')
    if abs(after.sum(axis=1)-1).max() > 2e-14:
        raise ValueError('physical mass fractions do not close')
    if np.any(source[:, -1] != 0):
        raise ValueError('inert metals have a nuclear source')
    regions = result['mixing_regions']
    if (not regions or regions[0][0] != 0 or regions[-1][1] != len(old)
            or any(b <= a for a, b in regions)
            or any(left[1] != right[0] for left, right in zip(regions, regions[1:]))):
        raise ValueError('mixing regions do not partition the star')
    delta = after-before; residual = delta-dt*source
    flux = np.zeros((len(new)-1, 7)); seen = set()
    for face, f in result['metal_boundary_fluxes']:
        if face in seen or not 0 <= face < len(new)-1 or len(f) != 6:
            raise ValueError('invalid or duplicate species boundary flux')
        seen.add(face)
        flux[face] = dt/mass*np.array([f[0], f[1], -sum(f), *f[2:]])
    if seen != {end-1 for _, end in regions[:-1]} or not np.isfinite(flux).all():
        raise ValueError('missing or invalid species boundary flux')
    region_error = 0.
    for begin, end in regions:
        r = w[begin:end]@residual[begin:end]/mass
        if begin: r -= flux[begin-1]
        if end < len(new): r += flux[end-1]
        region_error = max(region_error, float(abs(r).max()))
        if end > begin+1 and not np.all(new[begin:end, 5:11] == new[begin, 5:11]):
            raise ValueError('convective reservoir is not homogeneous')
    catalyst_before = float(w@old[:, 7:10].sum(axis=1))
    if catalyst_before <= 0:
        raise ValueError('this CN calculation requires a positive catalyst inventory')
    catalyst_error = float(abs(w@(new[:, 7:10].sum(axis=1)-old[:, 7:10].sum(axis=1)))/catalyst_before)
    # A fixed inert-group binding coefficient cancels exactly in the global
    # ledger. Check its conservation separately; do not count settling as
    # nuclear energy. Atomic masses include the neutralizing electrons.
    release = -float(w@delta[:, :6]@(ATOMIC/MASS-1))*LIGHT**2/dt
    nuclear = float(w@(energy[:, 3]+energy[:, 4])); photon = float(new[-1, 4])
    if not np.isfinite([release, nuclear, photon]).all() or nuclear < 0 or photon <= 0:
        raise ValueError('invalid power in nuclear energy audit')
    grav = -(energy[:, 0]-energy[:, 2]+energy[:, 1]*(1/new[:, 2]-1/old[:, 2]))/dt
    report = dict(region_species_error=region_error,
                  global_species_error=float(abs(w@residual/mass).max()),
                  catalyst_number_error=catalyst_error,
                  inert_mass_error=float(abs(w@delta[:, -1]/mass)),
                  mass_energy_error=(release-nuclear)/max(nuclear, photon),
                  first_law_error=float(w@(energy[:, 3]-energy[:, 5]+grav))/photon-1)
    limits = dict(region_species_error=2e-11, global_species_error=2e-11,
                  catalyst_number_error=2e-9, inert_mass_error=1e-13,
                  mass_energy_error=2e-6, first_law_error=2e-8)
    if any(not np.isfinite(v) or abs(v) > limits[k] for k, v in report.items()):
        raise ValueError('independent conservation check: '+str(report))
    return new, report, w[:, None]/mass*residual


def abundance_step_change(before, after):
    # Includes He4 and the inert group, so new metal coordinates cannot evade
    # the step proposal. The native abundance guard remains independent.
    return float(abs(physical(after)-physical(before)).max())


def compare(coarse, fine, composition_mode='pointwise', diagnostics=None,
            previous=None, intermediate=None, composition_relative_tolerance=1e-4,
            structure_relative_tolerance=1e-4, helium3_relative_tolerance=None):
    if not np.isfinite(composition_relative_tolerance) or not 0<composition_relative_tolerance<=.01:
        raise ValueError('composition relative tolerance must be positive and at most 0.01')
    if not np.isfinite(structure_relative_tolerance) or not 0<structure_relative_tolerance<=.01:
        raise ValueError('structure relative tolerance must be positive and at most 0.01')
    if helium3_relative_tolerance is not None and (
            not np.isfinite(helium3_relative_tolerance) or not 0<helium3_relative_tolerance<=.01):
        raise ValueError('helium-3 relative tolerance must be positive and at most 0.01')
    a, b = np.array(coarse['model']), np.array(fine['model'])
    pa, pb = physical(a), physical(b); w = weights(a[:, 0])
    mixed_a = [(x, y) for x, y in coarse['mixing_regions'] if y>x+1]
    mixed_b = [(x, y) for x, y in fine['mixing_regions'] if y>x+1]
    boundary = max((abs(x-y) for ra, rb in zip(mixed_a, mixed_b) for x, y in zip(ra, rb)), default=0)
    if len(mixed_a)!=len(mixed_b): boundary = len(a)
    species_error = abs(pa-pb)/(1e-8+composition_relative_tolerance*np.maximum(pa, pb))
    # A separately checked isotope budget changes only this local time-error
    # estimate. Global composition, fuel and conservation criteria retain
    # their original scales.
    if helium3_relative_tolerance is not None:
        species_error[:, 1] = abs(pa[:, 1]-pb[:, 1])/(1e-8+
            helium3_relative_tolerance*np.maximum(pa[:, 1],pb[:, 1]))
    distribution_error = float(w@abs(a[:, 5]-b[:, 5])/(w@b[:, 5])/1e-6)
    metrics = dict(structure=float(abs(np.log(a[:, 1:4]/b[:, 1:4])).max()/structure_relative_tolerance),
                   species=float(species_error.max()),
                   luminosity=float(abs(a[-1, 4]/b[-1, 4]-1)/1e-3),
                   hydrogen=distribution_error,
                   mixed_boundary=float(boundary))
    if composition_mode=='moving-boundary':
        if previous is None or intermediate is None:
            raise ValueError('moving-boundary comparison requires start and intermediate reservoirs')
        keep = np.ones(len(a), dtype=bool)
        paths = [[(x,y) for x,y in record['mixing_regions'] if y>x+1]
                 for record in [previous, intermediate, coarse, fine]]
        motion = 0
        topology_fallback = len({len(path) for path in paths})!=1
        if not topology_fallback:
            for reservoirs in zip(*paths):
                for boundaries in zip(*reservoirs):
                    motion = max(motion,max(boundaries)-min(boundaries))
            if motion<=1:
                for reservoirs in zip(*paths):
                    for boundaries in zip(*reservoirs):
                        keep[min(boundaries):max(boundaries)] = False
        # Endpoints may have the same boundary but cross the cell at different
        # times. Account for the cell traversed over the whole interval,
        # permitting at most one cell per boundary, not an arbitrary band.
        if diagnostics is not None:
            diagnostics.update(boundary_cells=np.flatnonzero(~keep).tolist(),
                               all_node_species_error_norm=metrics['species'],
                               hydrogen_distribution_error_norm=distribution_error,
                               topology_fallback=topology_fallback)
        metrics['species'] = float(species_error[keep].max())
        metrics['boundary_motion'] = float(motion)
        metrics['composition_l1'] = float(((w@abs(pa-pb))/
            (w.sum()*1e-8+composition_relative_tolerance*(w@np.maximum(pa,pb)))).max())
        # Fuel accuracy is the difference in total hydrogen. Redistribution
        # is constrained separately by composition_l1 and the local checks.
        metrics['hydrogen'] = (distribution_error if topology_fallback else
                              float(abs(w@(a[:,5]-b[:,5]))/(w@b[:,5])/1e-6))
    elif composition_mode!='pointwise':
        raise ValueError('unknown composition comparison')
    if diagnostics is not None:
        active_error = species_error.copy()
        if composition_mode=='moving-boundary': active_error[~keep] = -1
        node, species = np.unravel_index(active_error.argmax(), active_error.shape)
        diagnostics.update(maximum_species_error_node=int(node),
                           maximum_species_error_species=['H1','He3','He4','C12','C13','N14','inert'][species],
                           maximum_species_error_values=[float(pa[node,species]), float(pb[node,species])])
    return max(metrics.values()), metrics
