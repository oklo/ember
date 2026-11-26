import math

def same_restart(source, restored):
    persisted=('years','accepted','central_T_K','central_X','surface_X','surface_Z','luminosity_Lsun')
    if any(source[k]!=restored[k] for k in persisted):return False
    # The boundary solver reconstructs the photosphere from the saved interior.
    # Check both the small surface shift and fixed L proportional to R^2 Teff^4.
    r=restored['radius_Rsun']/source['radius_Rsun']
    t=restored['Teff_K']/source['Teff_K']
    return max(abs(r-1),abs(t-1))<1e-6 and math.isclose(r*r*t**4,1.,rel_tol=1e-11)


def audited_steps(history, attempts, report):
    """Match each accepted solver attempt to the two history rows it connects."""
    accepted = [row for row in attempts if row['accepted']]
    if not accepted or len(accepted) != len(history)-1 or len(accepted) != report['accepted_this_invocation']:
        return False
    return all(trial['converged'] and trial['audit_pass']
               and trial['start_years'] == left['years']
               and trial['step_years'] == right['step_years']
               and right['accepted'] == left['accepted']+1
               and right['years'] > left['years']
               for trial,left,right in zip(accepted,history[:-1],history[1:],strict=True))
