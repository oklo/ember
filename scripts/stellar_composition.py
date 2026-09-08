"""Fixed interior mixture shared by the atmosphere and EOS data importers."""


def interior_composition():
    # Same documented AAG21 CNO approximation as C++ solar_scaled(.7,.02).
    c, n, o = (10**(8.46-12)*12.011, 10**(7.83-12)*14.007, 10**(8.69-12)*15.999)
    total = c+n+o
    c12 = .0144*c/total/(1+1/89)
    return [.7, 0, .28, c12, c12/89, .0144*n/total, .0144*o/total, .0056]
