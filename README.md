calcPhiPsi
==========
Reference: PIPATH: An optimized algorithm for generating alpha-helical structures from PISEMA data, Asbury et al.
Calculate the phi/psi angles from PISEMA data, utilizing code from PIPATH.
Consider only the case where exact solutions can be found, i.e., all the gramian terms should be positive.
Pick the phi/psi angle combination that is closest to the value of phi/psi for an ideal helix, i.e., -60/-45.
