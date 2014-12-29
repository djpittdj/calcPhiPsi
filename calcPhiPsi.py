"""
Author: Jian Dai
"""
#!/usr/bin/python
from FSUNMRBase import NMRBase
from FSUPeptidePlane import Propagator
import re, math

class CalcPhiPsi:
	def __init__(self, s11, s22, s33, nuParallel, betaNH, infilename, outfilename):
		self.nmr = NMRBase()
		self.nmr.betaRad = self.nmr.toRad(betaNH)
		self.nmr.sigma11 = s11
		self.nmr.sigma22 = s22
		self.nmr.sigma33 = s33
		self.nmr.nuParallel = nuParallel
		self.propagator = Propagator(self.nmr)
		self.spectrum = []
		self.PhiPsi = []
		self.outfilename = outfilename

		# read in PISEMA data file
		infile = open(infilename , 'r')
		for line in infile.readlines():
			resid = int(line.split()[0])
			cs = float(line.split()[1])
			dipole = float(line.split()[2])
			# compute B in PAF for each observed data point
			sgn = 0
			if cs == 0.0 and dipole == 0.0 :
				ipos = []
				ipos.append(Vector(0.0, 0.0, 0.0))
			else:
				sgn = self.propagator.nmr.determineDipoleSign(dipole, cs, 2)
				[ipos, eps] = self.propagator.nmr.getBinPAF(cs, dipole, sgn, 1)
			data = [cs, dipole, ipos, sgn, eps, resid]
			self.spectrum.append(data)
		infile.close()

	def calcPhiPsiDeviation(self, phi, psi):
		ideal_phi = -60.0
		ideal_psi = -45.0
		dphi = math.fabs(phi - ideal_phi)
		dpsi = math.fabs(psi - ideal_psi)
		if dphi > 180.:
			dphi = 360. - dphi
		elif dphi < -180.:
			dphi = 360. + dphi
		if dpsi > 180.:
			dpsi = 360. - dpsi
		elif dpsi < -180.:
			dpsi = 360. + dpsi
		delta = math.sqrt(dphi**2.0 + dpsi**2.0)
		return delta

	def calcClosestPhiPsi(self, i, j):
		ipos = self.spectrum[i][2]
		jpos = self.spectrum[j][2]

		Tors = []
		Exact = []
		minDelta = 999.
		closestPhi = 0.0
		closestPsi = 0.0
		if (len(ipos) == 0 or len(jpos) == 0) :
			print "invalid sequence %i -> %i.\n" % (self.spectrum[i][-1], self.spectrum[j][-1])
		elif ipos[0].length() > 0. and jpos[0].length() > 0.:
			for di in range(0, len(ipos), 1):
				for dj in range(0, len(jpos), 1):
					[tors, exact, grams, eps] = self.propagator.computeTorsionsByChiralities(ipos[di], jpos[dj], self.nmr.betaRad, 0, 1)
					Tors.extend(tors)
					Exact.extend(exact)

			for k in range(0, len(Tors)):
				if Exact[k]:
					[phi, psi] = Tors[k]
					delta = self.calcPhiPsiDeviation(phi, psi)
					if delta < minDelta:
						minDelta = delta
						closestPhi = phi
						closestPsi = psi

		return [closestPhi, closestPsi]

	def calcTors(self):
		num_aa = len(self.spectrum)
		outfile = open(self.outfilename, 'w')
		for i in range(0, num_aa-1):
			j = i+1
			resid = self.spectrum[i][-1]
			[closestPhi, closestPsi] = self.calcClosestPhiPsi(i, j)
			(self.PhiPsi).append([closestPhi, closestPsi])
			print>>outfile, '%4i%4i%10.3f%10.3f' % (resid, resid+1, closestPhi, closestPsi)
		outfile.close()

sigma11 = 57.3
sigma22 = 81.2
sigma33 = 227.8
nuParallel = 10.735
betaNH = 17.0
dataFile = 'kdpf_pisema.dat'
calcPhiPsi = CalcPhiPsi(sigma11, sigma22, sigma33, nuParallel, betaNH, dataFile, 'phipsi.dat')
calcPhiPsi.calcTors()
