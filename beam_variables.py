beam_variables = dict()
beam_variables = {1.92 : {1.0 : [11207,11420], 2.0 : [11214, 11421] ,3.0 : [11218, 11422], 4.0 : [11226, 11423], 5.0 : [11228, 11425], 6.0: ''}, 
				3.7 : {1.75: [11594, 11595], 3.75: [11600], 7.0: [11605], 10.0: [11608]},
				7.906: {1.0: [11989], 2.0: [11990], 4.0: [11991], 5.0: [12047, 12057], 8.0: [11993], 10.0: [12049, 12052, 12058]},
				9.91: {2.0: [12423], 4.0: [12342, 12424], 8.0: [12316, 12340, 12373, 12374], 15.0: [12367]},
				}

##Dictionary to specify targets for specific run energies
targets = dict()
targets = {"LH2": [11207, 11214, 11218, 11226, 11228, 11420, 11421, 11422, 11423, 11425, 11432, 11989, 11990, 11991, 11993, 12052, 12367, 12340],
			"LD2": [11594, 11595, 11600, 11605, 11608, 12047, 12049, 12057, 12058, 12316, 12342, 12423, 12424, 12373, 12374],
			}

##Dictionary to speficy the Invariant Mass min and max histogram bin values for each run energy
bin_variables_w = dict()
bin_variables_w = {1.92: [0.4, 1.8],
				3.70: [0.7, 2.7],
				7.906: [1.0, 3.5],
				9.91: [0.7, 2.5],
				}

##Dictionary to specify the Q-Squared min and max histogram bin values for each energy
bin_variables_q = dict()
bin_variables_q = {1.92: [0.4, 1.0],
				3.70: [1.0, 5.0],
				7.906: [3.0, 12],
				9.91: [5.0, 14],
				}	
						

##Dictionary to specify the energy cut min and max for each run energy
cut_energy = dict()
cut_energy = {1.92: [0.6, 0.8],
			3.70: [0.78, 1.15],
			7.906: [0.9, 1.1],
			9.91: [0.7, 1.3],
			}

elastic_peak = dict()
elastic_peak = {1.92: [0.9, 1.05], 
			3.70: [1.4, 1.55],
			7.906: [1.1, 1.5],
			9.91: [1.15, 1.5],
			}


def BeamEnergyFromRun(run):
	for energy, current in beam_variables.items():
		for runs in current.values():
			for vals in runs:
				if(vals == run):
					return energy
	return "RUN NOT FOUND IN DATABASE"

def CurrentFromRun(run):
	for energy, current in beam_variables.items():
		for beam_current, runs in current.items():
			for this_run in runs:
				if(this_run == run):
					return beam_current
	return "####ERROR: RUN NOT FOUND IN DATABASE####"

# def TargetFromRun(run):
# 	for target, energy in targets.items():
# 		if (BeamEnergyFromRun(run) in energy):
# 			return target

def TargetFromRun(run):
	for targ, runs in targets.items():
		for vals in runs:
			if(vals == run):
				return targ
	

def WBinFromRun(run):
	E = BeamEnergyFromRun(run)
	if (E in bin_variables_w.keys()):
		bin_index = list(bin_variables_w.keys()).index(E)
		return list(bin_variables_w.values())[bin_index]
		

def qBinFromRun(run):
	E = BeamEnergyFromRun(run)
	if (E in bin_variables_q.keys()):
		bin_index = list(bin_variables_q.keys()).index(E)
		return list(bin_variables_q.values())[bin_index]


def CutFromRun(run):
	E = BeamEnergyFromRun(run)
	if (E in cut_energy.keys()):
		cut_index = list(cut_energy.keys()).index(E)
		return list(cut_energy.values())[cut_index]

def PeakMinMaxFromRun(run):
	E = BeamEnergyFromRun(run)
	if (E in elastic_peak.keys()):
		peak_index = list(elastic_peak.keys()).index(E)
		return list(elastic_peak.values())[peak_index]

