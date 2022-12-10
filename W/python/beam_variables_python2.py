beam_variables = dict()
#kin_var --> {Energy: {Current: [Run Numbers]}}
beam_variables = {1.92 : {1.0 : 11207, 2.0 : 11214,3.0 : 11218, 4.0 : 11226, 5.0 : 11228, 6.0: ''}, 
	  		5.56 : {1.0 : '', 2.0 : '',3.0 : '', 4.0 : '', 5.0 : '', 6.0: ''},
			3.70 : {1.0 : 11420, 2.0 : 11421,3.0 : 11422, 4.0 : 11423, 5.0 : 11424, 6.0: 11425},

			}

##Dictionary to specify targets for specific run energies
targets = dict()
targets = {	"LH2": [1.92, 3.70, 5.6],
			"Targ": [7.0],
			}

##Dictionary to speficy the Invariant Mass min and max histogram bin values for each run energy
bin_variables_w = dict()
bin_variables_w = {1.92: [0.4, 1.8],
				3.70: [0.7, 2.7],
				5.56: [],

				}

##Dictionary to specify the Q-Squared min and max histogram bin values for each energy
bin_variables_q = dict()
bin_variables_q = {1.92: [0.4, 1.0],
				3.70: [1.0, 5.0],
				5.56: [],

				}	

qlist = list(bin_variables_q)							

##Dictionary to specify the energy cut min and max for each run energy
cut_energy = dict()
cut_energy = {1.92: [0.6, 0.8],
			3.70: [0.78, 1.15],
			5.56: [],

			}


def BeamEnergyFromRun(run):
	for energy, current in beam_variables.items():
		if (run in current.values()):
			return energy

def CurrentFromRun(run):
	for energy, current, in beam_variables.items():
		if (run in current.values()):
			return current.keys()[current.values().index(run)]

def TargetFromRun(run):
	for target, energy in targets.items():
		if (BeamEnergyFromRun(run) in energy):
			return target

def WBinFromRun(run):
	E = BeamEnergyFromRun(run)
	if (E in bin_variables_w.keys()):
		return bin_variables_w.values()[bin_variables_w.keys().index(E)]

def qBinFromRun(run):
	E = BeamEnergyFromRun(run)
	if (E in bin_variables_q.keys()):
		return bin_variables_q.values()[bin_variables_q.keys().index(E)]


def CutFromRun(run):
	E = BeamEnergyFromRun(run)
	if (E in cut_energy.keys()):
		return cut_energy.values()[cut_energy.keys().index(E)]


