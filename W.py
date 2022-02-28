from ROOT import *
from beam_variables import *
import numpy as np

def getXfromBin(bin_val):
	return ((max_bin_w - min_bin_w)/bins)*(bin_val) + min_bin_w

def plot_W(runnum, rootFiles, save, name, qtrans, nEntries):
	global W, Mp, hW, hist_W, hW_max, num_events, energies, hW_stack, stack_can, bins, max_bin_w, min_bin_w

	W, hW, hW_max, num_events= [], [], [], []
	Mp = 0.938

	##DEFINE BOUNDS ON THE ELASTIC PEAK##
	bins = 200
	# el_peak_min_x = 0.96
	# el_peak_max_x = 1.1

	#Canvas and Information to put on plot
	legend = TLegend(0.1, 0.6, 0.38, 0.9)
	legend.SetHeader("Beam current and # of Events in Peak Region", "C")
	legend.SetTextSize(0.02)

	#Run info to put in plot
	run_info = TPaveText(0.55,0.75,0.89,0.89,"NDC")
	#run_info.AddText("Momentum transfer, Q^{2} = " + str(np.round(q_squared, 3)) + " (GeV/c)^{2}")
	run_info.AddText("Target: " + TargetFromRun(int(runnum[0])))
	run_info.AddText("N Entries: " + str(nEntries[0]))
	run_info.AddText("Required hits on track: 3+ ")
	run_info.AddText("Shower & Pre-shower cuts: " + str(CutFromRun(int(runnum[0]))[0]) + " < E/p < " + str(CutFromRun(int(runnum[0]))[1]))
	run_info.SetTextSize(.025)
	run_info.SetFillColor(0)

	lineColor = 1	

	if(len(runnum) == 1):
		stack_title = "BigBite Elastic Peak"
	elif(len(runnum) > 1):
		stack_title = "BigBite Elastic Peaks"

	stack_can = TCanvas("stack_can", stack_title)
	hW_stack = THStack("hW_stack", stack_title)


#####################################################################################
###******* WHERE THE WORK IS DONE -------- CALCULATING INVARIANT MASS W **********###
#####################################################################################


	for run in runnum:
		index = runnum.index(run)
		min_bin_w = WBinFromRun(int(run))[0]
		max_bin_w = WBinFromRun(int(run))[1]

		cut_min = CutFromRun(int(run))[0]
		cut_max = CutFromRun(int(run))[1]

		E = BeamEnergyFromRun(int(run))

		el_peak_min_x = PeakMinMaxFromRun(int(run))[0]
		el_peak_max_x = PeakMinMaxFromRun(int(run))[1]

		#rootFile = "/Users/john/UVa/SBS/analysis/rootFiles/gmn_replayed_" + run + "_stream0_seg0_0.root"
		#rootFile = "/Users/john/UVa/SBS/inv_mass/replayed_11-21-21/bbgem_replayed_" + runnum + "_stream0_seg0_3.root"


		tCh = TChain("T")
		tCh.Add(rootFiles[index])


		hist_W = TH1F('hist_W_' + run , 'Invariant Mass (Elastic/Inelastic Peaks)', bins, min_bin_w, max_bin_w)
		

		cut1 = TCut("bb.gem.track.nhits[0] > 3")
		cut2 = TCut("(bb.sh.e + bb.ps.e)/bb.tr.p > " + str(cut_min) + " && (bb.ps.e + bb.ps.e)/bb.tr.p < " + str(cut_max))
		cut3 = TCut("bb.ps.e>0.1")

		all_cuts = cut1 + cut2 + cut3

		###Variables: E = Beam Energy; bb.tr.p[0] = scattered electron energy; bb.tr.pz[0] = scattered electron energy in "z"
		###Mp = mass of proton

		eq = "(sqrt(%f^2 + 2*%f*(%f - bb.tr.p[0]) - 2*%f*bb.tr.p[0]*(1 - cos(acos(bb.tr.pz[0]/bb.tr.p[0])))))>> hist_W_" + run
		

		tCh.Draw((eq  % (Mp, Mp, E, E)), all_cuts)

		num_events.append(hist_W.Integral(hist_W.GetXaxis().FindBin(el_peak_min_x), hist_W.GetXaxis().FindBin(el_peak_max_x)))

		#PLOT DISPLAY STUFF#
		if(lineColor == 5):
			lineColor+=1
		hist_W.SetLineColor(lineColor)
		lineColor+=1

		#legend.AddEntry(hist_W, "Current: " + str(CurrentFromRun(int(run))) + " \mu A; " + str(int(num_events[index])) + ' ' + " events")
		legend.AddEntry(hist_W, "Current: " + str(CurrentFromRun(int(run))) + " \mu A; ")

		##Store calculated variables in arrays##
		hW_stack.Add(hist_W)
		hW.append(hist_W)
		hW_max.append(hW[index].GetMaximum())



		# x_max = hist_W.GetMaximumBin()
		# x_Wval = ((max_bin_w - min_bin_w)/bins)*(x_max) + min_bin_w
		# print("x_Wval: ", x_Wval)



		hist_W.GetXaxis().SetRange(80, 120)
		
		W.append(getXfromBin(hist_W.GetMaximumBin()))

		hist_W.GetXaxis().SetRange(0, bins)


	#Line info to show peak cuts
	line_min = TLine(el_peak_min_x, 0, el_peak_min_x, max(hW_max)*1.05)
	line_max = TLine(el_peak_max_x, 0, el_peak_max_x, max(hW_max)*1.05)
	line_W = TLine(np.mean(W), 0, np.mean(W), max(hW_max)*1.05)
	line_W.SetLineStyle(10)
	line_min.SetLineStyle(6)
	#line_min.SetLineColor(kGray)
	line_max.SetLineStyle(6)
	#line_max.SetLineColor(kGray)
	legend.AddEntry(line_W, "Average elastic peak, W = " + str(np.round(np.mean(W), 4)), "L")
	legend.AddEntry(line_max, "Region for elastic peak", "L")
	legend.SetTextSize(0.02)


	hW_stack.Draw("nostack")
	hW_stack.GetXaxis().SetLabelOffset(0.02)
	legend.Draw("same")
	line_min.Draw("same")
	line_max.Draw("same")
	line_W.Draw("same")
	run_info.Draw()
	hW_stack.GetXaxis().SetTitle("Invariant mass, W [GeV/c]")
	hW_stack.GetXaxis().SetTitleOffset(1.3)
	hW_stack.GetYaxis().SetTitle("Counts [N]")
	stack_can.Update()


	if(save):
		if(name == None):
			stack_can.Print("./elastic_peak_plot.pdf")
		else:
			stack_can.Print("./" + name + ".pdf")
	
	return 1

