The Pol2011 directory is organized as follows:
==============================================

Most macros include the file "../interface/commonVar.h" in which
common settings are defined: pT, |y| bins, kinematical cuts, bins in
cosTheta and phi, etc.

In order to use the macros the user need to create two
sub-directories:
mkdir Figures
mkdir RootFiles

1.) runData.C is the steering macro for PolData.C/h which loops over
the TTree produced from the JPsiAnalyzerPAT.cc macro, applies a series
of muon quality and kinematical cuts, and stores the TLorentzVector of
the two muons of the selected events, as well as a series of mass
histograms in an output file. The default name of the root output file is:

	   RootFiles/selEvents_data_Ups_1Aug2011.root

2.) runMassFit.C is the steering macro for upsilon_2StepFit.C which
loops over all pT / |y| bins and performs a fit to the invariant mass
region in the upsilon mass region. The fit is a 2-step fit in which
the BG shape and normalization is fixed from the L and R mass sideband
windows which is then imposed in the fit to the signal region, defined
inbetween the L and R sidebands. The signal consists of 3 Crystal ball
functions with common tail parameters (alpha and n) where only the
mass and width of the 1S are left free and the corresponding mass and
widths of the 2S and 3S are fixed by the respective mass ratios as
given in the PDG2010 tables. 
The CB parameters (alpha and n) are fixed from the pT integrated bin,
for each rapidity separately, and are imposed on the pT differential
bins. 

The fitted TF1 objects for the signal (3 CB functions) as well as the
BG function are stored in output files called

   RootFiles/data_Ups_rap1_pT1.root

etc. The defintions of the pT and rapidity bins are defined in
"../interface/commonVar.h".

This macro also saves the graphical representation of the fit in
corresponding pdf files.

2a.) In order to visualize the series of pdf files created in the
previous step the user can go into the subdirectory "../latex" and
execute the command 
pdflatex massFits.tex


3.) runCopyTreeEntries is the steering macro for
CopyTreeEntries.C. This macro loops again over the same pT and y bins
and adds to the previously created output file
(data_Ups_rap1_pT1.root) the TLorentzVectors of the 2 leptons,
specific to that pT and y bin. It furthermore, calculates the width of
the L and R sideband windows and projects corresponding events into
the L and R (cosTheta,phi) 2D histos.

3a.) the macro PlotCosThetaPhiBG.C plots the individual (cosTheta,
phi) distributions of the L and R mass sidebands. The individual
figures can then be assembled into one pdf file using the file 
../latex/cosThetaPhi_BG.tex

4.) runTrimEventContent.C is the steering macro for
TrimEventContent.C. This macro loops over the events in a given pT and
y bin and stores in the "data" TTree only those events which are
selected within a certain nSigma window around the signal polemass. It
also adds the Left and Right Sideband mass windows in a given fraction
and stores one output BG histogram, as a function of (cosTheta,
phi). Furthermore, it also projects the (cosTheta, phi) distribution
in a 2D histogram of the signal events; the fraction of BG in the
nSigma window around the signal is stored in a 1D histogram. Inputs to
the macro are the 
    *) fraction with which the L BG should be added to the right one
    *) nSigma within which the events should be projected
    *) which state is being considered: //[0]... 1S, [1]... 2S, [2]... 3S

