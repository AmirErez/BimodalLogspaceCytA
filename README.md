# BimodalLogspaceCytA
Bimodal distributions in log space: submission to Cytometry A 

Contents:

* Manuscript: Latex code for manuscript
* Data: Experimental data used in the manuscript
	* 20140920-OT1-dynamics	- Data used for Fig. 3 in the manuscript with MIFlowCyt document
	* FeinermanScience2008	- Previously published data used for Figs. 4,5,7 in the manuscript.
				  From Feinerman et al., Science Vol. 321, Issue 5892, pp. 1081-1084 Aug 2008
* Scripts: MATLAB scripts to make figures
	* LogspacePaperPlots.m 	- Makes all figures in the paper. When done, run CopyFigures.sh to copy the figures to Manuscript/
        * CopyFigures.sh		- Copies figures from Scripts/Plots/ to Manuscript/ directory, renaming them in the process.
	* ReadERKData_Amir.m	- Reads FeinermanScience2008 data and makes figures. Called from LogspacePaperPlots.m	
