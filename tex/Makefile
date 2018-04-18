# latex Makefile
LATEX=latex
BIBTEX=bibtex
DVIPS=dvips
PS2PDF=ps2pdf
TARGETS=column_derivation.ps

all: clean $(TARGETS) pdf

pdf: column_derivation.pdf

clean: 
	@rm -f *.aux *.bbl *.blg *.dvi *.log

column_derivation.dvi: column_derivation.tex
	${LATEX} column_derivation.tex

ps: column_derivation.ps

column_derivation.ps: column_derivation.dvi column_derivation.bib
	${BIBTEX} column_derivation
	${LATEX} column_derivation.tex
	${BIBTEX} column_derivation
	${LATEX} column_derivation.tex
	${LATEX} column_derivation.tex
	${DVIPS} column_derivation.dvi -o column_derivation.ps
	#${DVIPS} column_derivation.dvi -o column_derivation1.ps
	#awk -f awk.f column_derivation1.ps > column_derivation.ps
	#cp temp.ps column_derivation.ps 

column_derivation.pdf: column_derivation.ps
	${PS2PDF} column_derivation.ps

bw:
	@rm -f apjform_bw.aux apjform_bw.bbl apjform_bw.blg apjform_bw.dvi apjform_bw.log
	${LATEX} apjform_bw.tex
	${LATEX} apjform_bw.tex
	${DVIPS} apjform_bw.dvi -o apjform_bw.ps
	${PS2PDF} apjform_bw.ps

apj:
	${LATEX} apjform.tex
	${LATEX} apjform.tex
	${DVIPS} apjform.dvi -o apjform.ps
	${PS2PDF} apjform.ps

tar:
	cp OutflowOverlay2426.ps outflows_on_CO32.eps outflows_on_BGPS.eps outflows_on_CO.eps outflows_on_IRAS100.eps outflows_on_24UM.eps  S201_BGPSon24UM.eps S201_BGPSon8UM.eps AFGL4029_BGPSon24UM.eps AFGL4029_BGPSon8UM.eps LWCas_BGPSon24UM.eps LWCas_BGPSon8UM.eps W5S_BGPSon24UM.eps W5S_BGPSon8UM.eps W5SW_BGPSon24UM.eps W5SW_BGPSon8UM.eps W5W_BGPSon24UM.eps W5W_BGPSon8UM.eps W5NW_BGPSon24UM.eps W5NW_BGPSon8UM.eps OutflowHistograms.ps outflowoverlays/OutflowOverlay_8um_12.eps outflowoverlays/OutflowOverlay_24_12.eps outflowoverlays/OutflowOverlay_bgps_12.eps outflowplots/OutflowSpectra12.eps column_derivation_mnras/
	cp column_derivation.tex outflowsumtable_mnras.tex outflowtable_mnras.tex preface_mnras.tex regiontable_mnras.tex column_derivation_mnras/
	cp column_derivation.bib column_derivation_mnras/
	tar -czf column_derivation.tar.gz column_derivation_mnras
	tar -czf column_derivation_overlays.tar.gz OnlineFigures



