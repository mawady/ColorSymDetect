# Color-based Reflection Symmetry Detection

Unsupervised Reflection Symmetry Detection using Log-Gabor Filters and Textural & Color Histograms
Release Date: July 2017

Users of this software are encouraged to cite the following article:
Mohamed Elawady, Christophe Ducottet, Olivier Alata, Cécile Barat, and Philippe Colantoni. "Wavelet-based Reflection Symmetry Detection via Textural and Color Histograms." arXiv preprint arXiv:1707.02931 (2017). [https://arxiv.org/abs/1707.02931]

The sample photos are kindly provided by Aesthetic Visual Analysis (AVA) dataset.

Contact: Mohamed Elawady (mohamed [dot] elawady [at] univ-st-etienne [dot] fr)
http://github.com/mawady/ColorSymDetect


FILE:
main.m - code for computer symmetry on an image and display the candidates' info (axes, voting representation)
symBilOurCentLogGaborHSV.m - code for main function to compute symmetry on an image

OUTPUT FORMAT:
SymOcLgHSV :	(Cx6 Matrix) Information of symmetry candidates (sorted decendly)
				C --> Number of symmetry candidates
				Columns: x1,y1,x2,y2,score,normalized_score
				x1,y1,x2,y2 --> 2D coordinates of line endpoints representing a symmetry axis candidate
				score --> symmetry measure of a symmetry axis candidate
				normalized_score --> symmetry measure mormalized by maximum score				
voteMap	:		(MxN Matrix) Probability estimation of symmetry voting space (polar coordinate system)
				M --> number of bins for angular component
				N --> number of bins for displacement component
