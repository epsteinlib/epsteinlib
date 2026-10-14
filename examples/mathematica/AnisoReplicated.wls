#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2026 Jonathan K. Hartlaub <contact@jkhartlaub.com>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Anisotropic Epstein Zeta Functions Replicated*)


(* ::Text:: *)
(*This script generates the results for all numerical experiments presented in the arXiv preprint "Computation of anisotropic singular sums from high-order derivatives of Epstein zeta functions", arXiv:2609.28282, 2026. The grids are coarser than in the manuscript so that the evaluation of this script takes less than a minute on a standard laptop equipped with 12 cores.*)


(* ::Subsection::Closed:: *)
(*Setup Packages*)


(* ::Text:: *)
(*Please ensure that EpsteinLib is correctly installed and accessible before running the script. That is*)
(*	*)
(*	1. 	go to the project: cd epsteinlib*)
(*	2.	build library: meson setup build*)
(*	3.	compile library: meson compile -C build*)
(*	*)
(*To test the C library, run meson test -C build.*)


(* ::Text:: *)
(*Load the Mathematica wrapper of the Epstein zeta function, so that it can be used in parallel on all kernels.*)


SetDirectory[NotebookDirectory[]<>"../../mathematica"];
ParallelEvaluate[Quiet[<<"EpsteinZeta.wl"]];
<<"EpsteinZeta.wl";


(* ::Text:: *)
(*Install MaTeX for LaTeX figure annotations. If MaTeX is not already installed, un-comment the ResourceFunction line to install it.*)


(*ResourceFunction["MaTeXInstall"][];*)
Needs["MaTeX`"];


(* ::Text:: *)
(*Set options for pretty plots.*)


SetOptions[RegionPlot,TicksStyle->{FontSize->18},LabelStyle->{FontSize->18}];
SetOptions[MaTeX,"Preamble"->{"\\usepackage{xcolor}","\\usepackage{bm}"}];
SetOptions[MaTeX,FontSize->18];


(* ::Subsection::Closed:: *)
(*Auxiliary Functions*)


(* ::Text:: *)
(*Minimum of absolute and relative error.*)


Error[true_,est_]:=Module[{errAbs,errRel}
	,errAbs=Abs[true-SetPrecision[est,100]]
	;errRel=If[PossibleZeroQ[true]
	,If[PossibleZeroQ[est],0,Infinity]
	,errAbs/Abs[true]
	]
	;Min[errAbs,errRel]
]


(* ::Text:: *)
(*Colors for phases of complex values.*)


ArgColor[arg_]=Hue[(arg+Pi)/(2 Pi)];


(* ::Section:: *)
(*Figure 1*)


(* ::Subsection::Closed:: *)
(*Prepare*)


(* ::Text:: *)
(*Fix parameters for the (regularized) anisotropic Epstein zeta functions.*)


zfunc[y1_,y2_]:=EpsteinZetaAniso[3.5,IdentityMatrix@2,{1,1}/3,{y1,y2},{0,2}]
zfuncreg[y1_,y2_]:=EpsteinZetaAnisoReg[3.5,IdentityMatrix@2,{1,1}/3,{y1,y2},{0,2}]


(* ::Text:: *)
(*Generate anisotropic Epstein zeta and regularized anisotropic Epstein zeta values for fixed \[Nu], identity a, all-even \[alpha], all-equal x, and y in a 2D grid [-0.7, 0.7] x [-0.7, 0.7].*)


offset=0.01;xmin =- 0.7-offset;xmax = 0.7+offset; ymin = -0.7-offset; ymax = 0.7+offset;
stepsize = 0.07; (* set to 0.007 to match the resolution of the manuscript *)
zanitab = Flatten[ParallelTable[{x,y,Abs[zfunc[x,y]]},{x,xmin,xmax, stepsize},{y,ymin,ymax,stepsize}],1];
stepsize = 0.07; (* the regularized function is smooth, so the coarse stepsize suffices *)
zaniregtab = Flatten[ParallelTable[{x,y,Abs[zfuncreg[x,y]]},{x,xmin,xmax, stepsize},{y,ymin,ymax,stepsize}],1];


(* ::Subsection::Closed:: *)
(*Compose*)


plt=GraphicsRow[{
	Show[
		ListPlot3D[
			zanitab
			,AxesLabel->{MaTeX@"y_1",MaTeX@"y_2"}
			,PlotRange->{Automatic,Automatic,{-0.8,6.8}}
			,PlotStyle->{Opacity[0.8]}
			,ColorFunctionScaling->False
			,ColorFunction->Function[{x,y,z},ArgColor[Arg[zfunc[x,y]]]]
			,RegionFunction->Function[{x,y,z},Norm[{x,y}]>0.015]
			,Mesh->None
			,PlotLabel->MaTeX["|Z_{\\mathbb Z^2, 7/2,(0,2)}((1,1)/3,\\bm y)|"]
		]
		,ListPointPlot3D[
			{{0,0,Abs@zfunc[0,0]}}
			,PlotStyle->{Opacity[0.8],ArgColor[zfunc[0,0]],PointSize[0.015]}
		]
	]
	,Show[
		ListPlot3D[
			zaniregtab
			,AxesLabel->{MaTeX@"y_1",MaTeX@"y_2"}
			,PlotRange->{Automatic,Automatic,{1.5,6.5}}
			,ColorFunction->Function[{x,y,z},ArgColor[Arg[zfuncreg[x,y]]]]
			,PlotStyle->{Opacity[0.8]}
			,ColorFunctionScaling->False
			,Mesh->None
			,PlotLabel->MaTeX["|Z^{\\rm reg}_{\\mathbb Z^2, 7/2,(0,2)}((1,1)/3,\\bm y)|"]
		]
	]
	,BarLegend[
		{ArgColor[#]&,{-Pi-1/100,Pi+1/100}}
		,Ticks->{{-Pi,"-\[Pi]"},{-Pi/2,"-\[Pi]/2"},{0,"0"},{Pi/2,"\[Pi]/2"},{Pi,"\[Pi]"}}
	]
	}
]


(* ::Text:: *)
(**)


(* ::Subsection:: *)
(*Display*)


(* ::Text:: *)
(*Complex phase plots of the 2D anisotropic Epstein zeta function (left) and its regularized version (right) for \[Nu] = 7/2, \[CapitalLambda] = \:2124\.b2 and x =(1/3,1/3), where the height of the surface encodes magnitude and the color encodes the phase.*)


plt


(* ::Section:: *)
(*Table 1*)


(* ::Subsection::Closed:: *)
(*Prepare*)


(* ::Text:: *)
(*Set the reference associated with applying the poly-Laplacian.*)


EpsteinZetaAnisoPolyLaplaceAnalytic[\[Nu]_,a_,x_,y_,n_]:=EpsteinZeta[\[Nu]-2n,a,x,y]


(* ::Text:: *)
(*Rebuild the poly-Laplacian through sums of anisotropic Epstein zeta functions.*)


EpsteinZetaAnisoPolyLaplace[\[Nu]_,a_,x_,y_,n_]:=n!Total@Table[
	1/Apply[Times,\[Alpha]!]EpsteinZetaAniso[\[Nu],a,x,y,2\[Alpha]]
	,{\[Alpha],Select[Tuples[Range[0,n],Length[a]],Total[#]==n&]}
]


(* ::Text:: *)
(*Generate test values.*)


vectorStepSize[d_]:=ConstantArray[1/2,d] (* set to ConstantArray[1/5,d] to replicate the result from the manuscript *)
errs=Table[
	Table[
		Max@Flatten[ParallelTable[
			Table[
				Error[
					EpsteinZetaAnisoPolyLaplaceAnalytic[\[Nu],a,x,y,n]
					,EpsteinZetaAnisoPolyLaplace[\[Nu],a,x,y,n]
				]
				,{y,Select[Tuples[Range[ConstantArray[0,d],ConstantArray[1/2,d],vectorStepSize[d]]],Not[PossibleZeroQ[Norm[#]]&&\[Nu]==d+2n]&]}
			]
			,{\[Nu],0,1,.1}
			,{x,Tuples[Range[ConstantArray[0,d],ConstantArray[1/2,d],vectorStepSize[d]]]}
			,{a,{IdentityMatrix[d]}}
			]
		]
		,{n,6}
	]
	,{d,{3,2,1}}
]


(* ::Subsection::Closed:: *)
(*Compose*)


tbl=TableForm[
	Map[ScientificForm[#,3]&,errs,{2}]
	,TableHeadings->{
		MaTeX["d="<>ToString[#]]&/@{3,2,1}
		,MaTeX["n="<>ToString[#]]&/@Range[6]
	}
]


(* ::Subsection:: *)
(*Display*)


(* ::Text:: *)
(*Maximum error of sums of anisotropic Epstein zeta functions associated with applying the poly-Laplacian \[CapitalDelta]\:207f over a grid of values in \[Nu], x and y.*)


tbl


(* ::Section:: *)
(*Table 2*)


(* ::Subsection::Closed:: *)
(*Prepare*)


(* ::Text:: *)
(*Set the reference for non-mixed second-order derivatives for the square lattice and all-equal vector arguments.*)


EpsteinZetaAnisoAnalytic[
	\[Nu]_
	,a_/;PossibleZeroQ[Norm[a-IdentityMatrix[Length[a]]]]
	,x_/;PossibleZeroQ[Norm[x-ConstantArray[First[x],Length[x]]]]
	,y_/;PossibleZeroQ[Norm[y-ConstantArray[First[y],Length[y]]]]
	,\[Alpha]_/;AllTrue[\[Alpha],EvenQ]&&PossibleZeroQ[2-Total[\[Alpha]]]
]:=EpsteinZeta[\[Nu]-Total[\[Alpha]],a,x,y]/Length[a]


(* ::Text:: *)
(*Generate test values.*)


stepsize=1/2; (* set to 1/100 to match the resolution of the manuscript *)
errs=Table[
	{Median[#],Max[#]}&[Flatten[ParallelTable[
		Table[
			Error[
				EpsteinZetaAnisoAnalytic[\[Nu],IdentityMatrix@d,ConstantArray[x,d],ConstantArray[y,d],2UnitVector[d,1]]
				,EpsteinZetaAniso[\[Nu],IdentityMatrix@d,ConstantArray[x,d],ConstantArray[y,d],2UnitVector[d,1]]
			]
			,{y,Select[Range[0,1/2,stepsize],Not[PossibleZeroQ[#]&&\[Nu]==d+2 ]&]}
		]
		,{\[Nu],-10,10,1/10}
		,{x,0,1/2,stepsize}
		,{a,{IdentityMatrix[d]}}
	]]]
	,{d,Reverse@Range[6]}
]


(* ::Subsection::Closed:: *)
(*Compose*)


tbl=TableForm[
	Map[ScientificForm[#,3]&,Transpose[errs],{2}]
	,TableHeadings->{
		{MaTeX["E^{\\rm med}"],MaTeX["E^{\\rm max}"]}
		,MaTeX["d="<>ToString[#]]&/@Reverse@Range[6]
	}
]


(* ::Subsection:: *)
(*Display*)


(* ::Text:: *)
(*Median and maximum error of anisotropic Epstein zeta functions corresponding to non-mixed second-order derivatives for the square lattice and all-equal vector arguments over a grid of values in \[Nu], x and y.*)


tbl


(* ::Section:: *)
(*Figure 2*)


(* ::Subsection::Closed:: *)
(*Prepare*)


(* ::Text:: *)
(*Set the reference in terms of the Ramanujan tau L-function with exact trivial zeros.*)


EpsteinZetaAnisoAnalytic[
	\[Nu]_
	,IdentityMatrix[8]
	,ConstantArray[-1/2,8]
	,ConstantArray[1/2,8]
	,ConstantArray[1,8]
]:=If[IntegerQ[\[Nu]/2]&&\[Nu]<=0
	,0
	,2^(-\[Nu]/2) RamanujanTauL[\[Nu]/2]
];


(* ::Text:: *)
(*Generate test values.*)


stepsize = 1/2; (* set to 1/20 to match the resolution of the manuscript *)
errs=ParallelTable[
	{
		\[Nu]
		,N@Error[
			EpsteinZetaAnisoAnalytic[\[Nu],IdentityMatrix[8],ConstantArray[-1/2,8],ConstantArray[1/2,8],ConstantArray[1,8]]
			,EpsteinZetaAniso[\[Nu],IdentityMatrix[8],ConstantArray[-1/2,8],ConstantArray[1/2,8],ConstantArray[1,8]]
		]
	}
	,{\[Nu],-10-1/5,10+1/5,stepsize}
]


(* ::Subsection::Closed:: *)
(*Compose*)


plt=ListLogPlot[
	errs
	,PlotTheme->"Scientific"
	,PlotLabel->MaTeX["Z_{\\mathbb Z^8,\\nu,\\bm 1}(-\\bm 1/2,\\bm 1/2)"]
	,FrameLabel->{MaTeX["\\nu"],MaTeX["E"]}
	,PlotStyle->{PointSize[0.015],RGBColor[0.365248, 0.427802, 0.758297]}
	,PlotRange->{{-11-1/3,11+1/3},{2 10^(-20),5 10^(-13)}}
	,GridLines->{Range[-10,10,5],10^Range[-21,-12,2]}
	,AspectRatio->.5
	,ImageSize->400
]


(* ::Subsection:: *)
(*Display*)


(* ::Text:: *)
(*Error for a special case of the eight-dimensional anisotropic Epstein zeta function for the square lattice and all-one anisotropy in which case it is given by the Ramanujan tau L-function.*)


plt


(* ::Section:: *)
(*Figure X*)


(* ::Subsection:: *)
(*Prepare*)


(* ::Subsection:: *)
(*Compose*)


(* ::Subsection:: *)
(*Display*)
