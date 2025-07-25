#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
   SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Subsection:: *)
(*Setup*)


(* ::Text:: *)
(*Load Mathematica wrapper in every kernel. Remove the "Quiet" if unsure about the parallel evaluation.*)


BaseDirectory=FileNameJoin[{NotebookDirectory[], "..", "..", "mathematica"}];
SetDirectory[BaseDirectory];<<"EpsteinZeta.wl";
ParallelEvaluate[SetDirectory[BaseDirectory];Quiet[<<"EpsteinZeta.wl"]];


(* ::Text:: *)
(*Set 0^0 in every kernel.*)


Unprotect[Power];
Power[0|0.|0*I|0.*I,0|0.|0*I|0.*I]:=1;
Power[z_,0]/;Abs[z]==0:=1;
Protect[Power];


ParallelEvaluate[
Unprotect[Power];
Power[0|0.|0*I|0.*I,0|0.|0*I|0.*I]:=1;
Power[z_,0]/;Abs[z]==0:=1;
Protect[Power];
];


(* ::Text:: *)
(*Minimum of absolute and relative error.*)


Error[exact_,approx_]:=Min[Abs[exact-approx],If[PossibleZeroQ[exact],1,Abs[exact-approx]/Abs[exact]]];


(* ::Text:: *)
(*Install MaTeX for LATEX figure annotations. If MaTeX is not installed already un-comment the ResourceFunction line to install it.*)


(*ResourceFunction["MaTeXInstall"][];*)
Needs["MaTeX`"];


(* ::Text:: *)
(*Set options for pretty plots.*)


SetOptions[RegionPlot,TicksStyle->{FontSize->18},LabelStyle->{FontSize->18}];
SetOptions[MaTeX,"Preamble"->{"\\usepackage{xcolor}","\\usepackage{bm}"}];
SetOptions[MaTeX,FontSize->18];


(* ::Subsection:: *)
(*Define test problem*)


(* ::Text:: *)
(*Set constants.*)


\[Sigma]0 = 100; \[Delta]0 = 0.5; a0 = -500; b0 = 500; plotMin = -2*\[Sigma]0; plotMax = -plotMin; \[Nu]0 = 9/10;


(* ::Text:: *)
(*Define the algebraic singularity, the Gaussian, their fourier transforms and the derivatives of the Gaussian.*)


Singularity1D[y_,\[Nu]_]=(y^2)^(-\[Nu]/2);
Gaussian1D[y_,\[Sigma]_]=Exp[-Pi*(y/\[Sigma])^2];
SingularityFourier1D[z_,\[Nu]_]=FourierTransform[Singularity1D[y,\[Nu]],y,z,FourierParameters->{0, -2 Pi}];
GaussianFourier1D[z_,\[Sigma]_]=FourierTransform[Gaussian1D[y,\[Sigma]],y,z,FourierParameters->{0, -2 Pi}];
GaussianDerivatives1D[y_,\[Sigma]_,n_]:=Total[Table[(-Pi/\[Sigma]^2)^(n-k)Binomial[n,k](n-k)!/(n-2k)!(2y)^(n-2k)Gaussian1D[y,\[Sigma]] ,{k,0,n/2}]];


(* ::Text:: *)
(*Define sum and integral.*)


Sum1D[x_,\[Nu]_,\[Sigma]_] :=Total[Table[
	Piecewise[{{0,PossibleZeroQ[y-x]},{Singularity1D[y-x,\[Nu]]*Gaussian1D[y,\[Sigma]],Not[PossibleZeroQ[y-x]]}}]
	,{y,Ceiling[a0+\[Delta]0],Floor[b0+\[Delta]0]}
]]
Integral1D[x_,\[Nu]_,\[Sigma]_] = InverseFourierTransform[GaussianFourier1D[z,\[Sigma]]*SingularityFourier1D[z,\[Nu]],z,x,FourierParameters->{0, -2 Pi}];


(* ::Subsection:: *)
(*Singular Euler-Maclaurin expansion*)


(* ::Text:: *)
(*Define Euler-Maclaurin expansion.*)


SEM1D[x_,\[Nu]_,\[Sigma]_,n_] := Integral1D[x,\[Nu],\[Sigma]]+Total[Table[1/j!/(-2 Pi I)^j EpsteinZetaRegDer[\[Nu],{{1}},{0},{0},{j}]GaussianDerivatives1D[x,\[Sigma],j],{j,0,n,2}]]


(* ::Text:: *)
(*Show correction of integral and the corresponding error.*)


GraphicsRow[{
	Show[
	ListLinePlot[{
		{#,SEM1D[#,\[Nu]0,\[Sigma]0,0]}&/@Range[plotMin,plotMax,\[Sigma]0/20]
		,{#,Integral1D[#,\[Nu]0,\[Sigma]0]}&/@Range[plotMin,plotMax,\[Sigma]0/20]
		},PlotTheme->"Scientific"
		,PlotLegends->Placed[Map[MaTeX,{"\\text{SEM}","\\text{Integral}"}],{0.8,0.8}]
		,PlotRange->All
		,FrameTicks->{{{#,MaTeX[ToString[#]]}&/@Range[0,30,10],{#,None}&/@Range[0,30,10]},{{{-2*\[Sigma]0, MaTeX["-2\\sigma"]},{-\[Sigma]0, MaTeX["-\\sigma"]},{0,MaTeX["0"]},{\[Sigma]0, MaTeX["\\sigma"]},{2*\[Sigma]0, MaTeX["2\\sigma"]}},{#,None}&/@Range[-2\[Sigma]0,2\[Sigma]0,\[Sigma]0]}}
		,FrameLabel-> {{MaTeX["\\sum_{y=a+1}^{b}\\frac{e^{-\\pi(y/\\sigma)^2}}{|y-x|^{\\nu}}"],None},{MaTeX["x"],None}}
	],ListPlot[
		{#,Sum1D[#,\[Nu]0,\[Sigma]0]}&/@Range[plotMin,plotMax,\[Sigma]0/10]
		,PlotLegends->Placed[{MaTeX["\\text{Sum}"]},{0.8,0.6}]
		,PlotStyle->{Red}
	],ImageSize->500
	],ListLogPlot[Table[
		ParallelTable[{x,Error[Sum1D[x,\[Nu]0,\[Sigma]0],SEM1D[x,\[Nu]0,\[Sigma]0,n]]},{x,plotMin,plotMax,\[Sigma]0/20}]
		,{n,0,4,2}
		],PlotTheme->"Scientific"
		,PlotLegends->Placed[Map[MaTeX["\\text{Order: }"<>ToString[#]]&,2 Range[0,2]],{After,Center}]
		,PlotRange->All
		,FrameTicks->{{{10^#,MaTeX["10^{"<>ToString[#]<>"}"]}&/@Range[-16,-7,3],{10^#,None}&/@Range[-16,-7,3]},{{{-2*\[Sigma]0, MaTeX["-2\\sigma"]},{-\[Sigma]0, MaTeX["-\\sigma"]},{0,MaTeX["0"]},{\[Sigma]0, MaTeX["\\sigma"]},{2*\[Sigma]0, MaTeX["2\\sigma"]}},{#,None}&/@Range[-2\[Sigma]0,2\[Sigma]0,\[Sigma]0]}}
		,FrameLabel-> {{MaTeX["\\text{Error}"],None},{MaTeX["x"],None}}
		,ImageSize->600
	]
	},ImageSize->1150
]
