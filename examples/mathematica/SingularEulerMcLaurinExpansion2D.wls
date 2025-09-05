#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
   SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Singular Euler-Maclaurin Expansion in 2D*)


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
SetOptions[MaTeX,"Preamble"->{"\\usepackage{xcolor}","\\usepackage{bm}","\\usepackage{dsfont}"}];
SetOptions[MaTeX,FontSize->18];


(* ::Subsection:: *)
(*Define functions*)


(* ::Text:: *)
(*Compiled Sum		\!\(\**)
(*SubscriptBox["\[Sum]", *)
(*RowBox[{*)
(*RowBox[{"(", *)
(*RowBox[{*)
(*SubscriptBox["y", *)
(*RowBox[{"1", ","}]], *)
(*SubscriptBox["y", "2"]}], ")"}], "\[Element]", *)
(*SuperscriptBox[*)
(*TemplateBox[{},*)
(*"Integers"], "2"]}]]\( *)
(*\*SuperscriptBox[\(( *)
(*\*SuperscriptBox[\(( *)
(*\*SubscriptBox[\(y\), \(1\)] - *)
(*\*SubscriptBox[\(x\), \(1\)])\), \(2\)] + *)
(*\*SuperscriptBox[\(( *)
(*\*SubscriptBox[\(y\), \(2\)] - *)
(*\*SubscriptBox[\(x\), \(2\)])\), \(2\)])\), \(\(-\[Nu]\)/2\)] *)
(*\*SuperscriptBox[\(e\), \(\(-\(\[Pi]( *)
(*\*SubsuperscriptBox[\(y\), \(1\), \(2\)] + *)
(*\*SubsuperscriptBox[\(y\), \(2\), \(2\)])\)\)/*)
(*\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)]\)\)*)


Sum2DCompiled=Compile[{
	{x1, _Real}
	,{x2, _Real}
	,{\[Nu], _Real}
	,{\[Sigma], _Real}
  ,{lim, _Integer}
	},Module[
		{res}
		,res=0.
		;Do[
			res+=Piecewise[{
				{0,(y1-x1)^2+(y2-x2)^2==0}
				,{((y1-x1)^2+(y2-x2)^2)^(-\[Nu]/2)Exp[-Pi*(y1^2+y2^2)/\[Sigma]^2],True}
			}]
	                   ,{y1,-lim,lim}
			,{y2,-lim,lim}
		]
		;Return[res]
	]
];
Sum2D[{x1_,x2_},\[Nu]_,\[Sigma]_,lim_:100]:=Sum2DCompiled[N[x1],N[x2],N[\[Nu]],N[\[Sigma]],lim];


(* ::Text:: *)
(*Analytic expression for the integral		 \!\( *)
(*\*SubsuperscriptBox[\(\[Integral]\), \(-\[Infinity]\), \(\[Infinity]\)]\( *)
(*\*SuperscriptBox[\(( *)
(*\*SuperscriptBox[\(( *)
(*\*SubscriptBox[\(y\), \(1\)] - *)
(*\*SubscriptBox[\(x\), \(1\)])\), \(2\)] + *)
(*\*SuperscriptBox[\(( *)
(*\*SubscriptBox[\(y\), \(2\)] - *)
(*\*SubscriptBox[\(x\), \(2\)])\), \(2\)])\), \(\(-\[Nu]\)/2\)] *)
(*\*SuperscriptBox[\(e\), \(\(-\(\[Pi]( *)
(*\*SubsuperscriptBox[\(y\), \(1\), \(2\)] + *)
(*\*SubsuperscriptBox[\(y\), \(2\), \(2\)])\)\)/*)
(*\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)] \[DifferentialD]\ *)
(*\*SubscriptBox[\(y\), \(1\)] \[DifferentialD]\ *)
(*\*SubscriptBox[\(y\), \(2\)]\)\)*)


Integral2D[{x1_,x2_},\[Nu]_,\[Sigma]_]:=\[Pi]^(\[Nu]/2) (\[Sigma]^2)^(1-\[Nu]/2) Gamma[1-\[Nu]/2] LaguerreL[-(\[Nu]/2),-((\[Pi] (x1^2+x2^2))/\[Sigma]^2)]


(* ::Text:: *)
(*Analytic expression for the partial derivatives of the gaussian	 e^(-\[Pi](\!\( *)
(*\*SubsuperscriptBox[\(y\), \(1\), \(2\)] + *)
(*\*SubsuperscriptBox[\(y\), \(2\), \(2\)]\))/\[Sigma]^2)*)


Gaussian2D[{y1_,y2_},\[Sigma]_]=Exp[-((\[Pi] (y1^2+y2^2))/\[Sigma]^2)];
GaussianDerivatives2D[{y1_,y2_},\[Sigma]_,\[Alpha]_]:=Gaussian2D[{y1,y2},\[Sigma]]Total[Table[
	Apply[Times,(-Pi/\[Sigma]^2)^(\[Alpha]-\[Beta])Binomial[\[Alpha],\[Beta]](\[Alpha]-\[Beta])!/(\[Alpha]-2\[Beta])!(2{y1,y2})^(\[Alpha]-2\[Beta])]
		,{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]/2]]}
]]


(* ::Text:: *)
(*SEM Operator*)


SEM2D[{x1_,x2_},\[Nu]_,\[Sigma]_,n_:0] := Integral2D[{x1,x2},\[Nu],\[Sigma]]+Total[Table[
	Apply[Times,1/\[Beta]!/(-2 Pi I)^\[Beta]] Re@EpsteinZetaRegDer[\[Nu],IdentityMatrix[2],ConstantArray[0,2],ConstantArray[0,2],\[Beta]]GaussianDerivatives2D[{x1,x2},\[Sigma],\[Beta]]
		,{\[Beta],Tuples[Range[0,n,2],2]}
	]]


(* ::Subsection:: *)
(*Error Analysis*)


(* ::Text:: *)
(*Set test parameters.*)


\[Sigma]0 = 10; \[Nu]0 = 9/10;


(* ::Text:: *)
(*SEM and Error for chosen parameters.*)


GraphicsRow[{
	Show[
	ListLinePlot[{
		{#,SEM2D[{#,1},\[Nu]0,\[Sigma]0,0]}&/@Range[-2\[Sigma]0,2\[Sigma]0,\[Sigma]0/20]
		,{#,Integral2D[{#,1},\[Nu]0,\[Sigma]0]}&/@Range[-2\[Sigma]0,2\[Sigma]0,\[Sigma]0/20]
		},PlotTheme->"Scientific"
		,PlotLegends->Placed[Map[MaTeX,{"\\text{SEM}","\\text{Integral}"}],{0.8,0.8}]
		,PlotRange->All
		,GridLines->{{-2\[Sigma]0, -\[Sigma]0,0,\[Sigma]0,2\[Sigma]0},Range[0,30,10]}
		,FrameTicks->{{{#,MaTeX[ToString[#]]}&/@Range[0,30,10],{#,None}&/@Range[0,30,10]},{{{-2*\[Sigma]0, MaTeX["-2\\sigma"]},{-\[Sigma]0, MaTeX["-\\sigma"]},{0,MaTeX["0"]},{\[Sigma]0, MaTeX["\\sigma"]},{2*\[Sigma]0, MaTeX["2\\sigma"]}},{#,None}&/@Range[-2\[Sigma]0,2\[Sigma]0,\[Sigma]0]}}
		,FrameLabel-> {{MaTeX["\\sum_{\\bm y \\in \\mathds{Z}^d}\\frac{e^{-\\pi(\\bm y/\\sigma)^2}}{|\\bm y-(x,1)|^{\\nu}}"],None},{MaTeX["x"],None}}
	],ListPlot[
		{#,Sum2D[{#,1},\[Nu]0,\[Sigma]0]}&/@Range[-2\[Sigma]0,2\[Sigma]0,\[Sigma]0/10]
		,PlotLegends->Placed[{MaTeX["\\text{Sum}"]},{0.8,0.6}]
		,PlotStyle->{Red}
	],ImageSize->500
	],ListLogPlot[Table[
		ParallelTable[{x,Error[Sum2D[{x,1},\[Nu]0,\[Sigma]0],SEM2D[{x,1},\[Nu]0,\[Sigma]0,n]]},{x,-2\[Sigma]0,2\[Sigma]0,\[Sigma]0/10}]
		,{n,0,10,2}
		],PlotTheme->"Scientific"
		,PlotLegends->Placed[Map[MaTeX["\\text{Order: }"<>ToString[#]]&,2 Range[0,10]],{After,Center}]
		,PlotRange->All
	  ,Joined->True
		,GridLines->{{-2\[Sigma]0, -\[Sigma]0,0,\[Sigma]0,2\[Sigma]0},10^Range[-16,2,3]}
		,FrameTicks->{{{10^#,MaTeX["10^{"<>ToString[#]<>"}"]}&/@Range[-16,2,3],{10^#,None}&/@Range[-16,-7,3]},{{{-2*\[Sigma]0, MaTeX["-2\\sigma"]},{-\[Sigma]0, MaTeX["-\\sigma"]},{0,MaTeX["0"]},{\[Sigma]0, MaTeX["\\sigma"]},{2*\[Sigma]0, MaTeX["2\\sigma"]}},{#,None}&/@Range[-2\[Sigma]0,2\[Sigma]0,\[Sigma]0]}}
		,FrameLabel-> {{MaTeX["\\text{Error}"],None},{MaTeX["x"],None}}
		,ImageSize->600
	]
	},ImageSize->1150
]
