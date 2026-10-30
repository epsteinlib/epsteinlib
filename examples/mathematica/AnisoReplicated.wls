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
(*Set 0\:2070 = 1 in every kernel.*)


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
(*Install MaTeX for LaTeX figure annotations. If MaTeX is not already installed, un-comment the ResourceFunction line to install it.*)


(*ResourceFunction["MaTeXInstall"][];*)
Needs["MaTeX`"];


(* ::Text:: *)
(*Set options for pretty plots.*)


SetOptions[RegionPlot,TicksStyle->{FontSize->18},LabelStyle->{FontSize->18}];
SetOptions[MaTeX,"Preamble"->{"\\usepackage{xcolor}","\\usepackage{bm}"}];
SetOptions[MaTeX,FontSize->18];


(* ::Text:: *)
(*Darken Blue.*)


Unprotect[Blue];
Blue=RGBColor[0.365248, 0.427802, 0.758297];
Protect[Blue];


(* ::Subsection::Closed:: *)
(*Auxiliary Functions and Constants*)


(* ::Text:: *)
(*Minimum of absolute and relative error.*)


Options[Error]={InputPrecision->200,OutputPrecision->MachinePrecision};
Error[true_,est_,OptionsPattern[]]:=Module[{errAbs,errRel}
	,errAbs=Abs[true-SetPrecision[est,OptionValue[InputPrecision]]]
	;errRel=If[PossibleZeroQ[true]
	,If[PossibleZeroQ[est],0,Infinity]
	,errAbs/Abs[true]
	]
	;SetPrecision[Min[errAbs,errRel],OptionValue[OutputPrecision]]
]
SetAttributes[Error,Listable]


(* ::Text:: *)
(*Hexagonal lattice matrix.*)


AHex={{1,1/2},{0,Sqrt[3]/2}};


(* ::Text:: *)
(*Polynomial for upper Crandall derivatives.*)


polyp[y_,\[Alpha]_,\[Beta]_]:=Apply[Times,(-\[Pi])^(\[Alpha]-\[Beta])Binomial[\[Alpha],\[Beta]](\[Alpha]-\[Beta])!/(\[Alpha]-2\[Beta])!(2y)^(\[Alpha]-2\[Beta])]


(* ::Subsection::Closed:: *)
(*Compiled Direct Sum Representation*)


(* ::Text:: *)
(*Define the reference function in terms of a compiled direct sum representation.*)


Sum2DCompiled=Compile[{
		{\[Nu],_Real},{a11,_Real},{a12,_Real},{a21,_Real},{a22,_Real},{x1,_Real},{x2,_Real},{y1,_Real},{y2,_Real},{\[Alpha]1,_Integer},{\[Alpha]2,_Integer},{lim,_Integer}
	},Module[
		{res=0.+I*0.,c=0.+I*0.,y=0.+I*0.,t=0.+I*0.,term,diff1,diff2,l1,l2}
		,Do[
			l1=a11*z1+a12*z2
			;l2=a21*z1+a22*z2
			;diff1=l1-x1
			;diff2=l2-x2
			;term=Piecewise[{
				{0.,diff1==0&&diff2==0}
				,{Piecewise[{{1,\[Alpha]1==0},{(diff1)^\[Alpha]1,True}}]Piecewise[{{1,\[Alpha]2==0},{(diff2)^\[Alpha]2,True}}]Exp[-2*Pi*I(y1*l1+y2*l2)]/(diff1^2+diff2^2)^(\[Nu]/2),True}
			}]
			(*Kahan summation step*)
			;y=term-c ;t=res+y ;c=(t-res)-y ;res=t
			,{z1,-lim,lim}
			,{z2,-lim,lim}
		]
		;res
	]
];


(* ::Text:: *)
(*Access the direct sum.*)


Options[EpsteinZetaAnisoDirect]={SumLim->1000};
EpsteinZetaAnisoDirect[
	\[Nu]_
	,a_/;Dimensions[a]=={2,2}
	,x_/;Length[x]==2
	,y_/;Length[y]==2
	,\[Alpha]_/;Element[\[Alpha],NonNegativeIntegers]&&Length[\[Alpha]]==2
	,OptionsPattern[]
]:=Sum2DCompiled[N[\[Nu]],N[a[[1,1]]],N[a[[1,2]]],N[a[[2,1]]],N[a[[2,2]]],N[x[[1]]],N[x[[2]]],N[y[[1]]],N[y[[2]]],\[Alpha][[1]],\[Alpha][[2]],OptionValue[SumLim]];


(* ::Subsection::Closed:: *)
(*Harmonic Polynomials*)


(* ::Text:: *)
(*Define harmonic polynomials, so that Total@hTab = h_{k,\[Alpha]}(y).*)


q[k_,n_,m_,d_]:=Product[2n+d-2-4k-2j,{j,m+1,Floor[n/2]-k}]
\[Omega][k_,n_,d_]:=2^(-k)/q[k,n,0,d]/Product[2n+d-2k-2j,{j,1,k}]
c[k_,\[Alpha]_,\[Gamma]_,d_]:=Total@Table[
	(-2)^(-Total[\[Beta]])Total[\[Beta]]! Binomial[Total[\[Beta]]+k,k]q[k,Total[\[Alpha]],Total[\[Beta]],d]Apply[Times,FactorialPower[\[Gamma]-\[Beta],\[Alpha]+\[Beta]-\[Gamma]]/\[Beta]!Binomial[\[Alpha],\[Alpha]+\[Beta]-\[Gamma]]]
	,{\[Beta],Tuples[Range[MapThread[Max,{0\[Alpha],\[Gamma]-\[Alpha]}],\[Gamma]-\[Alpha]/2]]}
]
hTab[k_,\[Alpha]_][y_]:=\[Omega][k,Total[\[Alpha]],Length[\[Alpha]]]Table[
	c[k,\[Alpha],\[Gamma],Length[\[Alpha]]]Apply[Times,y^(2\[Gamma]-\[Alpha])]
	,{\[Gamma],Select[Tuples[Range[Ceiling[\[Alpha]/2],0\[Alpha]+Total[\[Alpha]]-k]],Total[#]==Total[\[Alpha]]-k&]}
]


(* ::Text:: *)
(*Relative backward error, condition number and relative error of one benchmark row {d,k,y,\[Alpha],...,hN,tpre,teval}. Rows where both the exact and the computed values vanish are dropped.*)


hErr[dat_]:=Module[{d=dat[[1]],tab,hExact,hAbs,hN}
	,tab=hTab[dat[[2]],dat[[3+d;;2+2d]]][Rationalize[dat[[3;;2+d]],0]]
	;hExact=Total@tab
	;hAbs=Total@Abs@tab
	;hN=SetPrecision[dat[[-3]],100]
	;If[PossibleZeroQ@hExact
		,If[PossibleZeroQ@hN,Nothing,{Infinity,Infinity,Infinity}]
		,N@{Abs[hN-hExact]/hAbs,hAbs/Abs[hExact],Abs[hN-hExact]/Abs[hExact]}
	]
]


(* ::Section:: *)
(*Figure 1*)


(* ::Subsection::Closed:: *)
(*Definitions*)


(* ::Text:: *)
(*Colors for phases of complex values.*)


ArgColor[arg_]=Hue[(arg+Pi)/(2 Pi)];


(* ::Text:: *)
(*Fix parameters for the (regularized) anisotropic Epstein zeta functions.*)


zfunc[y1_,y2_]:=EpsteinZetaAniso[3.5,IdentityMatrix@2,{1,1}/3,{y1,y2},{0,2}]
zfuncreg[y1_,y2_]:=EpsteinZetaAnisoReg[3.5,IdentityMatrix@2,{1,1}/3,{y1,y2},{0,2}]


(* ::Subsection::Closed:: *)
(*Computation*)


(* ::Text:: *)
(*Generate anisotropic Epstein zeta and regularized anisotropic Epstein zeta values for fixed \[Nu], identity a, all-even \[alpha], all-equal x, and y in a 2D grid [-0.7, 0.7] x [-0.7, 0.7].*)


offset=0.01; xmin=-0.7-offset; xmax=0.7+offset; ymin=-0.7-offset; ymax=0.7+offset;
stepsize=0.07; (* set to 0.007 to match the resolution of the manuscript *)
zanitab=Flatten[ParallelTable[{x,y,Abs[zfunc[x,y]]},{x,xmin,xmax,stepsize},{y,ymin,ymax,stepsize}],1];
stepsize=0.07; (* the regularized function is smooth, so the coarse stepsize suffices *)
zaniregtab=Flatten[ParallelTable[{x,y,Abs[zfuncreg[x,y]]},{x,xmin,xmax,stepsize},{y,ymin,ymax,stepsize}],1];


(* ::Subsection::Closed:: *)
(*Visualization*)


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


(* ::Subsection:: *)
(*Result*)


(* ::Text:: *)
(*Complex phase plots of the 2D anisotropic Epstein zeta function (left) and its regularized version (right) for \[Nu] = 7/2, \[CapitalLambda] = \:2124\.b2 and x =(1/3,1/3), where the height of the surface encodes magnitude and the color encodes the phase.*)


plt


(* ::Section:: *)
(*Table 1*)


(* ::Subsection::Closed:: *)
(*Definitions*)


(* ::Text:: *)
(*Set the reference associated with applying the poly-Laplacian.*)


EpsteinZetaAnisoPolyLaplaceAnalytic[\[Nu]_,a_,x_,y_,n_]:=EpsteinZeta[\[Nu]-2n,a,x,y]


(* ::Text:: *)
(*Rebuild the poly-Laplacian through sums of anisotropic Epstein zeta functions.*)


EpsteinZetaAnisoPolyLaplace[\[Nu]_,a_,x_,y_,n_]:=n!Total@Table[
	1/Apply[Times,\[Alpha]!]EpsteinZetaAniso[\[Nu],a,x,y,2\[Alpha]]
	,{\[Alpha],Select[Tuples[Range[0,n],Length[a]],Total[#]==n&]}
]


(* ::Subsection::Closed:: *)
(*Computation*)


(* ::Text:: *)
(*Generate test values.*)


vectorStepSize[d_]:=ConstantArray[1/2,d] (* set to ConstantArray[1/5,d] to match the resolution of the manuscript *)
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
(*Visualization*)


tbl=TableForm[
	Map[ScientificForm[#,3]&,errs,{2}]
	,TableHeadings->{
		MaTeX["d="<>ToString[#]]&/@{3,2,1}
		,MaTeX["n="<>ToString[#]]&/@Range[6]
	}
]


(* ::Subsection:: *)
(*Result*)


(* ::Text:: *)
(*Maximum error of sums of anisotropic Epstein zeta functions associated with applying the poly-Laplacian \[CapitalDelta]\:207f over a grid of values in \[Nu], x and y.*)


tbl


(* ::Section:: *)
(*Table 2*)


(* ::Subsection::Closed:: *)
(*Definitions*)


(* ::Text:: *)
(*Set the reference for non-mixed second-order derivatives for the square lattice and all-equal vector arguments.*)


EpsteinZetaAnisoAnalytic[
	\[Nu]_
	,a_/;PossibleZeroQ[Norm[a-IdentityMatrix[Length[a]]]]
	,x_/;PossibleZeroQ[Norm[x-ConstantArray[First[x],Length[x]]]]
	,y_/;PossibleZeroQ[Norm[y-ConstantArray[First[y],Length[y]]]]
	,\[Alpha]_/;AllTrue[\[Alpha],EvenQ]&&PossibleZeroQ[2-Total[\[Alpha]]]
]:=EpsteinZeta[\[Nu]-Total[\[Alpha]],a,x,y]/Length[a]


(* ::Subsection::Closed:: *)
(*Computation*)


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
(*Visualization*)


tbl=TableForm[
	Map[ScientificForm[#,3]&,Transpose[errs],{2}]
	,TableHeadings->{
		{MaTeX["E^{\\rm med}"],MaTeX["E^{\\rm max}"]}
		,MaTeX["d="<>ToString[#]]&/@Reverse@Range[6]
	}
]


(* ::Subsection:: *)
(*Result*)


(* ::Text:: *)
(*Median and maximum error of anisotropic Epstein zeta functions corresponding to non-mixed second-order derivatives for the square lattice and all-equal vector arguments over a grid of values in \[Nu], x and y.*)


tbl


(* ::Section:: *)
(*Figure 2*)


(* ::Subsection::Closed:: *)
(*Definitions*)


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


(* ::Subsection::Closed:: *)
(*Computation*)


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
(*Visualization*)


plt=ListLogPlot[
	errs
	,PlotTheme->"Scientific"
	,PlotLabel->MaTeX["Z_{\\mathbb Z^8,\\nu,\\bm 1}(-\\bm 1/2,\\bm 1/2)"]
	,FrameLabel->{MaTeX["\\nu"],MaTeX["E"]}
	,PlotStyle->{PointSize[0.015],Blue}
	,PlotRange->{{-11-1/3,11+1/3},{2 10^(-20),5 10^(-13)}}
	,GridLines->{Range[-10,10,5],10^Range[-21,-12,2]}
	,AspectRatio->.5
	,ImageSize->400
]


(* ::Subsection:: *)
(*Result*)


(* ::Text:: *)
(*Error for a special case of the eight-dimensional anisotropic Epstein zeta function for the square lattice and all-one anisotropy, in which case it is given by the Ramanujan tau L-function.*)


plt


(* ::Section:: *)
(*Figure 3*)


(* ::Subsection::Closed:: *)
(*Definitions*)


(* ::Text:: *)
(*Set test parameters.*)


lim0=1000;
a0={IdentityMatrix[2],IdentityMatrix[2],IdentityMatrix[2],AHex,AHex,AHex};
x0={{0,0},{0,1/2},{0,1/3},{0,0},{0,1/2},{0,1/3}};
y0={{0,0},{1/3,1/4},{1/5,1/6},{0,0},{1/3,1/4},{1/5,1/6}};
\[Alpha]0={{2,0},{1,1},{0,2},{2,0},{1,1},{0,2}};


(* ::Subsection::Closed:: *)
(*Computation*)


(* ::Text:: *)
(*Generate test values.*)


stepsize = 1/4; (* set to 1/50 to match the resolution of the manuscript *)
errs=Table[
	ParallelTable[
		{\[Nu]
		,Error[
			EpsteinZetaAnisoDirect[\[Nu],a0[[i]],x0[[i]],y0[[i]],\[Alpha]0[[i]]]
			,EpsteinZetaAniso[\[Nu],a0[[i]],x0[[i]],y0[[i]],\[Alpha]0[[i]]]
		]
		}
		,{\[Nu],9-1/10,13+1/10,stepsize}
	],{i,6}
]


(* ::Subsection::Closed:: *)
(*Visualization*)


(* ::Text:: *)
(*Set labels depicting the test values.*)


labels={
	"Z_{\\mathbb{Z}^2,\\nu,(2,0)}(\\bm 0,\\bm 0)"
	,"Z_{\\mathbb{Z}^2,\\nu,(1,1)}((0,{\\textstyle\\frac{1}{2}}),({\\textstyle\\frac{1}{3}},{\\textstyle\\frac{1}{4}}))"
	,"Z_{\\mathbb{Z}^2,\\nu,(0,2)}((0,{\\textstyle\\frac{1}{3}}),({\\textstyle\\frac{1}{5}},{\\textstyle\\frac{1}{6}}))"
	,"Z_{\\Lambda_{\\mathrm{hex}},\\nu,(2,0)}(\\bm 0,\\bm 0)"
	,"Z_{\\Lambda_{\\mathrm{hex}},\\nu,(1,1)}((0,{\\textstyle\\frac{1}{2}}),({\\textstyle\\frac{1}{3}},{\\textstyle\\frac{1}{4}}))"
	,"Z_{\\Lambda_{\\mathrm{hex}},\\nu,(0,2)}((0,{\\textstyle\\frac{1}{3}}),({\\textstyle\\frac{1}{5}},{\\textstyle\\frac{1}{6}}))"
}


plt=GraphicsGrid[
	ArrayReshape[
		Table[
			ListLogPlot[
			{First@#,Last@#}&/@errs[[i]]
			,PlotTheme->"Scientific"
			,PlotLabel->MaTeX[labels[[i]]]
			,FrameLabel->{MaTeX@"\\nu",MaTeX@"E"}
			,PlotStyle->{PointSize[0.025],Blue}
			,ImageSize->240
			]
			,{i,6}
		]
		,{2,3}
	]
	,ImageSize->1000
]


(* ::Subsection:: *)
(*Result*)


(* ::Text:: *)
(*Error for low-order anisotropic Epstein zeta functions for large \[Nu], where direct summation gives accurate reference values.*)


plt


(* ::Section:: *)
(*Figure 4*)


(* ::Subsection::Closed:: *)
(*Definitions*)


(* ::Text:: *)
(*Set test parameters.*)


lim0=1000;
a0={IdentityMatrix[2],IdentityMatrix[2],IdentityMatrix[2],AHex,AHex,AHex};
x0={{0,0},{0,1/2},{0,1/3},{0,0},{0,1/2},{0,1/3}};
y0={{0,0},{1/3,1/4},{1/5,1/6},{0,0},{1/3,1/4},{1/5,1/6}};
\[Alpha]0={{8,10},{9,10},{10,9},{8,10},{9,10},{10,9}};


(* ::Subsection::Closed:: *)
(*Computation*)


(* ::Text:: *)
(*Generate test values.*)


stepsize = 1/4; (* set to 1/50 to match the resolution of the manuscript *)
errs=Table[
	ParallelTable[
		{\[Nu]
		,Error[
			EpsteinZetaAnisoDirect[\[Nu],a0[[i]],x0[[i]],y0[[i]],\[Alpha]0[[i]]]
			,EpsteinZetaAniso[\[Nu],a0[[i]],x0[[i]],y0[[i]],\[Alpha]0[[i]]]
		]
		}
		,{\[Nu],25-1/10,29+1/10,stepsize}
	],{i,6}
]


(* ::Subsection::Closed:: *)
(*Visualization*)


(* ::Text:: *)
(*Set labels depicting the test values.*)


labels={
	"Z_{\\mathbb{Z}^2,\\nu,(8,10)}(\\bm 0,\\bm 0)"
	,"Z_{\\mathbb{Z}^2,\\nu,(9,10)}((0,{\\textstyle\\frac{1}{2}}),({\\textstyle\\frac{1}{3}},{\\textstyle\\frac{1}{4}}))"
	,"Z_{\\mathbb{Z}^2,\\nu,(10,9)}((0,{\\textstyle\\frac{1}{3}}),({\\textstyle\\frac{1}{5}},{\\textstyle\\frac{1}{6}}))"
	,"Z_{\\Lambda_{\\mathrm{hex}},\\nu,(8,10)}(\\bm 0,\\bm 0)"
	,"Z_{\\Lambda_{\\mathrm{hex}},\\nu,(9,10)}((0,{\\textstyle\\frac{1}{2}}),({\\textstyle\\frac{1}{3}},{\\textstyle\\frac{1}{4}}))"
	,"Z_{\\Lambda_{\\mathrm{hex}},\\nu,(10,9)}((0,{\\textstyle\\frac{1}{3}}),({\\textstyle\\frac{1}{5}},{\\textstyle\\frac{1}{6}}))"
};


plt=GraphicsGrid[
	ArrayReshape[
		Table[
			ListLogPlot[
			{First@#,Last@#}&/@errs[[i]]
			,PlotTheme->"Scientific"
			,PlotLabel->MaTeX[labels[[i]]]
			,FrameLabel->{MaTeX@"\\nu",MaTeX@"E"}
			,PlotStyle->{PointSize[0.025],Blue}
			,ImageSize->240
			]
			,{i,6}
		]
		,{2,3}
	]
	,ImageSize->1000
]


(* ::Subsection:: *)
(*Result*)


(* ::Text:: *)
(*Error for high-order anisotropic Epstein zeta functions for large \[Nu], where direct summation gives accurate reference values.*)


plt


(* ::Section:: *)
(*Figure 5*)


(* ::Subsection::Closed:: *)
(*Definitions*)


(* ::Text:: *)
(*Define the (derivatives of the) upper Crandall functions.*)


GUpper[\[Nu]_,y_]=Piecewise[{
	{-2/\[Nu],Norm[y]<10^(-64)}
	,{ExpIntegralE[1-\[Nu]/2,Pi y . y],Norm[y]<6}
	},0
];
GUpperDer[\[Nu]_,y_,\[Alpha]_]:=Total[Flatten[Table[
	Module[
		{pol}
		,pol=polyp[y,\[Alpha],\[Beta]]
		;If[Abs[pol]<10^-60
			,0
			,pol GUpper[\[Nu]+2Total[\[Alpha]-\[Beta]],y]
		]
	]
	,{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]/2]]}
]]]


(* ::Text:: *)
(*Define sums of the derivative Crandall formula in real and reciprocal space.*)


SumRealDer[\[Nu]_,A_,x_,y_,\[Alpha]_,sumLim_]:=Total@Table[
	Apply[Times, z^\[Alpha]]GUpper[\[Nu],z]Exp[-2Pi I y . (z+x)]
	,{z,(A . #-x)&/@Tuples[Range[-sumLim,sumLim],Length[A]]}
]
SumFourierDer[\[Nu]_,A_,x_,y_,\[Alpha]_,sumLim_]:=Total@Table[
	GUpperDer[Length[\[Alpha]]-\[Nu],k,\[Alpha]]Exp[-2Pi I x . k]
	,{k,(Inverse[Transpose[A]] . #+y)&/@Tuples[Range[-sumLim,sumLim],Length[A]]}
]


(* ::Text:: *)
(*Prototype of the anisotropic Epstein zeta function based on the Crandall derivative representation.*)


Options[epsteinZetaAniso]={InputPrecision->100,SumLim->6};
epsteinZetaAniso[\[Nu]Ext_?NumericQ,AExt_/;MatrixQ[AExt,NumericQ],xExt_/;VectorQ[xExt,NumericQ],yExt_/;VectorQ[yExt,NumericQ],\[Alpha]_/;VectorQ[\[Alpha],IntegerQ],OptionsPattern[]]:=Module[
	{EPS,prec,sumLim,\[Nu]N,d,AN,xN,yN,\[Mu],sumReal,sumReci,res}
	,EPS=10^-50
	;prec=OptionValue[InputPrecision]
	;sumLim=OptionValue[SumLim]
	;d=Length[\[Alpha]]
	;\[Nu]N=N[\[Nu]Ext,prec]
	;AN=N[AExt,prec]
	;\[Mu]=Abs[Det[AN]]^(1/Length[AN])
	;AN=1/\[Mu] AN
	;xN=1/\[Mu] N[xExt,prec]
	;yN=\[Mu] N[yExt,prec]
	;Do[
		If[
			AExt==IdentityMatrix[d]&&OddQ[\[Alpha][[i]]]&&(PossibleZeroQ[xN[[i]]]&&PossibleZeroQ[yN[[i]]])
			,Return[0]
		]
		,{i,Length[\[Alpha]]}
	]
	;If[Abs[\[Nu]N]<EPS&&Total[\[Alpha]]>0
		,res = 0
		,sumReal=SetPrecision[SumRealDer[\[Nu]N,AN,xN,yN,\[Alpha],sumLim],prec]
		;sumReci=SetPrecision[SumFourierDer[\[Nu]N,AN,xN,yN,\[Alpha],sumLim],prec]
		;res=SetPrecision[\[Mu]^(Total[\[Alpha]]-\[Nu]N)\[Pi]^(\[Nu]N/2)/Gamma[\[Nu]N/2](sumReal+N[-2 Pi I]^(-Total[\[Alpha]])/Abs[Det[AN]]sumReci),prec]
	]
	;Return[res]
]


(* ::Text:: *)
(*Set test parameters.*)


\[Nu]0=5/2;
x0:={ConstantArray[0,dim],UnitVector[dim,1]/2,UnitVector[dim,1]/2};
y0:={ConstantArray[0,dim],ConstantArray[0,dim],UnitVector[dim,1]/3};


(* ::Subsection::Closed:: *)
(*Computation*)


(* ::Text:: *)
(*Generate test values.*)


stabnrange=Range[0,30,3]; (* set to Range[0,61,1] to match the resolution of the manuscript *)
errs=Table[
	ParallelTable[
		Module[{\[Nu],a,x,y,ref}
			,\[Nu]=5/2
			;a=IdentityMatrix[d]
			;x={
				ConstantArray[0,d]
				,UnitVector[d,1]/2
				,UnitVector[d,1]/2
			}[[i]]
			;y={
				ConstantArray[0,d]
				,ConstantArray[0,d]
				,UnitVector[d,1]/3
			}[[i]]
			;\[Alpha]=n UnitVector[d,1]
			;ref=epsteinZetaAniso[\[Nu],a,x,y,\[Alpha],InputPrecision->100]
			;{
				Error[ref,EpsteinZetaAniso[\[Nu],a,x,y,\[Alpha]]]
				,Error[ref,epsteinZetaAniso[\[Nu],a,x,y,\[Alpha],InputPrecision->MachinePrecision]]
			}
		]
		,{n,stabnrange}
	]
	,{i,3}
	,{d,3}
];
Do[
	stabTabC[i]=Transpose@{stabnrange,#}&/@errs[[i,All,All,1]]
	;stabTabM[i]=Transpose@{stabnrange,#}&/@errs[[i,All,All,2]]
	,{i,3}
]


(* ::Text:: *)
(*Prepare fits of the derivative Crandall prototype errors.*)


fitStab=Table[
	Fit[{#[[1]],Log[#[[2]]]}&/@Select[stabTabM[i][[d]],#[[2]]>0&],{1,n},n]
	,{i,3}
	,{d,3}
];


(* ::Subsection::Closed:: *)
(*Visualization*)


(* ::Text:: *)
(*Set labels depicting the test values.*)


labels={
	"Z_{\\mathbb Z^d,5/2,n\\bm e_1}(\\bm 0,\\bm 0)"
	,"Z_{\\mathbb Z^d,5/2,n\\bm e_1}({\\textstyle\\frac{1}{2}}\\bm e_1,\\bm 0)"
	,"Z_{\\mathbb Z^d,5/2,n\\bm e_1}({\\textstyle\\frac{1}{2}}\\bm e_1,{\\textstyle\\frac{1}{3}}\\bm e_1)"
};


(* ::Text:: *)
(*Set legend depicting the colors of the dimensions.*)


legend=PointLegend[ColorData[97]/@Range[3],Table[MaTeX["d="<>ToString[d]],{d,3}]];


plt=GraphicsGrid[{
	Append[
		Table[Show[
				ListLogPlot[
					stabTabM[i]
					,PlotTheme->"Scientific"
					,PlotStyle->ColorData[97]/@Range[3]
					,PlotLabel->MaTeX[labels[[i]]]
					,FrameLabel->{MaTeX@"n",MaTeX@"E"}
					,PlotRange->All
					,Epilog->Inset[
						Column[Table[MaTeX["\\rho_{d="<>ToString[d]<>"}="<>ToString@NumberForm[Exp@Coefficient[fitStab[[i,d]],n],{10,2}]],{d,3}]]
						,Scaled[{0.05,0.95}]
						,{Left,Top}
					]
					,ImageSize->240
				]
				,LogPlot[
					Evaluate@Exp@fitStab[[i]]
					,{n,0,61}
					,PlotStyle->(Directive[Opacity[0.5],ColorData[97][#]]&/@Range[3])
				]
			]
			,{i,3}
		]
		,legend
	]
	,Append[
		Table[
			ListLogPlot[
				Map[{#[[1]],Max[10^-17,#[[2]]]}&,stabTabC[i],{2}]
				,PlotTheme->"Scientific"
				,PlotStyle->ColorData[97]/@Range[3]
				,PlotLabel->MaTeX[labels[[i]]]
				,FrameLabel->{MaTeX@"n",MaTeX@"E"}
				,PlotRange->All
				,ImageSize->240
			]
			,{i,3}
		]
		,legend
	]
	}
	,ImageSize->1200
]


(* ::Subsection:: *)
(*Result*)


(* ::Text:: *)
(*Error of a prototype of the derivative Crandall representation for the anisotropic Epstein zeta function (first row) and of the stable harmonic-decomposition Crandall representation implemented in EpsteinLib (second row).*)


plt


(* ::Section:: *)
(*Figure 6*)


(* ::Subsection::Closed:: *)
(*Definitions*)


(* ::Text:: *)
(*Set test parameters.*)


\[Nu]0=1/2;
taylorNMax=12; (* set to 61 to match the resolution of the manuscript *)
stepsize=1/4; (* set to 1/20 to match the resolution of the manuscript *)


(* ::Subsection::Closed:: *)
(*Computation*)


(* ::Text:: *)
(*Generate test values where, for each n, the Taylor coefficients of all-even \[Alpha] with |\[Alpha]|=n are added to the partial sum.*)


errs=Table[
	With[{a=IdentityMatrix[d],x=ConstantArray[0,d]}
		,With[{coeffs=Table[
				ParallelTable[
					{\[Alpha],(-2Pi I)^n/Apply[Times,\[Alpha]!]EpsteinZetaAniso[\[Nu]0,a,x,ConstantArray[0,d],\[Alpha]]}
					,{\[Alpha],Select[Tuples[Range[0,n,2],d],Total[#]==n&]}
				]
				,{n,0,taylorNMax}
			]}
			,Transpose@{
				Range[0,taylorNMax]
				,Max/@Transpose@ParallelTable[
					Module[{ref=EpsteinZetaReg[\[Nu]0,a,x,y],acc=0}
						,Table[
							acc+=Total[Last[#]Apply[Times,y^First[#]]&/@coeffs[[n+1]]]
							;Error[ref,acc]
							,{n,0,taylorNMax}
						]
					]
					,{y,N@Select[Tuples[Range[-1/2,1/2,stepsize],d],Norm[#]<=1/2&]}
				]
			}
		]
	]
	,{d,3}
];


(* ::Subsection::Closed:: *)
(*Visualization*)


plt=ListLogPlot[
	errs
	,PlotTheme->"Scientific"
	,PlotLabel->MaTeX["Z^{\\rm reg}_{\\mathbb Z^d,1/2}(\\bm 0,\\bm y)"]
	,FrameLabel->{MaTeX@"n",MaTeX["E^{\\rm max}"]}
	,PlotLegends->Table[MaTeX["d="<>ToString[d]],{d,3}]
	,PlotRange->All
	,ImageSize->400
]


(* ::Subsection:: *)
(*Result*)


(* ::Text:: *)
(*Maximum error of the reconstruction of the regularized Epstein zeta function through its Taylor series, whose coefficients are given in terms of anisotropic Epstein zeta functions, of order n around y = 0.*)


plt


(* ::Section:: *)
(*Table 3*)


(* ::Subsection::Closed:: *)
(*Computation*)


(* ::Text:: *)
(*Import the benchmark data for \[Alpha]=(n,2n), where the C routine forces double arithmetic.*)


dataStab=Table[(* remove the dummy YINCSTAB2D (1. / 4) line in benchmark_harmonics.c to match the resolution of the manuscript *)
	Import[NotebookDirectory[]<>"../../tests/csv/benchmark_harmonic_stab_"<>ToString[n]<>"_"<>ToString[2n]<>"_2D.csv"]
	,{n,0,10}
];


(* ::Text:: *)
(*Generate test values.*)


errsStab=Table[ParallelTable[hErr[dat],{dat,dataStab[[j]]}],{j,11}];


(* ::Subsection::Closed:: *)
(*Visualization*)


tbl=TableForm[
	Map[ScientificForm[N@#,2]&
		,Table[{
			Max@dataStab[[j,All,-2]]
			,Max@dataStab[[j,All,-1]]
			,Median@errsStab[[j,All,2]]
			,Max@errsStab[[j,All,2]]
			,Max@errsStab[[j,All,1]]
		},{j,11}]
		,{2}
	]
	,TableHeadings->{
		MaTeX["\\bm\\alpha=("<>ToString[#]<>","<>ToString[2#]<>")"]&/@Range[0,10]
		,MaTeX/@{"t_{\\rm pre}","t_{\\rm eval}^{\\rm max}","\\kappa_{\\rm med}","\\kappa_{\\rm max}","E_{\\rm bwd}^{\\rm max}"}
	}
];


(* ::Subsection:: *)
(*Result*)


(* ::Text:: *)
(*Evaluation times, condition numbers and backward errors of the double algorithm for two-dimensional harmonic polynomials compared with an exact implementation over a grid of values.*)


tbl


(* ::Section:: *)
(*Table 4*)


(* ::Subsection::Closed:: *)
(*Computation*)


(* ::Text:: *)
(*Import the benchmark data for all-equal \[Alpha] with |\[Alpha]|=24 precomputed by the adaptive algorithm.*)


dataHigh=Table[(* remove the dummy YINC[DIM_MAX] line in benchmark_harmonics.c to match the resolution of the manuscript *)
	Import[NotebookDirectory[]<>"../../tests/csv/benchmark_harmonic_"<>ToString[d]<>"D.csv"]
	,{d,4}
];


(* ::Text:: *)
(*Generate test values.*)


sampleStep=10; (* set to 1 to match the resolution of the manuscript *)
dataHigh=Table[
	Import[NotebookDirectory[]<>"../../tests/csv/benchmark_harmonic_"<>ToString[d]<>"D.csv"][[;;;;sampleStep]]
	,{d,4}
];
errsHigh=Table[ParallelTable[hErr[dat],{dat,dataHigh[[d]]}],{d,4}];


(* ::Subsection::Closed:: *)
(*Visualization*)


tbl=TableForm[
	Map[ScientificForm[N@#,2]&
		,Table[{
			Max@dataHigh[[d,All,-2]]
			,Max@dataHigh[[d,All,-1]]
			,100Count[errsHigh[[d,All,2]],_?(#>8&)]/Length[errsHigh[[d]]]
			,Max@errsHigh[[d,All,2]]
			,Max@errsHigh[[d,All,3]]
		},{d,4}]
		,{2}
	]
	,TableHeadings->{
		MaTeX["d="<>ToString[#]]&/@Range[4]
		,MaTeX/@{"t_{\\rm pre}","t_{\\rm eval}^{\\rm max}","\\%(\\kappa_h>8)","\\kappa_{\\rm max}","E_{\\rm rel}^{\\rm max}"}
	}
];


(* ::Subsection:: *)
(*Result*)


(* ::Text:: *)
(*Evaluation times, condition numbers and relative errors of the adaptive algorithm for one- to four-dimensional harmonic polynomials compared with an exact implementation over a grid of values.*)


tbl


(* ::Section:: *)
(*Figure 7*)


(* ::Subsection::Closed:: *)
(*Computation*)


(* ::Text:: *)
(*Import the runtimes of anisotropic Epstein zeta functions in EpsteinLib.*)


timing=Table[
	Import[NotebookDirectory[]<>"../../tests/csv/epsteinZetaAniso_timing_"<>ToString[d]<>"D.csv"]
	,{d,6}
];


(* ::Text:: *)
(*Set relative runtimes and obtain their fits.*)


relTiming=Table[{#[[1]],#[[2]]/timing[[d,1,2]]}&/@timing[[d,1;;21]],{d,3}];
fitTime=Table[Fit[relTiming[[d,1;;If[d>1,11,21]]],{1,n},n],{d,3}];


(* ::Subsection::Closed:: *)
(*Visualization*)


plt=GraphicsRow[{
	ListLogPlot[
		timing
		,PlotTheme->"Scientific"
		,PlotStyle->ColorData[97]/@Range[6]
		,FrameLabel->{MaTeX["|\\bm\\alpha|"],MaTeX["\\text{runtime }[s]"]}
		,PlotRange->{{-0.5,10.5},Automatic}
		,ImageSize->300
	]
	,Show[
		ListPlot[
			relTiming
			,PlotTheme->"Scientific"
			,PlotStyle->ColorData[97]/@Range[3]
			,FrameLabel->{MaTeX["|\\bm\\alpha|"],MaTeX["\\text{relative runtime}"]}
			,PlotRange->All
			,Epilog->Inset[
				Column[Table[MaTeX["a_{d="<>ToString[d]<>"}\\approx "<>ToString@NumberForm[Coefficient[fitTime[[d]],n],{10,2}]],{d,3}]]
				,Scaled[{0.05,0.95}]
				,{Left,Top}
			]
			,ImageSize->300
		]
		,Plot[
			Evaluate@fitTime
			,{n,0,20}
			,PlotStyle->(Directive[Opacity[0.5],ColorData[97][#]]&/@Range[3])
		]
	]
	,PointLegend[ColorData[97]/@Range[6],Table[MaTeX["d="<>ToString[d]],{d,6}]]
	}
	,ImageSize->800
]


(* ::Subsection:: *)
(*Result*)


(* ::Text:: *)
(*Evaluation times (left) and relative evaluation times (right) of anisotropic Epstein zeta functions as a function of the order of the anisotropy.*)


plt


(* ::Section:: *)
(*Figure 8*)


(* ::Subsection::Closed:: *)
(*Definitions*)


(* ::Text:: *)
(*fcc lattice matrix.*)


AFcc={{1/Sqrt[2],1/Sqrt[2],0},{1/Sqrt[2],0,1/Sqrt[2]},{0,1/Sqrt[2],1/Sqrt[2]}};


(* ::Text:: *)
(*Set test parameters.*)


\[Nu]0=5/2;
\[Sigma]0=10;
a0=AFcc;
lim0=35;
stepsize=1/2; (* set to 1/10 to match the resolution of the manuscript *)
x0range=Range[-1/10,10+1/10,stepsize];


(* ::Text:: *)
(*Define the (derivatives of the) Gaussian.*)


Gaussian3D[y_,\[Sigma]_]:=Exp[-Pi y . y/\[Sigma]^2];
GaussianDerivatives3D[y_,\[Sigma]_,\[Alpha]_]:=Gaussian3D[y,\[Sigma]]Total[Table[
	(1/\[Sigma]^2)^Total[\[Alpha]-\[Beta]]polyp[y,\[Alpha],\[Beta]]
	,{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]/2]]}
]]


(* ::Text:: *)
(*Define the reference direct sum.*)


Sum3D[\[Nu]_,a_,x_,\[Sigma]_,lim_]:=Module[{pts,d2,g}
	,pts=Transpose[N[Tuples[Range[-lim,lim],3]] . Transpose[N[a]]]
	;d2=Total[(pts-N[x])^2]
	;g=Exp[-Pi Total[pts^2]/\[Sigma]^2]
	;Total[g Unitize[d2]/(d2+1-Unitize[d2])^(\[Nu]/2),Method->"CompensatedSummation"]
]


(* ::Text:: *)
(*Integral and SEM for the singular Euler-Maclaurin expansion.*)


Integral3D[\[Nu]_,a_,x_,\[Sigma]_]:=2\[Pi]^((\[Nu]-1)/2)\[Sigma]^(3-\[Nu])Gamma[(3-\[Nu])/2]Hypergeometric1F1[\[Nu]/2,3/2,-\[Pi] x . x/\[Sigma]^2]/Abs[Det[a]]
SEM3D[\[Nu]_,a_,x_,\[Sigma]_,\[ScriptL]_]:=Integral3D[\[Nu],a,x,\[Sigma]]+Total[Table[
	If[PossibleZeroQ@Norm@x&&AnyTrue[\[Beta],OddQ]
		,0
		,Re@EpsteinZetaAniso[\[Nu],a,x,ConstantArray[0,3],\[Beta]]GaussianDerivatives3D[x,\[Sigma],\[Beta]]/Apply[Times,\[Beta]!]
	]
	,{\[Beta],Select[Tuples[Range[0,2\[ScriptL]],3],Total[#]<=2\[ScriptL]&]}
]]


(* ::Subsection::Closed:: *)
(*Computation*)


(* ::Text:: *)
(*Generate test values.*)


Sum3DTab=ParallelTable[{x,Sum3D[\[Nu]0,a0,{x,0,1/2},\[Sigma]0,lim0]},{x,x0range}];
Int3DTab=Table[{x,Integral3D[\[Nu]0,a0,{x,0,1/2},\[Sigma]0]},{x,First@x0range-1,Last@x0range+1,1/10}];
SEM3DTab=ParallelTable[{x,SEM3D[\[Nu]0,a0,{x,0,1/2},\[Sigma]0,0]},{x,First@x0range-1,Last@x0range+1,1/20}];
errs=Table[
	Transpose@{x0range,Error[Sum3DTab[[All,2]],ParallelTable[SEM3D[\[Nu]0,a0,{x,0,1/2},\[Sigma]0,\[ScriptL]],{x,x0range}]]}
	,{\[ScriptL],0,6}
];


(* ::Subsection::Closed:: *)
(*Visualization*)


plt=GraphicsRow[{
	Show[
		ListLinePlot[
			{SEM3DTab,Int3DTab}
			,PlotTheme->"Scientific"
			,PlotLegends->Placed[{
				MaTeX["\\mathcal D^{(0)}_{\\Lambda_{\\rm fcc}}g(\\bm x)+\\textstyle\\frac{1}{V_{\\Lambda_{\\rm fcc}}}\\int f_{\\bm x}"]
				,MaTeX["\\textstyle\\frac{1}{V_{\\Lambda_{\\rm fcc}}}\\int f_{\\bm x}"]
				},{0.7,0.8}
			]
			,PlotLabel->MaTeX["\\bm x=(x,0,1/2)"]
			,PlotRange->{{-.35,10.35},{-5,105}}
			,PlotStyle->{Black,Blue}
			,FrameLabel->{MaTeX@"x",MaTeX@"f_{\\bm x}(\\bm z)=\\frac{e^{-\\pi|\\bm z|^2/100}}{|\\bm z-\\bm x|^{5/2}}"}
			,ImageSize->360
		]
		,ListPlot[
			Sum3DTab
			,PlotStyle->Red
			,PlotLegends->Placed[{MaTeX["\\textstyle\\sum_{\\Lambda_{\\rm fcc}}f_{\\bm x}"]},{0.2,0.2}]
		]
	]
	,ListLogPlot[
		errs
		,PlotTheme->"Scientific"
		,PlotLabel->MaTeX["\\mathcal D^{(\\ell)}_{\\Lambda_{\\rm fcc}}g(x,0,1/2)"]
		,FrameLabel->{MaTeX@"x",MaTeX@"E"}
		,PlotLegends->Table[MaTeX["\\ell="<>ToString[\[ScriptL]]],{\[ScriptL],0,6}]
		,PlotRange->All
		,ImageSize->400
	]
	}
	,ImageSize->1000
]


(* ::Subsection:: *)
(*Result*)


(* ::Text:: *)
(*Singular Euler-Maclaurin expansion for the Gaussian and the fcc lattice (left) and its error for different orders (right).*)


plt
