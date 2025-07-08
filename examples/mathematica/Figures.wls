#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan.busse@dlr.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*EpsteinLib Figures*)


(* ::Text:: *)
(*This notebook generates the results for all figures presented in the paper "Computation and Properties of the Epstein Zeta Function With High-Performance Implementation in EpsteinLib". Evaluation of this script takes less than 3 minutes on a standard laptop equipped with 12 cores.*)


(* ::Subsection:: *)
(*Setup Packages*)


(* ::Text:: *)
(*Please ensure that EpsteinLib is correctly installed and accessible before running the notebook. That is*)
(*	*)
(*	1. 	go to the project: cd C*)
(*	2.	build library: meson setup build*)
(*	3.	compile library: meson compile -C build*)
(*	*)
(*To test the C library, run meson test -C build.*)


(* ::Text:: *)
(*Load the Mathematica wrapper of the Epstein zeta function,  such that it can be used in parallel on all Kernels.*)


SetDirectory[NotebookDirectory[]<>"../../mathematica"];
ParallelEvaluate[Quiet[<<"EpsteinZeta.wl"]];
<<"EpsteinZeta.wl";


(* ::Text:: *)
(*Install MaTeX for LATEX figure annotations. If MaTeX is not installed already un-comment the ResourceFunction line to install it.*)


(*ResourceFunction["MaTeXInstall"][];*)
Needs["MaTeX`"];


(* ::Text:: *)
(*Set options for pretty plots.*)


SetOptions[RegionPlot,TicksStyle->{FontSize->18},LabelStyle->{FontSize->18}];
SetOptions[MaTeX,"Preamble"->{"\\usepackage{xcolor}","\\usepackage{bm}"}];
SetOptions[MaTeX,FontSize->18];


(* ::Section:: *)
(*Figure 1*)


(* ::Text:: *)
(*Generate Epstein zeta and regularized Epstein zeta values for a 2D grid in the range [-0.7, 0.7] x [-0.7, 0.7], avoiding the origin.*)


offset=0.01; xmin=-0.7-offset; xmax=0.7+offset; ymin =-0.7-offset; ymax=0.7+offset; stepsize=0.01;
ztab3D[1]=Flatten[ParallelTable[{x,y,Abs[EpsteinZeta[1/2,{{1}},{y},{x}]]},{x,xmin,0-offset,stepsize},{y,ymin,0-offset,stepsize}],1];
ztab3D[2]=Flatten[ParallelTable[{x,y,Abs[EpsteinZeta[1/2,{{1}},{y},{x}]]},{x,xmin,0-offset,stepsize},{y,0+offset,ymax,stepsize}],1];
ztab3D[3]=Flatten[ParallelTable[{x,y,Abs[EpsteinZeta[1/2,{{1}},{y},{x}]]},{x,0+offset,xmax,stepsize},{y,ymin,0-offset,stepsize}],1];
ztab3D[4]=Flatten[ParallelTable[{x,y,Abs[EpsteinZeta[1/2,{{1}},{y},{x}]]},{x,0+offset,xmax,stepsize},{y,0+offset,ymax,stepsize}],1];
ztab3D[5]=Join[Flatten[ParallelTable[{x,y,Abs[EpsteinZeta[1/2,{{1}},{y},{x}]]},{x,xmin,0-offset, stepsize},{y,{0}}],1],Flatten[ParallelTable[{x,y,Abs[EpsteinZeta[1/2,{{1}},{y},{x}]]},{x,0+offset,xmax, stepsize},{y,{0}}],1]];
ztab3D[6]=Join[Flatten[ParallelTable[{x,y,Abs[EpsteinZeta[1/2,{{1}},{y},{x}]]},{x,{0}},{y,ymin,0-offset,stepsize}],1],Flatten[ParallelTable[{x,y,Abs[EpsteinZeta[1/2,{{1}},{y},{x}]]},{x,{0}},{y,0+offset,ymax,stepsize}],1]];
zregtab3D[1]=Flatten[ParallelTable[{x,y,Abs[EpsteinZetaReg[1/2,{{1}},{y},{x}]]},{y,xmin,0-offset,stepsize},{x,ymin,ymax,stepsize}],1];
zregtab3D[2]=Flatten[ParallelTable[{x,y,Abs[EpsteinZetaReg[1/2,{{1}},{y},{x}]]},{y,0+offset,xmax,stepsize},{x,ymin,ymax,stepsize}],1];


(* ::Text:: *)
(*Create a side-by-side 3D plot of the absolute values of the Epstein zeta function and its regularized version.*)


GraphicsRow[{Show[ListPlot3D[ztab3D[#],AxesLabel->{MaTeX@"y",MaTeX@"x"},ColorFunction->"Rainbow",PlotLabel->MaTeX@"|Z_{\\mathbb{Z},\\:\\! 1/2}\\begin{vmatrix}x\\\\ y\\end{vmatrix}|",PlotRange->{-0.5,8},Mesh->None]&/@Range[6],
Show[ListLinePlot3D[ztab3D[#],PlotStyle->Thickness[0.01],PlotRange->{-0.5,8}]&/@{5,6}]],
Show[ListPlot3D[zregtab3D[#],AxesLabel->{MaTeX@"y",MaTeX@"x"},ColorFunction->"Rainbow",PlotLabel->MaTeX@"|Z_{\\mathbb{Z},\\:\\! 1/2}^{\\mathrm{reg}}\\begin{vmatrix}x\\\\ y\\end{vmatrix}|",PlotRange->{-0.5,8},Mesh->None]&/@Range[2],
Show[ListLinePlot3D[ztab3D[#],PlotStyle->Thickness[0.01],ColorFunction->"Rainbow"]&/@{5,6}]]}]


(* ::Section:: *)
(*Figure 2*)


(* ::Text:: *)
(*Define regions for different computational methods of the Epstein zeta function based on input parameters nu and x.*)


\[Nu]minsmall=-2; \[Nu]minbig=-10; \[Nu]maxsmall=2; \[Nu]maxbig=40; xminsmall=0; xminbig=xminsmall; xmaxsmall=3; xmaxbig=40;
asymptotic=\[Nu]>12.5&&\[Nu]<1/0.3x&&\[Nu]>1/2.35x;
ptaylor=If[x<0.5,\[Nu]>Log[0.5]/Log[0.5x],\[Nu]>x]&&!asymptotic;
contfrac=x>1.5&&!ptaylor&&!asymptotic;
recursive=!contfrac&&\[Nu]<-1/2;
qtaylor=!ptaylor&&!contfrac&&!recursive;
regionssmall={ptaylor,contfrac,qtaylor,recursive};
regionsbig=Join[regionssmall,{asymptotic}];


(* ::Text:: *)
(*Create a side-by-side plot showing the regions where different computational methods are optimal for the upper gamma function, for both large and small scale parameter ranges.*)


GraphicsRow[{Show[RegionPlot[regionsbig,{\[Nu],\[Nu]minbig,\[Nu]maxbig},{x,xminbig,xmaxbig},PlotLegends->SwatchLegend[{Automatic,Automatic,Automatic,Automatic,Opacity[0.8,Gray]},If[TrueQ[#],None,MaTeX["\\text{"<>#<>"}"]]&/@{"P-Taylor","cont. frac.",True,True,"asymptotic"}], Epilog->{Dashed,Thick,Line[{{\[Nu]minsmall,xminsmall},{\[Nu]maxsmall,xminsmall},{\[Nu]maxsmall,xmaxsmall},{\[Nu]minsmall,xmaxsmall},{\[Nu]minsmall,xminsmall}}]}, FrameLabel->{{MaTeX["x"],""},{MaTeX["\\nu"],""}}, ImageSize->Medium, PlotPoints->100],Plot[If[\[Nu]<1/2,2/2^(1/\[Nu]),\[Nu]],{\[Nu],0,\[Nu]maxbig},PlotStyle->{Black,Thick}]],
Show[RegionPlot[regionssmall,{\[Nu],\[Nu]minsmall,\[Nu]maxsmall},{x,xminsmall,xmaxsmall},PlotLegends->(MaTeX["\\text{"<>#<>"}"]&/@{"P-Taylor","cont. frac.","Q-Taylor","recursive"}),FrameLabel->{{MaTeX["x"],""},{MaTeX["\\nu"],""}}, ImageSize->Medium, PlotPoints->100],Plot[If[\[Nu]<1/2,2/2^(1/\[Nu]),\[Nu]],{\[Nu],0,\[Nu]maxsmall}],Plot[If[\[Nu]<1/2,2/2^(1/\[Nu]),\[Nu]],{\[Nu],0,\[Nu]maxsmall},PlotStyle->{Black,Thick}]]},ImageSize->1000]


(* ::Section:: *)
(*Figure 3*)


(* ::Text:: *)
(*Load reference data for the upper incomplete gamma function from CSV files. Each line contains nu, x, Re(Gamma), and Im(Gamma) values.  After compiling the library with meson, reference values can be regenerated by *)
(**)
(*	1. 	got to project: cd C*)
(*	2.	 execute binary: ./build/benchmark/epsteinlib_benchmark_gamma*)


gammaDataBig=Import[SetDirectory[NotebookDirectory[]]<>"/../../benchmark/csv/gamma_big.csv"];
gammaDataSmall=Import[SetDirectory[NotebookDirectory[]]<>"/../../benchmark/csv/gamma_small.csv"];


(* ::Text:: *)
(*Calculate high-precision (25 digits) reference values for the upper incomplete gamma function using Mathematica's built-in Gamma function (takes around 20 seconds on a standard laptop with 12 kernels).*)


gammaMatBig=ParallelTable[SetPrecision[Gamma[Rationalize[arg[[1]]],Rationalize[arg[[2]]]],25],{arg,gammaDataBig}];
gammaMatSmall=ParallelTable[SetPrecision[Gamma[Rationalize[arg[[1]]],Rationalize[arg[[2]]]],25],{arg,gammaDataSmall}];


(* ::Text:: *)
(*Define a function to compute the minimum of absolute and relative errors, then calculate errors between C implementation results and high-precision Mathematica results.*)


SetAttributes[MinErr,Listable];
MinErr[exact_,approx_]:=Min[Abs[exact-approx],Abs[exact-approx]/Abs[exact]];


gammaErrorsBig=Flatten[#]&/@Transpose[{gammaDataBig[[All,1;;2]],MinErr[gammaMatBig,gammaDataBig[[All,3]]]}];
gammaErrorsSmall=Flatten[#]&/@Transpose[{gammaDataSmall[[All,1;;2]],MinErr[gammaMatSmall,gammaDataSmall[[All,3]]]}];


(* ::Text:: *)
(*Create density plots of the errors for both large and small parameter ranges of the upper incomplete gamma function (takes around 50 seconds on a standard laptop, due to large number of entries in gammaErrorsBig(Small)).*)


GraphicsRow[{ListDensityPlot[gammaErrorsBig,ColorFunction->"ThermometerColors",Frame->True,PlotRange->All,FrameLabel->{{MaTeX["x"],""},{MaTeX["\\nu"],""}},ImageSize->Medium,PlotLegends->Automatic,InterpolationOrder->0],
ListDensityPlot[gammaErrorsSmall,PlotTheme->"Scientific",ColorFunction->"ThermometerColors",Frame->True,PlotRange->All,FrameLabel->{{MaTeX["x"],""},{MaTeX["\\nu"],""}},ImageSize->Medium,PlotLegends->Automatic,InterpolationOrder->0]},ImageSize->950]


(* ::Section:: *)
(*Figure 4*)


(* ::Subsection:: *)
(*Analytic formulas for special cases of the (regularized) Epstein zeta function*)


(* ::Text:: *)
(*Define analytic formulas for special cases of the Epstein zeta function. These formulas represent closed-form expressions for specific parameters.*)


epsteinZeta[\[Nu]_,{{1}},{-(1/2)},{0}]:=\!\(TraditionalForm\`2\  \*
TemplateBox[{"\[Nu]", RowBox[{"1", "/", "2"}]},
"HurwitzZeta"]\)
epsteinZeta[\[Nu]_,{{1,0},{0,2}},{-1,-2},{0,0}]:=\!\(TraditionalForm\`2\  \((\(-
\*SuperscriptBox[\(2\), \(-
\*FractionBox[\(\[Nu]\), \(2\)]\)]\) +
\*SuperscriptBox[\(2\), \(1 - \[Nu]\)] + 1)\)\  \*
TemplateBox[{FractionBox["\[Nu]", "2"]},
"DirichletBeta"]\  \*
TemplateBox[{FractionBox["\[Nu]", "2"]},
"Zeta"]\)
epsteinZeta[\[Nu]_,{{1,1/2},{0,Sqrt[3]/2}},{0,0},{0,0}]:=\!\(TraditionalForm\`2
\*SuperscriptBox[\(3\), \(1 -
\*FractionBox[\(\[Nu]\), \(2\)]\)]\  \*
TemplateBox[{FractionBox["\[Nu]", "2"]},
"Zeta"]\  \((\*
TemplateBox[{FractionBox["\[Nu]", "2"], FractionBox["1", "3"]},
"Zeta2"] - \*
TemplateBox[{FractionBox["\[Nu]", "2"], FractionBox["2", "3"]},
"Zeta2"])\)\)
epsteinZeta[\[Nu]_,{{1,0,0},{0,1,0},{0,0,2}},{0,0,-(1/2)},{1/2,0,0}]:=\!\(TraditionalForm\`
\*SuperscriptBox[\(4\), \(\[Nu]/2\)]\  \*
TemplateBox[{RowBox[{"\[Nu]", "-", "1"}]},
"DirichletBeta"]\)
epsteinZeta[\[Nu]_,6 IdentityMatrix[3],{-1,-1,-1},{1/12,1/12,1/12}]:=\!\(TraditionalForm\`
\*SuperscriptBox[\(3\), \(-
\*FractionBox[\(\[Nu]\), \(2\)]\)]\  \*
TemplateBox[{RowBox[{"\[Nu]", "-", "1"}]},
"DirichletBeta"]\)
epsteinZeta[\[Nu]_,DiagonalMatrix[{2 Sqrt[2],4,2}],{0,-1,-1},{1/(4 Sqrt[2]),0,0}]:=\!\(TraditionalForm\`
\*SuperscriptBox[\(2\), \(1 -
\*FractionBox[\(\[Nu]\), \(2\)]\)]\  \*
TemplateBox[{RowBox[{"\[Nu]", "-", "1"}]},
"DirichletBeta"]\)
epsteinZeta[\[Nu]_,IdentityMatrix[4],{1/2,0,0,0},{0,0,0,0}]:=2^\[Nu] (DirichletBeta[\[Nu]/2] DirichletBeta[\[Nu]/2-1]+DirichletLambda[\[Nu]/2] DirichletLambda[\[Nu]/2-1])
epsteinZeta[\[Nu]_,IdentityMatrix[6],ConstantArray[0,6],{1/2,1/2,0,0,0,0}]:=\!\(TraditionalForm\`4\  \*
TemplateBox[{RowBox[{FractionBox["\[Nu]", "2"], "-", "2"}]},
"DirichletBeta"]\  \*
TemplateBox[{FractionBox["\[Nu]", "2"]},
"DirichletEta"]\)
epsteinZeta[\[Nu]_,IdentityMatrix[8],ConstantArray[0,8],ConstantArray[1/2,8]]:=\!\(TraditionalForm\`\(-16\)\  \*
TemplateBox[{RowBox[{FractionBox["\[Nu]", "2"], "-", "3"}]},
"DirichletEta"]\  \*
TemplateBox[{FractionBox["\[Nu]", "2"]},
"Zeta"]\)


(* ::Text:: *)
(*Define the regularized version of the Epstein zeta function. This includes the 'shat' singularity for handling the regularzation term and the main 'epsteinZetaReg' function that combines the regular Epstein zeta function with the regularzation term.*)


shat[\[Nu]_,A_, y_]:=With[{d=Length@A},If[IntegerQ[(\[Nu]-d)/2]\[And]\[Nu]>(d-1/2),With[{k=Round[(\[Nu]-d)/2]},Pi^(k+d/2)/Gamma[k+d/2](-1)^(k+1)/k! (Pi y . y)^k Log[Pi y . y]],Pi^(\[Nu]-d/2)*(y . y)^((\[Nu]-d)/2)*Gamma[(d-\[Nu])/2]/Gamma[\[Nu]/2]]];
epsteinZetaReg[\[Nu]_,A_,x_,y_]:=If[PossibleZeroQ@Norm[y],epsteinZeta[\[Nu],A,x,y],SetPrecision[Exp[2 Pi I x . y]epsteinZeta[\[Nu],A,x,y]-shat[\[Nu],A,y]/Abs[Det[A]],200]];


(* ::Subsection:: *)
(*Load Epstein zeta reference values*)


(* ::Text:: *)
(*List the problem parameters with it's figure name.*)


benchmarkParameters= {{{{1}},{-(1/2)},{0},"S_1","s1"}
	, {{{1,0},{0,2}},{-1,-2},{0,0},"S_2^{(1)}","s21"}
	, {{{1,1/2},{0,Sqrt[3]/2}},{0,0},{0,0},"S_2^{(2)}","s22"}
	, {{{1,0,0},{0,1,0},{0,0,2}},{0,0,-(1/2)},{1/2,0,0},"S_3^{(1)}","s31"}
	, {6 IdentityMatrix[3],{-1,-1,-1},{1/12,1/12,1/12},"S_3^{(2)}","s32"}
	, {DiagonalMatrix[{2 Sqrt[2],4,2}],{0,-1,-1},{1/(4 Sqrt[2]),0,0},"S_3^{(3)}","s33"}
	,{IdentityMatrix[4],{1/2,0,0,0},{0,0,0,0},"S_4","s4"}
	,{IdentityMatrix[6],ConstantArray[0,6],{1/2,1/2,0,0,0,0},"S_6","s6"}
	, {IdentityMatrix[8],ConstantArray[0,8],ConstantArray[1/2,8],"S_8","s8"}};


(* ::Text:: *)
(*Define a function to load reference values from CSV files. Each line in the CSV contains: nu value, real part of Epstein zeta (or regularized version), imaginary part of Epstein zeta (or regularized version), and evaluation time.  After compiling the library with meson, reference values can be regenerated by *)
(**)
(*	1. 	got to project: cd C*)
(*	2.	execute binary: ./build/benchmark/epsteinlib_benchmark_epstein*)


referenceValues[funcStr_,problemStr_]:=Import[SetDirectory[NotebookDirectory[]]<> "/../../benchmark/"<>"csv/"<>funcStr<>"_"<>problemStr<>".csv"];


(* ::Text:: *)
(*Load reference values for both the regular and regularized Epstein zeta functions. This includes nu parameters, complex Epstein zeta values, and regularized Epstein zeta values for each benchmark case.*)


Do[\[Nu]Parameters[str]=referenceValues["epsteinZeta",str][[All,1]],{str,benchmarkParameters[[All,5]]}];
Do[epsteinZetaCValues[str]=referenceValues["epsteinZeta",str][[All,2]]+I*referenceValues["epsteinZeta",str][[All,3]],{str,benchmarkParameters[[All,5]]}]
Do[epsteinZetaRegCValues[str]=referenceValues["epsteinZetaReg",str][[All,2]]+I*referenceValues["epsteinZetaReg",str][[All,3]],{str,benchmarkParameters[[All,5]]}]


(* ::Text:: *)
(*Define an error function that computes the minimum of absolute and relative errors.*)


error[ref_,comp_]:=Block[{errorAbs,errorRel},errorAbs=Abs[ref-comp];errorRel=If[PossibleZeroQ@ref,1,errorAbs/Abs[ref]];Min[errorAbs,errorRel]]


(* ::Text:: *)
(*Calculate the errors between the C implementation results and high-precision Mathematica results for the (regularized) Epstein zeta function. This comparison helps validate the accuracy and performance of the C implementation across various parameter configurations.*)


errors = Table[ParallelTable[{\[Nu]Parameters[p[[5]]][[i]],Quiet[error[epsteinZetaCValues[p[[5]]][[i]],SetPrecision[epsteinZeta[Rationalize[\[Nu]Parameters[p[[5]]][[i]],10^(-20)],p[[1]],p[[2]],p[[3]]],200]]]},{i,Length[\[Nu]Parameters[p[[5]]]]}],{p,benchmarkParameters}];
(* The errors of the regularized Epstein zeta function are similar to the non-regularized version *)
(* errors = Table[ParallelTable[{\[Nu]Parameters[p[[5]]][[i]],Quiet[error[epsteinZetaRegCValues[p[[5]]][[i]],SetPrecision[epsteinZetaReg[Rationalize[\[Nu]Parameters[p[[5]]][[i]],10^(-20)],p[[1]],p[[2]],p[[3]]],200]]]},{i,Length[\[Nu]Parameters[p[[5]]]]}],{p,benchmarkParameters}]; *)


(* ::Subsection:: *)
(*Plot Errors*)


(* ::Text:: *)
(*Generate a grid of log plots showing the errors for each benchmark case. These plots visualize the accuracy of the Epstein zeta function implementation across different parameter ranges and special cases.*)


GraphicsGrid[ArrayReshape[Table[ListLogPlot[errors[[i]], PlotTheme->"Scientific",PlotLabel->MaTeX@benchmarkParameters[[i]][[4]],FrameLabel->{MaTeX@"\\nu",MaTeX@"E"},PlotRange->{2 10^(-12),10^(-20)},PlotStyle->{PointSize[0.02],ColorData[97,1]}, GridLines->Automatic ],{i,Length[benchmarkParameters]}],{Length[benchmarkParameters]/3,3}],ImageSize->Full]


(* ::Section:: *)
(*Figure 5*)


(* ::Text:: *)
(*Prepare timing data for performance analysis of the Epstein zeta function. This includes loading execution times from reference values and organizing them for plotting.*)


timings = Table[Transpose[{\[Nu]Parameters[str],referenceValues["epsteinZeta",str][[All,4]]}],{str,benchmarkParameters[[All,5]]}];
(*timings of the regularized Epstein zeta function are similar to those of the non regularized version*)
(*timings = Table[Transpose[{\[Nu]Parameters[str],referenceValues["epsteinZetaReg",str][[All,4]]}],{str,benchmarkParameters[[All,5]]}];*)


GraphicsGrid[ArrayReshape[Table[ListLogPlot[timings[[i]],PlotTheme->"Scientific",PlotLabel->MaTeX@benchmarkParameters[[i]][[4]],FrameLabel->{MaTeX@"\\nu",MaTeX@"\\text{time } [s]"},PlotStyle->{PointSize[0.02],ColorData[97,1]}, GridLines->Automatic ],{i,Length[benchmarkParameters]}], {Length[benchmarkParameters]/3,3}], ImageSize->Full]


(* ::Section:: *)
(*Figure 6*)


(* ::Text:: *)
(*Set up the dispersion relation in 3D using the Epstein zeta function. This involves defining the dispersion relation function and generating a table of values for different parameters.*)


epsteinZetaDispersion[\[Nu]_,k_]:=Re[EpsteinZeta[\[Nu],IdentityMatrix[Length@k],ConstantArray[0,Length@k],k]];
epsteinZetaDispersionTab3D=Table[{q,epsteinZetaDispersion[\[Nu],{0,0,0}]-epsteinZetaDispersion[\[Nu],{q ,0,0}]},{\[Nu],{3.5,4,5}},{q,-1/2,1/2,1/100}];


Show[Plot[59\[Nu],{\[Nu],0,0.28},PlotStyle->{Black,Opacity[0.7]},Axes->False,Frame->True,PlotTheme->"Scientific",FrameLabel->{MaTeX["k_1"],MaTeX["\\hbar \\,\\omega(\\bm k)/(JS)"]},PlotLegends->False, Epilog->{Text[Style[MaTeX["\\sim k_1^2"]],{-0.41,35}],Text[MaTeX["\\sim\\sqrt{k_1}"],{0.4,35}], Text[MaTeX["\\sim k_1 "],{0.32,18}]}],
Plot[53Sqrt[\[Nu]],{\[Nu],0,0.37},PlotStyle->{Black,Opacity[0.7]}],Plot[220\[Nu]^2,{\[Nu],-0.38,0},PlotStyle->{Black,Opacity[0.7]}],
ListLinePlot[epsteinZetaDispersionTab3D,PlotLegends->Placed[(MaTeX["\\nu ="<>ToString[#]]&/@{3.5,4,5}),{Center,Top}]],PlotRange->{{-0.51,0.51},{-1,40}}]


(* ::Section:: *)
(*Figure 7*)


(* ::Text:: *)
(*Define the energy of the Casimir effect in 1D using both the Epstein zeta function and an asymptotic expression. This allows for comparison between the exact and approximate methods.*)


casimirEnergyEpstein[z_,L_:1]:=Pi EpsteinZeta[-1,{{1/z,0,0},{0,1/L,0},{0,0,1/L}},ConstantArray[0,3],ConstantArray[0,3]]
casimirEnergyAsymptotic[z_,L_:1]:=8(-Pi^2/720)(L/z)^3


(* ::Text:: *)
(*Implement a method to calculate the Casimir force by taking the derivative of the energy using finite differences. This includes defining helper functions for coefficient calculation and derivative approximation.*)


coeffH[derPoints_,derOrder_:1]:=LinearSolve[Table[If[i==1,1,(1/(i-1)!)*(j-((derPoints+1)/2))^(i-1)],{i,1,derPoints},{j,1,derPoints}],Table[If[i==derOrder+1,(derPoints+1)*derOrder,0],{i,1,derPoints}]]
derivative[f_,x_,derPoints_,derOrder_,stepsize_:1/10000]:=With[{coeffs=coeffH[derPoints,derOrder],deltax = ((stepsize^derOrder)*derOrder*(derPoints+1))},Sum[coeffs[[i]]*f[x+stepsize(i-((derPoints+1)/2))]/deltax,{i,1,derPoints}]]


(* ::Text:: *)
(*Force of the Casimir effect in 1D by taking finite differences of 9'th order.*)


casimirForceEpstein[z_]:=With[{derPoints=10,derOrder=1},-derivative[casimirEnergyEpstein,z,derPoints,derOrder]]
casimirForceAsymptotic[d_]:=With[{L=1},-8 Pi^2L^2/(240 d^4)]
casimirEnergyCurveFit[d_]:=8Pi Exp[-2\[Pi] /d]/d^3


(* ::Text:: *)
(*Generate plots to visualize the Casimir force calculations. This includes comparisons between the exact Epstein zeta method, asymptotic approximation, and the differences between these approaches.*)


GraphicsRow[{Show[ListPlot[Table[{d,casimirForceEpstein[d]},{d,0.25,0.65,0.01}],PlotTheme->"Scientific", PlotRange->{{0.23,0.68},Automatic},FrameLabel->{MaTeX["L"], MaTeX["F"]}, GridLines->Automatic, Frame -> True, Axes -> False, PlotStyle -> {Red,PointSize[0.02],Opacity[0]}, ImageSize->Medium],
Plot[casimirForceAsymptotic[d],{d,0.242,0.65},PlotStyle->{Black},PlotLegends->Placed[{MaTeX ["\\text{asymptotic}"],None},{Center,Right}]],
ListPlot[Table[{d,casimirForceEpstein[d]},{d,0.25,0.65,0.01}],PlotStyle->{Red,PointSize[0.02]},PlotLegends->Placed[{MaTeX["\\text{Epstein zeta}"],None},{Center,Right}]]],
Show[ListLogPlot[ Table[{d,Re[casimirForceEpstein[d]-casimirForceAsymptotic[d]]},{d,0.25,0.64,0.01}],PlotTheme->"Scientific",PlotRange->{{0.23,0.68},Automatic},FrameLabel ->{MaTeX["L"], MaTeX["F_{Z}-F_{\\rm as}"]}, GridLines->Automatic, Frame -> True, Axes -> False, PlotStyle -> {Opacity[0],PointSize[0.02]}, ImageSize->Medium],
LogPlot[casimirEnergyCurveFit[d],{d,0.2402,0.65},PlotStyle->{Black},PlotLegends->Placed[{MaTeX ["\\text{curve fit}"],None},{Center,Right}]],
ListLogPlot[Table[{d,Re[casimirForceEpstein[d]-casimirForceAsymptotic[d]]},{d,0.25,0.64,0.01}],PlotTheme->"Scientific",PlotStyle -> {Red,PointSize[0.02]},PlotLegends->Placed[{MaTeX["\\text{true difference}"],None},{Center,Right}]]]},ImageSize->800]


(* ::Section:: *)
(*Figure 8*)


(* ::Text:: *)
(*Calculate the Casimir energy in 3D using the Epstein zeta function. *)


casimirEnergy3D[d_]:=Pi EpsteinZeta[-1,DiagonalMatrix[1/d],ConstantArray[0,Length@d],ConstantArray[0,Length@d]]
CasimirTab3D= Flatten[ParallelTable[{i1,i2,Re[casimirEnergy3D[{i1,i2,1/(i1 i2)}]]},{i1,0.69,1.31,0.0155},{i2,0.69,1.31,0.0155}],1];


ListContourPlot[CasimirTab3D,PlotLegends->Automatic]
