#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for epsteinZetaDer in epsteinZeta.c*)


(* ::Subsection:: *)
(*Setup*)


(* ::Text:: *)
(*Load Mathematica wrapper in every kernel. Remove the "Quiet" if unsure about the parallel evaluation.*)


baseDir=FileNameJoin[{NotebookDirectory[], "..", "..", "mathematica"}];
SetDirectory[baseDir];<<EpsteinZeta.wl
ParallelEvaluate[SetDirectory[baseDir];Quiet[<<EpsteinZeta.wl]];


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


(* ::Subsection:: *)
(*Auxiliary functions*)


(* ::Text:: *)
(*Generate random rational numbers and random integers between 0 and 6.*)


randRational := RandomInteger[{-10^(16),10^(16)}]/10^(16)
randRationalStable := RandomChoice[{-1,1}](1/10+8/10 RandomInteger[{0,10^(16)}]/10^(16))
randDerivative := RandomInteger[{0,6}]


(* ::Text:: *)
(*Calculates the directional derivative (\[Del]\:22c5 v)\:207f f[x] with finite differences by Taylors method given a stencil size "size" and stencil width "h".*)


DirDerTaylorCoeffs[n_,size_]:=LinearSolve[
	Table[(1/(i-1)!)*(j-((size+1)/2))^(i-1),{i,size},{j,size}]
	,n(size+1)Map[KroneckerDelta[#,n+1]&,Range[size]]
]
DirectionalDerivative[v_,n_,size_:15,h_:1/1000][f_][x_]:=If[PossibleZeroQ[n]
	,f[x]
	,Total[DirDerTaylorCoeffs[n,size] Map[f[x+h v #]&,Range[size]-(size+1)/2]]/((h^n)*n(size+1))
]


(* ::Text:: *)
(*Set Zeta without derivatives.*)


SetZeta[\[Nu]_,a_,x_,y_]:=Exp[2 Pi I x . y]EpsteinZeta[\[Nu],a,x,y]


(* ::Subsection:: *)
(*Reference functions*)


(* ::Text:: *)
(*(Stable) analytic derivatives in 1D in x = y = 0.*)


SetZetaDerAnalytic[\[Nu]_,{{a_}},{0},{0},{n_}]:=Piecewise[{{Abs[1/a]^(\[Nu]-n)(-2Pi I)^n 2 Zeta[\[Nu]-n],EvenQ[n]},{0,OddQ[n]}}]


(* ::Text:: *)
(*(Unstable) analytic derivatives in 1D for x = 0.*)


SetZetaDerAnalytic[\[Nu]_,{{a_}},{0},{y_},{n_}]:=Abs[1/a]^(\[Nu]-n)Re@Derivative[n][2 LerchPhi[Exp[-2 Pi I #],\[Nu],0]&][a y]


(* ::Text:: *)
(*Derivatives of the upper Crandall function.*)


p[y_,\[Alpha]_,\[Beta]_]:=Apply[Times,(-\[Pi])^(\[Alpha]-\[Beta])Binomial[\[Alpha],\[Beta]](\[Alpha]-\[Beta])!/(\[Alpha]-2\[Beta])!(2y)^(\[Alpha]-2\[Beta])]
g[\[Nu]_,y_]=Piecewise[{{0,Norm[y]>=4},{ -2/\[Nu],Norm[y]<10^(-64)},{ExpIntegralE[1-\[Nu]/2,Pi y . y],Norm[y]>=10^(-64)&&Norm[y]<4}}];
gDer[\[Nu]_,y_,\[Alpha]_]:=Total[Flatten[Table[p[y,\[Alpha],\[Beta]]g[\[Nu]+2Total[\[Alpha]-\[Beta]],y],{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]/2]]}]]]


(* ::Text:: *)
(*Derivatives of sums in real and reciprocal space.*)


sumRealDer[\[Nu]_,A_,x_,y_,\[Alpha]_,sumLim_]:=Total@Table[Apply[Times,(-2Pi I z)^\[Alpha]]g[\[Nu],z]Exp[-2Pi I y . z],{z,(A . #-x)&/@Tuples[Range[-sumLim,sumLim],Length[A]]}]
sumFourierDer[\[Nu]_,A_,x_,y_,\[Alpha]_,sumLim_]:=Total@Table[gDer[Length[\[Alpha]]-\[Nu],k+y,\[Alpha]]Exp[-2Pi I x . k],{k,(Inverse[Transpose[A]] . #)&/@Tuples[Range[-sumLim,sumLim],Length[A]]}]


(* ::Text:: *)
(*Prototype for epsteinZetaDer without special cases such non-positive even integers nu.*)


setZetaDer[\[Nu]Ext_,AExt_,xExt_,yExt_,\[Alpha]_,sumLim_:10]:=Module[
	{prec,\[Nu]N,AN,xN,yN,\[Mu],sumReal,sumReci,res}
	,prec=60
	;\[Nu]N=N[\[Nu]Ext,prec]
	;AN=N[AExt,prec]
	;\[Mu]=Abs[Det[AN]]^(1/Length[AN])
	;AN=1/\[Mu] N[AExt,prec]
	;xN=1/\[Mu] N[xExt,prec]
	;yN=\[Mu] N[yExt,prec]
	;If[\[Nu]N<0.5&&PossibleZeroQ[\[Nu]N/2-Round[\[Nu]N/2]]
	,Print["No early return"];res = None
	,sumReal=SetPrecision[sumRealDer[\[Nu]N,AN,xN,yN,\[Alpha],sumLim],prec]
	;sumReci=SetPrecision[sumFourierDer[\[Nu]N,AN,xN,yN,\[Alpha],sumLim],prec]
	;res=SetPrecision[\[Mu]^(Total[\[Alpha]]-\[Nu]N)\[Pi]^(\[Nu]N/2)/Gamma[\[Nu]N/2](sumReal+1/Abs[Det[AN]]sumReci),prec]
	]
	;Return[res]
]


(* ::Subsection:: *)
(*Create Epstein zeta 1D analytic reference values in x = 0*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[100];
Do[
	\[Nu]1D[i]=10 randRational
	;a1D[i]={{1/2+i/100}}
	;x1D[i]={0}
	;y1D[i]=If[i<50,{0},{randRational}]
	;n1D[i]={randDerivative}
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


zDerTab1D=ParallelTable[
	With[{result=SetPrecision[SetZetaDerAnalytic[\[Nu]1D[i],a1D[i],x1D[i],y1D[i],n1D[i]],60]}
	,Flatten[{N[\[Nu]1D[i],30]
	,N[a1D[i],30]
	,N[x1D[i],30]
	,N[y1D[i],30]
	,n1D[i]
	,Re[result]
	,Im[result]
	}]],{i,range}
];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/setZetaDer_1D_Ref.csv",zDerTab1D];


(* ::Subsection:: *)
(*Create set zeta prototype reference values*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[100];
Do[
	\[Nu]0[i]=10 randRational
	;a0[i]={{1,1/2},{0,Sqrt@3/2}}
	;x0[i]=If[i<10,{0,0},{randRational,randRational}]
	;y0[i]={randRational,randRational}
	;\[Alpha]0[i]={randDerivative,randDerivative}
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


zDerTab=ParallelTable[
	With[{result=N[setZetaDer[\[Nu]0[i],a0[i],x0[i],y0[i],\[Alpha]0[i]],30]}
	,{N[\[Nu]0[i],30]
	,N[a0[i][[1,1]],30]
	,N[a0[i][[1,2]],30]
	,N[a0[i][[2,1]],30]
	,N[a0[i][[2,2]],30]
	,N[x0[i][[1]],30]
	,N[x0[i][[2]],30]
	,N[y0[i][[1]],30]
	,N[y0[i][[2]],30]
	,\[Alpha]0[i][[1]]
	,\[Alpha]0[i][[2]]
	,Re[result]
	,Im[result]
	}],{i,range}
];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/setZetaDer_prototype_Ref.csv",zDerTab];


(* ::Subsection:: *)
(*Create Laplacian set zeta reference values by finite differences*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[100];
Do[
	\[Nu]L[i]=10 randRational
	;aL[i]={{1,1/2},{0,Sqrt@3/2}}
	;xL[i]=If[i<10,{0,0},{randRational,randRational}]
	;yL[i]={randRationalStable,randRationalStable}
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


zLaplaceTab=ParallelTable[
	(* big stencil width and size works best for "low" precision reference values *)
	With[{result=DirectionalDerivative[{1,1},2,(*stencil size: *)25,(*stencil width: *)1/250][SetZeta[\[Nu]L[i],aL[i],xL[i],#]&][yL[i]]}
	,{N[\[Nu]L[i]]
	,N[aL[i][[1,1]]]
	,N[aL[i][[1,2]]]
	,N[aL[i][[2,1]]]
	,N[aL[i][[2,2]]]
	,N[xL[i][[1]]]
	,N[xL[i][[2]]]
	,N[yL[i][[1]]]
	,N[yL[i][[2]]]
	,Re[result]
	,Im[result]
	}],{i,range}
];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/setZetaDer_laplace_Ref.csv",zLaplaceTab];


(* ::Subsection:: *)
(*Create set zeta derivative reference values for odd derivatives*)


(* ::Text:: *)
(*We employ, that SetZetaDer[\[Nu],A,x,y,\[Alpha]] should be 0 for x\:1d62 = y\:1d62 = 0 and \[Alpha]\:1d62 is odd for some 1 <= i <= d.*)


(* does not yet work for {{1,1,0},{1,0,1},{0,1,1}} (here |Det(a)|=2) or derivatives > 5*)
range = Range[100];
Do[
	\[Nu]Odd[i]=10 randRational
	;aOdd[i]=IdentityMatrix@3
	;{xOdd[i],yOdd[i],\[Alpha]Odd[i]}=Transpose[RandomSample[{{0,0,2 RandomInteger[{0,2}]+1},{randRational,randRational,RandomInteger[{0,2}]},{randRational,randRational,RandomInteger[{0,2}]}}]]
	,{i,range}
]


zOddTab=ParallelTable[
	Flatten[{
	N[\[Nu]Odd[i]]
	,N[aOdd[i]]
	,N[xOdd[i]]
	,N[yOdd[i]]
	,\[Alpha]Odd[i]
	,0.
	,0.
	}],{i,range}
];


SetDirectory[NotebookDirectory[]];
Export["csv/setZetaDer_odd_Ref.csv",zOddTab];
