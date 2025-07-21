#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for epsteinZetaReGDer in epsteinZeta.c*)


(* ::Subsection:: *)
(*Setup*)


(* ::Text:: *)
(*Load Mathematica wrapper in every kernel. Remove the "Quiet" if unsure about the parallel evaluation.*)


1+1


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


(* ::Subsection:: *)
(*Reference functions*)


(* ::Text:: *)
(*Derivatives of the upper Crandall function.*)


p[y_,\[Alpha]_,\[Beta]_]:=Apply[Times,(-\[Pi])^(\[Alpha]-\[Beta])Binomial[\[Alpha],\[Beta]](\[Alpha]-\[Beta])!/(\[Alpha]-2\[Beta])!(2y)^(\[Alpha]-2\[Beta])]
GReg[\[Nu]_,z_]=Piecewise[{
	{-2/\[Nu],Norm[z]<10^(-64)&&Not[Element[-\[Nu]/2,NonNegativeIntegers]]}
	,{-Gamma[\[Nu]/2,0,Pi z . z]/(Pi z . z)^(\[Nu]/2),Norm[z]>=10^(-64)&&Not[Element[-\[Nu]/2,NonNegativeIntegers]]}}
];
GRegDer[\[Nu]_,y_,\[Alpha]_]:=Total[Flatten[Table[p[y,\[Alpha],\[Beta]]GReg[\[Nu]+2Total[\[Alpha]-\[Beta]],y],{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]/2]]}]]];
G[\[Nu]_,y_]=Piecewise[{{0,Norm[y]>=4},{ -2/\[Nu],Norm[y]<10^(-64)},{ExpIntegralE[1-\[Nu]/2,Pi y . y],Norm[y]>=10^(-64)&&Norm[y]<4}}];
GDer[\[Nu]_,y_,\[Alpha]_]:=Total[Flatten[Table[p[y,\[Alpha],\[Beta]]G[\[Nu]+2Total[\[Alpha]-\[Beta]],y],{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]/2]]}]]]


(* ::Text:: *)
(*Derivatives of sums in real and reciprocal space.*)


sumRealDer[\[Nu]_,A_,x_,y_,\[Alpha]_,sumLim_]:=Total@Table[Apply[Times,(-2Pi I z)^\[Alpha]]G[\[Nu],z]Exp[-2Pi I y . z],{z,(A . #-x)&/@Tuples[Range[-sumLim,sumLim],Length[A]]}]
sumFourierNonZeroDer[\[Nu]_,A_,x_,y_,\[Alpha]_,sumLim_]:=Total@Table[GDer[Length[\[Alpha]]-\[Nu],k+y,\[Alpha]]Exp[-2Pi I x . k],{k,(Inverse[Transpose[A]] . #)&/@Select[Tuples[Range[-sumLim,sumLim],Length[A]],Not[PossibleZeroQ[Norm[#]]]&]}]


(* ::Text:: *)
(*Prototype for epsteinZetaDer without special cases such non-positive even integers nu.*)


epsteinZetaRegDer[\[Nu]Ext_,AExt_,xExt_,yExt_,\[Alpha]_,sumLim_:10]:=Module[
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
	;sumReci=SetPrecision[sumFourierNonZeroDer[\[Nu]N,AN,xN,yN,\[Alpha],sumLim]+GRegDer[Length[AN]-\[Nu]N,yN,\[Alpha]],prec]
	;res=SetPrecision[\[Mu]^(Total[\[Alpha]]-\[Nu]N)\[Pi]^(\[Nu]N/2)/Gamma[\[Nu]N/2](sumReal+1/Abs[Det[AN]]sumReci),prec]
	]
	;Return[res]
]


(* ::Subsection:: *)
(*Create Epstein zeta regularized prototype reference values*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[100];
Do[
	\[Nu]0[i]=10 randRational
	;a0[i]=IdentityMatrix@2
	;x0[i]=If[i<10,{0,0},{randRational,randRational}]
	;y0[i]={randRational,randRational}
	;\[Alpha]0[i]={randDerivative,randDerivative}
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


zDerTab=ParallelTable[
	With[{result=N[epsteinZetaRegDer[\[Nu]0[i],a0[i],x0[i],y0[i],\[Alpha]0[i]],30]}
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
Export["csv/epsteinZetaRegDer_prototype_Ref.csv",zDerTab];
