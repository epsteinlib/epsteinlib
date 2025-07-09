#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for epsteinZetaDer in epsteinZeta.c*)


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
(*Setup reference functions*)


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
(*Random Values*)


randRational := RandomInteger[{-10^(16),10^(16)}]/10^(16)
randDerivative := RandomInteger[{0,6}]


(* ::Subsection:: *)
(*Create polynomial p reference values*)


(* ::Text:: *)
(*Setup test parameters*)


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
(*Generate reference values*)


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
(*Save reference values to file*)


SetDirectory[NotebookDirectory[]];
Export["csv/setZetaDer_Ref.csv",zDerTab];
