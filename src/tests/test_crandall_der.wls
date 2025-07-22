#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for epsteinZeta . c*)


(* ::Subsection:: *)
(*Setup reference functions*)


(* ::Text:: *)
(*Helper functions for regularized derivatives for \[Nu] \[Element] d - 2k*)


YDer[k_,z_,\[Alpha]_]:=k!Pi^k Apply[Times,\[Alpha]!] Total@Table[
	Piecewise[{
		{Apply[Times,Binomial[2\[Beta],\[Alpha]]z^(2\[Beta]-\[Alpha])/\[Beta]!],Total[\[Beta]]==k}
		,{0,Total[\[Beta]]!=k}
	}],{\[Beta],Tuples[Range[Ceiling[\[Alpha]/2],ConstantArray[k,Length[z]]]]}
]


(* ::Text:: *)
(*Derivatives of Crandall functions*)


p[y_,\[Alpha]_,\[Beta]_]:=Apply[Times,(-\[Pi])^(\[Alpha]-\[Beta])Binomial[\[Alpha],\[Beta]](\[Alpha]-\[Beta])!/(\[Alpha]-2\[Beta])!(2y)^(\[Alpha]-2\[Beta])]
G[\[Nu]_,z_]=Piecewise[{{-2/\[Nu],Norm[z]<10^(-64)},{ExpIntegralE[1-\[Nu]/2,Pi z . z],Norm[z]>=10^(-64)}}];
GReg[\[Nu]_,z_]=Piecewise[{
	{-2/\[Nu],Norm[z]<10^(-64)&&Not[Element[-\[Nu]/2,NonNegativeIntegers]]}
	,{-Gamma[\[Nu]/2,0,Pi z . z]/(Pi z . z)^(\[Nu]/2),Norm[z]>=10^(-64)&&Not[Element[-\[Nu]/2,NonNegativeIntegers]]}}
];
GRegDerAnalytic[\[Nu]_,z_,\[Alpha]_]:=Derivative[0,\[Alpha]][GReg][\[Nu],z]
GDerAnalytic[\[Nu]_,z_,\[Alpha]_]:=Derivative[0,\[Alpha]][G][\[Nu],z]


(* ::Subsection:: *)
(*Random Values*)


randRational := RandomInteger[{-10^(16),10^(16)}]/10^(16)
randDerivative := RandomInteger[{0,6}]
randDerivativeBeta[a_] := RandomInteger[{0,Floor[a/2]}]


(* ::Subsection:: *)
(*Create (\[Pi] y\.b2)^k derivatives reference values*)


(* ::Text:: *)
(*Setup test parameters*)


range = Range[100];
Do[
	k0[i]=RandomInteger[{0,10}];
	z0[i]={randRational,randRational,randRational};
	\[Alpha]0[i]={randDerivative,randDerivative,randDerivative};
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values*)


YDerTab=ParallelTable[Flatten[{k0[i],N[z0[i],30],\[Alpha]0[i],N[SetPrecision[Derivative[\[Alpha]0[i]][(Pi # . #)^k0[i]&][z0[i]],200],30]}],{i,range}];


(* ::Text:: *)
(*Save reference values to file*)


SetDirectory[NotebookDirectory[]];
Export["csv/polynomial_y_der_Ref.csv",YDerTab];


(* ::Subsection:: *)
(*Create polynomial p reference values*)


(* ::Text:: *)
(*Setup test parameters*)


range = Range[800];
Do[
	y0[i]=10{randRational,randRational,randRational};
	\[Alpha]0[i]={randDerivative,randDerivative,randDerivative};
	\[Beta]0[i]={randDerivativeBeta[\[Alpha]0[i][[1]]],randDerivativeBeta[\[Alpha]0[i][[2]]],randDerivativeBeta[\[Alpha]0[i][[3]]]};
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values*)


pDerTab=ParallelTable[{N[y0[i][[1]],30],N[y0[i][[2]],30],N[y0[i][[3]],30],\[Alpha]0[i][[1]],\[Alpha]0[i][[2]],\[Alpha]0[i][[3]],\[Beta]0[i][[1]],\[Beta]0[i][[2]],\[Beta]0[i][[3]],N[SetPrecision[p[y0[i],\[Alpha]0[i],\[Beta]0[i]],200],30]},{i,range}];


(* ::Text:: *)
(*Save reference values to file*)


SetDirectory[NotebookDirectory[]];
Export["csv/polynomial_p_Ref.csv",pDerTab];


(* ::Subsection:: *)
(*Create upper Crandall derivatives reference values*)


randRational := RandomInteger[{-10^(16),10^(16)}]/10^(16)
randDerivative := RandomInteger[{0,6}]


(* ::Text:: *)
(*Setup test parameters*)


range = Range[200];
Do[
	\[Nu]0[i]=With[{\[Nu]=10 randRational},If[PossibleZeroQ@(\[Nu]-2),1/2,\[Nu]]];
	z0[i]={randRational,randRational,randRational};
	\[Alpha]0[i]={randDerivative,randDerivative,randDerivative};
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values*)


crandallDerTab=ParallelTable[With[{g=GDerAnalytic[\[Nu]0[i],z0[i],\[Alpha]0[i]]},{N[\[Nu]0[i],30], N[z0[i][[1]],30],N[z0[i][[2]],30],N[z0[i][[3]],30],\[Alpha]0[i][[1]],\[Alpha]0[i][[2]],\[Alpha]0[i][[3]],N[Re[g],30],N[Im[g],30]}],{i,range}];


(* ::Text:: *)
(*Save reference values to file*)


SetDirectory[NotebookDirectory[]];
Export["csv/crandall_g_der_Ref.csv",crandallDerTab];


(* ::Subsection:: *)
(*Create regularized Crandall derivatives reference values*)


randRational := RandomInteger[{-10^(16),10^(16)}]/10^(16)
randDerivative := RandomInteger[{0,6}]


(* ::Text:: *)
(*Setup test parameters*)


range = Range[100];
Do[
	\[Nu]0[i]=With[{\[Nu]=10 randRational},If[PossibleZeroQ@(\[Nu]-2),1/2,\[Nu]]];
	z0[i]=If[i<10,{0,0,0},{randRational,randRational,randRational}];
	\[Alpha]0[i]={randDerivative,randDerivative,randDerivative};
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values*)


crandallgRegDerTab=ParallelTable[With[{g=GRegDer[\[Nu]0[i],z0[i],\[Alpha]0[i]]},{N[\[Nu]0[i],30], N[z0[i][[1]],30],N[z0[i][[2]],30],N[z0[i][[3]],30],\[Alpha]0[i][[1]],\[Alpha]0[i][[2]],\[Alpha]0[i][[3]],N[Re[g],30],N[Im[g],30]}],{i,range}];


(* ::Text:: *)
(*Save reference values to file*)


SetDirectory[NotebookDirectory[]];
Export["csv/crandall_gReg_der_Ref.csv",crandallgRegDerTab];
