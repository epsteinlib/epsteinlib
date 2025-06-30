#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for epsteinZeta . c*)


(* ::Subsection:: *)
(*Setup reference functions*)


p[y_,\[Alpha]_,\[Beta]_]:=Apply[Times,(-\[Pi])^(\[Alpha]-\[Beta])Apply[Binomial,Transpose[{\[Alpha],\[Beta]}],{1}](\[Alpha]-\[Beta])!/(\[Alpha]-2\[Beta])!(2y)^(\[Alpha]-2\[Beta])]


g[\[Nu]_,z_]=Piecewise[{{-2/\[Nu],PossibleZeroQ[Norm[z]]},{ExpIntegralE[1-\[Nu]/2,Pi z . z],Not[PossibleZeroQ[Norm[z]]]}}];
gDer[\[Nu]_,z_,\[Alpha]_]=Derivative[0,\[Alpha]][g][\[Nu],z];


(* ::Subsection:: *)
(*Random Values*)


randRational := RandomInteger[{-10^(16),10^(16)}]/10^(16)
randDerivative := RandomInteger[{0,6}]
randDerivativeBeta[a_] := RandomInteger[{0,Floor[a/2]}]


(* ::Subsection:: *)
(*Create polynomial p reference values*)


(* ::Text:: *)
(*Setup test parameters*)


range = Range[1000];
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
(*Create Crandall derivatives reference values*)


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


crandallDerTab=ParallelTable[With[{g=gDer[\[Nu]0[i],z0[i],\[Alpha]0[i]]},{N[\[Nu]0[i],30], N[z0[i][[1]],30],N[z0[i][[2]],30],N[z0[i][[3]],30],\[Alpha]0[i][[1]],\[Alpha]0[i][[2]],\[Alpha]0[i][[3]],N[Re[g],30],N[Im[g],30]}],{i,range}];


(* ::Text:: *)
(*Save reference values to file*)


SetDirectory[NotebookDirectory[]];
Export["csv/crandall_g_der_Ref.csv",crandallDerTab];
