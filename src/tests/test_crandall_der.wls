#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for epsteinZeta . c*)


(* ::Subsection:: *)
(*Setup reference functions*)


g[\[Nu]_,z_]=Piecewise[{{-2/\[Nu],PossibleZeroQ[Norm[z]]},{ExpIntegralE[1-\[Nu]/2,Pi z . z],Not[PossibleZeroQ[Norm[z]]]}}];
gDer[\[Nu]_,z_,\[Alpha]_]=Derivative[0,\[Alpha]][g][\[Nu],z];


(* ::Subsection:: *)
(*Create reference values*)


randRational := RandomInteger[{-10^(16),10^(16)}]/10^(16)
<<<<<<< Updated upstream
randDerivative := RandomInteger[{0,10}]
=======
randDerivative := RandomInteger[{0,5}]
>>>>>>> Stashed changes


(* ::Text:: *)
(*Setup test parameters*)


range = Range[200];
Do[
\[Nu]0[i]=With[{\[Nu]=10 randRational},If[PossibleZeroQ@(\[Nu]-2),1/2,\[Nu]]];
z0[i]={randRational,randRational,randRational};
\[Alpha]0[i]={randDerivative,randDerivative,randDerivative};
,{i,Last[range]}
]


(* ::Text:: *)
(*Generate reference values*)


crandallDerTab=ParallelTable[With[{g=gDer[\[Nu]0[i],z0[i],\[Alpha]0[i]]},{N[\[Nu]0[i],30], N[z0[i][[1]],30],N[z0[i][[2]],30],N[z0[i][[3]],30],\[Alpha]0[i][[1]],\[Alpha]0[i][[2]],\[Alpha]0[i][[3]],N[Re[g],30],N[Im[g],30]}],{i,range}];


(* ::Text:: *)
(*Save reference values to file*)


SetDirectory[NotebookDirectory[]];
Export["csv/crandall_g_der_Ref.csv",crandallDerTab];


(* ::Subsection:: *)
(*Deleteme*)


upperCrand[\[Nu]_,z_]:=Piecewise[{{0,Norm[z]>=4},{ -2/\[Nu],Norm[z]<10^(-64)},{ExpIntegralE[1-\[Nu]/2,Pi z . z],Norm[z]>=10^(-64)&&Norm[z]<4}}];
upperCrandDer[\[Nu]_,z_,\[Alpha]_]:=Total[
Table[
p[\[Alpha],\[Beta],z]
(*upperCrand[\[Nu]+2multNorm[\[Alpha]]-2multNorm[\[Beta]],z]*)
,{\[Beta],betaSetMult[\[Alpha]/2]}
]
]


crandallDerTab=ParallelTable[With[{g=upperCrandDer[\[Nu]0[i],z0[i],\[Alpha]0[i]]},{N[\[Nu]0[i],30], N[z0[i][[1]],30],N[z0[i][[2]],30],\[Alpha]0[i][[1]],\[Alpha]0[i][[2]],N[Re[g],30],N[Im[g],30]}],{i,range}];
