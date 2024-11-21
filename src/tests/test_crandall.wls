#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for epsteinZeta . c*)


(* ::Subsection:: *)
(*Setup reference functions*)


g[\[Nu]_,z_]:=If[PossibleZeroQ@Norm[z],-2/\[Nu],ExpIntegralE[1-\[Nu]/2,Pi z . z]];


(* ::Subsection:: *)
(*Create reference values*)


randRational := RandomInteger[{-10^(16), 10^(16)}]/10^(16)


(* ::Text:: *)
(*Setup test parameters*)


range = Range[200];
For[i = 1, i <= Last[range], i++,
\[Nu]0[i]=With[{\[Nu]=10 randRational},If[PossibleZeroQ@(\[Nu]-2),1/2,\[Nu]]];
z0[i]={randRational,randRational}]


(* ::Text:: *)
(*Generate reference values*)


crandallTab=N[ParallelTable[With[{g=g[\[Nu]0[i],z0[i]]},{\[Nu]0[i], z0[i][[1]],z0[i][[2]],Re[g],Im[g]}], {i,range}], 30];


(* ::Text:: *)
(*Save reference values to file*)


SetDirectory[NotebookDirectory[]];
Export["csv/crandall_g_Ref.csv",crandallTab];
