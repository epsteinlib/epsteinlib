#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for epsteinZeta.c*)


(* ::Subsection:: *)
(*Setup Epstein Zeta reference equations*)


(* ::Text:: *)
(*Crandalls Formula for \[Lambda] = 1 and x,y small.*)


g[\[Nu]_,z_]:=If[PossibleZeroQ@Norm[z],-2/\[Nu],ExpIntegralE[1-\[Nu]/2,Pi z . z]];
sumReal[\[Nu]_,A_,x_,y_,L_:8]:=Total[g[\[Nu],#-x]Exp[-2 Pi I*y . #]&/@(A . # & /@Tuples[Range[-L,L],Length@A])]
sumReciprocal[\[Nu]_,A_,x_,y_,L_:8]:=Total[g[Length@A-\[Nu],#]Exp[-2 Pi I x . #]&/@(Transpose[Inverse[A]] . #+y&/@Tuples[Range[-L,L],Length@A])];
epsteinZeta[\[Nu]_,A_,x_,y_]:=(1/Pi)^(-\[Nu]/2)/Gamma[\[Nu]/2](sumReal[\[Nu],A,x,y]+1/Abs[Det[A]] sumReciprocal[\[Nu],A,x,y]);


(* ::Text:: *)
(*Analytic formulas for special cases*)


epsteinZeta[\[Nu]_,IdentityMatrix[2],{0,0}, {-1/2,-1/2}] := -4DirichletEta[\[Nu]/2]DirichletBeta[\[Nu]/2]
epsteinZeta[\[Nu]_,IdentityMatrix[2],{-1,-1}, {1/2,1/2}] := -2^(2-\[Nu]/2) (-2+2^(\[Nu]/2)) DirichletBeta[\[Nu]/2] Zeta[\[Nu]/2]
epsteinZeta[\[Nu]_,IdentityMatrix[2],{-1,-1}, {1/2,0}] := 2^(2-\[Nu]) (-2+2^(\[Nu]/2)) DirichletBeta[\[Nu]/2] Zeta[\[Nu]/2]


(* ::Text:: *)
(*Epstein Zeta regularised as a difference*)


shat[\[Nu]_,A_, y_] := With[{d=Length@A},If[IntegerQ[(\[Nu]-d)/2]\[And]\[Nu]>(d-1/2),With[{k=Round[(\[Nu]-d)/2]},Pi^(k+d/2)/Gamma[k+d/2] (-1)^(k+1)/k! (Pi y . y)^k Log[Pi y . y]],Pi^(\[Nu]-d/2)*(y . y)^((\[Nu]-d)/2)*Gamma[(d-\[Nu])/2]/Gamma[\[Nu]/2]]];
epsteinZetaReg[\[Nu]_,A_,x_,y_]:=If[PossibleZeroQ@Norm[y],epsteinZeta[\[Nu],A,x,y],SetPrecision[Exp[2 Pi I x . y]epsteinZeta[\[Nu],A,x,y]-shat[\[Nu],A,y]/Abs[Det[A]],200]];


(* ::Subsection:: *)
(*Create 2D reference values*)


$MaxExtraPrecision=1000;


randRational := RandomInteger[{-10^(16), 10^(16)}]/10^(16)


(* ::Text:: *)
(*Setup test parameters*)


maxRange = 200;
a0[i_]:= {{1,0},{0,1}}; \[Nu]0[i_] := 9/10; x0[i_] := {0,0}; y0[i_]:= {0,0};
a0[1] = {{2,0},{0,1}};(*non identity matrix*)
x0[2] = {1,0}; (*x in lattice*)
y0[3] = {1,0}; (*y in lattice*)
\[Nu]0[4] = 2; y0[4]:={0.1,0};(*\[Nu] = dimension*)
For[i = 5, i < 20, i++,
a0[i] = {{1+i/35,0},{0,1-i/35}};
\[Nu]0[i] = (randRational);
x0[i] = {(randRational/5), (randRational/5)};
y0[i] = {(randRational/10),(randRational/10)}](*diagonal-random matrices*)
For[i = 20, i < 100, i++,
a0[i] = {{1+randRational/2,randRational/10},{0,1}};
\[Nu]0[i] = (randRational);
x0[i] = {(randRational/5), (randRational/5)};
y0[i] = {(randRational/10),(randRational/10)}](*random parameters*)
a0[i_]:=IdentityMatrix[2]; x0[i_]:={0,0}; y0[i_]:={-1/2,-1/2}(*analytic formula 1*)
For[i = 100, i < 120, i++,
\[Nu]0[i]=With[{\[Nu]=randRational},If[PossibleZeroQ@(\[Nu]-2),1/2,\[Nu]]]]
a0[i_]:=IdentityMatrix[2]; x0[i_]:={-1,-1}; y0[i_]:={1/2,1/2}(*analytic formula 2*)
For[i = 100, i < 120, i++,
\[Nu]0[i]=With[{\[Nu]=randRational},If[PossibleZeroQ@(\[Nu]-2),1/2,\[Nu]]]]
a0[i_]:=IdentityMatrix[2]; x0[i_]:={-1,-1}; y0[i_]:={1/2,0}(*analytic formula 3*)
For[i = 100, i < 120, i++,
\[Nu]0[i]=With[{\[Nu]=randRational},If[PossibleZeroQ@(\[Nu]-2),1/2,\[Nu]]]]


(* ::Text:: *)
(*Generate reference values*)


epsteinZetaRef=N[ParallelTable[With[{zeta=N[epsteinZeta[\[Nu]0[i],a0[i],x0[i],y0[i]],30]},{Re[\[Nu]0[i]], Im[\[Nu]0[i]], a0[i][[1,1]],a0[i][[1,2]], a0[i][[2,1]], a0[i][[2,2]], x0[i][[1]], x0[i][[2]], y0[i][[1]], y0[i][[2]], Re[zeta], Im[zeta]}], {i,Range[maxRange]}], 16];
epsteinZetaRegRef=N[ParallelTable[With[{zeta=N[epsteinZetaReg[\[Nu]0[i],a0[i],x0[i],y0[i]],30]},{Re[\[Nu]0[i]], Im[\[Nu]0[i]], a0[i][[1,1]],a0[i][[1,2]], a0[i][[2,1]], a0[i][[2,2]], x0[i][[1]], x0[i][[2]], y0[i][[1]], y0[i][[2]], Re[zeta], Im[zeta]}], {i,Range[maxRange]}], 16];


(* ::Text:: *)
(*Save reference values to file*)


SetDirectory[NotebookDirectory[]];
Export["csv/epsteinZeta_Ref.csv",epsteinZetaRef];
Export["csv/epsteinZetaReg_Ref.csv",epsteinZetaRegRef];
