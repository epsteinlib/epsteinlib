#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Epstein zeta Benchmark*)


(* ::Text:: *)
(*Analytic formulas for special cases of the Epstein zeta function*)


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
(*Set regularized Epstein zeta reference function*)


shat[\[Nu]_,A_, y_] := With[{d=Length@A},If[IntegerQ[(\[Nu]-d)/2]\[And]\[Nu]>(d-1/2),With[{k=Round[(\[Nu]-d)/2]},Pi^(k+d/2)/Gamma[k+d/2] (-1)^(k+1)/k! (Pi y . y)^k Log[Pi y . y]],Pi^(\[Nu]-d/2)*(y . y)^((\[Nu]-d)/2)*Gamma[(d-\[Nu])/2]/Gamma[\[Nu]/2]]];
epsteinZetaReg[\[Nu]_,A_,x_,y_]:=If[PossibleZeroQ@Norm[y],epsteinZeta[\[Nu],A,x,y],SetPrecision[Exp[2 Pi I x . y]epsteinZeta[\[Nu],A,x,y]-shat[\[Nu],A,y]/Abs[Det[A]],200]];


(* ::Text:: *)
(*Call Epstein zeta function*)


SetDirectory[NotebookDirectory[]<>"../../mathematica"];
<<"EpsteinZeta.wl"


(* ::Text:: *)
(*Load EpsteinZeta.wl in every local Kernel*)


ParallelEvaluate[Quiet[<<"EpsteinZeta.wl"]];


(* ::Text:: *)
(*Compute minimum of absolute and relative error*)


error[ref_,comp_]:=Block[{errorAbs,errorRel},errorAbs=Abs[ref-comp];errorRel=If[PossibleZeroQ@ref,1,errorAbs/Abs[ref]];Min[errorAbs,errorRel]]


(* ::Text:: *)
(*Set the parameters to compare functions*)


stepsize = 1/4;
\[Nu]Parameters = Range[-10+2^(-15),10,stepsize];
fixedParameters ={
 {{{1}},{-(1/2)},{0},"S_1"},
 {{{1,0},{0,2}},{-1,-2},{0,0},"S_2^{(1)}"},
 {{{1,1/2},{0,Sqrt[3]/2}},{0,0},{0,0},"S_2^{(2)}"},
 {{{1,0,0},{0,1,0},{0,0,2}},{0,0,-(1/2)},{1/2,0,0},"S_3^{(1)}"},
 {6 IdentityMatrix[3],{-1,-1,-1},{1/12,1/12,1/12},"S_3^{(2)}"},
 {DiagonalMatrix[{2 Sqrt[2],4,2}],{0,-1,-1},{1/(4 Sqrt[2]),0,0},"S_3^{(3)}"},
 {IdentityMatrix[4],{1/2,0,0,0},{0,0,0,0},"S_4"},
 {IdentityMatrix[6],ConstantArray[0,6],{1/2,1/2,0,0,0,0},"S_6"},
 {IdentityMatrix[8],ConstantArray[0,8],ConstantArray[1/2,8],"S_8"}
};


(* ::Text:: *)
(*Calculate errors of the Epstein zeta function (approximately 25 seconds for a step size of 1/4 on 12 kernels)*)


errors = Table[ParallelTable[{\[Nu],Quiet[error[EpsteinZeta[\[Nu],p[[1]],p[[2]],p[[3]]],SetPrecision[epsteinZeta[Rationalize[\[Nu],10^(-20)],p[[1]],p[[2]],p[[3]]],200]]]},{\[Nu],\[Nu]Parameters}],{p,fixedParameters}];


(* ::Text:: *)
(*The errors of the regularized Epstein zeta function are similar to those of the Epstein zeta function*)


(*
errors = Table[ParallelTable[{\[Nu],Quiet[error[EpsteinZetaReg[\[Nu],p[[1]],p[[2]],p[[3]]],SetPrecision[epsteinZetaReg[Rationalize[\[Nu],10^(-20)],p[[1]],p[[2]],p[[3]]],200]]]},{\[Nu],\[Nu]Parameters}],{p,fixedParameters}];
*)


(* ::Text:: *)
(*Install MaTeX for L AT EX figure annotations*)


(*ResourceFunction["MaTeXInstall"][];*)
Needs["MaTeX`"];
SetOptions[MaTeX,FontSize->18];


(* ::Text:: *)
(*Plot minimum of absolute and relative error*)


GraphicsGrid[
 ArrayReshape[Table[
  ListLogPlot[
   errors[[i]],
   PlotTheme->"Scientific",
   PlotLabel->MaTeX@fixedParameters[[i]][[4]],
   FrameLabel->{MaTeX@"\\nu",MaTeX@"E"},
   PlotRange->{10^(-18), Automatic}
  ],{i,Length[fixedParameters]}],
 {Length[fixedParameters]/3,3}],
ImageSize->Full]
