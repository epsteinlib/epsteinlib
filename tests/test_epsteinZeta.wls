#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for epsteinZeta.c*)


(* ::Subsection::Closed:: *)
(*Setup*)


(* ::Text:: *)
(*Condition number of an invertible matrix.*)


ConditionNumber[a_]:=Module[
	{sings,kappa}
	,sings=SingularValueList[a]
	;kappa=Max[sings]/Min[sings]
	;kappa
]


(* ::Subsection::Closed:: *)
(*Setup Epstein zeta reference equations*)


(* ::Text:: *)
(*Crandall function and sums in real and reciprocal space.*)


g[\[Nu]_,z_]:=If[PossibleZeroQ@Norm[z],-2/\[Nu],ExpIntegralE[1-\[Nu]/2,Pi z . z]];
SumReal[\[Nu]_,A_,x_,y_,L_:10]:=Total[g[\[Nu],#-x]Exp[-2 Pi I*y . #]&/@(A . # & /@Tuples[Range[-L,L],Length@A])]
SumFourier[\[Nu]_,A_,x_,y_,L_:10]:=Total[g[Length@A-\[Nu],#]Exp[-2 Pi I x . #]&/@(Transpose[Inverse[A]] . #+y&/@Tuples[Range[-L,L],Length@A])];


(* ::Text:: *)
(*Epstein zeta prototype.*)


Options[epsteinZeta]={InputPrecision->100,SumLim->10};
epsteinZeta[\[Nu]Ext_,AExt_,xExt_,yExt_,OptionsPattern[]]:=Module[
	{EPS,prec,sumLim,\[Nu]N,AN,xN,yN,\[Mu],sumReal,sumReci,res}
	,EPS=10^-50
	;prec=OptionValue[InputPrecision]
	;sumLim=OptionValue[SumLim]
	;\[Nu]N=N[\[Nu]Ext,prec]
	;AN=N[AExt,prec]
	;\[Mu]=Abs[Det[AN]]^(1/Length[AN])
	;AN=1/\[Mu] AN
	;xN=1/\[Mu] N[xExt,prec]
	;yN=\[Mu] N[yExt,prec]
	;If[Abs[\[Nu]N]<EPS&&Total[\[Alpha]]>0
		,res = 0
		,sumReal=SetPrecision[SumReal[\[Nu]N,AN,xN,yN,sumLim],prec]
		;sumReci=SetPrecision[SumFourier[\[Nu]N,AN,xN,yN,sumLim],prec]
		;res=SetPrecision[\[Mu]^(-\[Nu]N)\[Pi]^(\[Nu]N/2)/Gamma[\[Nu]N/2](sumReal+1/Abs[Det[AN]]sumReci),prec]
	]
	;Return[res]
]


(* ::Text:: *)
(*Analytic formulas for special cases*)


epsteinZeta[\[Nu]_,IdentityMatrix[2],{0,0}, {-1/2,-1/2}] := -4DirichletEta[\[Nu]/2]DirichletBeta[\[Nu]/2]
epsteinZeta[\[Nu]_,IdentityMatrix[2],{-1,-1}, {1/2,1/2}] := -2^(2-\[Nu]/2) (-2+2^(\[Nu]/2)) DirichletBeta[\[Nu]/2] Zeta[\[Nu]/2]
epsteinZeta[\[Nu]_,IdentityMatrix[2],{-1,-1}, {1/2,0}] := 2^(2-\[Nu]) (-2+2^(\[Nu]/2)) DirichletBeta[\[Nu]/2] Zeta[\[Nu]/2]


(* ::Text:: *)
(*Epstein zeta regularized as a difference*)


shat[\[Nu]_,A_, y_] := With[{d=Length@A},If[IntegerQ[(\[Nu]-d)/2]\[And]\[Nu]>(d-1/2),With[{k=Round[(\[Nu]-d)/2]},Pi^(k+d/2)/Gamma[k+d/2] (-1)^(k+1)/k! (Pi y . y)^k Log[Pi y . y]],Pi^(\[Nu]-d/2)*(y . y)^((\[Nu]-d)/2)*Gamma[(d-\[Nu])/2]/Gamma[\[Nu]/2]]];
epsteinZetaReg[\[Nu]_,A_,x_,y_]:=If[PossibleZeroQ@Norm[y],epsteinZeta[\[Nu],A,x,y],SetPrecision[Exp[2 Pi I x . y]epsteinZeta[\[Nu],A,x,y]-shat[\[Nu],A,y]/Abs[Det[A]],200]];


(* ::Subsection::Closed:: *)
(*Create 2D reference values*)


$MaxExtraPrecision=1000;


randRational := RandomInteger[{-10^(16), 10^(16)}]/10^(16)


(* ::Text:: *)
(*Setup test parameters.*)


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
(*Setup test parameters.*)


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
Export["csv/epsteinZeta_ref.csv",epsteinZetaRef];
Export["csv/epsteinZetaReg_ref.csv",epsteinZetaRegRef];


(* ::Subsection::Closed:: *)
(*Test random matrices*)


(* ::Text:: *)
(*Generate special 2D reference values.*)


matrixTestData2D={
	{{1,2},{0,1}}
	,{{1,0},{2,1}}
	,{{0,2},{1,0}}
	,{{0,1},{1,1}}
	,{{1,1},{1,0}}
	,{{0,2},{1,1}}
	,{{1,1},{2,0}}
	,{{0,1},{1,1}}
};


(* ::Text:: *)
(*Generate random 2D matrix reference values.*)


While[Length[matrixTestData2D]<75
	,temp=RandomInteger[{-10^(16), 10^(16)},{2,2}]/10^15
	;If[ConditionNumber[temp]<10
		,AppendTo[matrixTestData2D,temp]
	]
	;Clear[temp]
]


(* ::Text:: *)
(*Generate random 3D matrix reference values.*)


matrixTestData3D={};
While[Length[matrixTestData3D]<25
	,temp=RandomInteger[{-10^(16), 10^(16)},{3,3}]/10^15
	;If[ConditionNumber[temp]<10
		,AppendTo[matrixTestData3D,temp]
	]
	;Clear[temp]
]


(* ::Text:: *)
(*Combine matrix data.*)


matrixTestData=Join[matrixTestData2D,matrixTestData3D];


(* ::Text:: *)
(*Generate zeta reference values.*)


epsteinZetaRandomMatricesRef=ParallelTable[
	Module[{d,\[Nu],x,y,zeta}
	,d=Length[a]
	;\[Nu]=RandomInteger[{-10^(16), 10^(16)}]/10^16
	;x=a . RandomInteger[{-10^(16), 10^(16)},{d}]/10^16/2
	;y=Inverse[Transpose[a]] . RandomInteger[{-10^(16), 10^(16)},{d}]/10^16/2
	;zeta=epsteinZeta[\[Nu],a,x,y]
	;Flatten[Join[{d},N@{\[Nu],a,x,y,Re[zeta],Im[zeta]}]]
	]
	,{a,matrixTestData}
];


(* ::Text:: *)
(*Save reference values to file*)


SetDirectory[NotebookDirectory[]];
Export["csv/epsteinZeta_random_matrices_ref.csv",epsteinZetaRandomMatricesRef];


epsteinZetaRegRandomMatricesRef=Table[
	Module[{d,\[Nu],x,y,zeta}
	,d=Length[a]
	;\[Nu]=RandomInteger[{-10^(16), 10^(16)}]/10^16
	;x=a . RandomInteger[{-10^(16), 10^(16)},{d}]/10^16/2
	;y=Inverse[Transpose[a]] . RandomInteger[{-10^(16), 10^(16)},{d}]/10^16/2
	;zeta=epsteinZetaReg[\[Nu],a,x,y]
	;Flatten[Join[{d},N@{\[Nu],a,x,y,Re[zeta],Im[zeta]}]]
	]
	,{a,matrixTestData}
];


(* ::Text:: *)
(*Save reference values to file*)


SetDirectory[NotebookDirectory[]];
Export["csv/epsteinZetaReg_random_matrices_ref.csv",epsteinZetaRegRandomMatricesRef];
