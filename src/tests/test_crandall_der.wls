#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for epsteinZeta . c*)


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
(*Setup reference functions*)


(* ::Text:: *)
(*Harmonic polynomials.*)


cInner[n_,i_,k_,d_]=2^i Product[2n+d-2-4k-2j,{j,1,i}]
cOuter[n_,k_,d_]=2^-k Product[(2n+d-2ell)/((2n+d+2-4ell)(2n+d-4ell)),{ell,k}]
hInner[\[Alpha]_,\[Beta]_,\[Gamma]_,n_,k_,d_]:=(-1)^Total[\[Beta]]/cInner[n,Total[\[Beta]],k,d]Binomial[Total[\[Beta]]+k,k]Total[\[Beta]]!/Apply[Times,Factorial[\[Beta]]] Apply[Times,Binomial[\[Alpha],\[Alpha]+\[Beta]-\[Gamma]](\[Gamma]-\[Beta])!/(2\[Gamma]-\[Alpha]-2\[Beta])!]
h[k_][x_,\[Alpha]_]:=Module[
	{d,n}
	,d=Length[\[Alpha]]
	;n=Total[\[Alpha]]
	;cOuter[n,k,d]Total@Table[
	Apply[Times,Power[x,2\[Gamma]-\[Alpha]]]
	Total[Table[
	hInner[\[Alpha],\[Beta],\[Gamma],n,k,d]
	,{\[Beta],Tuples[Range[MapThread[Max,{0\[Alpha],\[Gamma]-\[Alpha]}],\[Gamma]-\[Alpha]/2]]}
	]]
	,{\[Gamma],Select[Tuples[Range[Ceiling[\[Alpha]/2],0\[Alpha]+Total[\[Alpha]]-k]],Total[#]==Total[\[Alpha]]-k&]}
	]
]


(* ::Text:: *)
(*Helper functions for  derivatives for \[Nu] \[Element] d - 2k.*)


p[z_,\[Alpha]_,\[Beta]_]:=Apply[Times,(-\[Pi])^(\[Alpha]-\[Beta])Binomial[\[Alpha],\[Beta]](\[Alpha]-\[Beta])!/(\[Alpha]-2\[Beta])!(2z)^(\[Alpha]-2\[Beta])]
l[z_,\[Alpha]_,\[Beta]_]:=-(-1)^Total[\[Alpha]-\[Beta]]Total[\[Alpha]-\[Beta]]!/Total[\[Alpha]-\[Beta]]Apply[Times,Binomial[\[Alpha],\[Beta]](\[Alpha]-\[Beta])!/(\[Alpha]-2\[Beta])!(2z)^(\[Alpha]-2\[Beta])]
YDer[k_/;Element[k,NonNegativeIntegers],z_,\[Alpha]_]:=Piecewise[{
	{(\[Pi] z . z)^k,Total[\[Alpha]]==0}
	,{k!Pi^k Apply[Times,\[Alpha]!] Total@Table[
			Piecewise[{
			{Apply[Times,Binomial[2\[Beta],\[Alpha]]z^(2\[Beta]-\[Alpha])/\[Beta]!],Total[\[Beta]]==k}
			,{0,Total[\[Beta]]!=k}
			}],{\[Beta],Tuples[Range[Ceiling[\[Alpha]/2],ConstantArray[k,Length[z]]]]}
		],True
	}
}]
LDer[z_,\[Alpha]_]:=Piecewise[{
	{Log[\[Pi] z . z],Total[\[Alpha]]==0}
	,{Total[Flatten[Table[
			l[z,\[Alpha],\[Beta]]/(z . z)^(Total[\[Alpha]-\[Beta]])
			,{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]/2]]}
		]]],True
	}
}]


(* ::Text:: *)
(*Derivatives of singularity for arguments nu = dim + 2k.*)


s[k_,{z1_,z2_,z3_}]=\[Pi]^(k+Length[{z1,z2,z3}]/2)/Gamma[k+Length[{z1,z2,z3}]/2](-1)^(k+1)/k!(\[Pi] {z1,z2,z3} . {z1,z2,z3})^k Log[Pi {z1,z2,z3} . {z1,z2,z3}]
sDer[k_/;Element[k,NonNegativeIntegers],z_,\[Alpha]_]:=(
	\[Pi]^(k+Length[z]/2)/Gamma[k+Length[z]/2](-1)^(k+1)/k!Total@Table[
			Apply[Times,Binomial[\[Alpha],\[Beta]]]YDer[k,z,\[Beta]]LDer[z,\[Alpha]-\[Beta]]
			,{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]]]}
	]
)


(* ::Text:: *)
(*Derivatives of Crandall functions.*)


G[\[Nu]_,z_]=Piecewise[{{-2/\[Nu],Norm[z]<10^(-64)},{ExpIntegralE[1-\[Nu]/2,Pi z . z],Norm[z]>=10^(-64)}}];
GReg[\[Nu]_,z_]=Piecewise[{
	{-2/\[Nu],Norm[z]<10^(-64)&&Not[Element[-\[Nu]/2,NonNegativeIntegers]]}
	,{-Gamma[\[Nu]/2,0,Pi z . z]/(Pi z . z)^(\[Nu]/2),Norm[z]>=10^(-64)&&Not[Element[-\[Nu]/2,NonNegativeIntegers]]}}
];
GRegDer[\[Nu]_,z_,\[Alpha]_]:=Piecewise[{
	{
		(-1)^(-\[Nu]/2)/Round[-\[Nu]/2]!(HarmonicNumber[Round[-\[Nu]/2]]-EulerGamma)YDer[Round[-\[Nu]/2],z,\[Alpha]]-Total[
			Table[Piecewise[{{0,n==Round[-\[Nu]/2]},{(-1)^n YDer[n,z,\[Alpha]]/(n-Round[-\[Nu]/2])/n!,n!=Round[-\[Nu]/2]}}],{n,0,100}]
		],(Norm[z]<=1/10)&&Element[-\[Nu]/2,NonNegativeIntegers]
	},{
		GDer[\[Nu],z,\[Alpha]]-\[Pi]^(-(Length[z]-\[Nu])/2)Gamma[(Length[z]-\[Nu])/2]sDer[Round[-\[Nu]/2],z,\[Alpha]]
		,(Norm[z]>1/10)&&Element[-\[Nu]/2,NonNegativeIntegers]
	},{
		Total[Flatten[Table[p[z,\[Alpha],\[Beta]]GReg[\[Nu]+2Total[\[Alpha]-\[Beta]],z],{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]/2]]}]]]
		,Not[Element[-\[Nu]/2,NonNegativeIntegers]]
	}
}]
GRegDerAnalytic[\[Nu]_,z_,\[Alpha]_]:=Derivative[0,\[Alpha]][GReg][\[Nu],z]
GDerAnalytic[\[Nu]_,z_,\[Alpha]_]:=Derivative[0,\[Alpha]][G][\[Nu],z]


(* ::Subsection:: *)
(*Random Values*)


randRational := RandomInteger[{-10^(16),10^(16)}]/10^(16)
randDerivative := RandomInteger[{0,6}]
randDerivativeBeta[a_] := RandomInteger[{0,Floor[a/2]}]


(* ::Subsection:: *)
(*Create polynomial p reference values*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[100];
Do[
	z0[i]=10{randRational,randRational,randRational};
	\[Alpha]0[i]={randDerivative,randDerivative,randDerivative};
	\[Beta]0[i]={randDerivativeBeta[\[Alpha]0[i][[1]]],randDerivativeBeta[\[Alpha]0[i][[2]]],randDerivativeBeta[\[Alpha]0[i][[3]]]};
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


pTab=ParallelTable[Flatten[Join[N[z0[i],30],\[Alpha]0[i],\[Beta]0[i],{N[SetPrecision[p[z0[i],\[Alpha]0[i],\[Beta]0[i]],200],30]}]],{i,range}];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/polynomial_p_Ref.csv",pTab];


(* ::Subsection:: *)
(*Create coefficients c reference values*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[999];
Do[
	n0[j]=RandomInteger[{0,20}]
	;i0[j]=RandomInteger[{0,Floor[n0[j]/2]}]
	;k0[j]=RandomInteger[{0,i0[j]}]
	;d0[j]=RandomInteger[{1,10}]
	,{j,range}
]


(* ::Text:: *)
(*Generate reference values.*)


cInnerTab=ParallelTable[{n0[j],i0[j],k0[j],d0[j],N@cInner[n0[j],i0[j],k0[j],d0[j]]},{j,range}];
cOuterTab=ParallelTable[{n0[j],k0[j],d0[j],N@cOuter[n0[j],k0[j],d0[j]]},{j,range}];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/coeffs_c_inner_Ref.csv",cInnerTab];
Export["csv/coeffs_c_outer_Ref.csv",cOuterTab];


(* ::Subsection:: *)
(*Create references values for inner sum of harmonic polynomials*)


range = Range[100];
j=1;
While[j<=Last[range]
	,
	\[Alpha]0[j]={RandomInteger[{0,60}],RandomInteger[{0,60}],RandomInteger[{0,60}]}
	;k0[j]=RandomInteger[{0,Floor[Total[\[Alpha]0[j]]/2]}]
	;tups=Select[Tuples[Range[Ceiling[\[Alpha]0[j]/2],0 \[Alpha]0[j]+Total[\[Alpha]0[j]]-k0[j]]],Total[#]==Total[\[Alpha]0[j]]-k0[j]&]
	;If[Length[tups]==0,Continue[]]
	;\[Gamma]0[j]=RandomChoice@tups
	;tups=Tuples[Range[MapThread[Max,{0\[Alpha]0[j],\[Gamma]0[j]-\[Alpha]0[j]}],\[Gamma]0[j]-\[Alpha]0[j]/2]]
	;If[Length[tups]==0,Continue[]]
	;\[Beta]0[j]=RandomChoice@tups
	;j++
]
Clear[j,tups]


(* ::Text:: *)
(*Generate reference values.*)


hInnerTab=ParallelTable[Flatten@{\[Alpha]0[j],\[Beta]0[j],\[Gamma]0[j],k0[j],N@hInner[\[Alpha]0[j],\[Beta]0[j],\[Gamma]0[j],Total[\[Alpha]0[j]],k0[j],3]},{j,range}];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/harmonic_h_inner_Ref.csv",hInnerTab];


(* ::Subsection:: *)
(*Create upper Crandall derivatives reference values*)


randRational := RandomInteger[{-10^(16),10^(16)}]/10^(16)
randDerivative := RandomInteger[{0,6}]


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[10];
Do[
	\[Nu]0[i]=With[{\[Nu]=10 randRational},If[PossibleZeroQ@(\[Nu]-2),1/2,\[Nu]]];
	z0[i]={randRational,randRational,randRational};
	\[Alpha]0[i]={randDerivative,randDerivative,randDerivative};
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


crandallDerTab=ParallelTable[With[{g=GDerAnalytic[\[Nu]0[i],z0[i],\[Alpha]0[i]]},{N[\[Nu]0[i],30], N[z0[i][[1]],30],N[z0[i][[2]],30],N[z0[i][[3]],30],\[Alpha]0[i][[1]],\[Alpha]0[i][[2]],\[Alpha]0[i][[3]],N[Re[g],30],N[Im[g],30]}],{i,range}];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/crandall_g_der_Ref.csv",crandallDerTab];


(* ::Subsection:: *)
(*Create polynomial l reference values*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[100];
Do[
	z0[i]=10{randRational,randRational,randRational};
	\[Alpha]0[i]={Max[randDerivative,1],randDerivative,randDerivative};
	\[Beta]0[i]={randDerivativeBeta[\[Alpha]0[i][[1]]],randDerivativeBeta[\[Alpha]0[i][[2]]],randDerivativeBeta[\[Alpha]0[i][[3]]]};
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


lTab=ParallelTable[Flatten[Join[N[z0[i],30],\[Alpha]0[i],\[Beta]0[i],{N[SetPrecision[l[z0[i],\[Alpha]0[i],\[Beta]0[i]],200],30]}]],{i,range}];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/polynomial_l_Ref.csv",lTab];


(* ::Subsection:: *)
(*Create (\[Pi] y\.b2)^k derivatives reference values*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[100];
Do[
	k0[i]=RandomInteger[{0,10}];
	z0[i]={randRational,randRational,randRational};
	\[Alpha]0[i]=If[i<10,{0,0,i-1},{randDerivative,randDerivative,randDerivative}];
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


YDerTab=ParallelTable[Flatten[{k0[i],N[z0[i],30],\[Alpha]0[i],N[SetPrecision[Derivative[\[Alpha]0[i]][(Pi # . #)^k0[i]&][z0[i]],200],30]}],{i,range}];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/polynomial_y_der_Ref.csv",YDerTab];


(* ::Subsection:: *)
(*Create Log[\[Pi] y\.b2] derivatives reference values*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[100];
Do[
	z0[i]={randRational,randRational,randRational};
	\[Alpha]0[i]=If[i<10,{0,0,i-1},{randDerivative,randDerivative,randDerivative}];
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


LDerTab=ParallelTable[Flatten[{N[z0[i],30],\[Alpha]0[i],N[SetPrecision[Derivative[\[Alpha]0[i]][Log[Pi # . #]&][z0[i]],200],30]}],{i,range}];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/log_l_der_Ref.csv",LDerTab];


(* ::Subsection:: *)
(*Create singularity derivatives reference values*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[100];
Do[
	k0[i]=RandomInteger[{0,10}];
	z0[i]={randRational,randRational,randRational};
	\[Alpha]0[i]=If[i<10,{0,0,i-1},{RandomInteger[{0,4}],RandomInteger[{0,4}],RandomInteger[{0,4}]}];
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


sDerTab=ParallelTable[Flatten[{k0[i],N[z0[i],30],\[Alpha]0[i],N[SetPrecision[Derivative[\[Alpha]0[i]][s[k0[i], #]&][z0[i]],200],30]}],{i,range}];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/singularity_s_der_Ref.csv",sDerTab];


(* ::Subsection:: *)
(*Create regularized Crandall derivatives reference values*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[100];
Do[
	\[Nu]0[i]=With[{\[Nu]=10 randRational},If[PossibleZeroQ@(\[Nu]-2),1/2,\[Nu]]];
	z0[i]=If[i<10,{0,0,0},{randRational,randRational,randRational}];
	\[Alpha]0[i]={randDerivative,randDerivative,randDerivative};
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


crandallgRegDerTab=ParallelTable[With[{g=GRegDer[\[Nu]0[i],z0[i],\[Alpha]0[i]]},{N[\[Nu]0[i],30], N[z0[i][[1]],30],N[z0[i][[2]],30],N[z0[i][[3]],30],\[Alpha]0[i][[1]],\[Alpha]0[i][[2]],\[Alpha]0[i][[3]],N[Re[g],30],N[Im[g],30]}],{i,range}];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/crandall_gReg_der_Ref.csv",crandallgRegDerTab];


(* ::Subsection:: *)
(*Create regularized Crandall function prototype reference values for nu = dim + 2*k*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[100];
Do[
	\[Nu]0[i]=-2 RandomInteger[{0,12}]
	;z0[i]={randRational,randRational}
	;\[Alpha]0[i]={randDerivative,randDerivative}
	,{i,range}
]


GRegKDerTab=ParallelTable[Flatten[Join[{\[Nu]0[i]},N[z0[i],16],\[Alpha]0[i],{N[SetPrecision[GRegDer[\[Nu]0[i],z0[i],\[Alpha]0[i]],200],30],i}]],{i,range}];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/crandall_gReg_der_d2k_prototype_Ref.csv",GRegKDerTab]
