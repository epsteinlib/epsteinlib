#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2025-2026 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for harmonics . c*)


(* ::Subsection::Closed:: *)
(*Setup*)


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
(*Reference functions*)


(* ::Text:: *)
(*Harmonic polynomials.*)


cInner[n_,i_,k_,d_]:= Product[2n+d-2-4k-2j,{j,i+1,n/2-k}]
cOuter[n_,k_,d_]:=2^-k/cInner[n,0,k,d] Product[(2n+d-2ell)/((2n+d+2-4ell)(2n+d-4ell)),{ell,k}]
hInner[\[Alpha]_,\[Beta]_,\[Gamma]_,n_,k_,d_]:=(-1/2)^Total[\[Beta]] cInner[n,Total[\[Beta]],k,d]Binomial[Total[\[Beta]]+k,k]Total[\[Beta]]!/Apply[Times,Factorial[\[Beta]]] Apply[Times,Binomial[\[Alpha],\[Alpha]+\[Beta]-\[Gamma]](\[Gamma]-\[Beta])!/(2\[Gamma]-\[Alpha]-2\[Beta])!]
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


(* ::Subsection::Closed:: *)
(*Random Values*)


randRational := RandomInteger[{-10^(16),10^(16)}]/10^(16)


(* ::Subsection::Closed:: *)
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
Export["csv/coeffs_c_inner_ref.csv",cInnerTab];
Export["csv/coeffs_c_outer_ref.csv",cOuterTab];


(* ::Subsection:: *)
(*Create references values for inner sum of harmonic polynomials*)


(* ::Text:: *)
(*Generate reference arguments.*)


range = Range[999];
j=1;
maxn=60;
While[j<=Last[range]
	,
	\[Alpha]0[j]={temp1=RandomInteger[{1,maxn}],temp2=RandomInteger[{0,maxn-temp1}],RandomInteger[{0,maxn-temp1-temp2}]}
	;k0[j]=RandomInteger[{0,Floor[Total[\[Alpha]0[j]]/2]}]
	;tups=Select[Tuples[Range[Ceiling[\[Alpha]0[j]/2],0 \[Alpha]0[j]+Total[\[Alpha]0[j]]-k0[j]]],Total[#]==Total[\[Alpha]0[j]]-k0[j]&]
	;If[Length[tups]==0,Continue[]]
	;\[Gamma]0[j]=RandomChoice@tups
	;tups=Tuples[Range[MapThread[Max,{0\[Alpha]0[j],\[Gamma]0[j]-\[Alpha]0[j]}],\[Gamma]0[j]-\[Alpha]0[j]/2]]
	;If[Length[tups]==0,Continue[]]
	;\[Beta]0[j]=RandomChoice@tups
	;z0[j]=2{randRational,randRational,randRational};
	;j++
]
Clear[j,tups]


(* ::Text:: *)
(*Generate reference values.*)


hInnerTab=ParallelTable[Flatten@{\[Alpha]0[j],\[Beta]0[j],\[Gamma]0[j],k0[j],N@hInner[\[Alpha]0[j],\[Beta]0[j],\[Gamma]0[j],Total[\[Alpha]0[j]],k0[j],3]},{j,range}];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/harmonic_h_inner_ref.csv",hInnerTab];


(* ::Subsection:: *)
(*Create 1D harmonic h reference values*)


range=Range[999];
Do[
	\[Alpha]0[j]={RandomInteger[{0,50}]}
	;k0[j]=RandomInteger[{0,Floor[First@\[Alpha]0[j]/2]}]
	;z0[j]={randRational}
,{j,range}
]


(* ::Text:: *)
(*Generate reference values.*)


h1DTab=ParallelTable[Flatten@{k0[j],N@z0[j],\[Alpha]0[j],N@h[k0[j]][z0[j],\[Alpha]0[j]]},{j,range}];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/harmonic_h_1D_ref.csv",h1DTab];


(* ::Subsection:: *)
(*Create harmonic h reference values*)


(* ::Text:: *)
(*Generate reference arguments.*)


range = Range[200];
j=1;
While[j<=Last[range]
	,
	\[Alpha]0[j]={RandomInteger[{0,20}],RandomInteger[{0,20}],RandomInteger[{0,20}]}
	;If[Total[\[Alpha]0[j]]>30,Continue[]]
	;k0[j]=RandomInteger[{0,Floor[Total[\[Alpha]0[j]]/2]}]
	;tups=Select[Tuples[Range[Ceiling[\[Alpha]0[j]/2],0 \[Alpha]0[j]+Total[\[Alpha]0[j]]-k0[j]]],Total[#]==Total[\[Alpha]0[j]]-k0[j]&]
	;If[Length[tups]==0,Continue[]]
	;\[Gamma]0[j]=RandomChoice@tups
	;tups=Tuples[Range[MapThread[Max,{0\[Alpha]0[j],\[Gamma]0[j]-\[Alpha]0[j]}],\[Gamma]0[j]-\[Alpha]0[j]/2]]
	;If[Length[tups]==0,Continue[]]
	;\[Beta]0[j]=RandomChoice@tups
	;z0[j]=2{randRational,randRational,randRational};
	;j++
]
Clear[j,tups]


(* ::Text:: *)
(*Generate reference values for z = (1,1,1), which focuses on the precision of harmonics.*)


h3DTabUnity=ParallelTable[
		With[
			{res=SetPrecision[h[k0[j]][{1,1,1},\[Alpha]0[j]],200]}
			,Flatten@{k0[j],N@{1,1,1},\[Alpha]0[j],N@res}
		]
	,{j,range}
];


(* ::Text:: *)
(*Save reference values to file.*)


Export[NotebookDirectory[]<>"csv/harmonic_h_3D_unity_ref.csv",h3DTabUnity];


(* ::Text:: *)
(*Generate reference values for random z.*)


h3DTabRandom=ParallelTable[
		With[
			{res=SetPrecision[h[k0[j]][z0[j],\[Alpha]0[j]],200]}
			,Flatten@{k0[j],N@z0[j],\[Alpha]0[j],N@res}
		]
	,{j,range}
];


(* ::Text:: *)
(*Save reference values to file.*)


Export[NotebookDirectory[]<>"csv/harmonic_h_3D_random_ref.csv",h3DTabRandom];
