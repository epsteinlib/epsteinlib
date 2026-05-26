#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for epsteinZeta . c*)


(* ::Text:: *)
(*This file introduces reference functions beyond the methods employed in EpsteinLib. In particular, there are multiple strategies for computing the derivatives of helper functions for the singularity. Most of the derivatives for helper functions are provided here based on straightforward Fa\[AGrave] di Bruno applications, as well as the harmonic approach. Here, we favor using the opposite method to compute high-precision reference values, than the one employed in EspteinLib for the particular function.*)


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
hHarmonic[k_][x_,\[Alpha]_]:=Module[
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
(*Harmonic method for the derivatives of for  \|01d42b \[RightTeeArrow] (\[Pi] \|01d42b \[CenterDot] \|01d42b)^(d+2 ell).*)


YDerHarmonic[ell_/;Element[ell,NonNegativeIntegers],z_,\[Alpha]_]:=Module[
	{n,d}
	,d=Length[\[Alpha]]
	;n=Total[\[Alpha]]
	;\[Pi]^ell 2^Total[\[Alpha]]Total[Table[
		hHarmonic[k][z,\[Alpha]]FactorialPower[ell,n-k]FactorialPower[ell+d/2-1,k](z . z)^(ell+k-n)
		,{k,0,Total[\[Alpha]]/2}
	]]
]


(* ::Text:: *)
(*Derivatives of singularity for arguments nu = dim + 2k.*)


sDer[k_/;Element[k,NonNegativeIntegers],z_,\[Alpha]_]:=(
	\[Pi]^(k+Length[z]/2)/Gamma[k+Length[z]/2](-1)^(k+1)/k!Total@Table[
			Apply[Times,Binomial[\[Alpha],\[Beta]]]YDer[k,z,\[Beta]]LDer[z,\[Alpha]-\[Beta]]
			,{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]]]}
	]
)


(* ::Text:: *)
(*YDerGeneral implements the \[Alpha]-derivatives of \|01d42b \[RightTeeArrow] (\[Pi] \|01d42b \[CenterDot] \|01d42b)^\[Nu] for real arguments \[Nu].*)


polalt[z_,\[Alpha]_,\[Beta]_]:=Apply[Times,Binomial[\[Alpha],\[Beta]](\[Alpha]-\[Beta])!/(\[Alpha]-2\[Beta])!(2z)^(\[Alpha]-2\[Beta])]
YDerGeneral[\[Nu]_,z_,\[Alpha]_]:=\[Pi]^\[Nu] Total[Flatten[Table[
	FactorialPower[\[Nu],Total[\[Alpha]]-Total[\[Beta]]]polalt[z,\[Alpha],\[Beta]]/(z . z)^(Total[\[Alpha]]-Total[\[Beta]]-\[Nu])
	,{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]/2]]}
]]]


(* ::Text:: *)
(*Derivatives of singularity for real arguments \[Nu].*)


sing[\[Nu]_,y_]:=Module[
	{ell,d}
	,d=Length[y]
	;ell=Max[0,Round[(\[Nu]-d)/2]]
	;Which[
		PossibleZeroQ[\[Nu]-(d+2 ell)]
			,\[Pi]^(ell+d/2)/Gamma[ell+d/2](-1)^(ell+1)/ell!(\[Pi] y . y)^ell Log[Pi y . y]
		,True
			,\[Pi]^(\[Nu]/2)/Gamma[\[Nu]/2]Gamma[(d-\[Nu])/2](\[Pi] y . y)^((\[Nu]-d)/2)
	]
]
sDerFull[\[Nu]_,z_,\[Alpha]_]:=Module[
	{ell,d}
	,d=Length[z]
	;ell=Max[0,Round[(\[Nu]-d)/2]]
	;Which[
		PossibleZeroQ[\[Nu]-(d+2 ell)]
			,sDer[ell,z,\[Alpha]]
		,True
			,\[Pi]^(\[Nu]/2)Gamma[(d-\[Nu])/2]/Gamma[\[Nu]/2]YDerGeneral[(\[Nu]-d)/2,z,\[Alpha]]
	]
]


(* ::Text:: *)
(*This is the harmonic method for the singularity derivatives as employed in EpsteinLib.*)


sDerHarmonic[\[Nu]_,z_,\[Alpha]_]:=Module[
	{n,k,d}
	,d=Length[\[Alpha]]
	;n=Total[\[Alpha]]
	;If[PossibleZeroQ[Norm[z]],Return[0]]
	;Sum[
		Module[
			{ell,m,c,q,res}
			,m=n-2k
			;ell=Max[m+k,Round[(\[Nu]-d)/2]]
			;hHarmonic[k][z,\[Alpha]]Which[
				PossibleZeroQ[\[Nu]-(d+2ell)]
					,2^k (-2 \[Pi]^2)^(m+k)/FactorialPower[\[Nu]/2-k-1,m]sing[\[Nu]-2(m+k),z]
				(1+(HarmonicNumber[ell]-HarmonicNumber[ell+k-(m+2k)]+PolyGamma[ell+d/2]-PolyGamma[ell+d/2-k])/Log[Pi z . z])
				,True
					,Gamma[(d-\[Nu])/2+k+m]/Gamma[\[Nu]/2-k](-1)^(k+m)2^(2k+m)Pi^(\[Nu]-d/2)(z . z)^((\[Nu]-d)/2-k-m)
			]
		]
		,{k,0,n/2}
	]
]


(* ::Text:: *)
(*Crandall functions.*)


GLower[\[Nu]_,z_]:=Piecewise[
	{{2/\[Nu],Norm[z]<10^(-64)}}
	,Gamma[\[Nu]/2,0,Pi z . z]/(Pi z . z)^(\[Nu]/2)
]
GUpper[\[Nu]_,z_]:=Piecewise[{
	{-2/\[Nu],Norm[z]<10^(-64)}
	,{ExpIntegralE[1-\[Nu]/2,Pi z . z],Norm[z]<4}
}]
GReg[\[Nu]_,z_]:=With[{\[Nu]Z=Round[\[Nu]]}
	,Piecewise[{{
			GUpper[\[Nu],z]-\[Pi]^(-(Length[z]-\[Nu])/2)Gamma[(Length[z]-\[Nu])/2]sDer[Round[-\[Nu]/2]
			,z,ConstantArray[0,Length[z]]],(Norm[z]>10^(-64))&&PossibleZeroQ[\[Nu]-\[Nu]Z]&&Element[-\[Nu]Z/2,NonNegativeIntegers]
		}}
		,-GLower[\[Nu],z]
	]
]


(* ::Text:: *)
(*Derivatives of Crandall functions.*)


GDer[\[Nu]_,y_,\[Alpha]_]:=Piecewise[{
	{GUpper[\[Nu],y],Total[\[Alpha]]==0}
	,{Total[Flatten[Table[p[y,\[Alpha],\[Beta]]GUpper[\[Nu]+2Total[\[Alpha]-\[Beta]],y],{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]/2]]}]]],Total[\[Alpha]]>0}
}]
GRegDer[\[Nu]_,z_,\[Alpha]_]:=Piecewise[{
	{
		GReg[\[Nu],z]
		,Total[\[Alpha]]==0
	}
	,{
		(-1)^(-\[Nu]/2)/Round[-\[Nu]/2]!(HarmonicNumber[Round[-\[Nu]/2]]-EulerGamma)YDer[Round[-\[Nu]/2],z,\[Alpha]]-Total[
			Table[Piecewise[{{0,n==Round[-\[Nu]/2]},{(-1)^n YDer[n,z,\[Alpha]]/(n-Round[-\[Nu]/2])/n!,n!=Round[-\[Nu]/2]}}],{n,0,100}]
		],Norm[z]<=1/10&&Element[-\[Nu]/2,NonNegativeIntegers]&&Total[\[Alpha]]>0
	},{
		GDer[\[Nu],z,\[Alpha]]-\[Pi]^(-(Length[z]-\[Nu])/2)Gamma[(Length[z]-\[Nu])/2]sDer[Round[-\[Nu]/2],z,\[Alpha]]
		,Norm[z]>1/10&&Element[-\[Nu]/2,NonNegativeIntegers]&&Total[\[Alpha]]>0
	},{
		Total[Flatten[Table[p[z,\[Alpha],\[Beta]]GReg[\[Nu]+2Total[\[Alpha]-\[Beta]],z],{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]/2]]}]]]
		,Not[Element[-\[Nu]/2,NonNegativeIntegers]]&&Total[\[Alpha]]>0
	}
}]
GRegDerAnalytic[\[Nu]_,z_,\[Alpha]_]:=Derivative[0,\[Alpha]][GReg][\[Nu],z]
GDerAnalytic[\[Nu]_,z_,\[Alpha]_]:=Derivative[0,\[Alpha]][GUpper][\[Nu],z]


(* ::Text:: *)
(*Regularized Crandall function.*)


GRegHarmonic[\[Nu]_,z_,k_,\[Alpha]_]:=Module[
	{ell,d,n,res}
	,d=Length[z]
	;ell=Max[0,Round[-\[Nu]/2]]
	;n=Total[\[Alpha]]
	;res=Which[
		ell>=Max[k,n-k]&&PossibleZeroQ[\[Nu]+2ell]
			,Which [
				z . z<1/100
					,(-\[Pi] z . z)^(ell+k-n)/(ell+k-n)!
					(HarmonicNumber[ell]+PolyGamma[ell+d/2]-PolyGamma[ell+d/2-k]-EulerGamma)
					-Total[Table[If[j==ell+k-n,0,(-\[Pi] z . z)^j/(j-(ell+k-n))/j!],{j,0,100}]]
				,True
					,GUpper[-2(k+ell-n),z]
					+(-\[Pi] z . z)^(ell+k-n)/(ell+k-n)!
					(Log[\[Pi] z . z]+HarmonicNumber[ell]-HarmonicNumber[ell+k-n]+PolyGamma[ell+Length[z]/2]-PolyGamma[ell+Length[z]/2-k])
			]
		,True
			,-GLower[\[Nu]+2(n-k),z]
	]
	;res
]


(* ::Subsection::Closed:: *)
(*Random Values*)


randRational := RandomInteger[{-10^(16),10^(16)}]/10^(16)
randDerivative := RandomInteger[{0,6}]
randDerivativeBeta[a_] := RandomInteger[{0,Floor[a/2]}]


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


(* ::Subsection::Closed:: *)
(*Create references values for inner sum of harmonic polynomials*)


(* ::Text:: *)
(*Generate reference arguments.*)


range = Range[100];
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


(* ::Subsection::Closed:: *)
(*Create 1D harmonic h reference values*)


range=Range[100];
Do[
	\[Alpha]0[j]={RandomInteger[{0,50}]}
	;k0[j]=RandomInteger[{0,Floor[First@\[Alpha]0[j]/2]}]
	;z0[j]={randRational}
,{j,range}
]


(* ::Text:: *)
(*Generate reference values.*)


h1DTab=ParallelTable[Flatten@{k0[j],N@z0[j],\[Alpha]0[j],N@hHarmonic[k0[j]][z0[j],\[Alpha]0[j]]},{j,range}];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/harmonic_h_1D_ref.csv",h1DTab];


(* ::Subsection::Closed:: *)
(*Create harmonic h reference values*)


(* ::Text:: *)
(*Generate reference arguments.*)


range = Range[100];
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
(*Generate reference values.*)


h3DTab=ParallelTable[Flatten@{k0[j],N@z0[j],\[Alpha]0[j],N@hHarmonic[k0[j]][z0[j],\[Alpha]0[j]]},{j,range}];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/harmonic_h_3D_ref.csv",h3DTab];


(* ::Subsection::Closed:: *)
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
Export["csv/crandall_g_der_ref.csv",crandallDerTab];


(* ::Subsection:: *)
(*Create (\[Pi] y\.b2)^k derivatives reference values for slow harmonic method*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[500];
Do[
	k0[i]=RandomInteger[{0,10}]
	;z0[i]={randRational,randRational,randRational}
	;\[Alpha]0[i]=If[i<10,{0,0,i-1},{RandomInteger[{0,10}],RandomInteger[{0,10}],RandomInteger[{0,10}]}]
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


YDerTabH=ParallelTable[Flatten[{k0[i],N[z0[i],30],\[Alpha]0[i],N[YDer[k0[i],z0[i],\[Alpha]0[i]]]}],{i,range}];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/polynomial_y_der_harmonic_ref.csv",YDerTabH];


(* ::Subsection:: *)
(*Create (\[Pi] y\.b2)^k derivatives reference values for fast, naive method*)


(* ::Text:: *)
(*Setup test parameters.  The harmonic method is out-scaled for these parameters and has an higher-error in the C implementation than the naive method used in YDer. The latter method for the computation of the derivatives of (\[Pi] y\.b2)^k (which is used in EpsteinLib)  has reasonable scaling in |\[Alpha]| and d.*)


range = Range[499];
Do[
	k0[i]=RandomInteger[{0,10}]
	;z0[i]={randRational,randRational,randRational}
	;\[Alpha]0[i]=If[i<10,{0,0,i-1},{RandomInteger[{10,50}],RandomInteger[{0,50}],RandomInteger[{0,50}]}]
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


YDerTab=Join[
	YDerTabH
	,ParallelTable[Flatten[{k0[i],N[z0[i],30],\[Alpha]0[i],N[YDer[k0[i],z0[i],\[Alpha]0[i]]]}],{i,range}]
];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/polynomial_y_der_ref.csv",YDerTab];


(* ::Subsection::Closed:: *)
(*Create singularity derivatives reference values for fast, low-precision naive method*)


(* ::Text:: *)
(*Setup test parameters.*)


range = 250;
Do[
	k0[i]=RandomInteger[{0,10}]
	;\[Nu]0[i]=randRational
	;z0[i]={randRational,randRational,randRational}
	;\[Alpha]0[i]=If[i<10,{0,0,i-1},{RandomInteger[{0,6}],RandomInteger[{0,6}],RandomInteger[{0,6}]}]
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values, both for special cases \[Nu] \[Element] d - 2k and real \[Nu].*)


sDerTab=Join[
	ParallelTable[Flatten[{N[3+2*k0[i]],N[z0[i]],\[Alpha]0[i],N[SetPrecision[sDerFull[3+2 *k0[i], z0[i],\[Alpha]0[i]],200]]}],{i,range}]
	,ParallelTable[Flatten[{N[\[Nu]0[i]],N[z0[i]],\[Alpha]0[i],N[SetPrecision[sDerFull[\[Nu]0[i], z0[i],\[Alpha]0[i]],200]]}],{i,range}]
];


(* ::Text:: *)
(*Safe reference values.*)


Export[NotebookDirectory[]<>"csv/singularity_s_der_ref.csv",sDerTab];


(* ::Subsection::Closed:: *)
(*Create singularity derivatives reference values for high-precision harmonic method*)


(* ::Text:: *)
(*We include the prior reference values and add high order derivatives on top. If one replaces "sDerTab " above with the here computed "sDerTabHarmonic" the tests for the naive method fail dramatically.*)


(* ::Text:: *)
(*Setup test parameters.*)


range = 250;
Do[
	k0[i]=RandomInteger[{0,20}]
	;\[Nu]0[i]=randRational
	;z0[i]={randRational,randRational,randRational}
	;\[Alpha]0[i]=If[i<10,{0,0,i-1},{RandomInteger[{0,10}],RandomInteger[{0,10}],RandomInteger[{0,10}]}]
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values, both for special cases \[Nu] \[Element] d - 2k and real \[Nu].*)


sDerTabHarmonic=Join[
	sDerTab
	,ParallelTable[Flatten[{N[3+2*k0[i]],N[z0[i]],\[Alpha]0[i],N[SetPrecision[sDerFull[3+2 *k0[i], z0[i],\[Alpha]0[i]],200]]}],{i,range-1}]
	,ParallelTable[Flatten[{N[\[Nu]0[i]],N[z0[i]],\[Alpha]0[i],N[SetPrecision[sDerFull[\[Nu]0[i], z0[i],\[Alpha]0[i]],200]]}],{i,range}]
];


(* ::Text:: *)
(*Safe reference values.*)


Export[NotebookDirectory[]<>"csv/singularity_s_der_harmonic_ref.csv",sDerTabHarmonic];


(* ::Subsection::Closed:: *)
(*Create regularized Crandall harmonic prototype reference values*)


(* ::Text:: *)
(*Generate random reference data.*)


randomGeneralData=Table[
	Module[
		{d,\[Nu],z,k,\[Alpha],res}
		,d=RandomChoice[Range[4]]
		;\[Nu]=RandomInteger[{-10^(16), 10^(16)}]/10^15
		;z=RandomChoice[
			{0.05,0.05,0.9}->{
				ConstantArray[0,d]
				,PadRight[{10^-14},d]
				,RandomInteger[{-10^(16), 10^(16)},{d}]/10^15
			}
		]
		;\[Alpha]=RandomInteger[{0,Floor[80/d]},{d}]
		;k=RandomInteger[{0,Floor[Total[\[Alpha]]/2]}]
		;res=SetPrecision[GRegHarmonic[\[Nu],z,k,\[Alpha]],200]
		;Flatten[{d,N@\[Nu],N@z,k,\[Alpha],N@res}]
	]
	,{500}
];


(* ::Text:: *)
(*Generate special reference data where nu an non-positive even integer.*)


randomIntegerData=Table[
	Module[
		{d,\[Nu],z,k,\[Alpha],res}
		,d=RandomChoice[Range[4]]
		;\[Nu]=2 RandomInteger[{-6,0}]
		;z=RandomChoice[
			{0.05,0.05,0.9}->{
				ConstantArray[0,d]
				,PadRight[{10^-14},d]
				,RandomInteger[{-10^(16), 10^(16)},{d}]/10^15
			}
		]
		;\[Alpha]=RandomInteger[{0,Floor[80/d]},{d}]
		;k=RandomInteger[{0,Floor[Total[\[Alpha]]/2]}]
		;res=SetPrecision[GRegHarmonic[\[Nu],z,k,\[Alpha]],200]
		;Flatten[{d,N@\[Nu],N@z,k,\[Alpha],N@res}]
	]
	,{i,499}
];
Export[NotebookDirectory[]<>"csv/crandall_gReg_harmonic_ref.csv",%];


(* ::Text:: *)
(*Combine tests data.*)


testData=Join[randomGeneralData,randomIntegerData];


(* ::Text:: *)
(*Save reference values to file*)


Export[NotebookDirectory[]<>"csv/crandall_gReg_harmonic_ref.csv",testData];
