#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for epsteinZetaRegDer in epsteinZeta.c*)


(* ::Subsection::Closed:: *)
(*Setup*)


(* ::Text:: *)
(*Load Mathematica wrapper in every kernel. Remove the "Quiet" if unsure about the parallel evaluation.*)


BaseDirectory=FileNameJoin[{NotebookDirectory[], "..", "..", "mathematica"}];
SetDirectory[BaseDirectory];<<EpsteinZeta.wl
ParallelEvaluate[SetDirectory[BaseDirectory];Quiet[<<EpsteinZeta.wl]];


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


(* ::Subsection::Closed:: *)
(*Auxiliary functions*)


(* ::Text:: *)
(*Bain path lattice matrix and hexagonal lattice matrix.*)


BMat[a_]=1/Sqrt[a+1]DiagonalMatrix[{Sqrt[a],1,1}] . {{1,1,0},{1,0,1},{0,1,1}};
AHex={{1,1/2},{0,Sqrt@3/2}};


(* ::Text:: *)
(*Generate random rational numbers and random integers between 0 and 6.*)


randRational := RandomInteger[{-10^(16),10^(16)}]/10^(16)
randRationalStable := RandomChoice[{-1,1}](1/10+8/10 RandomInteger[{0,10^(16)}]/10^(16))
randDerivative := RandomInteger[{0,6}]


(* ::Text:: *)
(*Calculates the directional derivative (\[Del]\:22c5 v)\:207f f[x] with finite differences by Taylors method given a stencil size "size" and stencil width "h".*)


DirDerTaylorCoeffs[n_,size_]:=LinearSolve[
	Table[(1/(i-1)!)*(j-((size+1)/2))^(i-1),{i,size},{j,size}]
	,n(size+1)Map[KroneckerDelta[#,n+1]&,Range[size]]
]
DirectionalDerivative[v_,n_,size_:15,h_:1/1000][f_][x_]:=If[PossibleZeroQ[n]
	,f[x]
	,Total[DirDerTaylorCoeffs[n,size] Map[f[x+h v #]&,Range[size]-(size+1)/2]]/((h^n)*n(size+1))
]


(* ::Subsection::Closed:: *)
(*Reference functions*)


(* ::Text:: *)
(*Helper functions for the derivatives.*)


p[z_,\[Alpha]_,\[Beta]_]:=Apply[Times,(-\[Pi])^(\[Alpha]-\[Beta])Binomial[\[Alpha],\[Beta]](\[Alpha]-\[Beta])!/(\[Alpha]-2\[Beta])!(2z)^(\[Alpha]-2\[Beta])]
l[z_,\[Alpha]_,\[Beta]_]:=-(-1)^Total[\[Alpha]-\[Beta]]Total[\[Alpha]-\[Beta]]!/Total[\[Alpha]-\[Beta]]Apply[Times,Binomial[\[Alpha],\[Beta]](\[Alpha]-\[Beta])!/(\[Alpha]-2\[Beta])!(2z)^(\[Alpha]-2\[Beta])]


(* ::Text:: *)
(*Derivatives of singularity for arguments nu = dim + 2k.*)


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
sDer[k_/;Element[k,NonNegativeIntegers],z_,\[Alpha]_]:=(
	\[Pi]^(k+Length[z]/2)/Gamma[k+Length[z]/2](-1)^(k+1)/k!Total@Table[
			Apply[Times,Binomial[\[Alpha],\[Beta]]]YDer[k,z,\[Beta]]LDer[z,\[Alpha]-\[Beta]]
			,{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]]]}
	]
)


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


(* ::Text:: *)
(*Derivatives of sums in real and reciprocal space.*)


sumRealDer[\[Nu]_,A_,x_,y_,\[Alpha]_,sumLim_]:=Total@Table[Apply[Times,(-2Pi I z)^\[Alpha]]GUpper[\[Nu],z]Exp[-2Pi I y . z],{z,(A . #-x)&/@Tuples[Range[-sumLim,sumLim],Length[A]]}]
sumFourierNonZeroDer[\[Nu]_,A_,x_,y_,\[Alpha]_,sumLim_]:=Total@Table[GDer[Length[\[Alpha]]-\[Nu],k+y,\[Alpha]]Exp[-2Pi I x . k],{k,(Inverse[Transpose[A]] . #)&/@Select[Tuples[Range[-sumLim,sumLim],Length[A]],Not[PossibleZeroQ[Norm[#]]]&]}]


(* ::Text:: *)
(*Prototype for epsteinZetaDer without special cases such non-positive even integers nu.*)


Options[epsteinZetaRegDer]={InputPrecision->120,SumLim->10};
epsteinZetaRegDer[\[Nu]Ext_,AExt_,xExt_,yExt_,\[Alpha]_,OptionsPattern[]]:=Module[
	{dim,prec,\[Nu]N,AN,xN,yN,\[Mu],sumReal,sumReci,res,k,cor,vol,sumLim}
	,prec=OptionValue[InputPrecision]
	;sumLim=OptionValue[SumLim]
	;dim = Length[AExt]
	;If[\[Nu]Ext>=dim&&PossibleZeroQ[(\[Nu]Ext-dim)/2-Round[(\[Nu]Ext-dim)/2]]
		,\[Nu]N=Round[\[Nu]Ext];If[PossibleZeroQ[Total[\[Alpha]]],Return[" bug for no derivatives in nu equals dim plus 2k"]]
		,\[Nu]N=N[\[Nu]Ext,prec]
	]
	;AN=N[AExt,prec]
	;vol=Abs[Det[AN]]
	;\[Mu]=Abs[Det[AN]]^(1/Length[AN])
	;AN=1/\[Mu] N[AExt,prec]
	;xN=1/\[Mu] N[xExt,prec]
	;yN=\[Mu] N[yExt,prec]
	;sumReal=SetPrecision[sumRealDer[\[Nu]N,AN,xN,yN,\[Alpha],sumLim],prec]
	;sumReci=SetPrecision[sumFourierNonZeroDer[\[Nu]N,AN,xN,yN,\[Alpha],sumLim]+GRegDer[Length[AN]-\[Nu]N,yN,\[Alpha]],prec]
	;res=SetPrecision[\[Mu]^(Total[\[Alpha]]-\[Nu]N)\[Pi]^(\[Nu]N/2)/Gamma[\[Nu]N/2](sumReal+1/Abs[Det[AN]]sumReci),prec]
	;k=Max[0,Round[(\[Nu]N-dim)/2]]
	;If[PossibleZeroQ[k-(\[Nu]N-dim)/2]
		,res+=Pi^(k+dim/2)/Gamma[k+dim/2](-1)^(k+1)/k!YDer[k,N[yExt],\[Alpha]]Log[\[Mu] \[Mu]]/Abs[Det[AExt]]
	]
	;Return[res]
]


(* ::Text:: *)
(*Harmonic polynomials*)


cInner[n_,i_,k_,d_]=2^i Product[2n+d-2-4k-2j,{j,1,i}];
cOuter[n_,k_,d_]=2^-k Product[(2n+d-2ell)/((2n+d+2-4ell)(2n+d-4ell)),{ell,k}];
hHarmonic[k_/;Element[k,NonNegativeIntegers]][x_List,\[Alpha]_List]:=Module[
	{res,d,n}
	,d=Length[\[Alpha]]
	;n=Total[\[Alpha]]
	;res=cOuter[n,k,d]Total@Table[
		Apply[Times,Power[x,2\[Gamma]-\[Alpha]]]
		Total[Table[
			(-1)^Total[\[Beta]]/cInner[n,Total[\[Beta]],k,d]
			 Binomial[Total[\[Beta]]+k,k]Total[\[Beta]]!/Apply[Times,Factorial[\[Beta]]]
			Apply[Times,Binomial[\[Alpha],\[Alpha]+\[Beta]-\[Gamma]](\[Gamma]-\[Beta])!/(2\[Gamma]-\[Alpha]-2\[Beta])!]
			,{\[Beta],Tuples[Range[MapThread[Max,{0\[Alpha],\[Gamma]-\[Alpha]}],\[Gamma]-\[Alpha]/2]]}
			]]
		,{\[Gamma],Select[Tuples[Range[Ceiling[\[Alpha]/2],0\[Alpha]+Total[\[Alpha]]-k]],Total[#]==Total[\[Alpha]]-k&]}
		]
	;res
]


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


(* ::Text:: *)
(*Sums for harmonic method in real and reciprocal space.*)


SumRealHarmonic[\[Nu]_,A_,x_,y_,\[Alpha]_,kIndex_,sumLim_]:=Total[Table[
	Module[{h}
		,h=hHarmonic[kIndex][z,\[Alpha]]
		;If[PossibleZeroQ[h],0,h GUpper[\[Nu],z]Exp[-2\[Pi] I y . z]]
	]
	,{z,(A . #-x)&/@Tuples[Range[-sumLim,sumLim],Length[A]]}
]]
SumFourierRegHarmonic[\[Nu]_,A_,x_,y_,\[Alpha]_,kIndex_,sumLim_]:=Total[Table[
	hHarmonic[kIndex][y+k,\[Alpha]]GUpper[Length[x]-\[Nu]+2Total[\[Alpha]]-4kIndex,k+y]Exp[-2\[Pi] I x . k]
	,{k,(Inverse[Transpose[A]] . #)&/@Select[Tuples[Range[-sumLim,sumLim],Length[A]],Not[PossibleZeroQ[Norm[#]]]&]}
]]


(* ::Text:: *)
(*Harmonic method for the regularized epstein Zeta derivatives.*)


Options[epsteinZetaRegDerHarmonic]={InputPrecision->100,SumLim->8};
epsteinZetaRegDerHarmonic[\[Nu]Ext_,AExt_,xExt_,yExt_,\[Alpha]_,OptionsPattern[]]:=Module[
	{EPS,prec,sumLim,\[Nu]True,AN,xN,yN,\[Mu],res,n,nc,m,l,d}
	,EPS=10^-50
	;prec=OptionValue[InputPrecision]
	;sumLim=OptionValue[SumLim]
	;d=Length[\[Alpha]]
	;n = Total[\[Alpha]]
	;\[Nu]True=N[\[Nu]Ext,prec]
	;AN=N[AExt,prec]
	;\[Mu]=Abs[Det[AN]]^(1/Length[AN])
	;AN=1/\[Mu] N[AExt,prec]
	;xN=1/\[Mu] N[xExt,prec]
	;yN=\[Mu] N[yExt,prec]
	;If[Abs[\[Nu]Ext]<EPS&&n>0
		,res = 0
		,res=(-2Pi)^n Total@Table[
			Module[
				{\[Nu]N,sumReal,sumReci,resk}
				,\[Nu]N=\[Nu]True-2k
				;m=Total[\[Alpha]]-2k
				;If[\[Nu]N<-1&&Abs[\[Nu]N/2-Round[\[Nu]N/2]]<EPS
					,0
					,If[Abs[\[Nu]N]<EPS
						,If[Norm[xN]>EPS
							,0
							,sumReal=-I^n hHarmonic[k][xN,\[Alpha]]
							;resk=\[Pi]^(\[Nu]N/2)sumReal
							;resk
						]
						,sumReal=I^n SumRealHarmonic[\[Nu]N,AN,xN,yN,\[Alpha],k,sumLim]
						;nc=Module[{h}
							,h=hHarmonic[k][yN,\[Alpha]]
							;If[PossibleZeroQ@h
								,0
								,h GRegHarmonic[d-\[Nu]True,yN,k,\[Alpha]]
							]
						]
						;sumReci=(-1)^k (SumFourierRegHarmonic[\[Nu]N,AN,xN,yN,\[Alpha],k,sumLim]+nc)
						;resk=\[Pi]^(\[Nu]N/2)/Gamma[\[Nu]N/2](sumReal+1/Abs[Det[AN]]sumReci)
						;resk
					]
				]
			]
		,{k,0,Floor[n/2]}
		]
	]
	;res*=\[Mu]^(Total[\[Alpha]]-\[Nu]True)
	;l=Max[0,Round[(\[Nu]True-d)/2]]
	;If[PossibleZeroQ[l-(\[Nu]True-d)/2]
		,res+=Pi^(l+d/2)/Gamma[l+d/2](-1)^(l+1)/l!YDer[l,N[yExt,prec],\[Alpha]]Log[\[Mu] \[Mu]]/Abs[Det[AExt]]
	]
	;res
]


(* ::Subsection::Closed:: *)
(*Create Epstein zeta regularized prototype reference values*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[987];
Do[
	\[Nu]0[i]=10 randRational
	;a0[i]=AHex
	;x0[i]=If[i<10,{0,0},{randRational,randRational}]
	;y0[i]={randRational,randRational}
	;\[Alpha]0[i]={randDerivative,randDerivative}
	,{i,range}
]


(* ::Text:: *)
(*Set special reference values.*)


zDerTabSpecial={
	With[{\[Nu]=1/2,a=IdentityMatrix@2,x={0,0},y={0,0},\[Alpha]={2,0}}
		,res=epsteinZetaRegDer[\[Nu],a,x,y,\[Alpha]]
		;Flatten@Join[N[{\[Nu],a,x,y}],{\[Alpha],N[Re[res]],N[Im[res]]}]
	]
	,With[{\[Nu]=1/2,a=IdentityMatrix@2,x={1/10,0},y={0,0},\[Alpha]={2,0}}
		,res=epsteinZetaRegDer[\[Nu],a,x,y,\[Alpha]]
		;Flatten@Join[N[{\[Nu],a,x,y}],{\[Alpha],N[Re[res]],N[Im[res]]}]
	]
	,With[{\[Nu]=1/2,a=IdentityMatrix@2,x={1/10,0},y={1/20,0},\[Alpha]={2,0}}
		,res=epsteinZetaRegDer[\[Nu],a,x,y,\[Alpha]]
		;Flatten@Join[N[{\[Nu],a,x,y}],{\[Alpha],N[Re[res]],N[Im[res]]}]
	]
	,With[{\[Nu]=1/2,a=AHex,x={1/10,0},y={1/20,0},\[Alpha]={2,0}}
		,res=epsteinZetaRegDer[\[Nu],a,x,y,\[Alpha]]
		;Flatten@Join[N[{\[Nu],a,x,y}],{\[Alpha],N[Re[res]],N[Im[res]]}]
	]
	,With[{\[Nu]=1/2,a=AHex,x={1+1/10,0},y={1/20,0},\[Alpha]={2,0}}
		,res=epsteinZetaRegDer[\[Nu],a,x,y,\[Alpha]]
		;Flatten@Join[N[{\[Nu],a,x,y}],{\[Alpha],N[Re[res]],N[Im[res]]}]
	]
	,With[{\[Nu]=1/2,a=AHex,x={1+1/10,0},y={1/20,2},\[Alpha]={2,0}}
		,res=epsteinZetaRegDer[\[Nu],a,x,y,\[Alpha]]
		;Flatten@Join[N[{\[Nu],a,x,y}],{\[Alpha],N[Re[res]],N[Im[res]]}]
	]
	,With[{\[Nu]=10+1/2,a=IdentityMatrix@2,x={0,0},y={1+1/10,0},\[Alpha]={2,0}}
			,res=epsteinZetaRegDer[\[Nu],a,x,y,\[Alpha]]
			;Flatten@Join[N[{\[Nu],a,x,y}],{\[Alpha],N[Re[res]],N[Im[res]]}]
	]
	,With[{\[Nu]=10+1/2,a=AHex,x={1+1/10,0},y={1/20,2},\[Alpha]={2,0}}
		,res=epsteinZetaRegDer[\[Nu],a,x,y,\[Alpha]]
		;Flatten@Join[N[{\[Nu],a,x,y}],{\[Alpha],N[Re[res]],N[Im[res]]}]
	]
		,With[{\[Nu]=1/2,a=AHex,x={1+1/10,0},y={1/20,2},\[Alpha]={2,0}}
		,res=epsteinZetaRegDer[\[Nu],a,x,y,\[Alpha]]
		;Flatten@Join[N[{\[Nu],a,x,y}],{\[Alpha],N[Re[res]],N[Im[res]]}]
	]
	,With[{\[Nu]=4,a=IdentityMatrix@2,x={0,0},y={0,0},\[Alpha]={2,0}}
		,res=epsteinZetaRegDer[\[Nu],a,x,y,\[Alpha]]
		;Flatten@Join[N[{\[Nu],a,x,y}],{\[Alpha],N[Re[res]],N[Im[res]]}]
	]
	,With[{\[Nu]=6,a=IdentityMatrix@2,x={0,0},y={0,0},\[Alpha]={2,0}}
		,res=epsteinZetaRegDer[\[Nu],a,x,y,\[Alpha]]
		;Flatten@Join[N[{\[Nu],a,x,y}],{\[Alpha],N[Re[res]],N[Im[res]]}]
	]
	,With[{\[Nu]=6,a=AHex,x={0,0},y={0,0},\[Alpha]={2,0}}
		,res=epsteinZetaRegDer[\[Nu],a,x,y,\[Alpha]]
		;Flatten@Join[N[{\[Nu],a,x,y}],{\[Alpha],N[Re[res]],N[Im[res]]}]
	]
};


(* ::Text:: *)
(*Generate random reference values.*)


zDerTabRandom=ParallelTable[
	With[{result=N[epsteinZetaRegDer[\[Nu]0[i],a0[i],x0[i],y0[i],\[Alpha]0[i]],30]}
	,Flatten[Join[N[{\[Nu]0[i],a0[i],x0[i],y0[i]},30],\[Alpha]0[i],{Re[result],Im[result]}]]
	],{i,range}
];


(* ::Text:: *)
(*Join reference values.*)


zDerTab=Join[zDerTabSpecial,zDerTabRandom];


(* ::Text:: *)
(*Save reference values to file.*)


Export[NotebookDirectory[]<>"csv/epsteinZetaRegDer_prototype_ref.csv",zDerTab];


(* ::Subsection::Closed:: *)
(*Create Epstein zeta regularized prototype reference values for nu = dim + 2k.*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[100];
Do[
	\[Nu]Int[i]=2+2 RandomInteger[{0,6}]
	;a0[i]={{1,1/2},{0,Sqrt@3/2}}
	;x0[i]=If[i<10,{0,0},{randRational,randRational}]
	;y0[i]={randRational,randRational}
	;\[Alpha]0[i]={Max[1,randDerivative],randDerivative}
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


zDerd2kTab=ParallelTable[
	With[{result=epsteinZetaRegDer[\[Nu]Int[i],a0[i],x0[i],y0[i],\[Alpha]0[i]]}
	,Flatten[Join[N[{\[Nu]Int[i],a0[i],x0[i],y0[i]},30],\[Alpha]0[i],N[{Re[result],Im[result]},30]]]
	],{i,range}
];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/epsteinZetaRegDer_d2k_prototype_ref.csv",zDerd2kTab];


(* ::Subsection::Closed:: *)
(*Create Epstein zeta regularized derivative reference values for Bain type lattice matrix*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[12];
Do[
	\[Nu]Int[i]=2+2 RandomInteger[{0,6}]
	;A0[i]=RandomReal[{1/10,1}]
	;x0[i]=If[i<10,{0,0,0},{randRational,randRational,randRational}]
	;y0[i]=If[i<5,{0,0,0},{randRational,randRational,randRational}]
	;\[Alpha]0[i]={Max[1,randDerivative],randDerivative,randDerivative}
	,{i,range}
]


range = Range[12];
Do[
	\[Nu]Int[i]=N[2+2 RandomInteger[{0,6}]]
	;a0[i]=N[BMat[RandomReal[{1/10,1}]]]
	;x0[i]=If[i<10,{0.,0.,0.},N@{randRational,randRational,randRational}]
	;y0[i]=If[i<5,{0.,0.,0.},N@{randRational,randRational,randRational}]
	;\[Alpha]0[i]={Max[1,randDerivative],randDerivative,randDerivative}
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


zDerBainTab=ParallelTable[
	With[{result=epsteinZetaRegDer[N[\[Nu]Int[i]],a0[i],x0[i],y0[i],\[Alpha]0[i]]}
	,Flatten[Join[{\[Nu]Int[i],a0[i],x0[i],y0[i]},\[Alpha]0[i],{If[PossibleZeroQ[Re[result]],0.,N[Re[result]]],N[Im[result]]}]]
	],{i,range}
];


SetDirectory[NotebookDirectory[]];
Export["csv/epsteinZetaRegDer_bain_prototype_ref.csv",zDerBainTab];


(* ::Subsection:: *)
(*Create Epstein zeta regularized 1D reference values*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[999];
SeedRandom[2]; (* TODO: cancellation error, seed 2 and up to 40 ders *)
Do[
	\[Nu]1D[i]=10 randRational
	;a1D[i]=1+{{RandomInteger[{-2^10,2^8}]/2^10/2}}
	;x1D[i]=If[i<5,{0},{1+randRational}]
	;y1D[i]=If[i<10,{0},{randRational}]
	;n1D[i]={RandomInteger[{0,20}]}
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


zDerTab1D=ParallelTable[
	With[{result=N[epsteinZetaRegDer[\[Nu]1D[i],a1D[i],x1D[i],y1D[i],n1D[i],InputPrecision->100,SumLim->8]]}
	,Flatten[{N[\[Nu]1D[i]]
	,N[a1D[i]]
	,N[x1D[i]]
	,N[y1D[i]]
	,n1D[i]
	,Re[result]
	,Im[result]
	}]],{i,range}
];


(* ::Text:: *)
(*Save reference values to file.*)


Export[NotebookDirectory[]<>"csv/epsteinZetaRegDer_1D_ref.csv",zDerTab1D];
