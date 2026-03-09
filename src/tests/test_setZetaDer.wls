#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for epsteinZetaDer in epsteinZeta.c*)


(* ::Subsection:: *)
(*Setup*)


(* ::Text:: *)
(*Load Mathematica wrapper in every kernel. Remove the "Quiet" if unsure about the parallel evaluation.*)


BaseDirectory=FileNameJoin[{NotebookDirectory[],"..", "..", "mathematica"}];
SetDirectory[BaseDirectory];<<"EpsteinZeta.wl";
ParallelEvaluate[SetDirectory[BaseDirectory];Quiet[<<"EpsteinZeta.wl"]];


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
(*Generate random rational numbers and random integers between 0 and 6.*)


randRational := RandomInteger[{-10^(16),10^(16)}]/10^(16)
randRationalStable := RandomChoice[{-1,1}](1/10+8/10 RandomInteger[{0,10^(16)}]/10^(16))
randDerivative := RandomInteger[{0,6}]


(* ::Text:: *)
(*Set Zeta without derivatives.*)


SetZeta[\[Nu]_,a_,x_,y_]:=Exp[2 Pi I x . y]EpsteinZeta[\[Nu],a,x,y]


(* ::Text:: *)
(*Hexagonal lattice.*)


AHex={{1,1/2},{0,Sqrt@3/2}};


(* ::Text:: *)
(*Bain-path lattice matrix.*)


BMat[a_] = 1/Sqrt[a + 1]*DiagonalMatrix[{Sqrt[a], 1, 1}] . {{1, 1, 0}, {1, 0, 1}, {0, 1, 1}};


(* ::Subsection::Closed:: *)
(*Reference functions*)


(* ::Text:: *)
(*(Stable) analytic even derivatives in 1D.*)


SetZetaDerAnalytic[\[Nu]_,{{a_}},{x_},{y_},{n_?EvenQ}]:=(-4 Pi^2)^Floor[n/2] SetZetaDer[\[Nu]-Floor[n],{{a}},{x},{y},{0}]


(* ::Text:: *)
(*Derivatives of the upper Crandall function.*)


p[y_,\[Alpha]_,\[Beta]_]:=Apply[Times,(-\[Pi])^(\[Alpha]-\[Beta])Binomial[\[Alpha],\[Beta]](\[Alpha]-\[Beta])!/(\[Alpha]-2\[Beta])!(2y)^(\[Alpha]-2\[Beta])]
g[\[Nu]_,y_]=Piecewise[{
	{-2/\[Nu],Norm[y]<10^(-64)}
	,{ExpIntegralE[1-\[Nu]/2,Pi y . y],Norm[y]<4}
	},0
];
gLower[\[Nu]_,y_]=Piecewise[{
	{2/\[Nu],Norm[y]<10^(-64)}
	,{Gamma[\[Nu]/2,0,Pi y . y]/(Pi y . y)^(\[Nu]/2),Norm[y]<4}
	}
];
gDer[\[Nu]_,y_,\[Alpha]_]:=Total[Flatten[Table[
	Module[
		{pol}
		,pol=p[y,\[Alpha],\[Beta]]
		;If[Abs[pol]<10^-60
			,0
			,pol g[\[Nu]+2Total[\[Alpha]-\[Beta]],y]
		]
	]
	,{\[Beta],Tuples[Range[ConstantArray[0,Length[\[Alpha]]],\[Alpha]/2]]}
]]]


(* ::Text:: *)
(*Derivatives of sums in real and reciprocal space.*)


SumRealDer[\[Nu]_,A_,x_,y_,\[Alpha]_,sumLim_]:=Total@Table[Apply[Times,(-2Pi I z)^\[Alpha]]g[\[Nu],z]Exp[-2Pi I y . z],{z,(A . #-x)&/@Tuples[Range[-sumLim,sumLim],Length[A]]}]
SumFourierDer[\[Nu]_,A_,x_,y_,\[Alpha]_,sumLim_]:=Total@Table[gDer[Length[\[Alpha]]-\[Nu],k+y,\[Alpha]]Exp[-2Pi I x . k],{k,(Inverse[Transpose[A]] . #)&/@Tuples[Range[-sumLim,sumLim],Length[A]]}]


(* ::Text:: *)
(*Prototype for epsteinZetaDer without special cases such non-positive even integers nu.*)


Options[setZetaDer]={InputPrecision->100,SumLim->10};
setZetaDer[\[Nu]Ext_,AExt_,xExt_,yExt_,\[Alpha]_,OptionsPattern[]]:=Module[
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
		,sumReal=SetPrecision[SumRealDer[\[Nu]N,AN,xN,yN,\[Alpha],sumLim],prec]
		;sumReci=SetPrecision[SumFourierDer[\[Nu]N,AN,xN,yN,\[Alpha],sumLim],prec]
		;res=SetPrecision[\[Mu]^(Total[\[Alpha]]-\[Nu]N)\[Pi]^(\[Nu]N/2)/Gamma[\[Nu]N/2](sumReal+1/Abs[Det[AN]]sumReci),prec]
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
(*Sum in real and reciprocal space.*)


SumRealHarmonic[\[Nu]_,A_,x_,y_,\[Alpha]_,kIndex_,sumLim_]:=Total[Table[
	Module[{h}
		,h=hHarmonic[kIndex][z,\[Alpha]]
		;If[PossibleZeroQ[h],0,h g[\[Nu],z]Exp[-2\[Pi] I y . z]]
	]
	,{z,(A . #-x)&/@Tuples[Range[-sumLim,sumLim],Length[A]]}
]]
SumFourierHarmonic[\[Nu]_,A_,x_,y_,\[Alpha]_,kIndex_,sumLim_]:=Total[Table[
	hHarmonic[kIndex][y+k,\[Alpha]]g[Length[x]-\[Nu]+2Total[\[Alpha]]-4kIndex,k+y]Exp[-2\[Pi] I x . k]
	,{k,(Inverse[Transpose[A]] . #)&/@Tuples[Range[-sumLim,sumLim],Length[A]]}
]]


(* ::Text:: *)
(*Prototype for epsteinZetaDer with harmonic functions without special cases such non-positive even integers nu.*)


Options[setZetaDerHarmonic]={InputPrecision->100,SumLim->10};
setZetaDerHarmonic[\[Nu]Ext_,AExt_,xExt_,yExt_,\[Alpha]_,OptionsPattern[]]:=Module[
	{EPS,prec,sumLim,\[Nu]True,AN,xN,yN,\[Mu],res,n}
	,EPS=10^-50
	;prec=OptionValue[InputPrecision]
	;sumLim=OptionValue[SumLim]
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
				;If[\[Nu]N<-1&&Abs[\[Nu]N/2-Round[\[Nu]N/2]]<EPS
					,0
					,If[Abs[\[Nu]N]<EPS
						,If[Norm[xN]>EPS
							,0
							,sumReal=-I^n hHarmonic[k][xN,\[Alpha]]
							;resk=\[Pi]^(\[Nu]N/2)sumReal
							;resk*=(-2Pi)^n
							;resk
						]
						,sumReal=I^n SumRealHarmonic[\[Nu]N,AN,xN,yN,\[Alpha],k,sumLim]
						;sumReci=(-1)^k SumFourierHarmonic[\[Nu]N,AN,xN,yN,\[Alpha],k,sumLim]
						;resk=\[Pi]^(\[Nu]N/2)/Gamma[\[Nu]N/2](sumReal+1/Abs[Det[AN]]sumReci)
						;resk
					]
				]
			]
		,{k,0,Floor[n/2]}
		]
	]
	;\[Mu]^(Total[\[Alpha]]-\[Nu]True)res
]


(* ::Text:: *)
(*Isolate singularities in real sum.*)


SumRealHarmonicHighExp[\[Nu]_,A_,x_,y_,\[Alpha]_,kIndex_,sumLim_]:=Total[Table[
	Module[{h,crandall}
		,h=hHarmonic[kIndex][z,\[Alpha]]
		;If[
			PossibleZeroQ[h]
			,0
			,h g[\[Nu],z] Exp[-2\[Pi] I y . z]
		]
	]
	,{z,(A . #-x)&/@Select[Tuples[Range[-sumLim,sumLim],Length[A]],Not[MemberQ[{0,1},Norm[#]]]&]}
]]+Total[Table[
	Module[{h,crandall}
		,h=hHarmonic[kIndex][z,\[Alpha]]
		;If[PossibleZeroQ[h]
			,0
			,h (-gLower[\[Nu],z])Exp[-2\[Pi] I y . z]
		]
	]
	,{z,(A . #-x)&/@Select[Tuples[Range[-sumLim,sumLim],Length[A]],MemberQ[{0,1},Norm[#]]&]}
]]


(* ::Text:: *)
(*Prototype for large nu values.*)


Options[setZetaDerHighExp]={InputPrecision->100,SumLim->10,EPSCANCELLATION->10^-16};
setZetaDerHighExp[\[Nu]Ext_,AExt_,xExt_,yExt_,\[Alpha]_,OptionsPattern[]]:=Module[
	{EPS,EPSCANC,prec,sumLim,\[Nu]True,AN,xN,yN,\[Mu],res,n}
	,EPS=10^-50
	;EPSCANC=OptionValue[EPSCANCELLATION]
	;prec=OptionValue[InputPrecision]
	;sumLim=OptionValue[SumLim]
	;n=Total[\[Alpha]]
	;\[Nu]True=N[\[Nu]Ext,prec]
	;AN=N[AExt,prec]
	;\[Mu]=Abs[Det[AN]]^(1/Length[AN])
	;AN=1/\[Mu] N[AExt,prec]
	;xN=1/\[Mu] N[xExt,prec]
	;yN=\[Mu] N[yExt,prec]
	;If[Abs[\[Nu]Ext]<EPS&&n>0
		,res=0
		,res=(-2Pi)^n Total@Table[
			Module[
				{\[Nu]N,sumReal,sumReci,resk}
				,\[Nu]N=\[Nu]True-2k
				;If[\[Nu]N<-1&&Abs[\[Nu]N/2-Round[\[Nu]N/2]]<EPS
					,0
					,If[Abs[\[Nu]N]<EPS
						,If[Norm[xN]>EPS
							,0
							,resk=-hHarmonic[k][xN,\[Alpha]]
							;resk*=(-2 Pi I)^n \[Pi]^(\[Nu]N/2)
							;resk
						]
						,sumReal=I^n SumRealHarmonicHighExp[\[Nu]N,AN,xN,yN,\[Alpha],k,sumLim]
						;sumReci=(-1)^k SumFourierHarmonic[\[Nu]N,AN,xN,yN,\[Alpha],k,sumLim]
						;resk=sumReal+1/Abs[Det[AN]]sumReci
						;resk*=\[Pi]^(\[Nu]N/2)/Gamma[\[Nu]N/2]
						;resk
					]
				]
			]
			,{k,0,Floor[n/2]}
		]
		;res+=(-2 Pi I)^n Total[Table[
			Module[{resk}
				,resk=Total[Table[
					If[(\[Nu]True-2k)<-1&&Abs[(\[Nu]True-2k)/2-Round[(\[Nu]True-2k)/2]]<EPS
						,0
						,If[Abs[\[Nu]True-2k]<EPS
							,0
							,Module[
								{h}
								,h=hHarmonic[k][z,\[Alpha]]
								;If[PossibleZeroQ[h]||PossibleZeroQ[Norm[z]]
									,0
									,h (z . z)^k
								]
							]
						]
					]
					,{k,0,Floor[n/2]}
				]]
				;If[Abs[resk]<EPSCANC
					,0
					,resk Exp[-2\[Pi] I yN . z](z . z)^(-\[Nu]True/2)
				]
			]
			,{z,(AN . #-xN)&/@Select[Tuples[Range[-sumLim,sumLim],Length[AN]],MemberQ[{0,1},Norm[#]]&]}
		]]
	]
	;\[Mu]^(Total[\[Alpha]]-\[Nu]True)res
]


(* ::Subsection::Closed:: *)
(*Create Epstein zeta 1D analytic reference values in x = 0*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[100];
SeedRandom[2]; (* TODO: cancellation error, seed 2 and up to 40 ders *)
Do[
	\[Nu]1D[i]=10 randRational
	;a1D[i]={{1/2+i/100}}
	;x1D[i]=If[i<5,{0},{randRational}]
	;y1D[i]=If[i<10,{0},{randRational}]
	;n1D[i]={RandomInteger[{0,40}]}
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


zDerTab1D=ParallelTable[
	With[{result=N[setZetaDer[\[Nu]1D[i],a1D[i],x1D[i],y1D[i],n1D[i],InputPrecision->100,SumLim->8]]}
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


SetDirectory[NotebookDirectory[]];
Export["csv/setZetaDer_1D_Ref.csv",zDerTab1D];


(* ::Subsection::Closed:: *)
(*Create set zeta 2D reference values*)


(* ::Text:: *)
(*Restrict x and y to their lattice cells.*)


randomElemLatticeVector[a_]:=a . Table[RandomInteger[{-10^(16)/2,10^(16)/2}]/10^(16),{Length[a]}]


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[186];
SeedRandom[42];
AHex={{1,1/2},{0,Sqrt[3]/2}};
Do[
	\[Nu]0[i]=10 randRational
	;a0[i]=AHex
	;x0[i]=If[i<10,{0,0},randomElemLatticeVector[AHex]]
	;y0[i]=randomElemLatticeVector[Inverse[Transpose[AHex]]]
	;\[Alpha]0[i]={RandomInteger[{0,15}],RandomInteger[{0,15}]}
	,{i,range}
]


(* ::Text:: *)
(*Generate random reference values.*)


zDerRandomTab=ParallelTable[
	Flatten@{N@\[Nu]0[i],N@a0[i],N@x0[i],N@y0[i],\[Alpha]0[i],{N@Re@#,N@Im@#}&/@{setZetaDer[\[Nu]0[i],a0[i],x0[i],y0[i],\[Alpha]0[i]]}}
	,{i,range}
];


(* ::Text:: *)
(*Append specific test values.*)


zDerSpecialTab={
	Flatten@{N@1/2,N@IdentityMatrix@2,{0.,0.},{0.,0.},{2,0},N@setZetaDer[1/2,IdentityMatrix@2,{0,0},{0,0},{2,0}],0.}
	,Flatten@{N@1/2,N@IdentityMatrix@2,{0.,0.},{.1,.2},{2,0},{N@Re@#,N@Im@#}&/@{setZetaDer[1/2,IdentityMatrix@2,{0,0},{1,2}/10,{2,0}]}}
	,Flatten@{N@1/2,N@IdentityMatrix@2,{.1,.2},{0.,0.},{2,0},{N@Re@#,N@Im@#}&/@{setZetaDer[1/2,IdentityMatrix@2,{1,2}/10,{0,0},{2,0}]}}
	,Flatten@{N@1/2,N@AHex,{.1,.2},{0.,0.},{2,0},{N@Re@#,N@Im@#}&/@{setZetaDer[1/2,AHex,{1,2}/10,{0,0},{2,0}]}}
	,Flatten@{N@1/2,N@AHex,{0.,0.},{.1,0.3},{2,0},{N@Re@#,N@Im@#}&/@{setZetaDer[1/2,AHex,{0.,0.},{1,3}/10,{2,0}]}}
	,Flatten@{N@1/2,N@AHex,{.1,.2},{.1,0.3},{2,0},{N@Re@#,N@Im@#}&/@{setZetaDer[1/2,AHex,{1,2}/10,{1,3}/10,{2,0}]}}
	,Flatten@{N@1/2,N@IdentityMatrix@2,{3.,4.},{0.,0.},{2,0},{N@Re@#,N@Im@#}&/@{setZetaDer[1/2,IdentityMatrix@2,{3,4},{0,0},{2,0}]}}
	,Flatten@{N@1/2,N@IdentityMatrix@2,{0.,0.},{3.,4.},{2,0},{N@Re@#,N@Im@#}&/@{setZetaDer[1/2,IdentityMatrix@2,{0,0},{3,4},{2,0}]}}
	,Flatten@{N@1/2,N@AHex,{1.1,3.2},{0.,0.},{2,0},{N@Re@#,N@Im@#}&/@{setZetaDer[1/2,AHex,{1,3}+{1,2}/10,{0,0},{2,0}]}}
	,Flatten@{N@1/2,N@AHex,{0.,0.},{1.1,3.2},{2,0},{N@Re@#,N@Im@#}&/@{setZetaDer[1/2,AHex,{0,0},{1,3}+{1,2}/10,{2,0}]}}
	,Flatten@{N@1/2,N@IdentityMatrix@2,1.+{.1,.2},{.0,.0},{2,0},{N@Re@#,N@Im@#}&/@{setZetaDer[1/2,IdentityMatrix@2,1+{1,2}/10,{0,0},{2,0}]}}
	,Flatten@{N@1/2,N@IdentityMatrix@2,1.+{.1,.2},{.1,.0},{2,0},{N@Re@#,N@Im@#}&/@{setZetaDer[1/2,IdentityMatrix@2,1+{1,2}/10,{1,0}/10,{2,0}]}}
	,Flatten@{N@1/2,N@IdentityMatrix@2,1.+{.1,.2},3+{.1,.0},{2,0},{N@Re@#,N@Im@#}&/@{setZetaDer[1/2,IdentityMatrix@2,1+{1,2}/10,3+{1,0}/10,{2,0}]}}
	,Flatten@{N@1/2,N@IdentityMatrix@2,1.+{.1,.2},3+{.1,.0},{2,4},{N@Re@#,N@Im@#}&/@{setZetaDer[1/2,IdentityMatrix@2,1+{1,2}/10,3+{1,0}/10,{2,4}]}}
};


(* ::Text:: *)
(*Add reference values for large nu.*)


zDerLargeNu =Flatten[
	ParallelTable[{
			Flatten@{N@{\[Nu],{{1,0},{0,1}},{0,0},{0,0}},\[Alpha],{N@Re@#,N@Im@#}&/@{setZetaDer[\[Nu],{{1,0},{0,1}},{0,0},{0,0},\[Alpha]]}}
			,Flatten@{N@{\[Nu],{{1,0},{0,1}},{0.,0.},{0,.4}},\[Alpha],{N@Re@#,N@Im@#}&/@{setZetaDer[\[Nu],{{1,0},{0,1}},{0,0},{0,4/10},\[Alpha]]}}
			,Flatten@{N@{\[Nu],{{1,0},{0,1}},{0.,0.4},{0,.0}},\[Alpha],{N@Re@#,N@Im@#}&/@{setZetaDer[\[Nu],{{1,0},{0,1}},{0,4/10},{0,0},\[Alpha]]}}
			,Flatten@{N@{\[Nu],{{1,0},{0,1}},{0.,0.2},{0,.4}},\[Alpha],{N@Re@#,N@Im@#}&/@{setZetaDer[\[Nu],{{1,0},{0,1}},{0,2/10},{0,4/10},\[Alpha]]}}
		}
		,{\[Nu],{1/2,30,83/3,11}}
		,{\[Alpha],{{2,0},{6,8},{6,9}}}
	]
,2];
(*
Note: we observe catastrophic errors for large nu and x and y not in elementary lattice cell.
This seems to be a floating point arithmetic error, where the small error of projecting x gets amplified.
To test this, add above:
	,Flatten@{N@{\[Nu],{{1,0},{0,1}},{0.,0.2},{0,1.3}},\[Alpha],{N@Re@#,N@Im@#}&/@{setZetaDer[\[Nu],{{1,0},{0,1}},{0,2/10},1+{0,3/10},\[Alpha]]}}
	,Flatten@{N@{\[Nu],{{1,0},{0,1}},{0.,1.2},{0,.3}},\[Alpha],{N@Re@#,N@Im@#}&/@{setZetaDer[\[Nu],{{1,0},{0,1}},1+{0,2/10},{0,3/10},\[Alpha]]}}
	,Flatten@{N@{\[Nu],{{1,0},{0,1}},{0.,1.2},{0,1.3}},\[Alpha],{N@Re@#,N@Im@#}&/@{setZetaDer[\[Nu],{{1,0},{0,1}},1+{0,2/10},1+{0,3/10},\[Alpha]]}}
*)


(* ::Text:: *)
(*Combine reference values.*)


zDerTab=Join[zDerRandomTab,zDerSpecialTab,zDerLargeNu];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/setZetaDer_2D_Ref.csv",zDerTab];


(* ::Subsection::Closed:: *)
(*Create set zeta derivative reference values for odd derivatives*)


(* ::Text:: *)
(*We employ, that SetZetaDer[\[Nu],A,x,y,\[Alpha]] should be 0 for x\:1d62 = y\:1d62 = 0 and \[Alpha]\:1d62 is odd for some 1 <= i <= d.*)


(* does not yet work for {{1,1,0},{1,0,1},{0,1,1}} (here |Det(a)|=2) or derivatives > 5*)
range = Range[100];
Do[
	\[Nu]Odd[i]=10 randRational
	;aOdd[i]=IdentityMatrix@3
	;{xOdd[i],yOdd[i],\[Alpha]Odd[i]}=Transpose[RandomSample[{{0,0,2 RandomInteger[{0,2}]+1},{randRational,randRational,RandomInteger[{0,2}]},{randRational,randRational,RandomInteger[{0,2}]}}]]
	,{i,range}
]


zOddTab=ParallelTable[
	Flatten[{
	N[\[Nu]Odd[i]]
	,N[aOdd[i]]
	,N[xOdd[i]]
	,N[yOdd[i]]
	,\[Alpha]Odd[i]
	,0.
	,0.
	}],{i,range}
];


SetDirectory[NotebookDirectory[]];
Export["csv/setZetaDer_odd_Ref.csv",zOddTab];


(* ::Subsection::Closed:: *)
(*Create Laplacian derivative set zeta reference values*)


(* ::Text:: *)
(*Setup test parameters.*)


range = Range[100];
Do[
	\[Nu]L[i]=10 randRational
	;aL[i]={{1,1/2},{0,Sqrt@3/2}}
	;xL[i]=If[i<10,{0,0},{randRational,randRational}]
	;yL[i]=aL[i] . {randRationalStable,randRationalStable}
	,{i,range}
]


(* ::Text:: *)
(*Generate reference values.*)


zLaplaceTab=ParallelTable[
	With[{result=(-2 Pi I)^2 SetZeta[\[Nu]L[i]-2,aL[i],xL[i],yL[i]]}
	,{N[\[Nu]L[i]]
	,N[aL[i][[1,1]]]
	,N[aL[i][[1,2]]]
	,N[aL[i][[2,1]]]
	,N[aL[i][[2,2]]]
	,N[xL[i][[1]]]
	,N[xL[i][[2]]]
	,N[yL[i][[1]]]
	,N[yL[i][[2]]]
	,Re[result]
	,Im[result]
	}],{i,range}
];


(* ::Text:: *)
(*Save reference values to file.*)


SetDirectory[NotebookDirectory[]];
Export["csv/setZetaDer_laplace_Ref.csv",zLaplaceTab];


(* ::Subsection::Closed:: *)
(*Check behavior in critical nu values*)


(* ::Text:: *)
(*Generate reverence values for nu = dim + |\[Alpha]| + even natural.*)


smallVal = 10^-200;
zetaSpecialNuBiggerDim = {
	(* 1D: nu = 1 + Total[alpha] + 2Total[alpha] even => nu odd *)
	Flatten@{1,5.,N@{{1}},{0.},{0.},{2},{N@Re@#,N@Im@#}  &/@ {setZetaDer[5,{{1}},{0},{0},{2}]}}
    ,Flatten@{1,7.,N@{{1}},{0.},{0.},{4},{N@Re@#,N@Im@#}  &/@ {setZetaDer[7,{{1}},{0},{0},{4}]}}
    ,Flatten@{1,5.,N@{{1}},{0.1},{0.},{2},{N@Re@#,N@Im@#}  &/@ {setZetaDer[5,{{1}},{1/10},{0},{2}]}}
    ,Flatten@{1,7.,N@{{1}},{0.1},{0.},{4},{N@Re@#,N@Im@#}  &/@ {setZetaDer[7,{{1}},{1/10},{0},{4}]}}
    ,Flatten@{1,5.,N@{{1}},{0.},{N@smallVal},{2},{N@Re@#,N@Im@#}  &/@ {setZetaDer[5,{{1}},{0},{smallVal},{2}]}}
    ,Flatten@{1,5.,N@{{1/2}},{0.},{0.},{2},{N@Re@#,N@Im@#}  &/@ {setZetaDer[5,{{1/2}},{0},{0},{2}]}}
    ,Flatten@{1,5.,N@{{2}},{0.},{0.},{2},{N@Re@#,N@Im@#}  &/@ {setZetaDer[5,{{2}},{0},{0},{2}]}}
    (* 2D: nu = 2 + Total[alpha] + 2Total[alpha] even => nu even *)
	,Flatten@{2,6.,N@IdentityMatrix@2,{0.,0.},{0.,0.},{2,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[6,IdentityMatrix@2,{0,0},{0,0},{2,0}]}}
    ,Flatten@{2,6.,N@IdentityMatrix@2,{0.,0.},{0.,0.},{1,1},{N@Re@#,N@Im@#}  &/@ {setZetaDer[6,IdentityMatrix@2,{0,0},{0,0},{1,1}]}}
    ,Flatten@{2,6.,N@IdentityMatrix@2,{0.,0.},{0.,0.},{2,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[6,IdentityMatrix@2,{0,0},{0,0},{2,0}]}}
    ,Flatten@{2,8.,N@IdentityMatrix@2,{0.,0.},{0.,0.},{4,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[8,IdentityMatrix@2,{0,0},{0,0},{4,0}]}}
    ,Flatten@{2,8.,N@IdentityMatrix@2,{0.,0.},{0.,0.},{2,2},{N@Re@#,N@Im@#}  &/@ {setZetaDer[8,IdentityMatrix@2,{0,0},{0,0},{2,2}]}}
    ,Flatten@{2,6.,N@IdentityMatrix@2,{0.1,0.2},{0.,0.},{2,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[6,IdentityMatrix@2,{1,2}/10,{0,0},{2,0}]}}
    ,Flatten@{2,6.,N@IdentityMatrix@2,{0.1,0.2},{0.,0.},{1,1},{N@Re@#,N@Im@#}  &/@ {setZetaDer[6,IdentityMatrix@2,{1,2}/10,{0,0},{1,1}]}}
    ,Flatten@{2,8.,N@IdentityMatrix@2,{0.1,0.2},{0.,0.},{2,2},{N@Re@#,N@Im@#}  &/@ {setZetaDer[8,IdentityMatrix@2,{1,2}/10,{0,0},{2,2}]}}
    ,Flatten@{2,6.,N@IdentityMatrix@2,{0.,0.},{0.,N@smallVal},{2,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[6,IdentityMatrix@2,{0,0},{0,smallVal},{2,0}]}}
    ,Flatten@{2,6.,N@IdentityMatrix@2,{0.,0.},{N@smallVal,0.},{2,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[6,IdentityMatrix@2,{0,0},{smallVal,0},{2,0}]}}
    ,Flatten@{2,6.,N@AHex,{0.,0.},{0.,0.},{2,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[6,AHex,{0,0},{0,0},{2,0}]}}
    ,Flatten@{2,6.,N@AHex,{0.,0.},{0.,0.},{1,1},{N@Re@#,N@Im@#}  &/@ {setZetaDer[6,AHex,{0,0},{0,0},{1,1}]}}
    ,Flatten@{2,8.,N@AHex,{0.,0.},{0.,0.},{2,2},{N@Re@#,N@Im@#}  &/@ {setZetaDer[8,AHex,{0,0},{0,0},{2,2}]}}
    ,Flatten@{2,6.,N@AHex,{0.1,0.2},{0.,0.},{2,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[6,AHex,{1,2}/10,{0,0},{2,0}]}}
    ,Flatten@{2,6.,N@AHex,{0.,0.},{0.,N@smallVal},{2,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[6,AHex,{0,0},{0,smallVal},{2,0}]}}
    ,Flatten@{2,8.,N@AHex,{0.,0.},{0.,0.},{4,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[8,AHex,{0,0},{0,0},{4,0}]}}
    ,Flatten@{2,8.,N@AHex,{0.,0.},{0.,0.},{3,1},{N@Re@#,N@Im@#}  &/@ {setZetaDer[8,AHex,{0,0},{0,0},{3,1}]}}
    ,Flatten@{2,10.,N@AHex,{0.,0.},{0.,0.},{4,2},{N@Re@#,N@Im@#}  &/@ {setZetaDer[10,AHex,{0,0},{0,0},{4,2}]}}
    ,Flatten@{2,8.,N@AHex,{0.1,0.2},{0.,0.},{4,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[8,AHex,{1,2}/10,{0,0},{4,0}]}}
    ,Flatten@{2,8.,N@AHex,{0.,0.},{0.,N@smallVal},{4,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[8,AHex,{0,0},{0,smallVal},{4,0}]}}
	(* 3D: nu = 3 + Total[alpha] + 2Total[alpha] even => nu odd *)
	,Flatten@{3,7.,N@IdentityMatrix@3,{0.,0.,0.},{0.,0.,0.},{2,0,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[7,IdentityMatrix@3,{0,0,0},{0,0,0},{2,0,0}]}}
    ,Flatten@{3,7.,N@IdentityMatrix@3,{0.,0.,0.},{0.,0.,0.},{1,1,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[7,IdentityMatrix@3,{0,0,0},{0,0,0},{1,1,0}]}}
    ,Flatten@{3,7.,N@IdentityMatrix@3,{0.,0.,0.},{0.,0.,0.},{2,0,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[7,IdentityMatrix@3,{0,0,0},{0,0,0},{2,0,0}]}}
    ,Flatten@{3,7.,N@IdentityMatrix@3,{0.,0.,0.},{0.,0.,0.},{2,0,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[7,IdentityMatrix@3,{0,0,0},{0,0,0},{2,0,0}]}}
    ,Flatten@{3,9.,N@IdentityMatrix@3,{0.,0.,0.},{0.,0.,0.},{4,0,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[9,IdentityMatrix@3,{0,0,0},{0,0,0},{4,0,0}]}}
    ,Flatten@{3,9.,N@IdentityMatrix@3,{0.,0.,0.},{0.,0.,0.},{2,2,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[9,IdentityMatrix@3,{0,0,0},{0,0,0},{2,2,0}]}}
    ,Flatten@{3,7.,N@IdentityMatrix@3,{0.1,0.2,0.3},{0.,0.,0.},{2,0,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[7,IdentityMatrix@3,{1,2,3}/10,{0,0,0},{2,0,0}]}}
    ,Flatten@{3,7.,N@IdentityMatrix@3,{0.1,0.2,0.3},{0.,0.,0.},{1,1,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[7,IdentityMatrix@3,{1,2,3}/10,{0,0,0},{1,1,0}]}}
    ,Flatten@{3,7.,N@IdentityMatrix@3,{0.,0.,0.},{0.,0.,N@smallVal},{2,0,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[7,IdentityMatrix@3,{0,0,0},{0,0,smallVal},{2,0,0}]}}
    ,Flatten@{3,7.,N@IdentityMatrix@3,{0.,0.,0.},{N@smallVal,0.,0.},{2,0,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[7,IdentityMatrix@3,{0,0,0},{smallVal,0,0},{2,0,0}]}}
    ,Flatten@{3,7.,N@BMat[1/10],{0.,0.,0.},{0.,0.,0.},{2,0,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[7,BMat[1/10],{0,0,0},{0,0,0},{2,0,0}]}}
    ,Flatten@{3,7.,N@BMat[1/2],{0.,0.,0.},{0.,0.,0.},{2,0,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[7,BMat[1/2],{0,0,0},{0,0,0},{2,0,0}]}}
    ,Flatten@{3,7.,N@BMat[1/2],{0.1,0.2,0.3},{0.,0.,0.},{2,0,0},{N@Re@#,N@Im@#}  &/@ {setZetaDer[7,BMat[1/2],{1,2,3}/10,{0,0,0},{2,0,0}]}}
};


(* ::Text:: *)
(*Generate nu = dim - even natural reference values.*)


zetaSpecialNuEqualsDimMinus2k=Flatten[
	ParallelTable[
		Flatten@{d,N@nu,N@IdentityMatrix@d,N@PadRight[{xoffset},d],N[PadRight[{yoffset},d]],n UnitVector[d,1],{N@Re@#,N@Im@#}  &/@ {setZetaDerHarmonic[nu,IdentityMatrix@d,PadRight[{xoffset},d],PadRight[{yoffset},d],n UnitVector[d,1]]}}
		,{d,{1,2}}
		,{yoffset,{0,smallVal,2^-20,1/10}}
		,{xoffset,{0,smallVal,2^-20,1/10}}
		,{n,{2,4}}
		,{nu,-10,2,1}
	]
,4];


(* ::Text:: *)
(*Combine reference values.*)


zetaSpecialTab = Join[zetaSpecialNuBiggerDim,zetaSpecialNuEqualsDimMinus2k];


(* ::Text:: *)
(*Export reference values.*)


SetDirectory[NotebookDirectory[]];
Export["csv/setZetaDer_special_exponents_Ref.csv", zetaSpecialTab];
