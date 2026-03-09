#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for tools . c*)


(* ::Subsection:: *)
(*Setup reference functions*)


(* ::Text:: *)
(*Condition number of an invertible matrix.*)


ConditionNumber[a_]:=Module[
	{sings,kappa}
	,sings=SingularValueList[a]
	;kappa=Max[sings]/Min[sings]
	;kappa
]


(* ::Subsection:: *)
(*Test Singular Value Decomposition*)


(* ::Text:: *)
(*Generate 2D matrix reference values.*)


matrixTestData2D={};
While[Length[matrixTestData2D]<800
	,temp=RandomReal[{-10,10},{2,2}]
	;If[ConditionNumber[temp]<10
		,AppendTo[matrixTestData2D,Flatten[Join[{2,temp,Transpose[Inverse[temp]]}]]]
	]
	;Clear[temp]
]


(* ::Text:: *)
(*Generate 3D matrix reference values.*)


matrixTestData3D={};
While[Length[matrixTestData3D]<100
	,temp=RandomReal[{-10,10},{3,3}]
	;If[ConditionNumber[temp]<10
		,AppendTo[matrixTestData3D,Flatten[Join[{3,temp,Transpose[Inverse[temp]]}]]]
	]
	;Clear[temp]
]


(* ::Text:: *)
(*Generate 4D matrix reference values.*)


matrixTestData4D={};
While[Length[matrixTestData4D]<99
	,temp=RandomReal[{-10,10},{4,4}]
	;If[ConditionNumber[temp]<10
		,AppendTo[matrixTestData4D,Flatten[Join[{4,temp,Transpose[Inverse[temp]]}]]]
	]
	;Clear[temp]
]


(* ::Text:: *)
(*Combine reference data.*)


matrixTestData=Join[matrixTestData2D,matrixTestData3D,matrixTestData4D];


(* ::Text:: *)
(*Safe test cases.*)


SetDirectory[NotebookDirectory[]];
Export["csv/sdv_Ref.csv", matrixTestData];
