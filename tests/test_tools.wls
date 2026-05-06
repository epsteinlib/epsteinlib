#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for tools . c*)


(* ::Subsection::Closed:: *)
(*Setup reference functions*)


(* ::Text:: *)
(*Condition number of an invertible matrix.*)


ConditionNumber[a_]:=Module[
	{sings,kappa}
	,sings=SingularValueList[a]
	;kappa=Max[sings]/Min[sings]
	;kappa
]


(* ::Text:: *)
(*Lattice matrix along the Bain path, a physically relevant test matrix in 3D.*)


BMat[a_]=1/Sqrt[a+1]*DiagonalMatrix[{Sqrt[a],1,1}] . {{1,1,0},{1,0,1},{0,1,1}};


(* ::Subsection:: *)
(*Test Singular Value Decomposition*)


(* ::Text:: *)
(*Test case entries.*)


matrixTestEntry[a_]:=Flatten[{Length[a],N@a,N@Inverse[Transpose[a]]}]


(* ::Text:: *)
(*Generate special 2D reference values.*)


matrixTestData2D={
	matrixTestEntry[{{1,2},{0,1}}]
	,matrixTestEntry[{{1,0},{2,1}}]
	,matrixTestEntry[{{0,2},{1,0}}]
	,matrixTestEntry[{{0,1},{1,1}}]
	,matrixTestEntry[{{1,1},{1,0}}]
	,matrixTestEntry[{{0,2},{1,1}}]
	,matrixTestEntry[{{1,1},{2,0}}]
	,matrixTestEntry[{{0,1},{1,1}}]
};


(* ::Text:: *)
(*Generate random 2D matrix reference values.*)


While[Length[matrixTestData2D]<500
	,temp=RandomReal[{-10,10},{2,2}]
	;If[ConditionNumber[temp]<10
		,AppendTo[matrixTestData2D,matrixTestEntry[temp]]
	]
	;Clear[temp]
]


(* ::Text:: *)
(*Generate special 3D matrix reference values.*)


matrixTestData3D=Join[
	{
		matrixTestEntry[{{0,2,3},{1,1,1},{4,5,6}}]
		,matrixTestEntry[{{0,0,1},{0,1,0},{1,0,0}}]
	},
	Table[
		matrixTestEntry[BMat[a]]
	,{a,1/10,1,1/10}
	]
];


(* ::Text:: *)
(*Generate random 3D matrix reference values.*)


While[Length[matrixTestData3D]<400
	,temp=RandomReal[{-10,10},{3,3}]
	;If[ConditionNumber[temp]<10
		,AppendTo[matrixTestData3D,matrixTestEntry[temp]]
	]
	;Clear[temp]
]


(* ::Text:: *)
(*Generate special 4D matrix reference value.*)


matrixTestData4D=Map[matrixTestEntry,Permutations[
	{{1,0,0,0},{0,1,0,0},{5,-7,1,0},{0,0,0,1}}
]];


(* ::Text:: *)
(*Generate 4D matrix reference values.*)


While[Length[matrixTestData4D]<99
	,temp=RandomReal[{-10,10},{4,4}]
	;If[ConditionNumber[temp]<10
		,AppendTo[matrixTestData4D,matrixTestEntry[temp]]
	]
	;Clear[temp]
]


(* ::Text:: *)
(*Combine reference values.*)


matrixTestData=Join[matrixTestData2D,matrixTestData3D,matrixTestData4D];


(* ::Text:: *)
(*Safe reference values.*)


Export[NotebookDirectory[]<>"csv/svd_ref.csv", matrixTestData];
