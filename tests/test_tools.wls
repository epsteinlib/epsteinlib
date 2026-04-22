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


(* ::Text:: *)
(*Maximum number of limbs, hard-coded here and in C as well.*)


HPDYADMAXLIMBS = 32;


(* ::Text:: *)
(*Count  Trailing  Zeros  in  64 - bit  integer  and  return  number  of  consecutive  zero  bits  starting  from  the  least  significant  bit .*)


ctz64Ref[x_] := If[x == 0, 64, IntegerExponent[x, 2]]


(* ::Text:: *)
(*Most  Significant  Bit  position  in  32-bit  integer.*)


msb32Ref[x_] := If[x == 0, 0, Floor[Log2[x]]]


(* ::Text:: *)
(*Convert unsigned 64-bit integer to high precision dyadic number (hpdyad) representation.*)
(*	Input: x (unsigned long long), s (sign: \[PlusMinus]1)*)
(*	Output: {x, s, limb[0], limb[1], n, exp2, sign}*)


toHpdyad[x_, s_] := Module[{tz, shifted, l0, l1, n, finalSign}
	,If[x == 0
		,Return[{x, s, 0, 0, 0, 0, 1}]
	]
	;tz = IntegerExponent[x, 2]
	;shifted = Quotient[x, 2^tz]
	;l0 = Mod[shifted, 2^32]
	;l1 = Quotient[shifted, 2^32]
	;n = If[l1 != 0, 2, 1]
	;finalSign = s
	(* Output: {x_input, sign_input, limb0, limb1, n, exp2, sign_output} *)
	;{x, s, l0, l1, n, tz, finalSign}
]


(* ::Text:: *)
(*Count leading zeros in 32-bit representation of x. Returns number of zero bits before the most significant 1 and returns 32 if x = 0.*)


clz32[x_]:=If[x == 0,32,31-Floor[Log[2,x]]]


(* ::Text:: *)
(*Shifts array of 32-bit limbs left by s bits. Handles carry propagation between adjacent limbs.Low s bits of limb[i-1] overflow into high positions of limb[i]. Result kept modulo 2^32 per limb.*)


shiftLimbsLeft[limbs_, s_] := Module[{n, result}
	,n = Length[limbs]
	;result = Table[0, {n}]
	;Do[
		result[[i]] = Mod[BitShiftLeft[limbs[[i]], s], 2^32] + If[i > 1, BitShiftRight[limbs[[i-1]], 32 - s], 0]
		,{i, n, 1, -1}
	]
	;result
]


(* ::Text:: *)
(*Normalizes hpdyad representation to canonical form. Removes leading zero limbs from high end. *)
(*	If result is zero: returns sign=+1, exp2=0, all limbs=0.*)
(* Otherwise: left-shifts entire mantissa so MSB of highest non-zero limb is 1. Adjusts exp2 to compensate for the shift.*)
(* 	Returns {new_sign, new_exp2, new_limbs}.*)


hpdyadNormalize[{sign_,exp2_,limbs_}]:=Module[
	{n,s,newLimbs,newExp2}
	,n=Length[limbs]
	;While[n>0&&limbs[[n]]==0,n=n-1]
	;If[n==0,Return[{1,0,Table[0,{HPDYADMAXLIMBS}]}]]
	;s=clz32[limbs[[n]]]
	;If[s==0
		,newLimbs=limbs
		;newExp2=exp2
		,newLimbs=shiftLimbsLeft[limbs,s]
		;newExp2=exp2-s
	]
	;newLimbs=Join[Take[newLimbs,n],Table[0,{HPDYADMAXLIMBS -n}]]
	(* Trim low-end zeros *)
	;While[n>0&&newLimbs[[1]]==0
		,newLimbs=Join[Drop[newLimbs,1],{0}]
		;n=n-1
		;newExp2=newExp2+32
	]
	;{sign,newExp2,newLimbs}
]


(* ::Text:: *)
(*Define multiplication test function. Takes two unsigned long long values and their signs, returns product as double. Simulates: hpdyad_set_ull \[RightArrow] hpdyad_mul \[RightArrow] hpdyad_to_double_trunc.*)


hpdyadMulTest[valA_,signA_,valB_,signB_]:=Module[{prod}
	,prod=signA*valA*signB*valB
	;If[prod==0
	,0.0
	,N[prod]
	]
]


(* ::Text:: *)
(*Convert arbitrary integer to full HPDYADMAXLIMBS-limb hpdyad representation. Returns {limb0, limb1, ..., limb7, exp2, sign}.*)


toHpdyadFull[x_,s_]:=Module[
	{num,den,denExp,tz,shifted,limbs}
	,If[x==0
		,Return[Join[Table[0,{HPDYADMAXLIMBS}],{0,1}]]
	]
	;num=Numerator[x]
	;den=Denominator[x]
	;denExp=IntegerExponent[den,2]
	;tz=IntegerExponent[num,2]-denExp
	;shifted=Quotient[num,2^IntegerExponent[num,2]]
	;limbs=Table[Mod[Quotient[shifted,2^(32*i)],2^32],{i,0,HPDYADMAXLIMBS-1}]
	;Flatten[{limbs,tz,s}]
]


(* ::Text:: *)
(*Simulate hpdyad_mul: multiply two high precision dyadic numbers (hpdyad's) and return normalized result.*)
(*	Input: valA, signA, valB, signB*)
(*	Output: {limb0-7, exp2, sign} after normalization.*)


hpdyadMulTest[valA_,signA_,valB_,signB_]:=Module[
	{prod,prodHpdyad,normalized,limbList}
	,If[valA==0||valB==0
	    ,Return[Join[Table[0,{HPDYADMAXLIMBS}],{0,1}]]
	]
	;prod=signA*valA*signB*valB
	;prodHpdyad=toHpdyadFull[Abs[prod],Sign[prod]]
	;limbList=Take[prodHpdyad,HPDYADMAXLIMBS]
	;normalized=hpdyadNormalize[{prodHpdyad[[HPDYADMAXLIMBS+2]],prodHpdyad[[HPDYADMAXLIMBS+1]],limbList}]
	;Flatten[{normalized[[3]],normalized[[2]],normalized[[1]]}]
]


(* ::Text:: *)
(*Arbitrary precision integer right-shift with sticky bit.*)


hpdyadShrStickyRef[{sign_, exp2_, limbs_}, bits_] := Module[
	{M, Mshifted, sticky, newLimbs, newExp2}
	,M = FromDigits[Reverse[limbs], 2^32]
	;If[M == 0
		,{1, 0, Table[0, {HPDYADMAXLIMBS}]}
		,sticky = If[Mod[M, 2^bits] != 0, 1, 0]
		;Mshifted = BitOr[Floor[M / 2^bits], sticky]
		;newExp2 = exp2 + bits
		;newLimbs = Reverse[IntegerDigits[Mshifted, 2^32, HPDYADMAXLIMBS]]
		;{sign, newExp2, newLimbs}
	]
]


(* ::Text:: *)
(*Convert  integer  to  hpdyad.*)


ratToHpdyad[n_] := Module[
	{raw, sign, absN}
	,absN = Abs[n]
	;sign = If[n == 0, 1, Sign[n]]
	;raw = toHpdyadFull[absN, sign]
	;hpdyadNormalize[{raw[[HPDYADMAXLIMBS+2]], raw[[HPDYADMAXLIMBS+1]], raw[[1;;HPDYADMAXLIMBS]]}]
]


(* ::Subsection::Closed:: *)
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


Export[NotebookDirectory[]<>"csv/svd_Ref.csv", matrixTestData];
