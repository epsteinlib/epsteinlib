#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*Unit tests for tools . c*)


(* ::Subsection:: *)
(*Setup r	eference functions*)


(* ::Text:: *)
(*Count  Trailing  Zeros  in  64 - bit  integer  and  return  number  of  consecutive  zero  bits  starting  from  the  least  significant  bit .*)


ctz64Ref[x_] := If[x == 0, 64, IntegerExponent[x, 2]]


(* ::Text:: *)
(*Most  Significant  Bit  position  in  32-bit  integer.*)


msb32Ref[x_] := If[x == 0, 0, Floor[Log2[x]]]


(* ::Text:: *)
(*Convert unsigned 64-bit integer to apint representation.*)
(*	Input: x (unsigned long long), s (sign: \[PlusMinus]1)*)
(*	Output: {x, s, limb[0], limb[1], n, exp2, sign}*)
(*   *)
(*   Algorithm:*)
(*   - Factor out trailing zeros: x = shifted * 2^tz*)
(*   - Split shifted into two 32-bit limbs (base 2^32, little-endian)*)
(*   - Set n = number of nonzero limbs (1 or 2)*)
(*   - Set exp2 = tz (power of 2 factored out)*)
(*   - Canonicalize zero: if x=0 then {n=0, exp2=0, sign=1}*)


toApint[x_, s_] := Module[{tz, shifted, l0, l1, n, finalSign}
	,If[x == 0,
	Return[{x, s, 0, 0, 0, 0, 1}]
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


(* ::Subsection:: *)
(*Test Count Trailing Zeros (64-bit) and Most Significant Bit (32-bit) functions*)


(* ::Text:: *)
(*Setup reference values.*)


SeedRandom[42];
range = Range[995];

Do[
	k0[i] = RandomInteger[{0, 63}]
	;maxM = 2^(64 - k0[i]) - 1
	;m0[i] = 2*RandomInteger[{0, Floor[maxM/2]}] + 1
	;x0[i] = m0[i]*2^k0[i]
	;msbM0[i] = If[m0[i] < 2^32, msb32Ref[m0[i]],32 + msb32Ref[BitShiftRight[m0[i], 32]] ]
	,{i, range}
]

specialCases = {
	{0, 64, 0}
	,{2^63, 63, 0}
	,{1, 0, 0}
	,{2^32 - 1, 0, 31}
};


(* ::Text:: *)
(*Generate  reference  values .*)


ctzmsbTab = ParallelTable[
	{x0[i], ctz64Ref[x0[i]], msbM0[i]}
	,{i, range}
];
ctzmsbTab = Join[specialCases, ctzmsbTab];


(* ::Text:: *)
(*Export reference values.*)


SetDirectory[NotebookDirectory[]];
Export["csv/ctzmsb_Ref.csv", ctzmsbTab];


(* ::Subsection:: *)
(*Test arbitrary precision integer apint*)


(* ::Text:: *)
(*Generate reference  values.*)


SeedRandom[42];
testCases = Join[
  (* Edge cases *)
  {toApint[0, 1], toApint[0, -1], toApint[1, 1], toApint[1, -1],
   toApint[2^32 - 1, 1], toApint[2^32, 1], toApint[2^64 - 1, 1]},

  (* Random cases *)
  ParallelTable[
    toApint[RandomInteger[{0, 2^64 - 1}], RandomChoice[{-1, 1}]],
    {1000}
  ]
];


(* ::Text:: *)
(*Export reference values.*)


Export["csv/apint_set_ull_Ref.csv", testCases]
