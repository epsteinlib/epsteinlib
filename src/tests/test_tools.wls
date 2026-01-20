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


toApint[x_, s_] := Module[{tz, shifted, l0, l1, n, finalSign}
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
(*Normalizes apint representation to canonical form.Removes leading zero limbs from high end. *)
(*	If result is zero: returns sign=+1, exp2=0, all limbs=0.*)
(* Otherwise: left-shifts entire mantissa so MSB of highest non-zero limb is 1. Adjusts exp2 to compensate for the shift.*)
(* 	Returns {new_sign, new_exp2, new_limbs}.*)


apintNormalize[{sign_, exp2_, limbs_}] := Module[
{n, s, newLimbs, newExp2}
,n = Length[limbs]
;While[n > 0 && limbs[[n]] == 0, n = n - 1]
;If[n == 0, Return[{1, 0, Table[0, {8}]}]]
;s = clz32[limbs[[n]]]
;If[s == 0
,newLimbs = limbs
;newExp2 = exp2
,newLimbs = shiftLimbsLeft[limbs, s]
;newExp2 = exp2 - s
]
;newLimbs = Join[Take[newLimbs, n], Table[0, {8 - n}]]
;{sign, newExp2, newLimbs}
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
    {992}
  ]
];


(* ::Text:: *)
(*Export reference values.*)


SetDirectory[NotebookDirectory[]];
Export["csv/apint_set_ull_Ref.csv", testCases];


(* ::Subsection:: *)
(*Test arbitrary precision integer normalization function*)


(* ::Text:: *)
(*Generate test inputs.*)


SeedRandom[42];
testCases = {
	(* 1. Zero *)
	{1, 0, {0, 0, 0, 0, 0, 0, 0, 0}}
	,{-1, 5, {0, 0, 0, 0, 0, 0, 0, 0}}
	,(* 2. Single limb, needs max shift (limb = 1) *)
	{1, 0, {1, 0, 0, 0, 0, 0, 0, 0}}
	,{-1, 10, {1, 0, 0, 0, 0, 0, 0, 0}}
	,(* 3. Already normalized (MSB set) *)
	{1, 0, {BitShiftLeft[1, 31], 0, 0, 0, 0, 0, 0, 0}}
	,{1, 5, {0, BitShiftLeft[1, 31], 0, 0, 0, 0, 0, 0}}
	,(* 4. Middle zeros *)
	{1, 0, {5, 0, 0, 7, 0, 0, 0, 0}}
	,{-1, -3, {100, 0, 200, 0, 0, 0, 0, 0}}
	,(* 5. Multiple leading zeros to trim *)
	{1, 0, {42, 0, 0, 0, 0, 0, 0, 0}}
	,{1, 0, {1234567, 9876543, 0, 0, 0, 0, 0, 0}}
	,(* 6. Needs various shift amounts *)
	{1, 0, {BitShiftLeft[1, 15], 0, 0, 0, 0, 0, 0, 0}}
	,{1, 0, {BitShiftLeft[1, 30], 0, 0, 0, 0, 0, 0, 0}}
	,{1, 0, {0, BitShiftLeft[1, 7], 0, 0, 0, 0, 0, 0}}
	,(* 7. Full mantissa with middle zeros *)
	{1, 100, {16^^FFFFFFFF, 0, 16^^12345678, 0, 16^^ABCDEF00, 0, 0, 16^^80000000}}
	,(* 8. Dangerous: shift causes carry between limbs *)
	{1, 0, {16^^FFFF0000, 16^^0000FFFF, 0, 0, 0, 0, 0, 0}}
	,{1, 0, {16^^80000001, 16^^00000001, 0, 0, 0, 0, 0, 0}}
	,(* 9. Large exp2 values *)
	{1, 1000000, {123, 0, 0, 0, 0, 0, 0, 0}}
	,{-1, -500000, {16^^7FFFFFFF, 0, 0, 0, 0, 0, 0, 0}}
	,(* 10. Random cases *)
	{1, RandomInteger[{-1000, 1000}], RandomInteger[{0, 16^^FFFFFFFF}, 8]}
	,{-1, RandomInteger[{-1000, 1000}], RandomInteger[{0, 16^^FFFFFFFF}, 8]}
	,{1, RandomInteger[{-1000, 1000}], RandomInteger[{0, 16^^FFFFFFFF}, 8]}
	,{-1, RandomInteger[{-1000, 1000}], RandomInteger[{0, 16^^FFFFFFFF}, 8]}
};


(* ::Text:: *)
(*Generate reference values.*)


results = Table[
	Module[
		{input, output}
		,input = testCases[[i]]
		;output = apintNormalize[input]
		;Join[
			{input[[1]]}
			,{input[[2]]}
			,input[[3]]
			,{output[[1]]}
			,{output[[2]]}
			,output[[3]]
		]
	]
	,{i, Length[testCases]}
];


(* ::Text:: *)
(*Export reference values.*)


SetDirectory[NotebookDirectory[]];
Export["csv/apint_normalize_Ref.csv", N[results, 30]]
