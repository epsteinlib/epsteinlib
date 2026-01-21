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


apintNormalize[{sign_,exp2_,limbs_}]:=Module[
	{n,s,newLimbs,newExp2}
	,n=Length[limbs]
	;While[n>0&&limbs[[n]]==0,n=n-1]
	;If[n==0,Return[{1,0,Table[0,{8}]}]]
	;s=clz32[limbs[[n]]]
	;If[s==0
		,newLimbs=limbs
		;newExp2=exp2
		,newLimbs=shiftLimbsLeft[limbs,s]
		;newExp2=exp2-s
	]
	;newLimbs=Join[Take[newLimbs,n],Table[0,{8-n}]]
	(* Trim low-end zeros *)
	;While[n>0&&newLimbs[[1]]==0
		,newLimbs=Join[Drop[newLimbs,1],{0}]
		;n=n-1
		;newExp2=newExp2+32
	]
	;{sign,newExp2,newLimbs}
]


(* ::Text:: *)
(*Define multiplication test function. Takes two unsigned long long values and their signs, returns product as double. Simulates: apint_set_ull \[RightArrow] apint_mul \[RightArrow] apint_to_double_trunc.*)


apintMulTest[valA_,signA_,valB_,signB_]:=Module[{prod}
	,prod=signA*valA*signB*valB
	;If[prod==0
	,0.0
	,N[prod]
	]
]


(* ::Text:: *)
(*Convert arbitrary integer to full 8-limb apint representation. Returns {limb0, limb1, ..., limb7, exp2, sign}.*)


toApintFull[x_,s_]:=Module[{tz,shifted,limbs}
	,If[x==0
		,Return[{0,0,0,0,0,0,0,0,0,1}]
	]
	;tz=IntegerExponent[x,2]
	;shifted=Quotient[x,2^tz]
	;limbs=Table[Mod[Quotient[shifted,2^(32*i)],2^32],{i,0,7}]
	;Flatten[{limbs,tz,s}]
]


(* ::Text:: *)
(*Simulate apint_mul: multiply two apints and return normalized result.*)
(*	Input: valA, signA, valB, signB*)
(*	Output: {limb0-7, exp2, sign} after normalization.*)


apintMulTest[valA_,signA_,valB_,signB_]:=Module[
	{prod,prodApint,normalized,limbList}
	,If[valA==0||valB==0
		,Return[{0,0,0,0,0,0,0,0,0,1}]
	]
	;prod=signA*valA*signB*valB
	;prodApint=toApintFull[Abs[prod],Sign[prod]]
	;limbList=Take[prodApint,8]
	;normalized=apintNormalize[{prodApint[[10]],prodApint[[9]],limbList}]
	;Flatten[{normalized[[3]],normalized[[2]],normalized[[1]]}]
]


(* ::Text:: *)
(*Arbitrary precision integer right-shift with sticky bit.*)


apintShrStickyRef[{sign_, exp2_, limbs_}, bits_] := Module[
	{M, Mshifted, sticky, newLimbs, newExp2}
	,M = FromDigits[Reverse[limbs], 2^32]
	;If[M == 0
		,{1, 0, Table[0, {8}]}
		,sticky = If[Mod[M, 2^bits] != 0, 1, 0]
		;Mshifted = BitOr[Floor[M / 2^bits], sticky]
		;newExp2 = exp2 + bits
		;newLimbs = Reverse[IntegerDigits[Mshifted, 2^32, 8]]
		;{sign, newExp2, newLimbs}
	]
]


(* ::Text:: *)
(*Convert  integer  to  apint.*)


intToApint[n_] := Module[{raw, sign, absN}
,absN = Abs[n]
;sign = If[n == 0, 1, Sign[n]]
;raw = toApintFull[absN, sign]
;apintNormalize[{raw[[10]], raw[[9]], raw[[1;;8]]}]
]


(* ::Subsection::Closed:: *)
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


(* ::Subsection::Closed:: *)
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
testCases = Join[{
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
	},(* 10. Random values *)
	Table[
		{1, RandomInteger[{-1000, 1000}], RandomInteger[{0, 16^^FFFFFFFF}, 8]}
		,{i,9810}
	]
];


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
Export["csv/apint_normalize_Ref.csv", results];


(* ::Subsection::Closed:: *)
(*Test exact multiplication*)


(* ::Text:: *)
(*Special cases for reference.*)


specialCases={
	{0,1,0,1}
	,{0,1,12345,1}
	,{12345,-1,0,1}
	,{1,1,1,1}
	,{1,-1,1,-1}
	,{1,1,999999999999,1}
	,{4294967295,1,1,1}
	,{4294967295,1,4294967295,1}
	,{4294967296,1,4294967296,1}
	,{18446744073709551615,1,1,1}
	,{18446744073709551615,1,2,1}
	,{18446744073709551615,-1,2,-1}
	,{9223372036854775807,1,2,1}
	,{4294967295,1,4294967296,1}
	,{281474976710656,1,65536,1}
	,{65536,1,65536,1}
	,{4096,1,4096,1}
	,{16777216,1,16777216,1}
	,{1099511627776,1,16777216,1}
	,{281474976710655,1,65537,1}
	,{549755813888,1,549755813888,1}
	,{4294967295,-1,4294967295,-1}
	,{1,-1,1,1}
	,{999999999999,-1,999999999999,-1}
	,{72057594037927936,1,256,1}
	,{562949953421312,1,32768,1}
	,{18014398509481984,1,1024,1}
	,{4503599627370496,1,4096,1}
	,{1125899906842624,1,16384,1}
	,{17179869184,1,1073741824,1}
};


(* ::Text:: *)
(*Random reference values.*)


SeedRandom[42];

randomCases=Table[
	{RandomInteger[{0,18446744073709551615}]
	,RandomChoice[{1,-1}]
	,RandomInteger[{0,18446744073709551615}]
	,RandomChoice[{1,-1}]
	}
	,{9696}
];

allCases=Join[specialCases,randomCases];


(* ::Text:: *)
(*Generate reference outputs,*)


apintMulTestTable=ParallelTable[
	Flatten[Join[
		allCases[[i]]
		,apintMulTest@@allCases[[i]]
	]]
	,{i,1,Length[allCases]}
];


apintMulTestTable[[1;;10]]


(* ::Text:: *)
(*Export reference values.*)


SetDirectory[NotebookDirectory[]];
Export["csv/apint_mul_Ref.csv",apintMulTestTable];


(* ::Subsection::Closed:: *)
(*Test right shift with sticky bits function*)


(* ::Text:: *)
(*Setup reference values.*)


SeedRandom[42];
nTests = 899;
nEdge = 100;


(* ::Text:: *)
(*Edge cases*)


edgeCases=Join[
	(*Zeroinput*)
	Table[{1,0,Table[0,{8}],RandomInteger[{0,300}]},{10}]
	(*Noshift(bits=0)*)
	,Table[Module[{n,limbs},n=RandomInteger[{1,8}];limbs=PadRight[Table[RandomInteger[{1,2^32-1}],{n}],8,0];{RandomChoice[{-1,1}],RandomInteger[{-100,100}],limbs,0}],{10}]
	(*Shiftwithinsinglelimb*)
	,Table[{RandomChoice[{-1,1}],RandomInteger[{-100,100}],PadRight[{RandomInteger[{1,2^32-1}]},8,0],RandomInteger[{1,31}]},{10}]
	(*Shiftexactlyonelimbboundary*)
	,Table[Module[{n,limbs},n=RandomInteger[{2,8}];limbs=PadRight[Table[RandomInteger[{0,2^32-1}],{n}],8,0];{RandomChoice[{-1,1}],RandomInteger[{-100,100}],limbs,32}],{10}]
	(*Shiftacrosslimbboundary*)
	,Table[Module[{n,limbs},n=RandomInteger[{2,8}];limbs=PadRight[Table[RandomInteger[{0,2^32-1}],{n}],8,0];{RandomChoice[{-1,1}],RandomInteger[{-100,100}],limbs,RandomInteger[{33,63}]}],{10}]
	(*Completeshift-out,nosticky*)
	,Table[{RandomChoice[{-1,1}],RandomInteger[{-100,100}],PadRight[{2^RandomInteger[{0,31}]},8,0],RandomInteger[{32,300}]},{10}]
	(*Completeshift-out,withsticky*)
	,Table[{RandomChoice[{-1,1}],RandomInteger[{-100,100}],PadRight[{RandomInteger[{1,2^32-1}]},8,0],RandomInteger[{32,300}]},{10}]
	(*Exactly256bitshift*)
	,Table[{RandomChoice[{-1,1}],RandomInteger[{-100,100}],Table[RandomInteger[{1,2^32-1}],{8}],256},{10}]
	(*Maxlimbs,smallshift*)
	,Table[{RandomChoice[{-1,1}],RandomInteger[{-100,100}],Table[RandomInteger[{1,2^32-1}],{8}],RandomInteger[{1,31}]},{10}]
	(*All0xFFFFFFFFlimbs*)
	,Table[Module[{n,limbs},n=RandomInteger[{1,8}];limbs=PadRight[Table[2^32-1,{n}],8,0];{RandomChoice[{-1,1}],RandomInteger[{-100,100}],limbs,RandomInteger[{1,100}]}],{10}]
];


(* ::Text:: *)
(*Random cases.*)


Do[
	sign0[i] = RandomChoice[{-1, 1}]
	;n0[i] = RandomInteger[{0, 8}]
	;exp0[i] = RandomInteger[{-1000, 1000}]
	;limbs0[i] = If[n0[i] == 0
		,Table[0, {8}]
		,Join[
			Table[RandomInteger[{0, 2^32 - 1}], {n0[i] - 1}]
			,{RandomInteger[{1, 2^32 - 1}]}
			,Table[0, {8 - n0[i]}]
		]
	]
;bits0[i] = RandomInteger[{0, 300}]
,{i, nTests}
]


(* ::Text:: *)
(*Generate reference values.*)


shrStickyTab=ParallelTable[
	Module[
		{apintIn,apintOut,sign,exp2,limbs,bits}
		,{sign,exp2,limbs,bits}=If[i<=Length[edgeCases]
			,edgeCases[[i]]
			,{sign0[i-Length[edgeCases]],exp0[i-Length[edgeCases]]
			,limbs0[i-Length[edgeCases]],bits0[i-Length[edgeCases]]}
		]
		;apintIn={sign,exp2,limbs}
		;apintOut=apintShrStickyRef[apintIn,bits]
		;Flatten[Join[
		{sign,exp2}
		,limbs
		,{bits}
		,{apintOut[[1]],apintOut[[2]]}
		,apintOut[[3]]
		]]
	]
	,{i,Length[edgeCases]+nTests}
];


(* ::Text:: *)
(*Export reference values.*)


SetDirectory[NotebookDirectory[]];
Export["csv/apint_shr_sticky_Ref.csv", shrStickyTab];


(* ::Subsection:: *)
(*Test arbitrary precision integer addition function*)


nTests = 904;
(* Zero operands *)
zeroEdgeCases = {
{0, 0}
,{0, 12345678}
,{87654321, 0}
,{0, -999999}
,{-111111, 0}
};

(* ============================================ *)
(* CATEGORY 2: Same sign, simple (10 cases)    *)
(* ============================================ *)
sameSignSimple = {
{100, 200}
,{-500, -300}
,{1000000, 2000000}
,{-9999999, -1111111}
,{2^32, 2^32}
,{-2^40, -2^40}
,{2^63, 2^63}
,{-2^100, -2^100}
,{123456789012345, 987654321098765}
,{-111111111111111, -222222222222222}
};

(* ============================================ *)
(* CATEGORY 3: Same sign with carries (10 cases) *)
(* ============================================ *)
sameSignCarry = {
{2^32-1, 1}
,{2^32-1, 2^32-1}
,{2^64-1, 1}
,{2^64-1, 2^64-1}
,{2^96-1, 2^96-1}
,{2^128-1, 2^128-1}
,{-(2^32-1), -1}
,{-(2^64-1), -(2^64-1)}
,{2^32-100, 200}
,{2^64-1000, 2000}
};

(* ============================================ *)
(* CATEGORY 4: Different signs, small d (15 cases) *)
(* ============================================ *)
diffSignSmallD = {
{1000, -300}
,{-1000, 300}
,{2^40, -2^38}
,{-2^40, 2^38}
,{2^100, -2^98}
,{2^64, -1}
,{2^64, -(2^32-1)}
,{2^96, -(2^96-2)}
,{1234567890, -123456789}
,{-9876543210, 987654321}
,{2^128, -2^126}
,{2^200, -2^198}
,{2^32+1, -1}
,{2^64+1, -1}
,{-2^100-1, 1}
};

(* ============================================ *)
(* CATEGORY 5: Different signs, |a| < |b|, small d (10 cases) *)
(* ============================================ *)
diffSignBLargerSmallD = {
{300, -1000}
,{-300, 1000}
,{2^38, -2^40}
,{-2^38, 2^40}
,{1, -2^10}
,{2^32-1, -2^34}
,{2^98, -2^100}
,{123456789, -1234567890}
,{-987654321, 9876543210}
,{2^126, -2^128}
};

(* ============================================ *)
(* CATEGORY 6: Equal magnitude opposite sign (10 cases) *)
(* ============================================ *)
equalMagnitude = {
{100, -100}
,{-100, 100}
,{2^32, -2^32}
,{-2^32, 2^32}
,{2^64, -2^64}
,{2^100, -2^100}
,{-2^128, 2^128}
,{123456789, -123456789}
,{2^200, -2^200}
,{999999999999, -999999999999}
};

(* ============================================ *)
(* CATEGORY 7: Borrow propagation (15 cases)   *)
(* ============================================ *)
borrowPropagation = {
{2^32, -1}
,{2^33, -1}
,{2^34, -1}
,{2^35, -1}
,{2^32, -(2^32-1)}
,{2^33, -(2^33-1)}
,{2^34, -(2^34-1)}
,{2^32*2, -3}
,{2^33*2, -3}
,{2^32+2^16, -2^16}
,{2^33+2^17, -2^17}
,{2^34+2^18, -2^18}
,{-2^32, 1}
,{-2^33, 1}
,{-2^34, 1}
};

(* ============================================ *)
(* CATEGORY 8: Small exp2 differences only (10 cases) *)
(* ============================================ *)
smallExp2Diff = {
{2^100, 1}
,{2^100, -1}
,{2^100, 2^99}
,{2^100, -2^99}
,{-2^100, 2^99}
,{2^128, 2^127}
,{2^128, -2^127}
,{1, 2^20}
,{-1, 2^20}
,{2^100, 2^110}
};

(* ============================================ *)
(* CATEGORY 9: Powers of 2 (edge cases) (10 cases) *)
(* ============================================ *)
powersOfTwo = {
{2^31, 2^31}
,{2^32, 2^32}
,{2^63, 2^63}
,{2^64, 2^64}
,{2^31, -2^31}
,{2^32, -2^32}
,{2^64, -2^64}
,{2^128, 2^128}
,{2^128, -2^128}
,{2^200, 2^200}
};

(* ============================================ *)
(* CATEGORY 10: Random with small d (9910 cases) *)
(* ============================================ *)
SeedRandom[42];
randomTests = Table[
	Module[
		{a, b, aApint, bApint, d, range, aExp2, bRange}
		,range = RandomChoice[{2^50, 2^100, 2^150, 2^200, 2^250}]
		;a = RandomInteger[{-range, range}]
		;aApint = intToApint[a]
		;aExp2 = aApint[[2]]
		;bRange = 2^(aExp2 + 31)
		;b = RandomInteger[{-bRange, bRange}]
		;{a, b}
	]
,{nTests}
];

(* ============================================ *)
(* COMBINE ALL TEST CASES                      *)
(* ============================================ *)
allEdgeCases = Join[
zeroEdgeCases
,sameSignSimple
,sameSignCarry
,diffSignSmallD
,diffSignBLargerSmallD
,equalMagnitude
,borrowPropagation
,smallExp2Diff
,powersOfTwo
];

allTests = Join[allEdgeCases, randomTests];


(* ::Text:: *)
(*Generate reference values.*)


testData = ParallelTable[
 Module[{a,b,aApint, bApint, refApint, ref}
,{a,b}={allTests[[i, 1]], allTests[[i, 2]]}
;ref = a + b
;aApint = intToApint[a]
;bApint = intToApint[b]
;refApint = intToApint[ref]
;Flatten[{aApint[[1]], aApint[[2]], aApint[[3]], bApint[[1]], bApint[[2]], bApint[[3]], refApint[[1]], refApint[[2]], refApint[[3]]}]
]
,{i, 1,Length[allTests]}
];


(* ::Text:: *)
(*Export test values.*)


SetDirectory[NotebookDirectory[]];
Export["csv/apint_add_Ref.csv", testData];
