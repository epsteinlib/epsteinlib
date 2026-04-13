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
(*Maximum number of limbs, hard-coded here and in C as well.*)


APINTMAXLIMBS = 32;


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
(*Normalizes hpdyad representation to canonical form. Removes leading zero limbs from high end. *)
(*	If result is zero: returns sign=+1, exp2=0, all limbs=0.*)
(* Otherwise: left-shifts entire mantissa so MSB of highest non-zero limb is 1. Adjusts exp2 to compensate for the shift.*)
(* 	Returns {new_sign, new_exp2, new_limbs}.*)


apintNormalize[{sign_,exp2_,limbs_}]:=Module[
	{n,s,newLimbs,newExp2}
	,n=Length[limbs]
	;While[n>0&&limbs[[n]]==0,n=n-1]
	;If[n==0,Return[{1,0,Table[0,{APINTMAXLIMBS}]}]]
	;s=clz32[limbs[[n]]]
	;If[s==0
		,newLimbs=limbs
		;newExp2=exp2
		,newLimbs=shiftLimbsLeft[limbs,s]
		;newExp2=exp2-s
	]
	;newLimbs=Join[Take[newLimbs,n],Table[0,{APINTMAXLIMBS -n}]]
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


apintMulTest[valA_,signA_,valB_,signB_]:=Module[{prod}
	,prod=signA*valA*signB*valB
	;If[prod==0
	,0.0
	,N[prod]
	]
]


(* ::Text:: *)
(*Convert arbitrary integer to full APINTMAXLIMBS-limb apint representation. Returns {limb0, limb1, ..., limb7, exp2, sign}.*)


toApintFull[x_,s_]:=Module[
	{num,den,denExp,tz,shifted,limbs}
	,If[x==0
		,Return[Join[Table[0,{APINTMAXLIMBS}],{0,1}]]
	]
	;num=Numerator[x]
	;den=Denominator[x]
	;denExp=IntegerExponent[den,2]
	;tz=IntegerExponent[num,2]-denExp
	;shifted=Quotient[num,2^IntegerExponent[num,2]]
	;limbs=Table[Mod[Quotient[shifted,2^(32*i)],2^32],{i,0,APINTMAXLIMBS-1}]
	;Flatten[{limbs,tz,s}]
]


(* ::Text:: *)
(*Simulate hpdyad_mul: multiply two high precision dyadic numbers (hpdyad's) and return normalized result.*)
(*	Input: valA, signA, valB, signB*)
(*	Output: {limb0-7, exp2, sign} after normalization.*)


apintMulTest[valA_,signA_,valB_,signB_]:=Module[
	{prod,prodApint,normalized,limbList}
	,If[valA==0||valB==0
	    ,Return[Join[Table[0,{APINTMAXLIMBS}],{0,1}]]
	]
	;prod=signA*valA*signB*valB
	;prodApint=toApintFull[Abs[prod],Sign[prod]]
	;limbList=Take[prodApint,APINTMAXLIMBS]
	;normalized=apintNormalize[{prodApint[[APINTMAXLIMBS+2]],prodApint[[APINTMAXLIMBS+1]],limbList}]
	;Flatten[{normalized[[3]],normalized[[2]],normalized[[1]]}]
]


(* ::Text:: *)
(*Arbitrary precision integer right-shift with sticky bit.*)


apintShrStickyRef[{sign_, exp2_, limbs_}, bits_] := Module[
	{M, Mshifted, sticky, newLimbs, newExp2}
	,M = FromDigits[Reverse[limbs], 2^32]
	;If[M == 0
		,{1, 0, Table[0, {APINTMAXLIMBS}]}
		,sticky = If[Mod[M, 2^bits] != 0, 1, 0]
		;Mshifted = BitOr[Floor[M / 2^bits], sticky]
		;newExp2 = exp2 + bits
		;newLimbs = Reverse[IntegerDigits[Mshifted, 2^32, APINTMAXLIMBS]]
		;{sign, newExp2, newLimbs}
	]
]


(* ::Text:: *)
(*Convert  integer  to  apint.*)


ratToApint[n_] := Module[
	{raw, sign, absN}
	,absN = Abs[n]
	;sign = If[n == 0, 1, Sign[n]]
	;raw = toApintFull[absN, sign]
	;apintNormalize[{raw[[APINTMAXLIMBS+2]], raw[[APINTMAXLIMBS+1]], raw[[1;;APINTMAXLIMBS]]}]
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


While[Length[matrixTestData2D]<800
	,temp=RandomReal[{-10,10},{2,2}]
	;If[ConditionNumber[temp]<10
		,AppendTo[matrixTestData2D,matrixTestEntry[temp]]
	]
	;Clear[temp]
]


(* ::Text:: *)
(*Generate special 3D matrix reference values.*)


matrixTestData3D={
	matrixTestEntry[{{0,2,3},{1,1,1},{4,5,6}}]
};


(* ::Text:: *)
(*Generate random 3D matrix reference values.*)


While[Length[matrixTestData3D]<100
	,temp=RandomReal[{-10,10},{3,3}]
	;If[ConditionNumber[temp]<10
		,AppendTo[matrixTestData3D,matrixTestEntry[temp]]
	]
	;Clear[temp]
]


(* ::Text:: *)
(*Generate 4D matrix reference values.*)


matrixTestData4D={};
While[Length[matrixTestData4D]<99
	,temp=RandomReal[{-10,10},{4,4}]
	;If[ConditionNumber[temp]<10
		,AppendTo[matrixTestData4D,matrixTestEntry[temp]]
	]
	;Clear[temp]
]


(* ::Text:: *)
(*Combine reference data.*)


matrixTestData=Join[matrixTestData2D,matrixTestData3D,matrixTestData4D];


(* ::Text:: *)
(*Safe reference values.*)


SetDirectory[NotebookDirectory[]];
Export["csv/sdv_Ref.csv", matrixTestData];


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
Export["csv/hpdyad_set_ull_Ref.csv", testCases];


(* ::Subsection::Closed:: *)
(*Test arbitrary precision integer normalization function*)


(* ::Text:: *)
(*Generate test inputs.*)


SeedRandom[42];
testCases = Join[{
	(* 1. Zero *)
	{1, 0, PadRight[{}, APINTMAXLIMBS, 0]}
	,{-1, 5, PadRight[{}, APINTMAXLIMBS, 0]}
	(* 2. Single limb *)
	,{1, 0, PadRight[{1}, APINTMAXLIMBS, 0]}
	,{-1, 10, PadRight[{1}, APINTMAXLIMBS, 0]}
	(* 3. Already normalized *)
	,{1, 0, PadRight[{BitShiftLeft[1, 31]}, APINTMAXLIMBS, 0]}
	,{1, 5, PadRight[{0, BitShiftLeft[1, 31]}, APINTMAXLIMBS, 0]}
	(* 4. Middle zeros *)
	,{1, 0, PadRight[{5, 0, 0, 7}, APINTMAXLIMBS, 0]}
	,{-1, -3, PadRight[{100, 0, 200}, APINTMAXLIMBS, 0]}
	(* 5. Multiple leading zeros *)
	,{1, 0, PadRight[{42}, APINTMAXLIMBS, 0]}
	,{1, 0, PadRight[{1234567, 9876543}, APINTMAXLIMBS, 0]}
	(* 6. Shift amounts *)
	,{1, 0, PadRight[{BitShiftLeft[1, 15]}, APINTMAXLIMBS, 0]}
	,{1, 0, PadRight[{BitShiftLeft[1, 30]}, APINTMAXLIMBS, 0]}
	,{1, 0, PadRight[{0, BitShiftLeft[1, 7]}, APINTMAXLIMBS, 0]}
	(* 7. Full mantissa *)
	,{1, 100, PadRight[{16^^FFFFFFFF, 0, 16^^12345678, 0, 16^^ABCDEF00, 0, 0, 16^^80000000}, APINTMAXLIMBS, 0]}
	(* 8. Dangerous *)
	,{1, 0, PadRight[{16^^FFFF0000, 16^^0000FFFF}, APINTMAXLIMBS, 0]}
	,{1, 0, PadRight[{16^^80000001, 16^^00000001}, APINTMAXLIMBS, 0]}
	(* 9. Large exp2 *)
	,{1, 1000000, PadRight[{123}, APINTMAXLIMBS, 0]}
	,{-1, -500000, PadRight[{16^^7FFFFFFF}, APINTMAXLIMBS, 0]}
	}
	,(* 10. Random values *)
		Table[
			{1, RandomInteger[{-1000, 1000}], RandomInteger[{0, 16^^FFFFFFFF}, APINTMAXLIMBS]}
			,{i,981}
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
Export["csv/hpdyad_normalize_Ref.csv", results];


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
	,{969}
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


(* ::Text:: *)
(*Export reference values.*)


SetDirectory[NotebookDirectory[]];
Export["csv/hpdyad_mul_Ref.csv",apintMulTestTable];


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
	Table[{1,0,Table[0,{APINTMAXLIMBS}],RandomInteger[{0,300}]},{10}]
	(*Noshift(bits=0)*)
	,Table[Module[{n,limbs},n=RandomInteger[{1,APINTMAXLIMBS}];limbs=PadRight[Table[RandomInteger[{1,2^32-1}],{n}],APINTMAXLIMBS,0];{RandomChoice[{-1,1}],RandomInteger[{-100,100}],limbs,0}],{10}]
	(*Shiftwithinsinglelimb*)
	,Table[{RandomChoice[{-1,1}],RandomInteger[{-100,100}],PadRight[{RandomInteger[{1,2^32-1}]},APINTMAXLIMBS,0],RandomInteger[{1,31}]},{10}]
	(*Shiftexactlyonelimbboundary*)
	,Table[Module[{n,limbs},n=RandomInteger[{2,APINTMAXLIMBS}];limbs=PadRight[Table[RandomInteger[{0,2^32-1}],{n}],APINTMAXLIMBS,0];{RandomChoice[{-1,1}],RandomInteger[{-100,100}],limbs,32}],{10}]
	(*Shiftacrosslimbboundary*)
	,Table[Module[{n,limbs},n=RandomInteger[{2,APINTMAXLIMBS}];limbs=PadRight[Table[RandomInteger[{0,2^32-1}],{n}],APINTMAXLIMBS,0];{RandomChoice[{-1,1}],RandomInteger[{-100,100}],limbs,RandomInteger[{33,63}]}],{10}]
	(*Completeshift-out,nosticky*)
	,Table[{RandomChoice[{-1,1}],RandomInteger[{-100,100}],PadRight[{2^RandomInteger[{0,31}]},APINTMAXLIMBS,0],RandomInteger[{32,300}]},{10}]
	(*Completeshift-out,withsticky*)
	,Table[{RandomChoice[{-1,1}],RandomInteger[{-100,100}],PadRight[{RandomInteger[{1,2^32-1}]},APINTMAXLIMBS,0],RandomInteger[{32,300}]},{10}]
	(*Exactly256bitshift*)
	,Table[{RandomChoice[{-1,1}],RandomInteger[{-100,100}],Table[RandomInteger[{1,2^32-1}],{APINTMAXLIMBS}],256},{10}]
	(*Maxlimbs,smallshift*)
	,Table[{RandomChoice[{-1,1}],RandomInteger[{-100,100}],Table[RandomInteger[{1,2^32-1}],{APINTMAXLIMBS}],RandomInteger[{1,31}]},{10}]
	(*All0xFFFFFFFFlimbs*)
	,Table[Module[{n,limbs},n=RandomInteger[{1,APINTMAXLIMBS}];limbs=PadRight[Table[2^32-1,{n}],APINTMAXLIMBS,0];{RandomChoice[{-1,1}],RandomInteger[{-100,100}],limbs,RandomInteger[{1,100}]}],{10}]
];


(* ::Text:: *)
(*Random cases.*)


Do[
	sign0[i] = RandomChoice[{-1, 1}]
	;n0[i] = RandomInteger[{0,APINTMAXLIMBS}]
	;exp0[i] = RandomInteger[{-1000, 1000}]
	;limbs0[i] = If[n0[i] == 0
		,Table[0, {APINTMAXLIMBS}]
		,Join[
			Table[RandomInteger[{0, 2^32 - 1}], {n0[i] - 1}]
			,{RandomInteger[{1, 2^32 - 1}]}
			,Table[0, {APINTMAXLIMBS - n0[i]}]
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
Export["csv/hpdyad_shr_sticky_Ref.csv", shrStickyTab];


(* ::Subsection::Closed:: *)
(*Test arbitrary precision integer addition function*)


(* ::Text:: *)
(*Setup reference values.*)


edgeCases=Join[
	(* Zero operands *)
	{{0,0},{0,12345678},{87654321,0},{0,-999999},{-111111,0}}
	(* Same sign, simple *)
	,{{100,200},{-500,-300},{1000000,2000000},{-9999999,-1111111},{2^32,2^32},{-2^40,-2^40},{2^63,2^63},{-2^100,-2^100},{123456789012345,987654321098765},{-111111111111111,-222222222222222}}
	(* Same sign with carries *)
	,{{2^32-1,1},{2^32-1,2^32-1},{2^64-1,1},{2^64-1,2^64-1},{2^96-1,2^96-1},{2^128-1,2^128-1},{-(2^32-1),-1},{-(2^64-1),-(2^64-1)},{2^32-100,200},{2^64-1000,2000}}
	(* Different signs, large d *)
	,{{2^500,-2^50},{-2^500,2^50},{2^950,-1},{-2^950,1},{2^900,-2^100},{2^800,-2^20},{-2^800,2^20},{2^500,-(2^500-1)},{2^700,-2^80},{-2^700,2^80},{2^900,-2^130},{2^940,-2^10},{-2^940,2^10},{2^500,-2^400},{2^400,-2^500}}
	(* Different signs, |a|<|b|, large d *)
	,{{2^50,-2^500},{-2^50,2^500},{1,-2^950},{-1,2^950},{2^100,-2^900},{2^20,-2^800},{-2^20,2^800},{2^80,-2^700},{-2^80,2^700},{2^130,-2^900}}
	(* Equal magnitude opposite sign *)
	,{{100,-100},{-100,100},{2^32,-2^32},{-2^32,2^32},{2^64,-2^64},{2^100,-2^100},{-2^128,2^128},{123456789,-123456789},{2^500,-2^500},{2^950,-2^950}}
	(* Borrow propagation *)
	,{{2^500,-1},{2^700,-1},{2^900,-1},{2^950,-1},{2^500,-(2^500-1)},{2^700,-(2^700-1)},{2^900,-(2^900-1)},{2^500*2,-3},{2^700*2,-3},{2^500+2^400,-2^400},{2^700+2^500,-2^500},{2^900+2^700,-2^700},{-2^500,1},{-2^700,1},{-2^900,1}}
	(* Large exp2 differences *)
	,{{2^950,1},{2^950,-1},{2^950,2^500},{2^950,-2^500},{-2^950,2^500},{2^940,2^50},{2^940,-2^50},{1,2^950},{-1,2^950},{2^500,2^950}}
	(* Powers of 2 *)
	,{{2^31,2^31},{2^500,2^500},{2^700,2^700},{2^900,2^900},{2^31,-2^31},{2^500,-2^500},{2^700,-2^700},{2^900,-2^900},{2^950,2^950},{2^950,-2^950}}
];


(* ::Text:: *)
(*Setup random reference values.*)


SeedRandom[42];
nTests=904;
randomTests=Table[
	Module[
		{a,b,range}
		,{
			range=RandomChoice[{2^500,2^700,2^900,2^940,2^950}]
			,a=RandomInteger[{-range,range}]
			,b=RandomInteger[{-range,range}]
			,{a,b}
		}
	]
	,{nTests}
];

allTests=Join[edgeCases,randomTests];


(* ::Text:: *)
(*Generate reference values.*)


testData=ParallelTable[
	Module[{a,b,aApint,bApint,refApint,ref}
		,{a,b}={allTests[[i,1]],allTests[[i,2]]}
		;ref=a+b
		;aApint=ratToApint[a]
		;bApint=ratToApint[b]
		;refApint=ratToApint[ref]
		;Flatten[{aApint[[1]],aApint[[2]],aApint[[3]],bApint[[1]],bApint[[2]],bApint[[3]],refApint[[1]],refApint[[2]],refApint[[3]]}]
	]
	,{i,Length[allTests]}
];


(* ::Text:: *)
(*Export test values.*)


SetDirectory[NotebookDirectory[]];
Export["csv/hpdyad_add_Ref.csv", testData];


(* ::Subsection::Closed:: *)
(*Test arbitrary precision rational addition function*)


(* ::Text:: *)
(*Setup edge case reference values.*)


fractionalEdgeCases = Join[
	(* Zero operands *)
	{{0,0},{0,1/2},{1/4,0},{0,-3/8},{-7/16,0}}
	(* Same sign fractions *)
	,{{1/2,1/4},{-1/2,-1/4},{3/4,1/8},{-7/8,-1/16},{1/2,1/2}}
	(* Mixed integer and fraction *)
	,{{1,1/2},{-1,-1/2},{100,1/4},{-100,-3/4},{2^32,1/2}}
	(* Large integer + small fraction *)
	,{{2^100,1/2},{2^500,1/4},{2^900,1/8},{-2^100,-1/2},{-2^500,-1/4}}
	(* Small fraction + large integer *)
	,{{1/2,2^100},{1/4,2^500},{1/8,2^900},{-1/2,-2^100},{-1/4,-2^500}}
	(* Opposite signs with fractions *)
	,{{1/2,-1/4},{-1/2,1/4},{3/4,-1/2},{-3/4,1/2},{1,-1/2}}
	(* Cancellation to fraction *)
	,{{1,-1/2},{2,-3/2},{3/2,-1},{-3/2,1},{5/4,-1/4}}
	(* Cancellation to zero *)
	,{{1/2,-1/2},{-1/4,1/4},{3/8,-3/8},{-7/16,7/16},{2^(-10),-2^(-10)}}
	(* Very small fractions *)
	,{{2^(-50),2^(-50)},{2^(-100),2^(-100)},{-2^(-50),-2^(-50)},{2^(-50),-2^(-51)}}
	(* Large exp2 difference with fractions *)
	,{{2^200,2^50},{2^50,2^200},{2^300,2^150},{2^150,2^300}}
	(* Negative exp2 both operands *)
	,{{2^(-10),2^(-20)},{2^(-50),2^(-100)},{-2^(-10),-2^(-20)},{2^(-10),-2^(-20)}}
	(* Borrow with fractions *)
	,{{1,-2^(-50)},{2^100,-2^(-50)},{2^500,-2^(-50)},{1/2,-2^(-60)}}
	(* Carry with fractions *)
	,{{1/2,1/2},{3/4,1/4},{7/8,1/8},{2^(-10)-2^(-50),2^(-50)}}
  ];


(* ::Text:: *)
(*Setup random reference values.*)


SeedRandom[44];
nTestsFrac=939;
randomFractionalTests=Table[
	Module[
		{a,b,expA,expB,mantA,mantB}
		,expA=RandomInteger[{-50,900}]
		;expB=RandomInteger[{expA-200,expA+200}]
		;mantA=RandomInteger[{1,2^64}]*RandomChoice[{1,-1}]
		;mantB=RandomInteger[{1,2^64}]*RandomChoice[{1,-1}]
		;a=mantA*2^expA
		;b=mantB*2^expB
		;{a,b}
	]
	,{nTestsFrac}
];

allFractionalTests=Join[fractionalEdgeCases,randomFractionalTests];


(* ::Text:: *)
(*Generate reference values.*)


fractionalTestData = ParallelTable[
	Module[{a,b,aApint,bApint,refApint,ref}
		,{a,b} = allFractionalTests[[i]]
		;ref = a + b
		;aApint = ratToApint[a]
		;bApint = ratToApint[b]
		;refApint = ratToApint[ref]
		;Flatten[{aApint[[1]],aApint[[2]],aApint[[3]],bApint[[1]],bApint[[2]],bApint[[3]],refApint[[1]],refApint[[2]],refApint[[3]]}]
	]
	,{i,Length[allFractionalTests]}
];


(* ::Text:: *)
(*Export reference values.*)


SetDirectory[NotebookDirectory[]];
Export["csv/hpdyad_add_fractional_Ref.csv", fractionalTestData];


(* ::Subsection::Closed:: *)
(*Test apint to double conversion with round-to-nearest-even*)


(* ::Text:: *)
(*Setup edge case reference values.*)


edgeCases=Join[
	(* Zero *)
	{0}
	(* Small exact integers *)
	,{1,-1,2,-2,100,-100,12345,-12345,2^31-1,-(2^31-1),2^32-1,2^32,2^32+1}
	(* Exact powers of 2 *)
	,{2^52,2^53-1,2^53,2^53+1,2^53+2,2^54,2^63,2^64,-2^52,-2^53,-(2^53+1),-2^64}
	(* 53-bit boundary: exactly representable *)
	,{2^53-1,-(2^53-1),2^52+2^51,-(2^52+2^51)}
	(* 54-bit: rounding needed *)
	,{2^54+1,2^54+2,2^54+3,2^54+4,-(2^54+1),-(2^54+2),-(2^54+3),-(2^54+4)}
	(* Round to even cases: 2^54 + k where k determines rounding *)
	,{2^54+0,2^54+1,2^54+2,2^54+3,2^54+4,2^54+5,2^54+6,2^54+7,2^54+8}
	(* Large values within range *)
	,{2^100,2^100+2^47,2^100+2^48,-(2^100),2^200,2^500,2^700,2^900,2^950,-2^950}
	(* Near max double exponent *)
	,{2^1023,2^1023+2^970,-(2^1023),-(2^1023+2^970)}
	(* Odd/even rounding at various scales *)
	,{3*2^52,3*2^52+1,3*2^52+2,5*2^60,5*2^60+2^8,7*2^100,7*2^100+2^48}
	(* Multi-limb with rounding *)
	,{2^64+2^11,2^64+2^11+1,2^64+2^12,2^96+2^43,2^96+2^43+1,2^96+2^44}
	(* Sticky bit cases: many low bits set *)
	,{2^64+2^11+2^5+1,2^100+2^47+2^30+2^10+1,2^200+2^147+2^100+2^50+2^10+1}
];


(* ::Text:: *)
(*Setup  random  reference  values .*)


SeedRandom[43];
nTests=922;
randomTests=Table[
Module[
	{bits,x,sign}
	,bits=RandomInteger[{1,1000}]
	;x=RandomInteger[{2^(bits-1),2^bits-1}]
	;sign=RandomChoice[{1,-1}]
	;sign*x
	]
	,{nTests}
];
allTests=Join[edgeCases,randomTests];


(* ::Text:: *)
(*Generate reference values.*)


testData=ParallelTable[
	Module[
		{x,xApint,refDouble}
		,x=allTests[[i]]
		;xApint=ratToApint[x]
		;refDouble=N[x]
		;Flatten[{xApint[[1]],xApint[[2]],xApint[[3]],refDouble}]
	]
	,{i,Length[allTests]}
];


(* ::Text:: *)
(*Export  test  values.*)


SetDirectory[NotebookDirectory[]];
Export["csv/hpdyad_to_double_Ref.csv",testData];
