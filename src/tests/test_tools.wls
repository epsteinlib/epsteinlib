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
(*Most  Significant  Bit  position  in  32-bit  integer*)


msb32Ref[x_] := If[x == 0, 0, Floor[Log2[x]]]


(* ::Subsection:: *)
(*Test Count Trailing Zeros (64-bit) and Most Significant Bit (32-bit) functions*)


(* ::Text:: *)
(*Setup reference values.*)


SeedRandom[42];
range = Range[995];

Do[
  k0[i] = RandomInteger[{0, 63}];
  maxM = 2^(64 - k0[i]) - 1;
  m0[i] = 2*RandomInteger[{0, Floor[maxM/2]}] + 1;
  x0[i] = m0[i]*2^k0[i];
  msbM0[i] = If[m0[i] < 2^32, msb32Ref[m0[i]],32 + msb32Ref[BitShiftRight[m0[i], 32]] ];
  , {i, range}
]

specialCases = {
  {0, 64, 0},
  {2^63, 63, 0},
  {1, 0, 0},
  {2^32 - 1, 0, 31}
};


(* ::Text:: *)
(*Generate  reference  values .*)


ctzmsbTab = ParallelTable[
  {x0[i], ctz64Ref[x0[i]], msbM0[i]},
  {i, range}
];
ctzmsbTab = Join[specialCases, ctzmsbTab];


(* ::Text:: *)
(*Export reference values.*)


SetDirectory[NotebookDirectory[]];
Export["csv/ctzmsb_Ref.csv", ctzmsbTab];
