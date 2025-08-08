(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
   SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


BeginPackage["EpsteinZeta`"];


EpsteinZeta::usage="EpsteinZeta[\[Nu],A,x,y] computes the Epstein zeta function sum_{z in Lambda} exp(- 2 Pi I y.z)/|z - x|^\[Nu]
using the algorithm in Crandall, R., Unified algorithms for polylogarithm, L-series, and zeta variants. Algorithmic Reflections: Selected Works. PSIpress (2012).";
EpsteinZetaReg::usage="EpsteinZetaReg[\[Nu],A,x,y] computes the regularized Epstein zeta function exp(2 Pi I x.y) sum_{z in Lambda} exp(- 2 Pi I y.(z - x))/|z - x|^\[Nu]  - s\:0302(y)/|det(\[CapitalLambda])|
as in Andreas A. Buchheit et al., Exact continuum representation of long-range interacting systems and emerging exotic phases in unconventional superconductors. Phys. Rev. Res. 5 (4 Oct. 2023), p. 043065.  using a modification of the algorithm in Crandall, R., Unified algorithms for polylogarithm, L-series, and zeta variants. Algorithmic Reflections: Selected Works. PSIpress (2012).";
SetZetaDer::usage="SetZetaDer[\[Nu],A,x,y,\[Alpha]] computes the partiall derivatives of the set zeta function sum_{z in (Lambda - x)} exp(- 2 Pi I y.z)/|z|^\[Nu] with respect to y and some multi-index \[Alpha].";
EpsteinZetaRegDer::usage="EpsteinZetaRegDer[\[Nu],A,x,y] computes the partiall derivatives of the regularized Epstein zeta function exp(2 Pi I x.y) sum_{z in Lambda} exp(- 2 Pi I y.(z - x))/|z - x|^\[Nu]  - s\:0302(y)/|det(\[CapitalLambda])| with respect to y and some multi-index \[Alpha].";
GBessel::usage"GBessel[\[Nu],k,r] computes the incomplete Bessel function 2 int_0^1 t^(-nu-1) exp(-pi k^2 / t^2) exp(-pi r^2 t^2) dt";


Begin["Private`"];


(* Find the library *)
currentDir = DirectoryName[$InputFileName];
parentDir = DirectoryName[currentDir];

subdirsLimited = FileNameJoin[{parentDir, #}] & /@ {"result/lib", "build/src"} // Select[DirectoryQ];

$LibraryPath = Join[{currentDir}, {parentDir}, subdirsLimited, $LibraryPath];

libBaseName = "libepstein";
libPath = FindLibrary[libBaseName];

If[libPath === $Failed,
  Print[StringForm[
    "Library '``' not found. To resolve this:\n" <>
    "1. Ensure the library file (``.so for Linux, ``.dylib for macOS or ``.dll for Windows) exists.\n" <>
    "2. Place the library file in one of these locations:\n" <>
    "   a) The directory of EpsteinZeta.wl:\n``\n" <>
    "   b) The project root directory of EpsteinZeta.wl:\n``\n" <>
    "   c) Any of the following subdirectories:\n``\n" <>
    "   d) Or in any directory listed in the default $LibraryPath:\n``\n" <>
    "3. If you've just added the file, you may need to restart the Mathematica kernel.\n" <>
    "4. If issues persist, check file permissions and ensure the file is readable.",
    libBaseName, libBaseName, libBaseName, libBaseName,
    Column[{currentDir}],
    Column[{parentDir}],
    Column[subdirsLimited],
    Column[Complement[$LibraryPath, {currentDir}, subdirsLimited]]
  ]]
]


(* Interface to C functions *)
foreignFunctionEpsteinZeta = ForeignFunctionLoad[libPath, "epstein_zeta_mathematica_call", {"RawPointer"::["CDouble"], "CDouble", "CInt", "RawPointer"::["CDouble"], "RawPointer"::["CDouble"], "RawPointer"::["CDouble"]} -> "CInt"]
If[Head[foreignFunctionEpsteinZeta] =!= ForeignFunction,
  Print["ForeignFunctionLoad for epstein_zeta_mathematica_call failed."]
]

foreignFunctionEpsteinZetaReg = ForeignFunctionLoad[libPath, "epstein_zeta_reg_mathematica_call", {"RawPointer"::["CDouble"], "CDouble", "CInt", "RawPointer"::["CDouble"], "RawPointer"::["CDouble"], "RawPointer"::["CDouble"]} -> "CInt"]
If[Head[foreignFunctionEpsteinZetaReg] =!= ForeignFunction,
  Print["ForeignFunctionLoad for epstein_zeta_reg_mathematica_call failed."]
]

foreignFunctionSetZetaDer = ForeignFunctionLoad[libPath, "set_zeta_der_mathematica_call", {"RawPointer"::["CDouble"], "CDouble", "CInt", "RawPointer"::["CDouble"], "RawPointer"::["CDouble"], "RawPointer"::["CDouble"], "RawPointer"::["CUnsignedInt"]} -> "CInt"]
If[Head[foreignFunctionEpsteinZetaReg] =!= ForeignFunction,
  Print["ForeignFunctionLoad for set_zeta_der_mathematica_call failed."]
]

foreignFunctionEpsteinZetaRegDer = ForeignFunctionLoad[libPath, "epstein_zeta_reg_der_mathematica_call", {"RawPointer"::["CDouble"], "CDouble", "CInt", "RawPointer"::["CDouble"], "RawPointer"::["CDouble"], "RawPointer"::["CDouble"], "RawPointer"::["CUnsignedInt"]} -> "CInt"]
If[Head[foreignFunctionEpsteinZetaRegDer] =!= ForeignFunction,
  Print["ForeignFunctionLoad for epstein_zeta_reg_der_mathematica_call failed."]
]

foreignFunctionGBessel = ForeignFunctionLoad[libPath, "incomplete_bessel_g_mathematica_call", {"RawPointer"::["CDouble"], "CDouble", "CInt", "RawPointer"::["CDouble"], "RawPointer"::["CDouble"]} -> "CInt"]
If[Head[foreignFunctionGBessel] =!= ForeignFunction,
  Print["ForeignFunctionLoad for incomplete_bessel_g_mathematica_call failed."]
]

(* For checking nan output of C function *)
NaNQ = ResourceFunction["NaNQ"];


(* Internal routines for C function access *)
epsteinZetaInternal[\[Nu]_, A_, x_, y_, function_, foreignFunction_] :=
Module[
  {d = Length[A], failed = False, aMemory, xMemory, yMemory, zetaMemory, epsteinZetaObject, realPart, imagPart},

  (* Dimensional error handling *)
  If[Length[x] != d, Message[function::dimerrx, d, Length[x], x, y, A]; failed = True];
  If[Length[y] != d, Message[function::dimerry, d, Length[y], x, y, A]; failed = True];
  If[failed, Return[$Failed]];

  aMemory = RawMemoryAllocate["CDouble", d * d];
  xMemory = RawMemoryAllocate["CDouble", d];
  yMemory = RawMemoryAllocate["CDouble", d];

  Table[RawMemoryWrite[xMemory, N[x[[i]]], i-1], {i, 1, d}];
  Table[RawMemoryWrite[yMemory, N[y[[i]]], i-1], {i, 1, d}];
  Table[RawMemoryWrite[aMemory, N[A[[i,j]]], d*(i-1)+(j-1)], {i, 1, d}, {j, 1, d}];

  zetaMemory = RawMemoryAllocate["CDouble", 2];
  epsteinZetaObject = foreignFunction[zetaMemory, N[\[Nu]], d, aMemory, xMemory, yMemory];

  realPart = RawMemoryRead[zetaMemory, 0];
  imagPart = RawMemoryRead[zetaMemory, 1];

  If[PossibleZeroQ@epsteinZetaObject,
    If[!NaNQ[N@realPart] && !NaNQ[N@imagPart],
      realPart + I*imagPart,
      ComplexInfinity
    ],
    Print["Error: Calculation failed"];
    Print["Input parameters: \[Nu] = ", \[Nu], ", A = ", A, ", x = ", x, ", y = ", y];
    Print["Dimension: ", d];
    "An Error occurred.";
    epsteinZetaObject
  ]
]

epsteinZetaDerivativesInternal[\[Nu]_, A_, x_, y_, \[Alpha]_, function_, foreignFunction_] :=
Module[
  {d = Length[A], failed = False, aMemory, xMemory, yMemory, \[Alpha]Memory, zetaMemory, epsteinZetaObject, realPart, imagPart},

  (* Dimensional error handling *)
  If[Length[x] != d, Message[function::dimerrx, d, Length[x], x, y, A, \[Alpha]]; failed = True];
  If[Length[y] != d, Message[function::dimerry, d, Length[y], x, y, A, \[Alpha]]; failed = True];
  If[Length[\[Alpha]] != d, Message[function::dimerr\[Alpha], d, Length[\[Alpha]], x, y, A, \[Alpha]]; failed = True];
  If[failed, Return[$Failed]];

  aMemory = RawMemoryAllocate["CDouble", d * d];
  xMemory = RawMemoryAllocate["CDouble", d];
  yMemory = RawMemoryAllocate["CDouble", d];
  \[Alpha]Memory = RawMemoryAllocate["CUnsignedInt", d];

  Table[RawMemoryWrite[xMemory, N[x[[i]]], i-1], {i, 1, d}];
  Table[RawMemoryWrite[yMemory, N[y[[i]]], i-1], {i, 1, d}];
  Table[RawMemoryWrite[aMemory, N[A[[i,j]]], d*(i-1)+(j-1)], {i, 1, d}, {j, 1, d}];
  Table[RawMemoryWrite[\[Alpha]Memory, \[Alpha][[i]], i-1], {i, 1, d}];

  zetaMemory = RawMemoryAllocate["CDouble", 2];
  epsteinZetaObject = foreignFunction[zetaMemory, N[\[Nu]], d, aMemory, xMemory, yMemory, \[Alpha]Memory];

  realPart = RawMemoryRead[zetaMemory, 0];
  imagPart = RawMemoryRead[zetaMemory, 1];

  If[PossibleZeroQ@epsteinZetaObject,
    If[!NaNQ[N@realPart] && !NaNQ[N@imagPart],
      realPart + I*imagPart,
      ComplexInfinity
    ],
    Print["Error: Calculation failed"];
    Print["Input parameters: \[Nu] = ", \[Nu], ", A = ", A, ", x = ", x, ", y = ", y, ", \[Alpha] = ", \[Alpha]];
    Print["Dimension: ", d];
    "An Error occurred.";
    epsteinZetaObject
  ]
]

GBesselInternal[\[Nu]_, k_, r_, function_, foreignFunction_] :=
Module[
  {d = Length[k], failed = False, rMemory, kMemory, gMemory, gObject, realPart, imagPart},

  (* Dimensional error handling *)
  If[Length[r] != d, Message[function::dimerrr, d, Length[r], k, r]; failed = True];
  If[failed, Return[$Failed]];

  kMemory = RawMemoryAllocate["CDouble", d];
  rMemory = RawMemoryAllocate["CDouble", d];

  Table[RawMemoryWrite[kMemory, N[k[[i]]], i-1], {i, 1, d}];
  Table[RawMemoryWrite[rMemory, N[r[[i]]], i-1], {i, 1, d}];

  gMemory = RawMemoryAllocate["CDouble", 2];
  gObject = foreignFunction[gMemory, N[\[Nu]], d, kMemory, rMemory];

  realPart = RawMemoryRead[gMemory, 0];
  imagPart = RawMemoryRead[gMemory, 1];

  If[PossibleZeroQ@gObject,
    If[!NaNQ[N@realPart] && !NaNQ[N@imagPart],
      realPart + I*imagPart,
      ComplexInfinity
    ],
    Print["Error: Calculation failed"];
    Print["Input parameters: \[Nu] = ", \[Nu], ", k = ", k, ", r = ", r];
    Print["Dimension: ", d];
    "An Error occurred.";
    gObject
  ]
]



(* Dimensional error handling messages *)
EpsteinZeta::dimerrx = "Input vector x = `3` has incorrect dimension. Expected dimension `1` (matching `1`\[Times]`1` matrix A = `5`), but got `2`."
EpsteinZeta::dimerry = "Input vector y = `4` has incorrect dimension. Expected dimension `1` (matching `1`\[Times]`1` matrix A = `5`), but got `2`."
EpsteinZetaReg::dimerrx = "Input vector x = `3` has incorrect dimension. Expected dimension `1` (matching `1`\[Times]`1` matrix A = `5`), but got `2`."
EpsteinZetaReg::dimerry = "Input vector y = `4` has incorrect dimension. Expected dimension `1` (matching `1`\[Times]`1` matrix A = `5`), but got `2`."
SetZetaDer::dimerrx = "Input vector x = `3` has incorrect dimension. Expected dimension `1` (matching `1`\[Times]`1` matrix A = `5`), but got `2`."
SetZetaDer::dimerry = "Input vector y = `4` has incorrect dimension. Expected dimension `1` (matching `1`\[Times]`1` matrix A = `5`), but got `2`."
SetZetaDer::dimerr\[Alpha] = "Input vector \[Alpha] = `6` has incorrect dimension. Expected dimension `1` (matching `1`\[Times]`1` matrix A = `5`), but got `2`."
EpsteinZetaRegDer::dimerrx = "Input vector x = `3` has incorrect dimension. Expected dimension `1` (matching `1`\[Times]`1` matrix A = `5`), but got `2`."
EpsteinZetaRegDer::dimerry = "Input vector y = `4` has incorrect dimension. Expected dimension `1` (matching `1`\[Times]`1` matrix A = `5`), but got `2`."
EpsteinZetaRegDer::dimerr\[Alpha] = "Input vector \[Alpha] = `6` has incorrect dimension. Expected dimension `1` (matching `1`\[Times]`1` matrix A = `5`), but got `2`."
GBessel::dimerrx = "Input vector r = `4` has incorrect dimension. Expected dimension `1` (matching k = `3`), but got `2`."


(* Define the public Epstein and set zeta functions *)
EpsteinZeta[\[Nu]_?NumericQ, A_/;MatrixQ[A] && AllTrue[Flatten[A], NumericQ], x_/;VectorQ[x] && AllTrue[x, NumericQ], y_/;VectorQ[y] && AllTrue[y, NumericQ]] := epsteinZetaInternal[\[Nu], A, x, y, EpsteinZeta, foreignFunctionEpsteinZeta]
EpsteinZetaReg[\[Nu]_?NumericQ, A_/;MatrixQ[A] && AllTrue[Flatten[A], NumericQ], x_/;VectorQ[x] && AllTrue[x, NumericQ], y_/;VectorQ[y] && AllTrue[y, NumericQ]] := epsteinZetaInternal[\[Nu], A, x, y, EpsteinZetaReg, foreignFunctionEpsteinZetaReg]
SetZetaDer[\[Nu]_?NumericQ, A_/;MatrixQ[A] && AllTrue[Flatten[A], NumericQ], x_/;VectorQ[x] && AllTrue[x, NumericQ], y_/;VectorQ[y] && AllTrue[y, NumericQ],  \[Alpha]_/;VectorQ[\[Alpha]] && AllTrue[\[Alpha], Element[#, NonNegativeIntegers]&]] := epsteinZetaDerivativesInternal[\[Nu], A, x, y, \[Alpha], SetZetaDer, foreignFunctionSetZetaDer]
EpsteinZetaRegDer[\[Nu]_?NumericQ, A_/;MatrixQ[A] && AllTrue[Flatten[A], NumericQ], x_/;VectorQ[x] && AllTrue[x, NumericQ], y_/;VectorQ[y] && AllTrue[y, NumericQ],  \[Alpha]_/;VectorQ[\[Alpha]] && AllTrue[\[Alpha], Element[#, NonNegativeIntegers]&]] := epsteinZetaDerivativesInternal[\[Nu], A, x, y, \[Alpha], EpsteinZetaRegDer, foreignFunctionEpsteinZetaRegDer]
GBessel[\[Nu]_?NumericQ,k_/;VectorQ[k] && AllTrue[k, NumericQ], r_/;VectorQ[r] && AllTrue[r, NumericQ]] := GBesselInternal[\[Nu], k, r, GBessel, foreignFunctionGBessel]


(* Check if package loaded successfully *)
epsteinLoad::info = "`1`";
If[libPath =!= $Failed &&
   Head[foreignFunctionEpsteinZeta] === ForeignFunction &&
   Head[foreignFunctionEpsteinZetaReg] === ForeignFunction &&
   Head[foreignFunctionSetZetaDer] === ForeignFunction &&
   Head[foreignFunctionEpsteinZetaRegDer] === ForeignFunction &&
   PossibleZeroQ[EpsteinZeta[-2, {{1}}, {1}, {0}]] &&
   PossibleZeroQ[EpsteinZetaReg[-2, {{1}}, {1}, {0}]]&&
   PossibleZeroQ[SetZetaDer[-2, {{1}}, {1}, {0}, {0}]]&&
   PossibleZeroQ[EpsteinZetaRegDer[-2, {{1}}, {1}, {0}, {0}]],
  Message[epsteinLoad::info,
   "The (regularized) Epstein zeta and the derivatives of the set zeta function can be called using:
	EpsteinZeta[\[Nu], A, x, y]
	EpsteinZetaReg[\[Nu], A, x, y]
	SetZetaDer[\[Nu], A, x, y, \[Alpha]]
	EpsteinZetaReg[\[Nu], A, x, y, \[Alpha]]

	\[Nu] is a real number
	A is a square matrix
	x and y are vectors of the same dimension as A
	\[Alpha] is a non-negative multi-index of the same dimension as A

The incomplete Bessel function can be called using
	GBessel[\[Nu],k,r]

	\[Nu] is a real number and
	k and r are vectors of the same dimension
"]
]


End[];


EndPackage[];
