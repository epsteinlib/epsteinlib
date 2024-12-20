(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
   SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


BeginPackage["EpsteinZeta`"];


EpsteinZeta::usage="epsteinZeta[\[Nu],A,x,y] computes the Epstein zeta function sum_{z in Lambda} exp(- 2 Pi I y.z)/|z - x|^\[Nu]
using the algorithm in Crandall, R., Unified algorithms for polylogarithm, L-series, and zeta variants. Algorithmic Reflections: Selected Works. PSIpress (2012).";
EpsteinZetaReg::usage="epsteinZetaReg[\[Nu],A,x,y] computes the regularised Epstein zeta function exp(2 Pi I x.y) sum_{z in Lambda} exp(- 2 Pi I y.(z - x))/|z - x|^\[Nu]  - s\:0302(y)/|det(\[CapitalLambda])|
as in Andreas A. Buchheit et al., Exact continuum representation of long-range interacting systems and emerging exotic phases in unconventional superconductors. Phys. Rev. Res. 5 (4 Oct. 2023), p. 043065.  using a modification of the algorithm in Crandall, R., Unified algorithms for polylogarithm, L-series, and zeta variants. Algorithmic Reflections: Selected Works. PSIpress (2012).";


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

(* For checking nan output of C function *)
NaNQ = ResourceFunction["NaNQ"];

(* Internal routine for C function access *)
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


(* Dimensional error handling messages *)
EpsteinZeta::dimerrx = "Input vector x = `3` has incorrect dimension. Expected dimension `1` (matching `1`×`1` matrix A = `5`), but got `2`."
EpsteinZeta::dimerry = "Input vector y = `4` has incorrect dimension. Expected dimension `1` (matching `1`×`1` matrix A = `5`), but got `2`."
EpsteinZetaReg::dimerrx = "Input vector x = `3` has incorrect dimension. Expected dimension `1` (matching `1`×`1` matrix A = `5`), but got `2`."
EpsteinZetaReg::dimerry = "Input vector y = `4` has incorrect dimension. Expected dimension `1` (matching `1`×`1` matrix A = `5`), but got `2`."


(* Define the public Epstein zeta functions *)
EpsteinZeta[\[Nu]_?NumericQ, A_/;MatrixQ[A] && AllTrue[Flatten[A], NumericQ], x_/;VectorQ[x] && AllTrue[x, NumericQ], y_/;VectorQ[y] && AllTrue[y, NumericQ]] := epsteinZetaInternal[\[Nu], A, x, y, EpsteinZeta, foreignFunctionEpsteinZeta]

EpsteinZetaReg[\[Nu]_?NumericQ, A_/;MatrixQ[A] && AllTrue[Flatten[A], NumericQ], x_/;VectorQ[x] && AllTrue[x, NumericQ], y_/;VectorQ[y] && AllTrue[y, NumericQ]] := epsteinZetaInternal[\[Nu], A, x, y, EpsteinZetaReg, foreignFunctionEpsteinZetaReg]


(* Check if package loaded successfully *)
If[libPath =!= $Failed &&
   Head[foreignFunctionEpsteinZeta] === ForeignFunction &&
   Head[foreignFunctionEpsteinZetaReg] === ForeignFunction &&
   PossibleZeroQ[EpsteinZeta[-2, {{1}}, {1}, {0}]] &&
   PossibleZeroQ[EpsteinZetaReg[-2, {{1}}, {1}, {0}]],
  Print["The (regularized) Epstein zeta function can be called using:
  EpsteinZeta[\[Nu], A, x, y]
  EpsteinZetaReg[\[Nu], A, x, y]
Where:
  \[Nu] is a real number
  A is a square matrix
  x and y are vectors of the same dimension as A"]
]

End[];


EndPackage[];
