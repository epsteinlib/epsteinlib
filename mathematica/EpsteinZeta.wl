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


(* Get the current file directory *)
currentDir = DirectoryName[$InputFileName];
parentDir = DirectoryName[currentDir];

(* Get relevant subdirectories*)
subdirsLimited = FileNameJoin[{parentDir, #}] & /@ {"result/lib", "build/src"} // Select[DirectoryQ];

(* Add current directory and limited subdirectories to $LibraryPath *)
$LibraryPath = Join[{currentDir}, {parentDir}, subdirsLimited, $LibraryPath];

(* Now try to find the library *)
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


foreignFunctionEpsteinZeta = ForeignFunctionLoad[libPath, "epstein_zeta_mathematica_call", {"RawPointer"::["CDouble"], "CDouble", "CInt", "RawPointer"::["CDouble"], "RawPointer"::["CDouble"], "RawPointer"::["CDouble"]} -> "CInt"]
If[Head[foreignFunctionEpsteinZeta] =!= ForeignFunction,
  Print["ForeignFunctionLoad for epstein_zeta_mathematica_call failed."]
]

foreignFunctionEpsteinZetaReg = ForeignFunctionLoad[libPath, "epstein_zeta_reg_mathematica_call", {"RawPointer"::["CDouble"], "CDouble", "CInt", "RawPointer"::["CDouble"], "RawPointer"::["CDouble"], "RawPointer"::["CDouble"]} -> "CInt"]
If[Head[foreignFunctionEpsteinZetaReg] =!= ForeignFunction,
  Print["ForeignFunctionLoad for epstein_zeta_reg_mathematica_call failed."]
]


(* Internal function to compute the Epstein zeta function *)
epsteinZetaCInternal[\[Nu]External_, aExternal_, xExternal_, yExternal_, foreignFunctionExternal_] :=
Module[
  {\[Nu]=\[Nu]External, a=aExternal, x = xExternal, y=yExternal, foreignFunction = foreignFunctionExternal,
  aMemory, xMemory, yMemory, dim, buffer, zetaMemory, epsteinZetaObject},

  dim = Length@a;

  aMemory = RawMemoryAllocate["CDouble", dim * dim];
  xMemory = RawMemoryAllocate["CDouble", dim];
  yMemory = RawMemoryAllocate["CDouble", dim];

  Table[RawMemoryWrite[xMemory, N[x[[i]]], i-1], {i, 1, dim}];
  Table[RawMemoryWrite[yMemory, N[y[[i]]], i-1], {i, 1, dim}];
  Table[RawMemoryWrite[aMemory, N[a[[i,j]]], dim*(i-1)+(j-1)], {i, 1, dim}, {j, 1, dim}];

  zetaMemory = RawMemoryAllocate["CDouble", 2];
  epsteinZetaObject = foreignFunction[zetaMemory, N[\[Nu]], dim, aMemory, xMemory, yMemory];

  If[PossibleZeroQ@epsteinZetaObject,
    RawMemoryRead[zetaMemory, 0] + I*RawMemoryRead[zetaMemory, 1],
    Print["Error: Calculation failed"];
    Print["Input parameters: \[Nu] = ", \[Nu], ", a = ", a, ", x = ", x, ", y = ", y];
    Print["Dimension: ", dim];
    "An Error occurred.";
    epsteinZetaObject
  ]
]


(* Define the public Epstein zeta functions *)
EpsteinZeta[\[Nu]_, A_, x_, y_] := epsteinZetaCInternal[\[Nu], A, x, y, foreignFunctionEpsteinZeta]
EpsteinZetaReg[\[Nu]_, A_, x_, y_] := epsteinZetaCInternal[\[Nu], A, x, y, foreignFunctionEpsteinZetaReg]


If[libPath =!= $Failed &&
   Head[foreignFunctionEpsteinZeta] === ForeignFunction &&
   Head[foreignFunctionEpsteinZetaReg] === ForeignFunction,
  Print["The (regularized) Epstein zeta function can be called using:"];
  Print["  EpsteinZeta[\[Nu], A, x, y]"];
  Print["  EpsteinZetaReg[\[Nu], A, x, y]"];
  Print["Where:"];
  Print["  \[Nu] is a real number"];
  Print["  A is a square matrix"];
  Print["  x and y are vectors of the same dimension as A"]
]


End[];


EndPackage[];
