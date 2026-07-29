(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
   SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


BeginPackage["EpsteinZeta`"];


EpsteinZeta::usage="epsteinZeta[\[Nu],A,x,y] computes the Epstein zeta function sum_{z in Lambda} exp(- 2 Pi I y.z)/|z - x|^\[Nu]
using the algorithm in Crandall, R., Unified algorithms for polylogarithm, L-series, and zeta variants. Algorithmic Reflections: Selected Works. PSIpress (2012).";
EpsteinZetaReg::usage="epsteinZetaReg[\[Nu],A,x,y] computes the regularized Epstein zeta function exp(2 Pi I x.y) sum_{z in Lambda} exp(- 2 Pi I y.(z - x))/|z - x|^\[Nu]  - s\:0302(y)/|det(\[CapitalLambda])|
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


(* True when the library wrote something that is not a usable machine
   double (NaN, inf, or an unconvertible read). *)
badReadQ[u_] := ! MachineNumberQ[u];


(* Internal routine for C function access *)
epsteinZetaInternal[\[Nu]_, A_, x_, y_, function_, foreignFunction_] :=
Module[
  {d = Length[A], aMemory, xMemory, yMemory, zetaMemory, status, realPart, imagPart},

  aMemory = RawMemoryAllocate["CDouble", d * d];
  xMemory = RawMemoryAllocate["CDouble", d];
  yMemory = RawMemoryAllocate["CDouble", d];

  Table[RawMemoryWrite[xMemory, N[x[[i]]], i-1], {i, 1, d}];
  Table[RawMemoryWrite[yMemory, N[y[[i]]], i-1], {i, 1, d}];
  Table[RawMemoryWrite[aMemory, N[A[[i,j]]], d*(i-1)+(j-1)], {i, 1, d}, {j, 1, d}];

  zetaMemory = RawMemoryAllocate["CDouble", 2];
  status = foreignFunction[zetaMemory, N[\[Nu]], d, aMemory, xMemory, yMemory];

  Which[
    status === $Failed,
      Message[function::marshal, \[Nu], A, x, y]; $Failed,
    status =!= 0,
      Message[function::cfail, status, \[Nu], A, x, y]; $Failed,
    True,
      realPart = RawMemoryRead[zetaMemory, 0];
      imagPart = RawMemoryRead[zetaMemory, 1];
      If[badReadQ[realPart] || badReadQ[imagPart],
        ComplexInfinity,
        realPart + I*imagPart]
   ]
]


(* Dimensional error handling messages *)
EpsteinZeta::cfail = "The library returned error status `1` for \[Nu] = `2`, A = `3`, x = `4`, y = `5`.";
EpsteinZeta::marshal = "The arguments \[Nu] = `1`, A = `2`, x = `3`, y = `4` could not be passed to the library. \
\[Nu] and every entry of A, x and y must be real and representable as a machine double.";
EpsteinZetaReg::cfail = EpsteinZeta::cfail;
EpsteinZetaReg::marshal = EpsteinZeta::marshal;


(* *)
numericVectorQ[v_] := VectorQ[v, NumericQ];
numericSquareMatrixQ[A_] := SquareMatrixQ[A] && MatrixQ[A, NumericQ];


(* Define the public Epstein zeta functions *)
EpsteinZeta[\[Nu]_?NumericQ, A_?numericSquareMatrixQ, x_?numericVectorQ, y_?numericVectorQ] /;
    Length[x] == Length[y] == Length[A] :=
  epsteinZetaInternal[\[Nu], A, x, y, EpsteinZeta, foreignFunctionEpsteinZeta]

EpsteinZetaReg[\[Nu]_?NumericQ, A_?numericSquareMatrixQ, x_?numericVectorQ, y_?numericVectorQ] /;
    Length[x] == Length[y] == Length[A] :=
  epsteinZetaInternal[\[Nu], A, x, y, EpsteinZetaReg, foreignFunctionEpsteinZetaReg]



epsteinLoad::info = "`1`";
epsteinLoad::selftest = "EpsteinZeta` loaded but failed its self-test: `1`. \
The library is returning incorrect results; do not trust its output.";

If[libPath =!= $Failed &&
   Head[foreignFunctionEpsteinZeta] === ForeignFunction &&
   Head[foreignFunctionEpsteinZetaReg] === ForeignFunction,

  With[{failed = Keys @ Select[<|
      "EpsteinZeta[-2,{{1}},{1},{0}] vanishes"    -> PossibleZeroQ[EpsteinZeta[-2, {{1}}, {1}, {0}]],
      "EpsteinZetaReg[-2,{{1}},{1},{0}] vanishes" -> PossibleZeroQ[EpsteinZetaReg[-2, {{1}}, {1}, {0}]],
      "EpsteinZeta has a pole at \[Nu] = d"       -> (EpsteinZeta[1, {{1}}, {0}, {0}] === ComplexInfinity),
      "EpsteinZetaReg is finite at \[Nu] = d"     -> MachineNumberQ[EpsteinZetaReg[1, {{1}}, {0}, {0}]]
    |>, Not]},

    If[failed === {},
      Message[epsteinLoad::info,
        "The (regularized) Epstein zeta function can be called using:
  EpsteinZeta[\[Nu], A, x, y]
  EpsteinZetaReg[\[Nu], A, x, y]
Where:
  \[Nu] is a real number
  A is a d\[Times]d square matrix
  x and y are vectors of length d"],
      Message[epsteinLoad::selftest, failed]]]
]

End[];


EndPackage[];
