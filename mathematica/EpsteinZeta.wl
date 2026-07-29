(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024-2026 Jonathan Busse <jonathan@jbusse.de>
   SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


BeginPackage["EpsteinZeta`"];


EpsteinZeta::usage="EpsteinZeta[\[Nu],A,x,y] computes the Epstein zeta function sum_{z in \[CapitalLambda]} exp(- 2 Pi I y.z)/|z - x|^\[Nu]
using the algorithm in Crandall, R., Unified algorithms for polylogarithm, L-series, and zeta variants. Algorithmic Reflections: Selected Works. PSIpress (2012).";
EpsteinZetaReg::usage="EpsteinZetaReg[\[Nu],A,x,y] computes the regularized Epstein zeta function exp(2 Pi I x.y) sum_{z in \[CapitalLambda]} exp(- 2 Pi I y.z)/|z - x|^\[Nu]  - s\:0302(y)/|det(\[CapitalLambda])|
as in Andreas A. Buchheit et al., Exact continuum representation of long-range interacting systems and emerging exotic phases in unconventional superconductors. Phys. Rev. Res. 5 (4 Oct. 2023), p. 043065.  using a modification of the algorithm in Crandall, R., Unified algorithms for polylogarithm, L-series, and zeta variants. Algorithmic Reflections: Selected Works. PSIpress (2012).";
EpsteinZetaAniso::usage="EpsteinZetaAniso[\[Nu],A,x,y,\[Alpha]] computes the anisotropic Epstein zeta function sum_{z in \[CapitalLambda]} (z-x)^\[Alpha] exp(- 2 Pi I y.z)/|z-x|^\[Nu].";
EpsteinZetaAnisoReg::usage="EpsteinZetaAnisoReg[\[Nu],A,x,y,\[Alpha]] computes the regularized anisotropic Epstein zeta function sum_{z in \[CapitalLambda]- x} z^\[Alpha] exp(- 2 Pi I y.z)/|z|^\[Nu]  - s\:0302(y)/((-2 \[Pi] I)^|\[Alpha]| |det(\[CapitalLambda])|)."


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

foreignFunctionEpsteinZetaAniso = ForeignFunctionLoad[libPath, "epstein_zeta_aniso_mathematica_call", {"RawPointer"::["CDouble"], "CDouble", "CInt", "RawPointer"::["CDouble"], "RawPointer"::["CDouble"], "RawPointer"::["CDouble"], "RawPointer"::["CUnsignedInt"]} -> "CInt"]
If[Head[foreignFunctionEpsteinZetaAniso] =!= ForeignFunction,
  Print["ForeignFunctionLoad for epstein_zeta_aniso_mathematica_call failed."]
]

foreignFunctionEpsteinZetaAnisoReg = ForeignFunctionLoad[libPath, "epstein_zeta_aniso_reg_mathematica_call", {"RawPointer"::["CDouble"], "CDouble", "CInt", "RawPointer"::["CDouble"], "RawPointer"::["CDouble"], "RawPointer"::["CDouble"], "RawPointer"::["CUnsignedInt"]} -> "CInt"]
If[Head[foreignFunctionEpsteinZetaAnisoReg] =!= ForeignFunction,
  Print["ForeignFunctionLoad for epstein_zeta_aniso_reg_mathematica_call failed."]
]

(* For checking nan output of C function *)
NaNQ = ResourceFunction["NaNQ"];


(* Internal routines for C function access *)
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

  realPart = RawMemoryRead[zetaMemory, 0];
  imagPart = RawMemoryRead[zetaMemory, 1];

  If[PossibleZeroQ@status,
    If[!NaNQ[N@realPart] && !NaNQ[N@imagPart],
      realPart + I*imagPart,
      ComplexInfinity
    ],
    Message[function::cfail, status, \[Nu], A, x, y];
    $Failed
  ]
]

epsteinZetaDerivativesInternal[\[Nu]_, A_, x_, y_, \[Alpha]_, function_, foreignFunction_] :=
Module[
  {d = Length[A], aMemory, xMemory, yMemory, \[Alpha]Memory, zetaMemory, status, realPart, imagPart},

  aMemory = RawMemoryAllocate["CDouble", d * d];
  xMemory = RawMemoryAllocate["CDouble", d];
  yMemory = RawMemoryAllocate["CDouble", d];
  \[Alpha]Memory = RawMemoryAllocate["CUnsignedInt", d];

  Table[RawMemoryWrite[xMemory, N[x[[i]]], i-1], {i, 1, d}];
  Table[RawMemoryWrite[yMemory, N[y[[i]]], i-1], {i, 1, d}];
  Table[RawMemoryWrite[aMemory, N[A[[i,j]]], d*(i-1)+(j-1)], {i, 1, d}, {j, 1, d}];
  Table[RawMemoryWrite[\[Alpha]Memory, \[Alpha][[i]], i-1], {i, 1, d}];

  zetaMemory = RawMemoryAllocate["CDouble", 2];
  status = foreignFunction[zetaMemory, N[\[Nu]], d, aMemory, xMemory, yMemory, \[Alpha]Memory];

  realPart = RawMemoryRead[zetaMemory, 0];
  imagPart = RawMemoryRead[zetaMemory, 1];

  If[PossibleZeroQ@status,
    If[!NaNQ[N@realPart] && !NaNQ[N@imagPart],
      realPart + I*imagPart,
      ComplexInfinity
    ],
    Message[function::cfail, status, \[Nu], A, x, y, \[Alpha]];
    $Failed
  ]
]


(* Error messages *)
EpsteinZeta::cfail = "The library call failed with status `1` for \[Nu] = `2`, A = `3`, x = `4`, y = `5`.";
EpsteinZetaReg::cfail = EpsteinZeta::cfail;
EpsteinZetaAniso::cfail = "The library call failed with status `1` for \[Nu] = `2`, A = `3`, x = `4`, y = `5`, \[Alpha] = `6`.";
EpsteinZetaAnisoReg::cfail = EpsteinZetaAniso::cfail;


(* Check arguments helpers *)
numericSquareMatrixQ[A_] := SquareMatrixQ[A] && MatrixQ[A, NumericQ];
numericVectorQ[v_] := VectorQ[v, NumericQ];
multiIndexQ[\[Alpha]_] := VectorQ[\[Alpha], IntegerQ[#] && NonNegative[#] &];


(* Define the public Epstein zeta functions *)
EpsteinZeta[\[Nu]_?NumericQ, A_?numericSquareMatrixQ, x_?numericVectorQ, y_?numericVectorQ] /;
    Length[x] == Length[y] == Length[A] :=
  epsteinZetaInternal[\[Nu], A, x, y, EpsteinZeta, foreignFunctionEpsteinZeta]

EpsteinZetaReg[\[Nu]_?NumericQ, A_?numericSquareMatrixQ, x_?numericVectorQ, y_?numericVectorQ] /;
    Length[x] == Length[y] == Length[A] :=
  epsteinZetaInternal[\[Nu], A, x, y, EpsteinZetaReg, foreignFunctionEpsteinZetaReg]

EpsteinZetaAniso[\[Nu]_?NumericQ, A_?numericSquareMatrixQ, x_?numericVectorQ, y_?numericVectorQ, \[Alpha]_?multiIndexQ] /;
    Length[x] == Length[y] == Length[\[Alpha]] == Length[A] :=
  epsteinZetaDerivativesInternal[\[Nu], A, x, y, \[Alpha], EpsteinZetaAniso, foreignFunctionEpsteinZetaAniso]

EpsteinZetaAnisoReg[\[Nu]_?NumericQ, A_?numericSquareMatrixQ, x_?numericVectorQ, y_?numericVectorQ, \[Alpha]_?multiIndexQ] /;
    Length[x] == Length[y] == Length[\[Alpha]] == Length[A] :=
  epsteinZetaDerivativesInternal[\[Nu], A, x, y, \[Alpha], EpsteinZetaAnisoReg, foreignFunctionEpsteinZetaAnisoReg]
(* Check if package loaded successfully *)
epsteinLoad::info = "`1`";
epsteinLoad::broken = "The Epstein zeta library at `1` loaded but failed its self-test. \
Results cannot be trusted! Check that the library is current and built for this platform.";

If[libPath =!= $Failed &&
   Head[foreignFunctionEpsteinZeta] === ForeignFunction &&
   Head[foreignFunctionEpsteinZetaReg] === ForeignFunction &&
   Head[foreignFunctionEpsteinZetaAniso] === ForeignFunction &&
   Head[foreignFunctionEpsteinZetaAnisoReg] === ForeignFunction,

  If[PossibleZeroQ[EpsteinZeta[-2, {{1}}, {1}, {0}]] &&
     PossibleZeroQ[EpsteinZetaReg[-2, {{1}}, {1}, {0}]] &&
     PossibleZeroQ[EpsteinZetaAniso[-2, {{1}}, {1}, {0}, {0}]] &&
     PossibleZeroQ[EpsteinZetaAnisoReg[-2, {{1}}, {1}, {0}, {0}]] &&
     EpsteinZeta[1, {{1}}, {0}, {0}] =!= EpsteinZetaReg[1, {{1}}, {0}, {0}] &&
     EpsteinZetaAniso[1, {{1}}, {0}, {0}, {0}] =!= EpsteinZetaAnisoReg[1, {{1}}, {0}, {0}, {0}],

    Message[epsteinLoad::info,
     "The (regularized) Epstein zeta and the (regularized) anisotropic Epstein zeta function can be called using:
	EpsteinZeta[\[Nu], A, x, y]
	EpsteinZetaReg[\[Nu], A, x, y]
	EpsteinZetaAniso[\[Nu], A, x, y, \[Alpha]]
	EpsteinZetaAnisoReg[\[Nu], A, x, y, \[Alpha]]

	\[Nu] is a real number
	A is a square matrix
	x and y are vectors of the same dimension as A
	\[Alpha] is a non-negative multi-index of the same dimension as A
"],

    Message[epsteinLoad::broken, libPath]]
]


End[];


EndPackage[];
