#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
   SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


SetDirectory[NotebookDirectory[]<>"../../mathematica/"];
<<"EpsteinZeta.wl"


madelungRef = -1.7475645946331821906362120355443974;


dim = 3;
A = IdentityMatrix[dim];
x = ConstantArray[0, dim];
y = ConstantArray[0.5, dim];
\[Nu] = 1.0;


madelung = Re[EpsteinZeta[\[Nu], A, x, y]];


Print["Madelung sum in 3 dimensions: ", NumberForm[madelung, 16]];
Print["Reference value:              ", NumberForm[madelungRef, 16]];
Print["Relative error:               +", ScientificForm[Abs[madelungRef - madelung]/Abs[madelungRef], 2]];
