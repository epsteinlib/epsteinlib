#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Text:: *)
(*Call Epstein zeta function*)


SetDirectory[NotebookDirectory[]<>"../../mathematica"];
<<"EpsteinZeta.wl"


(* ::Text:: *)
(*Load EpsteinZeta.wl in every local Kernel*)


ParallelEvaluate[Quiet[<<"EpsteinZeta.wl"]];


(* ::Text:: *)
(*Install MaTeX for L AT EX figure annotations*)


(*ResourceFunction["MaTeXInstall"][];*)
Needs["MaTeX`"];
SetOptions[MaTeX,FontSize->18];


(* ::Text:: *)
(*Taylor series around zero*)


TaylorSeriesInZero[partDerInZero_,var_,order_]:=Total[Table[Apply[Times,var^\[Alpha]/\[Alpha]!] partDerInZero[\[Alpha]],{\[Alpha],Select[Tuples[Range[0,order],Length[var]],Total[#]<=order &]}]]


(* ::Subsection:: *)
(*Example 1D*)


EpsteinZeta1D[y1_]:=Re@EpsteinZeta[3,IdentityMatrix[1],{0},{y1}]


EpsteinZetaTaylor1D[y1_,0]=TaylorSeriesInZero[Re@SetZetaDer[3,IdentityMatrix[1],{0},{0.1},#]&,{y1}-{0.1},0];
EpsteinZetaTaylor1D[y1_,1]=TaylorSeriesInZero[Re@SetZetaDer[3,IdentityMatrix[1],{0},{0.1},#]&,{y1}-{0.1},1];
EpsteinZetaTaylor1D[y1_,2]=TaylorSeriesInZero[Re@SetZetaDer[3,IdentityMatrix[1],{0},{0.1},#]&,{y1}-{0.1},2];
EpsteinZetaTaylor1D[y1_,2]=TaylorSeriesInZero[Re@SetZetaDer[3,IdentityMatrix[1],{0},{0.1},#]&,{y1}-{0.1},3]


zetaTab1D = ParallelTable[{y1,EpsteinZeta1D[y1]},{y1,-0.,0.5,0.005}];
zetaTaylorTab1D = Table[ParallelTable[{y1,EpsteinZetaTaylor1D[y1,order]},{y1,-0.,0.5,0.005}],{order,0,3}];


Show[
ListLinePlot[zetaTab1D,PlotTheme->"Scientific",PlotStyle->Black]
,ListLinePlot[zetaTaylorTab1D,PlotTheme->"Scientific",PlotStyle->Opacity[0.6],PlotLegends->Placed[MaTeX[#<>"\\text{'th order}"]&/@{"0","1","2"}, {0.2,0.3}]]
]


(* ::Subsection:: *)
(*Example 2D*)


EpsteinZeta2D[y1_,y2_]:=Re@EpsteinZeta[3,IdentityMatrix[2],{0,0},{y1,y2}]


EpsteinZetaTaylor2D[y1_,y2_,0]=TaylorSeriesInZero[Re@SetZetaDer[3,IdentityMatrix[2],{0,0},{0.1,0},#]&,{y1,y2}-{0.1,0},0];
EpsteinZetaTaylor2D[y1_,y2_,1]=TaylorSeriesInZero[Re@SetZetaDer[3,IdentityMatrix[2],{0,0},{0.1,0},#]&,{y1,y2}-{0.1,0},1];
EpsteinZetaTaylor2D[y1_,y2_,2]=TaylorSeriesInZero[Re@SetZetaDer[3,IdentityMatrix[2],{0,0},{0.1,0},#]&,{y1,y2}-{0.1,0},2];
EpsteinZetaTaylor2D[y1_,y2_,3]=TaylorSeriesInZero[Re@SetZetaDer[3,IdentityMatrix[2],{0,0},{0.1,0},#]&,{y1,y2}-{0.1,0},3]


zetaTab2D = Re@Flatten[ParallelTable[{y1,y2,EpsteinZeta2D[y1, y2]},{y1,-0.,0.2,0.005},{y2,-0.1,0.1,0.005}],1];
zetaTaylorTab2D = Re@Table[Flatten[Table[{y1,y2,EpsteinZetaTaylor2D[y1, y2,order]},{y1,-0.,0.2,0.005},{y2,-0.1,0.1,0.005}],1],{order,0,3}];


Show[
ListPlot3D[zetaTab2D,PlotTheme->"Scientific",ColorFunction->"ThermometerColors",Mesh->None]
,ListPlot3D[zetaTaylorTab2D,PlotStyle->Opacity[0.1],Mesh->None]
]
