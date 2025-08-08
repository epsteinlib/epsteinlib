#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
   SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*GBessel Benchmarks*)


(* ::Text:: *)
(*This notebook tests the implementation of the incomplete Bessel function.*)
(**)
(*	GBessel[nu,k] = 2 int_0^1 t^(-nu-1) exp(-pi k^2 / t^2) exp(-pi r^2 t^2) dt*)
(**)
(*Evaluation of this script takes less than 3 minutes on a standard laptop equipped with 12 cores.*)


(* ::Subsection::Closed:: *)
(*Setup Packages*)


(* ::Text:: *)
(*Please ensure that EpsteinLib is correctly installed and accessible before running the notebook. That is*)
(*	*)
(*	1. 	go to the project: cd epsteinlib*)
(*	2.	build library: meson setup build*)
(*	3.	compile library: meson compile -C build*)
(*	*)
(*To test the C library, run meson test -C build.*)


(* ::Text:: *)
(*Load the Mathematica wrapper of the Epstein zeta function,  such that it can be used in parallel on all Kernels.*)


SetDirectory[NotebookDirectory[]<>"../../mathematica"];
ParallelEvaluate[Quiet[<<"EpsteinZeta.wl"]];
<<"EpsteinZeta.wl";


(* ::Text:: *)
(*Install MaTeX for LATEX figure annotations. If MaTeX is not installed already un-comment the ResourceFunction line to install it.*)


(*ResourceFunction["MaTeXInstall"][];*)
Needs["MaTeX`"];


(* ::Text:: *)
(*Set options for pretty plots.*)


SetOptions[RegionPlot,TicksStyle->{FontSize->18},LabelStyle->{FontSize->18}];
SetOptions[MaTeX,"Preamble"->{"\\usepackage{xcolor}","\\usepackage{bm}"}];
SetOptions[MaTeX,FontSize->18];


(* ::Subsection::Closed:: *)
(*Utility functions*)


(* ::Input:: *)
(*Error[ref_,comp_]:=Block[{errorAbs,errorRel},errorAbs=Abs[ref-comp];errorRel=If[PossibleZeroQ@ref,1,errorAbs/Abs[ref]];Min[errorAbs,errorRel]];*)


(* ::Input:: *)
(*CharacterOverlay[int_,plot_]:=Overlay[{plot,MaTeX["("<>ToString[Alphabet[][[i]]]<>")"]},Alignment->{-1,1}]*)


(* ::Input:: *)
(*LogarithmicScaling[x_,min_,max_]:=Log[x/min]/Log[max/min]*)


(* ::Subsubsection:: *)
(*Incomplete Bessel Function*)


(* ::Subsection::Closed:: *)
(*Define reference function*)


(* ::Input:: *)
(*GBesselRef[nu_,k_,r_]:=NIntegrate[*)
(*2 t^(-nu)Exp[-Pi k . k/t^2]Exp[-Pi r . r t^2]/t*)
(*,{t,0,1}*)
(*,PrecisionGoal->40*)
(*,AccuracyGoal->40*)
(*,WorkingPrecision->40*)
(*,WorkingPrecision->60*)
(*]*)


(* ::Subsection:: *)
(*Quick Tests*)


1+1


(* ::Text:: *)
(*Small values*)


ref=GBesselRef[-2,{0},{0}]
val=GBessel[-2,{0},{0}]
Error[ref,val]


(* ::Text:: *)
(*Big values*)


ref=GBesselRef[-2,{5},{5}]
val=GBessel[-2,{5},{5}]
Error[ref,val]


(* ::Text:: *)
(*Lower Crandall*)


ref=GBesselRef[-2,{0},{1}]
val=GBessel[-2,{0},{1}]
Error[ref,val]


(* ::Text:: *)
(*Upper Crandall*)


ref=GBesselRef[-2,{1},{0}]
val=GBessel[-2,{1},{0}]
Error[ref,val]


(* ::Text:: *)
(*Taylor series*)


ref=GBesselRef[-2,{1/10},{2/10}]
val=GBessel[-2,{0.1},{0.2}]
Error[ref,val]

ref=GBesselRef[-2,{2/10},{1/10}]
val=GBessel[-2,{0.2},{0.1}]
Error[ref,val]


(* ::Text:: *)
(*Recursive algorithm*)


ref=GBesselRef[-2,{1},{2}]
val=GBessel[-2,{1},{2}]
Error[ref,val]

ref=GBesselRef[-2,{2},{1}]
val=GBessel[-2,{2},{1}]
Error[ref,val]


(* ::Subsection:: *)
(*Benchmark*)


(* ::Subsubsection:: *)
(*Generate Errors*)


(* ::Input:: *)
(*errorTab=Table[Flatten[*)
(*	ParallelTable[{k,r,Max[Error[Quiet@GBesselRef[nu,{k},{r}],N@GBessel[nu,{k},{r}]],10^(-18)]}*)
(*,{k,Range[1/100,4,1/4]}*)
(*,{r,Range[0,4,1/4]}*)
(*],1*)
(*],{nu,-Range[-2,4,2]}*)
(*];*)


(* ::Input:: *)
(*Min[N@#[[All,3]]]&/@errorTab*)
(*Max[#[[All,3]]]&/@errorTab*)


(* ::Subsubsection:: *)
(*Make Plot*)


(* ::Input:: *)
(*Min[N@#[[All,3]]]&/@errorTab*)
(*Max[#[[All,3]]]&/@errorTab*)


(* ::Input:: *)
(*besselPlot=GraphicsRow[*)
(*Join[*)
(*Table[*)
(*CharacterOverlay[i,*)
(*With[{tab=errorTab[[i]]}*)
(*,With[{min=Min@tab[[All,3]],max=10^(-15)}*)
(*,ListDensityPlot[tab*)
(*,PlotTheme->"Scientific"*)
(*,PlotRange->{{-1/20, 4+1/20}, {-1/20, 4+1/20},{min,max}}*)
(*,PlotLabel->MaTeX["\\nu="<>ToString[{2,0,-2,-4}[[i]]]]*)
(*,FrameLabel->{MaTeX@"|\\bm k|",MaTeX@"|\\bm r|"}*)
(*,PlotLegends->None*)
(*,ColorFunction->(ColorData["LakeColors"][LogarithmicScaling[#,min,max]]&)*)
(*,FrameTicks->{*)
(*{*)
(*Join[*)
(*{#,MaTeX[ToString@#],{0.014,0}}&/@Range[0,4,1]*)
(*,{#,,{0.008,0}}&/@Range[0,4,1/2]*)
(*]*)
(*,Join[*)
(*{#,,{0.014,0}}&/@Range[0,4,1]*)
(*,{#,,{0.008,0}}&/@Range[0,4,1/2]*)
(*]*)
(*}*)
(*,{*)
(*Join[*)
(*{#,MaTeX[ToString@#],{0.014,0}}&/@Range[0,4,1]*)
(*,{#,,{0.008,0}}&/@Range[0,4,1/2]*)
(*]*)
(*,Join[*)
(*{#,,{0.014,0}}&/@Range[0,4,1]*)
(*,{#,,{0.008,0}}&/@Range[0,4,1/2]*)
(*]*)
(*}*)
(*}*)
(*,ColorFunctionScaling->False*)
(*,InterpolationOrder->0*)
(*,ImageSize->200*)
(*]*)
(*]*)
(*]*)
(*],{i,Length@errorTab}*)
(*]*)
(*,{*)
(*GraphicsRow[*)
(*{*)
(*BarLegend[*)
(*{"LakeColors",{0,1}}*)
(*,Ticks->({LogarithmicScaling[10^#,10^(-18),10^(-15)],}&/@Range[-18,-14,1])*)
(*,ImageSize->215*)
(*]*)
(*,Graphics[{*)
(*White*)
(*,Style[Text[MaTeX["\\!10^{-15}"],{0.5,3/3}]]*)
(*,Style[Text[MaTeX["\\!10^{-16}"],{0.5,2/3}]]*)
(*,Style[Text[MaTeX["\\!10^{-17}"],{0.5,1/3}]]*)
(*,Style[Text[MaTeX["\\!10^{-18}"],{0.5,0}]]*)
(*,Style[Text[MaTeX["~~~~~~~~"],{0.5,0.00}]]*)
(*,Style[Text[MaTeX["~~~~~~~~"],{0.5,-0.05}]]*)
(*}*)
(*,AspectRatio->3.5*)
(*,Axes->False*)
(*,Frame->False  *)
(*,ImageSize->{70,280}*)
(*]*)
(*},Spacings->-10*)
(*]*)
(*}*)
(*],ImageSize->{1000,290},Spacings->20*)
(*]*)


(* ::Subsubsection:: *)
(*Show Plot*)


(* ::Input:: *)
(*besselPlot*)
