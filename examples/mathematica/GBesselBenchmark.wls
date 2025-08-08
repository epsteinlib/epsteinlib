#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
   SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
   SPDX-License-Identifier: AGPL-3.0-only *)


(* ::Section:: *)
(*GBessel Benchmarks*)


(* ::Subsection:: *)
(*Setup*)


(* ::Text:: *)
(*Load Library*)


SetDirectory[NotebookDirectory[]<>"../../mathematica/"];
<<"EpsteinZeta.wl"


(* ::Text:: *)
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


(* ::Input:: *)
(*GBessel[-2,{0},{0}]*)
(*GBesselRef[-2,{0},{0}]*)
(**)
(*GBessel[-2,{0},{1}]*)
(*GBesselRef[-2,{0},{1}]*)
(**)
(*GBessel[-2,{1},{0}]*)
(*GBesselRef[-2,{1},{0}]*)
(**)
(*GBessel[-2.,{0.2},{0.1}]*)
(*GBesselRef[-2,{2/10},{1/10}]*)
(**)
(*GBessel[-2,{0.1},{0.2}]*)
(*GBesselRef[-2,{1/10},{2/10}]*)
(**)
(*GBessel[-2,{2},{1}]*)
(*GBesselRef[-2,{2},{1}]*)
(**)
(*GBessel[-2,{1},{2}]*)
(*GBesselRef[-2,{1},{2}]*)
(**)
(*GBessel[-2,{1},{2}]*)
(*GBesselRef[-2,{1},{2}]*)
(**)
(*GBessel[-2,{2},{1}]*)
(*GBesselRef[-2,{2},{1}]*)


(* ::Subsection:: *)
(*Paper Benchmarks*)
