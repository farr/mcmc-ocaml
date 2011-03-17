(* evidence_tool.ml: Tool to perform direct integration evidence calcs. *)
(* Copyright (C) 2011 Will M. Farr <w-farr@northwestern.edu> *)

(* This program is free software: you can redistribute it and/or modify *)
(* it under the terms of the GNU General Public License as published by *)
(* the Free Software Foundation, either version 3 of the License, or *)
(* (at your option) any later version. *)

(* This program is distributed in the hope that it will be useful, *)
(* but WITHOUT ANY WARRANTY; without even the implied warranty of *)
(* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the *)
(* GNU General Public License for more details. *)

(* You should have received a copy of the GNU General Public License *)
(* along with this program.  If not, see <http://www.gnu.org/licenses/>. *)

(** Reads output of an MCMC on stdin, and prints three estimators of
    the evidence on stdout: the Harmonic mean evidence, the Lebesgue
    evidence, and the direct-integration evidence (in that order).
    The MCMC parameters from the input are treated as an array of
    floats (i.e. no conversion is performed from coordinates to
    params). *)

open Printf

let nbox = ref 64
let leb_eps = ref 0.1

let options = 
  [("-nbox", Arg.Set_int nbox,
    sprintf "the number of samples in each direct/Lebesgue integration box (default %d)" !nbox);
   ("-lebeps", Arg.Set_float leb_eps,
    sprintf "the truncation parameter eps in Lebesgue evidence (default %g)" !leb_eps)]

module Ev = Evidence.Make(
  struct
    type params = float array

    let to_coords (p : params) = p
  end)

let _ = 
  Arg.parse options (fun _ -> ()) "evidence_tool.{native,byte} OPTIONS ...";
  let samples = Read_write.read (fun pt -> pt) stdin in 
    printf "%g %g %g\n" 
      (Ev.evidence_harmonic_mean samples)
      (Ev.evidence_lebesgue ~n:(!nbox) ~eps:(!leb_eps) samples)
      (Ev.evidence_direct ~n:(!nbox) samples)
