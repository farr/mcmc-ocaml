(* gaussian_cauchy_efficiency.ml: Compute efficiency of RJMCMC proposal. *)
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

module Interp = Interpolate_pdf.Make(struct
  type point = float array

  let coord (pt : point) = pt
  let point (cd : float array) = cd
end)

module Ev = Evidence.Make(struct
  type params = float array
  let to_coords (x : params) = x
end)

let _ = Random.self_init ()

let _ = Printf.eprintf "WARNING: you must modify mcmc.ml to *not* jump between models according to the prior before running this code.\n%!"

let data = [|-1.44898; -0.0762953; 2.25525; -0.284584; 1.16297; 0.00864677; 
            0.211493; -1.03406; -0.304487; 0.868774; -0.489545; -0.179075; 
            -0.139734; 1.76451; -0.979861; -3.02674; -0.946947; 0.846734; 
            0.74059; 0.630729; -1.89231; -1.29979; 0.636769; 0.948439; -1.82115; 
            -0.712503; -0.58101; -0.642658; -1.27348; -1.02178; 1.61619; 2.62544; 
            -1.26247; 0.636684; -1.88267; -0.597691; 2.81632; 1.27734; -0.683581; 
            -0.0384763; -1.05507; -0.950249; -0.504831; 0.792793; -0.301302; 
            -0.405713; 2.12917; -1.52254; -0.500314; 1.56529; -1.05244; -0.39324; 
            -2.00543; -0.393424; 0.104712; 1.83405; -1.12503; 1.39026; 0.784583; 
            -0.542822; 0.696715; 0.554867; 0.118314; -0.404154; 0.123976; 
            0.492262; 0.921631; 1.59935; -0.0245473; -0.673963; -0.0106481; 
            -0.575168; -0.442549; -0.825369; 0.364721; -1.18684; 0.186002; 
            0.707357; 0.780187; -1.41614; 1.1693; 1.06944; 0.56919; -0.849411; 
            1.1442; 0.130583; -1.03244; 0.708557; -1.26797; -1.00726; -0.894115; 
            1.23684; 1.19339; -0.419484; -0.746722; 0.566607; -0.610874; 
            -0.427252; -0.581807; 1.04941|]

let pgaussian = 1.0 /. 5e8
let pcauchy = 1.0 -. pgaussian

let _ = Printf.eprintf "pg = %g, pc = %g\n%!" pgaussian pcauchy

let mumin = -1.0
let mumax = 1.0
let sigmamin = 0.5
let sigmamax = 1.5

let log_prior params = 
  match params with 
    | [| mu; sigma |] -> 
      if mu < mumin || mu > mumax || sigma < sigmamin || sigma > sigmamax then 
        neg_infinity
      else
        -0.693147 (* log 0.5 *)
    | _ -> raise (Invalid_argument "log_prior: bad params")

let gaussian_log_likelihood params = 
  match params with 
    | [|mu; sigma|] -> 
      let sum = ref 0.0 in 
        for i = 0 to Array.length data - 1 do 
          sum := !sum +. Stats.log_gaussian mu sigma data.(i)
        done;
        !sum +. 0.0
    | _ -> raise (Failure "gaussian_log_likelihood: bad params")

let cauchy_log_likelihood params = 
  match params with 
    | [|mu; sigma|] -> 
      let sum = ref 0.0 in 
        for i = 0 to Array.length data - 1 do 
          sum := !sum +. Stats.log_cauchy mu sigma data.(i)
        done;
        !sum +. 0.0
    | _ -> raise (Failure "cauchy_log_likelihood: bad params")

let jump params = 
  match params with 
    | [| mu; sigma |] -> 
      [|Mcmc.uniform_wrapping mumin mumax 0.1 mu;
        Mcmc.uniform_wrapping sigmamin sigmamax 0.1 sigma|]
    | _ -> raise (Invalid_argument "jump: bad params")

let log_jump_prob _ _ = 0.0

let _ = Mcmc.reset_counters ()
let gsamples = 
  Mcmc.mcmc_array ~nbin:10000 ~nskip:100 10000 gaussian_log_likelihood log_prior jump log_jump_prob [|0.0; 1.0|]
let (gaccept, greject) = Mcmc.get_counters ()

let _ = Mcmc.reset_counters ()
let csamples = 
  Mcmc.mcmc_array ~nbin:10000 ~nskip:100 10000 cauchy_log_likelihood log_prior jump log_jump_prob [|0.0; 1.0|]
let (caccept, creject) = Mcmc.get_counters ()

let _ = 
  let out = open_out "gaussian.dat" in
    Read_write.write (fun x -> x) out gsamples;
    close_out out;
    let out = open_out "cauchy.dat" in
      Read_write.write (fun x -> x) out csamples;
      close_out out

let _ = 
  let gev = Ev.evidence_harmonic_mean gsamples and 
      cev = Ev.evidence_harmonic_mean csamples in
    Printf.eprintf "Gaussian accept %d, reject %d\n" gaccept greject;
    Printf.eprintf "Cauchy accept %d, reject %d\n" caccept creject;
    Printf.eprintf "Gaussian harmonic mean ev = %g\n" gev;
    Printf.eprintf "Cauchy harmonic mean ev = %g\n" cev;
    Printf.eprintf "Odds Ratio = %g\n" (gev/.cev);
    Printf.eprintf "Bayes ratio = %g\n%!" (gev*.pgaussian/.(cev*.pcauchy))

let ginterp = 
  Interp.make (Array.map (fun {Mcmc.value = x} -> x) gsamples) [|mumin; sigmamin|] [|mumax; sigmamax|]

let cinterp = 
  Interp.make (Array.map (fun {Mcmc.value = x} -> x) csamples) [|mumin; sigmamin|] [|mumax; sigmamax|]

let do_it n = 
  let gjump params = Interp.draw_high_level n ginterp and 
      cjump params = Interp.draw_high_level n cinterp in
  let gprob _ params = log (Interp.jump_prob_high_level n ginterp () params) and 
      cprob _ params = log (Interp.jump_prob_high_level n cinterp () params) in 
  let _ = Mcmc.reset_counters () in
  let full_rj = 
    Mcmc.rjmcmc_array ~nbin:10000 ~nskip:100 10000
      (gaussian_log_likelihood, cauchy_log_likelihood)
      (log_prior, log_prior)
      (jump, jump)
      (log_jump_prob, log_jump_prob)
      (gjump, cjump)
      (gprob, cprob)
      (pgaussian, pcauchy)
      ([|0.0; 1.0|], [|0.0; 1.0|]) in
  let (naccept, nreject) = Mcmc.get_counters () and 
      er = Mcmc.rjmcmc_evidence_ratio full_rj in 
    Printf.printf "%d %g %d %d\n%!" n er naccept nreject

let _ = 
  let rec loop i = 
    if i >= 20000 then 
      ()
    else begin
      do_it i;
      loop (2*i)
    end in 
    loop 1
