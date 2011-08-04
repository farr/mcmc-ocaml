(*  nested.ml: Nested sampling routines.
    Copyright (C) 2011 Will M. Farr <w-farr@northwestern.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. *)

open Mcmc

(* Replaces the lowest-likelihood live point (the lpoints array is
   assumed to be sorted in order of increasing likelihood) with the
   given point.  Returns the now-dead point that has been replaced.
   The algorithm is basically insertion sort. *)
let replace_live_point lpoints newpt = 
  let deadpt = lpoints.(0) in 
    lpoints.(0) <- newpt;
    let n = Array.length lpoints in
    let rec loop i = 
      if i >= n - 1 then 
        deadpt
      else
        let p = lpoints.(i) and 
            pnext = lpoints.(i+1) in 
          if p.like_prior.log_likelihood > pnext.like_prior.log_likelihood then begin
            lpoints.(i) <- pnext;
            lpoints.(i+1) <- p;
            loop (i+1)
          end else
            (* The new point is in the correct position. *)
            deadpt in 
      loop 0

(* Provides a low and high estimate of the Lebesque integral given the
   list of sample points.  Each successive point is assumed to remove
   1/nlive of the remaining available prior mass.  A trapezoidal-rule
   estimate of the integral would be an average of low and high. *)
let estimate_integral nlive pts = 
  let volume_fraction = 1.0 /. (float_of_int nlive) in 
  let rec accum_loop low_accum high_accum pts live_vol = 
    match pts with 
      | [] -> (low_accum, high_accum) 
      | [pt] -> (low_accum, high_accum)
      | p1 :: (p2 :: rest as remaining) -> 
        let l1 = exp p1.like_prior.log_likelihood and 
            l2 = exp p2.like_prior.log_likelihood in 
        let vlarge = live_vol and 
            vsmall = live_vol*.(1.0 -. volume_fraction) in 
        let height = l2 -. l1 in 
          accum_loop (low_accum +. vsmall*.height) (high_accum +. vlarge*.height) remaining vsmall in 
    accum_loop 0.0 0.0 pts 1.0

let remaining_integral_negligable accumulated_integral remaining_volume livepts epsrel = 
  let n = Array.length livepts in 
  let live_estimate = remaining_volume*.(exp livepts.(n-1).like_prior.log_likelihood) in 
    live_estimate /. (accumulated_integral +. live_estimate) <= epsrel

let draw_new_live_point nmcmc mode_hop livepts log_likelihood log_prior = 
  let nlive = Array.length livepts in 
  let log_l_threshold = livepts.(0).like_prior.log_likelihood in 
  let jump_proposal = differential_evolution_proposal ~mode_hopping_frac:mode_hop livepts in 
  let mcmc_logl pt = 
    let ll = log_likelihood pt in 
      if ll >= log_l_threshold then 
        log_prior pt 
      else
        neg_infinity in 
  let mcmc_logp pt = 0.0 in (* Uniform "prior" for mcmc. *)
  let mcmc_logjp _ _ = 0.0 in (* Symmetric jump proposal. *)
  let next_pt = make_mcmc_sampler mcmc_logl mcmc_logp jump_proposal mcmc_logjp in 
  let pt = livepts.(Random.int nlive).value in 
  let current_point = ref {value = pt; like_prior = {log_likelihood = mcmc_logl pt; log_prior = mcmc_logp pt}} in
    for i = 0 to nmcmc - 1 do 
      current_point := next_pt !current_point
    done;
    let new_live = {value = !current_point.value; like_prior = {log_likelihood = log_likelihood !current_point.value;
                                                                log_prior = log_prior !current_point.value}} in 
      assert(new_live.like_prior.log_likelihood >= livepts.(0).like_prior.log_likelihood);
      new_live

(* Each point is associated with a particular remaining volume that
   runs from 1 down to 0.  All but the final nlive points are
   associated with a reduction in volume by a constant factor
   (1-1/nlive); the final nlive points divide the remaining volume
   equally, since they were live when the sampling stopped. *)
let evidence_error_and_weights nlive all_pts = 
  let vol_fraction = 1.0 /. (float_of_int nlive) in 
  let log_reduction_frac = log1p (~-.vol_fraction) in
  let n = Array.length all_pts in 
  let wts = Array.make n 0.0 in
  let low = ref 0.0 and 
      high = ref 0.0 in
  let ilive = n - nlive in (* ilive is the index where the points that were live at the end of the sampling start *)
    for i = 0 to ilive - 1 do 
      (* Each point here reduces the available volume by a constant
         factor. *)
      let li = exp all_pts.(i).like_prior.log_likelihood and 
          li1 = exp all_pts.(i+1).like_prior.log_likelihood in
      let dv = vol_fraction *. (exp ((float_of_int i) *. log_reduction_frac)) in
      let dlow = dv *. li and
          dhigh = dv *. li1 in 
        low := !low +. dlow;
        high := !high +. dhigh;
        wts.(i) <- wts.(i) +. 0.5*.dlow;
        wts.(i+1) <- wts.(i+1) +. 0.5*.dhigh
    done;
    let dv = vol_fraction *. (exp ((float_of_int (ilive-1)) *. log_reduction_frac)) in 
      for i = ilive to n - 1 do 
        let li1 = exp all_pts.(i-1).like_prior.log_likelihood and 
            li = exp all_pts.(i).like_prior.log_likelihood in
        let dlow = dv *. li1 and 
            dhigh = dv *. li in 
          low := !low +. dlow;
          high := !high +. dhigh;
          wts.(i-1) <- wts.(i-1) +. 0.5*.dlow;
          wts.(i) <- wts.(i) +. 0.5*.dhigh
      done;
      let ev = 0.5*.(!low +. !high) and 
          dev = !high -. !low in 
        for i = 0 to n - 1 do 
          wts.(i) <- wts.(i) /. ev
        done;
        (ev, dev, all_pts, wts)

      let nested_evidence ?(epsrel = 0.01) ?(nmcmc = 1000) ?(nlive = 1000) ?(mode_hopping_frac = 0.1) draw_prior log_likelihood log_prior = 
        let livepts = 
          Array.map 
            (fun pt -> {value = pt; like_prior = {log_likelihood = log_likelihood pt; log_prior = log_prior pt}})
            (Array.init nlive (fun _ -> draw_prior ())) in 
        let volume_reduction_factor = 1.0 -. 1.0 /. (float_of_int nlive) in 
          Array.fast_sort (fun p1 p2 -> Pervasives.compare p1.like_prior.log_likelihood p2.like_prior.log_likelihood) livepts;
          let rec nested_loop vol_remaining integral_estimate retired_pts livepts = 
            let new_pt = draw_new_live_point nmcmc mode_hopping_frac livepts log_likelihood log_prior in 
            let retired_pt = replace_live_point livepts new_pt in 
            let new_retired = retired_pt :: retired_pts in 
            let retired_l = exp retired_pt.like_prior.log_likelihood in 
            let new_vol = vol_remaining *. volume_reduction_factor in 
            let dv = vol_remaining -. new_vol in 
            let new_estimate = integral_estimate +. retired_l *. dv in 
              if remaining_integral_negligable new_estimate new_vol livepts epsrel then 
                evidence_error_and_weights nlive (Array.of_list (List.rev_append new_retired (Array.to_list livepts)))
              else
                nested_loop new_vol new_estimate new_retired livepts in 
            nested_loop 1.0 0.0 [] livepts

let total_error_estimate ev dev nlive = 
  let rel_error = 1.0 /. (sqrt (float_of_int nlive)) in 
    sqrt (dev*.dev +. rel_error*.rel_error*.ev*.ev)
