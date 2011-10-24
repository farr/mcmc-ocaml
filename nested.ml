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
open Stats

type 'a nested_output = float * float * 'a Mcmc.mcmc_sample array * float array

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

let remaining_integral_negligable log_accumulated_integral log_remaining_volume livepts epsrel = 
  let n = Array.length livepts in 
  let log_live_estimate = log_remaining_volume +. (livepts.(n-1).like_prior.log_likelihood) in 
    log_live_estimate -. (log_sum_logs log_accumulated_integral log_live_estimate) <= log epsrel

let draw_new_live_point to_array from_array nmcmc mode_hop livepts log_likelihood log_prior = 
  let nlive = Array.length livepts in 
  let log_l_threshold = livepts.(0).like_prior.log_likelihood in 
  let jump_proposal = differential_evolution_proposal ~mode_hopping_frac:mode_hop to_array from_array livepts in 
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
      if not (new_live.like_prior.log_likelihood >= livepts.(0).like_prior.log_likelihood) then
        raise (Failure (Printf.sprintf "Error in draw_new_live_point: new log(L) of %g larger than old log(L) of %g\n%!"
          new_live.like_prior.log_likelihood livepts.(0).like_prior.log_likelihood))
      else        
        new_live

(* Each point is associated with a particular remaining volume that
   runs from 1 down to 0.  All but the final nlive points are
   associated with a reduction in volume by a constant factor
   (1-1/nlive); the final nlive points divide the remaining volume
   equally, since they were live when the sampling stopped. *)
let evidence_error_and_weights nlive all_pts = 
  let vol_fraction = 1.0 /. (float_of_int nlive) in 
  let log_vol_fraction = log vol_fraction in 
  let log_reduction_frac = log1p (~-.vol_fraction) in
  let log_half = -0.69314718055994530942 in
  let n = Array.length all_pts in 
  let wts = Array.make n neg_infinity in
  let low = ref neg_infinity and 
      high = ref neg_infinity in
  let ilive = n - nlive in (* ilive is the index where the points that were live at the end of the sampling start *)
    for i = 0 to ilive - 1 do 
      (* Each point here reduces the available volume by a constant
         factor. *)
      let logli = all_pts.(i).like_prior.log_likelihood and 
          logli1 = all_pts.(i+1).like_prior.log_likelihood in
      let log_dv = log_vol_fraction +. ((float_of_int i) *. log_reduction_frac) in
      let log_dlow = log_dv +. logli and
          log_dhigh = log_dv +. logli1 in 
        low := log_sum_logs !low log_dlow;
        high := log_sum_logs !high log_dhigh;
        wts.(i) <- log_sum_logs wts.(i) (log_half +. log_dlow);
        wts.(i+1) <- log_sum_logs wts.(i+1) (log_half +. log_dhigh)
    done;
    let log_dv = log_vol_fraction +. ((float_of_int (ilive-1)) *. log_reduction_frac) in 
      for i = ilive to n - 1 do 
        let logli1 = all_pts.(i-1).like_prior.log_likelihood and 
            logli = all_pts.(i).like_prior.log_likelihood in
        let log_dlow = log_dv +. logli1 and 
            log_dhigh = log_dv +. logli in 
          low := log_sum_logs !low log_dlow;
          high := log_sum_logs !high log_dhigh;
          wts.(i-1) <- log_sum_logs wts.(i-1) (log_half +. log_dlow);
          wts.(i) <- log_sum_logs wts.(i) (log_half +. log_dhigh)
      done;
      let log_ev = log_half +. (log_sum_logs !low !high) and 
          log_dev = !high +. log1p (~-.(exp (!low -. !high))) in 
        for i = 0 to n - 1 do 
          wts.(i) <- wts.(i) -. log_ev
        done;
        (log_ev, log_dev, all_pts, wts)

let nested_evidence ?observer ?(epsrel = 0.01) ?(nmcmc = 1000) ?(nlive = 1000) ?(mode_hopping_frac = 0.1) to_float from_float draw_prior log_likelihood log_prior = 
  let observer = match observer with 
    | Some(o) -> o
    | None -> fun _ -> () in
  let livepts = 
    Array.map 
      (fun pt -> {value = pt; like_prior = {log_likelihood = log_likelihood pt; log_prior = log_prior pt}})
      (Array.init nlive (fun _ -> draw_prior ())) in 
  let vol_fraction = 1.0 /. (float_of_int nlive) in 
  let log_volume_reduction_factor = log1p (~-.vol_fraction) in 
    Array.fast_sort (fun p1 p2 -> Pervasives.compare p1.like_prior.log_likelihood p2.like_prior.log_likelihood) livepts;
    let rec nested_loop log_vol_remaining log_integral_estimate retired_pts livepts = 
      let new_pt = draw_new_live_point to_float from_float nmcmc mode_hopping_frac livepts log_likelihood log_prior in 
      let retired_pt = replace_live_point livepts new_pt in 
        ignore(observer retired_pt);
        let new_retired = retired_pt :: retired_pts in 
        let log_retired_l = retired_pt.like_prior.log_likelihood in 
        let log_new_vol = log_vol_remaining +. log_volume_reduction_factor in 
        let log_dv = log_vol_remaining +. vol_fraction in 
        let log_new_estimate = log_sum_logs log_integral_estimate (log_retired_l +. log_dv) in 
          if remaining_integral_negligable log_new_estimate log_new_vol livepts epsrel then 
            evidence_error_and_weights nlive (Array.of_list (List.rev_append new_retired (Array.to_list livepts)))
          else
            nested_loop log_new_vol log_new_estimate new_retired livepts in 
      nested_loop 0.0 neg_infinity [] livepts
        
let log_total_error_estimate log_ev log_dev nlive = 
  let log_rel_error2 = ~-. (log (float_of_int nlive)) in 
    0.5 *. (log_sum_logs (2.0*.log_dev) (log_rel_error2 +. 2.0*.log_ev))

let weight_binary_search_index x running_sums = 
  if x <= running_sums.(0) then 
    0
  else
    let rec bs_loop ilow ihigh = 
      if ihigh - ilow <= 1 then 
        ihigh
      else
        let imid = (ilow + ihigh) / 2 in 
          if x <= running_sums.(imid) then 
            bs_loop ilow imid
          else
            bs_loop imid ihigh in 
      bs_loop 0 (Array.length running_sums - 1)

let posterior_samples n (_, _, all_pts, log_wts) = 
  let npts = Array.length all_pts in 
    assert(Array.length log_wts = npts);
    let summed_weights = Array.make npts (exp log_wts.(0)) in 
      for i = 1 to npts - 1 do 
        summed_weights.(i) <- (exp log_wts.(i)) +. summed_weights.(i-1)
      done;
      Array.init n 
        (fun _ -> 
          let x = Random.float 1.0 in 
          let i = weight_binary_search_index x summed_weights in 
            all_pts.(i))
