(*  mcmc.ml: MCMC samplers, utility functions.
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

type like_prior = {
  log_likelihood : float;
  log_prior : float;
}

type 'a mcmc_sample = {
  value : 'a;
  like_prior : like_prior
}

let naccept = ref 0
let nreject = ref 0

let reset_counters () = 
  naccept := 0;
  nreject := 0

let get_counters () = 
  (!naccept, !nreject)

let make_mcmc_sampler log_likelihood log_prior jump_proposal log_jump_prob = 
  fun x -> 
    let start = x.value and 
        start_log_post = x.like_prior.log_likelihood +. x.like_prior.log_prior in 
    let proposed = jump_proposal start in 
    let proposed_like = log_likelihood proposed and 
        proposed_prior = log_prior proposed in 
    let proposed_log_posterior = proposed_like +. proposed_prior in 
    let log_forward_jump = log_jump_prob start proposed and 
        log_backward_jump = log_jump_prob proposed start in 
    let log_accept_prob = 
      proposed_log_posterior -. start_log_post +. log_backward_jump -. log_forward_jump in 
      if log (Random.float 1.0) < log_accept_prob then begin
        incr naccept;
        {value = proposed;
         like_prior = {log_likelihood = proposed_like; log_prior = proposed_prior}}
      end else begin
        incr nreject;
        x
      end

let mcmc_array ?(nbin = 0) ?(nskip = 1) n log_likelihood log_prior jump_proposal log_jump_prob start = 
  let current_sample = ref {value = start;
                  like_prior = {log_likelihood = log_likelihood start;
                                log_prior = log_prior start}} in
  let sample = make_mcmc_sampler log_likelihood log_prior jump_proposal log_jump_prob in 
    for i = 0 to nbin - 1 do 
      current_sample := sample !current_sample
    done;
    let samples = Array.make n !current_sample in 
      for i = 1 to (n - 1) * nskip do 
        current_sample := sample !current_sample;
        if i mod nskip = 0 then 
          samples.(i / nskip) <- !current_sample
      done;
      samples

let remove_repeat_samples eql samps = 
  let removed = ref [] in 
    for i = Array.length samps - 1 downto 1 do 
      if not (eql samps.(i).value samps.(i-1).value) then 
        removed := samps.(i) :: !removed
    done;
    removed := samps.(0) :: !removed;
    Array.of_list !removed

type ('a, 'b) rjmcmc_value = 
  | A of 'a
  | B of 'b

type ('a, 'b) rjmcmc_sample = ('a, 'b) rjmcmc_value mcmc_sample

let make_rjmcmc_sampler (lla, llb) (lpa, lpb) (jpa, jpb) (ljpa, ljpb) (jintoa, jintob) (ljpintoa, ljpintob) (pa,pb) = 
  assert(pa +. pb -. 1.0 < sqrt epsilon_float);
  let log_pa = log pa and log_pb = log pb in
  let jump_proposal = function 
    | A(a) -> 
        if Random.float 1.0 < pa then
          A(jpa a)
        else
          B(jintob a)
    | B(b) -> 
        if Random.float 1.0 < pb then 
          B(jpb b)
        else
          A(jintoa b) and 
      log_jump_prob x y = 
    match x,y with 
      | A(a), A(a') -> 
        log_pa +. ljpa a a'
      | A(a), B(b) -> 
        log_pb +. ljpintob a b
      | B(b), A(a) -> 
        log_pa +. ljpintoa b a
      | B(b), B(b') -> 
        log_pb +. ljpb b b' and 
      log_like = function 
        | A(a) -> lla a
        | B(b) -> llb b and 
      log_prior = function 
        | A(a) -> log_pa +. lpa a
        | B(b) -> log_pb +. lpb b in 
    make_mcmc_sampler log_like log_prior jump_proposal log_jump_prob

let rjmcmc_array ?(nbin = 0) ?(nskip = 1) n (lla, llb) (lpa, lpb) (jpa, jpb) (ljpa, ljpb) (jintoa, jintob) 
    (ljpintoa, ljpintob) (pa,pb) (a,b) = 
  let is_a = Random.float 1.0 < 0.5 in 
  let next_state = 
    make_rjmcmc_sampler (lla,llb) (lpa,lpb) (jpa,jpb) (ljpa,ljpb) (jintoa,jintob) (ljpintoa, ljpintob) (pa,pb) in 
  let value = if is_a then A(a) else B(b) and 
      log_like = if is_a then lla a else llb b and 
      log_prior = if is_a then lpa a +. (log pa) else lpb b +. (log pb) in 
  let current_state = ref {value = value; like_prior = {log_likelihood = log_like; log_prior = log_prior}} in
    for i = 0 to nbin - 1 do 
      current_state := next_state !current_state
    done;
    let states = Array.make n !current_state in 
      for i = 1 to (n-1) * nskip do 
        current_state := next_state !current_state;
        if i mod nskip = 0 then 
          states.(i/nskip) <- !current_state
      done;
      states

let rjmcmc_model_counts data = 
  let na = ref 0 and 
      nb = ref 0 in 
    for i = 0 to Array.length data - 1 do 
      match data.(i).value with 
        | A(_) -> incr na
        | B(_) -> incr nb
    done;
    (!na, !nb)

let rjmcmc_evidence_ratio samples = 
  let (na,nb) = rjmcmc_model_counts samples in 
    (float_of_int na) /. (float_of_int nb)

let log_sum_logs la lb = 
  if la = neg_infinity && lb = neg_infinity then 
    neg_infinity
  else if la > lb then 
    let lr = lb -. la in 
      la +. (log (1.0 +. (exp lr)))
  else
    let lr = la -. lb in 
      lb +. (log (1.0 +. (exp lr)))

let combine_jump_proposals props = 
  let ptot = List.fold_left (fun ptot (p,_,_) -> ptot +. p) 0.0 props in 
  let props = List.map (fun (p, jp, ljp) -> (p /. ptot, jp, ljp)) props in 
  let propose x = 
    let jump = let rec loop prob props = 
      match props with 
        | (p, jp, _) :: props -> 
            if prob < p then jp else loop (prob -. p) props
        | _ -> raise (Failure "combine_jump_proposals: internal error: no jump proposal to select") in 
      loop (Random.float 1.0) props in
      jump x in 
  let log_jp x y = 
    let log_jp = 
      List.fold_left 
        (fun log_jump (p,_,ljp) -> 
           let log_local = (log p) +. (ljp x y) in 
             log_sum_logs log_jump log_local)
        neg_infinity
        props in 
      log_jp in
    (propose, log_jp)

let uniform_wrapping xmin xmax dx x = 
  let delta_x = (Random.float 1.0 -. 0.5)*.dx in 
  let rec loop new_x = 
    if new_x < xmin then 
      loop (xmin +. (xmin -. new_x))
    else if new_x >= xmax then 
      loop (xmax -. (new_x -. xmax))
    else
      new_x in 
    loop (x +. delta_x)

let differential_evolution_proposal ?(mode_hopping_frac = 0.0) to_float from_float samples = 
  let pick_samples () = 
    let n = Array.length samples in 
    let i = Random.int n in 
    let j = let rec loop () = let jtry = Random.int n in if jtry = i then loop () else jtry in loop () in 
      (to_float samples.(i).value, to_float samples.(j).value) in 
    fun current -> 
      let z = to_float current in 
      let (x, y) = pick_samples () in 
      let ndim = Array.length z in 
      let d = 
        if mode_hopping_frac <> 0.0 && Random.float 1.0 < mode_hopping_frac then 
          1.0 
        else 
          Stats.draw_uniform 0.0 2.0 in          
      let z' = Array.make ndim 0.0 in 
        for i = 0 to ndim - 1 do 
          z'.(i) <- z.(i) +. d*.(y.(i) -. x.(i))
        done;
        from_float z'
