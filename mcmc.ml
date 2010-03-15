type like_prior = {
  log_likelihood : float;
  log_prior : float;
}

type 'a mcmc_sample = {
  value : 'a;
  like_prior : like_prior
}

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
      if log (Random.float 1.0) < log_accept_prob then 
        {value = proposed;
         like_prior = {log_likelihood = proposed_like; log_prior = proposed_prior}}
      else
        x

let mcmc_array n log_likelihood log_prior jump_proposal log_jump_prob start = 
  let samples = 
    Array.make n {value = start;
                  like_prior = {log_likelihood = log_likelihood start;
                                log_prior = log_prior start}} in 
  let sample = make_mcmc_sampler log_likelihood log_prior jump_proposal log_jump_prob in 
    for i = 1 to n - 1 do 
      samples.(i) <- sample samples.(i-1)
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
  let jump_proposal = function 
    | A(a) -> 
        if Random.float 1.0 < pa then 
          A(jpa a)
        else
          B(jintob ())
    | B(b) -> 
        if Random.float 1.0 < pb then 
          B(jpb b)
        else
          A(jintoa ()) and 
      log_jump_prob x y = 
    match x,y with 
      | A(a), A(a') -> 
          ljpa a a'
      | A(a), B(b) -> 
          (log pb) +. ljpintob b
      | B(b), A(a) -> 
          (log pa) +. ljpintoa a
      | B(b), B(b') -> 
          ljpb b b' and 
      log_like = function 
        | A(a) -> lla a
        | B(b) -> llb b and 
      log_prior = function 
        | A(a) -> (log pa) +. lpa a
        | B(b) -> (log pb) +. lpb b in 
    make_mcmc_sampler log_like log_prior jump_proposal log_jump_prob

let rjmcmc_array n (lla, llb) (lpa, lpb) (jpa, jpb) (ljpa, ljpb) (jintoa, jintob) 
    (ljpintoa, ljpintob) (pa,pb) (a,b) = 
  let is_a = Random.float 1.0 < pa in 
  let next_state = 
    make_rjmcmc_sampler (lla,llb) (lpa,lpb) (jpa,jpb) (ljpa,ljpb) (jintoa,jintob) (ljpintoa, ljpintob) (pa,pb) in 
  let value = if is_a then A(a) else B(b) and 
      log_like = if is_a then lla a else llb b and 
      log_prior = if is_a then lpa a +. (log pa) else lpb b +. (log pb) in 
  let states = Array.make n 
    {value = value; like_prior = {log_likelihood = log_like; log_prior = log_prior}} in 
    for i = 1 to n - 1 do 
      let last = states.(i-1) in 
        states.(i) <- next_state last
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

let log_sum_logs la lb = 
  if la > lb then 
    let lr = lb -. la in 
      la +. (log (1.0 +. (exp lr)))
  else
    let lr = la -. lb in 
      lb +. (log (1.0 +. (exp lr)))

let make_admixture_mcmc_sampler (lla, llb) (lpa, lpb) (jpa, jpb) (ljpa, ljpb) (pa, pb) (va, vb) = 
  let log_pa = log pa and 
      log_pb = log pb and
      log_va = log va and 
      log_vb = log vb in 
  let log_likelihood (lam,a,b) = 
    (* Likelihood includes priors, too, since not multiplicative. *)
    log_sum_logs
      ((log lam) +. (lla a) +. (lpa a) +. log_pa -. log_vb)
      ((log (1.0 -. lam)) +. (llb b) +. (lpb b) +. log_pb -. log_va) and 
      log_prior _ = 0.0 and 
      propose (lam,a,b) = 
    (Random.float 1.0, jpa a, jpb b) and 
      log_jump_prob (_,a,b) (_, a', b') = 
    (ljpa a a') +. (ljpb b b') in 
    make_mcmc_sampler log_likelihood log_prior propose log_jump_prob

let admixture_mcmc_array n (lla, llb) (lpa, lpb) (jpa, jpb) (ljpa, ljpb) (pa, pb) (va, vb) (a, b) = 
  let lam = Random.float 1.0 in 
  let start = 
    {value = (lam, a, b);
     like_prior = 
        {log_likelihood = 
            log_sum_logs 
              ((log lam) +. (lla a) +. (lpa a) +. (log pa))
              ((log (1.0 -. lam)) +. (llb b) +. (lpb b) +. (log pb));
         log_prior = 0.0}} in 
  let next = make_admixture_mcmc_sampler (lla,llb) (lpa,lpb) (jpa, jpb) (ljpa, ljpb) (pa,pb) (va,vb) in 
  let samps = Array.make n start in 
    for i = 1 to n - 1 do 
      let last = samps.(i-1) in 
        samps.(i) <- next last
    done;
    samps

let admixture_evidence_ratio_mcmc_array n lambdas = 
  let get_lam x = let (lam,_,_) = x.value in lam in 
  let mean_lam = Stats.meanf get_lam lambdas in
  let std_lam = Stats.stdf ~mean:mean_lam get_lam lambdas in 
  let r_mean = 1.0 /. (2.0 -. 3.0*.mean_lam) -. 1.0 in 
  let r_mean = if r_mean < 0.0 then 0.0 else r_mean in 
  let denom_root = 2.0 -. 3.0 *. mean_lam in 
  let dr_mean = 3.0 *. std_lam /. denom_root /. denom_root in 
  let log_prior r = 
    if r < 1.0 then 
      log 0.5
    else
      (log 0.5) -. (2.0)*.(log r) and 
      log_likelihood r = 
    let ll = ref 0.0 in 
      for i = 0 to Array.length lambdas - 1 do 
        let (lam, _, _) = lambdas.(i).value in 
          ll := !ll +. (log_sum_logs 
                          ((log lam) +. (log r)) 
                          (log (1.0 -. lam)))
      done;
      !ll +. 0.0 and 
      propose r = 
    let dr = (Random.float 1.0 -. 0.5)*.dr_mean in 
      r +. dr and 
      log_jump_prob _ _ = 0.0 in 
    mcmc_array n log_likelihood log_prior propose log_jump_prob r_mean
      
