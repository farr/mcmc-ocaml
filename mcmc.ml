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
  assert(pa +. pb -. 1.0 < sqrt epsilon_float);
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

let rjmcmc_evidence_ratio samples = 
  let (na,nb) = rjmcmc_model_counts samples in 
    (float_of_int na)/.(float_of_int nb)

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
  let delta_x = mod_float delta_x (xmax -. xmin) in 
  let xnew = x +. delta_x in 
    if xnew > xmax then 
      xmin +. (xnew -. xmax)
    else if xnew < xmin then 
      xmax -. (xmin -. xnew)
    else
      xnew

type ('a, 'b) memory = 
  | Left of 'a * 'b
  | Right of 'a * 'b

let make_memory_rjmcmc_sampler 
    (lla, llb) (lpa, lpb) (jpa, jpb) (ljpa, ljpb) (pa, pb) = 
  if abs_float (pa +. pb -. 1.0) > sqrt epsilon_float then
    raise (Invalid_argument "make_memory_rjmcmc_sampler: sum of model priors not 1.0");
  let log_pa = log pa and 
      log_pb = log pb in
  let log_like = function 
    | Left(a,b) -> 
        lla a
    | Right(a,b) -> 
        llb b in
  let log_prior = function 
    | Left(a,b) -> 
        lpa a +. log_pa
    | Right(a,b) ->
        lpb b +. log_pb in 
  let jump_proposal = function 
    | Left(a,b) -> 
        if Random.float 1.0 < pb then 
          (* Jump to a. *)
          Right(jpa a, jpb b)
        else
          Left(jpa a, b)
    | Right(a,b) -> 
        if Random.float 1.0 < pa then 
          Left(jpa a, jpb b)
        else
          Right(a, jpb b) in
  let log_jump_prob x y = 
    match x,y with 
      | Left(ax, _), Left(ay, _) -> ljpa ax ay +. log_pa
      | Left(ax, bx), Right(ay, by) -> 
          log_pb +. ljpa ax ay +. ljpb bx by
      | Right(ax, bx), Left(ay, by) -> 
          log_pa +. ljpa ax ay +. ljpb bx by
      | Right(_,bx), Right(_, by) -> 
          log_pb +. ljpb bx by in
    make_mcmc_sampler log_like log_prior jump_proposal log_jump_prob

let memory_rjmcmc_array ?(nskip = 1) n (lla, llb) (lpa, lpb) (jpa, jpb) (ljpa, ljpb) (pa, pb) (a, b) =
  let sampler = make_memory_rjmcmc_sampler (lla, llb) (lpa, lpb) (jpa, jpb) (ljpa, ljpb) (pa, pb) in 
  let state = if Random.float 1.0 < pa then 
    {value = Left(a,b);
     like_prior = {log_likelihood = lla a; log_prior = lpa a +. (log pa)}}
  else
    {value = Right(a,b);
     like_prior = {log_likelihood = llb b; log_prior = lpb b +. (log pb)}} in 
  let current_state = ref state in 
  let samples = Array.make n state in 
    for i = 1 to (n - 1)*nskip do 
      current_state := sampler !current_state;
      if i mod nskip = 0 then 
        samples.(i/nskip) <- !current_state
    done;
    samples

let memory_rjmcmc_model_counts samples = 
  let nl = ref 0 and 
      nr = ref 0 in 
    for i = 0 to Array.length samples - 1 do 
      match samples.(i) with 
        | {value = Left(_,_)} -> incr nl
        | {value = Right(_,_)} -> incr nr
    done;
    (!nl, !nr)

let memory_evidence_ratio samples = 
  let (nl, nr) = memory_rjmcmc_model_counts samples in 
    (float_of_int nl)/.(float_of_int nr)

let split_memory_array (pa,pb) samples = 
  if abs_float (pa +. pb -. 1.0) > sqrt epsilon_float then 
    raise (Invalid_argument "split_memory_array: priors do not sum to 1.0");
  let log_pa = log pa and 
      log_pb = log pb in 
  let left = ref [] and 
      right = ref [] in 
    for i = Array.length samples - 1 downto 0 do 
      match samples.(i) with 
        | {value = Left(a,_);
           like_prior = {log_likelihood = ll; log_prior = lp}} -> 
            left := {value = a; like_prior = {log_likelihood = ll; log_prior = lp -. log_pa}} :: !left
        | {value = Right(_,b);
           like_prior = {log_likelihood = ll; log_prior = lp}} -> 
            right := {value = b; like_prior = {log_likelihood = ll; log_prior = lp -. log_pb}} :: !right
    done;
    (Array.of_list !left,
     Array.of_list !right)
