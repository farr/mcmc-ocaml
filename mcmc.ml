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

let log_sum_logs la lb = 
  if la = neg_infinity && lb = neg_infinity then 
    neg_infinity
  else if la > lb then 
    let lr = lb -. la in 
      la +. (log (1.0 +. (exp lr)))
  else
    let lr = la -. lb in 
      lb +. (log (1.0 +. (exp lr)))

let make_admixture_mcmc_sampler (lla, llb) (lpa, lpb) (jpa, jpb) (ljpa, ljpb) (pa, pb) (va, vb) = 
  assert(abs_float (pa +. pb -. 1.0) < sqrt epsilon_float);
  let log_pa = log pa and 
      log_pb = log pb and
      log_va = log va and 
      log_vb = log vb in 
  let log_likelihood (lam,a,b) = 
    (* Likelihood includes priors, too, since not multiplicative. *)
    let lpa = ((log lam) +. (lla a) +. (lpa a) +. log_pa -. log_vb) and 
        lpb = ((log (1.0 -. lam)) +. (llb b) +. (lpb b) +. log_pb -. log_va) in
    let lp = log_sum_logs lpa lpb in 
      lp and
      log_prior _ = 0.0 and 
      propose (lam,a,b) = 
    (Random.float 1.0, jpa a, jpb b) and 
      log_jump_prob (_,a,b) (_, a', b') = 
    (ljpa a a') +. (ljpb b b') in 
    make_mcmc_sampler log_likelihood log_prior propose log_jump_prob

let admixture_mcmc_array n (lla, llb) (lpa, lpb) (jpa, jpb) (ljpa, ljpb) (pa, pb) (va, vb) (a, b) = 
  assert(abs_float (pa +. pb -. 1.0) < sqrt epsilon_float);
  let lam = Random.float 1.0 in 
  let start = 
    {value = (lam, a, b);
     like_prior = 
        {log_likelihood = 
            log_sum_logs 
              ((log lam) +. (lla a) +. (lpa a) +. (log pa) -. (log vb))
              ((log (1.0 -. lam)) +. (llb b) +. (lpb b) +. (log pb) -. (log va));
         log_prior = 0.0}} in 
  let next = make_admixture_mcmc_sampler (lla,llb) (lpa,lpb) (jpa, jpb) (ljpa, ljpb) (pa,pb) (va,vb) in 
  let samps = Array.make n start in 
    for i = 1 to n - 1 do 
      let last = samps.(i-1) in 
        samps.(i) <- next last
    done;
    samps

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

let rec bisect_root epsf epsabs epsrel f x0 x1 = 
  let rec bisect_loop fx0 fx1 x0 x1 = 
    let x = 0.5*.(x0 +. x1) and
        dx = abs_float (x1 -. x0) and 
        xmag = 0.5*.(abs_float x1 +. (abs_float x0)) in
      if dx <= epsabs +. epsrel*.xmag then 
        x
      else
        let fx = f x in 
          if abs_float fx < epsf then 
            x
          else if fx0*.fx < 0.0 then 
            bisect_loop fx0 fx x0 x
          else if fx*.fx1 < 0.0 then 
            bisect_loop fx fx1 x x1
          else
            raise (Failure "bisect_root (loop): lost bounds on root") in 
  let fx0 = f x0 and 
      fx1 = f x1 in 
    if fx0 *. fx1 > 0.0 then 
      raise (Failure "bisect_root: root not bounded by initial brackets")
    else if fx0 = 0.0 then 
      x0
    else if fx1 = 0.0 then 
      x1
    else
      bisect_loop fx0 fx1 x0 x1

let find_increasing_bracket f x0 x1 = 
  let rec inc_bracket_loop fx0 fx1 x1 = 
    if fx0 *. fx1 <= 0.0 then 
      x1
    else
      let x1 = 2.0*.x1 in 
      inc_bracket_loop fx0 (f x1) x1 in
    inc_bracket_loop (f x0) (f x1) x1

let dlog_like samples r = 
  let sum = ref 0.0 and 
      nint = Array.length samples in 
    for i = 0 to nint - 1 do 
      let {value = (lam,_,_)} = samples.(i) in 
      sum := !sum +. lam /. (r *. lam +. 1.0 -. lam)
    done;
    let n = float_of_int nint in 
      1.0 /. (r +. 1.0) -. !sum /. n

let mean_lambda_ratio samples = 
  let lam = Stats.meanf (fun x -> let {value = (lam,_,_)} = x in lam) samples in 
    if lam > (2.0/.3.0) then 
      infinity
    else if lam < 1.0 /. 3.0 then 
      0.0
    else
      1.0 /. (2.0 -. 3.0*.lam) -. 1.0

let max_like_admixture_ratio samples = 
  let epsabs = sqrt (epsilon_float) and 
      epsrel = sqrt (epsilon_float) and 
      epsf = 0.0 in
  let f r = dlog_like samples r in 
  let r1 = find_increasing_bracket f 0.0 1.0 in 
    bisect_root epsf epsabs epsrel f 0.0 r1

let dlog_prior_uniform r = 
  if r <= 1.0 then 0.0 else -2.0 /. r

let ddlog_prior_uniform r = 
  if r <= 1.0 then 0.0 else 2.0 /. (r *. r)

let max_posterior_admixture_ratio 
    ?(dlog_prior = dlog_prior_uniform)
    ?(ddlog_prior = ddlog_prior_uniform)
    samples = 
  let epsabs = sqrt epsilon_float and 
      epsrel = sqrt epsilon_float and
      epsf = sqrt epsilon_float in
  let f r = dlog_like samples r +. dlog_prior r in 
  let r1 = find_increasing_bracket f 0.0 1.0 in 
    bisect_root epsf epsabs epsrel f 0.0 r1

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

let should_accept_swap_with_higher_beta beta {like_prior = {log_likelihood = ll}} beta' {like_prior = {log_likelihood = oll}} = 
  let log_paccept = (beta'/.beta -. 1.0)*.ll +. (beta/.beta' -. 1.0)*.oll in
    (log (Random.float 1.0)) < log_paccept

let adjust_state current_beta other_beta ({like_prior = {log_likelihood = ll; log_prior = lp}} as other_state) = 
  {other_state with like_prior = {log_likelihood = ll*.(current_beta /. other_beta); log_prior = lp}}

let pt_beta () = 
  let rank = Mpi.comm_rank Mpi.comm_world and 
      np = Mpi.comm_size Mpi.comm_world in 
    (float_of_int (rank+1))/.(float_of_int np)

let pt_dbeta () = 
  let np = Mpi.comm_size Mpi.comm_world in 
    1.0 /. (float_of_int np)

let nswap = ref 0

let reset_nswap () = nswap := 0
let get_nswap () = !nswap

let maybe_exchange_temps beta ({like_prior = {log_likelihood = ll}} as state) = 
  let dbeta = pt_dbeta () in
  let rank = Mpi.comm_rank Mpi.comm_world and
      np = Mpi.comm_size Mpi.comm_world in
  let other = ref state and 
      state = ref state in
    (* 0 <--> 1, 2 <--> 3, ... *)
    if rank = np - 1 && np mod 2 = 1 then begin
      () (* Do Nothing. *)
    end else begin
      (* Evens send up first. *)
      if rank mod 2 = 0 then begin
        Mpi.send !state (rank+1) 0 Mpi.comm_world
      end else begin
        other := Mpi.receive (rank-1) 0 Mpi.comm_world
      end;
      (* Odds send down next. *)
      if rank mod 2 = 1 then begin
        Mpi.send !state (rank-1) 0 Mpi.comm_world
      end else begin
        other := Mpi.receive (rank+1) 0 Mpi.comm_world
      end;
      let should_swap = ref true in 
        (* Evens determine whether there should be a swap. *)
        if rank mod 2 = 0 then begin
          let other_beta = beta +. dbeta in 
            should_swap := should_accept_swap_with_higher_beta beta !state other_beta !other;
        end;
        (* Now evens tell odds about swap. *)
        if rank mod 2 = 0 then 
          Mpi.send_int (if !should_swap then 1 else 0) (rank+1) 0 Mpi.comm_world
        else
          should_swap := Mpi.receive_int (rank-1) 0 Mpi.comm_world = 1;
        if !should_swap then begin
          incr nswap;
          let other_beta = if rank mod 2 = 0 then beta +. dbeta else beta -. dbeta in
          state := adjust_state beta other_beta !other
        end
    end;
    (* 1 <--> 2, 3 <--> 4, ... *)
    if rank = 0 || (np mod 2 = 0 && rank = np - 1) then 
      () (* Do nothing. *)
    else begin
      (* Odds send up first. *)
      if rank mod 2 = 1 then 
        Mpi.send !state (rank+1) 0 Mpi.comm_world
      else
        other := Mpi.receive (rank-1) 0 Mpi.comm_world;
      (* Evens send down next. *)
      if rank mod 2 = 0 then 
        Mpi.send !state (rank-1) 0 Mpi.comm_world
      else
        other := Mpi.receive (rank+1) 0 Mpi.comm_world;
      let should_swap = ref true in 
        (* Odds determine whether there should be a swap *)
        if rank mod 2 = 1 then begin
          let other_beta = beta +. dbeta in 
            should_swap := should_accept_swap_with_higher_beta beta !state other_beta !other
        end;
        (* Now odds tell evens about swap. *)
        if rank mod 2 = 1 then 
          Mpi.send_int (if !should_swap then 1 else 0) (rank+1) 0 Mpi.comm_world
        else
          should_swap := Mpi.receive_int (rank-1) 0 Mpi.comm_world = 1;
        if !should_swap then 
          incr nswap;
          let other_beta = if rank mod 2 = 1 then beta +. dbeta else beta -. dbeta in 
            state := adjust_state beta other_beta !other
    end;
    !state

let make_pt_mcmc_sampler nswap log_like log_prior propose log_jp = 
  let count = ref 1 in 
  let beta = pt_beta () in
    fun ({value = v; like_prior = {log_likelihood = ll; log_prior = lp}} as state) -> 
      if !count mod nswap = 0 then begin
        incr count;
        maybe_exchange_temps beta state
      end else
        let new_v = propose v in 
        let new_ll = log_like new_v and 
            new_lp = log_prior new_v in
        let new_ll = beta *. new_ll in 
        let log_forward = log_jp v new_v and 
            log_backward = log_jp new_v v in 
        let logacceptp = new_ll +. new_lp -. ll -. lp +. log_backward -. log_forward in 
          incr count;
          if (log (Random.float 1.0)) < logacceptp then begin
            incr naccept;
            {value = new_v; like_prior = {log_likelihood = new_ll; log_prior = new_lp}}
          end else begin
            incr nreject;
            state
          end

let pt_mcmc_array ?(nskip = 1) n nswap log_like log_prior propose log_jp start = 
  let beta = pt_beta () in 
  let state = {value = start;
               like_prior = {log_likelihood = beta *. (log_like start);
                             log_prior = log_prior start}} in 
  let states = Array.make n state in 
  let current_state = ref state in 
  let next_state = make_pt_mcmc_sampler nswap log_like log_prior propose log_jp in
    for i = 1 to (n-1)*nskip do 
      current_state := next_state !current_state;
      if i mod nskip = 0 then 
        states.(i/nskip) <- !current_state
    done;
    states

let expected_log_like samples = 
  let n = Array.length samples and
      beta = pt_beta () in 
  let ll_sum = ref 0.0 in 
    for i = 0 to n - 1 do 
      let {like_prior = {log_likelihood = ll}} = samples.(i) in
      ll_sum := !ll_sum +. ll/.beta
    done;
    !ll_sum /. (float_of_int n)

let integrate_lls lls = 
  let db = pt_dbeta () in 
  let sum = ref 0.0 in 
  let n = Array.length lls in 
    for i = 0 to n - 1 do 
      sum := !sum +. db*.lls.(i)
    done;
    !sum +. 0.0

let thermodynamic_integrate samples = 
  let ll = expected_log_like samples in 
  let lls = if Mpi.comm_rank Mpi.comm_world = 0 then Array.make (Mpi.comm_size Mpi.comm_world) 0.0 else [| |] in 
    Mpi.gather_float ll lls 0 Mpi.comm_world;
  let ti_int = 
    if Mpi.comm_rank Mpi.comm_world = 0 then 
      integrate_lls lls
    else
      0.0 in 
    Mpi.broadcast_float ti_int 0 Mpi.comm_world
