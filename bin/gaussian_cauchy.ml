let ndim = ref 1
let nsamp = ref 10
let samples_fname = ref "samples.dat"
let params_fname = ref "params.dat"
let gauss_mcmc_fname = ref "gaussian.dat"
let cauchy_mcmc_fname = ref "cauchy.dat"
let rjmcmc_fname = ref "rjmcmc.dat"
let admcmc_fname = ref "admixture.dat"
let mu_min = ref (-1.0)
let mu_max = ref 1.0 
let sigma_min = ref 0.1
let sigma_max = ref 0.2
let nmcmc = ref 100000
let interp_reduce = ref 1
let gprior = ref 0.5
let cprior = ref 0.5
let randmcmc = ref false

let opts = 
  [("-ndim", Arg.Set_int ndim,
    Printf.sprintf "number of dimensions in sample (default %d)" !ndim);
   ("-nsamp", Arg.Set_int nsamp,
    Printf.sprintf "number of samples to draw from gaussian (default %d)" !nsamp);
   ("-samples", Arg.Set_string samples_fname,
    Printf.sprintf "filename for samples (default %s)" !samples_fname);
   ("-params", Arg.Set_string params_fname,
    Printf.sprintf "filename for params (default %s)" !params_fname);
   ("-gaussian-mcmc", Arg.Set_string gauss_mcmc_fname,
    Printf.sprintf "filename for gaussian mcmc (default %s)" !gauss_mcmc_fname);
   ("-cauchy-mcmc", Arg.Set_string cauchy_mcmc_fname,
    Printf.sprintf "filename for chaucy mcmc (default %s)" !cauchy_mcmc_fname);
   ("-rjmcmc", Arg.Set_string rjmcmc_fname,
    Printf.sprintf "filename for reverse-jump mcmc (default %s)" !rjmcmc_fname);
   ("-admcmc", Arg.Set_string admcmc_fname,
    Printf.sprintf "filename for admixture mcmc (default %s)" !admcmc_fname);
   ("-mu-min", Arg.Set_float mu_min,
    Printf.sprintf "minimum mean parameter (default %g)" !mu_min);
   ("-mu-max", Arg.Set_float mu_max,
    Printf.sprintf "maximum mean parameter (default %g)" !mu_max);
   ("-sigma-min", Arg.Set_float sigma_min,
    Printf.sprintf "minimum standard deviation parameter (default %g)" !sigma_min);
   ("-sigma-max", Arg.Set_float sigma_max,
    Printf.sprintf "maximum standard deviation parameter (default %g)" !sigma_max);
   ("-nmcmc", Arg.Set_int nmcmc,
    Printf.sprintf "number of mcmc samples (default %d)" !nmcmc);
   ("-seed", Arg.Int (fun s -> Random.init s),
    "RNG seed (default self_init)");
   ("-interpfrac", Arg.Set_int interp_reduce, 
    (Printf.sprintf "reduction factor in number of interpolating samples (default %d)" !interp_reduce));
   ("-gprior", Arg.Float (fun gp -> gprior := gp; cprior := 1.0 -. gp),
    (Printf.sprintf "gaussian model prior (will also set cauchy prior accordingly; default %g)" !gprior));
   ("-cprior", Arg.Float (fun cp -> cprior := cp; gprior := 1.0 -. cp),
    (Printf.sprintf "cauchy model prior (will also set gaussian prior accordingly; default %g)" !cprior));
   ("-randmcmc", Arg.Set randmcmc,
    "re-randomize the RNG after draws, before running any MCMC")]

let randomize () = 
  let input = open_in_bin "/dev/random" in 
    Random.init (input_binary_int input);
    close_in input

let random_between a b = 
  a +. Random.float (b-.a)

let random_array n a b = 
  Array.init n (fun _ -> random_between a b)

let print_array chan arr = 
  let n = Array.length arr in 
    for i = 0 to n - 2 do 
      Printf.fprintf chan "%g " arr.(i)
    done;
    Printf.fprintf chan "%g\n" arr.(n-1)

let make_params () = 
  let mus = random_array !ndim !mu_min !mu_max and 
      sigmas = random_array !ndim !sigma_min !sigma_max in 
  let out = open_out !params_fname in 
    Printf.fprintf out "# first line: means, second line: sigmas\n";
    print_array out mus;
    print_array out sigmas;
    close_out out;
    (mus, sigmas)

let draw_samples mus sigmas = 
  let samps = Array.make !nsamp [|0.0|] in 
    for i = 0 to !nsamp - 1 do 
      let samp = Array.make !ndim 0.0 in 
        for j = 0 to !ndim - 1 do 
          samp.(j) <- Stats.draw_gaussian mus.(j) sigmas.(j)
        done;
        samps.(i) <- samp
    done;
    let out = open_out !samples_fname in 
      Printf.fprintf out "# each line a gaussian sample\n";
      for i = 0 to !nsamp - 1 do 
        print_array out samps.(i)
      done;
      close_out out;
      samps

let jump_proposal sigmas state = 
  match state with 
    | [| mu; sigma |] -> 
        let new_mu = Array.copy mu and 
            new_sigma = Array.copy sigma in 
          for i = 0 to Array.length new_mu - 1 do 
            let dx = sigmas.(i) /. (sqrt (float_of_int !nsamp)) in 
              new_mu.(i) <- Mcmc.uniform_wrapping !mu_min !mu_max dx new_mu.(i);
              new_sigma.(i) <- Mcmc.uniform_wrapping !sigma_min !sigma_max dx new_sigma.(i)
          done;
          [|new_mu; new_sigma|]
    | _ -> raise (Failure "jump_proposal: bad state")

let log_jump_prob _ _ = 0.0

let log_prior state = 
  match state with 
    | [|mu; sigma|] -> 
        let i = ref 0 and 
            inb = ref true in 
          while (!i < !ndim && !inb) do 
            if mu.(!i) < !mu_min || mu.(!i) > !mu_max || sigma.(!i) < !sigma_min || sigma.(!i) > !sigma_max then 
              inb := false;
            incr i
          done;
          if !inb then 
            ~-.(float_of_int !ndim)*.((log (!mu_max -. !mu_min)) +. (log (!sigma_max -. !sigma_min)))
          else
            neg_infinity
    | _ -> raise (Invalid_argument "log_prior: bad state")

let log_like prob samples state = 
  match state with 
    | [|mu; sigma|] -> 
        let n = Array.length samples and 
            lp = ref 0.0 in 
          for i = 0 to n - 1 do 
            let samp = samples.(i) in 
              for j = 0 to !ndim - 1 do 
                lp := !lp +. prob mu.(j) sigma.(j) samp.(j)
              done
          done;
          !lp +. 0.0
    | _ -> raise (Invalid_argument "log_like: bad state")

let log_like_gaussian samples state = log_like Stats.log_gaussian samples state
let log_like_cauchy samples state = log_like Stats.log_cauchy samples state

module Points_coords = 
  struct
    type point = float array array

    let coord pt = 
      match pt with 
        | [|mu; sigma|] -> Array.append mu sigma
        | _ -> raise (Invalid_argument "coord: bad state")

    let point coords = 
      [|Array.sub coords 0 !ndim;
        Array.sub coords !ndim !ndim|]
  end

module Interp = Interpolate_pdf.Make(Points_coords)
module Ev = Evidence.Make(struct type params = Points_coords.point let to_coords = Points_coords.coord end)

let sample_std samples = 
  let std = Stats.multi_std samples and 
      nsamp = float_of_int (Array.length samples) in 
  let fac = 1.0 /. (sqrt (nsamp -. 1.0)) in 
    for i = 0 to !ndim - 1 do 
      std.(i) <- std.(i) *. fac
    done;
    std

let do_mcmc out log_like nmcmc_samp samples = 
  let sigmas = sample_std samples in
  let jump_proposal state = jump_proposal sigmas state in
  let log_like state = log_like samples state in 
  let samples = 
    Mcmc.mcmc_array nmcmc_samp log_like log_prior jump_proposal log_jump_prob 
      [|Stats.multi_mean samples;
        Stats.multi_std samples|] in 
    for i = 0 to nmcmc_samp - 1 do 
      let {Mcmc.value = samp} = samples.(i) in 
        match samp with 
          | [|mu; sigma|] -> print_array out (Array.append mu sigma)
          | _ -> raise (Failure "do_mcmc: bad sample")
    done;
    samples

let gaussian_mcmc n samples = 
  let out = open_out !gauss_mcmc_fname in 
  let samples = do_mcmc out log_like_gaussian n samples in 
    close_out out;
    samples

let cauchy_mcmc n samples = 
  let out = open_out !cauchy_mcmc_fname in 
  let samples = do_mcmc out log_like_cauchy n samples in 
    close_out out;
    samples

let param_bounds () = 
  (Array.append (Array.make !ndim !mu_min) (Array.make !ndim !sigma_min),
   Array.append (Array.make !ndim !mu_max) (Array.make !ndim !sigma_max))

let array_subsample reduction_factor arr = 
  let res = Array.make ((Array.length arr) / reduction_factor) arr.(0) in 
    for i = 0 to Array.length res - 1 do 
      res.(i) <- arr.(i*reduction_factor)
    done;
    res

let rjmcmc n gsamp csamp data = 
  let (low,high) = param_bounds () in 
  let gsamp = Array.map (fun {Mcmc.value = x} -> x) (array_subsample !interp_reduce gsamp) and 
      csamp = Array.map (fun {Mcmc.value = x} -> x) (array_subsample !interp_reduce csamp) in 
  let gpdf = Interp.make gsamp low high and 
      cpdf = Interp.make csamp low high and 
      mu = Stats.multi_mean data and 
      sigma = Stats.multi_std data in
  let log_like_gaussian state = log_like_gaussian data state and 
      log_like_cauchy state = log_like_cauchy data state in
  let sigmas = sample_std data in
  let jump_proposal state = jump_proposal sigmas state in
  let samples = 
    Mcmc.rjmcmc_array 
      n
      (log_like_gaussian, log_like_cauchy)
      (log_prior, log_prior)
      (jump_proposal, jump_proposal)
      (log_jump_prob, log_jump_prob)
      ((fun () -> Interp.draw gpdf),
       (fun () -> Interp.draw cpdf))
      ((fun gpt -> log (Interp.jump_prob gpdf () gpt)),
       (fun cpt -> log (Interp.jump_prob cpdf () cpt)))
      (!gprior, !cprior)
      ([|mu; sigma|], [|mu; sigma|]) in 
  let out = open_out !rjmcmc_fname in 
    for i = 0 to n - 1 do 
      let {Mcmc.value = samp} = samples.(i) in 
        match samp with 
          | Mcmc.A([|mu; sigma|]) -> 
              Printf.fprintf out "g ";
              print_array out (Array.append mu sigma)
          | Mcmc.B([|mu; sigma|]) -> 
              Printf.fprintf out "c ";
              print_array out (Array.append mu sigma)
          | _ -> raise (Failure "rjmcmc: bad sample")
    done;
    close_out out;
    samples

let admixture_mcmc n data = 
  let log_like_gaussian state = log_like_gaussian data state and 
      log_like_cauchy state = log_like_cauchy data state in 
  let mu = Stats.multi_mean data and 
      sigma = Stats.multi_std data in
  let vol = (!mu_max -. !mu_min)**(float_of_int !ndim)*.(!sigma_max -. !sigma_min)**(float_of_int !ndim) in
  let samples = 
    Mcmc.admixture_mcmc_array
      n
      (log_like_gaussian, log_like_cauchy)
      (log_prior, log_prior)
      (jump_proposal sigma, jump_proposal sigma)
      (log_jump_prob, log_jump_prob)
      (!gprior, !cprior)
      (vol, vol)
      ([|mu; sigma|], [|mu; sigma|]) in 
  let out = open_out !admcmc_fname in 
    for i = 0 to n - 1 do 
      let {Mcmc.value = (lam, gstate, cstate)} = samples.(i) in 
        match gstate, cstate with 
          | [|gmu; gsigma|], [|cmu; csigma|] ->
              Printf.fprintf out "%g " lam;
              print_array out (Array.concat [gmu; gsigma; cmu; csigma])
          | _ -> raise (Failure "admixture_mcmc: bad state")
    done;
    close_out out;
    samples

let accepted_frac na nr = 
  (float_of_int na) /. (float_of_int (na+nr))

let _ = 
  Random.self_init ();
  Arg.parse opts (fun _ -> ()) "gaussian_cauchy.{native,byte} OPT ...";
  let (mu,sigma) = make_params () in 
  let data = draw_samples mu sigma in 
  if !randmcmc then randomize ();
  Mcmc.reset_counters ();
  let gmcmc = gaussian_mcmc !nmcmc data in
  let (gna, gnr) = Mcmc.get_counters () in 
  Mcmc.reset_counters ();
  let cmcmc = cauchy_mcmc !nmcmc data in 
  let (cna, cnr) = Mcmc.get_counters () in
  Mcmc.reset_counters ();
  let rjmcmc = rjmcmc !nmcmc gmcmc cmcmc data in 
  let (rjna, rjnr) = Mcmc.get_counters () in
  Mcmc.reset_counters ();
  let admcmc = admixture_mcmc !nmcmc data in 
  let (adna, adnr) = Mcmc.get_counters () in
  let (ng, nc) = Mcmc.rjmcmc_model_counts rjmcmc in 
  let lam = Stats.meanf (fun {Mcmc.value = (lam,_,_)} -> lam) admcmc in 
  let gevhm = Ev.evidence_harmonic_mean gmcmc and 
      cevhm = Ev.evidence_harmonic_mean cmcmc in 
  let gevlb = Ev.evidence_lebesgue gmcmc and 
      cevlb = Ev.evidence_lebesgue cmcmc in 
  let gevdr = Ev.evidence_direct gmcmc and 
      cevdr = Ev.evidence_direct cmcmc in 
  let max_like_ratio = 
    try 
      Mcmc.max_like_admixture_ratio admcmc 
    with 
      | Assert_failure(_) -> 0.0 /. 0.0 in
    Printf.printf "Accepted %g in gaussian, %g in cauchy MCMC.\n" 
      (accepted_frac gna gnr) (accepted_frac cna cnr);
    Printf.printf "Reverse jump: n_gauss = %d, n_cauchy = %d, ratio = %g, accepted frac = %g\n"
      ng nc ((float_of_int ng)/.(float_of_int nc)) (accepted_frac rjna rjnr);
    Printf.printf "Mean lambda = %g, for ratio of %g, max_like ratio of %g, accepted frac = %g\n" 
      lam (Mcmc.mean_lambda_ratio admcmc) max_like_ratio (accepted_frac adna adnr);
    Printf.printf "Direct evidence: gauss = %g, cauchy = %g, ratio = %g\n"
      gevdr cevdr (gevdr/.cevdr);
    Printf.printf "Harmonic mean evidence: gauss = %g, cauchy = %g, ratio = %g\n" 
      gevhm cevhm (gevhm/.cevhm);
    Printf.printf "Lebesgue evidence: gauss = %g, cauchy = %g, ratio = %g\n"
      gevlb cevlb (gevlb/.cevlb)
