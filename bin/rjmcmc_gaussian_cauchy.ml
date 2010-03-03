let npts = ref 20
let sigma = ref 0.1
let mu = ref 0.5
let nsamp = ref 10000
let nrjsamp = ref 1000000

let options = 
  [("-seed", Arg.Int (fun i -> Random.init i), "set the RNG seed (default is to self_init)");
   ("-npts", Arg.Set_int npts, Printf.sprintf "number of sample points (default: %d)" !npts);
   ("-nsamp", Arg.Set_int nsamp, Printf.sprintf "number of samples in the MCMC's (default: %d)" !nsamp);
   ("-sigma", Arg.Set_float sigma, 
    Printf.sprintf "standard deviation of gaussian from which pts are drawn (default %g)" !sigma);
   ("-mu", Arg.Set_float mu,
    Printf.sprintf "mean of gaussian from which points are drawn (default %g)" !mu);
   ("-nrjsamp", Arg.Set_int nrjsamp, 
    Printf.sprintf "number of samples in reverse-jump MCMC (default %d)" !nrjsamp)]

let random_between a b = 
  a +. (b-.a)*.(Random.float 1.0)

module Interp = Interpolate_pdf.Make(
  struct
    type point = float array
    let coord (p : point) = p
    let point (c : float array) = c
  end)

let _ = 
  Random.self_init ();
  Arg.parse options (fun _ -> ()) "rjmcmc_gaussian_cauchy.{byte,native} OPTIONS ...";
  assert(!mu >= 0.0 && !mu <= 1.0);
  assert(!sigma >= 0.0 && !sigma <= 1.0);
  let pts = Array.init !npts 
    (fun _ -> Stats.draw_gaussian !mu !sigma) in 
    for i = 0 to Array.length pts - 1 do 
      Printf.printf "Found point: %g\n" pts.(i)
    done;
  let mu0 = Stats.mean pts in 
  let sigma0 = Stats.std ~mean:mu0 pts in 
  let log_prior = function 
    | [|mu; sigma|] -> 
        if 0.0 <= mu && mu <= 1.0 && 0.0 <= sigma && sigma <= 1.0 then 
          0.0
        else
          neg_infinity
    | _ -> raise (Invalid_argument "log_prior") in 
  let log_like_gaussian = function 
    | [|mu; sigma|] -> 
        let ll = ref 0.0 in 
          for i = 0 to Array.length pts - 1 do 
            ll := !ll +. Stats.log_gaussian mu sigma pts.(i)
          done;
          !ll
    | _ -> raise (Invalid_argument "log_like_gaussian") in
  let log_like_cauchy = function 
    | [|mu; sigma|] -> 
        let ll = ref 0.0 in 
          for i = 0 to Array.length pts - 1 do 
            ll := !ll +. Stats.log_cauchy mu sigma pts.(i)
          done;
          !ll
    | _ -> raise (Invalid_argument "log_like_cauchy") in 
  let jump = function 
    | [|mu; sigma|] -> 
        let pop_sig = sigma /. (sqrt (float_of_int (Array.length pts))) in 
        let mu' = random_between (mu -. pop_sig) (mu +. pop_sig) and 
            sigma' = random_between (sigma -. pop_sig) (sigma +. pop_sig) in 
          [| mu'; sigma' |]
    | _ -> raise (Invalid_argument "jump") in 
  let log_jump_prob x y = 
    match x,y with 
      | [|mu; sigma|], [|mu'; sigma'|] -> 
          let pop_sig = sigma /. (sqrt (float_of_int (Array.length pts))) in 
            if mu' >= mu -. pop_sig && mu' <= mu +. pop_sig &&
              sigma' >= sigma -. pop_sig && sigma' <= sigma +. pop_sig then 
                ~-.(log (pop_sig *. 2.0))
            else
              neg_infinity
      | _ -> raise (Invalid_argument "log_jump_probability") in 
  let gaussian_samples = Mcmc.mcmc_array !nsamp log_like_gaussian log_prior jump log_jump_prob [|mu0;sigma0|] and 
      cauchy_samples = Mcmc.mcmc_array !nsamp log_like_cauchy log_prior jump log_jump_prob [|mu0; sigma0|] in 
  let gaussian_posterior_interp = 
    Interp.make (Array.map (fun x -> x.Mcmc.value) gaussian_samples) [|0.0; 0.0|] [|1.0; 1.0|] and 
      cauchy_posterior_interp = 
    Interp.make (Array.map (fun x -> x.Mcmc.value) cauchy_samples) [|0.0; 0.0|] [|1.0; 1.0|] in 
  let gaussian_mean = Stats.multi_mean (Array.map (fun x -> x.Mcmc.value) gaussian_samples) and 
      cauchy_mean = Stats.multi_mean (Array.map (fun x -> x.Mcmc.value) cauchy_samples) in 
    Printf.printf "gaussian mean parameters: [|%g; %g|]\ncauchy mean parameters: [|%g; %g|]\nactual parameters: [|%g; %g|]\n"
      gaussian_mean.(0) gaussian_mean.(1) cauchy_mean.(0) cauchy_mean.(1) !mu !sigma;
    let next_sample = 
      Mcmc.make_rjmcmc_sampler
        (log_like_gaussian,log_like_cauchy)
        (log_prior,log_prior)
        (jump,jump)
        (log_jump_prob,log_jump_prob)
        ((fun () -> Interp.draw gaussian_posterior_interp), (fun () -> Interp.draw cauchy_posterior_interp))
        ((fun x -> log (Interp.jump_prob gaussian_posterior_interp () x)),
         (fun x -> log (Interp.jump_prob cauchy_posterior_interp () x)))
        (0.5,0.5) and 
        v0 = [|mu0; sigma0|] and
        ng = ref 0 and 
        nc = ref 0 in 
    let sample = ref {Mcmc.value = Mcmc.A(v0);
                      like_prior = {Mcmc.log_likelihood = log_like_gaussian v0;
                                    log_prior = (log 0.5) +. log_prior v0}} in
      for i = 0 to !nrjsamp - 1 do 
        (match (!sample).Mcmc.value with 
           | Mcmc.A(_) -> incr ng
           | Mcmc.B(_) -> incr nc);
        sample := next_sample !sample
      done;      
      Printf.printf "%d gaussian counts, %d cauchy counts; ratio = %g\n" !ng !nc ((float_of_int !ng)/.(float_of_int !nc))
