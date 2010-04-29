(** Testing the various methods for evidence determination. *)

let ndim = ref 1
let nsamp = ref 10
let nmcmc = ref 1000000
let nevid = ref 100
let nsubsamp = ref 10000

let opts = 
  [("-ndim", Arg.Set_int ndim,
    (Printf.sprintf "Number of dimensions to draw gaussian from (default %d)" !ndim));
   ("-nsamp", Arg.Set_int nsamp,
    (Printf.sprintf "Number of samples to draw from a gaussian (default %d)" !nsamp));
   ("-nmcmc", Arg.Set_int nmcmc,
    (Printf.sprintf "Number of MCMC samples to draw (default %d)" !nmcmc));
   ("-nevid", Arg.Set_int nevid,
    (Printf.sprintf "Number of evidence calculations to perform (default %d)" !nevid));
   ("-nsubsamp", Arg.Set_int nsubsamp,
    (Printf.sprintf "Number of sub-samples to draw from the MCMC for each evidence calculation (default %d)"
       !nsubsamp))]

let _ = Random.self_init ()
let _ = Arg.parse opts (fun _ -> ()) 
  "compare_evidence OPTION ..."

let mus = Array.init !ndim (fun _ -> Random.float 1.0)
let sigmas = Array.init !ndim (fun _ -> Random.float 1.0)

let samples = 
  Array.init !nsamp 
    (fun _ -> 
       Array.init !ndim
         (fun i -> 
            Stats.draw_gaussian mus.(i) sigmas.(i)))

let _ = 
  let out = open_out "samples.dat" in 
    Array.iter 
      (fun samp -> 
         Array.iter 
           (fun s -> Printf.fprintf out "%g " s)
           samp;
         Printf.fprintf out "\n")
      samples;
    close_out out

let log_likelihood (mu, sigma) = 
  let ll = ref 0.0 in 
    for i = 0 to !nsamp - 1 do 
      let samp = samples.(i) in 
        for j = 0 to !ndim - 1 do 
          let x = samp.(j) and 
              mu = mu.(j) and 
              sigma = sigma.(j) in 
            ll := !ll +. Stats.log_gaussian mu sigma x
        done
    done;
    !ll +. 0.0

let log_prior _ = 0.0

let jump_proposal (mu,sigma) = 
  let scale = 1.0 /. (sqrt (float_of_int !nsamp)) in
  let new_mu = Array.copy mu and 
      new_sigma = Array.copy sigma in 
    for i = 0 to !ndim - 1 do 
      let dx = scale *. sigmas.(i) in 
      new_mu.(i) <- Mcmc.uniform_wrapping 0.0 1.0 dx new_mu.(i);
      new_sigma.(i) <- Mcmc.uniform_wrapping 0.0 1.0 dx new_sigma.(i)
    done;
    (new_mu, new_sigma)

let log_jump_prob _ _ = 0.0

let mcmcs = Mcmc.mcmc_array !nmcmc log_likelihood log_prior jump_proposal log_jump_prob (mus, sigmas)

let _ = 
  let out = open_out "mcmc.dat" in
  Array.iter 
    (fun {Mcmc.value = (mu,sigma)} -> 
       Array.iter (fun mu -> Printf.fprintf out "%g " mu) mu;
       Array.iter (fun s -> Printf.fprintf out "%g " s) sigma;
       Printf.fprintf out "\n")
    mcmcs;
    close_out out

let generate_subsample () = 
  Array.init 
    !nsubsamp
    (fun _ -> 
       mcmcs.(Random.int (Array.length mcmcs)))

module Ev = Evidence.Make(
  struct
    type params = float array * float array

    let to_coords (mu,sigma) = 
      Array.append mu sigma
  end)

let _ = 
  let out = open_out "evidences.dat" in
    for i = 0 to !nevid - 1 do 
      let mcmc = generate_subsample () in 
        Printf.fprintf out "%g %g %g\n" 
          (Ev.evidence_harmonic_mean mcmc) (Ev.evidence_lebesgue mcmc) (Ev.evidence_direct mcmc)
    done;
    close_out out            
