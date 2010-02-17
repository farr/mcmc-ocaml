module Ev = Evidence.Make(
  struct 
    type params = float array
    let to_coords (p : params) = p
  end)

let pi = 4.0 *. (atan 1.0)

let random_between a b = 
  a +. (Random.float (b-.a))

let log_multi_gaussian mu sigma x = 
  let n = Array.length x in 
  let res = ref 0.0 in 
    for i = 0 to n - 1 do 
      let dx = mu.(i) -. x.(i) in 
        res := !res -. 0.5*.(log (2.0*.pi));
        res := !res -. (log sigma.(i));
        res := !res -. dx*.dx/.(2.0*.sigma.(i)*.sigma.(i))
    done;
    !res +. 0.0

let multi_gaussian_propose sigma x = 
  Array.mapi
    (fun i x -> 
       x +. (random_between (~-.(sigma.(i))) sigma.(i)))
    x

let nsamp = ref 1000000
let nbin = ref 64
let ndim = ref 10

let clspecs = 
  [("-nsamp", Arg.Set_int nsamp, (Printf.sprintf "number of samples in MCMC (default: %d)" !nsamp));
   ("-nbin", Arg.Set_int nbin, (Printf.sprintf "number of samples per unit volume in integration (default: %d)" !nbin));
   ("-ndim", Arg.Set_int ndim, (Printf.sprintf "number of dimensions in MCMC (default: %d)" !ndim))]

let _ = 
  Random.self_init ();
  Arg.parse clspecs (fun _ -> ()) "evidence_direct OPTION ...";
  let mu = Array.init !ndim (fun _ -> Random.float 1.0) and 
      sigma = Array.init !ndim (fun _ -> Random.float 1.0) in 
  let log_like x = log_multi_gaussian mu sigma x and 
      log_prior x = 0.0 and 
      jump_proposal x = multi_gaussian_propose sigma x and 
      log_jump_prob x y = 0.0 in 
  let samples = Mcmc.mcmc_array !nsamp log_like log_prior jump_proposal log_jump_prob mu in 
  let samples = Mcmc.remove_repeat_samples (=) samples in 
  let ev = Ev.evidence_direct ~n:(!nbin) samples in 
    Printf.printf "%d %g\n%!" (!nbin) ev
