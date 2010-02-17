open OUnit

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

let log_uniform_prior x = 
  let i = ref 0 and 
      outside = ref false and
      n = Array.length x in 
    while (!i < n && not !outside) do
      if x.(!i) < 0.0 || x.(!i) > 1.0 then 
        outside := true;
      incr i
    done;
    if !outside then 
      neg_infinity
    else
      0.0

let multi_gaussian_propose sigma x = 
  let n = Array.length x in 
  let ndim = float_of_int n in 
  let xnew = Array.make n 0.0 in 
    for i = 0 to n - 1 do 
      let dx = sigma.(i) /. ndim in 
        xnew.(i) <- x.(i) +. (random_between (~-.dx) dx)
    done;
    xnew

let assert_equal_float ?(epsabs = 1e-8) ?(epsrel = 1e-8) = 
  assert_equal ~printer:string_of_float ~cmp:(cmp_float ~epsabs:epsabs ~epsrel:epsrel)

let test_evidence_direct_2d () = 
  let mu = Array.init 2 (fun _ -> Random.float 1.0) and 
      sigma = Array.init 2 (fun _ -> Random.float 1.0 +. 1.0) in 
  let all_samples = 
    (Mcmc.mcmc_array 10000 (log_multi_gaussian mu sigma) (fun _ -> 0.0) 
       (multi_gaussian_propose sigma) (fun _ _ -> 0.0) mu) in 
  let samples = Mcmc.remove_repeat_samples (=) all_samples in 
  let ev = Ev.evidence_direct ~n:64 samples in 
    assert_equal_float ~epsabs:0.5 ~msg:"evidence not 1 for gaussian posterior" 1.0 ev

let test_evidence_harmonic_mean_2d () = 
  let mu = Array.init 2 (fun _ -> Random.float 1.0) and 
      sigma = Array.init 2 (fun _ -> Random.float 0.1) in 
  let all_samples = 
    (Mcmc.mcmc_array 100000 (log_multi_gaussian mu sigma) log_uniform_prior 
       (multi_gaussian_propose sigma) (fun _ _ -> 0.0) mu) in 
  let ev = Ev.evidence_harmonic_mean all_samples in 
    assert_equal_float ~epsabs:0.5 ~msg:"evidence not 1 for gaussian posterior" 1.0 ev  

let tests = "evidence.ml tests" >:::
  ["evidence_direct in 2D" >:: test_evidence_direct_2d]
