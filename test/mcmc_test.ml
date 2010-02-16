open OUnit
open Mcmc

let pi = 4.0 *. (atan 1.0)

let assert_equal_float ?(epsabs = 1e-8) ?(epsrel = 1e-8) = 
  assert_equal ~printer:string_of_float ~cmp:(cmp_float ~epsabs:epsabs ~epsrel:epsrel)

let log_gaussian mu sigma x = 
  let d = mu -. x and
      log_two_pi = (-0.5)*.(log (2.0*.pi)) and 
      log_sigma = ~-. (log sigma) in 
  let log_exp_d2 = ~-.d*.d/.(2.0*.sigma*.sigma) in 
    log_two_pi +. log_sigma +. log_exp_d2

let random_between a b = 
  a +. (b -. a)*.(Random.float 1.0)

let mean samples = 
  let sum = ref 0.0 and 
      n = Array.length samples in 
    for i = 0 to n - 1 do 
      sum := !sum +. samples.(i).value
    done;
    !sum /. (float_of_int n)

let std ?mu samples = 
  let mu = 
    match mu with 
      | None -> mean samples 
      | Some(x) -> x in 
  let sum = ref 0.0 and
      n = Array.length samples in 
    for i = 0 to n - 1 do 
      let dx = samples.(i).value -. mu in 
        sum := !sum +. dx*.dx
    done;
    sqrt (!sum /. (float_of_int n))

let test_gaussian_post_uniform_proposal () = 
  let nsamp = 100000 and 
      mu = Random.float 1.0 and 
      sigma = Random.float 1.0 +. 1.0 in 
  let propose x = x +. (random_between (~-.sigma) sigma) and 
      posterior x = log_gaussian mu sigma x and 
      prior x = 0.0 and
      jump_prob x y = 0.0 in 
  let mcmc_sample = make_mcmc_sampler posterior prior propose jump_prob in 
  let samples = 
    Array.make nsamp {value = mu; 
                      like_prior = {log_likelihood = posterior mu; log_prior = prior mu}} in 
    for i = 1 to nsamp - 1 do 
      samples.(i) <- mcmc_sample samples.(i-1)
    done;
    let mean = mean samples in 
    let std = std ~mu:mean samples in 
      assert_equal_float ~epsrel:1e-1 ~epsabs:0.0 ~msg:"means differ" mu mean;
      assert_equal_float ~epsrel:1e-1 ~epsabs:0.0 ~msg:"sigmas differ" sigma std
      
let test_gaussian_post_left_biased_proposal () = 
  let nsamp = 100000 and 
      mu = Random.float 1.0 and 
      sigma = Random.float 1.0 +. 1.0 in 
  let propose x = 
    if Random.float 1.0 < 0.75 then 
      x -. Random.float sigma
    else
      x +. Random.float sigma in 
  let posterior x = log_gaussian mu sigma x and 
      prior x = 0.0 and
      jump_prob x y = if x > y then (log 0.75) else (log 0.25) in 
  let mcmc_sample = make_mcmc_sampler posterior prior propose jump_prob in 
  let samples = 
    Array.make nsamp {value = mu; 
                      like_prior = {log_likelihood = posterior mu; log_prior = prior mu}} in 
    for i = 1 to nsamp - 1 do 
      samples.(i) <- mcmc_sample samples.(i-1)
    done;
    let mean = mean samples in
    let std = std ~mu:mean samples in 
      assert_equal_float ~epsrel:2e-1 ~epsabs:0.0 ~msg:"means differ" mu mean;
      assert_equal_float ~epsrel:2e-1 ~epsabs:0.0 ~msg:"sigmas differ" sigma std

let test_prior_like () = 
  let nsamp = 100000 and 
      mu = Random.float 1.0 and 
      sigma = Random.float 1.0 +. 1.0 in 
  let propose x = x +. random_between (~-.sigma) sigma and 
      prior x = 0.25 *. (log_gaussian mu sigma x) and 
      like x = 0.75 *. (log_gaussian mu sigma x) and 
      jump_prob x y = 0.0 in
  let samples = mcmc_array nsamp like prior propose jump_prob mu in 
  let mean = mean samples in 
  let std = std ~mu:mean samples in 
    assert_equal_float ~msg:"means differ" ~epsrel:2e-1 ~epsabs:0.0 mu mean;
    assert_equal_float ~msg:"sigmas differ" ~epsrel:2e-1 ~epsabs:0.0 sigma std

let test_remove_repeat () = 
  let nsamp = 1000 and 
      mu = Random.float 1.0 and 
      sigma = Random.float 1.0 +. 1.0 in 
  let propose x = x +. random_between (~-.sigma) sigma and 
      prior x = 0.25 *. (log_gaussian mu sigma x) and 
      like x = 0.75 *. (log_gaussian mu sigma x) and 
      jump_prob x y = 0.0 in
  let samples = mcmc_array nsamp like prior propose jump_prob mu in 
  let no_repeats = remove_repeat_samples (=) samples in 
    for i = 0 to Array.length no_repeats - 2 do 
      assert_bool "repeated sample!" (no_repeats.(i).value <> no_repeats.(i+1).value)
    done

let tests = "mcmc.ml tests" >:::
  ["gaussian posterior, uniform jump proposal" >:: test_gaussian_post_uniform_proposal;
   "gaussian posterior, left-biased jump proposal" >:: test_gaussian_post_left_biased_proposal;
   "prior*like = gaussian, uniform jump" >:: test_prior_like;
   "remove_repeat" >:: test_remove_repeat]
