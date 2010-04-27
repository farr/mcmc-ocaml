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
  let std_error = 10.0*.sigma/.(sqrt (float_of_int nsamp)) in
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
      assert_equal_float ~epsrel:0.0 ~epsabs:std_error ~msg:"means differ" mu mean;
      assert_equal_float ~epsrel:0.0 ~epsabs:std_error ~msg:"sigmas differ" sigma std
      
let test_gaussian_post_left_biased_proposal () = 
  let nsamp = 100000 and 
      mu = Random.float 1.0 and 
      sigma = Random.float 1.0 +. 1.0 in 
  let std_error = 20.0*.sigma/.(sqrt (float_of_int nsamp)) in
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
      assert_equal_float ~epsrel:0.0 ~epsabs:std_error ~msg:"means differ" mu mean;
      assert_equal_float ~epsrel:0.0 ~epsabs:std_error ~msg:"sigmas differ" sigma std

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

let test_rjmcmc_gaussians () = 
  let nsamp = 1000000 and 
      mu1 = Random.float 1.0 and 
      sigma1 = Random.float 1.0 and 
      mu2 = Random.float 1.0 and 
      sigma2 = Random.float 1.0 in
  let propose1 x = Stats.draw_gaussian mu1 sigma1 and 
      propose2 x = Stats.draw_gaussian mu2 sigma2 and 
      propose_into1 () = Stats.draw_gaussian mu1 sigma1 and 
      propose_into2 () = Stats.draw_gaussian mu2 sigma2 in 
  let log_jump1 _ y = Stats.log_gaussian mu1 sigma1 y and 
      log_jump2 _ y = Stats.log_gaussian mu2 sigma2 y and 
      log_jump_into1 x = Stats.log_gaussian mu1 sigma1 x and 
      log_jump_into2 x = Stats.log_gaussian mu2 sigma2 x in 
  let log_like1 x = 0.5 *. (log_gaussian mu1 sigma1 x) and 
      log_like2 x = 0.3 *. (log_gaussian mu2 sigma2 x) and 
      log_prior1 x = 0.5 *. (log_gaussian mu1 sigma1 x) and 
      log_prior2 x = 0.7 *. (log_gaussian mu2 sigma2 x) in 
  let p1 = 0.5 and p2 = 0.5 in 
  let samples = rjmcmc_array nsamp
    (log_like1,log_like2)
    (log_prior1,log_prior2)
    (propose1,propose2)
    (log_jump1,log_jump2)
    (propose_into1,propose_into2)
    (log_jump_into1,log_jump_into2)
    (p1,p2)
    (mu1, mu2) in 
  let (n1, n2) = rjmcmc_model_counts samples in 
  let p1 = (float_of_int n1) /. (float_of_int (n1+n2)) and 
      p2 = (float_of_int n2) /. (float_of_int (n1+n2)) in 
    assert_equal_float ~epsrel:0.1 0.5 p1;
    assert_equal_float ~epsrel:0.1 0.5 p2

let test_admixture_gaussian_cauchy () = 
  let mu = Random.float 1.0 and 
      sigma = Random.float 1.0 +. 1.0 in 
  let ndata = 50 and 
      nlambda = 10000 in
  let samples = Array.init ndata (fun _ -> Stats.draw_gaussian mu sigma) in 
  let log_gaussian_like ms = 
    match ms with 
      | [|mu; sigma|] -> 
          let ll = ref 0.0 in 
            for i = 0 to Array.length samples - 1 do 
              ll := !ll +. log_gaussian mu sigma samples.(i)
           done;
            !ll +. 0.0
      | _ -> raise (Invalid_argument "log_gaussian_like") in
  let log_cauchy_like ms = 
    match ms with 
      | [|mu; sigma|] -> 
          let ll = ref 0.0 in 
            for i = 0 to Array.length samples - 1 do 
              ll := !ll +. Stats.log_cauchy mu sigma samples.(i)
            done;
            !ll +. 0.0
      | _ -> raise (Invalid_argument "log_cauchy_like") in 
  let log_prior ms = 
    match ms with 
      | [|mu; sigma|] -> 
          if 0.0 <= mu && mu <= 1.0 && 1.0 <= sigma && sigma <= 2.0 then 0.0 else neg_infinity 
      | _ -> raise (Invalid_argument "log_prior") in
  let propose ms = 
    match ms with 
      | [|mu; sigma|] -> 
          let dmu = sigma /. (sqrt (float_of_int ndata)) and 
              dsigma = sigma /. (sqrt (float_of_int ndata)) in 
            [|random_between (mu -. dmu) (mu +. dmu);
              random_between (sigma -. dsigma) (sigma +. dsigma)|]
      | _ -> raise (Invalid_argument "propose") in 
  let log_jump_prob _ _ = 0.0 in 
  let lambda_samples = 
    admixture_mcmc_array nlambda (log_gaussian_like, log_cauchy_like) (log_prior, log_prior) 
      (propose, propose) (log_jump_prob, log_jump_prob) (0.5, 0.5) (1.0,1.0) 
      ([|mu; sigma|], [|mu; sigma|]) in 
  let avg_lambda = Stats.meanf (fun {Mcmc.value = (lam,_,_)} -> lam) lambda_samples in 
    assert_bool "lambda small, favors cauchy" (avg_lambda > 0.5)

let test_admixture_lambda_dist () = 
  let xmax = 10.0 and 
      xmin = -10.0 in 
  let nsamples = 1000000 in
  let ea = Random.float 1.0 and 
      eb = Random.float 1.0 in 
  let r = ea /. eb in
  let el = (1.0 +. 2.0 *. r) /. (3.0 +. 3.0*.r) and 
      el2 = 2.0/.(r +. 1.0)*.(r /. 4.0 +. 1.0/.12.0) in
  let sigma_l = sqrt (el2 -. el*.el) in
  let std_error = sigma_l /. (sqrt (float_of_int nsamples)) in
  let log_like1 x = Stats.log_gaussian 0.0 1.0 x +. (log ea) in
  let log_like2 x = Stats.log_gaussian 0.0 1.0 x +. (log eb) in 
  let log_prior x = 
    if xmin <= x && x <= xmax then 
      ~-.(log (xmax -. xmin))
    else
      neg_infinity in
  let propose x = 
    let dx = 0.5 in 
    let xpdx = x +. dx and 
        xmdx = x -. dx in 
    let xmax = min xmax xpdx and 
        xmin = max xmin xmdx in 
      random_between xmin xmax in
  let log_jump_prob x y = 
    let dx = 0.5 in 
    let xpdx = x +. dx and 
        xmdx = x -. dx in 
    let xmax = min xmax xpdx and 
        xmin = max xmin xmdx in 
      ~-.(log (xmax -. xmin)) in
  Mcmc.reset_counters ();
  let samples = 
    admixture_mcmc_array nsamples (log_like1, log_like2) (log_prior, log_prior) 
      (propose, propose) (log_jump_prob, log_jump_prob) (0.5, 0.5) (xmax -. xmin, xmax -. xmin) (0.0, 0.0) in 
  let lam = Stats.meanf (fun {Mcmc.value = (lam,_,_)} -> lam) samples in 
    assert_equal_float ~epsabs:(50.0*.std_error) ~msg:"mean lambda differs" el lam

let test_combine_jump_proposal () = 
  let nsamples = 1000000 in 
  let propose_left x =
    let xnew = x -. (Random.float 1.0) in 
      xnew and
      propose_right x = x +. (Random.float 1.0) in 
  let log_jp_left x y = 
    if y <= x && y >= x -. 1.0 then begin
      0.0 
    end else
      neg_infinity and
      log_jp_right x y = 
    if y >= x && y <= x +. 1.0 then 
      0.0
    else
      neg_infinity in
  let (propose, log_jp) = 
    Mcmc.combine_jump_proposals 
      [(1.0, propose_left, log_jp_left);
       (2.0, propose_right, log_jp_right)] in 
  let samps = Mcmc.mcmc_array nsamples (fun x -> Stats.log_gaussian 0.0 1.0 x) (fun x -> 0.0) propose log_jp 0.0 in 
  let mu = Stats.meanf (fun {Mcmc.value = x} -> x) samps and 
      sigma = Stats.stdf (fun {Mcmc.value = x} -> x) samps in 
    assert_equal_float ~epsabs:0.05 0.0 mu;
    assert_equal_float ~epsrel:1e-2 1.0 sigma

let test_max_like_admixture_ratio () = 
  let nsamples = 1000000 in 
  let r = 0.5 +. Random.float 1.0 in 
  let samples = 
    Array.init nsamples
      (fun _ -> 
         let rec lam_loop () = 
           let lam = Random.float 1.0 in 
           let plam = r *. lam +. 1.0 -. lam in 
           if Random.float (r +. 1.0) < plam then 
             {Mcmc.value = (lam, (), ());
              like_prior = {Mcmc.log_likelihood = 0.0; log_prior = 0.0}}
           else
             lam_loop () in 
           lam_loop ()) in 
  let mlr = Mcmc.max_like_admixture_ratio samples in 
    assert_equal_float ~epsrel:1e-2 r mlr

let tests = "mcmc.ml tests" >:::
  ["gaussian posterior, uniform jump proposal" >:: test_gaussian_post_uniform_proposal;
   "gaussian posterior, left-biased jump proposal" >:: test_gaussian_post_left_biased_proposal;
   "prior*like = gaussian, uniform jump" >:: test_prior_like;
   "remove_repeat" >:: test_remove_repeat;
   "rjmcmc on gaussian posteriors in 1-D" >:: test_rjmcmc_gaussians;
   "admixture gaussian vs cauchy test" >:: test_admixture_gaussian_cauchy;
   "admixture lambda distribution" >:: test_admixture_lambda_dist;
   "combine_jump_proposal" >:: test_combine_jump_proposal;
   "max_like_admixture_ratio" >:: test_max_like_admixture_ratio]
