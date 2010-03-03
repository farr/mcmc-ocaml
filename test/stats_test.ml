open Asserts
open OUnit
open Stats

let test_mean () = 
  let xs = [|0.0; 1.0; 2.0; 3.0|] in 
    assert_equal_float (6.0/.4.0) (mean xs)

let test_mean_random () = 
  let xs = Array.init 10000 (fun _ -> Random.float 1.0) in 
    assert_equal_float ~epsrel:1e-1 0.5 (mean xs)

let test_std () = 
  let xs = [|1.0; 2.0; 3.0; 4.0; 5.0|] in 
    assert_equal_float ((sqrt 10.0)/.2.0) (std xs)

let test_std_random () = 
  let xs = Array.init 10000 (fun _ -> Random.float 1.0) in 
    assert_equal_float ~epsrel:1e-1 (sqrt (1.0/.3.0 -. 0.25)) (std xs)

let pi = 4.0*.(atan 1.0)

let test_gaussian () = 
  let gaus mu sigma x = exp (log_gaussian mu sigma x) in 
  let mu = Random.float 1.0 and 
      sigma = Random.float 1.0 in 
  let g0 = (1.0/.((sqrt (2.0*.pi))*.sigma)) in 
    assert_equal_float g0 (gaus mu sigma mu);
    assert_equal_float (g0*.(exp (-0.5))) (gaus mu sigma (mu+.sigma))

let test_gaussian_random () = 
  let mu = Random.float 1.0 and 
      sigma = Random.float 1.0 in 
  let xs = Array.init 10000 (fun _ -> draw_gaussian mu sigma) in 
  let mean = mean xs and 
      std = std xs in 
    assert_equal_float ~epsabs:0.1 mu mean;
    assert_equal_float ~epsabs:0.1 sigma std

let test_multi_mean () = 
  let xs = [|[|0.0; 1.0|];
             [|2.0; 3.0|];
             [|4.0; -5.0|]|] in 
  let mu = multi_mean xs in 
    assert_equal_float 2.0 mu.(0);
    assert_equal_float (-1.0/.3.0) mu.(1)

let tests = "stats.ml tests" >:::
  ["mean" >:: test_mean;
   "randomized mean" >:: test_mean_random;
   "std" >:: test_std;
   "randomized std" >:: test_std_random;
   "gaussian dist" >:: test_gaussian;
   "randomized gaussian draws" >:: test_gaussian_random;
   "multi_mean" >:: test_multi_mean]
