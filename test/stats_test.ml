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

let test_multi_std () = 
  let xs = [|[|0.662891; 0.218155; 0.464706; 0.148477; 0.39616|];
             [|0.43397; 0.161041; 0.625332; 0.508765; 0.261084|];
             [|0.147267; 0.403388; 0.643601; 0.892214; 0.269893|]|] in 
  let std = multi_std xs in 
    assert_equal_float_array ~epsabs:1e-3 ~epsrel:1e-3 
      [|0.258351; 0.126692; 0.098436; 0.371928; 0.0755716|]
      std

let test_meanf_stdf () = 
  let xs = Array.init 1000 (fun _ -> Random.float 1.0) in 
    assert_equal_float (mean xs) (meanf (fun x -> x) xs);
    assert_equal_float (std ~mean:(mean xs) xs) (stdf ~mean:(mean xs) (fun x -> x) xs)

let fcompare (x : float) y = Pervasives.compare x y

let test_find_nth () = 
  for i = 0 to 100 do 
    let n = 10 + Random.int 25 in 
    let xs = Array.init n (fun _ -> Random.float 1.0) in 
    let sxs = Array.copy xs in 
    let cxs = Array.copy xs in 
      Array.fast_sort fcompare sxs;
      let i = Random.int n in
        assert_equal ~cmp:(=) ~printer:string_of_float sxs.(i) (find_nth ~copy:true i xs);
        assert_bool "xs changed by find_nth" (xs = cxs);
        assert_equal ~cmp:(=) ~printer:string_of_float sxs.(i) (find_nth ~copy:false i xs);
        assert_bool "xs not disordered" (not (xs = cxs))
  done

let icompare (x : int) y = Pervasives.compare x y

let test_find_nthf () = 
  for i = 0 to 100 do 
    let n = 10 + Random.int 25 in 
    let xs = Array.init n (fun _ -> Random.int 100) in
    let sxs = Array.copy xs and 
        cxs = Array.copy xs in 
      Array.fast_sort icompare sxs;
      let i = Random.int n in 
        assert_equal ~cmp:(=) ~printer:string_of_int sxs.(i) (find_nthf ~copy:true icompare i xs);
        assert_bool "xs changed by find_nthf" (xs = cxs);
        assert_equal ~cmp:(=) ~printer:string_of_int sxs.(i) (find_nthf ~copy:false icompare i xs);
        assert_bool "xs not disordered" (not (xs = cxs))
  done

let test_log_lognormal () = 
  let mu = 0.328077 and 
      sigma = 0.330877 and
      x = 0.0553941 in
  let ln = log_lognormal mu sigma x in
    assert_equal_float ~epsabs:1e-3 (-44.3128) ln

let test_slow_autocorrelation () = 
  let data = 
    [|-1.44898; -0.0762953; 2.25525; -0.284584; 1.16297; 0.00864677; 0.211493|] in 
  let ac = Stats.slow_autocorrelation 3 data in 
    assert_equal_float_array ~epsrel:1e-5 
      [| 1.000000000000000; -0.117926565125732; -0.043808623348460|]
      ac

let tests = "stats.ml tests" >:::
  ["mean" >:: test_mean;
   "randomized mean" >:: test_mean_random;
   "std" >:: test_std;
   "randomized std" >:: test_std_random;
   "gaussian dist" >:: test_gaussian;
   "randomized gaussian draws" >:: test_gaussian_random;
   "multi_mean" >:: test_multi_mean;
   "meanf and stdf" >:: test_meanf_stdf;
   "multi_std" >:: test_multi_std;
   "find_nth" >:: test_find_nth;
   "find_nthf" >:: test_find_nthf;
   "log_lognormal" >:: test_log_lognormal]
