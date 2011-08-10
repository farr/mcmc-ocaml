open OUnit
open Asserts

let mcmc_samples_close a b = 
  let n = Array.length a in 
    assert_equal ~printer:string_of_int ~msg:"lengths differ" n (Array.length b);
    for i = 0 to n - 1 do 
      let {Mcmc.value = va;
           like_prior = {Mcmc.log_likelihood = lla;
                         log_prior = lpa}} = a.(i) and 
          {Mcmc.value = vb;
           like_prior = {Mcmc.log_likelihood = llb;
                         log_prior = lpb}} = b.(i) in 
        assert_equal ~printer:string_of_float ~msg:"values differ" 
          ~cmp:(cmp_float ~epsabs:1e-3 ~epsrel:1e-3) va vb;
        assert_equal ~printer:string_of_float ~msg:"log_likelihoods differ" 
          ~cmp:(cmp_float ~epsabs:1e-3 ~epsrel:1e-3) lla llb;
        assert_equal ~printer:string_of_float ~msg:"values differ" 
          ~cmp:(cmp_float ~epsabs:1e-3 ~epsrel:1e-3) lpa lpb
    done

let test_read_write_inverses () = 
  let samples = 
    Mcmc.mcmc_array 1000
      (fun x -> Stats.log_gaussian 0.0 1.0 x)
      (fun _ -> 0.0)
      (fun x -> x +. Random.float 1.0 -. 0.5)
      (fun x y -> 0.0)
      0.0 in
  let (fname, ofile) = Filename.open_temp_file "mcmc_test" ".dat" in 
    (try 
       Read_write.write (fun x -> [| x |]) ofile samples
     with 
       | x -> close_out ofile; Sys.remove fname; raise x);
    close_out ofile;
    let ifile = open_in fname in 
    let rsamples = 
      try
        Read_write.read (fun x -> x.(0)) ifile
      with 
        | x -> close_in ifile; Sys.remove fname; raise x in 
      close_in ifile;
      Sys.remove fname;
      mcmc_samples_close samples rsamples
    
let test_nested_read_write () = 
  let (log_ev, log_dev, arr_samples, log_wts) as result = 
    Nested.nested_evidence 
      (fun () -> [|Random.float 1.0|])
      (fun x -> Stats.log_gaussian 0.5 0.1 x.(0))
      (fun x -> 0.0) in 
  let samples = Array.map (fun x -> {x with Mcmc.value = x.Mcmc.value.(0)}) arr_samples in 
  let (fname, ofile) = Filename.open_temp_file "nested_test" ".dat" in 
    (try
       Read_write.write_nested ofile result
     with 
       | x -> close_out ofile; Sys.remove fname; raise x);
    close_out ofile;
    let ifile = open_in fname in 
    let (rlog_ev, rlog_dev, rarr_samples, rlog_wts) = 
      try 
        Read_write.read_nested ifile
      with 
        | x -> close_in ifile; Sys.remove fname; raise x in 
    let rsamples = Array.map (fun x -> {x with Mcmc.value = x.Mcmc.value.(0)}) rarr_samples in 
      close_in ifile;
      Sys.remove fname;
      assert_equal_float ~epsrel:0.01 ~msg:"evidence disagrees" log_ev rlog_ev;
      assert_equal_float ~epsrel:0.01 ~msg:"delta evidence disagrees" log_dev rlog_dev;
      mcmc_samples_close samples rsamples;
      assert_equal_float_array ~epsrel:0.01 ~msg:"weights disagree" log_wts rlog_wts
    
let tests = "read_write.ml tests" >:::
  ["read and write are inverses" >:: test_read_write_inverses;
   "read and write nested samples" >:: test_nested_read_write]
