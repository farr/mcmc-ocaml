open OUnit

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
    

let tests = "read_write.ml tests" >:::
  ["read and write are inverses" >:: test_read_write_inverses]
