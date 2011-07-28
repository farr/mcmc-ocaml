open OUnit

let tests = "all denest tests" >:::
  [Earray_test.tests;
   Evidence_test.tests;
   Interpolate_pdf_test.tests;
   Kd_tree_test.tests;
   Mcmc_test.tests;
   Read_write_test.tests;
   Stats_test.tests]

let _ = 
  let inp = open_in "/dev/random" in 
  let seed = input_binary_int inp in 
    Random.init seed;
    close_in inp

let _ = 
  let results = run_test_tt_main tests in 
  let nfail = 
    List.fold_left 
      (fun nfail result -> 
         match result with 
           | RSuccess(_) -> nfail
           | _ -> nfail + 1)
      0
      results in 
    exit nfail
