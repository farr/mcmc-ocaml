open OUnit

let tests = "all denest tests" >:::
  ["evidence.ml tests" >: Evidence_test.tests;
   "interpolate_pdf.ml tests" >: Interpolate_pdf_test.tests;
   "kd_tree.ml tests" >: Kd_tree_test.tests;
   "mcmc.ml tests" >: Mcmc_test.tests;
   "stats.ml tests" >: Stats_test.tests]

let _ = 
  Random.self_init ();
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
