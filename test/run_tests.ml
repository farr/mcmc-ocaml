open OUnit

let tests = "all denest tests" >:::
  ["kd_tree.ml tests" >: Kd_tree_test.tests]

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
