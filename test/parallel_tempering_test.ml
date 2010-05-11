(** Tests the parallel-tempering code in mcmc.ml, which uses MPI to
    run a number of chains in parallel.  Execute with 

    mpiexec -n <number-of-processes> ./parallel_tempering_test.{native,byte}

    The posterior being simulated is a 1-D gaussian; each chain stores
    its samples into the file pt_test_<beta>_samples.dat, and also
    prints to stdout the thermodynamic integration calculation of the
    evidence.
*)

let _ = 
  let inp = open_in_bin "/dev/random" in 
  let seed = input_binary_int inp in 
    close_in inp;
    Random.init seed

let log_likelihood x = Stats.log_gaussian 0.0 1.0 x
let log_prior x = if abs_float x <= 10.0 then log 0.1 else 0.0

let jump_proposal x = Mcmc.uniform_wrapping (-10.0) 10.0 1.0 x
let log_jump_prob x y = 0.0

let samples = 
  Mcmc.pt_mcmc_array ~nskip:11 100000 101 log_likelihood log_prior jump_proposal log_jump_prob 0.0

let _ = 
  let out = open_out ("pt_test_" ^ (string_of_float (Mcmc.pt_beta ())) ^ "_samples.dat") in 
    Read_write.write (fun x -> [| x |]) out samples;
    close_out out

let _ = 
  let evid = Mcmc.thermodynamic_integrate samples in 
    Printf.printf "Process %d: log evidence %g (true value %g)\n%!" (Mpi.comm_rank Mpi.comm_world) evid (-2.99573)
