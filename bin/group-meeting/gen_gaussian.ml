let _ = Random.self_init ()

let rand_float_between a b = 
  let x = Random.float 1.0 in 
    a +. (b -. a)*.x

let mu = Array.init 2 (fun _ -> rand_float_between (-1.0) 1.0)
let sigma = 
  let s = Array.make_matrix 2 2 0.0 in 
    s.(0).(0) <- rand_float_between 1.0 2.0;
    s.(1).(1) <- rand_float_between 1.0 2.0;
    s.(0).(1) <- rand_float_between 0.0 1.0;
    s.(1).(0) <- s.(0).(1);
    s

let _ = 
  let out = open_out "params.dat" in 
    Printf.fprintf out 
      "%g %g %g %g %g\n"
      mu.(0) mu.(1)
      sigma.(0).(0) sigma.(1).(1) sigma.(0).(1);
    close_out out

let log_like x = 
  let sum = ref 0.0 in 
    for i = 0 to 1 do 
      for j = 0 to 1 do 
        let di = x.(i) -. mu.(i) and 
            dj = x.(j) -. mu.(j) in 
          sum := !sum +. sigma.(i).(j)*.di*.dj
      done
    done;
    Stats.log_gaussian 0.0 1.0 !sum

let jump_propose = function 
  | [|x; y|] -> 
    [| x +. rand_float_between (-1.0) 1.0;
       y +. rand_float_between (-1.0) 1.0|]
  | _ -> raise (Invalid_argument "jump_propose: bad state")

let _ = 
  let samples = Mcmc.mcmc_array ~nskip:100 100 log_like (fun _ -> 0.0) jump_propose (fun _ _ -> 0.0)
    mu in 
  let out = open_out "samples.dat" in 
    Read_write.write (fun x -> x) out samples;
    close_out out
