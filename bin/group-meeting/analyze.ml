let _ = Random.self_init ()

let rand_between a b = 
  let x = Random.float 1.0 in
    a +. (b-.a)*.x

let samples = 
  let inp = open_in "samples.dat" in 
  let msamp = Read_write.read (fun x -> x) inp in 
    Array.map (fun {Mcmc.value = x} -> x) msamp

let log_like log_p mu sigma = 
  let sum = ref 0.0 in 
    for i = 0 to Array.length samples - 1 do 
      let s = samples.(i) and 
          msum = ref 0.0 in
        for i = 0 to 1 do 
          for j = 0 to 1 do 
            let di = s.(i) -. mu.(i) and 
                dj = s.(j) -. mu.(j) in
              msum := !msum +. di*.sigma.(i).(j)*.dj
          done
        done;
        sum := !sum +. (log_p !msum)
    done;
    !sum

let glog_like (mu,sigma) = log_like (fun x -> Stats.log_gaussian 0.0 1.0 x) mu sigma
let clog_like (mu,sigma) = log_like (fun x -> Stats.log_cauchy 0.0 1.0 x) mu sigma

let in_bounds low high x = x <= high && x >= low

let log_prior =
  let lp = 2.0*.(log 0.5) in
    fun (mu, sigma) -> 
      if in_bounds (-1.0) 1.0 mu.(0) &&
        in_bounds (-1.0) 1.0 mu.(1) &&
        in_bounds 1.0 2.0 sigma.(0).(0) && 
        in_bounds 1.0 2.0 sigma.(1).(1) &&
        in_bounds 0.0 1.0 sigma.(1).(0) then 
        lp
      else
        neg_infinity

let jump_propose (mu,sigma) = 
  let delta = 1.0 /. (sqrt (float_of_int (Array.length samples))) in
  let mdelta = ~-.delta in
  let new_sigma = Array.map (fun row -> Array.map (fun x -> x +. (rand_between mdelta delta)) row) sigma in 
    new_sigma.(1).(0) <- new_sigma.(0).(1);
    (Array.map (fun x -> x +. (rand_between mdelta delta)) mu,
     new_sigma)

let state_to_array (mu, sigma) = [|mu.(0); mu.(1); sigma.(0).(0); sigma.(1).(1); sigma.(0).(1)|]
let array_to_state : float array -> float array * float array array = function 
  | [|mu0; mu1; sigma00; sigma11; sigma01|] -> 
    ([|mu0; mu1|],
     [|[|sigma00; sigma01|];
       [|sigma01; sigma11|]|])
  | _ -> raise (Invalid_argument "array_to_state: bad array")
     
(* let gmcmc = Mcmc.mcmc_array ~nbin:10000 ~nskip:100 100000 glog_like log_prior jump_propose (fun _ _ -> 0.0) ([|0.0; 0.0|], [|[|1.0; 0.0|]; [|0.0; 1.0|]|]) *)
(* let cmcmc = Mcmc.mcmc_array ~nbin:10000 ~nskip:100 100000 clog_like log_prior jump_propose (fun _ _ -> 0.0) ([|0.0; 0.0|], [|[|1.0; 0.0|]; [|0.0; 1.0|]|]) *)

let gmcmc = 
  let inp = open_in "gaussian.mcmc" in
  let samps = Read_write.read array_to_state inp in 
    close_in inp;
    samps

let cmcmc = 
  let inp = open_in "cauchy.mcmc" in
  let samps = Read_write.read array_to_state inp in 
    close_in inp;
    samps

(* let _ =  *)
(*   let out = open_out "gaussian.mcmc" in  *)
(*     Read_write.write state_to_array out gmcmc; *)
(*     close_out out; *)
(*     let out = open_out "cauchy.mcmc" in  *)
(*       Read_write.write state_to_array out cmcmc; *)
(*       close_out out *)

module Ev = Evidence.Make(
  struct
    type params = float array * float array array
    let to_coords = state_to_array
  end)

module Interp = Interpolate_pdf.Make(
  struct
    type point = float array * float array array
    let coord = state_to_array
    let point = array_to_state
  end)
  
let rjmcmc = 
  let gpdf = Interp.make (Array.map (fun {Mcmc.value = v} -> v) gmcmc) [|-1.0; -1.0; 1.0; 1.0; 0.0|]
    [|1.0; 1.0; 2.0; 2.0; 1.0|] and 
      cpdf = Interp.make (Array.map (fun {Mcmc.value = v} -> v) cmcmc) [|-1.0; -1.0; 1.0; 1.0; 0.0|]
    [|1.0; 1.0; 2.0; 2.0; 1.0|] in
  Mcmc.rjmcmc_array ~nbin:10000 ~nskip:100 100000 
    (glog_like, clog_like)
    (log_prior, log_prior)
    (Interp.draw gpdf, Interp.draw cpdf)
    ((fun _ x -> Interp.jump_prob gpdf () a), (fun _ b -> Interp.jump_prob cpdf () b))
    ((fun _ -> Interp.draw gpdf), (fun _ -> Interp.draw cpdf))
    ((fun _ a -> Interp.jump_prob gpdf () a),
     (fun _ b -> Interp.jump_prob cpdf () b))
    (0.5, 0.5)
    (gmcmc.(0).Mcmc.value, cmcmc.(0).Mcmc.value)

let _ = 
  Printf.printf "Evidence for gaussian: %g\n" (Ev.evidence_direct gmcmc);
  Printf.printf "Evidence for cauchy: %g\n" (Ev.evidence_direct cmcmc);
  Printf.printf "Harmonic means: %g %g\n" (Ev.evidence_harmonic_mean gmcmc) (Ev.evidence_harmonic_mean cmcmc);
  Printf.printf "Reverse Jump ratio: %g\n" (Mcmc.rjmcmc_evidence_ratio rjmcmc)
