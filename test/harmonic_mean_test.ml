(* Simulate a gaussian, ask for the evidence. *)

module Ev = Evidence.Make(struct
  type params = float
  let to_coords (x : params) = [|x|]
end)

let _ = Random.self_init ()

let log_likelihood x = 
  let d = x*.x/.2.0 in 
    ~-.d -. 0.91893853320467274178

let log_prior x = if x < -1.0 || x > 1.0 then neg_infinity else -0.69314718055994530942

let jump_proposal x = x +. (Random.float 1.0 -. 0.5)

let log_jump_prob _ _ = 0.0

let samples = Mcmc.mcmc_array ~nbin:10000 ~nskip:100 100000 log_likelihood log_prior jump_proposal log_jump_prob 0.0

let compare_float (x : float) y = Pervasives.compare x y

let _ = 
  let bsamples = Array.copy samples in 
  let bevs = Array.make 10000 0.0 in
  let ev = Ev.evidence_harmonic_mean samples in 
    for i = 0 to 9999 do 
      for j = 0 to Array.length bsamples - 1 do 
        bsamples.(j) <- samples.(Random.int (Array.length samples))
      done;
      bevs.(i) <- Ev.evidence_harmonic_mean bsamples
    done;
    Array.fast_sort compare_float bevs;
    Printf.printf "%g %g %g\n" ev bevs.(1000) bevs.(9000) (* Should be 0.34134474606854294859 *)
