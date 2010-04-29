(** Various algorithms for computing the Bayesian evidence from an
    MCMC sample. *)

(** Input type for the [Make] functor. *)
module type MCMC_OUT = sig
  (** Parameters in the MCMC. *)
  type params

  (** [to_coords params] returns a float array representing the
      parameters. *)
  val to_coords : params -> float array
end

(** Output type for the [Make] functor. *)
module type EVIDENCE = sig
  (** Parameters. *)
  type params 

  (** MCMC samples. *)
  type sample = params Mcmc.mcmc_sample

  (** Embedded kD tree. *)
  module Kd : (Kd_tree.KD_TREE with type o = sample)

  (** Construct a kD tree. *)
  val kd_tree_of_samples : sample array -> float array -> float array -> Kd.tree

  (** Directly integrate evidence within the rectangular region in
      parameter space between the two given arrays. *)
  val evidence_direct : ?n : int -> sample array -> float

  (** Integrate evidence using harmonic mean. *)
  val evidence_harmonic_mean : sample array -> float

  (** Integrate evidence using Lebesque integral of 1/L, within the
  rectangular parameter region between the given arrays.  *)
  val evidence_lebesgue : ?n : int -> ?eps : float -> sample array -> float
end

module Make(MO : MCMC_OUT) : EVIDENCE with type params = MO.params = struct
  type params = MO.params
  type sample = params Mcmc.mcmc_sample

  module Kd = Kd_tree.Make(
    struct
      type t = sample
      let coord ({Mcmc.value = v} : t) = MO.to_coords v
    end)

  let kd_tree_of_samples samps = Kd.tree_of_objects (Array.to_list samps)

  let rec length_at_least n = function 
    | [] -> n = 0
    | _ :: xs -> 
        (n <= 0) || (length_at_least (n-1) xs)

  let rec depth = function 
    | Kd.Empty -> 0
    | Kd.Cell(_, _, _, left, right) -> 
        1 + (max (depth left)
               (depth right))

  let rec collect_subvolumes nmax = function 
    | Kd.Empty -> []
    | Kd.Cell(objs, _, _, left, right) as c -> 
        if not (length_at_least nmax objs) then 
          [c]
        else
          List.rev_append (collect_subvolumes nmax left) (collect_subvolumes nmax right)

  let log_likelihood {Mcmc.like_prior = {Mcmc.log_likelihood = ll}} = ll

  let log_prior {Mcmc.like_prior = {Mcmc.log_prior = lp}} = lp

  let log_posterior s = (log_likelihood s) +. (log_prior s)

  let compare_inverse_like s1 s2 = Pervasives.compare (~-.(log_likelihood s1)) (~-.(log_likelihood s2))

  let evidence_harmonic_mean samples = 
    let linv = ref 0.0 and 
        n = Array.length samples in 
      for i = 0 to n - 1 do 
        linv := !linv +. 1.0/.(exp (log_likelihood samples.(i)))
      done;
      (float_of_int n)/.(!linv)
          
  let median_sample (f : sample -> float) samples = 
    let n = List.length samples and 
        ssamp = 
      List.fast_sort 
        (fun s1 s2 -> 
           Pervasives.compare (f s1) (f s2))
        samples in 
      if n mod 2 = 0 then 
        0.5*.((f (List.nth ssamp (n/2 - 1))) +.
                (f (List.nth ssamp (n/2))))
      else 
        f (List.nth ssamp (n/2))
      
  let compare_samples s1 s2 = 
    let {Mcmc.value = v1} = s1 and 
        {Mcmc.value = v2} = s2 in 
      Pervasives.compare (MO.to_coords v1) (MO.to_coords v2)

  let rev_remove_dups comp l = 
    let rec rev_remove_dups_loop removed remaining = 
      match remaining with 
        | [] -> removed
        | [x] -> x :: removed
        | x :: (y :: _ as ys) -> 
            if comp x y = 0 then 
              rev_remove_dups_loop removed ys
            else
              rev_remove_dups_loop (x :: removed) ys in
      rev_remove_dups_loop [] l

  let array_to_list_remove_dups samples = 
    let l = Array.to_list samples in 
    let lsort = List.sort compare_samples l in 
      rev_remove_dups compare_samples lsort

  let evidence_direct ?(n = 64) samples = 
    let lsamples = array_to_list_remove_dups samples in
    let (low,high) = Kd.bounds_of_objects lsamples in 
    let sub_vs = collect_subvolumes n (Kd.tree_of_objects lsamples low high) in 
      List.fold_left
        (fun integral c -> 
           match c with 
             | Kd.Cell(objs, _, _, _, _) ->
                 let (low,high) = Kd.bounds_of_objects objs in 
                 let vol = Kd.bounds_volume low high and 
                     lp = median_sample log_posterior objs in 
                   integral +. vol*.(exp lp)
             | _ -> raise (Invalid_argument "evidence_direct: bad cell in integral accumulation"))
        0.0
        sub_vs

  let collect_samples_up_to_eps eps samps = 
    let rec collect_samples_loop collected_samples = function 
      | [] -> List.rev collected_samples
      | [x] -> List.rev (x :: collected_samples)
      | x :: ((y :: _) as ys) -> 
          let ilx = exp (~-.(log_likelihood x)) and 
              ily = exp (~-.(log_likelihood y)) in 
          let delta = ily -. ilx in 
            assert(delta >= 0.0);
            if delta > eps then 
              List.rev (x :: collected_samples)
            else
              collect_samples_loop (x :: collected_samples) ys in 
      collect_samples_loop [] samps

  let evidence_lebesgue ?(n = 64) ?(eps = 0.1) samples = 
    let samples_list = array_to_list_remove_dups samples in
    let samps = collect_samples_up_to_eps eps (List.fast_sort compare_inverse_like samples_list) in 
    let (low,high) = Kd.bounds_of_objects samps in 
    let sub_vs = collect_subvolumes n (kd_tree_of_samples samples low high) in 
    let prior_mass = 
      List.fold_left
        (fun pm c -> 
           match c with 
             | Kd.Cell(objs, _, _, _, _) ->
                 let (low,high) = Kd.bounds_of_objects objs in 
                 let vol = Kd.bounds_volume low high and 
                     lp = median_sample log_prior objs in 
                   pm +. vol*.(exp lp)
             | _ -> raise (Failure "evidence_lebesque: bad cell in integral accumulation of prior mass"))
        0.0
        sub_vs and 
        inv_l = 
      List.fold_left
        (fun inv_l samp -> inv_l +. (exp (~-.(log_likelihood samp))))
        0.0 
        samps and 
        n = float_of_int (List.length samps) in 
      n*.prior_mass/.inv_l    
end
