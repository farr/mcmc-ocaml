(** MCMC computation. *)

(** Store probability information unboxed. *)
type like_prior = {
  log_likelihood : float;
  log_prior : float;
}

(** One MCMC sample. *)
type 'a mcmc_sample = {
  value : 'a;
  like_prior : like_prior;
}

(** [make_mcmc_sampler log_likelihood log_prior proposal log_jump_prob] *)
val make_mcmc_sampler : ('a -> float) -> ('a -> float) ->
  ('a -> 'a) -> ('a -> 'a -> float) -> 
  ('a mcmc_sample -> 'a mcmc_sample)

(** [mcmc_array n ll lp prop ljp x] *)
val mcmc_array : int -> ('a -> float) -> ('a -> float) -> 
  ('a -> 'a) -> ('a -> 'a -> float) -> 
  'a -> ('a mcmc_sample) array
