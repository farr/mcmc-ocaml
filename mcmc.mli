(** MCMC computation.  Uses the Metropolis algorithm to choose values
    in an arbitrary space distributed (asymptotically as the chain
    length goes to infinity) according to a given probability
    distribution.  The functions in this module use the native OCaml
    random number generator to generate the random numbers used in the
    probabilistic sampling from the chain. *)

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

(** [make_mcmc_sampler log_likelihood log_prior jump_proposal
    log_jump_prob] returns a function of a single sample that returns
    the next sample in the MCMC chain described by the given functions
    on a parameter space.  The Metropolis algorithm is used to
    generate the samples; this algorithm only depends on the
    probability ratios (difference of logs) of the posterior and the
    jump probability, so the probability distributions given to
    [make_mcmc_sampler] need not be normalized.  Additionally, the
    algorithm only depends on the ratio of the forward and backward
    jump probabilities: [(log_jump_prob source target)
    -. (log_jump_prob target source)].  Therefore, if your jump
    proposal is symmetric, i.e. [(log_jump_prob a b) = (log_jump_prob
    b a)] for all [a] and [b], then it need not be computed at all
    because it will cancel in the sampling. *)
val make_mcmc_sampler : ('a -> float) -> ('a -> float) ->
  ('a -> 'a) -> ('a -> 'a -> float) -> 
  ('a mcmc_sample -> 'a mcmc_sample)

(** [mcmc_array n_samples log_likelihood log_prior jump_proposal
    log_jump_prob start] construts an array of samples of length
    [n_samples] from the MCMC chain described by the given parameters
    (see {!Mcmc.make_mcmc_sampler}).  The first sample in the array is
    always [start]. *)
val mcmc_array : int -> ('a -> float) -> ('a -> float) -> 
  ('a -> 'a) -> ('a -> 'a -> float) -> 
  'a -> ('a mcmc_sample) array

(** Remove repeated samples from an array of MCMC samples.  Note that
    this {b changes the distribution of the samples}---the repeated
    samples are required for the probability distribution of the chain
    to be correct.  However, it can be useful for, e.g. integration
    (see {!Evidence.EVIDENCE.evidence_direct} and friends) to remove
    repeated samples from the chain. *)
val remove_repeat_samples : ('a -> 'a -> bool) -> ('a mcmc_sample) array -> 
  ('a mcmc_sample) array
