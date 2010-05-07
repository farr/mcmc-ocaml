(** MCMC computation.  Uses the Metropolis algorithm to choose values
    in an arbitrary space distributed (asymptotically as the chain
    length goes to infinity) according to a given probability
    distribution.  The functions in this module use the native OCaml
    random number generator to generate the random numbers used in the
    probabilistic sampling from the chain. *)

(** We store the number of accepted and rejected steps and provide
    access through the following two functions.  [reset_counters ()]
    resets the counters to 0, and [get_counters ()] returns the
    counter values: number accepted steps and number rejected
    steps. *)
val reset_counters : unit -> unit
val get_counters : unit -> int * int

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

(** Values from a reverse-jump mcmc between two parameter spaces, of
    type ['a] and ['b] respectively. *)
type ('a, 'b) rjmcmc_value = 
  | A of 'a 
  | B of 'b

(** Samples from a reverse-jump mcmc between two parameter spaces, of
    type ['a] and ['b], respectively. *)
type ('a, 'b) rjmcmc_sample = ('a, 'b) rjmcmc_value mcmc_sample

(** [make_rjmcmc_sampler log_likelihoods log_priors
    internal_jump_proposals log_internal_jump_probabilities
    transition_jump_proposals log_transition_jump_probabilities
    model_priors] produces a reverse-jump MCMC sampling function.  The
    arguments to [make_rjmcmc_sampler] are as follows: 

    - [log_likelihoods] is a pair of functions that computes the log
      of the likelihood of parameters in each model.

    - [log_priors] is a pair of functions that computes the log of the
      priors in each model.

    - [internal_jump_proposals] is a pair of functions that propose
      jumps {b within} each model.

    - [log_internal_jump_probabilities] is a pair of functions that
      compute the probability of an internal jump between the first
      and second given state in each model.

    - [transition_jump_proposals] is a pair of functions that each
    produce a proposed jump into the model (from the opposite
    model).  For now, we assume that the models are totally separate
    (i.e. that one is not a superset of the other), so the proposal
    does not depend on the current state (in the other model).  This
    is the most likely property of [make_rjmcmc_sampler] to change.

    - [log_transition_jump_probabilities] is a pair of functions that
    return the log of the probability to propose a transition jump
    to the given state. 

    - [model_priors] are the priors on the models under consideration
    (note that these two numbers should sum to [1.0]).

    The [log_likelihoods], [log_priors], [internal_jump_proposals],
    and [log_internal_jump_probabilities] procedures are exactly the
    same as the corresponding functions for a single parameter-space
    MCMC. 

    The produced MCMC sample procedure will propose jumps between
    models with probability proportional to the corresponding model
    prior.
*)
val make_rjmcmc_sampler : 
  ('a -> float) * ('b -> float) -> 
  ('a -> float) * ('b -> float) -> 
  ('a -> 'a) * ('b -> 'b) -> 
  ('a -> 'a -> float) * ('b -> 'b -> float) -> 
  (unit -> 'a) * (unit -> 'b) -> 
  ('a -> float) * ('b -> float) -> 
  float * float -> 
  ('a, 'b) rjmcmc_sample -> ('a, 'b) rjmcmc_sample

(** [rjmcmc_array n log_likelihoods log_priors internal_jump_proposals
    log_internal_jump_probabilities transition_jump_proposals
    log_transition_jump_probabilities model_priors initial_states]
    produces an array of reverse-jump MCMC samples from the posterior
    of the two-parameter-space model described by its arguments,
    beginning with one of the [initial_states] (which one is chosen
    randomly according to the model priors).  See
    {!Mcmc.make_rjmcmc_sampler} for a description of the arguments. *)
val rjmcmc_array : 
  int -> 
  ('a -> float) * ('b -> float) -> 
  ('a -> float) * ('b -> float) -> 
  ('a -> 'a) * ('b -> 'b) -> 
  ('a -> 'a -> float) * ('b -> 'b -> float) -> 
  (unit -> 'a) * (unit -> 'b) -> 
  ('a -> float) * ('b -> float) -> 
  float * float -> 
  'a * 'b -> ('a, 'b) rjmcmc_sample array 

(** [rjmcmc_model_counts samples] computes the number of times
    parameters from the two models in a reverse-jump MCMC appear in
    the sequence of samples.  The first number returned is the number
    of times [A(_)] appears in [samples]; the second is the number of
    times [B(_)] appears. *)
val rjmcmc_model_counts : ('a, 'b) rjmcmc_sample array -> (int * int)

(** [combine_jump_proposals \[(p1, jump1, log_jump_prob1); ...\]]
    combines the given list of jump proposals into a single jump
    proposal and jump probability which chooses [jump1] with
    probability proportional to [p1], etc.  The probabilities need not
    be normalized. *)
val combine_jump_proposals : 
  (float * ('a -> 'a) * ('a -> 'a -> float)) list -> 
  ('a -> 'a) * ('a -> 'a -> float)

(** [uniform_wrapping xmin xmax dx x] returns a uniform random number
    within a range of size [dx] about the point [x], wrapping the
    value if it becomes smaller than [xmin] or greater than [xmax].
    ([dx] must be smaller than the range [xmax -. xmin].) *)
val uniform_wrapping : float -> float -> float -> float -> float

(** Type of two-model samples with memory.  The [Left] and [Right] tag
    records the "currently active" state of the model. *)
type ('a, 'b) memory = 
    private
  | Left of 'a * 'b
  | Right of 'a * 'b

(** [make_memory_rjmcmc_sampler (log_like_left, log_like_right)
    (log_prior_left, log_prior_right) (jump_proposal_left,
    jump_proposal_right) (log_jump_prob_left, log_jump_prob_right)
    (log_model_prior_left, log_model_prior_right)] constructs an MCMC
    sampler for comparing two models.  The sampler records the
    "currently active" state (either [Left] or [Right]); the other
    component of the state is a "memory" of the "current" location in
    the other model.  When proposing jumps between models, we always
    propose within one jump of the current locations in each model;
    when proposing jumps within a model, we fix the inactive position.
    In this way, we get both samples of the PDF within each model and,
    by comparing the ratio of counts with each model active, an
    estimate of the evidence ratio between the models. *)
val make_memory_rjmcmc_sampler : 
  ('a -> float) * ('b -> float) -> 
  ('a -> float) * ('b -> float) -> 
  ('a -> 'a) * ('b -> 'b) -> 
  ('a -> 'a -> float) * ('b -> 'b -> float) -> 
  float * float -> 
  ('a, 'b) memory mcmc_sample -> ('a, 'b) memory mcmc_sample

(** Like {!Mcmc.mcmc_array}, but for {!Mcmc.make_memory_rjmcmc_sampler}. *)
val memory_rjmcmc_array : 
  int -> 
  ('a -> float) * ('b -> float) -> 
  ('a -> float) * ('b -> float) -> 
  ('a -> 'a) * ('b -> 'b) -> 
  ('a -> 'a -> float) * ('b -> 'b -> float) -> 
  float * float -> 
  'a * 'b -> ('a, 'b) memory mcmc_sample array

(** [memory_evidence_ratio samples] computes an estimate of the
    evidence ratio between the models composing [samples]. *)
val memory_evidence_ratio : ('a, 'b) memory mcmc_sample array -> float

(** [split_memory_array (pa,pb) samples] returns two arrays from the
    memory rjmcmc samples in [samples].  Each array could, in
    principle, be the output of an individual MCMC in each of the
    models.  [pa] and [pb] are the priors for the two models; these
    are needed to recreate the correct priors {b within} the models
    from the memory rjmcmc priors. *)
val split_memory_array : float * float -> ('a, 'b) memory mcmc_sample array -> 
  'a mcmc_sample array * 'b mcmc_sample array
