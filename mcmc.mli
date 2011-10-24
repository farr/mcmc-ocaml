(*  mcmc.mli: MCMC samplers and utilities.
    Copyright (C) 2011 Will M. Farr <w-farr@northwestern.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. *)

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

(** [mcmc_array ?nbin ?nskip n_samples log_likelihood log_prior
    jump_proposal log_jump_prob start] construts an array of samples
    of length [n_samples] from the MCMC chain described by the given
    parameters (see {!Mcmc.make_mcmc_sampler}).  [nbin] samples are
    initially discarded as burn-in.  If [nskip] is provided, it is the
    number of samples to produce before recording a sample in the
    array ([nskip] defaults to 1, which causes every state produced by
    the sampler to be recorded). *)
val mcmc_array : ?nbin : int -> ?nskip : int -> int -> ('a -> float) -> ('a -> float) -> 
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
    model).  

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
  ('b -> 'a) * ('a -> 'b) -> 
  ('b -> 'a -> float) * ('a -> 'b -> float) -> 
  float * float -> 
  ('a, 'b) rjmcmc_sample -> ('a, 'b) rjmcmc_sample

(** [rjmcmc_array ?nbin ?nskip n log_likelihoods log_priors
    internal_jump_proposals log_internal_jump_probabilities
    transition_jump_proposals log_transition_jump_probabilities
    model_priors initial_states] produces an array of reverse-jump
    MCMC samples from the posterior of the two-parameter-space model
    described by its arguments, beginning with one of the
    [initial_states] (which one is chosen randomly according to the
    model priors).  See {!Mcmc.make_rjmcmc_sampler} for a description
    of the arguments. *)
val rjmcmc_array : 
  ?nbin : int -> 
  ?nskip : int -> 
  int -> 
  ('a -> float) * ('b -> float) -> 
  ('a -> float) * ('b -> float) -> 
  ('a -> 'a) * ('b -> 'b) -> 
  ('a -> 'a -> float) * ('b -> 'b -> float) -> 
  ('b -> 'a) * ('a -> 'b) -> 
  ('b -> 'a -> float) * ('a -> 'b -> float) -> 
  float * float -> 
  'a * 'b -> ('a, 'b) rjmcmc_sample array 

(** [rjmcmc_model_counts samples] computes the number of times
    parameters from the two models in a reverse-jump MCMC appear in
    the sequence of samples.  The first number returned is the number
    of times [A(_)] appears in [samples]; the second is the number of
    times [B(_)] appears. *)
val rjmcmc_model_counts : ('a, 'b) rjmcmc_sample array -> (int * int)

(** [rjmcmc_evidence_ratio samples] returns the evidence ratio for the
    reversible-jump MCMC [samples]. *)
val rjmcmc_evidence_ratio : ('a, 'b) rjmcmc_sample array -> float

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

(** [differential_evolution_proposal ?mode_hopping_frac to_float
    from_float samples] returns a valid jump proposal that will use
    the given [samples] to construct a differential-evolution
    proposal.  The optional [to_float] and [from_float] arguments are
    used to convert the mcmc samples into and out of float arrays.
    The [mode_hopping_frac] argument controls how often the proposal
    operates in "mode-hopping" mode compared to how often it operates
    in normal mode (see below).

    A differential evolution proposal chooses two of the points from
    [samples], computes the vector in coordinate space from one to the
    other, and then uses this vector to propose a coordinate increment
    to the current point.  In the standard mode, the proposal jumps
    along this vector with a randomly-drawn magnitude distributed
    uniformly between 0 and 2 times the vector's magnitude; in
    mode-hopping mode, the proposal uses exacty the difference vector
    for the proposed coordinate increment.  By default the proposal
    operates in standard mode (i.e. [mode_hopping_frac] is set to
    [0.0]).

    The proposals generated by [differential_evolution_proposal] are
    symmetric, but it is difficult to calculate their exact jump
    probability.  Therefore, they should not be used as part of a
    combined jump proposal unless all parts are symmetric.
*)
val differential_evolution_proposal : 
  ?mode_hopping_frac : float -> ('a -> float array) -> (float array -> 'a) ->
  'a mcmc_sample array -> 
  'a -> 'a
