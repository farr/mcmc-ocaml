(*  nested.mli: Interface for nested sampling routines.
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

(** Nested sampling. *)

(** [nested_evidence ?epsrel ?nmcmc ?nlive ?mode_hopping_frac
    draw_prior log_likelihood log_prior] returns a tuple
    [(log_evidence, log_delta_evidence, sampled_points,
    log_point_weights)], where [log_evidence] is an estimate of the
    log of the evidence integral for the model with [log_likelihood]
    and [log_prior], [log_delta_evidence] is the log of an estimate of
    the error in the evidence estimate, [sampled_points] is a list of
    the points sampled in the nested sampling algorithm, and
    [log_weights] is an array of the log of the weight of the
    corresponding sample points in the evidence integral.  This last
    quantity is useful for computing statistics from the sample
    points, as each point should contribute proportional to its
    weight.

    The algorithm continues sampling from parameter space until it
    estimates that the (fractional) contribution to the evidence
    integral from the remaining live points is smaller than [epsrel].
    Each new live point is drawn from an MCMC using [nmcmc] samples
    starting from a random live point and using the live points to
    construct a differential-evolution proposal.  There are [nlive]
    live points in total.  [mode_hopping_frac] controls the fraction
    of DE jumps that are proposed in mode-hopping mode (see
    {!Mcmc.differential_evolution_proposal}).  [draw_prior] is used to
    produce the initial distribution of live points.  The [observer]
    argument is called each time a live point is retired with the
    retired point as argument.
*)
val nested_evidence : 
  ?observer : (float array Mcmc.mcmc_sample -> unit) ->
  ?epsrel : float -> 
  ?nmcmc : int -> 
  ?nlive : int -> 
  ?mode_hopping_frac : float -> 
  (unit -> float array) -> 
  (float array -> float) -> 
  (float array -> float) -> 
  float * float * (float array Mcmc.mcmc_sample array) * (float array)

(** [total_error_estimate log_evidence log_delta_evidence nlive]
    returns the quadrature-sum of the systematic and statistical error
    for the evidence calculated by nested sampling.  [log_evidence]
    and [log_delta_evidence] are returned by the nested sampling
    algorithm, while the statistical error is estimated as a relative
    error of order 1/sqrt([nlive]).  *)
val log_total_error_estimate : float -> float -> int -> float
