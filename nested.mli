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
    draw_prior log_likelihood log_prior] returns a tuple [(evidence,
    delta_evidence, sampled_points)], where [evidence] is an estimate
    of the evidence integral for the model with [log_likelihood] and
    [log_prior], [delta_evidence] is an estimate of the error in the
    evidence estimate, and [sampled_points] is a list of the points
    sampled in the nested sampling algorithm.  

    The algorithm continues sampling from parameter space until it
    estimates that the (fractional) contribution to the evidence
    integral from the remaining live points is smaller than [epsrel].
    Each new live point is drawn from an MCMC using [nmcmc] samples
    starting from a random live point and using the live points to
    construct a differential-evolution proposal.  There are [nlive]
    live points in total.  [mode_hopping_frac] controls the fraction
    of DE jumps that are proposed in mode-hopping mode (see
    {!Mcmc.differential_evolution_proposal}).  [draw_prior] is used to
    produce the initial distribution of live points.
*)
val nested_evidence : 
  ?epsrel : float -> 
  ?nmcmc : int -> 
  ?nlive : int -> 
  ?mode_hopping_frac : float -> 
  (unit -> float array) -> 
  (float array -> float) -> 
  (float array -> float) -> 
  float * float * (float array Mcmc.mcmc_sample list)
