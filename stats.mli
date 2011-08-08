(*  stats.mli: Statistical functions.
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

(** Useful statistics functions. *)

(** [mean xs] returns the mean value. *)
val mean : float array -> float

(** [meanf f xs] returns the mean value of [f] over the [xs]. *)
val meanf : ('a -> float) -> 'a array -> float

(** [multi_mean xs] returns the mean value of the data in each
    dimension. *)
val multi_mean : float array array -> float array

(** [std ?mean xs] returns the standard deviation of the [xs].  In the
    common case that the mean is already known it can be supplied in
    the [?mean] optional argument to eliminate redundant
    computation. *)
val std : ?mean : float -> float array -> float

(** [stdf ?mean f xs] returns the standard deviation of [f] over the
    [xs].  As in {!Stats.std}, a previously-computed mean value can be
    suppled as the [~mean] argument to make the calculation more
    efficient. *)
val stdf : ?mean : float -> ('a -> float) -> 'a array -> float

(** [multi_std ?mean xs] returns the standard deviation of the
    multi-dimensional samples [xs].  It is completely equivalent to
    [Array.init dim (fun i -> std (Array.map (fun x -> x.(i)) xs))],
    but more efficient. *)
val multi_std : ?mean : float array -> float array array -> float array

(** [draw_cauchy x0 gamma] draws a random number from the Cauchy
    distribution with mode [x0] and HWHM [gamma]. *)
val draw_cauchy : float -> float -> float

(** [log_cauchy x0 gamma x] returns the log of the Cauchy PDF with
    mode [x0] and HWHM [gamma] evaluated at [x] *)
val log_cauchy : float -> float -> float -> float

(** [draw_gaussian mu sigma] draws a random number from the Gaussian
    (normal) distribution with mean [mu] and standard deviation
    [sigma]. *)
val draw_gaussian : float -> float -> float

(** [log_gaussian mu sigma x] returns the log of the Gaussian PDF with
    mean [mu] and standard deviation [sigma] at [x]. *)
val log_gaussian : float -> float -> float -> float

(** [log_multi_gaussian mu sigma x] returns the log of a
    multi-dimensional Gaussian PDF with mean [mu] and standard deviation
    [sigma]. *)
val log_multi_gaussian : float array -> float array -> float array -> float

(** [log_lognormal mu sigma x] returns the log of the log-normal PDF.
    Note that [mu] and [sigma] are the corresponding parameters in the
    log-normal distribution, and are not the mean and standard
    deviation. *)
val log_lognormal : float -> float -> float -> float

(** [draw_uniform min max] returns a uniformly distributed random
    number between [min] and [max]. *)
val draw_uniform : float -> float -> float
  
(** [find_nth ?copy n xs] returns the [n]th element of [xs] in
    ascending order.  The [?copy] parameter governs whether a copy of
    [xs] is made before the search; if not, then the array [xs] will
    be disordered on return.  The disordering actually makes
    subsequent searches more efficient (in the limit that each element
    of [xs] is requested, the algorithm performs a net quicksort,
    leaving the array in sorted order). *)
val find_nth : ?copy : bool -> int -> float array -> float

(** [find_nthf ?copy compare n xs] returns the [n]th element of [xs]
    ordered ascending in the given comparison function.  The [?copy]
    parameter governs whether a copy of [xs] is made before the
    search; if not, then the array [xs] will be disordered on return.
    The disordering actually makes subsequent searches more efficient
    (in the limit that every element of [xs] is requested, the
    algorithm performs a net quicksort, leaving the array in sorted
    order). *)
val find_nthf : ?copy : bool -> ('a -> 'a -> int) -> int -> 'a array -> 'a

(** [slow_autocorrelation nslides x] returns an array [nslides]
    elements long.  The element at index zero of this array is the
    variance of the series [x], and successive elements are the sum of
    products of elements of [x] with itself shifted in time.  The
    function is called [slow_autocorrelation] because it does not use
    fast (i.e. FFT) methods, but rather slow time-slide-and-sum
    methods to compute the autocorrelation.  The total computation
    time is O([nslides][n]), where [n] is the length of [x].*)
val slow_autocorrelation : int -> float array -> float array

(** [log_sum_logs logx logy] returns the log of the sum of [x] and [y]
    given [logx = log x] and [logy = log y] accurately. *)
val log_sum_logs : float -> float -> float
