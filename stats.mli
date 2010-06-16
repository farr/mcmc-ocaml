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
