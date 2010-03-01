(** Useful statistics functions. *)

(** [mean xs] returns the mean value. *)
val mean : float array -> float

(** [std ?mean xs] returns the standard deviation of the xs.  In the
    common case that the mean is already known it can be supplied in
    the [?mean] optional argument to eliminate redundant
    computation. *)
val std : ?mean : float -> float array -> float

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
