(** Code for reading and writing MCMC output to channel in
    human-readable form.  An MCMC is written using a function that
    converts from the arbitrary MCMC values to arrays of floats
    describing coordinates for the value.  The output is the
    coordinates separated by spaces, followed by the log_likelihood
    and log_prior for that value, and then a newline.  Corresponding
    input can be read into an array of MCMC results using a function
    that converts from coordinates to MCMC values. *)

(** [write_sample to_coords chan sample] outputs the sample on a
    single line of [chan].  The value of [sample] is converted to a
    [float array] using [to_coords]; this float array is written to
    [chan], followed by the log_likelihood and log_prior of the
    sample, all separated by spaces. *)
val write_sample : ('a -> float array) -> out_channel -> 'a Mcmc.mcmc_sample -> unit

(** [write to_coords chan mcmc] writes the samples in the [mcmc] array
    to [chan] using [to_coords] to convert sample values to float
    coordinates. *)
val write : ('a -> float array) -> out_channel -> ('a Mcmc.mcmc_sample array) -> unit

(** [read from_coords chan] reads from [chan] until an EOF is
    encountered a sequence of lines.  Each line contains some
    coordinates for a sample value ([from_coords] converts between
    coordinates and sample values), a log_likelihood, and a log_prior
    separated by non-newline whitespace.  Each line is converted to a
    sample, and the entire sample is returned as the result of the
    [read] call. *)
val read : (float array -> 'a) -> in_channel -> 'a Mcmc.mcmc_sample array
