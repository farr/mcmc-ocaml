(*  read_write.mli: Read and write MCMC samples to/from files.
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

(** [read_sample from_coords chan] reads a single line from the given
    channel and interprets it as a single sample written with
    write_sample, returning the corresponding sample. *)
val read_sample : (float array -> 'a) -> in_channel -> 'a Mcmc.mcmc_sample

(** [read from_coords chan] reads from [chan] until an EOF is
    encountered a sequence of lines.  Each line contains some
    coordinates for a sample value ([from_coords] converts between
    coordinates and sample values), a log_likelihood, and a log_prior
    separated by non-newline whitespace.  Each line is converted to a
    sample, and the entire sample is returned as the result of the
    [read] call. *)
val read : (float array -> 'a) -> in_channel -> 'a Mcmc.mcmc_sample array
