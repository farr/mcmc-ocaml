(*  earray.mli: Extendable arrays.
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

(** Extendable vectors.  All single-element operations complete in
    amoritzed constant-time. *)

(** Extendable vectors. *)
type 'a earray

(** Retrieve values. *)
val get : 'a earray -> int -> 'a

(** Set values.  The index is allowed to point one past the end of the
    array; in this special case, [set ea i x] is equivalent to [append ea
    x]. *)
val set : 'a earray -> int -> 'a -> unit

(** Convert to regular array.  The exten*)
val to_array : 'a earray -> 'a array

(** Construct from a regular array.  The extendable array and the
    array argument do not share storage. *)
val of_array : 'a array -> 'a earray

(** Add an element to the end of the extendable vector.  (Completes in
    amortized constant time.) *)
val append : 'a earray -> 'a -> unit

(** Cut the array down to the given length, discarding any elements
    that may have been placed in slots beyond the new length. *)
val chop : 'a earray -> int -> unit

(** Construct an extendable array of the given size, with all slots
    filled with the initialization argument. *)
val make : int -> 'a -> 'a earray

(** Current number of elements in earray. *)
val length : 'a earray -> int

(** Copy. *)
val copy : 'a earray -> 'a earray
