(** Extendable vectors.  All single-element operations complete in
    amoritzed constant-time. *)

(** Extendable vectors. *)
type 'a earray

(** Retrieve values. *)
val get : 'a earray -> int -> 'a

(** Set values. *)
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
