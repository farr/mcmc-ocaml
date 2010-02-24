(** Code that uses a kD tree to construct a piecewise-constant PDF
    from a list of samples (presumable from an MCMC). *)

(** Abstract probability space that can map to and from R^n. *)
module type PROB_SPACE = sig
  (** Points in the space. *)
  type point 
    
  (** Coordinates of a point. *)
  val coord : point -> float array

  (** Construct a point from its coordinates. *)
  val point : float array -> point
end

(** Output signature for [Make]. *)
module type INTERPOLATE_PDF = sig
  (** Points in the abstract prob. space. *)
  type point 

  (** Abstract type for the interpolated PDF. *)
  type interp_pdf 

  (** Construct an interpolating PDF from the given list of samples
      (with the given bounds on the coordinates for the probability
      space). *)
  val make : point list -> float array -> float array -> interp_pdf

  (** Draw a point from the interpolated PDF. *)
  val draw : interp_pdf -> point

  (** [jump_prob interp_pdf source target] computes the probability
      that the output point of [draw] will be [target].  The [source]
      argument isn't used, but appears for compatibility with the jump
      probability functions expected by [Mcmc]. *)
  val jump_prob : interp_pdf -> 'a -> point -> float
end

module Make(S : PROB_SPACE) : INTERPOLATE_PDF with type point = S.point = struct
  type point = S.point 

  module Kdt = Kd_tree.Make_filling(
    struct
      type t = point

      let coord = S.coord
    end)

  (** Stores an array of the points and also the necessary tree structure. *)
  type interp_pdf = point array * Kdt.tree

  let random_between a b = 
    a +. (b-.a)*.(Random.float 1.0)

  let make samps low high = 
    (Array.of_list samps, 
     Kdt.tree_of_objects samps low high)

  let rec find_cell samp = function 
    | Kdt.Empty -> raise (Failure "find_cell: couldn't find")
    | Kdt.Cell(_,_,_,Kdt.Empty,Kdt.Empty) as c -> c
    | Kdt.Cell(_,low,high,left,right) -> 
        if Kdt.in_tree samp left then 
          find_cell samp left
        else
          find_cell samp right

  let draw (samps, tree) = 
    let samp = samps.(Random.int (Array.length samps)) in 
      match find_cell samp tree with 
        | Kdt.Empty -> raise (Failure "draw: internal error")
        | Kdt.Cell(_, low, high, _, _) -> 
            let coord = Array.copy low in 
              for i = 0 to Array.length low - 1 do 
                coord.(i) <- random_between low.(i) high.(i)
              done;
              S.point coord

  let pt_in_tree pt = function 
    | Kdt.Cell(_, low, high, _, _) when Kdt.in_bounds pt low high -> true
    | _ -> false

  let rec find_cell_pt pt = function 
    | Kdt.Empty -> raise (Failure "find_cell_pt: internal error")
    | Kdt.Cell(_, _, _, Kdt.Empty, Kdt.Empty) as c -> c
    | Kdt.Cell(_, _, _, left, right) -> 
        if pt_in_tree pt left then 
          find_cell_pt pt left
        else
          find_cell_pt pt right

  let jump_prob (samps, tree) _ pt = 
    let pt = S.coord pt in 
    let nf = float_of_int (Array.length samps) in 
      match find_cell_pt pt tree with 
        | Kdt.Empty -> raise (Failure "jump_prob: internal error")
        | Kdt.Cell(samps, _, _, _, _) as c -> 
            let v = Kdt.volume c and 
                ns = float_of_int (List.length samps) in 
              ns/.nf/.v

end
