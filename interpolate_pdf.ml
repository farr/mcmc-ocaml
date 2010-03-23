(** Code that uses a kD tree to construct a piecewise-constant PDF
    from an array of samples (presumably from an MCMC). *)

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
  val make : point array -> float array -> float array -> interp_pdf

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

  module Kdt = Kd_tree.Make(
    struct
      type t = point

      let coord = S.coord
    end)

  type interp_pdf = point array * Kdt.tree

  let random_between a b = 
    a +. (b-.a)*.(Random.float 1.0)

  let random_in_volume low high = 
    let n = Array.length low in 
    let x = Array.make n 0.0 in 
      for i = 0 to n - 1 do 
        x.(i) <- random_between low.(i) high.(i)
      done;
      x

  let in_bounds pt low high = 
    let i = ref 0 and 
        n = Array.length pt in 
      while (!i < n) && (pt.(!i) >= low.(!i)) && (pt.(!i) <= high.(!i)) do 
        incr i
      done;
      !i = n

  let in_tree pt = function 
    | Kdt.Empty -> false
    | Kdt.Cell(_, low, high, _, _) -> 
        in_bounds pt low high

  let rec find_cell pt = function 
    | Kdt.Empty -> raise (Invalid_argument "find_cell: empty tree")
    | Kdt.Cell(_, _, _, Kdt.Empty, Kdt.Empty) as c -> 
        c
    | Kdt.Cell(_, _, _, left, right) -> 
        if in_tree pt left then 
          find_cell pt left
        else
          find_cell pt right

  let make pts low high = 
    (pts, Kdt.tree_of_objects (Array.to_list pts) low high)

  let draw (pts, tree) = 
    let pt = pts.(Random.int (Array.length pts)) in 
    match find_cell (S.coord pt) tree with 
      | Kdt.Empty -> raise (Failure "draw: empty tree")
      | Kdt.Cell(_,low,high,_,_) -> 
          S.point (random_in_volume low high)

  let jump_prob (pts, tree) _ pt = 
    let n = float_of_int (Array.length pts) in 
      match find_cell (S.coord pt) tree with 
        | Kdt.Empty -> raise (Failure "jump_prob: empty tree")
        | Kdt.Cell(objs, low, high, _, _) -> 
            let nobjs = float_of_int (List.length objs) and 
                v = Kdt.bounds_volume low high in 
              nobjs /. (v *. n)
end
