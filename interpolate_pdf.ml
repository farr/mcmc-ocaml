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

  module Kdt = Kd_tree.Make_filling(
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

  let center_pt low high = 
    let n = Array.length low in 
    let x = Array.make n 0.0 in 
      for i = 0 to n - 1 do 
        x.(i) <- 0.5*.(low.(i) +. high.(i))
      done;
      x

  let rec collect_extra_points = function 
    | Kdt.Null -> []
    | Kdt.Empty(low, high) -> [S.point (center_pt low high)]
    | Kdt.Cell(_,_,_,left,right) -> 
        List.rev_append (collect_extra_points left) (collect_extra_points right)

  let make pts low high = 
    let lpts = Array.to_list pts in 
    let tree = Kdt.tree_of_objects lpts low high in 
      (Array.of_list (List.rev_append (collect_extra_points tree) lpts),
       tree)

  let rec find_tree pt = function 
    | Kdt.Null -> raise (Failure "find_volume: internal error")
    | Kdt.Empty(_, _) as t -> t
    | Kdt.Cell(_,_,_,Kdt.Null,Kdt.Null) as t -> 
        t
    | Kdt.Cell(_,_,_,left,right) -> 
        if Kdt.in_tree_volume pt left then 
          find_tree pt left
        else
          find_tree pt right

  let find_volume pt tree = 
    match find_tree pt tree with 
      | Kdt.Empty(low,high) -> (low,high)
      | Kdt.Cell(_,low,high,_,_) -> (low,high)
      | _ -> raise (Failure "find_volume: internal error")

  let draw (pts, tree) = 
    let n = Array.length pts in 
    let pt = pts.(Random.int n) in 
    let (low,high) = find_volume pt tree in 
      S.point (random_in_volume low high)

  let jump_prob (pts, tree) _ pt = 
    let n = Array.length pts in 
    let nf = float_of_int n in 
      match find_tree pt tree with 
        | Kdt.Empty(_,_) as t -> 1.0 /. nf /. (Kdt.volume t)
        | Kdt.Cell(objs, _, _, _, _) as t -> 
            (float_of_int (List.length objs)) /. nf /. (Kdt.volume t)
        | _ -> raise (Failure "jump_prob: internal error")

end
