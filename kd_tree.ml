(*  kd_tree.ml: kD-trees.
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

(** A k-D tree datastructure. *)

(** Input functor for the k-D structure. *)
module type COORDINATE_OBJECT = sig
  (** The objects. *)
  type t

    (** [coord obj] returns a float array of coordinates for the given
        object.  [coord] can be called many times in the construction
        of a k-D tree, so it should be as fast as possible. *)
  val coord : t -> float array
end

(** Output signature of the [Make] functor. *)
module type KD_TREE = sig
  (** Objects. *)
  type o 

  (** Trees.  They are either empty, or else contain a list of
      objects, bounds, and two sub-trees.  If the bounds are
      [(low,high)], then the objects' coordinates, [x], satisfy
      [low.(i) <= x.(i)] and [x.(i) <= high.(i)] along each dimension
      [i].  The left sub-tree contains objects whose coordinates are
      less than or equal to the median along the dimension of largest
      spread of the objects; the right sub-tree contains objects whose
      coordinates are greater than the median along this dimension. *)
  type tree = private 
              | Cell of o list * float array * float array * tree * tree
              | Empty

  (** Constructs a tree from the given objects and bounds. *)
  val tree_of_objects : o list -> float array -> float array -> tree

  (** Returns two bounds that enclose the given objects. *)
  val bounds_of_objects : o list -> float array * float array

  (** [bounds_volume low high] computes the volume enclosed by the
      given rectangular bounds. *)
  val bounds_volume : float array -> float array -> float

  (** [volume tree] returns the coordinate volume spanned by the
      tree. *)
  val volume : tree -> float
end

module Make(O : COORDINATE_OBJECT) : KD_TREE with type o = O.t = struct
  type o = O.t

  type tree = 
    | Cell of o list * float array * float array * tree * tree
    | Empty

  let rec find_ith compare i = function 
    | [] -> raise (Invalid_argument "find_ith: no objects")
    | [x] -> 
        if i = 0 then 
          x
        else
          raise (Failure "find_ith: internal error")
    | x :: xs as objs -> 
        if List.for_all (fun y -> compare x y = 0) xs then 
          x
        else
          let pvt = List.nth objs (Random.int (List.length objs)) in 
          let (lte, gt) = List.partition (fun obj -> compare obj pvt <= 0) objs in 
          let n_lte = List.length lte in 
            if i < n_lte then 
              find_ith compare i lte
            else
              find_ith compare (i-n_lte) gt

  let compare_along_dim i o1 o2 = 
    let c1 = O.coord o1 and 
        c2 = O.coord o2 in 
      Pervasives.compare c1.(i) c2.(i)

  let compare_coords o1 o2 = 
    Pervasives.compare (O.coord o1) (O.coord o2)

  let bounds_of_objects = function 
    | [] -> raise (Invalid_argument "bounds_of_objects: no objects")
    | [x] -> (O.coord x, O.coord x)
    | x :: xs -> 
        let low = Array.copy (O.coord x) in 
        let high = Array.copy low in 
          List.iter 
            (fun x -> 
               let c = O.coord x in 
                 for i = 0 to Array.length low - 1 do 
                   if c.(i) < low.(i) then low.(i) <- c.(i);
                   if c.(i) > high.(i) then high.(i) <- c.(i)
                 done)
            xs;
          (low,high)

  let split_bounds low high olow ohigh dim = 
    let x = 0.5*.((O.coord olow).(dim) +. (O.coord ohigh).(dim)) in 
    let new_low = Array.copy low and 
        new_high = Array.copy high in 
      new_low.(dim) <- x;
      new_high.(dim) <- x;
      (new_low, new_high)

  let longest_dim low high = 
    let dim = ref (-1) and 
        dx_max = ref neg_infinity in 
      for i = 0 to Array.length low - 1 do 
        let dx = high.(i) -. low.(i) in 
          if dx > !dx_max then begin
            dim := i;
            dx_max := dx
          end
      done;
      !dim

  let find_max comp = function 
    | [] -> raise (Invalid_argument "find_max: no objects")
    | [x] -> x
    | x :: xs -> 
        List.fold_left (fun max x -> if comp x max > 0 then x else max) x xs

  let find_min comp = function 
    | [] -> raise (Invalid_argument "find_min: no objects")
    | [x] -> x
    | x :: xs -> 
        List.fold_left (fun min x -> if comp x min < 0 then x else min) x xs

  let adjust_for_empty_split comp objs objs2 =
    match objs, objs2 with 
      | [], [] -> [], []
      | [], objs -> 
          let min = find_min comp objs in 
            List.partition (fun obj -> comp obj min <= 0) objs
      | objs, [] -> 
          let max = find_max comp objs in 
            List.partition (fun obj -> comp obj max < 0) objs 
      | objs, objs2 -> objs, objs2

  let rec tree_of_objects objs low high = 
    match objs with 
      | [] -> Empty
      | [x] as l -> Cell(l, low, high, Empty, Empty)
      | (x :: xs as l) when List.for_all (fun y -> compare_coords x y = 0) xs -> 
          Cell(l, low, high, Empty, Empty)
      | objs -> 
          let n = List.length objs in 
          let i = n / 2 in 
          let (l,h) = bounds_of_objects objs in 
          let dim = longest_dim l h in 
          let comp x y = compare_along_dim dim x y in 
          let pvt = find_ith comp i objs in 
          let (lte, gt) = List.partition (fun x -> compare_along_dim dim x pvt <= 0) objs in 
          let (lte, gt) = adjust_for_empty_split comp lte gt in 
          let lt_bound = find_max comp lte and 
              gt_bound = find_min comp gt in 
          let (new_low, new_high) = split_bounds low high lt_bound gt_bound dim in 
            Cell(objs, low, high, 
                 tree_of_objects lte low new_high,
                 tree_of_objects gt new_low high)

  let bounds_volume low high = 
    let v = ref 1.0 in 
      for i = 0 to Array.length low - 1 do 
        v := !v *. (high.(i) -. low.(i))
      done;
      !v +. 0.0
    
  let volume = function 
    | Empty -> 0.0
    | Cell(_, low, high, _, _) -> bounds_volume low high
end
