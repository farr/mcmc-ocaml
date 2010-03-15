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

  (** Constructs a tree from the given objects. *)
  val tree_of_objects : o list -> tree

  (** [volume tree] returns the coordinate volume spanned by the
      objects in the tree. *)
  val volume : tree -> float
end

(** Output signature of the [Make_filling] functor.  The difference
    between space-filling kD trees and ordinary kD trees is that the
    filling kD trees (and their sub-cells) have bounds that exactly
    partition a rectangular region of R^n; regular kD trees have
    bounds that are as small as possible to enclose the objects
    contained in them.  At every level in a space-filling kD tree, the
    union of the cell volumes is equal to the total volume at the top
    of the tree; a normal kD tree will be missing some volume at lower
    levels because it shrinks the cell bounds down to the bounding box
    for the objects in the cell. *)
module type FILLING_KD_TREE = sig
  (** Objects. *)
  type o

  (** Trees. *)
  type tree = private 
              | Cell of o list * float array * float array * tree * tree
              | Empty of float array * float array (** Empty cell, but finite volume. *)
              | Null (** Empty cell of zero volume. *)

(** Construct a tree with the given bounds.  All objects should be
    in the range [\[low, high)].*)
  val tree_of_objects : o list -> float array -> float array -> tree

  (** Tree volume. *)
  val volume : tree -> float

  (** Is the given object contained within the volume of the tree? *)
  val in_tree_volume : o -> tree -> bool

  (** [in_bounds x low high] tests whether [x] is in the range
      \[[low], [high]). *)
  val in_bounds : float array -> float array -> float array -> bool

(** [integrate f tree] integrates the function [f] on the objects in
    the tree over the volume discretized by [tree].  *)
  val integrate : (o -> float) -> tree -> float

end

module Make(O : COORDINATE_OBJECT) : KD_TREE with type o = O.t = struct
  type o = O.t

  type tree = | Cell of o list * float array * float array * tree * tree
              | Empty

  let bounds_of_objects = function
    | [] -> raise (Invalid_argument "bounds_of_objects: need at least one object")
    | o :: objs -> 
        let low = Array.copy (O.coord o) in 
        let high = Array.copy low in 
          List.fold_left
            (fun (low, high) o ->
               let coord = O.coord o in 
                 for i = 0 to Array.length coord - 1 do 
                   if low.(i) > coord.(i) then low.(i) <- coord.(i);
                   if high.(i) < coord.(i) then high.(i) <- coord.(i)
                 done;
                 (low, high))
            (low, high)
            objs

  let compare_all_coords o1 o2 = 
    let c1 = O.coord o1 and 
        c2 = O.coord o2 in 
      Pervasives.compare c1 c2

  let compare_along_dimension dim o1 o2 = 
    let c1 : float array = O.coord o1 and 
        c2 = O.coord o2 in 
      Pervasives.compare c1.(dim) c2.(dim)

  let lte compare a b = 
    (compare a b) <= 0

  let rec all_equal compare = function 
    | [] -> true
    | [x] -> true
    | x :: (y :: _ as rest) -> compare x y = 0 && all_equal compare rest

  let rec find_nth_sorted compare objs n = 
    match objs with 
      | [] -> raise (Failure "find_nth_sorted: no objects")
      | [obj] -> 
          if n = 0 then obj else raise (Failure "find_nth_sorted: internal error, n > 0")
      | (o :: objs as l) when all_equal compare l -> o
      | objs -> 
          let len = List.length objs in 
          let pvt = List.nth objs (Random.int len) in 
          let (lte, gt) = List.partition (fun obj -> lte compare obj pvt) objs in 
          let len_lte = List.length lte in 
            if n < len_lte then 
              find_nth_sorted compare lte n
            else
              find_nth_sorted compare gt (n-len_lte)

  let partition_elt compare = function 
    | [] -> raise (Invalid_argument "median: no objects")
    | objs -> 
        let n = List.length objs in 
        let i = if n mod 2 = 0 then n/2 - 1 else n/2 in 
          find_nth_sorted compare objs i

  let longest_dimension low high = 
    let dx_max = ref (-1.0/.0.0) and 
        dim = ref (-1) in 
      for i = 0 to Array.length low - 1 do 
        let dx = high.(i) -. low.(i) in 
          if dx > !dx_max then begin
            dx_max := dx;
            dim := i
          end
      done;
      !dim

  let rec tree_of_objects = function 
    | [] -> Empty
    | (o :: _) as objs when all_equal compare_all_coords objs -> 
        let low = Array.copy (O.coord o) in 
          Cell(objs, low, low, Empty, Empty)
    | objs -> 
        let (low,high) = bounds_of_objects objs in 
        let dim = longest_dimension low high in 
        let compare = compare_along_dimension dim in 
        let part_elt = partition_elt compare objs in 
        let (lte, gt) = List.partition (fun obj -> lte compare obj part_elt) objs in 
          Cell(objs, low, high, 
               tree_of_objects lte,
               tree_of_objects gt)

  let volume = function 
    | Empty -> 0.0
    | Cell(_, low, high, _, _) -> 
        let vol = ref 1.0 in 
          for i = 0 to Array.length low - 1 do 
            vol := !vol *. (high.(i) -. low.(i))
          done;
          !vol
end

(** [Make_filling] is like [Make] except that the trees that it
    produces fill a rectangular coordinate volume completely instead
    of having bounds that exactly enclose a list of objects. *)
module Make_filling(O : COORDINATE_OBJECT) : FILLING_KD_TREE with type o = O.t = struct
  type o = O.t

  type tree = 
    | Cell of o list * float array * float array * tree * tree
    | Empty of float array * float array
    | Null

  let sub_bounds dim low high = 
    let new_low = Array.copy low and 
        new_high = Array.copy high in 
    let x_split = (low.(dim) +. high.(dim))*.0.5 in 
      new_low.(dim) <- x_split;
      new_high.(dim) <- x_split;
      (new_low, new_high)

  let in_bounds x low high = 
    let n = Array.length x in 
    let i = ref 0 and 
        inb = ref true in 
      while (!inb && !i <= n-1) do 
        let x = x.(!i) in 
          if x < low.(!i) || x >= high.(!i) then inb := false;
          incr i
      done;
      !inb

  let object_in_bounds o low high = in_bounds (O.coord o) low high

  let in_tree_volume o = function 
    | Cell(_, low, high, _, _) -> 
        object_in_bounds o low high
    | Empty(low,high) -> object_in_bounds o low high
    | Null -> false

  let longest_dimension low high = 
    let dx = ref neg_infinity and 
        dim = ref (-1) in
      for i = 0 to Array.length low - 1 do 
        let dxi = high.(i) -. low.(i) in 
          if dxi > !dx then begin
            dx := dxi;
            dim := i
          end
      done;
      !dim

  let all_equal = function 
    | [] -> true
    | [x] -> true
    | x :: xs -> 
        let xc = O.coord x in 
          List.for_all (fun y -> Pervasives.compare xc (O.coord y) = 0) xs

  let rec tree_of_objects objs low high = 
    match objs with 
      | [] -> Empty(low, high)
      | objs when all_equal objs -> 
          Cell(objs, low, high, Null, Null)
      | objs -> 
          let dim = longest_dimension low high in 
          let (new_low, new_high) = sub_bounds dim low high in 
          let (left, right) = 
            List.partition (fun obj -> object_in_bounds obj low new_high) objs in 
            Cell(objs, low, high,
                 tree_of_objects left low new_high,
                 tree_of_objects right new_low high)

  let vol low high = 
    let vol = ref 1.0 in 
      for i = 0 to Array.length low - 1 do 
        let dx = high.(i) -. low.(i) in 
          vol := !vol *. dx
      done;
      !vol +. 0.0

  let volume = function 
    | Empty(low, high) -> vol low high
    | Null -> 0.0
    | Cell(_, low, high, _, _) -> vol low high      

  let rec length_between a b = function 
    | [] -> 
        a <= 0 && 0 <= b 
    | x :: xs -> length_between (a-1) (b-1) xs

  let rec integrate f = function 
    | Null -> 0.0
    | Empty(_,_) -> 0.0
    | Cell([obj], low, high, _, _) -> 
        (vol low high)*.(f obj)
    | Cell(_,_,_,left,right) -> 
        (integrate f left) +. (integrate f right)
end
