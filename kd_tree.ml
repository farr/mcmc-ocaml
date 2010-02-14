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
    let n1 = Array.length c1 and 
        n2 = Array.length c2 in 
    let c = Pervasives.compare n1 n2 in 
      if c <> 0 then c else 
        let rec loop i = 
          if i = n1 then 0 else 
            let x = c1.(i) and 
                y = c2.(i) in 
            let c = Pervasives.compare x y in 
              if c <> 0 then c else
                loop (i+1) in 
          loop n1

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
    | o :: _ as objs when all_equal compare_all_coords objs -> 
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

