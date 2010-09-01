type 'a earray = 
    {mutable elts : 'a array;
     mutable size : int}

let get ev i = 
  if i < ev.size then 
    ev.elts.(i)
  else
    raise (Invalid_argument "get: requested element out of range")

let append ev x = 
  let nelts = Array.length ev.elts in 
    if ev.size = nelts then begin
      ev.elts <- Array.append ev.elts (Array.make (ev.size+1) x); (* +1 ensures that we always have at least one spot. *)
    end else begin
      ev.elts.(ev.size) <- x;
    end;
    ev.size <- ev.size + 1

let set ev i x = 
  if i < ev.size then 
    ev.elts.(i) <- x
  else if i = ev.size then 
    append ev x
  else
    raise (Invalid_argument "set: index out of range")

let to_array ev =
  Array.sub ev.elts 0 ev.size

let chop ev size = 
  if size < 0 then 
    raise (Invalid_argument "chop: must have nonnegative size")
  else if size > ev.size then 
    raise (Invalid_argument "chop: cannot extend array beyond its current size")
  else begin
    ev.elts <- Array.sub ev.elts 0 size;
    ev.size <- size
  end

let make size init = 
  if size >= 0 then 
    {elts = Array.make size init;
     size = size}
  else
    raise (Invalid_argument "make: must have nonnegative size")

let length ev = ev.size

let of_array arr = 
  {elts = Array.copy arr;
   size = Array.length arr}

let copy ea = 
  {elts = Array.copy ea.elts;
   size = ea.size}
