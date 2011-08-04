module UArray = struct
  let fold_left2 fn start a1 a2 = 
    let n = min (Array.length a1) (Array.length a2) in
    let accum = ref start in 
      for i = 0 to n - 1 do 
        accum := fn !accum a1.(i) a2.(i)
      done;
      !accum
end
