let mean xs = 
  let sum = ref 0.0 and
      n = Array.length xs in 
    for i = 0 to n - 1 do 
      sum := !sum +. xs.(i)
    done;
    !sum /. (float_of_int n)

let meanf f xs = 
  let sum = ref 0.0 and 
      n = Array.length xs in 
    for i = 0 to n - 1 do 
      sum := !sum +. (f xs.(i))
    done;
    !sum /. (float_of_int n)

let std_mean = mean

let std ?mean xs = 
  let mean = match mean with | None -> std_mean xs | Some(mu) -> mu in 
  let var = ref 0.0 and
      n = Array.length xs in 
    for i = 0 to n - 1 do 
      let x = xs.(i) -. mean in 
        var := !var +. x*.x
    done;
    sqrt (!var /. (float_of_int (n-1)))

let stdf ?mean f xs = 
  let mean = match mean with | None -> meanf f xs | Some(mu) -> mu in 
  let n = Array.length xs and 
      sum = ref 0.0 in 
    for i = 0 to Array.length xs - 1 do 
      let x = f xs.(i) in 
      let dx = x -. mean in 
        sum := !sum +. dx*.dx
    done;
    sqrt (!sum /. (float_of_int (n-1)))

let pi = 4.0 *. (atan 1.0)

let multi_mean xs = 
  let mu = Array.make (Array.length xs.(0)) 0.0 and 
      n = float_of_int (Array.length xs) in 
    for i = 0 to Array.length xs - 1 do 
      let x = xs.(i) in 
        for j = 0 to Array.length x - 1 do 
          mu.(j) <- mu.(j) +. x.(j)
        done
    done;
    for i = 0 to Array.length mu - 1 do 
      mu.(i) <- mu.(i) /. n
    done;
    mu

let multi_std ?mean xs = 
  let mean = match mean with | None -> multi_mean xs | Some(mu) -> mu in 
  let n = Array.length xs and 
      dim = Array.length xs.(0) in
  let std = Array.make dim 0.0 in 
    for i = 0 to n - 1 do 
      let x = xs.(i) in 
        for j = 0 to dim - 1 do 
          let dx = x.(j) -. mean.(j) in 
            std.(j) <- std.(j) +. dx*.dx
        done
    done;
    for i = 0 to dim - 1 do 
      std.(i) <- sqrt (std.(i) /. (float_of_int (n-1)))
    done;
    std

let draw_cauchy x0 gamma = 
  let p = Random.float 1.0 in 
    x0 +. gamma*.(tan (pi *. (p -. 0.5)))
    
let log_cauchy x0 gamma x = 
  let dx = (x -. x0) /. gamma in 
  0.0 -. (log (pi *. gamma)) -. 
    (log (1.0 +. dx*.dx))

let log_gaussian mu sigma x = 
  let dx = (x -. mu)/.sigma in 
    (-0.5)*.(log (2.0 *. pi)) -. (log sigma) 
      -. 0.5*.dx*.dx

(* Method from Press, Teukolsky, Vetterling and Flannery.  Numerical
   Recipes.  Cambridge University Press, 2007 (third edition).
   Section 7.3.9; original algorithm by Leva.*)
let draw_gaussian mu sigma = 
  let rec loop () = 
    let u = Random.float 1.0 and 
        v = 1.7156 *. (Random.float 1.0 -. 0.5) in 
    let x = u -. 0.449871 and 
        y = abs_float v +. 0.386595 in 
    let q = x*.x +. y*.(0.19600*.y -. 0.25472*.x) in 
      if q > 0.27597 && (q > 0.27846 || v*.v > (-4.0)*.(log u)*.u*.u) then 
        loop ()
      else
        mu +. sigma*.v/.u in 
    loop ()

let all_equal_float (xs : float array) start endd = 
  let x0 = xs.(start) in 
  let rec ae_loop i = 
    if i >= endd then 
      true
    else if xs.(i) = x0 then 
      ae_loop (i+1)
    else
      false in 
    ae_loop (start+1)
      

let find_nth ?(copy = true) nth (xs : float array) = 
  if nth < 0 || nth >= (Array.length xs) then raise (Invalid_argument "find_nth: nth outside array bounds");
  let xs = if copy then Array.copy xs else xs in 
  let rec find_nth_loop start nth endd = 
    if endd - start = 1 then 
      if nth <> 0 then 
        raise (Failure "find_nth: nth index not in array bounds")
      else
        xs.(start)
    else if all_equal_float xs start endd then 
      xs.(start)
    else begin
      let part = xs.(start + (Random.int (endd-start))) in 
      let rec swap_loop low high = 
        let new_low = let rec new_low_loop l = 
                        if l >= endd then l else if xs.(l) <= part then new_low_loop (l+1) else l in new_low_loop low and 
            new_high = let rec new_high_loop h = 
                         if h < start then h else if xs.(h) > part then new_high_loop (h-1) else h in new_high_loop high in
          if new_low > new_high then new_low else if new_low >= endd then endd else if new_high < start then new_low else begin
            let tmp = xs.(new_low) in 
              xs.(new_low) <- xs.(new_high);
              xs.(new_high) <- tmp;
              swap_loop new_low new_high
          end in
      let ilow = swap_loop start (endd-1) in 
        if nth < (ilow - start) then 
          find_nth_loop start nth ilow
        else
          find_nth_loop ilow (nth - (ilow - start)) endd
    end in 
    find_nth_loop 0 nth (Array.length xs)

let all_equal compare xs start endd = 
  let x0 = xs.(start) in 
  let rec ae_loop i = 
    if i >= endd then 
      true
    else if compare x0 xs.(i) = 0 then 
      ae_loop (i+1)
    else
      false in 
    ae_loop (start+1)
      
let find_nthf ?(copy = true) compare nth xs = 
  if nth < 0 || nth >= Array.length xs then raise (Invalid_argument "find_nthf: nth outside array bounds");
  let xs = if copy then Array.copy xs else xs in 
  let rec find_nth_loop start nth endd = 
    if endd - start = 1 then 
      if nth <> 0 then 
        raise (Failure "find_nthf: nth index not in array bounds")
      else
        xs.(start)
    else if all_equal compare xs start endd then 
      xs.(start)
    else begin
      let part = xs.(start + (Random.int (endd - start))) in 
      let rec swap_loop low high = 
        let new_low = let rec new_low_loop l = 
                        if l >= endd then l else if compare xs.(l) part <= 0 then new_low_loop (l+1) else l in new_low_loop low and
            new_high = let rec new_high_loop h = 
                         if h < start then h else if compare xs.(h) part > 0 then new_high_loop (h-1) else h in new_high_loop high in
          if new_low > new_high then new_low else if new_low >= endd then endd else if new_high < start then new_low else begin
            let tmp = xs.(new_low) in 
              xs.(new_low) <- xs.(new_high);
              xs.(new_high) <- tmp;
              swap_loop new_low new_high
          end in
      let ilow = swap_loop start (endd - 1) in
        if nth < (ilow - start) then 
          find_nth_loop start nth ilow
        else
          find_nth_loop ilow (nth - (ilow - start)) endd
    end in 
    find_nth_loop 0 nth (Array.length xs)
