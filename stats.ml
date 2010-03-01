let mean xs = 
  let sum = ref 0.0 and
      n = Array.length xs in 
    for i = 0 to n - 1 do 
      sum := !sum +. xs.(i)
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

let pi = 4.0 *. (atan 1.0)

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
