open Asserts
open OUnit

module Interp = Interpolate_pdf.Make(
  struct
    type point = float array

    let coord (p : point) = p
    let point (c : float array) = c
  end)

open Interp

let random_point () = 
  Array.init 2 (fun _ -> Random.float 1.0)

let pi = 4.0 *. (atan 1.0)

let linear_pdf a b x = 
  2.0*.x/.(b*.b -. a*.a)

let draw_linear a b = 
  assert(a >= 0.0);
  assert(b > a);
  let y = Random.float 1.0 in 
    sqrt (b*.b*.y +. a*.a*.(1.0 -. y))

let linear_mean a b = 
  2.0/.3.0*.(b +. a*.a/.(a +. b))

let mean pts = 
  let n = Array.length pts in 
  let nf = float_of_int n in 
  let mu = Array.make 2 0.0 in 
    for i = 0 to n - 1 do 
      mu.(0) <- mu.(0) +. pts.(i).(0);
      mu.(1) <- mu.(1) +. pts.(i).(1)
    done;
    mu.(0) <- mu.(0) /. nf;
    mu.(1) <- mu.(1) /. nf;
    mu

let test_mean_distribution () = 
  let b0 = Random.float 2.0 and 
      b1 = Random.float 4.0 in 
  let pts = 
    Array.init 10000 
      (fun _ -> [|draw_linear 0.0 b0; draw_linear 0.0 b1|]) in 
  let pdf = make pts [|0.0; 0.0|] [|b0; b1|] in 
  let pts = Array.init 10000 (fun _ -> draw pdf) in 
  let mu = mean pts in 
    assert_equal_float ~epsrel:0.05 ~msg:"first component" (linear_mean 0.0 b0) mu.(0);
    assert_equal_float ~epsrel:0.05 ~msg:"second component" (linear_mean 0.0 b1) mu.(1)

let tests = "interpolate_pdf.ml tests" >:::
  ["mean of linear PDF in 2D" >:: test_mean_distribution]
