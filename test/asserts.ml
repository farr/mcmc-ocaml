open OUnit

let assert_equal_float ?(epsabs = 1e-8) ?(epsrel = 1e-8) ?(msg = "") = 
  assert_equal ~msg:msg ~cmp:(cmp_float ~epsabs:epsabs ~epsrel:epsrel) ~printer:string_of_float

let assert_equal_float_array ?(epsabs = 1e-8) ?(epsrel = 1e-8) ?(msg = "assert_equal_float_array") x y =
  let n = Array.length x in 
    assert_bool (Printf.sprintf "%s: array lengths differ" msg) (Array.length y = n);
    for i = 0 to n - 1 do 
      assert_equal_float ~epsabs:epsabs ~epsrel:epsrel ~msg:msg x.(i) y.(i)
    done
