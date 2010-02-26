open OUnit

let assert_equal_float ?(epsabs = 1e-8) ?(epsrel = 1e-8) = 
  assert_equal ~cmp:(cmp_float ~epsabs:epsabs ~epsrel:epsrel) ~printer:string_of_float
