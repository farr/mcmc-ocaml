open OUnit

let random_array nmin nmax = 
  let n = nmin + (Random.int (nmax-nmin)) in 
  let a = Array.make n 0.0 in 
    for i = 0 to n - 1 do 
      a.(i) <- Random.float 1.0
    done;
    a

let test_to_of_array () = 
  let a = random_array 10 50 in 
  let ea = Earray.of_array a in 
    assert_equal (Earray.to_array ea) a

let test_append () = 
  let a = random_array 10 50 in 
  let b = random_array 10 50 in 
  let ea = Earray.of_array a in 
  let ab = Array.append a b in 
    for i = 0 to Array.length b - 1 do 
      Earray.append ea b.(i)
    done;
    assert_equal (Earray.to_array ea) ab

let test_get_set () = 
  let a = random_array 10 50 in 
  let ea = Earray.of_array a in 
  let i = Random.int (Array.length a) in 
    assert_equal (Earray.get ea i) (Array.get a i);
    let x = Random.float 1.0 in 
      Array.set a i x;
      Earray.set ea i x;
      assert_equal (Earray.to_array ea) a

let test_chop_length () = 
  let a = random_array 10 50 in 
  let nsmall = (Array.length a)/2 in 
  let asmall = Array.sub a 0 nsmall in 
  let ea = Earray.of_array a in 
    Earray.chop ea nsmall;
    assert_equal (Earray.length ea) nsmall;
    assert_equal (Earray.to_array ea) asmall

let test_make () = 
  assert_equal (Earray.to_array (Earray.make 25 1.3)) (Array.make 25 1.3)

let tests = "earray.ml tests" >:::
  ["to_array, of_array" >:: test_to_of_array;
   "test_append" >:: test_append;
   "get set tests" >:: test_get_set;
   "chop and length" >:: test_chop_length;
   "make" >:: test_make]
