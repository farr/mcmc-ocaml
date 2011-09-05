open OUnit
open Ellipse
open Lacaml.Impl.D
open Asserts
open Stats
open Mcmc

let mat_mul a b = 
  let n = Array.length a and 
      m = Array.length b.(0) and 
      l = Array.length a.(0) in 
  let res = Array.make_matrix n m 0.0 in 
    for i = 0 to n - 1 do 
      for j = 0 to m - 1 do 
        for k = 0 to l - 1 do 
          res.(i).(j) <- res.(i).(j) +. a.(i).(k)*.b.(k).(j)
        done
      done
    done;
    res

let random_symmetric_matrix n = 
  let mat = Array.make_matrix n n 0.0 in 
    for i = 0 to n - 1 do 
      mat.(i).(i) <- draw_uniform (-1.0) 1.0;
      for j = i+1 to n-1 do 
        let x = draw_uniform (-1.0) 1.0 in 
          mat.(i).(j) <- x;
          mat.(j).(i) <- x
      done
    done;
    mat

let random_point n = 
  let arr = Array.make n 0.0 in 
    for i = 0 to n - 1 do 
      arr.(i) <- draw_uniform (-1.0) 1.0
    done;
    {value = arr;
     like_prior = {log_likelihood = neg_infinity; log_prior = neg_infinity}}

let random_points ndim npts = 
  Array.init npts (fun _ -> random_point ndim)

let transpose mat = 
  let n = Array.length mat and 
      m = Array.length mat.(0) in 
  let mt = Array.make_matrix m n 0.0 in 
    for i = 0 to n - 1 do 
      for j = 0 to m - 1 do 
        mt.(j).(i) <- mat.(i).(j)
      done
    done;
    mt

let test_eigensystem () = 
  for i = 0 to 100 do 
    let m = random_symmetric_matrix 5 in 
    let (evals, evecs) = eigensystem m in 
    let evmat = Array.make_matrix 5 5 0.0 in 
      for i = 0 to 4 do 
        evmat.(i).(i) <- evals.(i)
      done;
      let m_reconstr = mat_mul evecs (mat_mul evmat (transpose evecs)) in 
        assert_equal_float_matrix ~msg:"reconstruction differs" m m_reconstr
  done

let test_enclose () = 
  for i = 0 to 100 do 
    let pts = random_points 5 1000 in 
    let ell = enclosing_ellipse 2.0 (fun x -> x) pts in 
      for i = 0 to 999 do 
        let pt = pts.(i).value in 
        let r = elliptical_range ell pt in 
          assert_bool "point outside ellipse" (r < 1.0)
      done
  done

let tests = "ellipse.ml tests" >:::
  ["reconstructed eigensystem" >:: test_eigensystem;
   "enclose test" >:: test_enclose]
