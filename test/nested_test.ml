(*  nested_test.ml: Test suite for the nested module.
    Copyright (C) 2011 Will M. Farr <w-farr@northwestern.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. *)

open Asserts
open Mcmc
open Nested
open OUnit
open Stats

let test_single_gaussian () = 
  let log_prior x = 
    if x.(0) < 1.0 && x.(0) > 0.0 && x.(1) < 1.0 && x.(1) > 0.0 then 
      0.0
    else
      neg_infinity in 
  let log_likelihood x = 
    let g1 = log_gaussian 0.5 0.1 x.(0) and 
        g2 = log_gaussian 0.5 0.1 x.(1) in 
      g1 +. g2 in 
  let draw_prior () = 
    [|draw_uniform 0.0 1.0; draw_uniform 0.0 1.0|] in 
  let (ev, dev, _, _) = nested_evidence draw_prior log_likelihood log_prior in 
    assert_equal_float ~msg:"bad Gaussian evidence" ~epsabs:0.1 1.0 ev;
    assert_bool "error estimate too large" ((total_error_estimate ev dev 1000) < 0.1)

let test_four_gaussians () = 
  let mu1 = [|0.25; 0.25|] and 
      mu2 = [|0.25; 0.75|] and 
      mu3 = [|0.75; 0.25|] and 
      mu4 = [|0.75; 0.75|] in 
  let sigma = [|0.05; 0.05|] in 
  let log_prior x = 
    if x.(0) < 1.0 && x.(0) > 0.0 && x.(1) < 1.0 && x.(1) > 0.0 then 
      0.0
    else
      neg_infinity in 
  let log_likelihood x = 
    let g1 = log_multi_gaussian mu1 sigma x and 
        g2 = log_multi_gaussian mu2 sigma x and 
        g3 = log_multi_gaussian mu3 sigma x and 
        g4 = log_multi_gaussian mu4 sigma x in 
      log ((exp g1) +. (exp g2) +. (exp g3) +. (exp g4)) in 
  let draw_prior () = 
    [|draw_uniform 0.0 1.0; draw_uniform 0.0 1.0|] in 
  let (ev, dev, _, _) = nested_evidence draw_prior log_likelihood log_prior in 
  let err = total_error_estimate ev dev 1000 in 
    assert_equal_float ~msg:"bad Gaussian evidence" ~epsabs:0.5 4.0 ev;
    assert_bool (Printf.sprintf "error estimate too large: %g" err)  (err < 0.5)

let test_single_gaussian_weights () = 
  let log_prior x = 
    if x.(0) < 1.0 && x.(0) > 0.0 && x.(1) < 1.0 && x.(1) > 0.0 then 
      0.0
    else
      neg_infinity in 
  let log_likelihood x = 
    let g1 = log_gaussian 0.5 0.1 x.(0) and 
        g2 = log_gaussian 0.5 0.1 x.(1) in 
      g1 +. g2 in 
  let draw_prior () = 
    [|draw_uniform 0.0 1.0; draw_uniform 0.0 1.0|] in 
  let nlive = 1000 in 
  let (ev, dev, pts, wts) = nested_evidence ~nlive:nlive draw_prior log_likelihood log_prior in 
  let mean = Util.UArray.fold_left2 (fun sum pt wt -> sum +. wt*.pt.value.(0)) 0.0 pts wts in 
  let wt_sum = Array.fold_left (+.) 0.0 wts in 
    assert_equal_float ~msg:"Bad weight sum" 1.0 wt_sum;
    assert_equal_float ~msg:"Bad mean" ~epsabs:0.1 0.5 mean

let tests = "nested.ml tests" >:::
  ["single Gaussian test" >:: test_single_gaussian;
   "four Gaussian test" >:: test_four_gaussians;
   "single Gaussian weights" >:: test_single_gaussian_weights]
