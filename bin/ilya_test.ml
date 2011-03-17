(* ilya_test.ml: Test program for admixture MCMC. *)
(* Copyright (C) 2011 Will M. Farr <w-farr@northwestern.edu> *)

(* This program is free software: you can redistribute it and/or modify *)
(* it under the terms of the GNU General Public License as published by *)
(* the Free Software Foundation, either version 3 of the License, or *)
(* (at your option) any later version. *)

(* This program is distributed in the hope that it will be useful, *)
(* but WITHOUT ANY WARRANTY; without even the implied warranty of *)
(* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the *)
(* GNU General Public License for more details. *)

(* You should have received a copy of the GNU General Public License *)
(* along with this program.  If not, see <http://www.gnu.org/licenses/>. *)

(* A is uniform on 0 to 1, B is sharply peaked. *)
type state =
  | A of float * float
  | B of float * float

let _ = Random.self_init ()

let bmin = 0.4
let bmax = 0.6

let log_likelihood = function 
  | A(x,_) -> if 0.0 <= x && x <= 1.0 then 0.0 else neg_infinity
  | B(_,x) -> if bmin <= x && x <= bmax then ~-.(log (bmax -. bmin)) else neg_infinity

let log_prior = function
  | A(x,y) | B(x,y) -> if 0.0 <= x && x <= 1.0 && 0.0 <= y && y <= 1.0 then 0.0 else neg_infinity

let p_switch = 0.1

let jump_proposal = function 
  | A(x,y) -> 
      if Random.float 1.0 < p_switch then 
        B(Random.float 1.0, Random.float 1.0)
      else
        A(Random.float 1.0, y)
  | B(x,y) -> 
      if Random.float 1.0 < p_switch then 
        A(Random.float 1.0, Random.float 1.0)
      else
        B(x,Random.float 1.0)

let log_jp _ _ = 0.0 (* Symmetric? *)

let samples = 
  Mcmc.mcmc_array
    1000000
    log_likelihood
    log_prior
    jump_proposal
    log_jp
    (B(0.5, 0.5))

let count_ab samples = 
  Array.fold_left
    (fun (na,nb) state -> 
       let {Mcmc.value = v} = state in 
       match v with 
         | A(_,_) -> (na+1, nb)
         | B(_,_) -> (na, nb+1))
    (0,0)
    samples

let _ = 
  let (na,nb) = count_ab samples in 
    Printf.printf "Number in a = %d, number in b = %d (evidence ratio = %g)\n"
      na
      nb
      ((float_of_int na)/.(float_of_int nb))

let alog_like x = if 0.0 <= x && x <= 1.0 then 0.0 else neg_infinity
let blog_like x = if bmin <= x && x <= bmax then ~-.(log (bmax-.bmin)) else neg_infinity
let log_prior _ = 0.0
let propose _ = Random.float 1.0
let log_jp _ _ = 0.0

let ad_samples = 
  Mcmc.admixture_mcmc_array
    1000000
    (alog_like, blog_like)
    (log_prior, log_prior)
    (propose, propose)
    (log_jp, log_jp)
    (0.5, 0.5)
    (1.0, 1.0)
    (0.5, 0.1)

let _ = 
  let mean_lam = Stats.meanf (fun samp -> let {Mcmc.value = (lam,_,_)} = samp in lam) ad_samples in 
    Printf.printf "\n\nMean lambda = %g\n" mean_lam;
    Printf.printf "Mean lambda ratio = %g\n" (1.0 /. (2.0 -. 3.0*.mean_lam) -. 1.0);
    Printf.printf "Max-likelihood ratio estimator = %g\n" (Mcmc.max_like_admixture_ratio ad_samples)
