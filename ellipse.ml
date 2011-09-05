(* ellipse.ml: Recursive ellipse decomposition of a domain of points.
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

open Mcmc
open Stats
open Lacaml.Impl.D

type ellipse = 
    {center : float array;
     axes : float array;
     orientation : float array array}

type 'a ellipse_tree_data = 
    {pts : 'a mcmc_sample array;
     left : 'a ellipse_tree;
     right : 'a ellipse_tree;
     ellipse : ellipse;
     circumcircle : (float array) * float}
and 'a ellipse_tree = 
  | Empty 
  | Tree of 'a ellipse_tree_data
  
let center to_coord pts = 
  let ndim = Array.length (to_coord pts.(0).value) and 
      n = Array.length pts in 
  let nf = float_of_int n in 
  let cen = Array.make ndim 0.0 in 
    for i = 0 to n - 1 do 
      let pt = to_coord pts.(i).value in 
      for j = 0 to ndim - 1 do 
        cen.(j) <- cen.(j) +. pt.(j)/.nf
      done
    done;
    cen

let sigma2 to_coord mu pts = 
  let ndim = Array.length (to_coord pts.(0).value) and 
      n = Array.length pts in 
  let nf = float_of_int n in 
  let sigma = Array.make_matrix ndim ndim 0.0 in 
    for i = 0 to n - 1 do 
      let pt = to_coord pts.(i).value in 
        for j = 0 to ndim - 1 do 
          let dxj = pt.(j) -. mu.(j) in
            sigma.(j).(j) <- sigma.(j).(j) +. dxj*.dxj/.nf;
            for k = j + 1 to ndim - 1 do 
              let dxk = pt.(k) -. mu.(k) in 
              let delta = dxj*.dxk/.nf in
                sigma.(j).(k) <- sigma.(j).(k) +. delta;
                sigma.(k).(j) <- sigma.(k).(j) +. delta
            done
        done
    done;
    sigma

let eigensystem sigma = 
  let mat = Mat.of_array sigma in 
  let (_, w, z, _) = syevr ~vectors:true mat in 
    (Vec.to_array w,  (Mat.to_array z))

let elliptical_range {center = c; axes = a; orientation = ori} pt = 
  let r = ref 0.0 in 
    for j = 0 to Array.length ori.(0) - 1 do 
      let d = ref 0.0 in 
        for i = 0 to Array.length ori - 1 do 
          d := !d +. (pt.(i) -. c.(i))*.ori.(i).(j)
        done;
        r := !r +. !d *. !d /. a.(j)
    done;
    !r +. 0.0

let max_elliptical_range to_coord ell pts = 
  Array.fold_left 
    (fun max_r pt ->
      let pt = to_coord pt.value in 
        max max_r (elliptical_range ell pt))
    neg_infinity
    pts

let rescale_ellipse sf to_coord ell pts = 
  let r_max = max_elliptical_range to_coord ell pts in 
  let dim_sf = sf ** (1.0 /. (float_of_int (Array.length ell.axes))) in 
    {ell with axes = Array.map (fun x -> x*.dim_sf*.r_max) ell.axes}

let enclosing_ellipse sf to_coord pts = 
  let cen = center to_coord pts in 
  let sigma2 = sigma2 to_coord cen pts in 
  let (evals, evecs) = eigensystem sigma2 in 
  let ell0 = {center = cen; axes = evals; orientation = evecs} in 
    rescale_ellipse sf to_coord ell0 pts

let distance p1 p2 = 
  let r = ref 0.0 in 
    for i = 0 to Array.length p1 - 1 do 
      let dx = p1.(i) -. p2.(i) in 
        r := !r +. dx*.dx
    done;
    sqrt !r

let union_circumcircles ((c1, r1) as cc1)  ((c2, r2) as cc2) =
  let r12 = distance c1 c2 in 
    if r12 +. r2 < r1 then 
      cc1 
    else if r12 +. r1 < r2 then 
      cc2 
    else begin 
      let ndim = Array.length c1 in 
      let rnew = r1 +. r12 +. r2 in 
      let cnew = Array.make ndim 0.0 in 
      let mag = 0.5*.(r2 +. r12 -. r1)/.r12 in 
        for i = 0 to ndim - 1 do 
          cnew.(i) <- c1.(i) +. mag*.(c2.(i) -. c1.(i))
        done;
        (cnew,rnew)
    end

let in_circumcircle pt (c1,r1) = 
  let r = distance pt c1 in 
    r < r1

let ellipse_circumcircle {center = c; axes = a} = 
  (c, Array.fold_left max neg_infinity a)

let widest_dimension {axes = a} = 
  let imax = ref (-1) and 
      rmax = ref neg_infinity in
    for i = 0 to Array.length a - 1 do 
      if a.(i) > !rmax then begin
        rmax := a.(i);
        imax := i
      end
    done;
    !imax

let array_partition pred arr = 
  let (l1, l2) = List.partition pred (Array.to_list arr) in 
    (Array.of_list l1, Array.of_list l2)

let rec ellipse_tree sf to_coord pts = 
  let ell = enclosing_ellipse sf to_coord pts in 
  let ndim = Array.length ell.axes in 
    assert(Array.length pts >= ndim + 1);
    let split = widest_dimension ell in 
    let (left_pts, right_pts) = 
      array_partition (fun pt -> (to_coord pt.value).(split) < ell.center.(split)) pts in 
    let left = if Array.length left_pts < ndim + 1 then Empty else ellipse_tree sf to_coord left_pts and 
        right = if Array.length right_pts < ndim + 1 then Empty else ellipse_tree sf to_coord right_pts in 
    let enclosing_circumcircle = 
      match left, right with 
        | Empty, Empty -> 
          ellipse_circumcircle ell
        | Tree({circumcircle = c}), Empty | Empty, Tree({circumcircle = c}) -> 
          union_circumcircles c (ellipse_circumcircle ell)
        | Tree({circumcircle = cl}), Tree({circumcircle = cr}) -> 
          union_circumcircles (union_circumcircles cl cr) (ellipse_circumcircle ell) in 
      Tree({pts = pts;
            left = left;
            right = right;
            ellipse = ell;
            circumcircle = enclosing_circumcircle})


