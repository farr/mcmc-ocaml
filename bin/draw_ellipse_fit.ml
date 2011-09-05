open Ellipse 
open Stats
open Graphics
open Mcmc

let pi = 3.1415926535897932385

let random_correlation_matrix () = 
  let a = draw_uniform 0.0 0.25 and 
      b = draw_uniform 0.0 0.25 and 
      theta = draw_uniform 0.0 (2.0*.pi) in 
  let c = cos theta and 
      s = sin theta in 
  let c2 = c*.c and 
      s2 = s*.s in 
    [|[|a*.c2 +. b*.s2; (a -. b)*.c*.s|];
      [|(a -. b)*.c*.s; b*.c2 +. a*.s2|]|]

let correlate_samples corr samps = 
  Array.map
    (fun samp -> 
      match corr,samp with 
        | [|[|c00; c01|];
            [|c10; c11|]|], [|x; y|] -> 
          [|x*.c00 +. y*.c01; x*.c01 +. y*.c11|]
        | _ -> raise (Failure "correlate_samples: Bad form"))
    samps

let offset_samples mu samps = 
  Array.map 
    (fun samp -> 
      match mu, samp with 
        | [|mux; muy|], [|x; y|] -> 
          [|x +. mux; y +. muy|]
        | _ -> raise (Failure "offset_samples: bad form"))
    samps

let draw_pt x y = 
  let nx = float_of_int (size_x ()) and 
      ny = float_of_int (size_y ()) in 
  let ix = int_of_float (x *. nx +. 0.5) and 
      iy = int_of_float (y *. ny +. 0.5) in 
    plot ix iy

let draw_ellipse_pt {center = cen; axes = a; orientation = o} theta = 
  let c = cos theta and 
      s = sin theta in 
  let x = cen.(0) +. c*.(sqrt a.(0))*.o.(0).(0) +. s*.(sqrt a.(1))*.o.(0).(1) and 
      y = cen.(1) +. c*.(sqrt a.(0))*.o.(1).(0) +. s*.(sqrt a.(1))*.o.(1).(1) in 
    draw_pt x y

let draw_ellipse ell npts = 
  for i = 0 to npts - 1 do 
    draw_ellipse_pt ell ((float_of_int i)*.2.0*.pi/.(float_of_int npts))
  done

let _ = 
  Random.self_init ();
  let s1 = random_correlation_matrix () and s2 = random_correlation_matrix () in 
  let mu1 = [|0.25; 0.25|] and 
      mu2 = [|0.75; 0.75|] in 
  let samps1 = 
    offset_samples mu1 (correlate_samples s1 (Array.init 250 (fun _ -> [|draw_gaussian 0.0 1.0; draw_gaussian 0.0 1.0|]))) and 
      samps2 = 
    offset_samples mu2 (correlate_samples s2 (Array.init 250 (fun _ -> [|draw_gaussian 0.0 1.0; draw_gaussian 0.0 1.0|]))) in 
  let samps = Array.append samps1 samps2 in 
  let samps = Array.map (fun pt -> {value = pt; like_prior = {log_likelihood = neg_infinity; log_prior = neg_infinity}}) samps in
  let tree = ellipse_tree 1.0 (fun x -> x) samps in
    Graphics.open_graph " 500x500";
    Array.iter (fun pt -> draw_pt pt.value.(0) pt.value.(1)) samps;
    let rec loop = function 
      | Empty -> ()
      | Tree({ellipse = ell; left = left; right = right}) -> 
        loop left;
        loop right;
        draw_ellipse ell 1000 in 
      loop tree;
      let _ = read_key () in 
        ()
