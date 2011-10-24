(*  read_write.ml: Reading and writing MCMC samples to/from files.
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

let write_sample to_coords chan 
    {value = v;
     like_prior = {log_likelihood = ll;
                   log_prior = lp}} = 
  Array.iter (fun x -> Printf.fprintf chan "%g " x) (to_coords v);
  Printf.fprintf chan "%g %g\n" ll lp 

let write to_coords chan mcmcs = 
  Array.iter 
    (fun samp -> write_sample to_coords chan samp)
    mcmcs

let read_one_line from_coords string = 
  let ea = Earray.make 0 0.0 in 
  let buf = Scanf.Scanning.from_string string in
    (try 
       while true do 
         Scanf.bscanf buf " %g " (fun x -> Earray.append ea x)
       done
     with 
       | End_of_file -> ());
    let pts = Earray.to_array ea in 
    let n = Array.length pts in 
      {value = from_coords (Array.sub pts 0 (n-2));
       like_prior = {log_likelihood = pts.(n-2);
                     log_prior = pts.(n-1)}}

let read_sample from_coords chan = 
  let line = input_line chan in 
    read_one_line from_coords line
      
let read from_coords chan = 
  let ea = Earray.of_array [||] in 
    try 
      while true do 
        Earray.append ea (read_sample from_coords chan)
      done;
      Earray.to_array ea
    with 
      | End_of_file -> Earray.to_array ea

let write_nested to_coords chan (log_ev, log_dev, all_pts, wts) = 
  Printf.fprintf chan "%g %g\n" log_ev log_dev;
  Array.iteri (fun i sample -> 
    let wt = wts.(i) in 
      Array.iter (fun x -> Printf.fprintf chan "%g " x) (to_coords sample.value);
      Printf.fprintf chan "%g %g " sample.like_prior.log_likelihood sample.like_prior.log_prior;
      Printf.fprintf chan "%g\n" wt)
    all_pts

let read_nested_log_ev_dev string = 
  let buf = Scanf.Scanning.from_string string in 
    Scanf.bscanf buf " %g %g " (fun x y -> (x,y))

let read_nested_sample from_coord string = 
  let buf = Scanf.Scanning.from_string string in 
  let ea = Earray.of_array [||] in 
    (try 
       while true do 
         Earray.append ea (Scanf.bscanf buf " %g " (fun x -> x))
       done;
       ()
     with 
       | End_of_file -> ());
    let pt = Earray.to_array ea in 
    let n = Array.length pt in 
      ({value = from_coord (Array.sub pt 0 (n-3));
        like_prior = {log_likelihood = pt.(n-3);
                      log_prior = pt.(n-2)}},
       pt.(n-1))

let read_nested from_coord chan = 
  let (log_ev, log_dev) = read_nested_log_ev_dev (input_line chan) in 
  let ea = Earray.of_array [||] in 
    (try 
       while true do
         Earray.append ea (read_nested_sample from_coord (input_line chan))
       done;
       ()
     with 
       | End_of_file -> ());
    let pts_and_wts = Earray.to_array ea in 
      (log_ev, log_dev, Array.map fst pts_and_wts, Array.map snd pts_and_wts)
