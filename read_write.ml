let write_sample to_coords chan 
    {Mcmc.value = v;
     like_prior = {Mcmc.log_likelihood = ll;
                   Mcmc.log_prior = lp}} = 
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
      {Mcmc.value = from_coords (Array.sub pts 0 (n-2));
       like_prior = {Mcmc.log_likelihood = pts.(n-2);
                     Mcmc.log_prior = pts.(n-1)}}

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
