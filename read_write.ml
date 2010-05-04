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

let read from_coords chan = 
  let coords_ll_lp = Parse_floats.parse_floats Lex_floats.lex_floats (Lexing.from_channel chan) in 
    Array.of_list
      (List.rev 
         (List.rev_map 
            (fun v -> 
               let n = Array.length v in 
                 {Mcmc.value = from_coords (Array.sub v 0 (n - 2));
                  like_prior = {Mcmc.log_likelihood = v.(n-2);
                                log_prior = v.(n-1)}})
            coords_ll_lp))
    
