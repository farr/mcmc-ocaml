(** This is a corrected test of Rodriguez's "extreme" example from
    arXiv:0709.1067.  We draw parameters (a,b) from N(0,100), and then
    conduct Nsamp repetitions of an experiment that consists in
    drawing (x,y) from

    x = exp(-50(a^2 + b^2)) + N(0,1)
    y = exp(-50((a-c)^2 + (b-c)^2)) + N(0,1)

    where c = 0.1 is assumed known.  We wish to compare the
    performance of two different priors used to recover the parameters
    (a,b).  The first is the "true" prior, independent N(0,100); the
    second is the "uninformative" prior, proportional to the sqrt of
    the Fisher Information matrix for the likelihood:

    pi(a,b) \propto |a-b| exp(-100(a^2+b^2) + 10(a+b))

    The program will run some number Ntrials of such recoveries,
    reporting the standard deviation on a and b after Nsamp draws of
    (x,y) pairs for each prior used in the analysis.
*)

open Printf
open Stats

let nsamp = ref 100
let ntrials = ref 100

let pwidth = ref 1.0 

let options = 
  [("-nsamp", Arg.Set_int nsamp, sprintf "N the number of (x,y) samples to draw for each experiment (default %d)" !nsamp);
   ("-ntrials", Arg.Set_int ntrials, sprintf "N the number of prior comparison trials to conduct (default %d)" !ntrials);
   ("-pwidth", Arg.Set_float pwidth, sprintf "sigma the width of the prior for (a,b) draw in each trial (default %g)" !pwidth)]

let mus ab = 
  let a = ab.(0) and 
      b = ab.(1) in 
  let amc = a -. 0.1 and
      bmc = b -. 0.1 in 
    [|exp(-50.0*.(a*.a +. b*.b)); 
      exp(-50.0*.(amc*.amc +. bmc*.bmc))|]

let drawxy ab = 
  match mus ab with 
    | [|mua; mub|] -> 
      let dx = draw_gaussian 0.0 1.0 and 
          dy = draw_gaussian 0.0 1.0 in 
        [|mua +. dx; mub +. dy|]
    | _ -> raise (Failure "drawxy: malformed ab argument")

let log_pi0 ab = 
  match ab with 
    | [|a; b|] -> (log_gaussian 0.0 !pwidth a) +. (log_gaussian 0.0 !pwidth b)
    | _ -> raise (Failure "log_pi0: malformed ab argument")

let log_pi1 ab = 
  match ab with 
    | [|a; b|] -> 
      let d = abs_float (a -. b) in 
        d *. (exp (-100.0*.(a*.a +. b*.b) +. 10.0*.(a+.b)))
    | _ -> raise (Failure "log_pi1: malformed ab argument")

let log_likelihood xys ab = 
  match mus ab with 
    | [|mua; mub|] -> 
      let logl = ref 0.0 in 
        for i = 0 to Array.length xys - 1 do 
          let xy = xys.(i) in 
          let dx = mua -. xy.(0) and 
              dy = mub -. xy.(1) in 
            logl := !logl +. (log_gaussian 0.0 1.0 dx) +. (log_gaussian 0.0 1.0 dy)
        done;
        !logl +. 0.0
    | _ -> raise (Failure "log_likelihood: malformed ab argument")

let sqr x = x *. x

let rec mus_to_as mus abigger = 
  let sign = if abigger then -1.0 else 1.0 in 
  let log_ma = log mus.(0) and 
      log_mb = log mus.(1) in 
  let disc = sign*.(sqrt (~-. (sqr log_mb) +. 2.0*.(log_ma -. 1.0)*.log_mb -. (sqr (log_ma +. 1.0)))) in
    [|0.05*.(1.0 -. log_ma +. log_mb -. disc);
      0.05*.(1.0 -. log_ma +. log_mb +. disc)|]

(* We compute the means from the given ab, and then jump with width ~
   1/sqrt(nsamp) in the means, and transform back to ab. *)
let jump_proposal ab = 
  match mus ab with 
    | [|mua; mub|] -> 
      let sigma = 1.0 /. (sqrt (float_of_int !nsamp)) in
      let dma = draw_uniform (~-.sigma) sigma and 
          dmb = draw_uniform (~-.sigma) sigma in 
      let new_mua = mua +. dma and 
          new_mub = mub +. dmb in 
        mus_to_as [|new_mua; new_mub|] (Random.float 1.0 > 0.5)
    | _ -> raise (Failure "jump_proposal: malformed ab argument")

(* This requires some explanation: since we jump uniformly in the
   means, we have 

   p(ab' | ab) da' db' = p(mu' | mu) dmu1' dmu2'

   Where the ab's and mu's satisfy the relations above.  But p(mu' |
   mu) is a constant (i.e. independent of mu' and mu), so

   p(ab' | ab) da' db' \propto dmu1' dmu2'

   or
   
   p(ab' | ab) \propto |d(mu1', mu2')/d(a',b')|

   where |d(x,y)/d(u,v)| stands for the absolute value of the
   determinant of the Jacobian of the map u,v |-> x(u,v), y(u,v).
   That is the quantity below:
*)
let log_jump_prob _ ab' = 
  match ab' with 
    | [|a'; b'|] ->
      let d = abs_float (a' -. b') in 
        (* Flat jumps in mu => p(a,b) da db = p(m1, m2) dm1 dm2 / det(d(m1,m2)/d(a,b)) *)
        (log d) +. ((-1.0) +. 10.0*.(1.0 -. 10.0*.a')*.a' +. 10.0*.(1.0 -. 10.0*.b')*.b')
    | _ -> raise (Failure "log_jump_prob: malformed ab argument")

let _ = 
  Arg.parse (Arg.align options) (fun _ -> ()) "correct_rodriguez_test.{byte,native} OPTIONS ...";
  Random.self_init ();
  let samples = 
    Mcmc.mcmc_array 
      ~nbin:1000
      ~nskip:1000
      10000
      (fun _ -> 0.0)
      log_pi0 
      jump_proposal
      log_jump_prob
      [|0.01; -0.01|] in 
  let out = open_out "prior.dat" in
    Read_write.write (fun x -> x) out samples;
    close_out out
    
