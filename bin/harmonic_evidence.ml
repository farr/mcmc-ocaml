open Printf

module Ev = Evidence.Make(struct
  type params = float array
  let to_coords x = x
end)

let nbstrap = ref 10000
let ifile = ref ""

let options = 
  [("-nbstrap", Arg.Set_int nbstrap,
    sprintf "number of bootstrap samples to use to estimate error (default %d)" !nbstrap);
   ("-i", Arg.Set_string ifile,
    sprintf "input filename");
   ("-seed", Arg.Int (fun s -> Random.init s),
    "seed the RNG used for bootstrap")]

let _ = 
  Random.self_init ();
  Arg.parse options (fun _ -> ()) "harmonic_evidence.{byte,native} OPTIONS ...";
  let inp = open_in !ifile in 
  let samples = Read_write.read (fun x -> x) inp in 
    close_in inp;
    let bsamples = Array.copy samples in 
    let ev = Ev.evidence_harmonic_mean samples in 
    let evs = Array.init !nbstrap
      (fun _ -> 
        let n = Array.length bsamples in 
          for i = 0 to n - 1 do 
            bsamples.(i) <- samples.(Random.int n)
          done;
          Ev.evidence_harmonic_mean bsamples) in
      Array.fast_sort (fun (x : float) y -> Pervasives.compare x y) evs;
      let ilow = (!nbstrap) / 20 and 
          ihigh = (!nbstrap * 19) / 20 in 
        printf "    Best        10%%         90%%    \n%10g %10g %10g\n" ev evs.(ilow) evs.(ihigh)
        
