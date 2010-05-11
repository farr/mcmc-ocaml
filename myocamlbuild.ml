open Ocamlbuild_plugin

let oUnit_dir = "/Users/farr/Documents/code/oUnit"
let ocamlmpi_dir = "/Users/farr/Documents/code/ocamlmpi"

let _ = dispatch begin function 
  | After_rules -> 
      ocaml_lib "mcmc";
      ocaml_lib ~extern:true ~dir:oUnit_dir "oUnit";
      ocaml_lib ~extern:true ~dir:ocamlmpi_dir "mpi"
  | _ -> ()
end
