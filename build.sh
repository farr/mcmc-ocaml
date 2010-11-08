#!/bin/bash

set -e

ocb () {
    ocamlbuild $*
}

rule () {
    case $1 in
        all*) ocb mcmc.cma mcmc.cmxa mcmc.docdir/index.html;;
        doc*) ocb mcmc.docdir/index.html;;
        lib*) ocb mcmc.cma mcmc.cmxa;;
        tes*) ocb run_tests.native && ./run_tests.native;;
        clean*) ocb -clean;;
        *) echo "Unknown action";;
    esac
}

if [ $# -eq 0 ]; then
    rule all
else 
    while [ $# -gt 0 ]; do
        rule $1;
        shift;
    done
fi