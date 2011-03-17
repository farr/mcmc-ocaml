(*  lex_floats.mll: Floating point lexer.
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

{ open Parse_floats }

let digit = ['0'-'9']
let uint = digit+
let pint = '+' uint
let nint = '-' uint
let int = uint | pint | nint

let exponent = 'e' | 'E'
let dot = '.' | ','

let eol = '\n'

let wspace = [' ' '\t' '\r']+

let float = int | (int dot int) | (int dot int exponent int) | (int exponent int)

let inf = "inf"

let neg_inf = "-inf"

let nan = "nan"

rule lex_floats = parse 
    eol { EOL }
  | wspace { lex_floats lexbuf }
  | float as x { FLOAT(float_of_string x) }
  | inf { FLOAT(infinity) }
  | neg_inf { FLOAT(neg_infinity) }
  | nan { FLOAT( 0.0 /. 0.0 ) }
  | eof { EOF }
  | _ { lex_floats lexbuf }
