(* config.ml *)
(* module to manage the parameters input to the program, and to display help *)

   

(***********************************************************)
(* Structure for the parameters of the trmerge program     *)
(***********************************************************)
type 'a context_t =
	{ 
	  mutable file: string;	          	(* input file *)
 	  mutable maxnbseq: int;                (* maximum number of genomic sequences in the input files *)
	  mutable strmode: int;                 (* whether and how strand is taken into account *)
	  mutable outfile: string;              (* output file *)
	  mutable verbose: bool;                (* if we want the output to be displayed on stdout *)
	  mutable ucsc: bool;                   (* ucsc option, prints values of the 9th field with two double quotes 
						   and semicolon to output in ucsc *)
	  mutable sorted: bool;                 (* boolean for whether file is already sorted according to strand, chr, start, end 
						   (sort -k7,7 -k1,1 -k4,4n -k5,5n) *)
	} 



(********************************************************)
(* trmerge context, these are the default parameters    *)
(********************************************************)
let context = 
  {	
    file = "";
    maxnbseq=500;
    strmode=0; 
    outfile="";
    verbose=false;
    ucsc=false;
    sorted=false;
  };;


let usage = 
" 
                     ********************************************          
                     *   compmerge - version v1.4 (June 2018)   *
                     *             Sarah Djebali                *
                     ********************************************

Usage: "^(Filename.basename (Sys.argv.(0)))^" file [options] 

Takes as input the gtf file and merges its transcripts, based on their intron structure
if they are spliced, and based on a simple stranded overlap otherwise.
Outputs a gtf file of merged transcripts. Only works with stranded transcripts.
More precisely, this program merges spliced transcripts that:
- overlap on the same strand by at least one bp,
- have all their overlapping introns that are actually the same,
- have no bp which are exonic in one and intronic in the other,
- have at least one intronic bp in common.

** file must be provided in gtf format (with gene_id before transcript_id).
** [options] can be:
   -o outfile:    outfile is the name of the gtf file the user wants the output to be provided in. 
                  -> default is file_merged.gtf.

   -v:            provides the output result in the standard output rather than in an output file.
                  -> default is unset.
   
   -ucsc:         format the output file in a way that complies with the ucsc browser 
		  (in order to directly load the file in the ucsc browser)
		  (namely adds two double quotes and one semi-colon to each value of a (key,value) pair).
                  -> default is unset.

   -s nbseq:      nbseq is an upper bound for the number of sequences you have in your input gff files.
                  -> default is 500.

   -so:           compmerge does not require the input file to be sorted, however this option enables to skip
                  the input file sorting in case this file is already sorted according to strand, chromosome, start, end. 
                  Note: this sorting could be performed outside compmerge using the unix sort command: 
                  sort -k7,7 -k1,1 -k4,4n -k5,5n file > sortedfile
** Please report any bug to sarahqd@gmail.com        
"




(***********************************************************************)
(* Read the arguments from the command line and updates the context    *)
(***********************************************************************)
let read_commandline () =
  let u = try 
      context.file <- Sys.argv.(1);
    with
      | Invalid_argument s -> Common.print_error usage
  in
  
  (* we start reading the arguments from the 3rd one since the two first are compulsory and are the input files *)
  let argnum = ref 1 and ok = ref true in

  (* This function returns the next argument. mustbe says whether this argument has to exist.
     In general mustbe is true since the arguments go by pairs *)
  let getarg mustbe =
    incr argnum; 
    if !argnum < (Array.length Sys.argv) then 
      Sys.argv.(!argnum)
    else if mustbe then 
      raise Not_found 
    else
      ""
  in
    (* Actually reading each of the arguments *)
    try 
      while (!ok) do
	match (getarg false) with
	  | "-st" -> context.strmode <- int_of_string (getarg true)  
	  | "-o"  -> context.outfile <- getarg true
	  | "-v"  -> context.verbose <- true
	  | "-s"  -> context.maxnbseq <- int_of_string (getarg true)
	  | "-ucsc" -> context.ucsc <- true
	  | "-so" -> context.sorted <- true
	  | "-h"
	  | "--help" -> Printf.printf "%s\n" usage; exit 1; 
	  | ""	-> ok := false
	  | s	-> failwith s
      done;
      Common.print_log "# Command line has been read\n";
    with
      | Not_found -> Common.print_error ("Missing parameter at the end of the command line\n"^usage^"\n");
      | Failure "int_of_string" -> Common.print_error ("Syntax error (incorrect integer value)\n"^usage^"\n");
      | Failure s -> Common.print_error ("Syntax error ("^s^")\n"^usage^"\n");;


