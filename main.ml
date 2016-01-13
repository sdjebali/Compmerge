(* main.ml *)

open Common
open Collection
open Config
open Feature
open Input
open Printf



(* This function takes as input an array of mergeable spliced transcripts ordered from 5' to 3', 
   and outputs a list of list of transcripts to merge later on in a single object (see main function).
   The main rules for compatibility of two spliced transcripts to be merged are:
   - all overlapping introns of the two tr are the same
   AND
   -there is no exon of one tr that overalps the intron of the other one
   Those two rules are implemented in the transcript module in feature.ml.
*)
let comp_and_merge trarr =
  let lgarr = Array.length trarr in
    if (lgarr==0) then
      []
    else
      let i = ref 1 and llindices_tr_to_merge = ref [] in  (* llindices_tr_to_merge represents the list of list of indices in a of the tr to be merged *)
      let h = Hashtbl.create lgarr and found = ref false in                  (* h is a hashtable where each key corresponds to a final list of indices of tr to be merged, 
						       in order words to a final merged tr *)
	Hashtbl.add h 0 0;  (* add then first indice of the array, corresponding to the first transcript, here exists since lgarr not null *)
	while (!i<lgarr) do
	  found:= false;
	  Hashtbl.iter      (* iter over all the keys = over all the list of indices of tr to be merged*)
	    (fun k v -> 
	      if ((not (!found)) && (List.for_all (Transcript.compatible_spliced_tr (trarr.(!i))) (List.map (fun i -> trarr.(i)) (Hashtbl.find_all h k)))) then 
		begin
		  Hashtbl.add h k (!i);
		  found:=true;
		end 
	    ) 
	    h;
	  if (not (!found)) then
	    Hashtbl.add h (!i) (!i);
	  incr i;
	done;
	Hashtbl.iter
	  (fun k v -> let elt_to_add = Hashtbl.find_all h k in 
			if (((!llindices_tr_to_merge)=[]) || ((List.hd (!llindices_tr_to_merge))<>elt_to_add)) then
			  llindices_tr_to_merge:=(Hashtbl.find_all h k)::(!llindices_tr_to_merge)
	  )
	  h;
	List.map (List.map (fun i -> trarr.(i))) (!llindices_tr_to_merge);;



(**********************************************  MAIN FUNCTION *************************************)

(* this is the function that executes the outest actions *)
let comp_and_merge_main () =
  (* Read the arguments on the command line *)
  read_commandline ();
  let u = Common.print_log (("# Input file is ")^(context.file)^("\n")) in
  
  (* Sort input file according to strand and then transcript, and put the result in an intermediate file 
     tmp_sorted_file placed by default in the /tmp directory since a /tmp directory always 
     exists on all systems, but it may be interesting to add an option to enable the user
     to specify his own tmp directory (in case of huge file for example). This argument
     would then need to be passed to the make_temp_file_name function *)
  let tmp_sorted_file = Common.make_sorted_temp_file_if_user_wants_stranded true context.sorted context.file in
  let u= Common.print_log ("# I have treated (sorted and put in temp file) the input file according to what the user wants\n") in

   (* Open input channel *)
  let inchan = open_in tmp_sorted_file in
  
  (* Open output channel *)
  let outchan = 
    if (context.verbose) then
      stdout
    else
      begin
	try
	  (open_out context.outfile) 
	with
	  | Sys_error s -> open_out ((context.file)^("_merged.gtf"))
      end
  in
 

  (* 1. Read from input channel and make list of Transcript objects containing themselves exon objects in their exon lists.
     note: those tr are already ordered according to strand first and then from 5' to 3' and have their exons already 
     ordered from 5' to 3'. *)
  let trlist = read_gff_into_tr_list inchan in
  let u = Common.print_log (("# I have ")^(string_of_int (List.length trlist))^(" initial transcripts\n")) in
  

  (* 2. Divide transcripts into spliced and unspliced and make segseqs out of those.
     note that this division does not affect the tr ordering meaning that within each segseq
     tr are still ordered according to strand first and then from 5' to 3' *)
  let monoextr_list = List.filter Transcript.monoex trlist and splicedtr_list = List.filter Transcript.splicedtr trlist in
  let u = Common.print_log (("# I have ")^(string_of_int (List.length monoextr_list))^(" mono-exonic transcripts and ")^(string_of_int (List.length splicedtr_list))^(" spliced transcripts\n")) in
  let monoex_segseq = SegSeq.make2 (Array.of_list monoextr_list) and splicedtr_segseq = SegSeq.make2 (Array.of_list splicedtr_list) in
  

  (* 3. Make stranded clusters of monoex transcripts, using simple overlap rule but that depends on strand, 
     (which explains why we need to cut by strand first). *)
  let monoex_segseq_cut_by_strand = cut_according_to_strand_tr monoex_segseq in    (* cuts according to strand, modifying the input segseq *)
  let larr_monoex = SegSeq.elements monoex_segseq_cut_by_strand in                 (* list of arrays where each array corresponds to a strand *)
  let lsegseq_monoex = List.map (fun a -> cut_according_to_position_tr (SegSeq.make2 a)) larr_monoex in (* for each strand we have a segseq with the positions 
													   where to cut to get stranded monoex tr clusters *)
  let u = Common.print_log (("# I have ")^(string_of_int (List.length lsegseq_monoex))^(" strand(s) for the monoex tr after their comparison and merging\n")) in

  
  (* 4. Make stranded clusters of spliced transcripts, using simple overlap rule but that depends on strand, 
     (which explains why we need to cut by strand first), and then apply the comp_and_merge algorithm to each of 
     those clusters. The three first steps are the same as for the monoex tr but then need to apply comp_and_merge. *)
  let splicedtr_segseq_cut_by_strand = cut_according_to_strand_tr splicedtr_segseq in        (* cuts according to strand, modifying the input segseq *)
  let larr_splicedtr = SegSeq.elements splicedtr_segseq_cut_by_strand in      (* list of arrays where each array corresponds to a strand *)
  let u = Common.print_log (("# I have ")^(string_of_int (List.length larr_splicedtr))^(" strand(s) for the spliced tr after cutting by strand\n")) in

  let lsegseq_splicedtr = List.map (fun a -> cut_according_to_position_tr (SegSeq.make2 a)) larr_splicedtr in (* for each strand makes a different cut 
														 according to position, and returns a segseq *)
  let llarr_splicedtr = List.map SegSeq.elements lsegseq_splicedtr in  (* a list of list of arrays, for each strand we have a list of arrays 
									  where each array is an array of transcripts that are mergeable and
									  should be submitted to the comp_and_merge algorithm *)
  let u = Common.print_log (("# I have ")^(string_of_int (List.length llarr_splicedtr))^(" strand(s) for the spliced tr after cutting by strand and position\n")) in
 

  let larr_splicedtr = List.flatten llarr_splicedtr in   (* a list of arrays of transcripts, where each array is a set of tr that are mergeable *)
  let lll_splicedtr = List.map comp_and_merge larr_splicedtr in  (* a list of list of list of transcripts where the outer list has as many elements as mergeable clusters
								    and where the inside lists contain the actual tr to merge *)
 
  let u = Common.print_log (("# I have performed the comp_and_merge action for spliced transcripts\n")) in
  let ll_splicedtr = List.flatten lll_splicedtr in  (* a list of list of transcripts where the inner lists gather the tr to merge, this is the final interesting
						       list for spliced tr *)           
  (* let u = Common.print_log (("# I have ")^(string_of_int (List.length ll_splicedtr))^(" final merged spliced transcripts\n")) in *)
           
  (* 5. Output merged monoex and spliced transcripts in output file, providing new ids and recording initial tr ids in the merged tr *)
  let aarr_for_merged_monoextr = Array.of_list (List.flatten (List.map SegSeq.elements lsegseq_monoex)) in
  let alist_for_merged_splicedtr = Array.of_list ll_splicedtr in 
  let u = Common.print_log (("# I made the two arrays to be printed\n")) in
  
				
  let lmergedtr_monoex = ref [] and lmergedtr_spliced = ref [] in
  let u1 = Array.iteri (fun i a -> (lmergedtr_monoex:=(arrtr_to_mergedtr i a)::(!lmergedtr_monoex))) aarr_for_merged_monoextr in
  let nbmonoextr = List.length (!lmergedtr_monoex) in
  let u = Common.print_log (("# I made the actual merged transcript objects to be printed for the mono-exonic class\n")) in
  let u2 = Array.iteri (fun i a -> (lmergedtr_spliced:=(listtr_to_mergedtr (nbmonoextr+i) a)::(!lmergedtr_spliced))) alist_for_merged_splicedtr in
  let u = Common.print_log (("# I made the actual merged transcript objects to be printed for the spliced class\n")) in
  
    List.iter (Transcript.print_gtf outchan context.ucsc) (!lmergedtr_monoex);
    List.iter (Transcript.print_gtf outchan context.ucsc) (!lmergedtr_spliced);
    let u = Common.remove_sorted_temp_file_if_user_wants context.sorted tmp_sorted_file in
    let u = Common.print_log ("# I have removed the temporary sorted file if it has been created\n") in
    let u = Common.print_log ("# compmerge did its work !\n") in
      flush outchan;;   (* note: tr and exons are printed in the reverse order but is it a problem? *)


comp_and_merge_main ();;

