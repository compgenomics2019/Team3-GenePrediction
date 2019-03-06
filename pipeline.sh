#!/usr/bin/bash

get_input() {
    # Function to parse arguments
    # Specifying usage message
    usage="Usage: pipeline.bash -i <input directory> -o <output directory> -[v]
              Listeria gene prediction pipeline. The options available are:
                        -i : Directory for genome sequences [required]
                        -o : Output directory [required]
                        -h : Print usage instructions"

  # Getopts block, will take in the arguments as inputs and assign them to variables
        while getopts "i:o:vh" option; do
                case $option in
                        i) input_directory=$OPTARG;;
                        o) output_directory=$OPTARG;;
                        v) v=1;;
                        h) echo "$usage"
                              exit 0;;
                       \?) echo "Invalid option."
                          "$usage"
                              exit 1;;
                esac
        done

  #Check for presence of required arguments
  if [ ! "$input_directory" ] || [ ! "$output_directory" ]
  then
    echo "ERROR: Required arguments missing!"
    echo "$usage"
    exit 1
  fi

  #Check if input directory is a directory

  if [ ! -d $input_directory ]
  then echo "ERROR: Not a valid directory"
  echo "$usage"
  exit 1
  fi

  #Check if output file is already present, give option to rewrite.
	if [ -d $output_directory ]
        then
		echo "Output directory already exists, would you like to overwrite? Reply with y/n"
		read answer
		case $answer in
			y) echo "Overwriting folder $output in subsequent steps";;
			n) echo "Folder overwrite denied, exiting pipeline"
				exit 1;;
			\?) echo "Incorrect option specified, exiting pipeline"
				exit 1;;
		esac
	fi

}

make_temp(){

  #Make temp Directory
  mkdir -p temp

  #Make output directory
  mkdir -p $output_directory

  #Parse input directory to get list of genomes
  ls $input_directory | xargs -L1 bash -c 'a=${0%.*};ext=${0#*.};echo $a >> temp/genomes_list.txt;echo $ext >>temp/genomes_list.txt'

  #Export input directory and to be used within xargs commands
  export input_directory

}

run_prodigal(){

  #Make directory within temp for prodigal genes
  mkdir temp/prodigal_results

  #Run prodigal on genes
  cat temp/genomes_list.txt | xargs -L2 bash -c 'prodigal -i $input_directory/"$0"."$1" -o temp/prodigal_results/"$0".gff -f gff'
}

run_genemark(){

  #Make directory within temp for genemark genes
  mkdir temp/genemark_results

  #Run genemark on genes
  cat temp/genomes_list.txt | xargs -L2 bash -c 'gms2.pl --seq $input_directory/"$0"."$1" --genome-type bacteria --output temp/genemark_results/"$0".gff --fnn temp/genemark_results/$0_nucl.fna --faa  temp/genemark_results/$0_prot.faa --format gff '

  #Move all results to ouput Directory
  mv temp/genemark_results/* $output_directory

}

run_rna (){

  #Make directory within temp for merged results
  mkdir temp/final_results

  #Make directory within temp for trascan genemark_results
  mkdir temp/aragorn

  #Run trnascan on genes
  cat temp/genomes_list.txt | xargs -L2 bash -c 'aragorn -fo -o temp/aragorn/"$0"_aragorn.fna $input_directory/"$0"."$1"'

  #Make directory within temp for rnammer
  mkdir temp/rnammer

  #Make
  cat temp/genomes_list.txt | xargs -L2 bash -c 'rnammer -S bac -m lsu,ssu,tsu -multi -gff temp/rnammer/"$0"_rnammer.gff -f temp/rnammer/"$0"_rnammer.fna < $input_directory/"$0"."$1"'

  cat temp/genomes_list.txt | xargs -L2 bash -c 'cat temp/rnammer/"$0"_rnammer.fna temp/aragorn/"$0"_aragorn.fna > temp/final_results/"$0"_ncRNA.fna'

  #Move RNA results to final directory
  mv temp/final_results/* $output_directory/
}

run_bedtools(){

  #Make directory within temp for merged results
  mkdir temp/merged_results

  #Find the genes that overlap between the prodigal and genemark predictions (with a minimum fraction of overlap of 0.9)
  cat temp/genomes_list.txt | xargs -L2 bash -c 'bedtools intersect -a temp/prodigal_results/"$0".gff -b temp/genemark_results/"$0".gff -r -f 0.90  > temp/merged_results/"$0"_prodigal_genemark.bed'

  #Find the genes predicted by prodigal that are not in genemark
  cat temp/genomes_list.txt | xargs -L2 bash -c 'bedtools intersect -a temp/prodigal_results/"$0".gff -b temp/genemark_results/"$0".gff -v  > temp/merged_results/"$0"_prodigal_only.bed'

  #Find the genes predicted by genemark that are not in prodigal
  cat temp/genomes_list.txt | xargs -L2 bash -c 'bedtools intersect -a temp/genemark_results/"$0".gff -b temp/prodigal_results/"$0".gff -v  > temp/merged_results/"$0"_genemark_only.bed'

  #Merge all 3 DNA prediction files (prodigal_only, genemark_only, and prodigal_genemark)
  cat temp/genomes_list.txt | xargs -L2 bash -c 'cat temp/merged_results/"$0"_prodigal_genemark.bed temp/merged_results/"$0"_prodigal_only.bed temp/merged_results/"$0"_genemark_only.bed | sort -n -k 4 > temp/final_results/"$0"_cds.bed'

  #Get fasta file for the coding regions
  cat temp/genomes_list.txt | xargs -L2 bash -c 'bedtools getfasta -fi $input_directory/"$0"."$1" -bed temp/final_results/"$0"_cds.bed -fo temp/final_results/"$0"_cds.fna'

  #Move all results to temp directory
  mv temp/final_results/* $output_directory
}

main(){
  get_input "$@"
  make_temp
  #run_prodigal "$input_directory"
  run_genemark "$input_directory"
  run_rna "$input_directory"

  #run_bedtools
}

main "$@"
