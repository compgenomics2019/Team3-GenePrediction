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
        while getopts "i:o:h" option; do
                case $option in
                        i) input_directory=$OPTARG;;
                        o) output_directory=$OPTARG;;
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

  #Parse input directory to get list of genomes
  ls $input_directory | xargs -L1 bash -c 'a=${0%.*};ext=${0#*.};echo $a >> temp/genomes_list.txt;echo $ext >>temp/genomes_list.txt'

  #Export input directory and to be used within xargs commands
  export input_directory

  #Make output directory
  mkdir -p $output_directory

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

  cat temp/genomes_list.txt | xargs -L2 bash -c 'cat temp/rnammer/"$0"_rnammer.fna temp/aragorn/"$0"_aragorn.fna > temp/final_results/"$0"_RNA.fna'

  #Move RNA results to final directory
  mv temp/final_results/* $output_directory/
}

main(){
  get_input "$@"
  make_temp
  run_genemark "$input_directory"
  run_rna "$input_directory"
}

main "$@"

