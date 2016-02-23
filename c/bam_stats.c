/*########LICENCE#########
* PCAP - NGS reference implementations and helper code for the ICGC/TCGA Pan-Cancer Analysis Project
* Copyright (C) 2014-2016 ICGC PanCancer Project
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not see:
*   http://www.gnu.org/licenses/gpl-2.0.html
*#########LICENCE#########*/


#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include "dbg.h"
#include "bam_access.h"
#include "bam_stats_output.h"

#include "khash.h"

static char *input_file = NULL;
static char *output_file = NULL;
static char *ref_file = NULL;
static int rna = 0;
int grps_size = 0;
stats_rd_t*** grp_stats;

int check_exist(char *fname){
	FILE *fp;
	if((fp = fopen(fname,"r"))){
		fclose(fp);
		return 1;
	}
	return 0;
}

void print_version (int exit_code){
  printf ("%s\n",VERSION);
	exit(exit_code);
}

void print_usage (int exit_code){

	printf ("Usage: bam_stats -i file -o file [-p plots] [-r reference.fa.fai] [-h] [-v]\n\n");
  printf ("-i --input     File path to read in.\n");
  printf ("-o --output    File path to output.\n\n");
	printf ("Optional:\n");
	printf ("-r --ref-file  File path to reference index (.fai) file.\n");
	printf ("               NB. If cram format is supplied via -b and the reference listed in the cram header can't be found bam_stats may fail to work correctly.\n");
	printf ("-a --rna       Uses the RNA method of calculating insert size (ignores anything outside ± ('sd'*standard_dev) of the mean in calculating a new mean)\n");

	printf ("Other:\n");
	printf ("-h --help      Display this usage information.\n");
	printf ("-v --version   Prints the version number.\n\n");
  exit(exit_code);
}

void options(int argc, char *argv[]){

  ref_file = NULL;

	const struct option long_opts[] =
	  {
             	{"version", no_argument, 0, 'v'},
             	{"help",no_argument,0,'h'},
              {"input",required_argument,0,'i'},
              {"ref-file",required_argument,0,'r'},
              {"output",required_argument,0,'o'},
              {"rna",no_argument,0, 'a'},
              { NULL, 0, NULL, 0}

   }; //End of declaring opts

   int index = 0;
   int iarg = 0;

   //Iterate through options
   while((iarg = getopt_long(argc, argv, "i:o:r:vha", long_opts, &index)) != -1){
   	switch(iarg){
   		case 'i':
        input_file = optarg;
        break;

   		case 'o':
				output_file = optarg;
   			break;

   		case 'r':
   		  ref_file = optarg;
   		  break;

   		case 'a':
        rna = 1;
        break;

   		case 'h':
        print_usage(0);
        break;

      case 'v':
        print_version(0);
        break;

			case '?':
        print_usage (1);
        break;

      default:
      	print_usage (1);

   	}; // End of args switch statement

   }//End of iteration through options

   //Do some checking to ensure required arguments were passed and are accessible files
   if (input_file==NULL || strcmp(input_file,"/dev/stdin")==0) {
    input_file = "-";   // htslib recognises this as a special case
   }
   if (strcmp(input_file,"-") != 0) {
     if(check_exist(input_file) != 1){
   	  printf("Input file (-i) %s does not exist.\n",input_file);
   	  print_usage(1);
     }
   }
   if (output_file==NULL || strcmp(output_file,"/dev/stdout")==0) {
    output_file = "-";   // we recognise this as a special case
   }
   if(ref_file){
     if(check_exist(ref_file) != 1){
      printf("Reference fasta index file (-r) %s does not exist.\n",ref_file);
      print_usage(1);
     }
   }

   return;
}

int main(int argc, char *argv[]){
	options(argc, argv);
	htsFile *input = NULL;
	bam_hdr_t *head = NULL;
  rg_info_t **grps = NULL;
  //Open bam file as object
  input = hts_open(input_file,"r");
  check(input != NULL, "Error opening hts file for reading '%s'.",input_file);

  //Set reference index file
  if(ref_file){
    hts_set_fai_filename(input, ref_file);
  }else{
    if(input->format.format == cram) log_warn("No reference file provided for a cram input file, if the reference described in the cram header can't be located bam_stats may fail.");
  }

  //Read header from bam file
  head = sam_hdr_read(input);
  check(head != NULL, "Error reading header from opened hts file '%s'.",input_file);


  grps = bam_access_parse_header(head, &grps_size, &grp_stats);
  check(grps != NULL, "Error fetching read groups from header.");

  //Process every read in bam file.
  int check = bam_access_process_reads(input,head,grps, grps_size, &grp_stats, rna);
  check(check==0,"Error processing reads in bam file.");

  int res = bam_stats_output_print_results(grps,grps_size,grp_stats,input_file,output_file);
  check(res==0,"Error writing bam_stats output to file.");

  bam_hdr_destroy(head);
  hts_close(input);

  return 0;

  error:
    if(grps) free(grps);
    if(head) bam_hdr_destroy(head);
    if(input) hts_close(input);
    return 1;
}

