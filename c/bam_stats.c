/*       LICENCE
* PCAP - NGS reference implementations and helper code for the ICGC/TCGA Pan-Cancer Analysis Project
* Copyright (C) 2014 ICGC PanCancer Project
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
*/


#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include "dbg.h"
#include <bam_access.h>

#include "khash.h"

static char *input_file;
static char *output_file;
static char *ref_file;
int grps_size = 0;
stats_rd_t*** grp_stats;
static char *bas_header = "bam_filename\tsample\tplatform\tplatform_unit\tlibrary\treadgroup\tread_length_r1\tread_length_r2\t#_mapped_bases\t#_mapped_bases_r1\t#_mapped_bases_r2\t#_divergent_bases\t#_divergent_bases_r1\t#_divergent_bases_r2\t#_total_reads\t#_total_reads_r1\t#_total_reads_r2\t#_mapped_reads\t#_mapped_reads_r1\t#_mapped_reads_r2\t#_mapped_reads_properly_paired\t#_gc_bases_r1\t#_gc_bases_r2\tmean_insert_size\tinsert_size_sd\tmedian_insert_size\t#_duplicate_reads\n";
static char *rg_line_pattern = "%s\t%s\t%s\t%s\t%s\t%s\t%ld\t%ld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%.3f\t%.3f\t%.3f\t%lld\n";


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
	return;
}

void print_usage (int exit_code){

	printf ("Usage: bam_stats -i file -o file [-p plots] [-r reference.fa.fai] [-h] [-v]\n\n");
  printf ("-i --input     File path to read in.\n");
  printf ("-o --output    File path to output.\n\n");
	printf ("Optional:\n");
	printf ("-r --ref-file  File path to reference index (.fai) file.\n");
	printf ("               NB. If cram format is supplied via -b and the reference listed in the cram header can't be found bam_stats may fail to work correctly.\n");
	printf ("-p --plots     Folder to contain quality score plots.\n\n");
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
              {"plots",required_argument,0,'p'},
              { NULL, 0, NULL, 0}

   }; //End of declaring opts

   int index = 0;
   int iarg = 0;

   //Iterate through options
   while((iarg = getopt_long(argc, argv, "i:o:r:pvh", long_opts, &index)) != -1){
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

   		case 'p':
   		  print_usage(1);
   			//plots = optarg;
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
   if(check_exist(input_file) != 1){
   	printf("Input file (-i) %s does not exist.\n",input_file);
   	print_usage(1);
   }
   if(ref_file){
     if(check_exist(ref_file) != 1){
      printf("Reference fasta index file (-r) %s does not exist.\n",ref_file);
      print_usage(1);
     }
   }

   return;
}

int calculate_mean_sd_median_insert_size(khash_t(ins) *inserts,double *mean, double *sd, double *median){

    uint64_t pp_mean = 0;
    uint64_t tt_mean = 0;
    uint32_t key;
    uint64_t val;

    kh_foreach(inserts,key,val,
        { pp_mean += key * val;
          tt_mean += val;
        });

    /*
    int i = 0;
    for(i=0;i<200000;i++){
      if(inserts[i]>0){
        pp_mean += (i+1) * inserts[i];
        tt_mean += inserts[i];
      }
    }
    */

    if(tt_mean){//Calculate mean , median, sd
      *mean = (double) ((double)pp_mean/(double)tt_mean);

      float midpoint = (float)((float)tt_mean / (float)2);
      float midpoint2 = (float)((float)tt_mean / (float)2 + (float)1);
      uint32_t insert = 0;
      uint32_t prev_insert = 0;
      uint64_t running_total = 0;

      kh_foreach(inserts,key,val,
            { insert = key;
              running_total += val;
              if(running_total >= midpoint) break;
              prev_insert = key;
            });

      /*
      int j=0;

      for(j=0;j<200000;j++){
        if(inserts[j]>0){
          insert = j+1;
          running_total += inserts[j];
          if(running_total >= midpoint) break ;
          prev_insert = j+1;
        }
      }*/

      if(tt_mean %2 == 0 && ( running_total - midpoint2 >= insert )){
        //warn "Thinks is even AND split between bins ";
        *median = (((double)insert + (double)prev_insert) / (double)2);
      }else{
        //warn "Thinks is odd or NOT split between bins";
        *median = (double)(insert);
      }

      //We have mean and median so calculate the SD
      uint64_t pp_sd = 0;
      uint64_t tt_sd = 0;

      kh_foreach(inserts,key,val,
            { double diff = (double)(key) - (*mean);
              pp_sd += (diff * diff) * (double)val;
              tt_sd += val;
            });

      /*
      int k=0;
      for(k=0;k<200000;k++){
        if(inserts[k]>0){
          double diff = (double)(k+1) - (*mean);
          pp_sd += (diff * diff) * inserts[k];
          tt_sd += inserts[k];
        }
      }*/
      kh_destroy(ins, inserts);

      if(tt_sd){
        double variance = fabs((double)((double)pp_sd / (double)tt_sd));
        *sd = sqrt(variance);
      }else{
        *sd = 0;
      }

    } //End of if we have data to calculate from.
  return 0;
}

int print_results(rg_info_t **grps){
  FILE *out = fopen(output_file,"w");
  check(out != NULL,"Error trying to open output file %s for writing.",output_file);

  int chk = fprintf(out,"%s",bas_header);
  check(chk==strlen(bas_header),"Error writing bas_header to output file.");

  //Iterate through each RG
  int i=0;
  for(i=0;i<grps_size;i++){
    if(grp_stats[i][0]->count==0 && grp_stats[i][1]->count==0) continue; // Skip empty read groups stats
    uint64_t unmapped_r1 = grp_stats[i][0]->umap;
    uint64_t unmapped_r2 = grp_stats[i][1]->umap;
    uint64_t unmapped =  grp_stats[i][0]->umap + grp_stats[i][1]->umap;

    uint64_t total_reads_r1 = grp_stats[i][0]->count;
    uint64_t total_reads_r2 = grp_stats[i][1]->count;
    uint64_t total_reads = grp_stats[i][0]->count + grp_stats[i][1]->count;

    uint64_t mapped_reads = total_reads - unmapped;

    uint64_t gc_r1 = grp_stats[i][0]->gc;
    uint64_t gc_r2 = grp_stats[i][1]->gc;

    uint64_t mapped_reads_r1 = 0;
    uint64_t mapped_reads_r2 = 0;
    uint64_t proper_pairs = 0;
    uint64_t mapped_bases = 0;
    uint64_t mapped_bases_r1 = 0;
    uint64_t mapped_bases_r2 = 0;
    uint64_t divergent_bases = 0;
    uint64_t divergent_bases_r1 = 0;
    uint64_t divergent_bases_r2 = 0;
    uint64_t dup_reads = 0;

    double mean_insert_size = 0;
    double insert_size_sd = 0;
    double median_insert_size = 0;

    if(mapped_reads>0){
      //Only need group one as they are pairs
      proper_pairs = grp_stats[i][0]->proper;
      mapped_reads_r1 = total_reads_r1 - unmapped_r1;
      mapped_reads_r2 = total_reads_r2 - unmapped_r2;

      mapped_bases_r1 =  grp_stats[i][0]->mapped_bases;
      mapped_bases_r2 =  grp_stats[i][1]->mapped_bases;
      mapped_bases =  grp_stats[i][0]->mapped_bases +  grp_stats[i][1]->mapped_bases;

      divergent_bases_r1 = grp_stats[i][0]->divergent;
      divergent_bases_r2 = grp_stats[i][1]->divergent;
      divergent_bases = grp_stats[i][0]->divergent + grp_stats[i][1]->divergent;

      calculate_mean_sd_median_insert_size(grp_stats[i][0]->inserts,&mean_insert_size,&insert_size_sd,&median_insert_size);
      dup_reads = grp_stats[i][0]->dups + grp_stats[i][1]->dups;
    }

    uint32_t read_length_r1 = grp_stats[i][0]->length;
    uint32_t read_length_r2 = grp_stats[i][1]->length;

    char *file = basename(input_file);

    chk = fprintf(out,rg_line_pattern,
                      file,
                      grps[i]->sample,
                      grps[i]->platform,
                      grps[i]->platform_unit,
                      grps[i]->lib,
                      grps[i]->id,
                      read_length_r1,
                      read_length_r2,
                      mapped_bases,
                      mapped_bases_r1,
                      mapped_bases_r2,
                      divergent_bases,
                      divergent_bases_r1,
                      divergent_bases_r2,
                      total_reads,
                      total_reads_r1,
                      total_reads_r2,
                      mapped_reads,
                      mapped_reads_r1,
                      mapped_reads_r2,
                      proper_pairs,
                      gc_r1,
                      gc_r2,
                      mean_insert_size,
                      insert_size_sd,
                      median_insert_size,
                      dup_reads);

      check(chk>0,"Error writing bas line to output file.");
      fflush(out);
  }

  fclose(out);
  return 0;

error:
  if(out) fclose(out);
  return -1;

}

int main(int argc, char *argv[]){
	options(argc, argv);
	htsFile *input;
	bam_hdr_t *head;

  //Open bam file as object
  input = hts_open(input_file,"r");
  //Set reference index file
  if(ref_file){
    hts_set_fai_filename(input, ref_file);
  }else{
    if(input->format.format == cram) log_warn("No reference file provided for a cram input file, if the reference described in the cram header can't be located bam_stats may fail.");
  }

  //Read header from bam file
  head = sam_hdr_read(input);
  rg_info_t **grps = parse_header(head, &grps_size, &grp_stats);
  check(grps != NULL, "Error fetching read groups from header.");

  //Process every read in bam file.
  int check = process_reads(input,head,grps, grps_size, &grp_stats);
  check(check==0,"Error processing reads in bam file.");

  int res = print_results(grps);
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

