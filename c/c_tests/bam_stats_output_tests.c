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

#include	<unistd.h>
#include "minunit.h"
#include "bam_stats_output.h"

char err[100];
char *exp_file = "../t/data/Stats.c.bam.bas";
char *input_file = "../t/data/Stats.bam";

int compare_files(char *file1, char *file2){
  FILE *fp1;
  FILE *fp2;
  char c1[500], c2[500];
  int cmp;
  fp1 = fopen(file1, "r");
  fp2 = fopen(file2, "r");

  if((fp1 == NULL) || (fp2 == NULL)){
    return -1;
  }else{
    while((fgets(c1 , 500, fp1) != NULL) && (fgets(c2 , 500, fp2) != NULL)){
      if((cmp = strcmp(c1, c2)) != 0){
        fprintf(stderr,"Expected: '%s'\nGot: '%s'",c1,c2);
        fclose(fp1);
        fclose(fp2);
        return -1;
      }
    }
  }
  fclose(fp1);
  fclose(fp2);
  return 0;
}

char *bam_stats_output_print_results_test_file(){
  rg_info_t **grps = NULL;
  int grps_size;
  stats_rd_t*** grp_stats;
  char *output_file = NULL;
  int res = bam_stats_output_print_results(grps, grps_size,grp_stats,input_file,output_file);
  if(res != -1){
    sprintf(err,"Should have encountered an error for bad output file.");
    return err;
  }
  output_file = "/cantwritehere";
  res = bam_stats_output_print_results(grps, grps_size,grp_stats,input_file,output_file);
  if(res != -1){
    sprintf(err,"Should have encountered an error for bad output file.");
    return err;
  }


  /****Check for output to file****/
  htsFile *input = NULL;
	bam_hdr_t *head = NULL;
  //Open bam file as object
  input = hts_open(input_file,"r");
  if(input==NULL){
    sprintf(err,"Error opening hts file for reading '%s'.\n",input_file);
    return err;
  }

  //Read header from bam file
  head = sam_hdr_read(input);
  if(head == NULL){
    sprintf(err,"Error reading header from opened hts file '%s'.\n",input_file);
    return err;
  }

  grps = bam_access_parse_header(head, &grps_size, &grp_stats);
  if(grps == NULL){
    sprintf(err,"Error fetching read groups from header.\n");
    return err;
  }
  //Process every read in bam file.
  int check = bam_access_process_reads(input,head,grps, grps_size, &grp_stats, 0);
  if(check!=0){
    sprintf(err,"Error processing reads in bam file.\n");
    return err;
  }
  output_file  = "../t/data/test_out.bam.bas";
  res = bam_stats_output_print_results(grps, grps_size,grp_stats,input_file,output_file);

  //Check test file is equal to expected
  if(compare_files(exp_file,output_file) != 0){
    sprintf(err,"Two files expected %s && got %s were not equal in content\n",exp_file,output_file);
  }

  //Delete test file
  int un = unlink(output_file);
  if(un != 0){
    sprintf(err,"Failed to delete tmp output file %s.\n",output_file);
  }
  return NULL;
}

char *bam_stats_output_print_results_test_stdout(){
  rg_info_t **grps = NULL;
  int grps_size;
  stats_rd_t*** grp_stats;
  char *output_file = NULL;
  int res = bam_stats_output_print_results(grps, grps_size,grp_stats,input_file,output_file);
  if(res != -1){
    sprintf(err,"Should have encountered an error for bad output file.");
    return err;
  }
  output_file = "/cantwritehere";
  res = bam_stats_output_print_results(grps, grps_size,grp_stats,input_file,output_file);
  if(res != -1){
    sprintf(err,"Should have encountered an error for bad output file.");
    return err;
  }


  /****Check for output to file****/
  htsFile *input = NULL;
	bam_hdr_t *head = NULL;
  //Open bam file as object
  input = hts_open(input_file,"r");
  if(input==NULL){
    sprintf(err,"Error opening hts file for reading '%s'.\n",input_file);
    return err;
  }

  //Read header from bam file
  head = sam_hdr_read(input);
  if(head == NULL){
    sprintf(err,"Error reading header from opened hts file '%s'.\n",input_file);
    return err;
  }

  grps = bam_access_parse_header(head, &grps_size, &grp_stats);
  if(grps == NULL){
    sprintf(err,"Error fetching read groups from header.\n");
    return err;
  }
  //Process every read in bam file.
  int check = bam_access_process_reads(input,head,grps, grps_size, &grp_stats, 0);
  if(check!=0){
    sprintf(err,"Error processing reads in bam file.\n");
    return err;
  }
  output_file  = "../t/data/test_out.bam.bas";

  /****Check for output to stdout****/
  //redirect stdout to file
  FILE *frp = freopen(output_file, "w", stdout);
  if(frp == NULL){
    sprintf(err,"Error reassigning stdout to file %s\n",output_file);
  }
  res = bam_stats_output_print_results(grps, grps_size,grp_stats,input_file,output_file);
  frp = freopen("/dev/stdout", "w", stdout);
  //Check test file is equal to expected
  if(compare_files(exp_file,output_file) != 0){
    sprintf(err,"Two files expected %s && got %s were not equal in content\n",exp_file,output_file);
  }
  //Delete test file
  int un = unlink(output_file);
  if(un != 0){
    sprintf(err,"Failed to delete tmp output file %s.\n",output_file);
  }
  return NULL;
}

char *all_tests() {
   mu_suite_start();
   mu_run_test(bam_stats_output_print_results_test_file);
   mu_run_test(bam_stats_output_print_results_test_stdout);
   return NULL;
}

RUN_TESTS(all_tests);