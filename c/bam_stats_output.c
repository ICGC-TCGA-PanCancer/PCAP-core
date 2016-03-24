/*       LICENCE
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
*/

#include <libgen.h>
#include <inttypes.h>
#include "bam_stats_output.h"
#include "bam_stats_calcs.h"

static char *bas_header = "bam_filename\tsample\tplatform\tplatform_unit\tlibrary\treadgroup\tread_length_r1\tread_length_r2\t#_mapped_bases\t#_mapped_bases_r1\t#_mapped_bases_r2\t#_divergent_bases\t#_divergent_bases_r1\t#_divergent_bases_r2\t#_total_reads\t#_total_reads_r1\t#_total_reads_r2\t#_mapped_reads\t#_mapped_reads_r1\t#_mapped_reads_r2\t#_mapped_reads_properly_paired\t#_gc_bases_r1\t#_gc_bases_r2\tmean_insert_size\tinsert_size_sd\tmedian_insert_size\t#_duplicate_reads\t#_mapped_pairs\t#_inter_chr_pairs\n";
static char *rg_line_pattern = "%s\t%s\t%s\t%s\t%s\t%s\t%"PRIu32"\t%"PRIu32"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%.3f\t%.3f\t%.3f\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n";


int bam_stats_output_print_results(rg_info_t **grps,int grps_size,stats_rd_t*** grp_stats,char *input_file,char *output_file){
  FILE *out = NULL;
  check(output_file != NULL, "Output file was NULL");
  if (strcmp(output_file,"-")==0) {
    out = stdout;
  } else {
    out = fopen(output_file,"w");
  }
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
    uint64_t mapped_pairs = 0;
    uint64_t inter_chr_pairs = 0;

    double mean_insert_size = 0;
    double insert_size_sd = 0;
    double median_insert_size = 0;

    if(mapped_reads>0){
      //Only need group one as they are pairs
      proper_pairs = grp_stats[i][0]->proper;
      mapped_pairs = grp_stats[i][0]->mapped_pairs;
      inter_chr_pairs = grp_stats[i][0]->inter_chr_pairs;
      mapped_reads_r1 = total_reads_r1 - unmapped_r1;
      mapped_reads_r2 = total_reads_r2 - unmapped_r2;

      mapped_bases_r1 =  grp_stats[i][0]->mapped_bases;
      mapped_bases_r2 =  grp_stats[i][1]->mapped_bases;
      mapped_bases =  grp_stats[i][0]->mapped_bases +  grp_stats[i][1]->mapped_bases;

      divergent_bases_r1 = grp_stats[i][0]->divergent;
      divergent_bases_r2 = grp_stats[i][1]->divergent;
      divergent_bases = grp_stats[i][0]->divergent + grp_stats[i][1]->divergent;

      bam_stats_calcs_calculate_mean_sd_median_insert_size(grp_stats[i][0]->inserts,&mean_insert_size,&insert_size_sd,&median_insert_size);
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
                      dup_reads,
                      mapped_pairs,
                      inter_chr_pairs);

      check(chk>0,"Error writing bas line to output file.");
      fflush(out);
  }

  if (out != stdout) fclose(out);
  return 0;

error:
  if(out && out != stdout) fclose(out);
  return -1;

}
