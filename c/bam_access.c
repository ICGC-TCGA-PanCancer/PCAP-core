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

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "bam_access.h"
#include "bam_stats_calcs.h"

int get_rg_index_from_rg_store(rg_info_t **grps, char *rg, int grps_size){
  int i=0;
  for(i=0; i<grps_size; i++){
    if (strncmp(rg,grps[i]->id,strlen(rg))==0) return i;
  }
  return -1;
}

void parse_rg_line(char *tmp_line, rg_info_t *group) {
//Now tokenise tmp_line on \t and read in
  char *tag = strtok(tmp_line,"\t");
  assert(strcmp(tag,"@RG")==0);
  group->id = strdup("\0");
  group->sample = strdup("\0");
  group->platform = strdup("\0");
  group->platform_unit = strdup("\0");
  group->lib = strdup("\0");
  tag = strtok(NULL,"\t");
  while(tag != NULL){
    assert(tag[2]==':');
    tag[2]=0;
    char *val = tag+3;
    if (strcmp("ID",tag)==0) group->id = strdup(val);
    if (strlen(val)==0) val = ".";
    if (strcmp("SM",tag)==0) group->sample = strdup(val);
    if (strcmp("PL",tag)==0) group->platform = strdup(val);
    if (strcmp("PU",tag)==0) group->platform_unit = strdup(val);
    if (strcmp("LB",tag)==0) group->lib = strdup(val);
    tag = strtok(NULL,"\t");
  }//End of iterating through tags in this RG tmp_line
  return;
}

rg_info_t **bam_access_parse_header(bam_hdr_t *head, int *grps_size, stats_rd_t ****grp_stats){
  assert(head != NULL);
  char *line = NULL;
  rg_info_t **groups;
  int size = 0;
  char *head_txt = head->text;
  char *head_bac = strdup(head_txt);
  //First pass counts read groups
  line = strtok(head_txt,"\n");
  while(line != NULL){
		//Check for a read group line
		if(strncmp(line,"@RG",3)==0){
      size++;
    }
    line = strtok(NULL,"\n");
  }
  if(size>0){
    //We now have the number of read groups, assign the RG id to each.
    groups = (rg_info_t**) malloc(sizeof(rg_info_t *) * size);
    check_mem(groups);
    char *ptr = NULL;
    line = strtok_r(head_bac,"\n",&ptr);
    int idx = 0;
    while(line != NULL){
      //Check for a read group line
      if(strncmp(line,"@RG",3)==0){
        groups[idx] = (rg_info_t *) malloc(sizeof(rg_info_t));
        check_mem(groups[idx]);
        char * tmp = strdup(line);
        parse_rg_line(tmp,groups[idx]);
        free (tmp);

        check((groups[idx]->id != NULL),"Error recognising ID from RG line. NULL found.");
        check((groups[idx]->id[0]!='\0'),"Error recognising ID from RG line. Empty string.");
        check((groups[idx]->sample != NULL),"Error recognising SM from RG line.");
        if(groups[idx]->sample[0] == '\0') groups[idx]->sample = strdup(".");
        check(groups[idx]->platform != NULL,"Error recognising PL from RG line.");
        if(groups[idx]->platform[0] == '\0') groups[idx]->platform = strdup(".");
        check(groups[idx]->lib != NULL,"Error recognising LB from RG line.");
        if(groups[idx]->lib[0] == '\0') groups[idx]->lib = strdup(".");
        check(groups[idx]->platform_unit != NULL,"Error recognising PU from RG line.");
        if(groups[idx]->platform_unit[0] == '\0') groups[idx]->platform_unit = strdup(".");
        idx++;
      }//End of iteration through header lines.
      line = strtok_r(NULL,"\n",&ptr);
    }
    if(line) free(line);
	}else{ //Deal with a possible lack of @RG lines.
    groups = malloc(sizeof(rg_info_t*) * 1);
    check_mem(groups);
    groups[0]->id = strdup(".");
    groups[0]->sample = strdup(".");
    groups[0]->platform = strdup(".");
    groups[0]->platform_unit = strdup(".");
    groups[0]->lib = strdup(".");
    size = 1;
	}
	*grp_stats = (stats_rd_t***) malloc(sizeof(stats_rd_t**) * (size));
  check_mem(*grp_stats);
  int j=0;
  for(j=0; j<size; j++){
    (*grp_stats)[j] = (stats_rd_t **) malloc((sizeof(stats_rd_t*)*2));
    check_mem((*grp_stats)[j]);
    (*grp_stats)[j][0] = (stats_rd_t *) malloc(sizeof(stats_rd_t));//Setup read one stats store
    check_mem((*grp_stats)[j][0]);
    (*grp_stats)[j][0]->length= 0;
    (*grp_stats)[j][0]->count= 0;
    (*grp_stats)[j][0]->dups= 0;
    (*grp_stats)[j][0]->gc= 0;
    (*grp_stats)[j][0]->umap= 0;
    (*grp_stats)[j][0]->divergent= 0;
    (*grp_stats)[j][0]->mapped_bases= 0;
    (*grp_stats)[j][0]->proper= 0;
    (*grp_stats)[j][0]->inserts = kh_init(ins);
    (*grp_stats)[j][1] = (stats_rd_t *)malloc(sizeof(stats_rd_t));//Setup read two stats store
    check_mem((*grp_stats)[j][1]);
    (*grp_stats)[j][1]->length= 0;
    (*grp_stats)[j][1]->count= 0;
    (*grp_stats)[j][1]->dups= 0;
    (*grp_stats)[j][1]->gc= 0;
    (*grp_stats)[j][1]->umap= 0;
    (*grp_stats)[j][1]->divergent= 0;
    (*grp_stats)[j][1]->mapped_bases= 0;
    (*grp_stats)[j][1]->proper= 0;
  }
  *grps_size = size;
	return groups;

error:
  if(groups) free(groups);
  if(line) free(line);
  if(grp_stats) free(grp_stats);
  return NULL;
}

int bam_access_process_reads(htsFile *input, bam_hdr_t *head, rg_info_t **grps, int grps_size, stats_rd_t ****grp_stats, int rna){
  assert(input != NULL);
  assert(head != NULL);
  assert(grps != NULL);

  bam1_t *b;
  //Iterate through each read in bam file.
  b = bam_init1();
  int ret;
  while((ret = sam_read1(input, head, b)) >= 0){
    if (b->core.flag & BAM_FSECONDARY && rna == 0) continue; //skip secondary hits so no double counts
    if (b->core.flag & BAM_FQCFAIL) continue; // skip vendor fail as generally aren't considered
    if (b->core.flag & BAM_FSUPPLEMENTARY) continue; // skip supplimentary

    uint8_t read = 1; //second read
    if (b->core.flag & BAM_FREAD1) read = 0; //first read

    char *rg = bam_aux2Z(bam_aux_get(b,"RG"));
    if(rg == NULL || strlen(rg)==0){
      rg = ".";
    }

    int rg_index = get_rg_index_from_rg_store(grps,rg,grps_size);
    check(rg_index>=0, "Error assigning @RG ID index for ID:%s.", rg);
    check(rg_index<grps_size, "Error assigning @RG ID index for ID:%s.", rg);

    // grp_stats[rg_index][read]; Stats for this RG/read order combination
    if((*grp_stats)[rg_index][read]->length == 0) (*grp_stats)[rg_index][read]->length = b->core.l_qseq;

    (*grp_stats)[rg_index][read]->count++;
    if(b->core.flag & BAM_FDUP) (*grp_stats)[rg_index][read]->dups++;

    //Get the count of GCs in the sequence.
    int i=0;
    for(i=0;i<b->core.l_qseq;i++){
      uint8_t base = bam_seqi(bam_get_seq(b),i);
      if(base==4||base==2) (*grp_stats)[rg_index][read]->gc++; //Check for G/C
    }

    //Count unmapped and go to next read as anything after this is for mapped only.
    if(b->core.flag & BAM_FUNMAP){
      (*grp_stats)[rg_index][read]->umap++;
      continue;
    }

    // everything after this point must require reads are mapped

    // Divergence calculation: Collect stats that will allow us to calculate the the number of bases that diverge from the reference.
    //                         This requires collecting the value from the NM tag and the mapped proportion of the query string.
    uint8_t *nm = 0;
    nm = bam_aux_get(b,"NM");

    if(nm){
      uint32_t nm_val = bam_aux2i(nm);
      if(nm_val>0){
        (*grp_stats)[rg_index][read]->divergent += nm_val;
        (*grp_stats)[rg_index][read]->mapped_bases += bam_access_get_mapped_base_count_from_cigar(b);
      }else{
        (*grp_stats)[rg_index][read]->mapped_bases += (bam_endpos(b) - b->core.pos) + 1;
      }
    }else{
      (*grp_stats)[rg_index][read]->mapped_bases += bam_access_get_mapped_base_count_from_cigar(b);
    }

    // stats that only assess read 1
    if(b->core.flag & BAM_FREAD1) {
      // Count all the pairs where both ends are not unmapped
      // already tested if this read is mapped above
      if(!(b->core.flag & BAM_FMUNMAP)) {
        (*grp_stats)[rg_index][read]->mapped_pairs++;

        // Insert size can only be calculated based on reads that are on same chr
        // so it is more sensible to generate the distribution based on $PROPER-pairs.
        // only assess read 1 as size is a factor of the pair
        if(b->core.flag & BAM_FPROPER_PAIR){
          (*grp_stats)[rg_index][read]->proper++;
          uint32_t ins = b->core.isize;
          int res;
          khint_t k;
          k = kh_put(ins,(*grp_stats)[rg_index][read]->inserts,abs(ins),&res);
          if(res){
            kh_value((*grp_stats)[rg_index][read]->inserts,k) = 1;
          }else{
            kh_value((*grp_stats)[rg_index][read]->inserts,k) = kh_value((*grp_stats)[rg_index][read]->inserts,k)+1;
          }
        }
        else if(b->core.tid != b->core.mtid) {
          // here count the reads where the chr are different
          (*grp_stats)[rg_index][read]->inter_chr_pairs++;
        }
      }
    }

  }
  bam_destroy1(b);
  return 0;
  error:
    if(b) bam_destroy1(b);
    return -1;
}

uint64_t bam_access_get_mapped_base_count_from_cigar(bam1_t *b){
#define _cop(c) ((c)&BAM_CIGAR_MASK)
#define _cln(c) ((c)>>BAM_CIGAR_SHIFT)
  assert(b != NULL);
  uint64_t count = 0;
  uint32_t *cigar = bam_get_cigar(b);
  //iterate through each cigar operation
  int i=0;
  for(i=0;i<b->core.n_cigar;i++){
    int op = _cop(cigar[i]);
    int l = _cln(cigar[i]);
    if(op == BAM_CINS|| op == BAM_CEQUAL || op == BAM_CMATCH || op == BAM_CDIFF ){ //Match,Insert,Missmatch,Equal all sum to generate the length of mapped bases
      count += l;
    }
  }
  return count;
}

