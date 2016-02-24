/*########LICENCE#########
* PCAP - NGS reference implementations and helper code for the ICGC/TCGA Pan-Cancer Analysis Project
* Copyright (C) 2014,2015,2016 ICGC PanCancer Project
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

#include <inttypes.h>
#include "minunit.h"
#include "bam_access.h"

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
    size_t l, m;
    char *s;
} kstring_t;
#endif

static inline int kputsn(const char *p, int l, kstring_t *s)
{
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	memcpy(s->s + s->l, p, l);
	s->l += l;
	s->s[s->l] = 0;
	return l;
}


static inline int kputs(const char *p, kstring_t *s)
{
	return kputsn(p, strlen(p), s);
}

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

char *test_bam = "../t/data/Stats.bam";
char *exp_plat = "GAII";
char *exp_sample = "PD1234a";
char *exp_lib = "PD1234a 140546_1054";
char *exp_rgid_1 = "29976";
char *exp_platform_unit_1 = "5178_6";
char *exp_rgid_2 = "29978";
char *exp_platform_unit_2 = "5085_6";
uint32_t exp_rd_length = 20;
uint64_t exp_rg1_rd1_tot_count = 7;
uint64_t exp_rg1_rd2_tot_count = 3;
uint64_t exp_rg2_rd1_tot_count = 3;
uint64_t exp_rg2_rd2_tot_count = 3;
uint64_t exp_rg1_rd1_dups = 4;
uint64_t exp_rg1_rd2_dups = 0;
uint64_t exp_rg2_rd1_dups = 0;
uint64_t exp_rg2_rd2_dups = 0;
uint64_t exp_rg1_rd1_gc = 63;
uint64_t exp_rg1_rd2_gc = 27;
uint64_t exp_rg2_rd1_gc = 27;
uint64_t exp_rg2_rd2_gc = 27;
uint64_t exp_rg1_rd1_umap = 1;
uint64_t exp_rg1_rd2_umap = 1;
uint64_t exp_rg2_rd1_umap = 2;
uint64_t exp_rg2_rd2_umap = 2;
uint64_t exp_rg1_rd1_divergent = 30;
uint64_t exp_rg1_rd2_divergent = 18;
uint64_t exp_rg2_rd1_divergent = 11;
uint64_t exp_rg2_rd2_divergent = 12;
uint64_t exp_rg1_rd1_mapped_bases = 115;
uint64_t exp_rg1_rd2_mapped_bases = 37;
uint64_t exp_rg2_rd1_mapped_bases = 20;
uint64_t exp_rg2_rd2_mapped_bases = 20;
uint64_t exp_rg1_rd1_proper = 6;
uint64_t exp_rg1_rd2_proper = 0;
uint64_t exp_rg2_rd1_proper = 1;
uint64_t exp_rg2_rd2_proper = 0;

char err[100];

char *test_bam_access_parse_header(){
  htsFile *input;
	bam_hdr_t *head;
	int grps_size = 0;
	stats_rd_t*** grp_stats;
  //Open bam file as object
  input = hts_open(test_bam,"r");
  if (input == NULL){
    sprintf(err,"Error opening bam file %s\n",test_bam);
    return err;
  }
  //Read header from bam file
  head = sam_hdr_read(input);
  if (head == NULL){
    sprintf(err,"Error reading header from bam file %s\n",test_bam);
    return err;
  }
  //Run header parsing
  rg_info_t **grps = bam_access_parse_header(head, &grps_size, &grp_stats);
  if(grps_size != 2){
    sprintf(err,"Didn't read two read groups from test bam: %d\n",grps_size);
    return err;
  }
  //Check rg 1
  if(strcmp(grps[0]->id,exp_rgid_1)!=0){
    sprintf(err,"Read group %d expected id %s but got %s.\n",1,exp_rgid_1,grps[0]->id);
    return err;
  }
  if(strcmp(grps[0]->platform_unit,exp_platform_unit_1)!=0){
    sprintf(err,"Read group %d expected platform_unit %s but got %s.\n",2,exp_platform_unit_1,grps[0]->platform_unit);
    return err;
  }
  //Check rg 2
  if(strcmp(grps[1]->id,exp_rgid_2)!=0){
    sprintf(err,"Read group %d expected id %s but got %s.\n",2,exp_rgid_2,grps[1]->id);
    return err;
  }
  if(strcmp(grps[1]->platform_unit,exp_platform_unit_2)!=0){
    sprintf(err,"Read group %d expected platform_unit %s but got %s.\n",2,exp_platform_unit_2,grps[1]->platform_unit);
    return err;
  }

  //Commmon to both RGs
  int i=0;
  for (i=0;i<2;i++){
    if(strcmp(grps[i]->platform,exp_plat)!=0){
      sprintf(err,"Read group %d expected platform %s but got %s.\n",(i+1),exp_plat,grps[i]->platform);
      return err;
    }
    if(strcmp(grps[i]->sample,exp_sample)!=0){
      sprintf(err,"Read group %d expected sample %s but got %s.\n",(i+1),exp_sample,grps[i]->sample);
      return err;
    }
    if(strcmp(grps[i]->lib,exp_lib)!= 0){
      sprintf(err,"Read group %d expected lib %s but got %s.\n",(i+1),exp_lib,grps[i]->lib);
      return err;
    }
  }
	return NULL;
}

char *test_bam_access_get_mapped_base_count_from_cigar(){
  htsFile *input;
  bam_hdr_t *head;
  kstring_t str = {0,0,0};
  char *sample_sam = "IL29_5178:2:54:17473:17010	579	1	9993	0	5S10M5S	=	9993	100	CTCTTCCGATCTTTAGGGTT	;\?;\?\?>>>>F<BBDEBEEFF	RG:Z:29976	NM:i:1";
	kputs(sample_sam,&str);
  //Open bam file as object
  input = hts_open(test_bam,"r");
  if (input == NULL){
    sprintf(err,"Error opening bam file %s\n",test_bam);
    return err;
  }
  //Read header from bam file
  head = sam_hdr_read(input);
  if (head == NULL){
    sprintf(err,"Error reading header from bam file %s\n",test_bam);
    return err;
  }
  //bam_access_get_mapped_base_count_from_cigar();
  bam1_t *b = bam_init1();
  int ret = sam_parse1(&str,head,b);
  if(ret<0){
    sprintf(err,"Error reading sam record 1 from bam file\n");
    return err;
  }
  int count = bam_access_get_mapped_base_count_from_cigar(b);
  if(count != 10){
    sprintf(err,"Error checking cigar. Expected 10 mapped bases got %d\n",count);
    bam_destroy1(b);
    return err;
  }

  kstring_t str_two = {0,0,0};
  char *sample_sam_two = "IL29_5178:2:54:17473:17010	579	1	9993	0	20M	=	9993	100	CTCTTCCGATCTTTAGGGTT	;\?;\?\?>>>>F<BBDEBEEFF	RG:Z:29976	NM:i:1";
  bam_destroy1(b);
  b = bam_init1();
  kputs(sample_sam_two,&str_two);
  ret = sam_parse1(&str_two,head,b);
  if(ret<0){
    sprintf(err,"Error reading sam record 2 from bam file '%d'\n",ret);
    return err;
  }
  count = bam_access_get_mapped_base_count_from_cigar(b);
  if(count != 20){
    sprintf(err,"Error checking cigar. Expected 20 mapped bases got %d\n",count);
    bam_destroy1(b);
    return err;
  }
  bam_destroy1(b);

	return NULL;
}

char *test_bam_access_process_reads_no_rna(){
  htsFile *input;
	bam_hdr_t *head;
	int grps_size = 0;
	stats_rd_t*** grp_stats;
  //Open bam file as object
  input = hts_open(test_bam,"r");
  if (input == NULL){
    sprintf(err,"Error opening bam file %s\n",test_bam);
    return err;
  }
  //Read header from bam file
  head = sam_hdr_read(input);
  if (head == NULL){
    sprintf(err,"Error reading header from bam file %s\n",test_bam);
    return err;
  }
  //Run header parsing
  rg_info_t **grps = bam_access_parse_header(head, &grps_size, &grp_stats);
  if(grps_size != 2){
    sprintf(err,"Didn't read two read groups from test bam: %d\n",grps_size);
    return err;
  }
  int check = bam_access_process_reads(input, head, grps, grps_size, &grp_stats, 0);
  if(check!=0){
    sprintf(err,"Error processing reads in bam file, non RNA.\n");
  }

  //Commmon to both RGs
  int i=0;
  for (i=0;i<2;i++){
    //Read lengths
    if(grp_stats[i][0]->length!=exp_rd_length){
      sprintf(err,"Read group %d read_2 length expected %"PRIu32" but got %"PRIu32".\n",(i+1),exp_rd_length,grp_stats[i][0]->length);
      return err;
    }
    if(grp_stats[i][1]->length!=exp_rd_length){
      sprintf(err,"Read group %d read_2 length expected %"PRIu32" but got %"PRIu32".\n",(i+1),exp_rd_length,grp_stats[i][1]->length);
      return err;
    }
  }

  /******RG1 checks*****/
  //Read count 1
  if(grp_stats[0][0]->count!=exp_rg1_rd1_tot_count){
    sprintf(err,"RG 1, read_1 count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd1_tot_count,grp_stats[0][0]->count);
    return err;
  }
  //Read count 2
  if(grp_stats[0][1]->count!=exp_rg1_rd2_tot_count){
    sprintf(err,"RG 1, read_2 count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd2_tot_count,grp_stats[0][1]->count);
    return err;
  }
  //Duplicate reads 1
  if(grp_stats[0][0]->dups!=exp_rg1_rd1_dups){
    sprintf(err,"RG 1, read_1 duplicate count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd1_dups,grp_stats[0][0]->dups);
    return err;
  }
  //Duplicate reads 2
  if(grp_stats[0][1]->dups!=exp_rg1_rd2_dups){
    sprintf(err,"RG 1, read_2 duplicate count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd2_dups,grp_stats[0][1]->dups);
    return err;
  }
  //GC 1
  if(grp_stats[0][0]->gc!=exp_rg1_rd1_gc){
    sprintf(err,"RG 1, read_1 gc count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd1_gc,grp_stats[0][0]->gc);
    return err;
  }
  //GC 2
  if(grp_stats[0][1]->gc!=exp_rg1_rd2_gc){
    sprintf(err,"RG 1, read_2 gc count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd2_gc,grp_stats[0][1]->gc);
    return err;
  }
  //Unmapped 1
  if(grp_stats[0][0]->umap!=exp_rg1_rd1_umap){
    sprintf(err,"RG 1, read_1 umap count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd1_umap,grp_stats[0][0]->umap);
    return err;
  }
  //Unmapped 2
  if(grp_stats[0][1]->umap!=exp_rg1_rd2_umap){
    sprintf(err,"RG 1, read_2 umap count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd2_umap,grp_stats[0][1]->umap);
    return err;
  }
  //divergent 1
  if(grp_stats[0][0]->divergent!=exp_rg1_rd1_divergent){
    sprintf(err,"RG 1, read_1 divergent count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd1_divergent,grp_stats[0][0]->divergent);
    return err;
  }
  //divergent 2
  if(grp_stats[0][1]->divergent!=exp_rg1_rd2_divergent){
    sprintf(err,"RG 1, read_2 divergent count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd2_divergent,grp_stats[0][1]->divergent);
    return err;
  }
  //mapped_bases 1
  if(grp_stats[0][0]->mapped_bases!=exp_rg1_rd1_mapped_bases){
    sprintf(err,"RG 1, read_1 mapped_bases count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd1_mapped_bases,grp_stats[0][0]->mapped_bases);
    return err;
  }
  //mapped_bases 2
  if(grp_stats[0][1]->mapped_bases!=exp_rg1_rd2_mapped_bases){
    sprintf(err,"RG 1, read_2 mapped_bases count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd2_mapped_bases,grp_stats[0][1]->mapped_bases);
    return err;
  }
  //proper 1
  if(grp_stats[0][0]->proper!=exp_rg1_rd1_proper){
    sprintf(err,"RG 1, read_1 proper count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd1_proper,grp_stats[0][0]->proper);
    return err;
  }
  //proper 2
  if(grp_stats[0][1]->proper!=exp_rg1_rd2_proper){
    sprintf(err,"RG 1, read_2 proper count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd2_proper,grp_stats[0][1]->proper);
    return err;
  }




  /******RG2 checks*****/
  if(grp_stats[1][0]->count!=exp_rg2_rd1_tot_count){
    sprintf(err,"RG 2, read_1 count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd1_tot_count,grp_stats[1][0]->count);
    return err;
  }
  if(grp_stats[1][1]->count!=exp_rg2_rd2_tot_count){
    sprintf(err,"RG 2, read_2 count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd2_tot_count,grp_stats[1][1]->count);
    return err;
  }
  //Duplicate reads 1
  if(grp_stats[1][0]->dups!=exp_rg2_rd1_dups){
    sprintf(err,"RG 1, read_1 duplicate count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd1_dups,grp_stats[1][0]->dups);
    return err;
  }
  //Duplicate reads 2
  if(grp_stats[1][1]->dups!=exp_rg2_rd2_dups){
    sprintf(err,"RG 1, read_2 duplicate count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd2_dups,grp_stats[1][1]->dups);
    return err;
  }
  //GC 1
  if(grp_stats[1][0]->gc!=exp_rg2_rd1_gc){
    sprintf(err,"RG 1, read_1 gc count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd1_gc,grp_stats[1][0]->gc);
    return err;
  }
  //GC 2
  if(grp_stats[1][1]->gc!=exp_rg2_rd2_gc){
    sprintf(err,"RG 1, read_2 gc count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd2_gc,grp_stats[1][1]->gc);
    return err;
  }
  //Unmapped 1
  if(grp_stats[1][0]->umap!=exp_rg2_rd1_umap){
    sprintf(err,"RG 2, read_1 umap count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd1_umap,grp_stats[1][0]->umap);
    return err;
  }
  //Unmapped 2
  if(grp_stats[1][1]->umap!=exp_rg2_rd2_umap){
    sprintf(err,"RG 2, read_2 umap count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd2_umap,grp_stats[1][1]->umap);
    return err;
  }
  //divergent 1
  if(grp_stats[1][0]->divergent!=exp_rg2_rd1_divergent){
    sprintf(err,"RG 2, read_1 divergent count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd1_divergent,grp_stats[1][0]->divergent);
    return err;
  }
  //divergent 2
  if(grp_stats[1][1]->divergent!=exp_rg2_rd2_divergent){
    sprintf(err,"RG 2, read_2 divergent count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd2_divergent,grp_stats[1][1]->divergent);
    return err;
  }
  //mapped_bases 1
  if(grp_stats[1][0]->mapped_bases!=exp_rg2_rd1_mapped_bases){
    sprintf(err,"RG 2, read_1 mapped_bases count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd1_mapped_bases,grp_stats[1][0]->mapped_bases);
    return err;
  }
  //mapped_bases 2
  if(grp_stats[1][1]->mapped_bases!=exp_rg2_rd2_mapped_bases){
    sprintf(err,"RG 2, read_2 mapped_bases count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd2_mapped_bases,grp_stats[1][1]->mapped_bases);
    return err;
  }
  //proper 1
  if(grp_stats[1][0]->proper!=exp_rg2_rd1_proper){
    sprintf(err,"RG 2, read_1 proper count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd1_proper,grp_stats[1][0]->proper);
    return err;
  }
  //proper 2
  if(grp_stats[1][1]->proper!=exp_rg2_rd2_proper){
    sprintf(err,"RG 2, read_2 proper count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd2_proper,grp_stats[1][1]->proper);
    return err;
  }
	return NULL;
}

char *test_bam_access_process_reads_rna(){ // rna flag in this method includes secondary reads
  htsFile *input;
	bam_hdr_t *head;

	int grps_size = 0;
	stats_rd_t*** grp_stats;
  //Open bam file as object
  input = hts_open(test_bam,"r");
  if (input == NULL){
    sprintf(err,"Error opening bam file %s\n",test_bam);
    return err;
  }
  //Read header from bam file
  head = sam_hdr_read(input);
  if (head == NULL){
    sprintf(err,"Error reading header from bam file %s\n",test_bam);
    return err;
  }
  //Run header parsing
  rg_info_t **grps = bam_access_parse_header(head, &grps_size, &grp_stats);
  if(grps_size != 2){
    sprintf(err,"Didn't read two read groups from test bam: %d\n",grps_size);
    return err;
  }
  int check = bam_access_process_reads(input, head, grps, grps_size, &grp_stats, 0);
  if(check!=0){
    sprintf(err,"Error processing reads in bam file, non RNA.\n");
  }

  //Commmon to both RGs
  int i=0;
  for (i=0;i<2;i++){
    //Read lengths
    if(grp_stats[i][0]->length!=exp_rd_length){
      sprintf(err,"Read group %d read_2 length expected %"PRIu32" but got %"PRIu32".\n",(i+1),exp_rd_length,grp_stats[i][0]->length);
      return err;
    }
    if(grp_stats[i][1]->length!=exp_rd_length){
      sprintf(err,"Read group %d read_2 length expected %"PRIu32" but got %"PRIu32".\n",(i+1),exp_rd_length,grp_stats[i][1]->length);
      return err;
    }
  }

  /******RG1 checks*****/
  //Read count 1
  if(grp_stats[0][0]->count!=exp_rg1_rd1_tot_count){
    sprintf(err,"RG 1, read_1 count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd1_tot_count,grp_stats[0][0]->count);
    return err;
  }
  //Read count 2
  if(grp_stats[0][1]->count!=exp_rg1_rd2_tot_count){
    sprintf(err,"RG 1, read_2 count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd2_tot_count,grp_stats[0][1]->count);
    return err;
  }
  //Duplicate reads 1
  if(grp_stats[0][0]->dups!=exp_rg1_rd1_dups){
    sprintf(err,"RG 1, read_1 duplicate count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd1_dups,grp_stats[0][0]->dups);
    return err;
  }
  //Duplicate reads 2
  if(grp_stats[0][1]->dups!=exp_rg1_rd2_dups){
    sprintf(err,"RG 1, read_2 duplicate count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd2_dups,grp_stats[0][1]->dups);
    return err;
  }
  //GC 1
  if(grp_stats[0][0]->gc!=exp_rg1_rd1_gc){
    sprintf(err,"RG 1, read_1 gc count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd1_gc,grp_stats[0][0]->gc);
    return err;
  }
  //GC 2
  if(grp_stats[0][1]->gc!=exp_rg1_rd2_gc){
    sprintf(err,"RG 1, read_2 gc count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd2_gc,grp_stats[0][1]->gc);
    return err;
  }
  //Unmapped 1
  if(grp_stats[0][0]->umap!=exp_rg1_rd1_umap){
    sprintf(err,"RG 1, read_1 umap count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd1_umap,grp_stats[0][0]->umap);
    return err;
  }
  //Unmapped 2
  if(grp_stats[0][1]->umap!=exp_rg1_rd2_umap){
    sprintf(err,"RG 1, read_2 umap count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd2_umap,grp_stats[0][1]->umap);
    return err;
  }
  //divergent 1
  if(grp_stats[0][0]->divergent!=exp_rg1_rd1_divergent){
    sprintf(err,"RG 1, read_1 divergent count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd1_divergent,grp_stats[0][0]->divergent);
    return err;
  }
  //divergent 2
  if(grp_stats[0][1]->divergent!=exp_rg1_rd2_divergent){
    sprintf(err,"RG 1, read_2 divergent count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd2_divergent,grp_stats[0][1]->divergent);
    return err;
  }
  //mapped_bases 1
  if(grp_stats[0][0]->mapped_bases!=exp_rg1_rd1_mapped_bases){
    sprintf(err,"RG 1, read_1 mapped_bases count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd1_mapped_bases,grp_stats[0][0]->mapped_bases);
    return err;
  }
  //mapped_bases 2
  if(grp_stats[0][1]->mapped_bases!=exp_rg1_rd2_mapped_bases){
    sprintf(err,"RG 1, read_2 mapped_bases count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd2_mapped_bases,grp_stats[0][1]->mapped_bases);
    return err;
  }
  //proper 1
  if(grp_stats[0][0]->proper!=exp_rg1_rd1_proper){
    sprintf(err,"RG 1, read_1 proper count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd1_proper,grp_stats[0][0]->proper);
    return err;
  }
  //proper 2
  if(grp_stats[0][1]->proper!=exp_rg1_rd2_proper){
    sprintf(err,"RG 1, read_2 proper count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg1_rd2_proper,grp_stats[0][1]->proper);
    return err;
  }




  /******RG2 checks*****/
  if(grp_stats[1][0]->count!=exp_rg2_rd1_tot_count){
    sprintf(err,"RG 2, read_1 count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd1_tot_count,grp_stats[1][0]->count);
    return err;
  }
  if(grp_stats[1][1]->count!=exp_rg2_rd2_tot_count){
    sprintf(err,"RG 2, read_2 count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd2_tot_count,grp_stats[1][1]->count);
    return err;
  }
  //Duplicate reads 1
  if(grp_stats[1][0]->dups!=exp_rg2_rd1_dups){
    sprintf(err,"RG 1, read_1 duplicate count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd1_dups,grp_stats[1][0]->dups);
    return err;
  }
  //Duplicate reads 2
  if(grp_stats[1][1]->dups!=exp_rg2_rd2_dups){
    sprintf(err,"RG 1, read_2 duplicate count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd2_dups,grp_stats[1][1]->dups);
    return err;
  }
  //GC 1
  if(grp_stats[1][0]->gc!=exp_rg2_rd1_gc){
    sprintf(err,"RG 1, read_1 gc count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd1_gc,grp_stats[1][0]->gc);
    return err;
  }
  //GC 2
  if(grp_stats[1][1]->gc!=exp_rg2_rd2_gc){
    sprintf(err,"RG 1, read_2 gc count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd2_gc,grp_stats[1][1]->gc);
    return err;
  }
  //Unmapped 1
  if(grp_stats[1][0]->umap!=exp_rg2_rd1_umap){
    sprintf(err,"RG 2, read_1 umap count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd1_umap,grp_stats[1][0]->umap);
    return err;
  }
  //Unmapped 2
  if(grp_stats[1][1]->umap!=exp_rg2_rd2_umap){
    sprintf(err,"RG 2, read_2 umap count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd2_umap,grp_stats[1][1]->umap);
    return err;
  }
  //divergent 1
  if(grp_stats[1][0]->divergent!=exp_rg2_rd1_divergent){
    sprintf(err,"RG 2, read_1 divergent count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd1_divergent,grp_stats[1][0]->divergent);
    return err;
  }
  //divergent 2
  if(grp_stats[1][1]->divergent!=exp_rg2_rd2_divergent){
    sprintf(err,"RG 2, read_2 divergent count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd2_divergent,grp_stats[1][1]->divergent);
    return err;
  }
  //mapped_bases 1
  if(grp_stats[1][0]->mapped_bases!=exp_rg2_rd1_mapped_bases){
    sprintf(err,"RG 2, read_1 mapped_bases count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd1_mapped_bases,grp_stats[1][0]->mapped_bases);
    return err;
  }
  //mapped_bases 2
  if(grp_stats[1][1]->mapped_bases!=exp_rg2_rd2_mapped_bases){
    sprintf(err,"RG 2, read_2 mapped_bases count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd2_mapped_bases,grp_stats[1][1]->mapped_bases);
    return err;
  }
  //proper 1
  if(grp_stats[1][0]->proper!=exp_rg2_rd1_proper){
    sprintf(err,"RG 2, read_1 proper count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd1_proper,grp_stats[1][0]->proper);
    return err;
  }
  //proper 2
  if(grp_stats[1][1]->proper!=exp_rg2_rd2_proper){
    sprintf(err,"RG 2, read_2 proper count incorrect. Expected %"PRIu64" but got %"PRIu64"\n",exp_rg2_rd2_proper,grp_stats[1][1]->proper);
    return err;
  }

  //Run it again but with RNA set ()


	return NULL;
}


char *all_tests() {
   mu_suite_start();
   mu_run_test(test_bam_access_parse_header);
   mu_run_test(test_bam_access_get_mapped_base_count_from_cigar);
   mu_run_test(test_bam_access_process_reads_no_rna);
   mu_run_test(test_bam_access_process_reads_rna);
   return NULL;
}

RUN_TESTS(all_tests);
