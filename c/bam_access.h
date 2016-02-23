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

#ifndef __bam_access_h__
#define __bam_access_h__

#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include "htslib/sam.h"
#include "dbg.h"
#include "khash.h"

KHASH_MAP_INIT_INT(ins,uint64_t)
//KHASH_INIT2(ins,, khint32_t, uint64_t, 1, kh_int_hash_func, kh_int_hash_equal)

typedef struct {
  uint32_t length;
	uint64_t count;
	uint64_t dups;
  uint64_t gc;
  uint64_t umap;
  uint64_t divergent;
  uint64_t mapped_bases;
  uint64_t proper;
  //list of counts of possible insert sizes....
  khash_t(ins) *inserts; //from counts of insert size from 0-200000 bp. See how this works, might need amore dynamic means for some data.
  //FQP is not included as we're not covering quality plots yet.
} stats_rd_t;

typedef struct{
  char *id;
  char *platform;
  char *platform_unit;
  char *lib;
  char *sample;
} rg_info_t;

rg_info_t **bam_access_parse_header(bam_hdr_t *head, int *grps_size, stats_rd_t ****grp_stats);

int bam_access_process_reads(htsFile *input, bam_hdr_t *head, rg_info_t **grps, int grps_size, stats_rd_t ****grp_stats, int rna);

uint64_t bam_access_get_mapped_base_count_from_cigar(bam1_t *b);

#endif
