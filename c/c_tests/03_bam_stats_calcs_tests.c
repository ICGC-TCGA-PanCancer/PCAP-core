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

#include <inttypes.h>
#include "minunit.h"
#include "bam_stats_calcs.h"
#include "khash.h"

KHASH_MAP_INIT_INT(insert,uint64_t)
//KHASH_INIT2(ins,, khint32_t, uint64_t, 1, kh_int_hash_func, kh_int_hash_equal)

double exp_mean = 150;
double exp_sd = 50;
double exp_median = 150;
char err[100];

char *bam_stats_calcs_calculate_mean_sd_median_insert_size_test(){
  int res;
  khash_t(ins) *inserts;
  inserts = kh_init(ins);
  khint_t k;
  k = kh_put(ins,inserts,abs(200),&res);
  kh_value(inserts,k) = 50;
  k = kh_put(ins,inserts,abs(100),&res);
  kh_value(inserts,k) = 50;

  double mean;
  double sd;
  double median;

  int check = bam_stats_calcs_calculate_mean_sd_median_insert_size(inserts, &mean, &sd, &median);
  if(check != 0) {
    sprintf(err,"Calculation failed to complete\n");
    return err;
  }

  if(mean != exp_mean){
    sprintf(err,"Mean from calculation %f is not as expected %f\n",mean,exp_mean);
    return err;
  }
  if(sd != exp_sd){
    sprintf(err,"SD from calculation %f is not as expected %f\n",sd,exp_sd);
    return err;
  }
  if(median != exp_median){
    sprintf(err,"median from calculation %f is not as expected %f\n",median,exp_median);
    return err;
  }



  return NULL;
}

char *all_tests() {
   mu_suite_start();
   mu_run_test(bam_stats_calcs_calculate_mean_sd_median_insert_size_test);
   return NULL;
}

RUN_TESTS(all_tests);
