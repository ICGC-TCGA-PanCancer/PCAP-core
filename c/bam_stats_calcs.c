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

#include <math.h>
#include <stdlib.h>
#include <inttypes.h>
#include "bam_access.h"

int compare( const void* a, const void* b)
{
     uint64_t int_a = * ( (uint64_t*) a );
     uint64_t int_b = * ( (uint64_t*) b );
     if ( int_a == int_b ) return 0;
     else if ( int_a < int_b ) return -1;
     else return 1;
}

int bam_stats_calcs_calculate_mean_sd_median_insert_size(khash_t(ins) *inserts,double *mean, double *sd, double *median){

    uint64_t pp_mean = 0;
    uint64_t tt_mean = 0;
    uint64_t key;
    uint64_t val;

    uint64_t *insert_bins;

    insert_bins = malloc(sizeof(uint64_t) * kh_size(inserts));

    int i=0;
    kh_foreach(inserts,key,val,
        { pp_mean += key * val;
          tt_mean += val;
          insert_bins[i] = key;
          i++;
        });

    if(tt_mean){//Calculate mean , median, sd
      *mean = (double) ((double)pp_mean/(double)tt_mean);


      //sort the array of insert sizes.
      qsort( insert_bins, kh_size(inserts), sizeof(uint64_t), compare );

      int midpoint2 = (int)(tt_mean / 2);
      int midpoint = (midpoint2 + 1);
      uint64_t insert = 0;
      uint64_t prev_insert = 0;
      uint64_t running_total = 0;
      uint64_t current_bin_count = 0;


      int j=0;
      for(j=0; j<=kh_size(inserts); j++){
        insert = insert_bins[j]; // The insert size...
        khint_t k;
        k = kh_get(ins,inserts,insert);
        uint64_t val = kh_val(inserts,k);
        running_total += val;
        current_bin_count = val;
        if(running_total >= midpoint) break;
        prev_insert = insert_bins[j];
      }

      if(tt_mean %2 == 0 && ( running_total - midpoint2 >= current_bin_count )){
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



      if(tt_sd){
        double variance = fabs((double)((double)pp_sd / (double)tt_sd));
        *sd = sqrt(variance);
      }else{
        *sd = 0;
      }

      kh_destroy(ins, inserts);

    } //End of if we have data to calculate from.
  return 0;
}