#include <getopt.h>
#include <inttypes.h>
#include "htslib/sam.h"
#include "khash.h"
#include "dbg.h"

KHASH_MAP_INIT_INT(posn,int32_t)
KHASH_MAP_INIT_INT(chrom,khash_t(posn))

char *bam_a_loc = NULL;
char *bam_b_loc = NULL;
char *ref_file = NULL;
int skip_z = 0;
int count_flag_diff = 0;


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

	printf ("Usage: diff_bams -a bam_a.bam -b bam_b.bsm [-r reference.fa] [] [] [-h] [-v]\n\n");
	printf ("Required:\n");
	printf ("-a --bam_a          The first BAM|CRAM file.\n");
	printf ("-b --bam_b          The second BAM|CRAM file.\n\n");
	printf ("Other:\n");
  printf ("-r --ref            Required for CRAM, genome.fa with co-located fai.\n");
  printf ("-c --count          Count flag differences.\n");
  printf ("-s --skip           Don't include reads with MAPQ=0 in comparison.\n\n");
  printf ("-h --help           Display this usage information.\n");
	printf ("-v --version        Prints the version number.\n\n");
  exit(exit_code);
}

void options(int argc, char *argv[]){

  const struct option long_opts[] =
    {
              {"version", no_argument, 0, 'v'},
              {"help",no_argument,0,'h'},
              {"bam_a",required_argument,0,'a'},
              {"ref",required_argument,0,'r'},
              {"bam_b",required_argument,0,'b'},
              {"skip",no_argument,0,'s'},
              {"count",no_argument,0,'c'},
              { NULL, 0, NULL, 0}

   }; //End of declaring opts

   int index = 0;
   int iarg = 0;

     //Iterate through options
   while((iarg = getopt_long(argc, argv, "a:b:r:scvh", long_opts, &index)) != -1){
    switch(iarg){
      case 's':
        skip_z = 1;
        break;

      case 'c':
        count_flag_diff = 1;
        break;

      case 'r':
        ref_file = optarg;
        break;

      case 'a':
        bam_a_loc = optarg;
        break;

      case 'b':
        bam_b_loc = optarg;
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
  if(ref_file != NULL){
    if(check_exist(ref_file) != 1){
      fprintf(stderr,"Reference fasta file (-r) %s does not exist.\n",ref_file);
      print_usage(1);
    }
  }

  if(strcmp(bam_a_loc,bam_b_loc)==0){
    fprintf(stderr,"bam_a '%s' and bam_b '%s' cannot be the same file\n",bam_a_loc,bam_b_loc);
    print_usage(1);
  }

  if(check_exist(bam_a_loc) != 1){
    fprintf(stderr,"First BAM|CRAM file (-a) %s does not exist.\n",bam_a_loc);
    print_usage(1);
  }

  if(check_exist(bam_b_loc) != 1){
    fprintf(stderr,"Second BAM|CRAM file (-b) %s does not exist.\n",bam_b_loc);
    print_usage(1);
  }

  if((strcmp(bam_a_loc + strlen(bam_a_loc) - 4,".cram")==0 || strcmp(bam_b_loc + strlen(bam_b_loc) - 4,".cram")==0) && ref_file == NULL) {
    fprintf(stderr,"Option '-r' must be defined if any CRAM files are provided\n");
    print_usage(1);
  }

   return;
}

int main(int argc, char *argv[]){
  htsFile *htsa = NULL;
  htsFile *htsb = NULL;
  bam_hdr_t *heada = NULL;
  bam_hdr_t *headb = NULL;
  khash_t(chrom) *chr_hash = NULL;
  bam1_t *reada = NULL;
  bam1_t *readb = NULL;
  options(argc, argv);
  //Open bam file a
  htsa = hts_open(bam_a_loc,"r");
  check(htsa != NULL, "Error opening hts file 'a' for reading '%s'.",bam_a_loc);
  //Open bam file b
  htsb = hts_open(bam_b_loc,"r");
  check(htsb != NULL, "Error opening hts file 'b' for reading '%s'.",bam_b_loc);

  if(ref_file){
    hts_set_fai_filename(htsa, ref_file);
    hts_set_fai_filename(htsb, ref_file);
  }

  heada = sam_hdr_read(htsa);
  check(heada != NULL, "Error reading header from opened hts file 'a' '%s'.",bam_a_loc);
  headb = sam_hdr_read(htsb);
  check(headb != NULL, "Error reading header from opened hts file 'b' '%s'.",bam_b_loc);

  if(heada->n_targets != headb->n_targets ){
    sentinel("Reference sequence count is different\n");
  }
  fprintf(stdout,"Reference sequence count passed\n");

  int i=0;
  for(i=0;i<heada->n_targets;i++){
    if(strcmp(heada->target_name[i],headb->target_name[i])!=0){
      sentinel("Reference sequences in different order\n");
    }
  }
  fprintf(stdout,"Reference sequence order passed\n");

  uint64_t count = 0;
  uint64_t flag_diffs = 0;
  uint64_t last_coord = 0;
  int chka = 0;
  int chkb = 0;
  reada = bam_init1();
  readb = bam_init1();
  while(1){
    count++;
    //Check the individual reads
    chka = sam_read1(htsa,heada,reada);
    chkb = sam_read1(htsb,headb,readb);

    if(skip_z==1){
      //Move to next non 0 MQ read
      while(reada->core.qual == 0 && chka >= 0) {
        chka = sam_read1(htsa,heada,reada);
      }
      while(readb->core.qual == 0 && chkb >= 0) {
        chkb = sam_read1(htsb,headb,readb);
      }
    }

    if(chka<0 && chkb<0){
      break;
    }
    if((chka>=0 && chkb <0) || (chkb>=0 && chka<0)){
      sentinel("Files have different number of records\n");
    }

    if(reada->core.tid != readb->core.tid || reada->core.pos != readb->core.pos || strcmp(bam_get_qname(reada),bam_get_qname(readb))!=0){
      sentinel("Files differ at record %"PRIu64" (qname) a=%s b=%s\n",count,bam_get_qname(reada),bam_get_qname(readb));
    }

    if(reada->core.flag != readb->core.flag){
      if(count_flag_diff==1){
        if(chr_hash == NULL){ chr_hash = kh_init(chrom);}
        flag_diffs++;
        int32_t pos = reada->core.pos;
        if((pos - last_coord) >=0 ){ //Checking for different chr
          int res;
          khiter_t k = kh_get(chrom, chr_hash, reada->core.tid); // query the hash table
          khash_t(posn) *pos_hash = NULL;
          if(k == kh_end(chr_hash)){
            k = kh_put(chrom,chr_hash,reada->core.tid,&res);
            pos_hash = kh_init(posn);
            khiter_t kpos = kh_put(posn,pos_hash,pos,&res);
            kh_value(pos_hash,kpos) = 1;
            kh_value(chr_hash,k) = *pos_hash;
          }else{
            *pos_hash = kh_value(chr_hash,k);
            khiter_t kpos = kh_get(posn, pos_hash, pos); // query the pos hash table
            if(kpos == kh_end(pos_hash)){
              kpos = kh_put(posn,pos_hash,pos,&res);
              kh_value(pos_hash,kpos) = 1;
            }else{
              kh_value(pos_hash,kpos) = kh_value(pos_hash,kpos)+1;
            }
            kh_value(chr_hash,k) = *pos_hash;
          }
        }
        last_coord = pos;
      }else{
        sentinel("Files differ at record %"PRIu64" (flags) a=%s b=%s\n",count,bam_get_qname(reada),bam_get_qname(readb));
      }
    }//End of if flags don't match
    if(count % 5000000 == 0) {
      fprintf(stdout,"Matching records: %"PRIu64"",count);
      if(count_flag_diff){
        fprintf(stdout,"\t(flag mismatch: %"PRIu64")",flag_diffs);
      }
      fprintf(stdout,"\r");
    }
  }//End of looping through all reads

  fprintf(stdout,"Matching records: %"PRIu64"\n",count);
  if(count_flag_diff && chr_hash != NULL){
    fprintf(stdout,"Flag mismatches: %"PRIu64"\n",flag_diffs);
    fprintf(stdout,"Locations of flag differences:\n");
    fprintf(stdout,"#Chr\tPos\tCount\n");
    khiter_t k;
    for (k = kh_begin(chr_hash); k != kh_end(chr_hash); ++k) {
      if (kh_exist(chr_hash, k)) {
        const int32_t keyp = kh_key(chr_hash,k);
        khash_t(posn) pos_hash = kh_value(chr_hash, k);
        khiter_t kpos;
        for (kpos = kh_begin(&pos_hash); kpos != kh_end(&pos_hash); ++kpos) {
          if (kh_exist(&pos_hash, kpos)) {
            const int32_t key = kh_key(&pos_hash,kpos);
            int32_t counter = kh_value(&pos_hash,kpos);
            fprintf(stdout,"%s\t%"PRIi32"\t%"PRIi32"\n",heada->target_name[keyp], key, counter);
            kh_destroy(posn,&pos_hash);
          }
        }

      }
    }
    kh_destroy(chrom,chr_hash);
  }

  bam_destroy1(reada);
  bam_destroy1(readb);
  bam_hdr_destroy(heada);
  bam_hdr_destroy(headb);
  hts_close(htsa);
  hts_close(htsb);
  return 0;
error:
  if(count_flag_diff && chr_hash != NULL){
    khiter_t k;
    for (k = kh_begin(chr_hash); k != kh_end(chr_hash); ++k) {
      if (kh_exist(chr_hash, k)) {
        const int32_t keyp = kh_key(chr_hash,k);
        khash_t(posn) pos_hash = kh_value(chr_hash, k);
        kh_destroy(posn,&pos_hash);
      }
    }
    kh_destroy(chrom,chr_hash);
  }
  if(reada) bam_destroy1(reada);
  if(readb) bam_destroy1(readb);
  if(heada) bam_hdr_destroy(heada);
  if(headb) bam_hdr_destroy(headb);
  if(htsa) hts_close(htsa);
  if(htsb) hts_close(htsb);
  return 1;
}