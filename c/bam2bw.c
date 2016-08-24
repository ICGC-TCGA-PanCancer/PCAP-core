#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include "bigWig.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "dbg.h"

char *out_file = "output.bam.bw";
char *input_file = NULL;
char *region_store = NULL;
char **our_region_list = NULL;
uint32_t region_list_count = 0;
int is_regions_file = 0;
uint8_t base_bit = 0;
int filter = 0;
char base = 0;
uint32_t single = 1;

typedef struct {
  uint32_t ltid;
  int      lstart,lcoverage,lpos,beg,end;
  htsFile *in;
  hts_idx_t *idx;
	bam_hdr_t *head;
	bigWigFile_t *out;
} tmpstruct_t;

void print_version (int exit_code){
  printf ("%s\n",VERSION);
	exit(exit_code);
}

int line_count (char *file_path){
  FILE *f = fopen(file_path,"r");
  int line_count = 0;
  check(f != NULL, "Error opening file '%s' to count lines.",file_path);
  char rd[512];
	while(fgets(rd, 512, f) != NULL){
    line_count++;
  }
  fclose(f);
  return line_count;
error:
  if(f) fclose(f);
  return -1;
}

int check_exist(char *fname){
	FILE *fp;
	if((fp = fopen(fname,"r"))){
		fclose(fp);
		return 1;
	}
	return 0;
}

void print_usage (int exit_code){
	printf("Usage: bam2bw -i input.[b|cr]am -o output.bw\n\n");
	printf("-i  --input [file]                                Path to the input [b|cr]am file.\n");
	printf("-F  --filter [int]                                SAM flags to filter. [default: %d]\n",filter);
	printf("-o  --outfile [file]                              Path to the output .bw file produced. [default:'%s']\n\n",out_file);
	printf("Optional: \n");
	printf("-b  --base [char]                                 To produce a single base coverage bigwig file [ACTG]. Must be used in conjunction with -r|--region\n");
	printf("-r  --region [file]                               A samtools style region (contig:start-stop) or a bed file of regions over which to produce the bigwig file\n\n");
  printf ("Other:\n");
	printf("-h  --help                                        Display this usage information.\n");
	printf("-v  --version                                     Prints the version number.\n\n");

  exit(exit_code);
}

int check_region_string(char *region){
  if(check_exist(region)==1){
    is_regions_file = 1;
    return 1;
  }else{
    //Parse this as a region
    int beg = 0;
    int end = 0;
    const char *res = hts_parse_reg(region, &beg, &end);
    if(res) return 1;
  }
  return -1;
}

int get_int_length(int input){
  return (input == 0 ? 1 : (int)(log10(input)+1));
}

int get_no_of_SQ_lines(char *bam_loc){
  htsFile *bam =  NULL;
  char *line = NULL;
  bam_hdr_t *header = NULL;
	bam = hts_open(bam_loc, "r");
	check(bam != 0,"Bam file %s failed to open to read header.",bam_loc);
	header = sam_hdr_read(bam);
	char *head_txt = header->text;
  line = strtok(head_txt,"\n");
  int count = 0;
	while(line != NULL){
		//Check for a read group line
		if(strncmp(line,"@SQ",3)==0){
		  count++;
		}//End of if this is an SQ line
		line = strtok(NULL,"\n");
	}
	free(line);
	hts_close(bam);
	bam_hdr_destroy(header);
  return count;
error:
  if(line) free(line);
	if(bam) hts_close(bam);
	if(header) bam_hdr_destroy(header);
  return -1;
}

int parse_SQ_line(char *line, char *name, int *length){
  char *tag = NULL;
  char *tmp = NULL;
  char *ptr = NULL;
  int nom = 0;
  int len = 0;
  tag = strtok_r(line,"\t",&ptr);
  while(tag != NULL){
    int chk=0;
    tmp = malloc(sizeof(char) * 512);
    check_mem(tmp);
    chk = sscanf(tag,"SN:%[^\t\n]",tmp);
    if(chk>0){
      strcpy(name,tmp);
      tag = strtok_r(NULL,"\t",&ptr);
      nom = 1;
      continue;
    }
    chk = sscanf(tag,"LN:%d",length);
    if(chk>0){
      len = 1;
      tag = strtok_r(NULL,"\t",&ptr);
      continue;
    }
    tag = strtok_r(NULL,"\t",&ptr);
  }
  if(tmp) free(tmp);
  if(tag) free(tag);

  if(nom && len) {
    return 1;
  }else{
    return -1;
  }
error:
  if(tmp) free(tmp);
  if(tag) free(tag);
  return -1;
}

int build_chromList_from_bam(chromList_t *chromList ,char *bam_loc){
  htsFile *bam =  NULL;
  char *line = NULL;
  char *tmp = NULL;
  bam_hdr_t *header = NULL;
  char *ptr;
	bam = hts_open(bam_loc, "r");
	check(bam != 0,"Bam file %s failed to open to read header.",bam_loc);
	header = sam_hdr_read(bam);
	char *head_txt = header->text;
  line = strtok(head_txt,"\n");
  int i = 0;
	while(line != NULL){
		//Check for a read group line
		if(strncmp(line,"@SQ",3)==0){
      char *tag = strtok_r(line,"\t",&ptr);
      while(tag != NULL){
        int chk=0;
        tmp = malloc(sizeof(char) * 512);
        check_mem(tmp);
        chk = sscanf(tag,"SN:%[^\t\n]",tmp);
        if(chk>0){
          chromList->chrom[i] = malloc(sizeof(char) * (strlen(tmp)+1));
          check_mem(chromList->chrom[i]);
          strcpy(chromList->chrom[i],tmp);
          tag = strtok_r(NULL,"\t",&ptr);
          free(tmp);
          continue;
        }
        free(tmp);
        uint32_t tmpint = 0;
        chk = sscanf(tag,"LN:%" SCNu32 "",&tmpint);
        if(chk>0){
          chromList->len[i] = tmpint;
          tag = strtok_r(NULL,"\t",&ptr);
          continue;
        }
        tag = strtok_r(NULL,"\t",&ptr);
      }//End of SQ line tags
      i++;
		}//End of if this is an SQ line
		line = strtok(NULL,"\n");
	}
	free(line);
	bam_hdr_destroy(header);
	hts_close(bam);
  return 1;
error:
  if(line) free(line);
	if(bam) hts_close(bam);
	if(header) bam_hdr_destroy(header);
	if(tmp) free (tmp);
  return -1;
}

void setup_options(int argc, char *argv[]){
	const struct option long_opts[] =
	{
             	{"input", required_argument, 0, 'i'},
             	{"filter", required_argument, 0, 'F'},
             	{"outfile",required_argument, 0, 'o'},
             	{"base",required_argument, 0, 'b'},
             	{"region",required_argument, 0, 'r'},
             	{"help", no_argument, 0, 'h'},
             	{"version", no_argument, 0, 'v'},
             	{ NULL, 0, NULL, 0}

   }; //End of declaring opts

   int index = 0;
   int iarg = 0;

   //Iterate through options
   while((iarg = getopt_long(argc, argv, "F:i:o:b:r:hv",long_opts, &index)) != -1){
    switch(iarg){
      case 'F':
        if(sscanf(optarg, "%i", &filter) != 1){
      		fprintf(stderr,"Error parsing -F|--filter argument '%s'. Should be an integer > 0",optarg);
      		print_usage(1);
      	}
        break;
   		case 'i':
				input_file = optarg;
				if(check_exist(input_file) != 1){
          fprintf(stderr,"Input bam file %s does not appear to exist.\n",input_file);
          print_usage(1);
        }
   			break;
   		case 'o':
				out_file = optarg;
   			break;
			case 'h':
				print_usage (0);
				break;
			case 'v':
				print_version (0);
				break;
			case 'r':
			  region_store = optarg;
			  //First check for a region format
			  int res = check_region_string(region_store);
        if(res<0){
          fprintf(stderr,"Region %s is not in correct format or not an existing bed file.\n",region_store);
          print_usage(1);
        }
			  break;
			case 'b':
			  if(sscanf(optarg, "%c", &base) != 1){
      		fprintf(stderr,"Error parsing -b|--base argument '%s'. Should be a single character > 0",optarg);
      		print_usage(1);
      	}
			  if(base != 'A' && base != 'C' && base != 'G' && base != 'T'){
          fprintf(stderr,"Single base %c provided must match [ACGT].\n",base);
          print_usage(1);
			  }else{
          switch (base){
            case 'A':
              base_bit = 1;
              break;
            case 'C':
              base_bit = 2;
              break;
            case 'G':
              base_bit = 4;
              break;
            case 'T':
              base_bit = 8;
              break;
          }
			  }
			  break;
			case '?':
        print_usage (1);
        break;
      default:
      	print_usage (0);
   	}; // End of args switch statement

   }//End of iteration through options

  if(input_file==NULL){
    fprintf(stderr,"Required option -i|--input-bam not defined.\n");
    print_usage(1);
  }

  if(base_bit != 0 && region_store==NULL){
    fprintf(stderr,"Option -r|--region must be used with the -b|--base option.\n");
    print_usage(1);
  }

  return;
}

// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t position, int n, const bam_pileup1_t *pl, void *data)
{
  tmpstruct_t *tmp = (tmpstruct_t*)data;
  int pos          = (int)position;
  int coverage     = n;
  int i;
  for (i=0;i<n;i++)
    if (pl[i].is_del ||
          (base_bit>0 && bam_seqi(bam_get_seq(pl[i].b), pl[i].qpos)!=base_bit)) coverage--;
  if (tmp->ltid != tid || tmp->lcoverage != coverage || pos > tmp->lpos+1) {
    if (tmp->lpos > 0 && tmp->lcoverage > 0){
      uint32_t start =  tmp->lstart;
      uint32_t stop = tmp->lpos+1;
      float cvg = tmp->lcoverage;
      bwAddIntervals(tmp->out,&tmp->head->target_name[tmp->ltid],&start,&stop,&cvg,single);
    }
    tmp->ltid       = tid;
    tmp->lstart     = pos;
    tmp->lcoverage  = coverage;
  }
  tmp->lpos = pos;
  return 0;
}

int main(int argc, char *argv[]){
	setup_options(argc, argv);
	tmpstruct_t tmp;
	int no_of_regions = 0;
	int sq_lines = get_no_of_SQ_lines(input_file);
	FILE *fp_bed = NULL;
	if(region_store){
    if(is_regions_file==1){
      //If we have a bedfile
      no_of_regions = line_count(region_store);
      check(no_of_regions>0,"Error counting entries in region file %s\n",region_store);
      our_region_list = calloc(no_of_regions, sizeof(char *));
      check_mem(region_list);
      int i=0;
      fp_bed = fopen(region_store,"r");
        char line[512];
        while(fgets(line,512,fp_bed)){
          char *contig = malloc(sizeof(char)*256);
          check_mem(contig);
          int beg,end;
          int chk = sscanf(line,"%s\t%d\t%d",contig,&beg,&end);
          check(chk==3,"Error reading line '%s' from regions bed file.",line);
          char *region = malloc(sizeof(char) * (strlen(contig)+get_int_length(beg)+get_int_length(end))+3);
          check_mem(region);
          chk = sprintf(region,"%s:%d-%d",contig,beg+1,end); //+1 to beginning as this is bed
          check(chk==((strlen(contig)+get_int_length(beg)+get_int_length(end))+2),"Error creating region line from bed entry '%s'.",line);
          our_region_list[i] = region;
          i++;
        }
      fclose(fp_bed);
    }else{
      //If we have a single region....
      no_of_regions = 1;
      our_region_list = calloc(no_of_regions, sizeof(char *));
      check_mem(region_list);
      our_region_list[0] = region_store;
    }
	}else{
    //If not...
    //Build a list of regions from the header of the bamfile
    no_of_regions = sq_lines;
    our_region_list = calloc(no_of_regions, sizeof(char *));
    htsFile *bam =  NULL;
    char *line = NULL;
    bam = hts_open(input_file, "r");
      check(bam != 0,"Bam file %s failed to open to read header.",input_file);
      bam_hdr_t *header = sam_hdr_read(bam);
      char *head_txt = header->text;
      line = strtok(head_txt,"\n");
      int i=0;
      while(line != NULL){
        //Check for a read group line
        if(strncmp(line,"@SQ",3)==0){
          char *contig = malloc(sizeof(char) * 256);
          int end,chk;
          chk = parse_SQ_line(line,contig,&end);
          check(chk==1,"Error parsing SQ line %s.\n",line);
          char *region = malloc(sizeof(char) * (strlen(contig)+get_int_length(end))+4);
          check_mem(region);
          chk = sprintf(region,"%s:%d-%d",contig,1,end);
          check(chk==((strlen(contig)+get_int_length(end))+3),"Error creating region line from bed entry '%s'.",line);
          our_region_list[i] = region;
          i++;
        }//End of if this is an SQ line
        line = strtok(NULL,"\n");
      }
      bam_hdr_destroy(header);
    free(line);
    hts_close(bam);
	}

  check(no_of_regions>0,"Error evaluating regions to parse. None found.");

  //Now create a bigwigfile
  //Open file as a bw file
  chromList_t *chromList = NULL;
  bigWigFile_t *bw_out = NULL;
  bam_plp_t buf = NULL;
	bam1_t *b = NULL;
  hts_itr_t *iter = NULL;
  //Generate the list of chromosomes using the bam header.
  chromList = calloc(1, sizeof(chromList_t));
  check_mem(chromList);
  chromList->nKeys = sq_lines;
  chromList->chrom = calloc(sq_lines,sizeof(char *));
  check_mem(chromList->chrom);
  chromList->len = calloc(sq_lines,sizeof(uint32_t));
  check_mem(chromList->len);
  //Iterate through the header of the bam file and get all contigs.
  int res = build_chromList_from_bam(chromList,input_file);
  check(res==1,"Error building chromList from bam header.");

  //Open output file
  bw_out = bwOpen(out_file, NULL, "w");
  check(bw_out!=NULL, "Error opening %s for writing.",out_file);
  //Initialise bw
  res = bwInit(1<<17);
  check(res==0,"Received an error in bwInit");
  //10 zoom levels
  res = bwCreateHdr(bw_out, 10);
  check(res==0,"Received error %d in bwCreateHeader",res);
  //Create the chromosome lists
  bw_out->cl = bwCreateChromList(chromList->chrom, chromList->len, chromList->nKeys);
  check(bw_out->cl != NULL, "Error generating chromlist as bw.");
  res = bwWriteHdr(bw_out);
  check(res==0,"Error %d writing bw header.",res);

  tmp.beg = 0; tmp.end = 0x7fffffff;
  tmp.lstart    = 0;
  tmp.lcoverage = 0;
  tmp.ltid      = 0;
  tmp.lpos      = 0;
  tmp.out       = bw_out;

  //Open input file
  tmp.in = hts_open(input_file, "r");
  check(tmp.in!=0,"Failed to open [CR|B]AM file '%s'.", input_file);
  tmp.head = sam_hdr_read(tmp.in);
  //Now we generate the bw info
  int i=0;
	for(i=0;i<no_of_regions;i++){
	  b = bam_init1();
    buf = bam_plp_init(0,0);

    tmp.idx = sam_index_load(tmp.in,input_file);
    check(tmp.idx!=0,"Fail to open [CR|B]AM index for file '%s'.", input_file);
    iter = sam_itr_querys(tmp.idx, tmp.head, our_region_list[i]);

     //TODO bam to bedgraph but instead put into BW.
    int result;
    while ((result = sam_itr_next(tmp.in, iter, b)) >= 0) {
      if((b->core.flag & filter)>0) continue; //Skip if this is a filtered read
      int ret, n_plp, tid, pos;
      const bam_pileup1_t *plp;
      ret = bam_plp_push(buf, b);
      if (ret < 0) break;
      while ((plp = bam_plp_next(buf, &tid, &pos, &n_plp)) != 0){
          pileup_func(tid, pos, n_plp, plp, &tmp);
      }
    }
    bam_plp_push(buf,0);

    //Check we've written everything...
    int n_plp, tid, pos;
    const bam_pileup1_t *plp;
    while ((plp = bam_plp_next(buf, &tid, &pos, &n_plp)) != 0){
      pileup_func(tid, pos, n_plp, plp, &tmp);
    }
    uint32_t start =  tmp.lstart;
    uint32_t stop = tmp.lpos+1;
    float cvg = tmp.lcoverage;
    bwAddIntervals(tmp.out,&tmp.head->target_name[tmp.ltid],&start,&stop,&cvg,single);

    bam_plp_destroy(buf);
    bam_destroy1(b);
    sam_itr_destroy(iter);
	}

  bwClose(bw_out);
  bwCleanup();
  hts_close(tmp.in);
  int clean=0;
  for(clean=0; clean<sq_lines; clean++){
    if(chromList->chrom){
      free(chromList->chrom[clean]);
    }
  }
  free(chromList->len);
  free(chromList);
  hts_idx_destroy(tmp.idx);
  bam_hdr_destroy(tmp.head);
  free(our_region_list);
  return 0;

error:
  if(our_region_list) free(our_region_list);
  if(fp_bed) fclose(fp_bed);
  return -1;
}
