# scan_tcga_files_for_beta.awk:
#
# input file: tab separated file of
# chromosome  start	    stop gene_name
#
# and
# scan through tcga files for matching positions.  take mean of
# matching Beta values and write out.  Also generate a header to
# identify each column.
#
# use barcode file to generate actual file names so that
# primary & metastatic scans can be distinguished
#
# run with:
#
# awk -f scan_tcga_files_for_beta.awk <chr_strt_stop_gene_file>
#
# work with different file list file by adding:
#
# tcgafilelst="m_files.txt"
#
# after '-f scan_3lst_tcga_files.awk'
#
# use different margin <n> with:
#
# margin=<n>
#
# Control column headers with
#  fullheader=1 to get the whole thing
#
# else only the unique patient ID is given
#
# Peter Stockwell, 29-Oct-2015
#
BEGIN{tcga_dir= \
"DNA_methylation/JHU_USC__HumanMethylation450/Level_3/";
tcga_tail=".txt";
margin = 500;
headed=0;
barcode_file="";
unseen=0;
skip_lines=1;
genename_field=0;
fullheader = 0;
}
NR>skip_lines{
if (barcode_file=="")
  {
  printf("No barcode file: need barcode_file=\"<barcodefile>\" in command\n");
  exit(1);
  }
if (!headed)
  {
  printf("#chr\tstart\tend");
  if (genename_field > 0)
    printf("\tGene");
  while (getline fline < barcode_file > 0)
    {
    fname = fileforbarcode(fline);
    if (fname != "")
      {
# create a header for this in the output
      if (fullheader)
        printf("\t%s",fline);
      else
        {
        ns = split(fline,fs,"-");
        printf("\t%s",fs[ns-4]);
        }
      }
    }
  printf("\n");
  headed = 1;
  close(barcode_file);
  }
printf("%s\t%s\t%s",$1,$2,$3);
if (genename_field > 0)
  printf("\t%s",$genename_field);
while (getline fline < barcode_file > 0)
  {
  fname = fileforbarcode(fline);
  if (fname != "")
    {
# found the file, so OK
    nmatch = 0;
    tot = 0.0;
    cmd = sprintf("awk '$4==\"%s\"&&%d<=$5&&$5<=%d&&$2!=\"NA\"' %s",$1,$2-margin,$3+margin,
                    fname);
#print cmd;
    while (cmd | getline ret > 0)
      {
#print ret;
      nrs = split(ret,rs);
      tot += rs[2] + 0.0;
      nmatch++;
      }
    if (nmatch > 0)
      printf("\t%.4f",tot/nmatch);
    else
      printf("\t-");
    close(cmd);
    }
  }
close(barcode_file)
printf("\n");
}
END{
# list unseen files if unseen=1 in cmd line
if(unseen==1)
  for (indx in notseen)
    printf("#Notfound: %s\n",indx);
}  

function fileforbarcode(bcode){
# look in appropriate subdir for a file matching bcode.
# return first matching entry.
# note unseen barcodes if unseen == 1.
# return "" for not found

lscmd=sprintf("ls %s/*%s%s 2> /dev/null",tcga_dir,bcode,tcga_tail);
if (lscmd | getline retval <= 0)
  {
  retval = "";
  if (unseen==1)
    notseen[bcode] = 1;
  }
close(lscmd);
return(retval);
}
