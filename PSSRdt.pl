#!/usr/bin/perl -w

# Author: Ruizheng Tian
# Program name: PSSRdt.pl-Polymorphic SSRs digging tool

###_______________________________________________________________________________
###
### Program name: PSSRdt.pl
### Author:       Ruizheng Tian
### Release date: 2019-9-1(version 1.0)
### Organization: Northwest A&F University

###
###
## _______________________________________________________________________________
## _______________________________________________________________________________
##
##
## DESCRIPTION: Program for the screening and research of 
##  			the information of polymorphic SSR/STR loci
##
## SYNTAX:   PSSRdt.pl Parameter <FASTA file>
##
##  	  <FASTAfile>    Single file in FASTA format ,and the fasta file was a data collection
##		  from the multiple transcriptome data of same specie .
##
##
## EXAMPLE: EXAMPLE: PSSRdt.pl  10  Species1.fasta

##
## _______________________________________________________________________________
##


#^^^ DECLARATION ^^^#

# Please check for arguments and input file ,if no destination file. #


use strict;
use Bio::SeqIO;
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^##
#use Bio::Tools::TandemRepeatsFinder;
#Bio::perl needs to be pre-installed.
#The usual installation process is as follows:
# Ubuntu:
# 	1.cpanm installation
#		wget http://xrl.us/cpanm -O /usr/bin/cpanm;
#		chmod +x /usr/bin/cpanm
# 	2.Bio::perl installation
#		cpanm --mirror http://mirrors.163.com/cpan --mirror-only Bio::Perl
# Windows(Installing strawberry perl first):
# 	1. Open cmd, input "cpan";
# 	2. input  install one bioperl version, sush as "install CJFIELDS/BioPerl-1.007001.tar.gz";
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^##


# Canceling the following "#" of three lines  , you can get the run process record in log.
#open(OUT,">findStr.log");
#select OUT;
# print("id,start_index,max_end_index,start_bp,repeated_element,repeated_times,end_bp\n");

my %STRS;
my $para = $ARGV[0];
my $seqio_obj = Bio::SeqIO->new(-file => $ARGV[1], -format => "fasta" );
while (my $seq_obj = $seqio_obj->next_seq ) {
	my @agcts=split('',$seq_obj->seq); #Converting records to character arrays.
	my $agct_length= length $seq_obj->seq;
	my $start_index=0;
	for(my $start_index=0;$start_index<$agct_length-1;){
		#Finding the STR from the leftmost element of the sequence array until the sequence array end.
		my @end_indexs=(0,0,0,0,0,0);
		#Recording the termination of the subscript shown in the search process, a total of 6 cycles.
		#So an array of 6 in length is required to record the termination subscripts when the periods are 1-6.
		my @repeated_times=(0,0,0,0,0,0); # Recording the number of repetition cycles when the periods are 1, 3, 4, 5, and 6 respectively
		for(my $t=1;$t<=6;$t=$t+1){ # Searching the length and the number of cycles of STRs when the number of cycles is 1, 3, 4, 5, and 6 respectively
			my $breaks=0;
			for(my $j=$start_index;$breaks==0&&$j<$agct_length-1;$j=$j+$t){
				my $k;
				for($k=0;$k<$t;$k++){
				# Comparing all the elements in a cycle. Only when a cycle of all the elements are the same, repeating completely a cycle.		
					if($j+$k+$t>=$agct_length-1 || ($agcts[$j+$k]  ne $agcts[$j+$k+$t])){
						$breaks=1;
						last;
					}
				}
				if($breaks==0){
					$repeated_times[$t-1]=$repeated_times[$t-1]+1;
					$end_indexs[$t-1]=$j+$k+$t;
				}
			}
		}
		my $max_end_index=$start_index;
		my $max_repeat_peroid=0;
		#Finding the longest STR during this search
		for(my $i=1;$i<=6;$i=$i+1){
			if($end_indexs[$i-1]>$max_end_index){
				$max_end_index=$end_indexs[$i-1];
				$max_repeat_peroid=$i;
			}
		}
		if((($max_end_index-$start_index>9 and $repeated_times[$max_repeat_peroid-1] >5 and $max_repeat_peroid!=2) or ($max_repeat_peroid==2 and $repeated_times[$max_repeat_peroid-1] >5))){
			#If it finds the STRs , the STRs will be put into the final result.
			my $repeated_times = $repeated_times[$max_repeat_peroid-1]+1;
			my $start_bp = "";
			if($start_index!=0){
				$start_bp =substr($seq_obj->seq,$start_index-$para>0?$start_index-$para:0,$start_index>$para?$para:$start_index);
			}
			my $repeated_element = substr($seq_obj->seq,$start_index,$max_repeat_peroid);
			my $end_bp = substr($seq_obj->seq,$max_end_index,$para);
			my $str = "$start_bp\_$repeated_element\_$end_bp";
			if($start_index>$para and $max_end_index+$para <$agct_length){			
				if( exists $STRS{$str}) {
					push @{$STRS{$str}}, $repeated_times[$max_repeat_peroid-1]+1;
				} else {
					my @nums=($repeated_times[$max_repeat_peroid-1]+1);
					$STRS{$str}=\@nums;
				}
			}
	#########		printf("%s,\t%d,\t%d,\t%s,\t%s,\t%d,\t%s\n",$seq_obj->id,$start_index,$max_end_index,$start_bp,$repeated_element,$repeated_times,$end_bp);
			$start_index=$max_end_index-$max_repeat_peroid>$start_index?$max_end_index-$max_repeat_peroid-1:$max_end_index-1;	
		}else{
			$start_index=$start_index+1;
		}
	}
}

open(OUT,">SSRs.total.details");
select OUT;
print ("Loci\tTotal_Sequence\tRepeat_times\n");
foreach my $key (keys %STRS){
	my $count = @{$STRS{$key}};
	print ($key,"	",$count ,"	");
	foreach my $num (@{$STRS{$key}}){
		print ($num,";");
	}
	print ("\n");
}

open(OUT,">SSRs.details");
select OUT;
#Judging the result file and only leaving the third column.
#Only reserving the result lines that existed differences among the digital in third column.
#for instance:  xxxxxxxx_ac_xxxxxxxx,  3   10 ,10, 10,  Three numbers "10" are the same, deleting this line.
#		   but  xxxxxxxx_ac_xxxxxxxx, 3   10,12,10   Keep this result line.
print ("Loci\tTotal_Sequence\tRepeat_times\n");
foreach my $key (keys %STRS){
	my $count = @{$STRS{$key}};
	my $should_print=0;
	for(my $i=1;$i<$count;$i++){
		if($STRS{$key}[$i]!=$STRS{$key}[0]){
			$should_print=1;
		}
	}
	if($should_print==1){
		print ($key,"	",$count ,"	");
		foreach my $num (@{$STRS{$key}}){
			print ($num,";");
		}
		print ("\n");
	}
}

