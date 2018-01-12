#! /usr/bin/perl -w
use strict;

###testing github

### read in all the .cntTable files
my %files; #hash mapping treatment filename to control filename
my %geneValues; #hash containing all the gene counts
my %teValues; #hash containing all the TE counts
my %l1Values; #hash containing all the L1 counts
my %h2bValues;
my $lineCount; #variable that tracks the line number in current file
my %index; #hash that indexes treatment filename to a unqiue number
my $count=0; #counter

my %H2Bindex = ("NM_001002916" => 0, "NM_170610" => 1,
"NM_021062" => 2, "NM_003526" => 3,
"NM_021063" => 4, "NM_138720" => 5,
"NM_003523" => 6, "NM_003522" => 7,
"NM_003518" => 8, "NM_003524" => 9,
"NM_003525" => 10, "NM_021058" => 11,
"NM_001312653" => 12, "NM_080593" => 13,
"NM_003519" => 14, "NM_003521" => 15,
"NM_003520" => 16, "NM_003527" => 17,
"NM_003528" => 18);

foreach my $filename (@ARGV){
    open(INPUT,$filename) or die "unable to open $filename\n";
    # if ($filename =~ m"^rmsk") {
   	    $lineCount=0;
   	    my @a;
	    while(my $line = <INPUT>){
	        if ($lineCount==0) {
	        	chomp $line;
	        	@a=split(/\t/,$line);
	        	$a[1]=~s/\.T$//g;
	        	$a[2]=~s/\.C$//g;
	        	$files{$a[1]}=$a[2]; #treatment ---> control
	        	$teValues{$a[1]}=[]; #initializing all the arrays inside the hashes
	        	$geneValues{$a[1]}=[];
	        	$l1Values{$a[1]}=[];
	        	$h2bValues{$a[1]}=[];
	        	$teValues{$a[2]}=[];
	        	$geneValues{$a[2]}=[];
	        	$l1Values{$a[2]}=[];
	        	$h2bValues{$a[2]}=[];
	        	unless (defined $index{$a[1]}) {
	       			$index{$a[1]}=$count;
	       			$count++; 		
	        	}
	        	unless (defined $index{$a[2]}) {
	        		$index{$a[2]}=$count;
	        		$count++;
	        	}
	        }
	        else{
	        	chomp $line;
	            my @b=split(/\t/,$line);
	            if ($b[0] =~ m":") {
	            	push @{$teValues{$a[1]}},$b[1];
	            	push @{$teValues{$a[2]}},$b[2];
	            	$teValues{$a[2]}=[$b[2]];
	            } 
	            else {
	            	push @{$geneValues{$a[1]}},$b[1];
	            	push @{$geneValues{$a[2]}},$b[2];
	            	foreach my $key (keys %H2Bindex){
	            		if($b[0]=~m/$key/){
	            			@{$h2bValues{$a[1]}}[$H2Bindex{$key}]=$b[1];
	            			@{$h2bValues{$a[2]}}[$H2Bindex{$key}]=$b[2];
	            		}
	            	}
	            }
	        }
	        $lineCount++;
	    }
    # } else {
    # 	$lineCount=0;
    # 	my @a;
	   #  while(my $line = <INPUT>){
	   #      if ($lineCount==0) {
	   #      	chomp $line;
	   #      	@a=split(/\t/,$line);
	   #      	$a[1]=~s/\.T$//g;
	   #      	$a[2]=~s/\.C$//g;
	   #      	$files{$a[1]}=$a[2]; #treatment ---> control
	   #      }
	   #      else{
	   #      	chomp $line;
	   #          my @b=split(/\t/,$line);
	   #          if ($b[0] =~ m":") {
	   #          	if (defined $l1Values{$a[1]}) {
	   #          		push @{$l1Values{$a[1]}},$b[1];
	   #          	} else {
	   #          		$l1Values{$a[1]}=[$b[1]];
	   #          	}
	   #          }
	   #      }
	   #      $lineCount++;
	   #  }	
    # }
}

### delete all the lines which have 0 values
# find the lengths of the 3 arrays i.e. gene, te and l1
# loop through 


my %sums;
foreach my $key1 (keys %files) {
	foreach my $key2 (@{$geneValues{$key1}}) {
		$sums{$key1}+=$key2;
	}
	foreach my $key2 (@{$teValues{$key1}}) {
		$sums{$key1}+=$key2;
	}
	foreach my $key2 (@{$geneValues{$files{$key1}}}) {
		$sums{$files{$key1}}+=$key2;
	}
	foreach my $key2 (@{$teValues{$files{$key1}}}) {
		$sums{$files{$key1}}+=$key2;
	}
}

my $mean=0;
foreach my $key (keys %sums) {
	print "$key ---> $sums{$key}\n";
	# foreach my $key1 (@{$h2bValues{$key}}){  # part that prints out the H2B expression
	# 	print "$key1\n";
	# }
	$mean+=$sums{$key};
}
$mean/=scalar(keys %sums);

print "mean = $mean\n";

### problem, how to deal with removing all 0 lines???????????????????????
### normalize based on sum of all reads
### print-out each set of treatment-control to a temp file and call R script for it which will make the plots
### delete all temp files created and done