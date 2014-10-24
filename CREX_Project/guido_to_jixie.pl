#!/usr/bin/perl

use warnings;
use diagnostics;
use strict;

open INFILE, "<prex_septumfield.dat";
open OUTFILE,">prex_septumfield.jixie.dat";
for my $i (0..7){
    my $input = <INFILE>;
    print OUTFILE $input;
}

my @matrix;

for my $i (0..14){
    for my $j (0..40){
	for my $k (0..40){
	    $matrix[$i][$j][$k] = <INFILE>;
	}
    }
}

close INFILE;

for my $i (0..40){
    for my $k (0..14){
	for my $j (0..40){
	    if( $matrix[$k][$j][$i] =~ /(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)/ ){
		if($2 >= 0){
		    print OUTFILE "$1 $3 $2 $4 $6 $5\n";
		}
	    }
	}
    }
}

close OUTFILE;
