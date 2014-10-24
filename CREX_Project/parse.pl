#!/usr/bin/perl

open INFILE, "<rtsep102008.map";
for my $i (0..5){
    my $input = <INFILE>;
    #print "$input";
}
my @slurp = <INFILE>;
close INFILE;

my @separ;
for ( @slurp ){
    my @numbers = ( $_ =~ /(\-?\d\.\d+(?:E(?:\+|\-|\W)?\d+)?)/g );
    #print "@numbers\n";
    for ( @numbers ){
	push( @separ, $_ );
    }
}


open OUTFILE, ">parse.dat";
my $index = 1;
for ( @separ ) {
    print OUTFILE "$_ ";
    if ( $index % 3 == 0 && $index != 0 ) {
	print OUTFILE "\n";
    }
    $index++;
}
close OUTFILE;
