#!/usr/bin/perl
use warnings;
use strict;

local $/ ="\n>";

my %id;
while(<>) {
    my ($prod,$protein_id,$state);
    $state = '';
    my $i =0;
    for my $line ( split(/\n/,$_) ) {
	next if $line =~ /^>/ && $line =~ /^Feature/;
	my @row = split(/\t/,$line);
#	warn $i++, " ", $line," ", scalar @row, "\n";

	if( @row == 3 ) {
	    if( $state eq 'CDS' ) {
		$id{$protein_id} = $prod;
#		print join("\t", $protein_id, $prod), "\n";
	    }
	    $state = pop @row;
	} elsif( @row == 5 ) {
	    for ( 1..3 ) { shift @row }
	    my ($qualifier, $value)  = @row;
	    if( $qualifier eq 'protein_id') {
		(undef,undef,$protein_id) = split(/\|/,$value);
	    } elsif( $qualifier eq 'product') {
		$prod = $value;
	    }
	}	    
    }
    # fencepost
#    print join("\t", $protein_id, $prod), "\n" if defined $protein_id;
}

for my $acc ( sort keys %id ) {
    print join("\t", $acc, $id{$acc}),"\n";
}
	
