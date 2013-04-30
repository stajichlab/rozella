#!/usr/bin/perl -w
use strict;

my $dir = shift || ".";

my %counts;
my %term2info;
my %totals;
opendir(IN,$dir) || die "$dir: $!";
for my $d (readdir(IN) ) {
    next unless $d =~ /_GOSlim_count/;
    my ($nm) = split(/_/,$d);
    opendir(DIR,$d) || die $!;
    for my $file (readdir(DIR) ) {
	next unless $file =~ /(\d+)\.(mf|cc|bp)\.txt$/;
	my ($ont) = $2;
	open(my $fh => "$d/$file") || die $!;
	while(<$fh>) { 
	    chomp;	    
	    my ($term,$desc,$count) = split(/\t/,$_);
	    next if $desc eq 'molecular_function' || 
		$desc eq 'biological_process' || 
		$desc eq 'cellular_component';
	    # skip top level catch-all

	    $counts{$ont}->{$nm}->{$term} += $count;
	    $totals{$ont}->{$nm} += $count;
	    $term2info{$ont}->{$term} = $desc unless exists $term2info{$term};
	}
    }
}

my @onts = sort keys %counts;
for my $ont ( @onts ) {
    open(my $fh => ">$ont.counts.tab") || die $!;
    open(my $ggfh => ">$ont.ggplot.tab") || die $!;
    my @nms = sort keys %{$counts{$ont}};
#    print $fh join("\t", qw(GO_ID GO_TERM), @nms), "\n";
    print $fh join("\t", qw(GO_TERM), @nms), "\n";
    print $ggfh join("\t", qw(GO_ID GO_Term Count Percentage Samples)), "\n";
    for my $term ( sort keys %{$term2info{$ont}}) {
        next if $term eq 'GO:0030533';
	print $fh join("\t", #$term, 
		       $term2info{$ont}->{$term}, 
		       map { sprintf("%.2f",
				     100 * 
				     ($counts{$ont}->{$_}->{$term} || 0) / 
				     $totals{$ont}->{$_}) } @nms),"\n";
	for my $nm (@nms ) {
	    print $ggfh join("\t",$term, $term2info{$ont}->{$term}, 
			     $counts{$ont}->{$nm}->{$term} || 0,
			     sprintf("%.2f",100 *   ($counts{$ont}->{$nm}->{$term} || 0) / 
				     $totals{$ont}->{$nm}),
			     $nm),"\n";
	}
    }
}
