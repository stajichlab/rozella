#!/usr/bin/perl -w
use strict;

=head1 NAME

gff3_to_gbktbl - Convert GFF3 files into TBL files for GenBank

=head1 USAGE

gff3_to_gbktbl -p product_names -s skip > out.tbl

=head1 DESCRIPTION

=head1 AUTHOR

Jason Stajich - jason.stajich[at]ucr.edu

=cut

use Bio::DB::SeqFeature::Store;
use GO::AppHandle;
use Getopt::Long;

my %GO_components = ('F' => 'go_function',
		     'C' => 'go_component',
		     'P' => 'go_process');
my $godbname = 'GO';
my $dbhost = '';
my $dbuser = '';
my $dbpass = '';
my $gbuser = $ENV{'USER'};
my $gbpass = $ENV{'MYSQL_PASS'};
my %abinit_versions = ('snap' => 'SNAP:2010-07-28',
			 'genemark' => 'GenemarkHMM:2.3e',
			 'augustus' => 'Augustus:2.6.1');

my $gffdb = '';

my $product_names;
my $debug = 0;
my $assoc_file; # GO Association file
my $interpro; # InterPro output
my $hypothetical = 'hypothetical protein';
my ($skip_names,$more_products,$output);
GetOptions(
    'p|product:s'       => \$product_names,
    'mp|moreproducts:s' => \$more_products,
    'a|association:s'   => \$assoc_file,
    'ipr:s'             => \$interpro,
    'v|debug!'          => \$debug,
    's|skip:s'          => \$skip_names,
    'o|output:s'        => \$output,
    'host|dbhost:s'     => \$dbhost,
    'dbuser:s'          => \$dbuser,  # for GO mysqldb
    'dbpass:s'          => \$dbpass,  # for GO mysqldb
    'godb:s'            => \$godbname,
    'gffdb:s'           => \$gffdb, # for features in Bio:DB::SeqFeature database
    'gbuser:s'          => \$gbuser, # if you have different credentials for GO db
    'gbpass:s'          => \$gbpass, # for GFFDB  
    );

my $ofh;
if( $output ) {
    open($ofh => ">$output") || die "cannot open $output: $!";
} else {
    $ofh = \*STDOUT;
}

my %skip_these;
if( $skip_names && -f $skip_names ){
    open(my $fh => $skip_names) || die "$skip_names: $!";
    while(<$fh>) {
	chomp;
	$skip_these{$_}++;
	if( s/-R[A-Z]$// ) {
	    $skip_these{$_}++;
	}
    }
}

my %GOterms;
if( $assoc_file ) {
    open(my $fh => $assoc_file) || die $!;
    while(<$fh>) {
	chomp;
	my ($ctr,$mRNA,$gene, undef,$go, $source, $evidence,undef, $component) = split(/\t/,$_);
	$GOterms{$mRNA}->{$component}->{$go} = [ $source, $evidence];
    }
}

my $go_apph = GO::AppHandle->connect(-dbname=> $godbname,
				     -dbhost => $dbhost,
				     -dbuser => $dbuser,
				     -dbauth => $dbpass,);

my %prodnames;

if( $product_names ) {
    open(my $fh => $product_names) || die $!;
    while(<$fh>) {
	chomp;
	my ($gene_name, $hit, $pid, $evalue, $product, $ec) = split(/\t/,$_);
#	$gene_name =~ s/-R[A-Z]$//;
	$prodnames{$gene_name} = [ $hit, $pid, $evalue, $product, $ec];
    }
}

my %more_prodnames;
if( $more_products ) {
    open(my $fh => $more_products) || die $!;
    while(<$fh>) {
	chomp;
	my ($gene_name, $product, $domain_count, $ec_num, 
	    $kegg, $metacyc, $reactome, $unipath) = split(/\t/,$_);
#	$gene_name =~ s/-R[A-Z]$//;
	$more_prodnames{$gene_name} = [ $product, $domain_count, $ec_num,
					$kegg, $metacyc, $reactome, $unipath];
    }
}
my %iprdomains;
if( $interpro ) {
    open(my $fh => $interpro) || die $!;
    while(<$fh>) {
	chomp;
	my ($locus,$md5,$len,$analysis,$domain,$desc,$qstart,$qend,$score,$status,
	    $date, $ipracc, $iprdesc, $go, $pathway) = split(/\t/,$_);
	if( $ipracc ) {
	    $desc =~ s/domain profile\.//;
	    next if ! length $desc ;
	    $iprdomains{$locus}->{$ipracc} = [$iprdesc, $analysis, $domain, $desc];
	}
    }
}

my $dbh = Bio::DB::SeqFeature::Store->new
    (-adaptor   => 'DBI::mysql',
     -dsn       => "dbi:mysql:host=$dbhost;dbname=$gffdb",
     -user      => $gbuser,
     -password  => $gbpass);

my $iter = $dbh->get_seq_stream(-type => 'scaffold:chromosome');
my @ids;

while( my $chrom = $iter->next_seq ) {
    push @ids, $chrom->display_id;
}
my %seen_products;
for my $chrom ( @ids ) {
    my $segment = $dbh->segment($chrom);
    my @genes = $segment->features(-type => ["gene"]);

    next if ! @genes;
    print $ofh ">Feature $chrom\n";
    for my $gene (@genes) {
	my ($gstart,$gend) = ($gene->start,$gene->end);
	($gstart,$gend) = ($gend,$gstart) if $gene->strand < 0;
	my $locus_name = $gene->display_id;
	next if( $skip_these{$locus_name} );
	print $ofh join("\t",$gstart,$gend, "gene"), "\n";
	print $ofh join("\t",'','','','locus_tag',$locus_name),"\n";
# These aren't great Notes
#	my ($note) = $gene->get_tag_values('Note');
#	print $ofh join("\t",'','','',"note", $note),"\n" if $note !~ /Protein of unknown/;
	for my $mRNA (sort { $a->strand * $a->start <=> $b->strand * $b->start} 
		      $gene->get_SeqFeatures(-type => 'mRNA') ) {
	    
	    my $first = 0;
	    for my $exon ( sort { $a->strand * $a->start <=> 
				      $b->strand * $b->start}
			   $mRNA->get_SeqFeatures(-type => 'exon') ) {
		my ($start,$end) = ($exon->start,$exon->end);
		($start,$end) = ($end,$start) if $exon->strand < 0;
		print $ofh join("\t",$start, $end, ! $first++ ? 'mRNA' : ''), "\n";
	    }
	    my $mRNAName = $mRNA->display_id;
	    print $ofh join("\t",'','','',"transcript_id","$mRNAName"),"\n";
#	    print $ofh join("\t",'','','',"protein_id","gnl|StajichUCR|$locus_name"),"\n";
# These aren't great Notes
#	    ($note) = $mRNA->get_tag_values('Note');
#	    print $ofh join("\t",'','','',"note", $note),"\n" if $note !~ /Protein of unknown/;
	    my $product = $more_prodnames{$mRNAName}->[0] || 'hypothetical protein';
	    if( exists $prodnames{$mRNAName} ) {
		$product = $prodnames{$mRNAName}->[3];
	    }
	    print $ofh join("\t", '','','', 'product', $product),"\n";
	    $first = 0;
	    for my $exon ( sort { $a->strand * $a->start <=> 
				      $b->strand * $b->start} 
			   $mRNA->get_SeqFeatures(-type => 'CDS') ) {
		my ($start,$end) = ($exon->start,$exon->end);
		($start,$end) = ($end,$start) if $exon->strand < 0;
		print $ofh join("\t",$start, $end,  ! $first++ ? 'CDS' : ''), "\n";
	    }
#	    print $ofh join("\t",'','','',"protein_id","gnl|StajichUCR|$locus_name"),"\n";	    
	    print $ofh join("\t",'','','',"transcript_id","$mRNAName"),"\n";	    

	    print $ofh join("\t", '','','', 'product', $product),"\n";
	    my (%xref,%ec_numbers);
#	    my @names = qw(KEGG MetaCyc Reactome UniPath);
#	    my $i = 3;
#	    for my $n ( @names ) {
#		if( $more_prodnames{$mRNAName}->[$i] ) {
#		    push @{$xref{$n}}, ( map { sprintf("%s:%s",
#						       $n, $_) } 
#					 split(',',$more_prodnames{$mRNAName}->[$i]));
#		}
#	    }
	    if( $more_prodnames{$mRNAName}->[2]) {
		for my $ec ( split(',',$more_prodnames{$mRNAName}->[2]) ) {
		    $ec_numbers{$ec}++;
		}
	    }
	    
	    if( exists $prodnames{$mRNAName} ) {
		if( $prodnames{$mRNAName}->[4] ) {
		    $ec_numbers{$prodnames{$mRNAName}->[4]}++;
		}
		my (undef,$spacc) = split(/\|/,$prodnames{$mRNAName}->[0]);
		print $ofh join("\t", '','','','inference', 
			   sprintf("similar to AA sequence:UniProt:%s",$spacc)), "\n";
	    }
	    for my $ec ( sort keys %ec_numbers ) {
		print $ofh join("\t",'','','',"EC_number",$ec),"\n"; 
	    }

	    if( my ($alias) = $mRNA->get_tag_values('Alias') ) {
		my $prog;
		if( $alias =~ /(snap|genemark|augustus)/ ) {
		    $prog = $1;
		    if( exists $abinit_versions{$prog} ) {
			print $ofh join("\t", '','','','inference', 
				   sprintf("ab initio prediction:%s",
					   $abinit_versions{$prog})),"\n";
		    } else {
			warn("no prog version stored for $alias\n");
		    }		
		} else {
		    warn("no prog version stored for $alias\n");
		}
	    }

	    if( exists $GOterms{$mRNAName} ) {
		my %seen;
		while( my ($comp_key,$component) = each %GO_components ) {
		    if( exists $GOterms{$mRNAName}->{$comp_key} ) {
			for my $term_acc ( keys %{$GOterms{$mRNAName}->{$comp_key}} ) {
			    my ($term_source, $term_ev) = 
				@{$GOterms{$mRNAName}->{$comp_key}->{$term_acc}};
			    
#			    my $go_terms = $go_graph->term_query({acc=>$term_acc});
			    my $t = $go_apph->get_term({acc=>$term_acc});
	#		    for my $t ( @{$go_terms || []} ) {
			    next if $seen{$t->acc}++;
			    my $term_desc = $t->name;
			    my $t_acc = $t->acc;
			    print $ofh join("\t", '','','',
				       $component, 
				       join('|', $term_desc,$t_acc,'',$term_ev)),"\n";
#			}
			}
		    }
		}
	    }
	    if( exists $iprdomains{$mRNAName} ) {
		for my $ipracc ( keys %{$iprdomains{$mRNAName}} ) {
		    push @{$xref{'InterPro'}}, 
		    sprintf("InterPro:%s",$ipracc);
		}
	    }
	    for my $xreftype ( sort keys %xref ) {
		for my $xref ( @{$xref{$xreftype}} ) {
		    print $ofh join("\t", '','','', 'db_xref',$xref),"\n";
		}
	    }
	}
	last if $debug;
    }
}
