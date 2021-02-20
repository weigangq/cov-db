#!/usr/bin/env perl

# 2/8/2021
# Parse TCS gml output
# Ref URL: https://savage.net.au/Ron/html/Fancy.Matching.of.Delimited.Text.html

use strict;
use Data::Dumper;
use Text::Balanced::Marpa;
use warnings;

my $parser = Text::Balanced::Marpa -> new(
                open  => ['['],
                close => [']'],
);

my $string = '';
while(<>) {
    chomp;
    $_ =~ s/\R//g; # remove generalized newline
    $string .= $_;
}

my @nodes;
my @edges;
my $result = $parser->parse(text => \$string);

die "parsing failed\n" if $result;
my($attributes);
my($indent);
my($text);
my $in_node = 0;
my $in_edge = 0;
my $nodeDepth;
my $nd;
my $edg;
my $numInNodes = 0;

for my $node ($parser -> tree -> traverse){
    next if ($node -> is_root);
    $attributes = $node -> meta;
    $text       = $$attributes{text};
    $text       =~ s/^\s+//;
    $text       =~ s/\s+$//;
    $nodeDepth = $node->depth();

    if ($nodeDepth == 2 && $text =~ /node/) {
	$nd = {
	    freq => 1,
	    out_wt => 0,
	    seq_str => 'NA',
	    iso => undef,
	    label => undef,
	    id => undef,
	};
	$in_node = 1;
    } 

    if ($nodeDepth == 2 && $text eq ']') {
	if ($in_node) {
	    push @nodes, $nd if $nd->{label} ne 'NA';
	    $in_node = 0;
	}

	if ($in_edge) {
	    push @edges, $edg;
	    $in_edge = 0;
	}
    } 

    if ($text =~ /id\s+(\d+)\s+label\s+\"(.*)\"/) {
	$nd->{id} = $1;
	my $lab = $2; # " "  for InNodes; "" for unlabeled nodes (not sure what)
	if ($lab) {
	    $lab =~ s/\s//g;
	    $nd->{label} = $lab ? $lab : 'InNode' . $numInNodes++;
#	    $numInNodes ++ unless $lab;
	} else {
	    $nd->{label} = "NA";
	}
    }

    if ($text =~ /outgroup\s+weight\s*=\s+(0\.\d+).+frequency=\s*(\d+)\s+(.+)\"\s+.+Sequence\s+.Sequence\s+=([ATCG]+)/) {
	$nd->{freq} = $2;
	$nd->{out_wt} = $1;
	$nd->{seq_str} = $4;
	$nd->{iso} = [ split(/\s+/, $3) ];
#	$nd->{iso} = $3;
    } 

    if ($nodeDepth == 2 && $text =~ /edge/) {
	$edg = {
	    node_from => undef,
	    node_to => undef,
	    nt_from => undef,
	    nt_to => undef,
	    label => undef,
	};
	$in_edge = 1;
    }

    if ($text =~ /source\s+(\d+)\s+target\s+(\d+)/) {
	$edg->{node_from} = $1;
	$edg->{node_to} = $2;
    }

    if ($text =~ /Changes\s+\"([ATCG])\s+([ATCG])\"/) {
	$edg->{nt_from} = $1;
	$edg->{nt_to} = $2;
    }

    if ($text =~ /linestyle\s+.+label\s+\"(\S+)\"/) {
	$edg->{label} = $1;
    }

=a
    $indent     = $node -> depth - 1;
    if (length($text)) {	
	print "Depth ", $indent + 1, ":", "\t" x $indent, "$text\n";
    }  else {
	print "Depth ", $indent + 1, ":\n";
    }
=cut
}

#print Dumper(\@nodes, \@edges);
warn "Num of nodes\t", scalar(@nodes), "\nNum of innodes\t", $numInNodes, "\nNum of edges\t", scalar(@edges), "\n";

my %nodes;
foreach my $nd (@nodes) { $nodes{$nd->{id}} = $nd }

print join "\t", ('edge_id', 'id_from', 'id_to', 'num_iso_from', 'num_iso_to', 'label_from', 'label_to', 'seq_from', 'seq_to', 'nt_from', 'nt_to');
print "\n";
foreach my $edg (@edges) {
    print join "\t", ($edg->{label},
		      $edg->{node_from},
		      $edg->{node_to},
		      $nodes{$edg->{node_from}}->{freq},
		      $nodes{$edg->{node_to}}->{freq},
		      $nodes{$edg->{node_from}}->{label},
		      $nodes{$edg->{node_to}}->{label},
		      $nodes{$edg->{node_from}}->{seq_str},
		      $nodes{$edg->{node_to}}->{seq_str},
		      $edg->{nt_from},
		      $edg->{nt_to}
    );
    print "\n";
}

exit;


