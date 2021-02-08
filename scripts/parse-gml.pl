#!/usr/bin/env perl

# 2/8/2021
# Parse TCS gml output

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
	    freq => undef,
	    out_wt => undef,
	    seq_str => undef,
	    iso => undef,
	    label => undef,
	    id => undef,
	};
	$in_node = 1;
    } 

    if ($nodeDepth == 2 && $text eq ']') {
	if ($in_node) {
	    push @nodes, $nd;
	    $in_node = 0;
	}

	if ($in_edge) {
	    push @edges, $edg;
	    $in_edge = 0;
	}
    } 

    if ($text =~ /id\s+(\d+)\s+label\s+\"(.*)\"/) {
	$nd->{id} = $1;
	my $lab = $2;
	$lab =~ s/\s//g;
	$nd->{label} = $lab ? $lab : 'NA';
	$numInNodes ++ unless $lab;
    } 

    if ($text =~ /outgroup\s+weight\s*=\s+(0\.\d+).+frequency=\s*(\d+)\s+(\S.+\d+)\s+.+Sequence\s+.Sequence\s+=([ATCG]+)/) {
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

    if ($text =~ /linestyle\s+.+label\s+\"(\d+)\"/) {
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
exit;


