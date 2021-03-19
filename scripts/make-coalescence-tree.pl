#!/usr/bin/env perl
# Produce a coalescent tree from simulated seqs in FASTA output
#  e.g., Post-A.fas
#
use Bio::SeqIO;
use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;
use Bio::Tree::Tree;
use Bio::Tree::Node;
use Bio::TreeIO;

my%opts;
&getopts('Vt', \%opts);
my $allele = 1;
#my $last_gen = $opts{g} || 200;
my (@lineages, @coales_events, @gens);
my $tag;
my $last_gen;
#my $sample_size = 0;
while(<>) {
    chomp;
    my @data = split;
#    $data[1] =~ /^Seq_(\d+)_\d+$/;
#    my $gen = $1;
    $last_gen = $data[1];
#    next unless $gen == $last_gen;
    $tag = $data[0];
    my @anc = split /\|/, $data[3];
 #   $sample_size = scalar @anc;
    push @lineages, { # samples, 0  through 19
	allele=>$allele++, # 1 through sample size (e.g., 20)
	ancestors=>\@anc, # 0 through gen size (e.g, 499)
    }
}

#
# make last generation resolved: this solves the missing-allele problem, but is probably wrong
# the question is why alleles don't resolve even at the last generation
#
#my @resolve_last = (0 .. $sample_size - 1);
#for (my $i = 0; $i <= $#lineages; $i++) {
#    push @{$lineages[$i]->{ancestors}}, $resolve_last[$i];
#}


if ($opts{V}) {
    for (my $i = 1; $i <= $last_gen; $i++) {
	print $i, "\t";
	foreach (@lineages) {
	    print $_->{ancestors}->[$i-1], "\t";
	}
	print "\n";
    }
}

my $coal_ct = 0;
my $br_len = 0;
#my @coals;
my $root=Bio::Tree::Node->new(-id=>'root');

&split_lineage(0, \@lineages, $root);
my $tree = Bio::Tree::Tree->new(-id=>$tag, -root=>$root, -nodelete=>1);
#print $tree->as_text("newick"), "\n" if $opts{t};
print $tree->as_text("newick"), "\n";

=a
my $first_coal=shift @coals;
if ($opts{B}) { # before first coalescence
    foreach (@lineages) {
	my @anc = @{$_->{ancestors}};
	for (my $i=0; $i<$first_coal; $i++) {
	    my $pa = "L";
	    my @subs = @anc[0..$i];
	    my $name = join "_", @subs;
	    $pa .= $name;
	    my $ch = $pa."_".$anc[$i+1];
	    print $pa, "\t", $ch, "\n";
	}
    }
}

if ($opts{A}) { # after first coalescence
    foreach (@lineages) {
	my @anc = @{$_->{ancestors}};
	for (my $i=$first_coal; $i<$last_gen; $i++) {
	    my $pa = "L";
	    my @subs = @anc[0..$i];
	    my $name = join "_", @subs;
	    $pa .= $name;
	    my $ch = $pa."_".$anc[$i+1];
	    print $pa, "\t", $ch, "\n";
	}
    }
}
=cut

exit;

# assuming coalescence at the root
# Throws error if not coalesced
sub split_lineage {
    my $begin = shift;
    my $ref = shift;
    my $parent = shift;
    my @set = @$ref;
    for (my $i=$begin; $i <= $last_gen; $i++) {
	my %seen;
	foreach (@set) {   $seen{$_->{ancestors}->[$i]}++;  } # collect all at gen i
	my @mems = keys %seen;
	$br_len ++;
	next if @mems == 1; # single lineage at gen i
	$coal_ct++;
#	push @coals, $i; # e.g, i=288 at gen 289
	print "Coalescent event $coal_ct at ", $i, "\n" if $opts{V}; 
	my $nd = Bio::Tree::Node->new(-id=>"co_".$coal_ct, -branch_length => $br_len);
	$parent->add_Descendent($nd);
	$br_len = 0;
	
	foreach my $pid (@mems) { # e.g., 41, 44, for each id at a split generation
	    my @subset =  grep {$_->{ancestors}->[$i] == $pid} @set;
	    if (@subset > 1) { # two or more alleles
		&split_lineage($i+1, \@subset, $nd);
	    } else { # reached an OTU
		my $ln = shift @subset;
		my $otu = Bio::Tree::Node->new(-id=>"S_".$ln->{allele}, -branch_length=>$last_gen-$i);
		$nd->add_Descendent($otu);
	    }
	}
	return
#	last; # stop at $i; don't proceed further
    } 
}


exit;
