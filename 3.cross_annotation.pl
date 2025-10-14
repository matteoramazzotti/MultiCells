use strict;
use warnings;

my $file = $ARGV[0];
my @genes;
my $genes = {};
my $header;


open(IN,$file);
while (<IN>) {
	chomp;
	
	if ($_ =~ /GeneID/) {
		$header = $_;
		next;
	}
	my @d = split("\t",$_);
	push(@genes,$d[0]);
	$genes->{$d[0]}->{symbol} = "";
	$genes->{$d[0]}->{line} = join("\t","MYGENEGOESHERE",@d[1..$#d]);
}
close(IN);

my $seen;
my $db_file = "dbs/hsa.gene_info";
open(IN,$db_file);
while (<IN>) {
	chomp;
	next if ($_ =~ /GeneID/);
	my @d = split("\t",$_);
	$seen -> {$d[2]} = 0 if (!$seen -> {$d[2]});
	next if ($seen -> {$d[2]} > 1);
	if ($genes->{$d[1]}) {
		$seen -> {$d[2]}++;
		if ($seen -> {$d[2]} > 1) {
			delete $genes->{$d[1]};
			next;
		}
		$genes->{$d[1]}->{symbol} = $d[2];
		$genes->{$d[1]}->{line} =~ s/MYGENEGOESHERE/$d[2]/;
	}
}
close(IN);

print STDOUT $header,"\n";

foreach my $g (@genes) {
	next if (!$genes->{$g});
	if ($genes->{$g}->{line} =~ /MYGENEGOESHERE/ ) {
		$genes->{$g}->{line} =~ s/MYGENEGOESHERE/$g/
	}
	print STDOUT $genes->{$g}->{line},"\n";
}