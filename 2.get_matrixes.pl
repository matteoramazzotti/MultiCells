use strict;
use warnings;
use v5.10;
use LWP::UserAgent;
use Data::Dumper;
use IPC::Run 'run';
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use File::Path qw(make_path);
use FindBin qw($Bin);

use Getopt::Long;

Getopt::Long::Configure ("bundling","no_ignore_case");

GetOptions( 
	'help|h' => \my $help,
	'verbose|v+' => \my $verbose,
	'gsm_list|l=s' => \my $gsm_sel_file,
	'gse_list|L=s' => \my $gse_sel_file,
	'out_folder|O=s' => \my $out_folder
);



if (!$gse_sel_file) {
	$gse_sel_file = "selected_gse.txt";
}
if (!$gsm_sel_file) {
	$gsm_sel_file = "selected_gsm.txt";
}
if (!$out_folder) {
	$out_folder = "$Bin";
}

my @gses;
my @gsms;
open(IN,$gse_sel_file);
while(<IN>) {
	chomp;
	next if ($_ eq "");
	push(@gses,$_);
}
close(IN);
open(IN,$gsm_sel_file);
while(<IN>) {
	chomp;
	next if ($_ eq "");
	push(@gsms,$_);
}
close(IN);

setup_dir("${out_folder}/raw_counts");
setup_dir("${out_folder}/raw_counts/mat");

my $header;
my @matches;
foreach my $gse (@gses) {

	my $matrix = get_matrix(gse=>$gse);
	chomp $matrix;
	save_matrix(path=>"${out_folder}/raw_counts/mat/${gse}_raw.txt",matrix=>$matrix);
	my $matrix_tr = transpose_matrix(matrix=>$matrix);
	$matrix_tr = remove_empty_columns(matrix => $matrix_tr);
	chomp $matrix_tr;
	if (!$header){
		$header = get_header(matrix=>$matrix_tr);
	}
	my $match = match_gsm(matrix=>$matrix_tr,gsms =>\@gsms);
	chomp $match;
	push(@matches,$match);
	sleep(5);
}
my $global_match = join("\n",$header,@matches);
my $global_match_tr = transpose_matrix(matrix=>$global_match);
# $global_match_tr =~ s/ /\t/g;
$global_match_tr = remove_empty_columns(matrix => $global_match_tr);

open(OUT,">","${out_folder}/raw_counts/matrix.tsv");
print OUT $global_match_tr;
close(OUT);

sub get_matrix {
	my %args = (
		gse => undef,
		@_
	);
	my $gse = $args{gse};
	my $ua = LWP::UserAgent->new;
	$ua->timeout(10);
	$ua->env_proxy;
	# my $can_accept = HTTP::Message::decodable;
	my $url = "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=${gse}&format=file&file=${gse}_raw_counts_GRCh38.p13_NCBI.tsv.gz";
	say $url;
	my $text;
	my $response = $ua->get($url);
	gunzip \$response->content => \$text or die "Gunzip failed: $GunzipError\n";
	return($text);
}

sub save_matrix {
	my %args = (
		matrix => "",
		path => ".",
		@_
	);

	open(OUT,">",$args{path});
	print OUT $args{matrix};
	close(OUT);
}

# sub transpose_matrix {
# 	my %args = (
# 		matrix => undef,
# 		@_
# 	);
# 	my $matrix = $args{matrix};
# 	my $tr_awk_line = '{ 
#     for (i=1; i<=NF; i++)  {
#       a[NR,i] = $i
#     }} NF>p { p = NF } END {
# 		 for(j=1; j<=p; j++) { 
# 			str=a[1,j]
# 		 for(i=2; i<=NR; i++){ str=str" "a[i,j]; } print str } 
# 	}';

# 	my @cmd1 = ["awk", $tr_awk_line];
# 	run @cmd1,"<",\$matrix, ">", \my $stdout;
# 	return($stdout);
# }

sub transpose_matrix {
	my %args = ( matrix => undef, @_ );
	my @rows = map { [ split(/\t/, $_, -1) ] } split(/\n/, $args{matrix});

	my $max = 0;
	for my $r (@rows) {
		$max = @$r if @$r > $max;
	}

	my @out;
	for my $col (0 .. $max-1) {
		push @out, join("\t", map { $_->[$col] // "" } @rows);
	}

	return join("\n", @out) . "\n";
}
sub remove_empty_columns {
	my (%args) = @_;
	my $matrix = $args{matrix};

	my @rows = map { [ split(/\t/, $_, -1) ] } split(/\n/, $matrix);
	return $matrix if !@rows;

	my $cols = 0;
	for my $r (@rows) {
		$cols = @$r if @$r > $cols;
	}

	my @keep;
	COL: for my $c (0 .. $cols-1) {
		for my $r (@rows) {
			next if $r->[$c] eq "";
			next if $r->[$c] =~ /^\s*$/;
			$keep[$c] = 1;
			next COL;
		}
		$keep[$c] = 0;
	}

	my @clean;
	for my $r (@rows) {
		my @new = map { $r->[$_] } grep { $keep[$_] } (0 .. $cols-1);
		push @clean, join("\t", @new);
	}

	return join("\n", @clean) . "\n";
}
sub get_header {
	my %args = (
		matrix => undef,
		@_
	);
	my $matrix = $args{matrix};
	my @tmp = split("\n",$matrix);
	return($tmp[0]);
}

sub match_gsm {
	my %args = (
		matrix => undef,
		gsms => undef,
		@_
	);
	my $matrix = $args{matrix};
	my @gsms = @{$args{gsms}};
	my $regex = join("|",@gsms);
	my @cmd1 = ["grep","--perl-regex", $regex];
	run @cmd1,"<",\$matrix, ">", \my $stdout;
	return($stdout);
}

sub setup_dir{
	my $dir = $_[0];
	unless (-d $dir) {
		print "Directory $dir does not exist. Creating $dir...\n";
		eval {
			make_path($dir);
			print "Directory $dir created successfully.\n";
		};
		if ($@) {
			die "Failed to create directory $dir: $@\n";
		}
	} else {
		print "Directory $dir already exists, proceding to query\n";
	}
}