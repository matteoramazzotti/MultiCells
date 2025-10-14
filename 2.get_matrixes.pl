use strict;
use warnings;
use v5.10;
use LWP::UserAgent;
use Data::Dumper;
use IPC::Run 'run';
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);


my $gse_sel_file = "selected_gse.txt";
my $gsm_sel_file = "selected_gsm.txt";
my @gses;
my @gsms;
open(IN,$gse_sel_file);
while(<IN>) {
	chomp;
	push(@gses,$_);
}
close(IN);
open(IN,$gsm_sel_file);
while(<IN>) {
	chomp;
	push(@gsms,$_);
}
close(IN);

my $header;
my @matches;
foreach my $gse (@gses) {

	my $matrix = get_matrix(gse=>$gse);
	chomp $matrix;
	save_matrix(path=>"raw_counts/mat/${gse}_raw.txt",matrix=>$matrix);
	my $matrix_tr = transpose_matrix(matrix=>$matrix);
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
$global_match_tr =~ s/ /\t/g;

open(OUT,">","raw_counts/matrix_2.tsv");
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

sub transpose_matrix {
	my %args = (
		matrix => undef,
		@_
	);
	my $matrix = $args{matrix};
	my $tr_awk_line = '{ 
    for (i=1; i<=NF; i++)  {
      a[NR,i] = $i
    }} NF>p { p = NF } END {
		 for(j=1; j<=p; j++) { 
			str=a[1,j]
		 for(i=2; i<=NR; i++){ str=str" "a[i,j]; } print str } 
	}';

	my @cmd1 = ["awk", $tr_awk_line];
	run @cmd1,"<",\$matrix, ">", \my $stdout;
	return($stdout);
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