use strict;
use warnings;
use v5.10;
use Data::Dumper;

use XML::Simple;
use Encode;
use Encode::Guess;

use File::Path qw(make_path);
use FindBin qw($Bin);
use IPC::Run 'run';

use Getopt::Long;

Getopt::Long::Configure ("bundling","no_ignore_case");

GetOptions( 
	'help|h' => \my $help,
	'verbose|v+' => \my $verbose,
	'query|q=s' => \my $query,
	'query-file|Q=s' => \my $query_file,
	'query-prefix|p=s' => \my $query_prefix, # string (when no quey file is provided) or location in table (when tsv file is provided) or none (query number)
	'max-result|m=i' => \my $max_res, # max number of series
	'check-counts|c' => \my $check_counts,
	'out_folder|O=s' => \my $out_folder
);

if (!$query && !$query_file) {
	say STDERR "No query and no query file provided";
	die;
}
if (!$out_folder) {
	$out_folder = "$Bin";
}



my @queries;

my $n = 0;

if ($query) {
	my $item = {
		query => $query
	};
	if ($query_prefix) {
		$item->{id} = $query_prefix;
	} else {
		$item->{id} = "query_$n";
	}
	push(@queries,$item);
}

if ($query_file) {
	open(IN,$query_file);
	while(<IN>) {
		chomp;
		next if ($_=~/^#/);
		$n++;
		my @data = split("\t",$_);
		my $item = {};
		$item->{query} = $data[0];
		if (scalar @data > 1) {
			$item->{id} = $data[1];
		} else {
			$item->{id} = "query_$n";
		}
		push(@queries,$item);
	}
	close(IN);
}


setup_dir("$out_folder/indexes");


foreach my $qObj (@queries) {
	open(OUT,">","$out_folder/indexes/query_$qObj->{id}.tsv");
	open(LOG,">","$out_folder/indexes/query_$qObj->{id}.log");
	say LOG "Performing: $qObj->{query}";
	say STDERR "Performing: $qObj->{query}";
	my $q = $qObj->{query};
	
	my $stdout = queryForSeries(query=>$q);

	if (!$stdout) {
		say LOG "Skipped: $q\n";
		say STDERR "Skipped: $q\n";
		sleep(10);
		next;
	}
	my $xml = XML::Simple->new();
	my $data = $xml->XMLin($stdout);
	my @main_struct;
	if (ref $data->{DocumentSummary} eq "HASH") {
		push(@main_struct,$data->{DocumentSummary});
	}
	if (ref $data->{DocumentSummary} eq "ARRAY") {
		@main_struct = @{$data->{DocumentSummary}};
	}
	my $n_of_ser = 0;
	foreach my $item (@main_struct) {
		print LOG "\tSeries: $item->{Accession}\n";
		
		#test if counts are available
		# my $base_url = url="www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=$item->{Accession}&format=file&file=$item->{Accession}_raw_counts_GRCh38.p13_NCBI.tsv.gz";
		my $counts_available = "true";
		if ($check_counts) {
			say STDERR "$item->{Accession}: Checking for counts...";
			$counts_available = test_url_for_counts(series=>$item->{Accession});
			sleep(10);
		}
		my @plat = split(";",$item->{GPL});
		@plat = map {"GPL$_"} @plat;
		

		my @samples = ();
		my $gsm_d = {};

		my @main;
		if (ref $item->{Samples}->{Sample} eq "HASH") {
			push(@main,$item->{Samples}->{Sample});
		}
		if (ref $item->{Samples}->{Sample} eq "ARRAY") {
			@main = @{$item->{Samples}->{Sample}};
		}

		
		foreach my $sample (@main) {
			push(@samples,$sample->{Accession});
			$gsm_d->{$sample->{Accession}} = {
				sample => $sample->{Accession},
				title => $sample->{Title},
				series => $item->{Accession},
				project => $item->{BioProject},
				platforms => join(";",@plat),
				organism => $item->{taxon},
				counts_available => $counts_available,
				gds_type => $item -> {gdsType}
			}
		}
		
		my $stdout2 = queryForSamples(samples=>\@samples);
		sleep(10);
		# if (!$stdout2) {
		# 	say LOG "Skipped: ",join(" ",@samples);
		# 	say STDERR "Skipped: ",join(" ",@samples);
		# 	foreach my $sample (@samples) {
		# 		my @line = map {$_ ? $_ : ""} ($sample->{Accession},$gsm_d->{$sample->{Accession}}->{srr},$gsm_d->{$sample->{Accession}}->{title},$gsm_d->{$sample->{Accession}}->{series},$gsm_d->{$sample->{Accession}}->{project},$gsm_d->{$sample->{Accession}}->{platforms},$gsm_d->{$sample->{Accession}}->{organism});
		# 		say OUT join("\t",@line);
		# 		# say OUT join("\t",$sample,"",$gsm_d->{$sample}->{title},$gsm_d->{$sample}->{series},$gsm_d->{$sample}->{project},$gsm_d->{$sample}->{platforms},$gsm_d->{$sample}->{organism});
		# 	}
		# 	next;
		# }
		my @p_details = split("\n",$stdout2);

		if (scalar(@p_details) == 0){
			say LOG "Skipped: ",join(" ",@samples);
			say STDERR "Skipped: ",join(" ",@samples);
			foreach my $sample (@main) {
				my @line = map {$_ ? $_ : ""} makeLine(sample=>$gsm_d->{$sample->{Accession}});
				say OUT join("\t",@line);
				# say OUT join("\t",$sample->{Accession},"",$gsm_d->{$sample->{Accession}}->{title},$gsm_d->{$sample->{Accession}}->{series},$gsm_d->{$sample->{Accession}}->{project},$gsm_d->{$sample->{Accession}}->{platforms},$gsm_d->{$sample->{Accession}}->{organism});
			}
			next;
		}
		foreach my $pd (@p_details) {
			next if ($pd =~ "BioProject");
			my @pdata = split(",",$pd);
			
			my $gsm = $pdata[29] ? $pdata[29] : $pdata[11] ? $pdata[11] : "";
			
			my $strategy = $pdata[12] ? $pdata[12] : "";
			$gsm_d->{$gsm}->{strategy} = $strategy;

			my $platform_name = $pdata[19] ? $pdata[19] : "";
			$gsm_d->{$gsm}->{platform_name} = $platform_name;

			if ($gsm_d->{$gsm}->{srr}) {
				$gsm_d->{$gsm}->{srr} .= ";$pdata[0]";
			} else {
				$gsm_d->{$gsm}->{srr} = $pdata[0] ? $pdata[0] : "";
			}

		}
		foreach my $sample (@main) {
			my @line = map {$_ ? $_ : ""} makeLine(sample=>$gsm_d->{$sample->{Accession}});
			say OUT join("\t",@line);
		}
		$n_of_ser++;
		if ($max_res) {
			last if ($n_of_ser == $max_res);
		}
	}
	close(OUT);
	close(LOG);
	if (scalar @queries > 1) {
		sleep(10);
	}
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


sub test_url_for_counts {
	use LWP::UserAgent;
	my %args = (
		@_
	);

	my $gse = $args{series};
	my $url = "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=${gse}&format=file&file=${gse}_raw_counts_GRCh38.p13_NCBI.tsv.gz";
	my $ua = LWP::UserAgent->new;
	$ua->timeout(10);
	$ua->env_proxy;

	my $response = $ua->head($url);

	if ($response->is_success) {
		return("true");
	} else {
		return("false");
	}
}

sub queryForSeries{
	my %args = (
		@_
	);
	my $q = $args{query};
	my @cmd1 = ["esearch","-db", "gds","-query", $q];
	my @cmd2 = ["efetch","-format", "docsum"];

	run @cmd1,"|",@cmd2, ">", \my $stdout, "2>", \my $stderr;

	my $enc = guess_encoding($stdout, qw/utf-8 iso-8859-1/);
	if (ref($enc)) {
		$stdout = $enc->decode($stdout);
	} else {
		$stdout = Encode::decode('iso-8859-1', $stdout);
	}
	return($stdout);
}

sub queryForSamples{
	my %args = (
		@_
	);
	my @samples = @{$args{samples}};
	my $q_samples = join(" OR ",@samples);
	my @cmd1 = ["esearch","-db", "sra","-query", $q_samples];
	my @cmd2 = ["efetch","-format", "runinfo"];
	run @cmd1,"|",@cmd2, ">", \my $stdout2,"2>", \my $stderr2;
	my $enc = guess_encoding($stdout2, qw/utf-8 iso-8859-1/);
	if (ref($enc)) {
		$stdout2 = $enc->decode($stdout2);
	} else {
		$stdout2 = Encode::decode('iso-8859-1', $stdout2);
	}
	return($stdout2);
}

sub makeLine {
	my %args = (
		@_
	);

	my $sample = $args{sample}; 

	my @arr = (
		$sample->{title},
		$sample->{series},
		$sample->{sample},
		$sample->{counts_available},
		$sample->{gds_type},
		$sample->{project},
		$sample->{organism},
		$sample->{platform_name},
		$sample->{platforms},
		$sample->{srr}
	);
	return(@arr);
}