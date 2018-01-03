#!/usr/bin/perl
# biotang, July 2014
use strict;
use warnings;
use File::Basename;
use 5.010001;

#############################################
#  Usage
#############################################
my $this_scp = basename $0;
die "Usage: $this_scp in.sam out.sam\n" if(not @ARGV);

#############################################
#  Command options
#############################################
my $file_sam = shift @ARGV;
my $file_out = shift @ARGV;

#############################################
#  Main
#############################################
my $runtime = time;
my ($num, $all) = &mark_uniq_mapped($file_sam, $file_out);
my $format_all = &num_readable($all);
my $format_unique = &num_readable($num);
my $format_unique_rate = 100*$num/$all;

#
printf STDERR "total paired-end reads:\t%s\n", $format_all;
printf STDERR "unique mapped reads:\t%s\n", $format_unique;
printf STDERR "unique mapped rate:\t%.2f%%\n", $format_unique_rate;

$runtime = time - $runtime;
print STDERR "runtime: $runtime sec\n";

#############################################
#  Main sub
#############################################
sub mark_uniq_mapped{
	my $file_sam = shift @_;
	my $file_out = shift @_;
	#
	my $uniq_num  = 0; #
	my $uniq_read = 0;
	my $total_num = 0; #
	my $ck = "";
	my @pair = ();
	my @line = ();
	#
	open FILE, $file_sam or die "";
	open OUT, ">$file_out" or die "";
	while(my $line=<FILE>){
		chomp $line;
		# header lines
		if($line=~m/^@/){
			print OUT "$line\n";
			next;
		}
		# read sam lines
		my %hit = &read_one_sam_line($line);
		
		#
		if($hit{'QNAME'} ne $ck){
			# process
			$total_num = $total_num + 2; # only for paired-end
			if($ck){
				my ($mk1, $mk2, $map1, $map2) = &process_one_read_pair(@pair);
				if($mk1==1 or $mk2==1){
					$uniq_num = $uniq_num + 2;
					$uniq_read++ if($map1==1);
					$uniq_read++ if($map2==1);
				}
				#
				my $str = $mk1.$mk2;
				foreach my $outline (@line){
					print OUT "$outline\tXU:Z:$str\n";
				}
			}
			# clean
			$ck = $hit{'QNAME'};
			@pair = ();
			@line = ();
		}
		push @pair, \%hit;
		push @line, $line;
	}
	# process the last one
	$total_num = $total_num + 2; # only for paired-end
	if($ck){
		my ($mk1, $mk2, $map1, $map2) = &process_one_read_pair(@pair);
		if($mk1==1 or $mk2==1){
			$uniq_num = $uniq_num + 2;
			$uniq_read++ if($map1==1);
			$uniq_read++ if($map2==1);
		}
		#
		my $str = $mk1.$mk2;
		foreach my $outline (@line){
			print OUT "$outline\tXU:Z:$str\n";
		}
	}
	#
	close FILE;
	close OUT;
	#
	return($uniq_read, $total_num);
}

#
sub process_one_read_pair{
	my $as1  = 0;
	my $xs1  = 0;
	my $as2  = 0;
	my $xs2  = 0;
	my $map1 = 0;# 0 mean the first read of pair is not mapped
	my $map2 = 0;# 0 mean the second read of pair is not mapped
	foreach my $hit (@_){
		my $bin_flag = &convert_to_bin_flag($$hit{'FLAG'});
		next if(&is_read_unmapped($bin_flag) );
		next if(&is_second_alignment($bin_flag) );
		if(&is_first_segment($bin_flag)){
			$as1  = $$hit{'AS'};
			$xs1  = $$hit{'XS'};
			$map1 = 1;
		} elsif(&is_second_segment($bin_flag)){
			$as2  = $$hit{'AS'};
			$xs2  = $$hit{'XS'};
			$map2 = 1;
		}
	}
	my $uniq_mk1 = 0;# 0 mean the first read of pair is not unique mapped
	my $uniq_mk2 = 0;# 0 mean the second read of pair is not unique mapped
	$uniq_mk1 = 1 if($as1>$xs1);
	$uniq_mk2 = 1 if($as2>$xs2);
	return($uniq_mk1, $uniq_mk2, $map1, $map2);
}

#
sub convert_to_bin_flag{
	return sprintf "%b", $_[0];
} 

#
sub is_read_unmapped{
	my $mk = 0;# no
	$mk=1 if($_[0]=~m/1[01]{2}$/);# yes
	return $mk;
}

#
sub is_first_segment{
	my $mk = 0; # it's the first alignment, not the second alignment
	$mk = 1 if($_[0]=~m/1[01]{6}$/);
	return $mk;
}

#
sub is_second_segment{
	my $mk = 0; # it's the first alignment, not the second alignment
	$mk = 1 if($_[0]=~m/10[01]{6}$/);
	return $mk;
}

#
sub is_second_alignment{
	my $mk = 0; # it's the first alignment, not the second alignment
	$mk = 1 if($_[0]=~m/1[01]{8}$/);
	return $mk;
}

#
sub read_one_sam_line{
	my $line = shift @_;
	my %hit  = ();
	my @col = split /\t/, $line;
	$hit{'QNAME'} = shift @col;
	$hit{'FLAG'}  = shift @col;
	$hit{'RNAME'} = shift @col;
	$hit{'POS'}   = shift @col;
	$hit{'MAPQ'}  = shift @col;
	$hit{'CIGAR'} = shift @col;
	$hit{'RNEXT'} = shift @col;
	$hit{'PNEXT'} = shift @col;
	$hit{'TLEN'}  = shift @col;
	$hit{'SEQ'}   = shift @col;
	$hit{'QUAL'}  = shift @col;
	foreach my $col (@col){
		my ($tag,$type,$value) = split /:/, $col;
		$hit{$tag} = $value;
	}
	return %hit;
}

sub num_readable{
	my $num = shift @_;
	$num =~ s/(?<=[0-9])(?=(?:[0-9]{3})+\z)/,/g;
	return $num;
}
