#!/usr/bin/perl
# jointgene.com
# biotang, 2017
use strict;
use warnings;
use File::Basename;

# usage
my $this_scp = basename $0;
die "Usage: perl $this_scp in.sam out.sam\n" if(not @ARGV);

# command arguments
my $file_sam = shift @ARGV;
my $file_out = shift @ARGV;

# 
#############################################
#  Main
#############################################
my $runtime = time;
my ($num, $all) = &collect_uniq_mapped($file_sam, $file_out);
my $format_all = &num_readable($all);
my $format_unique = &num_readable($num);
my $format_unique_rate = 100*$num/$all;
#
print  STDERR "total paired-end reads:\t$format_all\n";
print  STDERR "unique mapped reads:\t$format_unique\n";
printf STDERR "unique mapped rate:\t%.2f%%\n", $format_unique_rate;
$runtime = time - $runtime;
print STDERR "runtime: $runtime sec\n";

#############################################
#  Main function
#############################################
sub collect_uniq_mapped{
	my $file_sam = shift @_;
	my $file_out = shift @_;
	#
	my $total_num = 0;
	my $uniq_num  = 0;
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
		print OUT "$line\n" if(exists $hit{'XU'} and $hit{'XU'} ne "00");
		#
		my $bin_flag = &convert_to_bin_flag($hit{'FLAG'});
		if(not &is_second_alignment($bin_flag) ){
			$total_num++;
			next if(&is_read_unmapped($bin_flag) );
			$uniq_num++ if(exists $hit{'XU'} and $hit{'XU'} ne "00");
		}
	}
	#
	close FILE;
	close OUT;
	#
	return($uniq_num, $total_num);
}

#############################################
#  Functions
#############################################
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

sub convert_to_bin_flag{
	return sprintf "%b", $_[0];
} 

sub is_read_unmapped{
	my $mk = 0;# no
	$mk=1 if($_[0]=~m/1[01]{2}$/);# yes
	return $mk;
}

sub is_second_alignment{
	my $mk = 0; # it's the first alignment, not the second alignment
	$mk = 1 if($_[0]=~m/1[01]{8}$/);
	return $mk;
}

sub num_readable{
	my $num = shift @_;
	$num =~ s/(?<=[0-9])(?=(?:[0-9]{3})+\z)/,/g;
	return $num;
}

