#!/usr/bin/perl

use strict;

my $min_qv = shift @ARGV;
my $pass_percentage = shift @ARGV;
my $inputfile = shift @ARGV;
my $outputfile = shift @ARGV;
my $phredfile = "pipeline/phred.txt";
print "minimum_QV $min_qv\n";
print "pass_QVcheck_percentage $pass_percentage\n";
print "Input_file $inputfile\n";
print "Output_file $outputfile\n";


my %phred = ();
open (IN, "pipeline/phred.txt") || die $!;
while (my $line = <IN>) {
	chomp $line;
	my @vals = split(/\t/, $line);
	$phred{$vals[0]} = $vals[1];
}
close (IN);

#foreach my $i (keys %phred) {
#	print "$i\t$phred{$i}\n";
#}

my $line1 = "";
my $line2 = "";
my $line3 = "";
my $line4 = "";
my $lineid = 0;
open (OUTPUT, ">$outputfile");
open (IN, "$inputfile") || die $!;
while (my $line = <IN>) {
	chomp $line;
	$lineid += 1;
	if ($lineid == 1) {
		$line1 = $line;
	}
	elsif ($lineid == 2) {
		$line2 = $line;
	}
	elsif ($lineid == 3) {
		$line3 = $line;
	}
	elsif ($lineid == 4) {
		$line4 = $line;
		my @base = split(//, $line2);
		my @qv = split(//, $line4);
		my $length = @qv;
		my $pass_base = 0;
		for (my $i = 0; $i < $length; $i += 1) {
			my $qv = $qv[$i];
			my $phredscore = $phred{$qv} - 33;
			#print "$qv\t$phred{$qv}\n";
			if ($phredscore >= $min_qv) {
				$pass_base += 1;
			}
		}
		my $percentage = 0;
		if ($pass_base == 0) {
		}
		else {
			$percentage = ($pass_base / $length) * 100;
		}
		if ($percentage >= $pass_percentage) {
			print OUTPUT "$line1\n$line2\n$line3\n$line4\n";
		}
		$line1 = "";
		$line2 = "";
		$line3 = "";
		$line4 = "";
		$lineid = 0;
	}
	else {
		print "Error_line\n";
	}
}
close (IN);
