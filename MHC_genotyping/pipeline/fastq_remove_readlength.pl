#!/usr/bin/perl

use strict;

my $min_length = shift @ARGV;
my $inputfile = shift @ARGV;
my $outputfile = shift @ARGV;
print "minimum_length $min_length\n";
print "Input_file $inputfile\n";
print "Output_file $outputfile\n";

my $line1 = "";
my $line2 = "";
my $line3 = "";
my $line4 = "";
my $lineid = 0;
my $length = "";
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
		$length = length($line);
	}
	elsif ($lineid == 3) {
		$line3 = $line;
	}
	elsif ($lineid == 4) {
		$line4 = $line;
		if ($length >= $min_length) {
			print OUTPUT "$line1\n$line2\n$line3\n$line4\n";
		}
		$line1 = "";
		$line2 = "";
		$line3 = "";
		$line4 = "";
		$lineid = 0;
		$length = "";
	}
	else {
		print "Error_line\n";
	}
}
close (IN);
