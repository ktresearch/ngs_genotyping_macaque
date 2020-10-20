#!/usr/bin/perl

use strict;

my $inputfile = shift @ARGV;
my $outputfile = shift @ARGV;

open (IN, "$inputfile") || die $!;
open (OUTPUT, ">$outputfile");
my $check = 0;
while (my $line = <IN>) {	
	chomp $line;
	$check += 1;
	if ($check == 1) {
		print OUTPUT ">$line\:Length_";
	}
	elsif ($check == 2) {
		my $length = length($line);		
		print OUTPUT "$length\n$line\n";
	}
	elsif ($check == 3) {
		next;
	}
	elsif ($check == 4) {
		$check = 0;
		next;
	}
	else {
		print "error\n";
	}
}
close (IN);
close (OUTPUT);
