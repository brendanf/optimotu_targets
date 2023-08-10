#!/usr/bin/perl

$seq="";
while (<>) {
    if (/^>/) {
	print $id . $seq . "\n" if ($seq ne "");
	$id = $_;
	$seq = "";
    }
    else {
	s/[a-z\s]//g;
	$seq = $seq . $_;
    }
}

print $id . $seq . "\n" if ($seq ne "");
