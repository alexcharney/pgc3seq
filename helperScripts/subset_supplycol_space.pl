#!/usr/local/perl


my %subset = ();

my $list = $ARGV[0];
my $bigfile = $ARGV[1];
my $col =  $ARGV[2];
#my $output = $ARGV[3]; 

open(list, $list); 
while(<list>){
    chomp($_);
    $subset{$_} = 1;
}

#open(out, ">", $output);
open(big, $bigfile);
while(<big>){
    chomp($_);
    my @line = split(" ", $_);
    if(exists $subset{$line[$col]}){ 
#        print out "$_\n";
        print "$_\n";
    } 
}
