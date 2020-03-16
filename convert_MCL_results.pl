
open (my $inFile, '<', 'SW.txt') or die $!;
while (<$inFile>) {
  push(@array,split /\n+/);
  #print "@array\n";
}
close ($inFile);
open OUT, ">MCL-groups.txt";
$group==0;
my @headers = split('\s', $array[0]);
for ($n=1;$n<=666;$n++){
	@gec=();
	@gec=split('\s',$array[$n]);
	$group++;
	foreach $el (0..$#gec)
		{
			if ($gec[$el]==1)
			{
				print OUT "$headers[$el]\t $group\n";
			}
			else{}
		}
}