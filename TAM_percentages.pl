
use strict;
use warnings;
# parameters
my $pvalue_Control = 0.005;
my $percentage_Control = 0.2;
my $SumClust=0;
my $SumFunc=0;
my $SumFam=0;
my $SumHM=0;
my $SumTissue=0;
my $SumAll =0;
my $output_file="output.txt";
open(OUTPUT,">".$output_file);
my $directory = "C:/strawberry/perl/bin/4txtdir/";
opendir DIR, $directory ;
my @files = grep{ /\.txt/} readdir (DIR) ;
closedir DIR;
my $groupnumber=$#files+1;
foreach my $file (@files)
{
	my (@fields,@mydata);
	
    open (IN, $directory.$file);
	my $headline= <IN>;
    while (<IN>)


    {
        chomp;
        @fields = split (/\,/,$_) ;
		
		if ($#fields==7)
		{
			push(@mydata,@fields);
		}
		else 
		{
			my $val2 = $#fields - 7;
			splice @fields,1, $val2;
		}
    }
    close(IN);
	my (@category, @term , @count, @percent, @fold ,@pvalue, @bonferroni, @FDR);
	my $hitClu=0;
    my $hitFunc=0;
    my $hitFam=0; 
	my $hitHM =0;
	my $hitTissue =0;	
	for (my $i=0; $i <= $#mydata; $i+=8) { push (@category, $mydata[$i]) }; 
	for (my $i=1; $i <= $#mydata; $i+=8) { push (@term, $mydata[$i]) }; 
	for (my $i=2; $i <= $#mydata; $i+=8) { push (@count, $mydata[$i]) }; 
	for (my $i=3; $i <= $#mydata; $i+=8) { push (@percent, $mydata[$i]) }; 
	for (my $i=4; $i <= $#mydata; $i+=8) { push (@fold, $mydata[$i]) }; 
	for (my $i=5; $i <= $#mydata; $i+=8) { push (@pvalue, $mydata[$i]) }; 
	for (my $i=6; $i <= $#mydata; $i+=8) { push (@bonferroni, $mydata[$i]) }; 
	for (my $i=7; $i <= $#mydata; $i+=8) { push (@FDR, $mydata[$i]) }; 
	for my $e(0..$#category)
	{
		if ($category[$e] eq "Cluster")
		{
			if ( $count[$e]!=1 && $pvalue[$e] < $pvalue_Control && ($percent[$e]*$count[$e])>$percentage_Control)
			{
				$hitClu=1;
			}
		}
		if ($category[$e] eq "Function")
		{
			if ( $count[$e]!=1 && $pvalue[$e]<$pvalue_Control && ($percent[$e]*$count[$e])>$percentage_Control)
			{
				$hitFunc=1;
			}	
		}
		if ($category[$e] eq "Family")
		{
			if ( $count[$e]!=1 && $pvalue[$e]<$pvalue_Control && ($percent[$e]*$count[$e])>$percentage_Control)
			{
				$hitFam=1;
			}
		}
		if ($category[$e] eq "HMDD")
		{
			if ( $count[$e]!=1 && $pvalue[$e] < $pvalue_Control && ($percent[$e]*$count[$e]) > $percentage_Control)
			{
				$hitHM=1;
			}
		}
		if ($category[$e] eq "TissueSpecific")
		{
			if ( $count[$e]!=1 && $pvalue[$e] < $pvalue_Control && ($percent[$e]*$count[$e]) > $percentage_Control)
			{
				$hitTissue=1;
			}
		}
	}
		if ($hitClu == 1){ $SumClust++;}
		if ($hitFunc == 1){ $SumFunc++;}
		if ($hitFam == 1){ $SumFam++;}
		if ($hitHM == 1){ $SumHM++;}
		if ($hitTissue == 1){ $SumTissue++;}
		if ($hitHM == 1 || $hitClu == 1 || $hitFunc == 1 || $hitFam == 1 || $hitTissue == 1 ){ $SumAll++;}
}


my $CluP = ($SumClust*100/$groupnumber);
my $FuncP = ($SumFunc*100/$groupnumber);
my $FamP = ($SumFam*100/$groupnumber);
my $HMP = ($SumHM*100/$groupnumber);
my $TissueP= ($SumTissue*100/$groupnumber);
my $AllP = ($SumAll*100/$groupnumber);


print OUTPUT "Clust Per = ". $CluP. "\n";
print OUTPUT "Function Per =" .$FuncP."\n";
print OUTPUT "Family Per = ".$FamP."\n";
print OUTPUT "HMDD Per = ".$HMP."\n";
print OUTPUT "Tissue Per = ".$TissueP."\n";
print OUTPUT "All Per = ".$AllP."\n";
