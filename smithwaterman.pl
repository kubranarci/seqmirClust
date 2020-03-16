
use strict;
use warnings;
my $output_file="output.txt";
open(OUTPUT,">".$output_file);
my @seq_array;
my @sequences;
my @headers;


#read the text file, separate sequences and headers into 2 arrays
open (my $inFile, '<', 'hsa-fasta.txt') or die $!;
while (<$inFile>) 
{
  push(@seq_array,split /\s+/);
}
close ($inFile);
for (my$i=0; $i <= $#seq_array-1; $i+=2) { push (@headers, $seq_array[$i]) }; 
for (my$i=1; $i <= $#seq_array+1; $i+=2) { push (@sequences, $seq_array[$i]) };
print OUTPUT2, @headers;
my @blastmatrix=();
for (my $i=0; $i< $#sequences+1;$i++)
{
	for (my $j=0; $j< $#sequences+1;$j++)
	{
	$blastmatrix[$i][$j]='';
	}
}
# call blast and fill the matrix blast with the scores
for (my $e=0;$e<$#sequences+1;$e++) 
{ 
    for (my $k=0;$k<$#sequences+1;$k++)
	{
	$blastmatrix[$e][$k]= makeblast($sequences[$e],$sequences[$k]);
	}
}
###subroutine blast############
my $seq1;
my $seq2;


sub makeblast
{
	$seq1=shift;
	$seq2=shift;
	# scoring scheme
	my $MATCH    =  1; # +1 for letters that match
	my $MISMATCH = -1; # -1 for letters that mismatch
	my $GAP      = -1; # -1 for any gap


	# initialization
	my @matrix;
	$matrix[0][0]{score}   = 0;
	$matrix[0][0]{pointer} = "none";
	for(my $j = 1; $j <= length($seq1); $j++) 
	{
		$matrix[0][$j]{score}   = 0 * $j;
		$matrix[0][$j]{pointer} = "none";
	}
	for (my $i = 1; $i <= length($seq2); $i++) 
	{
		$matrix[$i][0]{score}   = 0 * $i;
		$matrix[$i][0]{pointer} = "none";
	}
	# fill
	my $max_i     = 0;
	my $max_j     = 0;
	my $max_score = 0;
	for(my $i = 1; $i <= length($seq2); $i++) 
	{
		for(my $j = 1; $j <= length($seq1); $j++) 
		{
			my ($diagonal_score, $left_score, $up_score);
			# calculate match score
			my $letter1 = substr($seq1, $j-1, 1);
			my $letter2 = substr($seq2, $i-1, 1);      
			if ($letter1 eq $letter2) 
			{
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
			}
			else 
			{
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
			}
			# calculate gap scores
			$up_score   = $matrix[$i-1][$j]{score} + $GAP;
			$left_score = $matrix[$i][$j-1]{score} + $GAP;
			if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) 
			{
				$matrix[$i][$j]{score}   = 0;
				$matrix[$i][$j]{pointer} = "none";
				next; # terminate this iteration of the loop
			}
			# choose best score
			if ($diagonal_score >= $up_score) 
			{
				if ($diagonal_score >= $left_score) 
				{
					$matrix[$i][$j]{score}   = $diagonal_score;
					$matrix[$i][$j]{pointer} = "diagonal";
				}
				else 
				{
					$matrix[$i][$j]{score}   = $left_score;
					$matrix[$i][$j]{pointer} = "left";
				}
			} 
			else 
			{
				if ($up_score >= $left_score) 
				{
					$matrix[$i][$j]{score}   = $up_score;
					$matrix[$i][$j]{pointer} = "up";
				}
				else 
				{
					$matrix[$i][$j]{score}   = $left_score;
					$matrix[$i][$j]{pointer} = "left";
				}
			}
			# set maximum score
			if ($matrix[$i][$j]{score} > $max_score) 
			{
				$max_i     = $i;
				$max_j     = $j;
				$max_score = $matrix[$i][$j]{score};
			}
		}
	}
	# trace-back
	my $align1 = "";
	my $align2 = "";
	# start at last cell of matrix
	my $j = length($seq1);
	my $i = length($seq2);
	while (1) 
	{
		last if $matrix[$i][$j]{pointer} eq "none"; # ends at first cell of matrix
		if ($matrix[$i][$j]{pointer} eq "diagonal") 
		{
			$align1 .= substr($seq1, $j-1, 1);
			$align2 .= substr($seq2, $i-1, 1);
			$i--;
			$j--;
		}
		elsif ($matrix[$i][$j]{pointer} eq "left") 
		{
			$align1 .= substr($seq1, $j-1, 1);
			$align2 .= "-";
			$j--;
		}
		elsif ($matrix[$i][$j]{pointer} eq "up") 
		{
			$align1 .= "-";
			$align2 .= substr($seq2, $i-1, 1);
			$i--;
		}    
	}
	my $ak = length($seq1);
	my $kara = length($seq2);
	return $matrix[$kara][$ak]{score};
}
#write the distance/blast matrix on outfile 
$"="\t";
print OUTPUT "matrix\t@headers\n";
for(my $i = 0; $i < $#sequences+1; $i++)
{
	# $#array_2d gives the highest index from the array
	print OUTPUT $headers[$i]."\t";
	for(my $j = 0; $j < $#sequences+1 ; $j++)
	{
		print OUTPUT "$blastmatrix[$i][$j]\t" ;
	}
	print OUTPUT "\n"; 
}
print "$#blastmatrix\n";