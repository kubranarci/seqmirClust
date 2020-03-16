
use strict;
use warnings;
my @seq_array;
my @sequences;
my @headers;
my $output_file="output.txt";
open(OUTPUT,">".$output_file);
open (my $inFile, '<', 'hsa-fasta.txt') or die $!;
while (<$inFile>) 
{
  push(@seq_array,split /\s+/);
}
close ($inFile);
for (my$i=1; $i < $#seq_array+1; $i+=2) { push (@sequences, $seq_array[$i]) };
for (my$i=0; $i < $#seq_array; $i+=2) { push (@headers, $seq_array[$i]) }; 
my @kmer_array= qw/AAA	AAT	AAC	AAG	ATA	ATT	ATC	ATG	ACA	ACT	ACC	ACG	AGA	AGT	AGC
	AGG	TAA	TAT	TAC	TAG	TTA	TTT	TTC	TTG	TCA	TCT	TCC	TCG	TGA	TGT	TGC	TGG	TAA	CAT	
	CAC	CAG	CTA	CTT	CTC	CTG	CCA	CCT	CCC	CCG	CGA	CGT	CGC	CGG	GAA	GAT	GAC	GAG	GTA	
	GTT	GTC	GTG	GCA	GCT	GCC	GCG	GGA	GGT	GGC	GGG/;
my @kmer_matrix;
$kmer_matrix[0][0]{score}= 0;


for my $i (0..$#sequences){
for my $j (0..$#kmer_array){


	if ( substr($sequences[$i],0)=~/$kmer_array[$j]/)
	{
		$kmer_matrix[$i][$j]{score}=1;
	}
	else 
	{
		$kmer_matrix[$i][$j]{score}=0;
	}
}
}
 no warnings 'uninitialized';
#write the distance/blast matrix on outfile 
$"="\t";
print OUTPUT "matrix\t@kmer_array\n";
for(my $i = 0; $i < $#headers+1; $i++)
{
   print OUTPUT "$headers[$i]\t";
   for(my $j = 0; $j < 64 ; $j++)
   {
      print OUTPUT "$kmer_matrix[$i][$j]{score}\t" ;
   }
   print OUTPUT "\n";
