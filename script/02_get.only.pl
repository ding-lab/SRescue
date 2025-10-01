$inBEDPE=$ARGV[0];
$inVCF=$ARGV[1];
$outdir=$ARGV[2];
open(file,"$inBEDPE");
open(res1,">$outdir/sronly.tsv");
open(res2,">$outdir/lronly.tsv");
open(res3,">$outdir/shared.tsv");
while(<file>)
{
	chomp;
	@token=split(/\t/,$_);
	# print $token[9]."\t".$token[10]."\n";
	if($token[9] ne "FALSE" and $token[10] eq "FALSE")
	{
		print res1 $token[0]."\t".$token[1]."\t".$token[7]."\t".$token[9]."\t".$token[6]."\t".$token[3]."\t".$token[4]."\n";
		$sronly{$token[0]."\t".$token[1]}=$token[7];
		$sronly{$token[3]."\t".$token[4]}=$token[7];
	}
	if($token[10] ne "FALSE" and $token[9] eq "FALSE")
	{
		print res2 $token[0]."\t".$token[1]."\t".$token[7]."\t".$token[10]."\t".$token[6]."\t".$token[3]."\t".$token[4]."\n";
		$lronly{$token[0]."\t".$token[1]}=$token[7];
		$lronly{$token[3]."\t".$token[4]}=$token[7];
	}
	if($token[10] ne "FALSE" and $token[9] ne "FALSE")
	{
		print res3 $token[0]."\t".$token[1]."\t".$token[7]."\t".$token[9].";".$token[10]."\t".$token[6]."\t".$token[3]."\t".$token[4]."\n";
		$shared{$token[0]."\t".$token[1]}=$token[7];
		$shared{$token[3]."\t".$token[4]}=$token[7];
	}
}
$file=$inVCF;
#$file=~s/bedpe$/vcf/g;	# specify paths explicitly
# print $file."\n";
open(file,"$file");
open(res4,">$outdir/sronly.vcf");
open(res5,">$outdir/lronly.vcf");
open(res6,">$outdir/shared.vcf");
while(<file>)
{
	chomp;
	if(grep(/^chr/,$_))
	{
		@token=split(/\t/,$_);
		# print $_."\n";
		if(exists($sronly{$token[0]."\t".$token[1]}))
		{
			# print res1 $token[0]."\t".$token[1]."\t".$sronly{$token[0]."\t".$token[1]}."\t".$token[2]."\n";
			print res4 $_."\n";
		}
		if(exists($shared{$token[0]."\t".$token[1]}) or exists($shared{$token[3]."\t".$token[4]}))
		{
			# print res3 $token[0]."\t".$token[1]."\t".$shared{$token[0]."\t".$token[1]}."\t".$token[2]."\n";
			print res6 $_."\n";
		}
		if(exists($lronly{$token[0]."\t".$token[1]}))
		{
			# print res2 $token[0]."\t".$token[1]."\t".$lronly{$token[0]."\t".$token[1]}."\t".$token[2]."\n";
			print res5 $_."\n";
		}
	}else{
		print res4 $_."\n";
		print res5 $_."\n";
		print res6 $_."\n";
	}
}
close(res1);
close(res2);
close(res3);
close(res4);
close(res5);
close(res6);
close(file);
